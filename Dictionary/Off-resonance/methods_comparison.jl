using KomaMRI, MAT, Suppressor, LinearAlgebra, CUDA, Distributions, Statistics, CSV, DataFrames

T2prime_list_ms = [25.0, 50.0, 75.0, 100.0]
phantom_length = 10.0
num_points = 2001
seq_file = "sequences/mpf_001_PhantomStudy_short_124.seq"

N_list = [3, 5, 7, 9, 11, 15, 21, 31]
SUBSET_SIZE = 12
MAX_SPINS_PER_CALL = 225_000

OUT_BASENAME = "offres_convergence"
OUT_CSV_ABS = OUT_BASENAME * "_ABS.csv"

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx_tbl = read(f_idx, "idx"); close(f_idx)
sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx_tbl)]
sampled_pairs_s = [(Float64(T1)/1000, Float64(T2)/1000) for (T1, T2) in sampled_pairs_ms]

function pick_subset(pairs::Vector{<:Tuple}, k::Int)
    n = length(pairs)
    k >= n && return pairs
    idxs = unique(round.(Int, range(1, n; length=k)))
    return pairs[idxs]
end

pairs_subset = pick_subset(sampled_pairs_s, SUBSET_SIZE)
keys_subset = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_subset]

phantom_length_m = phantom_length / 1000
pos = collect(range(-phantom_length_m/2, phantom_length_m/2, length=num_points))
zN = zeros(Float64, num_points)
xvec, yvec, zvec = (pos, zN, zN)

sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true
seq = read_seq(seq_file)

function build_phantom_single_pair_freq(T1::Real, T2::Real, Δω_values::AbstractVector,
                                        xvec::AbstractVector, yvec::AbstractVector, zvec::AbstractVector)
    npos = length(xvec)
    nfreq = length(Δω_values)
    N = npos * nfreq
    x = Vector{Float64}(undef, N)
    y = Vector{Float64}(undef, N)
    z = Vector{Float64}(undef, N)
    T1v = Vector{Float64}(undef, N)
    T2v = Vector{Float64}(undef, N)
    Δwv = Vector{Float64}(undef, N)
    idx = 1
    T1f = Float64(T1); T2f = Float64(T2)
    @inbounds for Δω in Δω_values
        Δωf = Float64(Δω)
        for i in 1:npos
            x[idx] = xvec[i]; y[idx] = yvec[i]; z[idx] = zvec[i]
            T1v[idx] = T1f; T2v[idx] = T2f; Δwv[idx] = Δωf
            idx += 1
        end
    end
    return Phantom{Float64}(x=x, y=y, z=z, T1=T1v, T2=T2v, Δw=Δwv)
end

function gauss_legendre_cauchy(γ_hz::Float64; N::Int)
    βtail = 0.5 ./ sqrt.(1 .- (2 .* (1:N-1)) .^ (-2))
    T = SymTridiagonal(zeros(N), βtail)
    evals, evecs = eigen(T)
    û = evals
    ŵ = 2 .* (evecs[1, :]).^2
    f_hz = γ_hz .* tan.(0.5π .* û)
    Δω = 2π .* f_hz
    w_norm = 0.5 .* ŵ
    return Δω, w_norm
end

function make_delta_omega_quantile(N::Int, T2prime_s::Float64)
    γ_hz = 1 / (2π * T2prime_s)
    lorentz = Cauchy(0.0, γ_hz)
    u = (collect(1:N) .- 0.5) ./ N
    Δf = quantile.(Ref(lorentz), u)
    Δω = 2π .* Δf
    w = fill(1.0/N, N)
    return Δω, w
end

function make_delta_omega(method::Symbol, N::Int, T2prime_s::Float64)
    if method === :quantile
        return make_delta_omega_quantile(N, T2prime_s)
    elseif method === :gauss
        γ_hz = 1 / (2π * T2prime_s)
        return gauss_legendre_cauchy(γ_hz; N=N)
    end
end

rel_diff_per_entry(a::AbstractVector, b::AbstractVector) = (norm(b) == 0 ? 0.0 : norm(a .- b) / norm(b))

function simulate_for(method::Symbol, N::Int; T2prime_ms::Float64)
    T2prime_s = T2prime_ms / 1000
    Δω_values, w_freq = make_delta_omega(method, N, T2prime_s)
    Npos = num_points
    safe_freq_batch = max(1, min(N, fld(MAX_SPINS_PER_CALL, Npos)))
    ranges = (safe_freq_batch >= N) ? [1:N] :
             [((i-1)*safe_freq_batch+1):min(i*safe_freq_batch, N) for i in 1:cld(N, safe_freq_batch)]
    out = Dict{Tuple{Int,Int}, Vector{ComplexF32}}()
    for (pair, key) in zip(pairs_subset, keys_subset)
        acc = ComplexF32[]
        for fr in ranges
            Δω_chunk = @view Δω_values[fr]
            w_chunk = @view w_freq[fr]
            big_phantom = build_phantom_single_pair_freq(pair[1], pair[2], Δω_chunk, xvec, yvec, zvec)
            sig_all = @suppress simulate(big_phantom, seq, sys; sim_params=sim_params)
            sig_clean = dropdims(sig_all; dims=(3,4))
            w_spins_c = repeat(w_chunk, inner=Npos) ./ Npos
            contrib = sig_clean * w_spins_c
            if isempty(acc)
                acc = ComplexF32.(contrib)
            else
                acc .+= contrib
            end
        end
        out[key] = acc
    end
    return out
end

function abs_error_vs_ref(method::Symbol, N::Int, ref_dict; T2prime_ms::Float64)
    cand = simulate_for(method, N; T2prime_ms=T2prime_ms)
    rels = [rel_diff_per_entry(cand[key], ref_dict[key]) for key in keys_subset]
    return (method=String(Symbol(method)), N=N, T2prime_ms=T2prime_ms,
            median_pct=100*median(rels), p95_pct=100*quantile(rels, 0.95), max_pct=100*maximum(rels))
end

const METHODS = [:gauss, :quantile]
const N_REF = 201

for T2p in T2prime_list_ms
    ref_dict = simulate_for(:gauss, N_REF; T2prime_ms=T2p)
    for m in METHODS, N in N_list
        row = abs_error_vs_ref(m, N, ref_dict; T2prime_ms=T2p)
        df = DataFrame([row])
        if !isfile(OUT_CSV_ABS)
            CSV.write(OUT_CSV_ABS, df)
        else
            CSV.write(OUT_CSV_ABS, df; append=true)
        end
    end
end

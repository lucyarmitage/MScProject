using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra, CUDA, Distributions

T2prime_ms = parse(Float64, ARGS[1])
batch_size_pairs = parse(Int, ARGS[2])
T2prime_s = T2prime_ms / 1000.0

phantom_length = 8
num_points = 4001
sliceOrientation = 1
N_freq_samples = 15
coverage = 0.99
seq_file = "sequences/mpf_001_PhantomStudy_124.seq"

phantom_length_m = phantom_length / 1000
pos = collect(range(-phantom_length_m/2, phantom_length_m/2, length=num_points))
zN = zeros(Float64, num_points)
xvec, yvec, zvec = if sliceOrientation == 1
    (pos, zN, zN)
elseif sliceOrientation == 2
    (zN, zN, pos)
elseif sliceOrientation == 3
    (zN, pos, zN)
else
    error("sliceOrientation must be 1, 2, or 3")
end

function build_combined_phantom(pairs_chunk::AbstractVector{<:Tuple{<:Real,<:Real}}, Δω_values::AbstractVector{<:Real}, xvec::Vector{Float64}, yvec::Vector{Float64}, zvec::Vector{Float64})
    npos = length(xvec)
    nfreq = length(Δω_values)
    N = length(pairs_chunk) * npos * nfreq
    x = Vector{Float64}(undef, N)
    y = Vector{Float64}(undef, N)
    z = Vector{Float64}(undef, N)
    T1v = Vector{Float64}(undef, N)
    T2v = Vector{Float64}(undef, N)
    Δwv = Vector{Float64}(undef, N)
    idx = 1
    for (T1, T2) in pairs_chunk
        T1 = Float64(T1); T2 = Float64(T2)
        for Δω in Δω_values
            @inbounds for i in 1:npos
                x[idx] = xvec[i]
                y[idx] = yvec[i]
                z[idx] = zvec[i]
                T1v[idx] = T1
                T2v[idx] = T2
                Δwv[idx] = Δω
                idx += 1
            end
        end
    end
    return Phantom{Float64}(x=x, y=y, z=z, T1=T1v, T2=T2v, Δw=Δwv)
end

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx"); close(f_idx)
sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s = [(T1/1000, T2/1000) for (T1, T2) in sampled_pairs_ms]

sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true
seq = read_seq(seq_file)

function gauss_legendre_cauchy(γ_hz::Float64; N::Int, coverage::Float64)
    α = (1 - coverage) / 2
    a = π*(α - 0.5)
    b = -a
    β = [0.0; 0.5 ./ sqrt.(1 .- (2 .* (1:N-1)) .^ (-2))]
    T = SymTridiagonal(zeros(N), β[2:end])
    evals, evecs = eigen(T)
    t̂ = evals
    ŵ = 2 .* (evecs[1, :]) .^ 2
    t = 0.5*(b-a) .* t̂ .+ 0.5*(b+a)
    w_interval = 0.5*(b-a) .* ŵ
    w_t = w_interval ./ (b - a)   # integrates to ~1 over [a,b]
    f_hz = γ_hz .* tan.(t)
    Δω = 2π .* f_hz    # rad/s
    return Δω, w_t
end

γ_hz = 1 / (2π * T2prime_s)
Δω_values, w_freq = gauss_legendre_cauchy(γ_hz; N=N_freq_samples, coverage=coverage)

out_file = "dict/prostate_$(phantom_length)mm_$(num_points)_$(T2prime_ms)ms_short.mat"
batch_results = Dict{Tuple{Int, Int}, Vector{ComplexF32}}()

for chunk_start in 1:batch_size_pairs:length(sampled_pairs_s)
    chunk_end = min(chunk_start + batch_size_pairs - 1, length(sampled_pairs_s))
    pairs_chunk = sampled_pairs_s[chunk_start:chunk_end]
    keys_chunk = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_chunk]
    big_phantom = build_combined_phantom(pairs_chunk, Δω_values, xvec, yvec, zvec)
    sig_all = @suppress simulate(big_phantom, seq, sys; sim_params=sim_params)
    sig_clean_all = dropdims(sig_all; dims=(3,4))
    spins_per_pair = num_points * length(Δω_values)
    w_spins = repeat(w_freq, inner=num_points) ./ num_points  # frequency weights × uniform space average; sums to 1
    start_idx = 1
    for key in keys_chunk
        sig_this = sig_clean_all[:, start_idx:start_idx + spins_per_pair - 1] 
        batch_results[key] = ComplexF32.(sig_this * w_spins)                    # weighted complex average over Δω and space
        start_idx += spins_per_pair
    end
    println("T2'=$(T2prime_ms) ms | Processed chunk $chunk_start:$chunk_end")
end

files_keys = sort(collect(keys(batch_results)))
num_entries = length(files_keys)
timepoints = size(first(values(batch_results)), 1)
bloch_matrix = zeros(ComplexF32, timepoints, num_entries)
idx_bloch = zeros(Float32, num_entries, 2)
for (j, key) in enumerate(files_keys)
    bloch_matrix[:, j] .= batch_results[key]
    idx_bloch[j, :] .= key
end
mkpath(dirname(out_file))
matwrite(out_file, Dict("dict0" => bloch_matrix, "idx" => idx_bloch))
println("Saved → $out_file")

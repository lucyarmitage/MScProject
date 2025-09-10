
using KomaMRI, MAT, Suppressor, FileIO, LinearAlgebra, CUDA

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
    @inbounds for (T1_raw, T2_raw) in pairs_chunk
        T1 = Float64(T1_raw)
        T2 = Float64(T2_raw)
        for Δω in Δω_values
            for i in 1:npos
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

function main()
    T2prime_ms = parse(Float64, ARGS[1])
    T2prime_s = T2prime_ms / 1000.0

    phantom_length = 8
    num_points = 2001
    sliceOrientation = 1
    N_freq_samples = 9
    seq_file = "sequences/mpf_001_PhantomStudy_short_124.seq"
    batch_size_pairs = 112

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

    f_idx = matopen("D_IDX_SP_Phantom2025.mat")
    idx_tbl = read(f_idx, "idx"); close(f_idx)
    sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx_tbl)]
    sampled_pairs_s = [(T1/1000, T2/1000) for (T1, T2) in sampled_pairs_ms]
    total_entries = length(sampled_pairs_s)

    sys = Scanner()
    sim_params = KomaMRICore.default_sim_params()
    sim_params["return_type"] = "mat"
    sim_params["sim_method"] = BlochDict()
    sim_params["gpu"] = true
    seq = read_seq(seq_file)

    γ_hz = 1 / (2π * T2prime_s)
    Δω_values, w_freq = gauss_legendre_cauchy(γ_hz; N=N_freq_samples)

    out_file = "dict/$(phantom_length)mm_$(num_points)_$(T2prime_ms)ms_short.mat"
    mkpath(dirname(out_file))

    bloch_matrix = Array{ComplexF32}(undef, 0, 0)
    idx_bloch = zeros(Float32, total_entries, 2)
    spins_per_pair = num_points * length(Δω_values)

    for chunk_start in 1:batch_size_pairs:total_entries
        chunk_end = min(chunk_start + batch_size_pairs - 1, total_entries)
        pairs_chunk = sampled_pairs_s[Int(chunk_start):Int(chunk_end)]
        keys_chunk = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_chunk]
        big_phantom = build_combined_phantom(pairs_chunk, Δω_values, xvec, yvec, zvec)
        sig_all = @suppress simulate(big_phantom, seq, sys; sim_params=sim_params)
        sig_clean_all = Array(dropdims(sig_all; dims=(3,4)))
        if isempty(bloch_matrix)
            T = size(sig_clean_all, 1)
            bloch_matrix = zeros(ComplexF32, T, total_entries)
        end
        w_spins = repeat(w_freq, inner=num_points) ./ num_points
        start_idx = 1
        for (j, key) in enumerate(keys_chunk)
            col = chunk_start + j - 1
            sig_this = @view sig_clean_all[:, start_idx:start_idx + spins_per_pair - 1]
            bloch_matrix[:, col] .= ComplexF32.(sig_this * w_spins)
            idx_bloch[col, 1] = Float32(key[1])
            idx_bloch[col, 2] = Float32(key[2])
            start_idx += spins_per_pair
        end
    end

    matwrite(out_file, Dict("dict0"=>bloch_matrix, "idx"=>idx_bloch))
end

main()

using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra, Distributions, CUDA, Printf

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx")
close(f_idx)
sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s  = [(T1 / 1000, T2 / 1000) for (T1, T2) in sampled_pairs_ms]

sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true
seq = read_seq("mpf_001_PhantomStudy_short_124.seq")

x_pos = collect(range(-7.5e-3, 7.5e-3, length=151))
batch_size_pairs = 250
timepoints = 1000

function build_batch_phantom(pairs_chunk, Δf_values, weights)
    big_phantom = Phantom{Float64}(name = "batch", x = Float64[], y = Float64[], z = Float64[], T1 = Float64[], T2 = Float64[])
    spins_per_pair = Int[]
    for (T1, T2) in pairs_chunk
        total_spins = 0
        for (Δf, w) in zip(Δf_values, weights)
            Δω = 2π * Δf
            aux = Phantom{Float64}(
                x = x_pos, y = zeros(length(x_pos)), z = zeros(length(x_pos)),
                T1 = fill(T1, length(x_pos)), T2 = fill(T2, length(x_pos)),
                Δw = fill(Δω, length(x_pos)), ρ = fill(w, length(x_pos))
            )
            big_phantom += aux
            total_spins += length(x_pos)
        end
        push!(spins_per_pair, total_spins)
    end
    return big_phantom, spins_per_pair
end

for γ in [1.0, 2.0, 3.0, 4.0]
    gamma_label = replace(@sprintf("%.1f", γ), "." => "_")
    N_samples = 15
    cutoff = 3 * γ
    Δf_values = collect(range(-cutoff, cutoff; length=N_samples))
    lorentz_pdf(f) = (γ / π) ./ (f.^2 .+ γ^2)
    weights = lorentz_pdf.(Δf_values)
    weights ./= sum(weights)

    out_folder = joinpath("/vol/bitbucket/la724/progress", "15mm_$gamma_label")
    mkpath(out_folder)
    batch_results = Dict{Tuple{Int, Int}, Vector{ComplexF32}}()

    for chunk_start in 1:batch_size_pairs:length(sampled_pairs_s)
        chunk_end = min(chunk_start + batch_size_pairs - 1, length(sampled_pairs_s))
        pairs_chunk = sampled_pairs_s[chunk_start:chunk_end]
        keys_chunk = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_chunk]

        if all(isfile(joinpath(out_folder, "signal_$(k[1])_$(k[2]).jld2")) for k in keys_chunk)
            continue
        end

        big_phantom, spins_per_pair = build_batch_phantom(pairs_chunk, Δf_values, weights)
        sig_all = nothing
        @suppress sig_all = simulate(big_phantom, seq, sys; sim_params=sim_params)
        sig_clean_all = dropdims(sig_all; dims=(3,4))

        start_idx = 1
        for (i, key) in enumerate(keys_chunk)
            stop_idx = start_idx + spins_per_pair[i] - 1
            sig_this = sig_clean_all[:, start_idx:stop_idx]
            sig_sum = vec(sum(sig_this, dims=2))
            batch_results[key] = sig_sum
            start_idx = stop_idx + 1
        end

        for (key, sig) in batch_results
            filepath = joinpath(out_folder, "signal_$(key[1])_$(key[2]).jld2")
            @save filepath key signal_mag=sig
        end
        empty!(batch_results)
    end

    files = readdir(out_folder)
    num_entries = length(files)
    bloch_matrix = zeros(ComplexF32, timepoints, num_entries)
    idx_bloch = zeros(ComplexF32, num_entries, 2)

    for (j, file) in enumerate(files)
        filepath = joinpath(out_folder, file)
        @load filepath key signal_mag
        bloch_matrix[:, j] .= signal_mag
        idx_bloch[j, :] .= key
    end

    mkpath("dict")
    out_file = "dict/dict_15mm_$gamma_label.mat"
    matwrite(out_file, Dict("dict0" => bloch_matrix, "idx" => idx_bloch))
end
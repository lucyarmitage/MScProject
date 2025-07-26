using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra, Printf

seq_folder = joinpath(@__DIR__, "b1_sequences")
out_folder = joinpath("progress", "progress_10mm_101_short_b1")
out_file = "blochdict_10mm_101_short_b1.mat"
mkpath(out_folder)

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx")
close(f_idx)

sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s  = [(T1 / 1000, T2 / 1000) for (T1, T2) in sampled_pairs_ms]

sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = false

x_pos = collect(range(-5e-3, 5e-3, length=101))

batch = Dict{Tuple{Int, Int, Float64}, Vector{ComplexF32}}()
batch_size = 50

seq_files = filter(f -> endswith(f, ".seq"), readdir(seq_folder))

for seq_file in seq_files
    m = match(r"b1_(\d+\.\d+)", seq_file)
    b1 = m === nothing ? 1.0 : parse(Float64, m.captures[1])
    rounded_b1 = round(b1; digits=2)

    println("Processing B1 = $rounded_b1 from $seq_file")
    flush(stdout)

    seq = read_seq(joinpath(seq_folder, seq_file))

    for (i, (T1, T2)) in enumerate(sampled_pairs_s)
        key = (Int(round(T1 * 1000)), Int(round(T2 * 1000)), rounded_b1)
        filepath = joinpath(out_folder, "signal_$(key[1])_$(key[2])_$(replace(string(rounded_b1), "." => "p")).jld2")

        if isfile(filepath)
            if i % 50 == 0
                println("Skipping existing signal for $key (iteration $i)")
                flush(stdout)
            end
            continue
        end

        obj = Phantom{Float64}(
            x = x_pos,
            y = zeros(length(x_pos)),
            z = zeros(length(x_pos)),
            T1 = fill(T1, length(x_pos)),
            T2 = fill(T2, length(x_pos))
        )

        sig = nothing
        @suppress sig = simulate(obj, seq, sys; sim_params=sim_params)

        sig_clean = dropdims(sig; dims=(3,4))
        signal_mag = vec(sum(sig_clean, dims=2))

        batch[key] = signal_mag

        if length(batch) >= batch_size || i == length(sampled_pairs_s)
            for (key, sig) in batch
                b1_str = replace(string(key[3]), "." => "p")
                filepath = joinpath(out_folder, "signal_$(key[1])_$(key[2])_$(b1_str).jld2")
                @save filepath key signal_mag=sig
            end
            println("Saved batch at iteration $i (size: $(length(batch))) for B1 = $rounded_b1")
            flush(stdout)
            empty!(batch)
        end
    end
end

files = filter(f -> endswith(f, ".jld2"), readdir(out_folder))
num_entries = length(files)
timepoints = 1000

bloch_matrix = zeros(ComplexF32, timepoints, num_entries)
idx_bloch = zeros(Float64, num_entries, 3)  # [T1, T2, B1]

for (j, file) in enumerate(files)
    filepath = joinpath(out_folder, file)
    @load filepath key signal_mag
    bloch_matrix[:, j] .= signal_mag
    idx_bloch[j, 1] = key[1]
    idx_bloch[j, 2] = key[2]
    idx_bloch[j, 3] = key[3]
end

matwrite(out_file, Dict("dict0" => bloch_matrix, "idx" => idx_bloch))
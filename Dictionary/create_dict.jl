using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra
# CUDA.set_runtime_version!(v"12.0")
# println("CUDA available: ", CUDA.has_cuda())
# println("GPU name: ", CUDA.name(CUDA.device()))
# CUDA.allowscalar(false)

out_folder = joinpath("progress", "progress_20mm_101")
out_file   = "blochdict_20mm_101.mat"

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx")
close(f_idx)

f_dict = matopen("D_Phantom2025_invEff096_SPinf_norefocusingTEadj_576InvTime_1000RF_10mm_101iso_0.mat")
dict_epg = read(f_dict, "dict0")
close(f_dict)

sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s  = [(T1 / 1000, T2 / 1000) for (T1, T2) in sampled_pairs_ms]

sys = Scanner()

sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict(save_Mz=true)
sim_params["gpu"] = true

seq = read_seq("mpf_001_new_short1.seq")

x_pos = collect(range(-10e-3, 10e-3, length=101))

for (i, (T1, T2)) in enumerate(sampled_pairs_s)
    key = (Int(round(T1 * 1000)), Int(round(T2 * 1000)))
    filepath = joinpath(out_folder, "signal_$(key[1])_$(key[2]).jld2")

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


    sig_clean = dropdims(sig, dims=4)
    sig_channel1 = sig_clean[:, :, 1]
    signal_mag = vec(sum(sig_channel1, dims=2))

    @save filepath key signal_mag

    if i % 50 == 0 || i == length(sampled_pairs_s)
        println("Processed $i / $(length(sampled_pairs_s)) entries")
        flush(stdout)
    end
end

files = readdir(out_folder)
num_entries = length(files)
timepoints = 1000

bloch_matrix = zeros(ComplexF32, timepoints, num_entries)
idx_bloch = zeros(ComplexF32, num_entries, 2)

for (j, file) in enumerate(files)
    filepath = joinpath(out_folder, file)
    @load filepath key signal_mag
    bloch_matrix[:, j] .= signal_mag
    idx_bloch[j, :] .= key
end

matwrite(out_file, Dict("dict0" => bloch_matrix, "idx" => idx_bloch))
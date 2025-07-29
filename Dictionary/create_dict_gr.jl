using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra

out_folder = joinpath("progress", "progress_10mm_101_gr_short")
out_file = "blochdict_10mm_101_gr_short.mat"
mkpath(out_folder)

batch = Dict{Tuple{Int, Int}, Vector{ComplexF32}}()
batch_size = 50 

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
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = false

seq = read_seq("sequences/mpf_001_PhantomStudy_short_124.seq")

N = 101
x_min, x_max = -5e-3, 5e-3
ϕ = (√5 - 1) / 2  # inverse of golden ratio ≈0.618
x_unit = [mod(i * ϕ, 1) for i in 0:N-1] # quasi-random in [0, 1]
x_pos = x_min .+ x_unit .* (x_max - x_min)  # scale

println(x_unit)
println(x_pos)


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


    sig_clean = dropdims(sig; dims=(3,4))
    signal_mag = vec(sum(sig_clean, dims=2))

    batch[key] = signal_mag

    if length(batch) >= batch_size || i == length(sampled_pairs_s)
        for (key, sig) in batch
            filepath = joinpath(out_folder, "signal_$(key[1])_$(key[2]).jld2")
            @save filepath key signal_mag=sig
        end
        println("Saved batch at iteration $i (size: $(length(batch)))")
        flush(stdout)
        empty!(batch)
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
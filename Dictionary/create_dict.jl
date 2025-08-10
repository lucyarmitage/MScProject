using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra, CUDA

phantom_length = 7        # mm
num_points = 11           # num sample points
batch_size_pairs = 800      # phantoms per batch
sliceOrientation = 1      # 1=coronal, 2=transverse, 3=sagittal

seq = read_seq("sequences/mpf_001_PhantomStudy_short_124.seq")

run_name = "$(phantom_length)mm_$(num_points)_short"

phantom_length_m = phantom_length / 1000
pos = collect(range(-phantom_length_m/2, phantom_length_m/2, length=num_points))

zN = zeros(Float64, num_points)
xvec, yvec, zvec = if sliceOrientation == 1
    (pos, zN, zN)  
elseif sliceOrientation == 2
    (zN, zN, pos)       
elseif sliceOrientation == 3
    (zN, pos, zN)   
end

function build_combined_phantom(pairs_chunk::AbstractVector{<:Tuple{<:Real,<:Real}})
    N = length(pairs_chunk) * num_points
    x   = Vector{Float64}(undef, N)
    y   = Vector{Float64}(undef, N)
    z   = Vector{Float64}(undef, N)
    T1v = Vector{Float64}(undef, N)
    T2v = Vector{Float64}(undef, N)

    for (phantom_num, (T1, T2)) in enumerate(pairs_chunk)
        spins = ((phantom_num-1)*num_points + 1):(phantom_num*num_points)
        @views copyto!(x[spins], xvec)
        @views copyto!(y[spins], yvec)
        @views copyto!(z[spins], zvec)
        @views fill!(T1v[spins], Float64(T1))
        @views fill!(T2v[spins], Float64(T2))
    end

    return Phantom{Float64}(x=x, y=y, z=z, T1=T1v, T2=T2v)
end

out_folder = joinpath("/vol/bitbucket/la724/progress", "test_$run_name")
out_file = "dict/test_$run_name.mat"
mkpath(out_folder)

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx"); close(f_idx)

sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s = [(T1 / 1000, T2 / 1000) for (T1, T2) in sampled_pairs_ms]

sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true

batch_results = Dict{Tuple{Int, Int}, Vector{ComplexF32}}()

for chunk_start in 1:batch_size_pairs:length(sampled_pairs_s)
    chunk_end = min(chunk_start + batch_size_pairs - 1, length(sampled_pairs_s))
    pairs_chunk = sampled_pairs_s[chunk_start:chunk_end]
    keys_chunk = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_chunk]

    if all(isfile(joinpath(out_folder, "signal_$(k[1])_$(k[2]).jld2")) for k in keys_chunk)
        println("Skipping batch $chunk_start:$chunk_end")
        continue
    end

    big_phantom = build_combined_phantom(pairs_chunk)

    @suppress sig_all = simulate(big_phantom, seq, sys; sim_params=sim_params)
    sig_clean_all = dropdims(sig_all; dims=(3,4))

    start_idx = 1
    for key in keys_chunk
        sig_this = sig_clean_all[:, start_idx:start_idx + num_points - 1]
        batch_results[key] = vec(sum(sig_this, dims=2))
        start_idx += num_points
    end

    for (key, sig) in batch_results
        filepath = joinpath(out_folder, "signal_$(key[1])_$(key[2]).jld2")
        @save filepath key signal_mag=sig
    end
    println("Saved chunk $chunk_start:$chunk_end")
    empty!(batch_results)
end

files = readdir(out_folder)
num_entries = length(files)
timepoints = 1000

bloch_matrix = zeros(ComplexF32, timepoints, num_entries)
idx_bloch = zeros(Float32, num_entries, 2)

for (j, file) in enumerate(files)
    filepath = joinpath(out_folder, file)
    @load filepath key signal_mag
    bloch_matrix[:, j] .= signal_mag
    idx_bloch[j, :] .= key
end

matwrite(out_file, Dict("dict0" => bloch_matrix, "idx" => idx_bloch))

using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra, CUDA

# Constructs a combined phantom for a batch of (T1, T2) pairs.
function build_combined_phantom(pairs_chunk::AbstractVector{<:Tuple{<:Real,<:Real}}, num_points, xvec, yvec, zvec)
    N = length(pairs_chunk) * num_points
    x = Vector{Float64}(undef, N)
    y = Vector{Float64}(undef, N)
    z = Vector{Float64}(undef, N)
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

# Defaults
phantom_length = 8.0
num_points = 2001
sliceOrientation = 1      # 1=coronal, 2=transverse, 3=sagittal
seq_file = "sequences/mpf_001_PhantomStudy_short_124.seq"
timepoints = 1000

# Parse ARGS
# 1: phantom_length (mm, Float64)
# 2: num_points (Int)
# 3: sliceOrientation (Int: 1|2|3)
# 4: timepoints (Int)
# 5: seq_file (String)
# 6: out_file (String)

if length(ARGS) >= 1
    phantom_length = parse(Float64, ARGS[1])
end
if length(ARGS) >= 2
    num_points = parse(Int, ARGS[2])
end
if length(ARGS) >= 3
    sliceOrientation = parse(Int, ARGS[3])
end
if length(ARGS) >= 4
    timepoints = parse(Int, ARGS[4])
end
if length(ARGS) >= 5
    seq_file = ARGS[5]
end
if length(ARGS) >= 6
    out_file = ARGS[6]
else 
    out_file = "dict/dict_$(phantom_length)mm_$(num_points)_short.mat"
end

batch_size_pairs = Int(floor(240000 / num_points))    # Assuming 16 GB

# Single phantom
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
    error("Invalid sliceOrientation: $sliceOrientation")
end

# (T1,T2) pairs 
f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx"); close(f_idx)
sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s  = [(T1 / 1000, T2 / 1000) for (T1, T2) in sampled_pairs_ms]

# Simulation parameters 
seq = read_seq(seq_file)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true


# Batch simulation 
batch_results = Dict{Tuple{Int, Int}, Vector{ComplexF32}}()
for batch_start in 1:batch_size_pairs:length(sampled_pairs_s)
    # Select (T1,T2) pairs for batch 
    batch_end = min(batch_start + batch_size_pairs - 1, length(sampled_pairs_s))
    pairs_batch = sampled_pairs_s[batch_start:batch_end]
    keys_batch  = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_batch]

    # Build combined phantom
    big_phantom = build_combined_phantom(pairs_batch, num_points, xvec, yvec, zvec)

    # Simulate whole batch
    sig_all = simulate(big_phantom, seq, sys; sim_params=sim_params)
    sig_clean_all = dropdims(sig_all; dims=(3,4))

    # Sum over spins to integrate slice profile.
    start_idx = 1
    for key in keys_batch
        sig_this = sig_clean_all[:, start_idx:start_idx + num_points - 1]
        batch_results[key] = vec(sum(sig_this, dims=2))
        start_idx += num_points
    end

    println("Processed batch $batch_start:$batch_end")
end

# Collect (T1_ms,T2_ms) keys
files_keys = collect(keys(batch_results))
num_entries = length(files_keys)

bloch_matrix = zeros(ComplexF32, timepoints, num_entries)  # Nt × N complex dictionary
idx_bloch = zeros(Float32, num_entries, 2)    # N × 2 table of [T1_ms, T2_ms]

# Fill columns with signal for each (T1,T2)
for (j, key) in enumerate(files_keys)
    bloch_matrix[:, j] .= batch_results[key]
    idx_bloch[j, :] .= key
end

# Save as .mat
matwrite(out_file, Dict("dict0" => bloch_matrix, "idx" => idx_bloch))
println("Saved dictionary to $out_file")
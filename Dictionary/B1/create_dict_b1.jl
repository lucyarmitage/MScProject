using KomaMRI, MAT, Suppressor, JLD2, FileIO, LinearAlgebra, CUDA

function build_combined_phantom(pairs_chunk::AbstractVector{<:Tuple{<:Real,<:Real}}, num_points, xvec, yvec, zvec)
    N = length(pairs_chunk) * num_points
    x  = Vector{Float64}(undef, N)
    y  = Vector{Float64}(undef, N)
    z  = Vector{Float64}(undef, N)
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

phantom_length = 8
num_points = 5001
sliceOrientation = 1
seq_file = "sequences/mpf_001_PhantomStudy_short_124.seq"
timepoints = 1000
batch_size_pairs = 120

if length(ARGS) >= 1
    b1_scale = parse(Float64, ARGS[1])
end
if length(ARGS) >= 2
    phantom_length = parse(Float64, ARGS[2])
end
if length(ARGS) >= 3
    num_points = parse(Int, ARGS[3])
end
if length(ARGS) >= 4
    sliceOrientation = parse(Int, ARGS[4])
end
if length(ARGS) >= 5
    timepoints = parse(Int, ARGS[5])
end
if length(ARGS) >= 6
    seq_file = ARGS[6]
end

out_file = if length(ARGS) >= 7
    ARGS[7]
else
    "dict_b1/dict_$(b1_scale)_$(phantom_length)mm_$(num_points)_short.mat"
end

mkpath(dirname(out_file))

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

f_idx = matopen("D_IDX_SP_Phantom2025.mat")
idx = read(f_idx, "idx"); close(f_idx)
sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx)]
sampled_pairs_s  = [(T1 / 1000, T2 / 1000) for (T1, T2) in sampled_pairs_ms]

seq_nom = read_seq(seq_file)
seq     = (b1_scale + 0im) * seq_nom      # B1 SCALING
sys     = Scanner()

sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"]  = BlochDict()
sim_params["gpu"] = true 

batch_results = Dict{Tuple{Int, Int}, Vector{ComplexF32}}()
for batch_start in 1:batch_size_pairs:length(sampled_pairs_s)
    batch_end   = min(batch_start + batch_size_pairs - 1, length(sampled_pairs_s))
    pairs_batch = sampled_pairs_s[batch_start:batch_end]
    keys_batch  = [(Int(round(T1*1000)), Int(round(T2*1000))) for (T1, T2) in pairs_batch]

    big_phantom = build_combined_phantom(pairs_batch, num_points, xvec, yvec, zvec)

    sig_all = simulate(big_phantom, seq, sys; sim_params=sim_params)
    sig_clean_all = dropdims(sig_all; dims=(3,4))

    start_idx = 1
    for key in keys_batch
        sig_this = sig_clean_all[:, start_idx:start_idx + num_points - 1]
        batch_results[key] = vec(sum(sig_this, dims=2))
        start_idx += num_points
    end
end

files_keys = sort(collect(keys(batch_results)))
num_entries = length(files_keys)

Nt_sim = length(batch_results[files_keys[1]])
Nt_use = min(timepoints, Nt_sim)

bloch_matrix = zeros(ComplexF32, Nt_use, num_entries)
idx_bloch    = zeros(Float32, num_entries, 3)

for (j, key) in enumerate(files_keys)
    sig = batch_results[key]
    bloch_matrix[:, j] .= sig[1:Nt_use]
    idx_bloch[j, 1] = key[1]
    idx_bloch[j, 2] = key[2] 
    idx_bloch[j, 3] = Float32(b1_scale)
end

matwrite(out_file, Dict(
    "dict0"       => bloch_matrix,
    "idx"         => idx_bloch,         
))
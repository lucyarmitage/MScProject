using KomaMRI, MAT, Suppressor, FileIO, LinearAlgebra, CUDA, Printf, Dates

phantom_length = 10.0
num_points = 100

batch_sizes = vcat([1, 5, 10, 50, 100], 250:250:10000)
num_repeats = 8

seq_path = "sequences/mpf_001_PhantomStudy_short_124.seq"
idx_mat_path = "D_IDX_SP_Phantom2025.mat"
out_dir = "gpu_batch_bench"; isdir(out_dir) || mkpath(out_dir)
csv_path = joinpath(out_dir, "batch_bench_timings.csv")

function gpu_seconds(f::Function)
    CUDA.synchronize()
    result = nothing
    t = @elapsed begin
        result = f()
        CUDA.synchronize()
    end
    return (t, result)
end

function build_combined_phantom(pairs_chunk::AbstractVector{<:Tuple{<:Real,<:Real}},
                                num_points::Int,
                                xvec::AbstractVector{<:Real},
                                yvec::AbstractVector{<:Real},
                                zvec::AbstractVector{<:Real})
    N = length(pairs_chunk) * num_points
    x = Vector{Float64}(undef, N)
    y = Vector{Float64}(undef, N)
    z = Vector{Float64}(undef, N)
    T1v = Vector{Float64}(undef, N)
    T2v = Vector{Float64}(undef, N)
    @inbounds for (p, (T1, T2)) in enumerate(pairs_chunk)
        rng = ((p-1)*num_points + 1):(p*num_points)
        @views copyto!(x[rng], xvec); @views copyto!(y[rng], yvec); @views copyto!(z[rng], zvec)
        @views fill!(T1v[rng], Float64(T1)); @views fill!(T2v[rng], Float64(T2))
    end
    return Phantom{Float64}(x=x, y=y, z=z, T1=T1v, T2=T2v)
end

pct(part, total) = 100 * part / max(total, eps())

function chunk_pairs_cyclic(pairs::AbstractVector{<:Tuple{<:Real,<:Real}}, bs::Int)
    n = length(pairs)
    if bs <= n
        return pairs[1:bs]
    else
        q, r = divrem(bs, n)
        return vcat([pairs for _ in 1:q]..., pairs[1:r])
    end
end

const CSV_HEADER = join((
    "batch_size","repeat","timepoints","spins","samples",
    "t_gpu_s","t_build_s","t_split_s","t_total_s",
    "pct_gpu","pct_overhead","gpu_name","started_at","finished_at","checksum",
    "mem_total_B","mem_free_before_B","mem_used_before_B",
    "mem_free_after_B","mem_used_after_B","mem_used_after_pct"
), ",")

function ensure_csv_header(path::AbstractString)
    if !isfile(path)
        open(path, "w") do io
            println(io, CSV_HEADER)
        end
    end
end

function append_csv_row(path::AbstractString, row::NamedTuple)
    open(path, "a") do io
        println(io, join((
            row.batch_size, row.repeat, row.timepoints, row.spins, row.samples,
            @sprintf("%.6f", row.t_gpu_s), @sprintf("%.6f", row.t_build_s),
            @sprintf("%.6f", row.t_split_s), @sprintf("%.6f", row.t_total_s),
            @sprintf("%.2f", row.pct_gpu), @sprintf("%.2f", row.pct_overhead),
            row.gpu_name, row.started_at, row.finished_at,
            @sprintf("%.6e", row.checksum),
            row.mem_total_B, row.mem_free_before_B, row.mem_used_before_B,
            row.mem_free_after_B, row.mem_used_after_B, @sprintf("%.2f", row.mem_used_after_pct)
        ), ","))
    end
end

mat = matopen(idx_mat_path)
idx_tbl = read(mat, "idx"); close(mat)
sampled_pairs_ms = [(row[1], row[2]) for row in eachrow(idx_tbl)]
sampled_pairs_s = [(Float64(T1)/1000, Float64(T2)/1000) for (T1, T2) in sampled_pairs_ms]

phantom_length_m = phantom_length / 1000
pos = collect(range(-phantom_length_m/2, phantom_length_m/2, length=num_points))
zN = zeros(Float64, num_points)
xvec, yvec, zvec = (pos, zN, zN)

seq = read_seq(seq_path)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true

warm_bs = min(max(1, first(batch_sizes)), length(sampled_pairs_s))
warm_pairs = sampled_pairs_s[1:warm_bs]
warm_phantom = build_combined_phantom(warm_pairs, num_points, xvec, yvec, zvec)
@suppress _ = simulate(warm_phantom, seq, sys; sim_params=sim_params)
CUDA.synchronize()

ensure_csv_header(csv_path)
gpu_name = CUDA.name(CUDA.device())

for bs in batch_sizes
    pairs_chunk = chunk_pairs_cyclic(sampled_pairs_s, bs)
    for rep in 1:num_repeats
        t_total0 = time()
        started_at = Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM:SS")
        t0 = time()
        big_phantom = build_combined_phantom(pairs_chunk, num_points, xvec, yvec, zvec)
        t_build = time() - t0
        mem_total_B = CUDA.total_memory()
        mem_free_before_B = CUDA.free_memory()
        mem_used_before_B = mem_total_B - mem_free_before_B
        t_gpu, sig_all = gpu_seconds() do
            @suppress simulate(big_phantom, seq, sys; sim_params=sim_params)
        end
        sig = dropdims(sig_all; dims=(3,4))
        tp = size(sig, 1)
        total_spins = bs * num_points
        mem_free_after_B = CUDA.free_memory()
        mem_used_after_B = mem_total_B - mem_free_after_B
        mem_used_after_pct = 100 * float(mem_used_after_B) / max(float(mem_total_B), eps())
        t1 = time()
        start_idx = 1
        cs = 0.0
        @inbounds for j in 1:bs
            sig_this = @view sig[:, start_idx:start_idx + num_points - 1]
            fp = vec(sum(sig_this, dims=2))
            cs += sum(abs2, fp)
            start_idx += num_points
        end
        t_split = time() - t1
        t_total = time() - t_total0
        samples = bs * tp
        overhead = t_total - t_gpu
        finished_at = Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM:SS")
        row_nt = (; batch_size=bs, repeat=rep, timepoints=tp, spins=total_spins, samples=samples,
                  t_gpu_s=t_gpu, t_build_s=t_build, t_split_s=t_split, t_total_s=t_total,
                  pct_gpu=pct(t_gpu, t_total), pct_overhead=pct(overhead, t_total),
                  gpu_name, started_at, finished_at, checksum=cs,
                  mem_total_B=mem_total_B, mem_free_before_B=mem_free_before_B, mem_used_before_B=mem_used_before_B,
                  mem_free_after_B=mem_free_after_B, mem_used_after_B=mem_used_after_B, mem_used_after_pct=mem_used_after_pct)
        append_csv_row(csv_path, row_nt)
        @printf("APPEND B=%-4d rep=%d tp=%d total=%.3fs gpu=%.3fs (%.1f%%) build=%.3fs split=%.3fs checksum=%.3e\n",
                bs, rep, tp, t_total, t_gpu, pct(t_gpu, t_total), t_build, t_split, cs)
        sig_all = nothing; sig = nothing; big_phantom = nothing
        GC.gc(); CUDA.reclaim()
    end
end

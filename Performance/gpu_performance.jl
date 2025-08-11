using KomaMRI, Suppressor, CUDA, CSV, DataFrames, Statistics

const GPU_MONITOR_STOP = Ref(false)

function start_gpu_monitor(log_path; interval::Real=0.1, device::Int=0)
    GPU_MONITOR_STOP[] = false
    if !isfile(log_path) || filesize(log_path) == 0
        open(log_path, "w") do io
            println(io, "ts,util_gpu,util_mem,mem_used_MB")
        end
    end
    @async begin
        while !GPU_MONITOR_STOP[]
            out = read(`nvidia-smi -i $device --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used --format=csv,noheader,nounits`, String)
            line = isempty(out) ? "" : first(split(chomp(out), '\n'))
            if !isempty(line)
                open(log_path, "a") do io
                    println(io, replace(line, ", " => ","))
                end
            end
            sleep(interval)
        end
    end
    nothing
end


stop_gpu_monitor() = (GPU_MONITOR_STOP[] = true; sleep(0.2))


function summarise_gpu_log(log_path)
    df = CSV.read(log_path, DataFrame)
    (maximum(df.mem_used_MB), mean(df.util_gpu), nrow(df))
end


function gpu_seconds(f::Function)
    if isdefined(CUDA, :Event) && isdefined(CUDA, :elapsed_time)
        CUDA.synchronize()
        ev_start = CUDA.Event()
        ev_stop = CUDA.Event()
        CUDA.record(ev_start)
        result = f()
        CUDA.record(ev_stop)
        CUDA.synchronize(ev_stop)
        t_ms = CUDA.elapsed_time(ev_start, ev_stop)
        return (t_ms/1000.0, result)
    else
        CUDA.synchronize()
        t0 = time_ns()
        result = f()
        CUDA.synchronize()
        t_s = (time_ns() - t0) / 1e9
        return (t_s, result)
    end
end


function build_combined_phantom(pairs_chunk::AbstractVector{<:Tuple{<:Real,<:Real}})
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


batch_sizes = [1, 5, 10, 50, 100, 200, 400, 600, 800, 1000, 1250, 1500]
num_repeats = 8
monitor_dt = 0.1        
device_id = 0
log_dir = "gpu_bench_logs"
mkpath(log_dir)

phantom_length_m = 10e-3
num_points = 101
pos = collect(range(-phantom_length_m/2, phantom_length_m/2, length=num_points))
zN = zeros(Float64, num_points)
xvec, yvec, zvec = (pos, zN, zN)

# Synthetic pairs 
sampled_pairs_s = [(1.0, 0.1) for _ in 1:maximum(batch_sizes)]

seq = read_seq("sequences/mpf_001_PhantomStudy_short_124.seq")
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "mat"
sim_params["sim_method"] = BlochDict()
sim_params["gpu"] = true



## WARM UP
warmup_pairs   = sampled_pairs_s[1:10]
warmup_phantom = build_combined_phantom(warmup_pairs)
@suppress _ = simulate(warmup_phantom, seq, sys; sim_params=sim_params)
CUDA.synchronize()

# Data frames to fill
runs_df = DataFrame(
    batch_size = Int[],
    repeat = Int[],
    spins  = Int[],
    timepoints = Int[],
    samples = Int[],
    time_s = Float64[],
    spins_per_s = Float64[],
    samples_per_s = Float64[],
    peak_mem_MB  = Float64[],
    avg_util_gpu = Float64[],
    util_samples = Int[],
)

for bs in batch_sizes
    pairs_chunk = sampled_pairs_s[1:bs]

    for rep in 1:num_repeats
        ph = build_combined_phantom(pairs_chunk)

        logfile = joinpath(log_dir, "gpu_bs_$(bs)_rep$(rep).csv")
        start_gpu_monitor(logfile; interval=monitor_dt, device=device_id)

        GC.gc(); CUDA.reclaim(); CUDA.synchronize()

        t_s, sig_all = gpu_seconds() do
            @suppress simulate(ph, seq, sys; sim_params=sim_params)
        end

        stop_gpu_monitor()

        timepoints = 1000
        spins = bs * num_points
        samples = timepoints * spins

        peak_mem, avg_util, nsamples = summarise_gpu_log(logfile)

        push!(runs_df, (
            bs, rep, spins, timepoints, samples, t_s,
            spins / t_s, samples / t_s, peak_mem, avg_util, nsamples
        ))
    end
end

runs_csv = joinpath(log_dir, "gpu_bench_runs.csv")
CSV.write(runs_csv, runs_df)

summary_df = combine(groupby(runs_df, :batch_size)) do g
    (;
        n_repeats = nrow(g),
        spins = first(g.spins),
        timepoints = first(g.timepoints),
        mean_time_s = mean(g.time_s),
        median_time_s = median(g.time_s),
        std_time_s = std(g.time_s),
        best_time_s = minimum(g.time_s),
        worst_time_s = maximum(g.time_s),
        mean_spins_per_s = mean(g.spins_per_s),
        median_spins_per_s = median(g.spins_per_s),
        mean_samples_per_s = mean(g.samples_per_s),
        median_samples_per_s = median(g.samples_per_s),
        max_peak_mem_MB = maximum(g.peak_mem_MB),
        mean_avg_util_gpu = mean(g.avg_util_gpu),
    )
end

summary_csv = joinpath(log_dir, "gpu_bench_summary.csv")
CSV.write(summary_csv, summary_df)

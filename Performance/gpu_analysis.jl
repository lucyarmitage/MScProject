using CSV, DataFrames, Statistics, StatsPlots, CategoricalArrays, Printf

runs = CSV.read("gpu_bench_runs.csv", DataFrame)

g = groupby(runs, :batch_size)

summary = combine(g,
    :time_s => (x -> length(x)) => :n_repeats,
    :time_s => median => :median_time_s,
    :time_s => mean => :mean_time_s,
    :time_s => std => :sd_time_s,
    :samples_per_s => median => :median_samples_per_s,
    :samples_per_s => mean => :mean_samples_per_s,
    :samples_per_s => std => :sd_samples_per_s,
    :avg_util_gpu => mean => :mean_avg_util_gpu,
    :avg_util_gpu => std => :sd_avg_util_gpu,
    :peak_mem_MB => mean => :mean_peak_mem_MB,
    :peak_mem_MB => std => :sd_peak_mem_MB
)

mkpath("gpu_plots")
CSV.write("gpu_plots/summary_min.csv", summary)

# Throughput vs batch size (mean ± SD)
p1 = scatter(summary.batch_size, summary.mean_samples_per_s;
    yerror = summary.sd_samples_per_s,
    xlabel = "Batch size",
    ylabel = "Mean samples/sec",
    title  = "Throughput vs Batch size",
    legend = false)
savefig(p1, "gpu_plots/throughput_vs_batch.png")

# Mean GPU utilisation vs batch size (mean ± SD)
p2 = scatter(summary.batch_size, summary.mean_avg_util_gpu;
    yerror = summary.sd_avg_util_gpu,
    xlabel = "Batch size", ylabel = "Mean GPU util (%)",
    title  = "GPU utilisation vs Batch size",
    legend = false)
savefig(p2, "gpu_plots/gpu_util_vs_batch.png")

# Peak GPU memory vs batch size (mean ± SD)
p3 = scatter(summary.batch_size, summary.mean_peak_mem_MB;
    yerror = summary.sd_peak_mem_MB,
    xlabel = "Batch size", ylabel = "Peak memory (MB)",
    title  = "Peak GPU memory vs Batch size",
    legend = false)
savefig(p3, "gpu_plots/peak_mem_vs_batch.png")
plot(p3)

# Per-phantom time vs batch size (mean ± 95% CI)
pp = combine(groupby(runs, :batch_size)) do g
    x  = g.time_s ./ g.batch_size
    m  = mean(x)
    sd = coalesce(std(x), 0.0)
    ci = 1.96 * sd / sqrt(max(1, nrow(g)))
    (; mean_pp_time_s = m, ci95_pp = ci)
end
p5 = scatter(pp.batch_size, pp.mean_pp_time_s;
    yerror = pp.ci95_pp,
    xlabel = "Batch size", ylabel = "time per phantom (s)",
    title  = "Per-phantom time vs batch size",
    legend = false)
savefig(p5, "gpu_plots/per_phantom_time.png")
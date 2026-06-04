## Per-walker diagnostic for an Emcee chain.
##
## Goal: figure out why a non-trivial fraction of samples sit at very low logp
## (the unmixedness artefact that derived_time_distribution had to filter).
## Three possibilities to distinguish:
##   (a) burn-in — all walkers start bad, settle later; cut early iterations.
##   (b) some walkers permanently stuck at low logp throughout; drop those.
##   (c) walkers cross prior boundaries and stay there.
##
## Usage:
##     julia --project=. diagnose_walkers.jl                              # combined chain, default
##     julia --project=. diagnose_walkers.jl results/emcee-chain-kr81.jld2

cd(@__DIR__)
import Pkg; Pkg.activate(".")

using JLD2
using DataFrames
using CairoMakie
using Dates
using Statistics
using Turing

include("BasalMixingModel.jl")

path = length(ARGS) >= 1 ? ARGS[1] : "results/emcee-chain-combined.jld2"
results = load_ensemble_results(path)
tag = splitext(basename(path))[1]   # e.g. "emcee-chain-combined"

ch = results.chain
n_iter, n_walkers = size(ch)
df = DataFrame(ch)
df.logjoint = vec(ch[:logjoint])
println("Chain: $path")
println("  iters=$n_iter, walkers=$n_walkers, total=$(nrow(df))")
println("  iter col: $(extrema(df.iter))   walker col: $(extrema(df.chain))")

# Reshape logjoint and t_0 into (n_iter × n_walkers) matrices.
# Sort by (chain, iter) so reshape gives one walker per column.
sort!(df, [:chain, :iter])
logp_mat = reshape(df.logjoint, n_iter, n_walkers)
t0_mat   = reshape(df.t_0,      n_iter, n_walkers)

mean_logp_per_walker = vec(mean(logp_mat; dims=1))
median_logp_per_walker = vec(median(logp_mat; dims=1))
max_logp = maximum(df.logjoint)

println("\nPer-walker mean logp (sorted):")
order = sortperm(mean_logp_per_walker)
for w in order
    println("  walker $(lpad(w, 2)): mean=$(round(mean_logp_per_walker[w], digits=2))  " *
            "median=$(round(median_logp_per_walker[w], digits=2))  " *
            "min=$(round(minimum(logp_mat[:, w]), digits=2))  " *
            "max=$(round(maximum(logp_mat[:, w]), digits=2))")
end

# Decide which walkers are "stuck": mean logp more than 5 below the best walker's mean.
best_walker_mean = maximum(mean_logp_per_walker)
stuck_walkers = findall(mean_logp_per_walker .< best_walker_mean - 5)
healthy_walkers = setdiff(1:n_walkers, stuck_walkers)
println("\nBest walker mean logp: $(round(best_walker_mean, digits=2))")
println("Stuck walkers (mean > 5 below best): $stuck_walkers  ($(length(stuck_walkers))/$n_walkers)")
println("Healthy walkers: $(length(healthy_walkers))/$n_walkers")

# Per-walker fraction-of-samples-near-prior-bound (-700 ± 5).
t0_high = vec(mean(abs.(t0_mat .- (-700.0)) .< 5.0; dims=1))
println("\nPer-walker fraction of samples at t_0 ∈ [-705, -695]:")
for w in order
    f = t0_high[w]
    if f > 0.05
        println("  walker $(lpad(w, 2)): $(round(f*100, digits=1)) % at the wall")
    end
end

## Trace figure ##
fig = Figure(size=(1100, 850))

# Panel A: logp vs iteration, one line per walker
axA = Axis(fig[1, 1], xlabel="iteration", ylabel="log p",
           title="logp trace per walker  ($tag)")
ylims!(axA, (-60, 5))
for w in 1:n_walkers
    col = w in stuck_walkers ? (:black, 0.7) : (:steelblue, 0.4)
    lw  = w in stuck_walkers ? 1.5 : 0.8
    lines!(axA, 1:n_iter, logp_mat[:, w]; color=col, linewidth=lw)
end
hlines!(axA, [max_logp]; color=:red, linestyle=:dash, linewidth=1)
text!(axA, "MAP = $(round(max_logp, digits=2))",
      position=(n_iter*0.7, max_logp+1), color=:red, fontsize=10)

# Panel B: t_0 vs iteration, one line per walker
axB = Axis(fig[1, 2], xlabel="iteration", ylabel="t_0 (kyr)",
           title="t_0 trace per walker")
ylims!(axB, (-3100, -600))
for w in 1:n_walkers
    col = w in stuck_walkers ? (:black, 0.7) : (:steelblue, 0.4)
    lw  = w in stuck_walkers ? 1.5 : 0.8
    lines!(axB, 1:n_iter, t0_mat[:, w]; color=col, linewidth=lw)
end
hlines!(axB, [-700.0, -3000.0]; color=:grey50, linestyle=:dot, linewidth=1)

# Panel C: histogram of mean-logp-per-walker
axC = Axis(fig[2, 1], xlabel="walker mean logp", ylabel="count of walkers",
           title="distribution of per-walker mean logp")
hist!(axC, mean_logp_per_walker; bins=20, color=(:steelblue, 0.7))
vlines!(axC, [best_walker_mean - 5]; color=:red, linestyle=:dash,
        label="stuck cutoff")
axislegend(axC; position=:rt)

# Panel D: t_0 vs logp scatter, coloured by walker class
axD = Axis(fig[2, 2], xlabel="t_0 (kyr)", ylabel="log p",
           title="t_0 vs logp, coloured by walker class")
ylims!(axD, (-60, 5))
for w in stuck_walkers
    scatter!(axD, t0_mat[:, w], logp_mat[:, w];
             color=(:black, 0.4), markersize=3)
end
for w in healthy_walkers
    scatter!(axD, t0_mat[:, w], logp_mat[:, w];
             color=(:steelblue, 0.4), markersize=3)
end

mkpath("plots")
outpath = joinpath("plots", string(Dates.today())*"_walker-diagnostic-$tag.png")
mysave(outpath, fig)
println("\nSaved diagnostic plot: $outpath")

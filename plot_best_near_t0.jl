## Plot a "best-like" run from the :combined chain restricted to a t_0 band.
##
## Usage:
##     julia --project=. plot_best_near_t0.jl                  # default band -750..-700
##     julia --project=. plot_best_near_t0.jl -750 -650        # custom band
##
## Picks argmax(logjoint) within the band, re-runs the forward model at those
## MAP-in-band parameters, and saves a standard four-panel figure.

cd(@__DIR__)
import Pkg; Pkg.activate(".")

using JLD2
using DataFrames
using CairoMakie
using Dates
using Distributions
using Turing

include("BasalMixingModel.jl")

t_lo = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : -750.0
t_hi = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : -700.0

results = load_ensemble_results("results/emcee-chain-combined.jld2")
df = DataFrame(results.chain)
df.logjoint = vec(results.chain[:logjoint])

mask = (df.t_0 .>= t_lo) .& (df.t_0 .<= t_hi)
n_in = count(mask)
println("Samples with t_0 in [$t_lo, $t_hi]: $n_in / $(nrow(df))")
n_in > 0 || error("No samples in the band — widen it.")

dfb = df[mask, :]
best_idx = argmax(dfb.logjoint)

p_best = (
    delta   = dfb.delta[best_idx],
    m_clean = dfb.m_clean[best_idx],
    f_dirty = dfb.f_dirty[best_idx],
    t_old   = dfb.t_old[best_idx],
    F_ar40  = dfb.F_ar40[best_idx],
)
t_0_pick = dfb.t_0[best_idx]
logp_pick = dfb.logjoint[best_idx]

println("Picked sample:")
println("  t_0     = $(round(t_0_pick; digits=2)) kyr")
println("  logp    = $(round(logp_pick; digits=2))")
println("  delta   = $(round(p_best.delta; digits=4))")
println("  m_clean = $(round(p_best.m_clean; digits=5))")
println("  f_dirty = $(round(p_best.f_dirty; digits=4))")
println("  t_old   = $(round(p_best.t_old; digits=2))")
println("  F_ar40  = $(round(p_best.F_ar40; digits=5))")

b = BasalMixingModel(depth=results.depth,
                     k81_obs_depths=results.k81.depth,
                     dar40_obs_depths=results.dar40.depth)
t_elapsed = -t_0_pick
RunBasalMixingModel!(p_best, b, (results.k81, results.dar40); dt=0.1, t1=t_elapsed)
b.k81.time_min = t_elapsed
b.dar40.time_min = t_elapsed
b.joint.time_min = t_elapsed

fig = plot_BasalMixingModelRun(b; k81_obs=results.k81, dar40_obs=results.dar40,
                               t_max=t_elapsed)

setup_suffix = isempty(results.setup) ? "" : "_$(results.setup)"
band_tag = "t0_$(round(Int, t_lo))to$(round(Int, t_hi))"
outpath = joinpath("plots",
                   string(Dates.today())*"_mixingmodel-best-in-band-$band_tag$setup_suffix.png")
mysave(outpath, fig)
println("Saved: $outpath")

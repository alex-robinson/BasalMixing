## Plotting for the basal-mixing ensemble run.
##
## Usage (interactive, after running run_basalmixing_ensemble.jl in the same
## session — `chain`, `k81`, `dar40`, `depth`, `priors` already live in `Main`):
##
##     include("plot_basalmixing_ensemble.jl")
##     plot_ensemble("results/emcee-chain-combined.jld2")
##
## Usage (overlay kr81-only on combined, the standard comparison figure):
##
##     plot_ensemble_comparison("results/emcee-chain-combined.jld2",
##                              "results/emcee-chain-kr81.jld2")
##
## Usage (one-shot from the shell, single chain):
##
##     julia --project=. plot_basalmixing_ensemble.jl results/emcee-chain-combined.jld2

## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using JLD2
using DataFrames
using CairoMakie
using Dates
using Distributions
using Turing   # loads FlexiChains so JLD2 can reconstruct the chain type on read

include("BasalMixingModel.jl")

function plot_prior_line!(ax, prior::Distribution; color=:red, kwargs...)
    lo, hi = quantile(prior, 0.001), quantile(prior, 0.999)
    x = range(lo, hi, length=300)
    lines!(ax, x, pdf.(prior, x); color=color, kwargs...)
end

"""
    logp_mask(df, logp_window) -> BitVector

Pragmatic filter for unmixed-chain visualisations: keep only samples with
`logjoint >= max(logjoint) - logp_window`. 10 logp below MAP is a posterior
weight of ~exp(-10) ≈ 4.5e-5 — negligible. Anything below that contributes
nothing to the posterior in a well-mixed chain, so if many samples accumulate
there they're an unmixedness artefact (e.g. Emcee walkers stuck at a prior
boundary) and should be excluded from the histogram / pdf estimates.

The raw logp scatter plot is left UNfiltered — that's where you want to *see*
the unmixed structure as a diagnostic. Set `logp_window=Inf` to disable.
"""
function logp_mask(df, logp_window::Float64)
    isfinite(logp_window) || return trues(nrow(df))
    return df.logjoint .>= maximum(df.logjoint) - logp_window
end

"""
    map_best_run(chain, k81, dar40, depth) -> (b, t_elapsed_map, p_best, df, best_idx)

Re-run the forward model at the chain's MAP parameters. The sampled-time model
only ever evaluated the likelihood at the single endpoint t = |t_0|, so the
re-run integrates exactly 0 → |t_0| and the trajectory beyond that is left
unrendered by `plot_BasalMixingModelRun` (it crops at `t_max`).
"""
function map_best_run(chain, k81, dar40, depth)
    df = DataFrame(chain)
    df.logjoint = vec(chain[:logjoint])
    best_idx = argmax(df.logjoint)

    p_best = (
        delta   = df.delta[best_idx],
        m_clean = df.m_clean[best_idx],
        f_dirty = df.f_dirty[best_idx],
        t_old   = df.t_old[best_idx],
        F_ar40  = df.F_ar40[best_idx],
    )
    b = BasalMixingModel(depth=depth, k81_obs_depths=k81.depth, dar40_obs_depths=dar40.depth)
    t_elapsed_map = -df.t_0[best_idx]
    RunBasalMixingModel!(p_best, b, (k81, dar40); dt=0.1, t1=t_elapsed_map)
    # Override argmin-RMSE: the inference time is the endpoint, by construction.
    b.k81.time_min   = t_elapsed_map
    b.dar40.time_min = t_elapsed_map
    b.joint.time_min = t_elapsed_map
    return (; b, t_elapsed_map, p_best, df, best_idx)
end

"""
    plot_ensemble(; chain, k81, dar40, depth, priors, ...; secondary=nothing)
    plot_ensemble(results::NamedTuple; secondary=nothing, kwargs...)
    plot_ensemble(path::String; kwargs...)

Render the three standard ensemble figures and return them as a `NamedTuple
(best, logp, hist)`:

- *best*: depth profiles of mixing rate, ⁸¹Kr closed-system age, and δ⁴⁰Ar
  at the MAP parameters (re-runs the forward model with `dt=0.1`).
- *logp*: per-parameter scatter of log-joint vs. parameter value, MAP highlighted.
- *hist*: marginal posterior histograms overlaid with each parameter's prior.

If `secondary` is provided (a NamedTuple from `load_ensemble_results`), each
figure additionally overlays the secondary chain's MAP run / scatter /
histograms in a dashed black style. Use this to compare e.g. the combined
likelihood (primary) against kr81-only (secondary).
"""
function plot_ensemble(;
    chain,
    k81,
    dar40,
    depth,
    priors,
    setup::AbstractString="",
    sampler_choice::Symbol=:unknown,
    likelihood::Symbol=:combined,
    secondary=nothing,
    logp_window::Float64=10.0,
    outdir::AbstractString="plots",
    save_figures::Bool=true,
)
    params = [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0]
    labels = ["delta (m)", "m_clean (m/yr)", "f_dirty",
              "t_old (kyr)", "F_ar40 (cc/m²/kyr)", "t_0 (kyr)"]

    primary = map_best_run(chain, k81, dar40, depth)
    sec     = secondary === nothing ? nothing :
              map_best_run(secondary.chain, secondary.k81, secondary.dar40, secondary.depth)

    primary_label = "primary (:$likelihood)"
    secondary_label = secondary === nothing ? nothing :
                      "secondary (:$(get(secondary, :likelihood, :unknown)))"

    ## Best-fit depth profiles ##
    overlay = sec === nothing ? nothing : (; b=sec.b, t_max=sec.t_elapsed_map, label=secondary_label)
    fig_best = plot_BasalMixingModelRun(primary.b;
                                        k81_obs=k81, dar40_obs=dar40,
                                        t_max=primary.t_elapsed_map,
                                        overlay=overlay)

    ## log-posterior scatter ##
    fig_logp = Figure(size=(900, 700))
    ipar = 0
    for (param, label) in zip(params, labels)
        string(param) in names(primary.df) || continue
        ipar += 1
        row, col = divrem(ipar-1, 2)
        ax = Axis(fig_logp[row+1, col+1], xlabel=label, ylabel="log p")
        ylims!(ax, (-50, 0))
        if sec !== nothing && string(param) in names(sec.df)
            scatter!(ax, sec.df[!, param], sec.df.logjoint;
                     alpha=0.4, markersize=5, color=:black,
                     label=secondary_label)
            scatter!(ax, [sec.df[sec.best_idx, param]], [sec.df.logjoint[sec.best_idx]];
                     color=:black, marker=:diamond, markersize=10)
        end
        scatter!(ax, primary.df[!, param], primary.df.logjoint;
                 alpha=0.6, markersize=6, color=:steelblue,
                 label=primary_label)
        scatter!(ax, [primary.df[primary.best_idx, param]], [primary.df.logjoint[primary.best_idx]];
                 color=:red, markersize=12, label="MAP")
    end

    ## Marginal histograms vs priors ##
    # Filter to samples within `logp_window` of MAP — see logp_mask docstring.
    mask_p = logp_mask(primary.df, logp_window)
    mask_s = sec === nothing ? nothing : logp_mask(sec.df, logp_window)
    n_kept_p = count(mask_p)
    n_kept_s = mask_s === nothing ? 0 : count(mask_s)
    println("Histogram filter: primary kept $n_kept_p/$(nrow(primary.df)) " *
            "(logp_window=$logp_window)" *
            (sec === nothing ? "" : "; secondary kept $n_kept_s/$(nrow(sec.df))"))

    fig_hist = Figure(size=(900, 700))
    ipar = 0
    for (param, label) in zip(params, labels)
        string(param) in names(primary.df) || continue
        ipar += 1
        row, col = divrem(ipar-1, 2)
        ax = Axis(fig_hist[row+1, col+1], xlabel=label, ylabel="density")
        hist!(ax, primary.df[mask_p, param]; bins=20, normalization=:pdf,
              color=(:steelblue, 0.7), label=primary_label)
        if sec !== nothing && string(param) in names(sec.df)
            # Use stephist (outline only) so the secondary doesn't visually
            # bury the primary fill.
            stephist!(ax, sec.df[mask_s, param]; bins=20, normalization=:pdf,
                      color=:black, linewidth=2, linestyle=:dash,
                      label=secondary_label)
        end
        if haskey(priors, param) && priors[param] isa Distribution
            plot_prior_line!(ax, priors[param];
                             label="Prior", linewidth=2, color=:grey50)
        end
        vlines!(ax, [primary.df[primary.best_idx, param]]; color=:red, linewidth=2, label="MAP")
    end

    if save_figures
        mkpath(outdir)
        prefix = joinpath(outdir, string(Dates.today())*"_")
        suffix = isempty(setup) ? "" : "_$setup"
        tag    = secondary === nothing ? "$(likelihood)" :
                 "$(likelihood)-vs-$(get(secondary, :likelihood, :secondary))"
        mysave(prefix*"mixingmodel-ens-best-$tag$suffix.png", fig_best)
        mysave(prefix*"mixingmodel-ens-logp-$tag$suffix.png", fig_logp)
        mysave(prefix*"mixingmodel-ens-hist-$tag$suffix.png", fig_hist)
    end

    return (best=fig_best, logp=fig_logp, hist=fig_hist)
end

plot_ensemble(results::NamedTuple; secondary=nothing, kwargs...) =
    plot_ensemble(;
        chain=results.chain,
        k81=results.k81,
        dar40=results.dar40,
        depth=results.depth,
        priors=results.priors,
        setup=get(results, :setup, ""),
        sampler_choice=get(results, :sampler_choice, :unknown),
        likelihood=get(results, :likelihood, :combined),
        secondary=secondary,
        kwargs...,
    )

plot_ensemble(path::AbstractString; kwargs...) =
    plot_ensemble(load_ensemble_results(String(path)); kwargs...)

"""
    plot_ensemble_comparison(primary_path, secondary_path; kwargs...)

Convenience wrapper: load two chains (e.g. combined vs kr81-only) and produce
the standard three figures with the secondary overlaid as a dashed black line
on top of the primary. Returns the same `(best, logp, hist)` NamedTuple as
`plot_ensemble`.
"""
function plot_ensemble_comparison(primary_path::AbstractString,
                                  secondary_path::AbstractString;
                                  kwargs...)
    primary   = load_ensemble_results(String(primary_path))
    secondary = load_ensemble_results(String(secondary_path))
    return plot_ensemble(primary; secondary=secondary, kwargs...)
end

"""
    derived_time_distribution(results::NamedTuple;
                              t_min=-3000.0, t_max=0.0,
                              logp_window=10.0) -> (t_grid, p_t, map_t, n_kept, n_total)

Empirical posterior on `t_0` (signed kyr BP) from the chain. The sampled-time
model puts `t_0` directly in the chain, so this is just a histogram on the
1-kyr grid spanning the prior support.

`logp_window` filters out samples with `logjoint < max(logjoint) - logp_window`
(see `logp_mask`) — needed when the Emcee chain has walkers stuck at a prior
boundary that contribute zero posterior weight but inflate the raw histogram.

`map_t` is the `t_0` of the `argmax(logjoint)` sample (NOT the histogram mode),
which is the right summary when the histogram is concentrated by unmixedness
rather than by likelihood.
"""
function derived_time_distribution(results::NamedTuple;
                                   t_min::Float64=-3000.0,
                                   t_max::Float64=0.0,
                                   logp_window::Float64=10.0)
    df = DataFrame(results.chain)
    df.logjoint = vec(results.chain[:logjoint])
    mask = logp_mask(df, logp_window)
    n_kept = count(mask)
    t_grid = collect(t_min:1.0:t_max)
    counts = zeros(Float64, length(t_grid))
    for v in df.t_0[mask]
        k = clamp(round(Int, v - t_grid[1]) + 1, 1, length(t_grid))
        counts[k] += 1.0
    end
    total = sum(counts)
    p_t = total > 0 ? counts ./ total : counts
    # MAP from argmax(logjoint), not histogram mode.
    map_t = df.t_0[argmax(df.logjoint)]
    return (t_grid=t_grid, p_t=p_t, map_t=map_t, n_kept=n_kept, n_total=nrow(df))
end

"""
    plot_derived_time(results; secondary=nothing, ...) -> (fig, t_grid, p_t, map_t)

Plot the sampled `t_0` posterior. If `secondary` is provided, overlay the
secondary chain's pdf as a dashed black line.
"""
function plot_derived_time(results::NamedTuple;
                           secondary=nothing,
                           logp_window::Float64=10.0,
                           outdir::AbstractString="plots",
                           save_figures::Bool=true)
    d     = derived_time_distribution(results; logp_window=logp_window)
    setup = get(results, :setup, "")
    likelihood = get(results, :likelihood, :combined)
    println("Derived-time filter: primary kept $(d.n_kept)/$(d.n_total) " *
            "(logp_window=$logp_window)")

    fig = Figure(size=(720, 420))
    ax  = Axis(fig[1, 1];
               xlabel = "t_0 (kyr)",
               ylabel = "density (1/kyr)",
               title  = "Posterior on t_0  (primary :$likelihood, logp_window=$logp_window)")
    lines!(ax, d.t_grid, d.p_t; linewidth=2, color=:steelblue,
           label="primary (:$likelihood)")
    vlines!(ax, [d.map_t]; color=:red, linewidth=2,
            label="MAP = $(round(Int, d.map_t)) kyr")

    if secondary !== nothing
        d2 = derived_time_distribution(secondary; logp_window=logp_window)
        sec_lik = get(secondary, :likelihood, :unknown)
        println("Derived-time filter: secondary kept $(d2.n_kept)/$(d2.n_total)")
        lines!(ax, d2.t_grid, d2.p_t; linewidth=2, color=:black,
               linestyle=:dash, label="secondary (:$sec_lik)")
    end

    axislegend(ax; position=:rt)

    if save_figures
        mkpath(outdir)
        prefix = joinpath(outdir, string(Dates.today())*"_")
        suffix = isempty(setup) ? "" : "_$setup"
        tag    = secondary === nothing ? "$(likelihood)" :
                 "$(likelihood)-vs-$(get(secondary, :likelihood, :secondary))"
        mysave(prefix*"derived-time-pred-$tag$suffix.png", fig)
    end
    return (fig=fig, t_grid=d.t_grid, p_t=d.p_t, map_t=d.map_t)
end

"""
    plot_derived_time_comparison(primary_path, secondary_path; kwargs...)

Convenience wrapper mirroring `plot_ensemble_comparison` for the t_0 pdf.
"""
function plot_derived_time_comparison(primary_path::AbstractString,
                                      secondary_path::AbstractString;
                                      kwargs...)
    primary   = load_ensemble_results(String(primary_path))
    secondary = load_ensemble_results(String(secondary_path))
    return plot_derived_time(primary; secondary=secondary, kwargs...)
end

# Allow `julia plot_basalmixing_ensemble.jl [path]` for a single chain, OR
# `julia plot_basalmixing_ensemble.jl primary secondary` for an overlay.
# Including this file from the REPL never auto-plots — call plot_ensemble*
# functions yourself.
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) >= 2
        plot_ensemble_comparison(ARGS[1], ARGS[2])
        plot_derived_time_comparison(ARGS[1], ARGS[2])
    else
        path = length(ARGS) >= 1 ? ARGS[1] : "results/emcee-chain-combined.jld2"
        plot_ensemble(path)
        plot_derived_time(load_ensemble_results(path))
    end
end

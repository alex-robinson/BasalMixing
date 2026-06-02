## Plotting for the basal-mixing ensemble run.
##
## Usage (interactive, after running run_basalmixing_ensemble.jl in the same
## session — `chain`, `k81`, `dar40`, `depth`, `priors` already live in `Main`):
##
##     include("plot_basalmixing_ensemble.jl")
##     plot_ensemble()                       # auto-pull from Main
##     # or explicitly:
##     plot_ensemble(; chain, k81, dar40, depth, priors)
##
## Usage (load a previously saved chain):
##
##     include("plot_basalmixing_ensemble.jl")
##     plot_ensemble("results/emcee-chain.jld2")
##
## Usage (one-shot from the shell):
##
##     julia --project=. plot_basalmixing_ensemble.jl results/emcee-chain.jld2

# If invoked as a script (julia plot_basalmixing_ensemble.jl ...), activate the
# project. When `include`d inside an existing session, leave the active env alone.
if abspath(PROGRAM_FILE) == @__FILE__
    cd(@__DIR__)
    import Pkg; Pkg.activate(".")
end

using JLD2
using DataFrames
using CairoMakie
using Dates
using Distributions
using Turing   # loads FlexiChains so JLD2 can reconstruct the chain type on read

include("BasalMixingModel.jl")

function plot_prior_line!(ax, prior::Distribution; color=:red, kwargs...)
    # Works for any Distributions.jl type
    lo, hi = quantile(prior, 0.001), quantile(prior, 0.999)
    x = range(lo, hi, length=300)
    lines!(ax, x, pdf.(prior, x); color=color, kwargs...)
end

"""
    load_ensemble_results(path::String) -> NamedTuple

Read a JLD2 snapshot written by `run_basalmixing_ensemble.jl`. Returns a
NamedTuple with `(chain, k81, dar40, depth, setup, priors, sampler_choice)`.
"""
function load_ensemble_results(path::String)
    isfile(path) || error("ensemble results file not found: $path")
    JLD2.jldopen(path, "r") do f
        chain = f["chain"]
        k81 = f["k81"]
        dar40 = f["dar40"]
        depth = f["depth"]
        setup = haskey(f, "setup") ? f["setup"] : ""
        priors = f["priors"]
        sampler_choice = haskey(f, "sampler_choice") ? f["sampler_choice"] : :unknown
        return (; chain, k81, dar40, depth, setup, priors, sampler_choice)
    end
end

"""
    plot_ensemble(path::String; kwargs...)
    plot_ensemble(; chain, k81, dar40, depth, priors,
                    setup="", sampler_choice=:unknown,
                    outdir="plots", save_figures=true)

Render the three standard ensemble figures and return them as a `NamedTuple
(best, logp, hist)`:

- *best*: depth profiles of mixing rate, ⁸¹Kr closed-system age, and δ⁴⁰Ar
  at the MAP parameters (re-runs the forward model with `dt=0.1, sampling=false`).
- *logp*: per-parameter scatter of log-joint vs. parameter value, MAP highlighted.
- *hist*: marginal posterior histograms overlaid with each parameter's prior.

If `save_figures=true`, the figures are written under `outdir/` with the
current date prefix.

If called with a single string argument, it's interpreted as a JLD2 path; the
file is loaded with `load_ensemble_results` and the resulting NamedTuple is
forwarded to the keyword form.
"""
function plot_ensemble(;
    chain,
    k81,
    dar40,
    depth,
    priors,
    setup::AbstractString="",
    sampler_choice::Symbol=:unknown,
    outdir::AbstractString="plots",
    save_figures::Bool=true,
)
    df = DataFrame(chain)
    df.logjoint = vec(chain[:logjoint])

    params = [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :time_pred]
    labels = ["delta (m)", "m_clean (m/yr)", "f_dirty",
              "t_old (kyr)", "F_ar40 (cc/m²/kyr)", "time_pred (kyr)"]
    best_idx = argmax(df.logjoint)

    ## Best-fit depth profiles ##
    p_best = (
        delta   = df.delta[best_idx],
        m_clean = df.m_clean[best_idx],
        f_dirty = df.f_dirty[best_idx],
        t_old   = df.t_old[best_idx],
        F_ar40  = df.F_ar40[best_idx],
    )
    b = BasalMixingModel(depth=depth, k81_obs_depths=k81.depth, dar40_obs_depths=dar40.depth)
    RunBasalMixingModel!(p_best, b, (k81, dar40); dt=0.1, sampling=false)
    fig_best = plot_BasalMixingModelRun(b; k81_obs=k81, dar40_obs=dar40)

    ## log-posterior scatter ##
    fig_logp = Figure(size=(900, 700))
    ipar = 0
    for (param, label) in zip(params, labels)
        string(param) in names(df) || continue
        ipar += 1
        row, col = divrem(ipar-1, 2)
        ax = Axis(fig_logp[row+1, col+1], xlabel=label, ylabel="log p")
        ylims!(ax, (-50, 0))
        scatter!(ax, df[!, param], df.logjoint;
                 alpha=0.6, markersize=6, color=:steelblue)
        scatter!(ax, [df[best_idx, param]], [df.logjoint[best_idx]];
                 color=:red, markersize=12, label="MAP")
    end

    ## Marginal histograms vs priors ##
    fig_hist = Figure(size=(900, 700))
    ipar = 0
    for (param, label) in zip(params, labels)
        string(param) in names(df) || continue
        ipar += 1
        row, col = divrem(ipar-1, 2)
        ax = Axis(fig_hist[row+1, col+1], xlabel=label, ylabel="density")
        hist!(ax, df[!, param]; bins=20, normalization=:pdf, color=(:steelblue, 0.7))
        if haskey(priors, param) && priors[param] isa Distribution
            plot_prior_line!(ax, priors[param];
                             label="Prior", linewidth=2, color=:grey50)
        end
        vlines!(ax, [df[best_idx, param]]; color=:red, linewidth=2, label="MAP")
    end

    if save_figures
        mkpath(outdir)
        prefix = joinpath(outdir, string(Dates.today())*"_")
        suffix = isempty(setup) ? "" : "_$setup"
        mysave(prefix*"mixingmodel-ens-best$suffix.png", fig_best)
        mysave(prefix*"mixingmodel-ens-logp$suffix.png", fig_logp)
        mysave(prefix*"mixingmodel-ens-hist$suffix.png", fig_hist)
    end

    return (best=fig_best, logp=fig_logp, hist=fig_hist)
end

# Convenience: forward a NamedTuple (e.g. from load_ensemble_results) to the kwarg form.
plot_ensemble(results::NamedTuple; kwargs...) =
    plot_ensemble(;
        chain=results.chain,
        k81=results.k81,
        dar40=results.dar40,
        depth=results.depth,
        priors=results.priors,
        setup=get(results, :setup, ""),
        sampler_choice=get(results, :sampler_choice, :unknown),
        kwargs...,
    )

# Convenience: a path string loads from JLD2 then plots.
plot_ensemble(path::AbstractString; kwargs...) =
    plot_ensemble(load_ensemble_results(String(path)); kwargs...)

"""
    plot_ensemble_from_workspace(; kwargs...)

Pick up `chain`, `k81`, `dar40`, `depth`, `priors` (and optionally `setup`,
`sampler_choice`) from `Main` and forward to `plot_ensemble`. Use this after
running `run_basalmixing_ensemble.jl` interactively:

    include("run_basalmixing_ensemble.jl")
    include("plot_basalmixing_ensemble.jl")
    plot_ensemble_from_workspace()
"""
function plot_ensemble_from_workspace(; kwargs...)
    for sym in (:chain, :k81, :dar40, :depth, :priors)
        isdefined(Main, sym) ||
            error("plot_ensemble_from_workspace: `$(sym)` is not defined in Main; " *
                  "run run_basalmixing_ensemble.jl first, or pass a JLD2 path to plot_ensemble.")
    end
    plot_ensemble(;
        chain  = Main.chain,
        k81    = Main.k81,
        dar40  = Main.dar40,
        depth  = Main.depth,
        priors = Main.priors,
        setup  = isdefined(Main, :setup) ? Main.setup : "",
        sampler_choice = isdefined(Main, :sampler_choice) ? Main.sampler_choice : :unknown,
        kwargs...,
    )
end

# Allow `julia plot_basalmixing_ensemble.jl path/to/chain.jld2`
if abspath(PROGRAM_FILE) == @__FILE__
    path = length(ARGS) >= 1 ? ARGS[1] : "results/emcee-chain.jld2"
    plot_ensemble(path)
end

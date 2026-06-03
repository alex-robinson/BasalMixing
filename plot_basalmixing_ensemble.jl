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
    # Works for any Distributions.jl type
    lo, hi = quantile(prior, 0.001), quantile(prior, 0.999)
    x = range(lo, hi, length=300)
    lines!(ax, x, pdf.(prior, x); color=color, kwargs...)
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
    model_kind::Symbol=:profile,
    outdir::AbstractString="plots",
    save_figures::Bool=true,
)
    df = DataFrame(chain)
    df.logjoint = vec(chain[:logjoint])

    params = [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0]
    labels = ["delta (m)", "m_clean (m/yr)", "f_dirty",
              "t_old (kyr)", "F_ar40 (cc/m²/kyr)", "t_0 (kyr)"]
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
    # Profile run fills b.states (grey time slices) and b.k81/b.dar40 (right panel's
    # dense time series). The predictions at t = sampled MAP t_0 are the same
    # values, just on this dense grid — what changes between :profile and :sampled
    # is only which time we annotate as "best".
    RunBasalMixingModel!(p_best, b, (k81, dar40); dt=0.1, sampling=false)
    if model_kind === :sampled
        # Chain stores t_0 in signed kyr BP; b.*.time_min holds positive
        # elapsed time, so flip the sign before annotating the trajectory.
        t_elapsed_map = -df.t_0[best_idx]
        b.k81.time_min   = t_elapsed_map
        b.dar40.time_min = t_elapsed_map
        b.joint.time_min = t_elapsed_map
    end
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
        model_kind=get(results, :model_kind, :profile),
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
        model_kind = isdefined(Main, :model_kind) ? Main.model_kind : :profile,
        kwargs...,
    )
end

"""
    derived_time_distribution(results::NamedTuple;
                              n_samples::Int=500,
                              t1::Float64=3000.0, dt::Float64=0.2)

Reconstruct the posterior on `t_0` (signed kyr BP, negative = past) from a
chain that doesn't sample it (`model_kind = :marginal`). Implements

  p(t | data) ≈ (1/N) Σ_s L(t | θ_s) / Z(θ_s),  with Z(θ_s) = ∫ L(t | θ_s) dt,

by re-running the forward model at each of `n_samples` chain draws (stride-
thinned across all walkers and iterations), reading the per-slice log-
likelihood from `b.k81.rmse` and `b.dar40.rmse`, normalising each curve to be
a probability density, and averaging. The forward model integrates in positive
elapsed time internally; the returned `t_grid` is sign-flipped to signed kyr
BP so the output is consistent with `:profile` / `:sampled` chains.

For chains that *do* sample `t_0` directly (`:profile`, `:sampled`), the t_0
column itself is the empirical posterior — no re-running needed — and this
function falls back to that.

Returns a NamedTuple `(t_grid, p_t, map_t)` where `t_grid` is the signed kyr
BP grid (negative values, ordered increasing), `p_t` is the density at each
grid point (1/kyr), and `map_t` is the maximum-a-posteriori `t_0`.
"""
function derived_time_distribution(results::NamedTuple;
                                   n_samples::Int=500,
                                   t1::Float64=3000.0,
                                   dt::Float64=0.2)
    chain    = results.chain
    k81      = results.k81
    dar40    = results.dar40
    depth    = results.depth
    model_kind = get(results, :model_kind, :profile)

    df = DataFrame(chain)
    n_total = nrow(df)

    # Fast paths for chains that already carry t_0.
    if model_kind in (:profile, :sampled) && "t_0" in names(df)
        # Histogram in signed kyr BP on the same 1-kyr grid the marginal uses.
        t_grid_full = collect(-t1:1.0:0.0)
        h_counts = zeros(Float64, length(t_grid_full))
        for v in df.t_0
            k = clamp(round(Int, v - t_grid_full[1]) + 1, 1, length(t_grid_full))
            h_counts[k] += 1.0
        end
        p_t = h_counts ./ (sum(h_counts) * 1.0)  # 1/kyr (Δt = 1 kyr)
        map_t = t_grid_full[argmax(p_t)]
        return (t_grid=t_grid_full, p_t=p_t, map_t=map_t)
    end

    # :marginal — average normalised L(t|θ_s) curves over chain draws.
    idx = unique(round.(Int, range(1, n_total; length=min(n_samples, n_total))))
    bs  = [BasalMixingModel(depth=depth, k81_obs_depths=k81.depth, dar40_obs_depths=dar40.depth)
           for _ in 1:Threads.maxthreadid()]

    n_pred       = length(bs[1].k81.time)
    t_grid       = collect(bs[1].k81.time)
    densities    = zeros(n_pred, length(idx))
    n_obs_k81    = nrow(k81)
    n_obs_dar40  = nrow(dar40)

    Threads.@threads for j in eachindex(idx)
        b = bs[Threads.threadid()]
        i = idx[j]
        p = (delta   = df.delta[i],
             m_clean = df.m_clean[i],
             f_dirty = df.f_dirty[i],
             t_old   = df.t_old[i],
             F_ar40  = df.F_ar40[i])
        ok = RunBasalMixingModel!(p, b, (k81, dar40); t1=t1, dt=dt, sampling=false)
        if !ok
            densities[:, j] .= 0.0
            continue
        end
        # Per-slice log-likelihood from RMSE
        log_L = Vector{Float64}(undef, n_pred)
        @inbounds for k in 1:n_pred
            log_L[k] = -0.5 * (n_obs_k81*b.k81.rmse[k]^2 + n_obs_dar40*b.dar40.rmse[k]^2)
        end
        log_Z = logsumexp(log_L)
        if isfinite(log_Z)
            @. densities[:, j] = exp(log_L - log_Z)  # per-slice mass; Δt=1 kyr -> density in 1/kyr
        end
    end

    p_t_elapsed = vec(mean(densities; dims=2))
    # Convert positive-elapsed grid to signed kyr BP, monotonically increasing.
    t_grid_bp   = reverse(-t_grid)
    p_t         = reverse(p_t_elapsed)
    map_t       = t_grid_bp[argmax(p_t)]
    return (t_grid=t_grid_bp, p_t=p_t, map_t=map_t)
end

"""
    plot_derived_time(results; outdir="plots", save_figures=true, kwargs...)

Plot the derived (or sampled) `t_0` posterior. Returns a NamedTuple
`(fig, t_grid, p_t, map_t)`.
"""
function plot_derived_time(results::NamedTuple;
                           outdir::AbstractString="plots",
                           save_figures::Bool=true,
                           kwargs...)
    d = derived_time_distribution(results; kwargs...)
    setup       = get(results, :setup, "")
    model_kind  = get(results, :model_kind, :profile)
    label_kind  = model_kind === :marginal ? "derived" : "sampled"

    fig = Figure(size=(720, 420))
    ax  = Axis(fig[1, 1];
               xlabel = "t_0 (kyr)",
               ylabel = "density (1/kyr)",
               title  = "$(label_kind) posterior on t_0  (model_kind=:$model_kind)")
    lines!(ax, d.t_grid, d.p_t; linewidth=2, color=:steelblue)
    vlines!(ax, [d.map_t]; color=:red, linewidth=2,
            label="MAP = $(round(Int, d.map_t)) kyr")
    axislegend(ax; position=:rt)

    if save_figures
        mkpath(outdir)
        prefix = joinpath(outdir, string(Dates.today())*"_")
        suffix = isempty(setup) ? "" : "_$setup"
        mysave(prefix*"derived-time-pred$suffix.png", fig)
    end
    return (fig=fig, t_grid=d.t_grid, p_t=d.p_t, map_t=d.map_t)
end

# Allow `julia plot_basalmixing_ensemble.jl path/to/chain.jld2`
if abspath(PROGRAM_FILE) == @__FILE__
    path = length(ARGS) >= 1 ? ARGS[1] : "results/emcee-chain.jld2"
    plot_ensemble(path)
end

path = "results/emcee-chain.jld2"
plot_ensemble(path)
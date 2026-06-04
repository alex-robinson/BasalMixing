## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise

using Turing
const AMH = Turing.Inference.AdvancedMH   # AdvancedMH is a Turing dep; reach it via Turing
using LinearAlgebra
using Random

using Dates
using CSV
using DataFrames
using JLD2

Random.seed!(42)

# Include code for basal mixing model
include("BasalMixingModel.jl")

"""
    tuned_mh(priors; var_linked=0.1, variances=NamedTuple())

Build a Turing `MH` sampler with a `LinkedRW(σ²)` proposal for every prior in
`priors` that is a `Distribution`. `LinkedRW` performs the random walk in
unconstrained (linked) space, so bounded priors (Uniform, truncated) stay in
support automatically. `σ²` is the proposal variance *in linked space*.

For Uniform(a,b) the link is logit((x-a)/(b-a)); for truncated lower-bounded
distributions the link is roughly log(x-lower). A `var_linked` of 0.05-0.25 is
a reasonable starting range; tune toward ~25-40% acceptance.

Override individual parameters via `variances=(:delta=0.04, ...)`.
"""
function tuned_mh(priors; var_linked::Float64=0.1, variances::NamedTuple=NamedTuple())
    props = Pair{Symbol,Any}[]
    LinkedRW = Turing.Inference.LinkedRW
    for name in keys(priors)
        prior = getfield(priors, name)
        prior isa Distribution || continue
        v = haskey(variances, name) ? variances[name] : var_linked
        push!(props, name => LinkedRW(v))
    end
    return MH(props...)
end

"""
    acceptance_rate(chain)

Mean acceptance rate. Uses the chain's `:accepted` flag if present and
populated (FlexiChain default in Turing >= 0.42); otherwise falls back to
counting non-repeated rows in the raw sample matrix. Returns a `Float64`,
never `missing` — some samplers (notably `Emcee`) store `Vector{Missing}`
for `:accepted` and we don't want that to crash downstream printing.
"""
function acceptance_rate(chain)
    try
        a = chain[:accepted]
        v = collect(skipmissing(vec(a)))
        if !isempty(v)
            return mean(v)
        end
    catch
        # fall through to row-diff fallback below
    end
    θ = Array(chain)
    if ndims(θ) >= 2 && size(θ, 1) > 1
        return mean(any(diff(θ, dims=1) .!= 0, dims=2))
    end
    return NaN
end

priors = (
    delta       = Uniform(0.3, 2.0),    # 1 m
    m_clean     = truncated(Normal(0.03, 0.002), lower=0.0), # 0.03 m/kyr
    f_dirty     = Uniform(4.0, 6.5),                        # 0.18 m/kyr / 0.03 m/kyr => f=6x
    t_old       = truncated(Normal(250.0,25.0), lower=0.0),  # 250 kyr
    F_ar40      = Uniform(0.004,0.007), #Normal(0.075,0.01),   # 0.075 m^3 / kyr
    σ_k81 = 30.0,
    σ_dar40 = 0.03,

    # Inception time t_0 in kyr (negative = past). Magnitude |t_0| is the
    # forward-model run length 0 → -t_0. Bounds [-3000, -700] kyr come from
    # the geological window for first GrIS glaciation (~2700 kyr) and the
    # ⁸¹Kr-implied minimum summit ice cover (~800 kyr).
    t_0         = Uniform(-3000.0, -700.0),
)

"""
    basal_mixing(k81_age_obs, dar40_obs, bs, dat, dt, priors, likelihood)

Sampled-time Bayesian model. Draws `t_0` from `priors.t_0` (signed kyr,
negative = past) and integrates the forward model exactly once over duration
|t_0|, then compares the predicted endpoint state to the observations.

`likelihood` selects which observation channel(s) constrain the posterior:
- `:combined` — both ⁸¹Kr ages and δ⁴⁰Ar (the default standard run)
- `:kr81`     — only ⁸¹Kr ages
- `:ar40`     — only δ⁴⁰Ar
"""
@model function basal_mixing(k81_age_obs, dar40_obs, bs, dat, dt, priors, likelihood::Symbol)

    ## Set priors ##
    delta     = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean   = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty   = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old     = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40    = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40
    σ_k81     = priors.σ_k81   isa Distribution ? σ_k81   ~ priors.σ_k81   : priors.σ_k81
    σ_dar40   = priors.σ_dar40 isa Distribution ? σ_dar40 ~ priors.σ_dar40 : priors.σ_dar40
    t_0       = priors.t_0     isa Distribution ? t_0     ~ priors.t_0     : priors.t_0

    # Pick this thread's private model state
    b = bs[Threads.threadid()]

    # Run the forward model exactly once over positive elapsed time 0 → -t_0.
    duration = -t_0
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40)
    success = RunBasalMixingModelToTime!(p, b, duration, dat; dt=dt)

    if !success
        success = RunBasalMixingModelToTime!(p, b, duration, dat; dt=dt*0.5)
    end

    if success
        k81_age_pred := b.k81.dat[:, 1]
        dar40_pred   := b.dar40.dat[:, 1]
    else
        # Predictions far from any plausible obs → very low likelihood
        k81_age_pred := fill(1e8, length(k81_age_obs))
        dar40_pred   := fill(1e8, length(dar40_obs))
    end

    if likelihood === :combined
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    elseif likelihood === :kr81
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
    elseif likelihood === :ar40
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    else
        error("Unknown likelihood: $likelihood (expected :combined, :kr81, or :ar40)")
    end

    return
end

## SCRIPT ##

depth, setup = generate_depths("highdirty";step=0.25)
(k81, dar40) = load_basalmixing_data(depth=depth)

# Allocate one BasalMixingModel per thread to avoid concurrent mutation of shared buffers
# when MCMCThreads() runs chains in parallel. Size by maxthreadid() — Julia >= 1.9 may
# return threadids beyond nthreads(:default) for interactive threads, and tasks can
# migrate between threads, so the indexed slot must always exist.
bs = [BasalMixingModel(depth=depth, k81_obs_depths=k81.depth, dar40_obs_depths=dar40.depth)
      for _ in 1:Threads.maxthreadid()]

## Sampler selection ##
# :emcee       -> Affine-invariant ensemble sampler (Goodman & Weare). Robust to
#                 correlations, no per-parameter tuning, no gradients required.
#                 Default.
# :mh_ram      -> Robust Adaptive Metropolis (Vihola 2012): online-adapts the
#                 full proposal covariance toward target acceptance α.
# :mh_tuned    -> MH with per-parameter LinkedRW proposals.
# :mh_default  -> Turing's default MH (uses priors as proposals; poor mixing here).
sampler_choice = :emcee

n_samples = 10_000   # used for the MH variants; per-chain count
n_chains  = 4

# Optional per-parameter linked-space variance overrides for :mh_tuned.
mh_variances = NamedTuple()

## Likelihood configurations to run back-to-back.
# Each entry produces a separate JLD2 snapshot under results/.
likelihoods = (:combined, :kr81)

results = Dict{Symbol,Any}()

for likelihood in likelihoods

    println("\n========================================")
    println("Running likelihood = :$likelihood")
    println("========================================")

    model = basal_mixing(k81.age, dar40.dar40, bs, (k81, dar40), 0.2, priors, likelihood)

    if sampler_choice === :mh_default
        spl = MH()
        chain = sample(model, spl, MCMCThreads(), n_samples, n_chains)
    elseif sampler_choice === :mh_tuned
        spl = tuned_mh(priors; var_linked=0.1, variances=mh_variances)
        chain = sample(model, spl, MCMCThreads(), n_samples, n_chains)
    elseif sampler_choice === :mh_ram
        spl = externalsampler(
            AMH.RobustAdaptiveMetropolis(; α=0.234, γ=2/3);
            adtype=Turing.ADTypes.AutoFiniteDiff(),
        )
        chain = sample(model, spl, MCMCThreads(), n_samples, n_chains)
    elseif sampler_choice === :emcee
        n_walkers = 32
        n_iter    = 500
        spl = Emcee(n_walkers, 2.0)
        chain = sample(model, spl, n_iter)
    else
        error("Unknown sampler_choice: $sampler_choice")
    end

    println("Sampler: $sampler_choice  |  acceptance ≈ ",
            round(acceptance_rate(chain), digits=3))

    results_path = "results/emcee-chain-$(likelihood).jld2"
    save_ensemble_results(results_path;
                          chain, k81, dar40, depth, setup, priors,
                          sampler_choice, likelihood)

    df = DataFrame(chain)
    df.logjoint = vec(chain[:logjoint])
    best_idx = argmax(df.logjoint)
    map_cols = filter(c -> string(c) in names(df),
                      [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0])
    @info "MAP for :$likelihood (joint logp = $(round(maximum(df.logjoint); digits=2)))" df[best_idx, map_cols]
    for p in map_cols
        v = df[!, p]
        qs = quantile(v, [0.025, 0.5, 0.975])
        @info "$(p): q025=$(round(qs[1]; sigdigits=4)) median=$(round(qs[2]; sigdigits=4)) q975=$(round(qs[3]; sigdigits=4)) std=$(round(std(v); sigdigits=4))"
    end

    results[likelihood] = (; chain, results_path)
end

## Plotting (overlays kr81-only on combined):
##
##     include("plot_basalmixing_ensemble.jl")
##     plot_ensemble_comparison("results/emcee-chain-combined.jld2",
##                              "results/emcee-chain-kr81.jld2")

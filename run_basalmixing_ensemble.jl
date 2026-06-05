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
    # Widened delta and f_dirty (2026-06-05) — the original (0.3, 2.0) and
    # (4.0, 6.5) bounds were clipping the posterior: both Emcee and NUTS
    # piled up on their upper edges. New bounds keep the priors physically
    # plausible (higher values are unphysical) while letting the typical
    # set sit in the interior, which NUTS needs for sane step-size adaptation.
    delta       = Uniform(0.3, 3.0),    # transition-zone thickness (m)
    m_clean     = truncated(Normal(0.03, 0.002), lower=0.0), # 0.03 m/kyr
    f_dirty     = Uniform(3.0, 8.0),    # dirty/clean mixing-rate ratio
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
    basal_mixing_nuts(k81_age_obs, dar40_obs, geom, priors, likelihood)

NUTS-compatible Bayesian model. Same structure as `basal_mixing`, but uses
the ODE-based forward integrator `run_basal_mixing_ode` (Tsit5 with adaptive
stepping) so the likelihood is differentiable w.r.t. all parameters
including `t_0`. `geom` is a `BasalMixingGeometry` NamedTuple — pass the one
allocated for the depth grid; no per-thread buffer pool is needed because
the ODE allocates parametrically per call.

The Euler-based `basal_mixing` below is kept for the Emcee path. Both share
the same prior block and likelihood-gating logic.
"""
@model function basal_mixing_nuts(k81_age_obs, dar40_obs, geom, priors, likelihood::Symbol)

    ## Set priors ##
    delta     = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean   = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty   = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old     = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40    = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40
    σ_k81     = priors.σ_k81   isa Distribution ? σ_k81   ~ priors.σ_k81   : priors.σ_k81
    σ_dar40   = priors.σ_dar40 isa Distribution ? σ_dar40 ~ priors.σ_dar40 : priors.σ_dar40
    t_0       = priors.t_0     isa Distribution ? t_0     ~ priors.t_0     : priors.t_0

    duration = -t_0
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty,
         t_old=t_old, F_ar40=F_ar40, t_target=duration)
    # Loose ODE tolerance — 1e-4 still ~100× tighter than σ_k81=30 kyr and
    # σ_dar40=0.03 ‰. The probe with this setting agreed with Euler to
    # within max |Δage| = 0.14 kyr, max |Δdar40| = 4.5e-4 ‰.
    k81_age_pred_v, dar40_pred_v, ok = run_basal_mixing_ode(p, geom;
                                                            reltol=1e-4, abstol=1e-7)

    if ok
        k81_age_pred := k81_age_pred_v
        dar40_pred   := dar40_pred_v
    else
        # ODE failed (very rare under sensible params); drive the likelihood
        # to ~-Inf with predictions far from observations.
        k81_age_pred := fill(typeof(k81_age_pred_v[1])(1e8), length(k81_age_obs))
        dar40_pred   := fill(typeof(dar40_pred_v[1])(1e8),   length(dar40_obs))
    end

    if likelihood === :combined
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    elseif likelihood === :kr81
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
    elseif likelihood === :ar40
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    else
        error("Unknown likelihood: $likelihood")
    end

    return
end

"""
    basal_mixing(k81_age_obs, dar40_obs, bs, dat, dt, priors, likelihood)

Emcee-compatible (gradient-free) Bayesian model. Uses the hand-rolled Euler
integrator `RunBasalMixingModelToTime!` with a per-thread buffer pool.

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
# :nuts        -> No-U-Turn Sampler via Turing's `NUTS`. Requires the
#                 ODE-based forward model (`run_basal_mixing_ode`) and
#                 ForwardDiff. The right choice when you want a true Bayesian
#                 posterior with no stuck walkers.
# :emcee       -> Affine-invariant ensemble sampler (Goodman & Weare). Robust to
#                 correlations, no per-parameter tuning, no gradients required.
#                 Watch for boundary-stuck walkers — see derived_time histogram.
# :mh_ram      -> Robust Adaptive Metropolis (Vihola 2012): online-adapts the
#                 full proposal covariance toward target acceptance α.
# :mh_tuned    -> MH with per-parameter LinkedRW proposals.
# :mh_default  -> Turing's default MH (uses priors as proposals; poor mixing here).
sampler_choice = :nuts

n_samples = 10_000   # used for the MH variants; per-chain count
n_chains  = 4

# NUTS-specific tuning: adaptation and post-adapt sample counts (per chain).
# Tuned (2026-06-05) so each chain fits inside the harness ~25 min kill window:
# probe data shows ~6.4 s/draw post-smoothing, so 250 draws ≈ 26 min single-chain.
# 4 chains × MCMCThreads in parallel → 600 post-adapt samples total at ~25 min wall.
n_adapts_nuts = 100
n_samples_nuts = 150
target_accept_nuts = 0.65

# Precompute the ODE-path geometry (Float64-only depth/thickness data).
# Cheap to allocate; reused by every NUTS likelihood eval.
geom = BasalMixingGeometry(depth, 3040.0, 3053.44, k81.depth, dar40.depth)

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

    if sampler_choice === :nuts
        model = basal_mixing_nuts(k81.age, dar40.dar40, geom, priors, likelihood)
    else
        model = basal_mixing(k81.age, dar40.dar40, bs, (k81, dar40), 0.2, priors, likelihood)
    end

    if sampler_choice === :nuts
        # max_depth=8 (256 leapfrogs/draw worst case) caps adaptation cost.
        # The AD-timing diagnostic shows ~60 ms/gradient → ≤15 s/draw at depth 8.
        # Default max_depth=10 lets adaptation drift to 60 s/draw, which is what
        # the killed smoke test was doing.
        spl = NUTS(n_adapts_nuts, target_accept_nuts;
                   max_depth=8,
                   adtype=Turing.ADTypes.AutoForwardDiff())
        chain = sample(model, spl, MCMCThreads(), n_samples_nuts, n_chains;
                       progress=true)
    elseif sampler_choice === :mh_default
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

    results_path = "results/$(sampler_choice)-chain-$(likelihood).jld2"
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

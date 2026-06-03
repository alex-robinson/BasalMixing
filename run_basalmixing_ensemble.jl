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
        # FlexiChain returns a DimArray; some samplers populate it with `missing`.
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
    #σ_k81           = Exponential(30.0),    # 30 kyr
    #t_old = 250.0,
    σ_k81 = 30.0,
    σ_dar40 = 0.03,

    # Inception time t_0 in kyr (negative = past). Magnitude |t_0| is the
    # forward-model run length 0 → -t_0. Bounds [-3000, -700] kyr come from
    # the geological window for first GrIS glaciation (~2700 kyr) and the
    # ⁸¹Kr-implied minimum summit ice cover (~800 kyr).
    t_0         = Uniform(-3000.0, -700.0),
)

@model function basal_mixing(k81_age_obs, dar40_obs, bs, dat, dt, priors)

    ## Set priors ##
    delta   = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old   = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40  = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40
    σ_k81   = priors.σ_k81   isa Distribution ? σ_k81   ~ priors.σ_k81   : priors.σ_k81
    σ_dar40 = priors.σ_dar40 isa Distribution ? σ_dar40 ~ priors.σ_dar40 : priors.σ_dar40

    # Pick this thread's private model state so concurrent chains don't trash each other's buffers.
    b = bs[Threads.threadid()]

    # Run the model
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40)
    success = RunBasalMixingModel!(p, b, dat; dt=dt, sampling=true)

    if !success
        # Try one more time with a smaller timestep
        success = RunBasalMixingModel!(p, b, dat; dt=dt*0.5, sampling=true)
    end

    if success
        # Extract the best-fit ages at the optimal time. b.joint.time_min is
        # positive elapsed kyr; the chain reports t_0 in signed kyr BP.
        kmin = b.joint.kmin
        t_0 := -b.joint.time_min
        k81_age_pred  := b.k81.dat[:,kmin]
        dar40_pred := b.dar40.dat[:,kmin]
    else
        # Assign a very high age so that Likelihood is very low
        t_0 := 0.0
        k81_age_pred  := fill(1e8, length(k81_age_obs))
        dar40_pred := fill(1e8, length(dar40_obs))
    end

    # Likelihoods:
    #   - observed k81 ages ~ Normal(predicted, σ_k81)
    #   - observed dar40 values ~ Normal(predicted, σ_dar40)

    k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
    dar40_obs ~ MvNormal(dar40_pred, σ_dar40 * I)

    return
end

"""
    basal_mixing_sampled(k81_age_obs, dar40_obs, bs, dat, dt, priors)

Sampled-time companion to `basal_mixing`: instead of running the forward model
to `t1=3000` and selecting `argmin(rmse)` (a profile likelihood), this model
draws `t_0` from `priors.t_0` (signed kyr, negative = past) and integrates the
forward model exactly once, from t_0 to t_end=0 (present day). Internally the
integrator runs on positive elapsed time `0 → -t_0`; the dynamics are
autonomous so this relabeling is exact. The likelihood compares the
observations to the predicted final state at t_end=0.

This is the proper Bayesian setup — the posterior on `t_0` is an actual
marginal, not a deterministic side-effect of point estimation.
"""
@model function basal_mixing_sampled(k81_age_obs, dar40_obs, bs, dat, dt, priors)

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

    # Run the forward model exactly once, integrator runs 0 → -t_0 in positive
    # elapsed time. RunBasalMixingModelToTime! stores the (positive) duration
    # into b.*.time_min; downstream plotting flips it back to signed kyr BP.
    duration = -t_0
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40)
    success = RunBasalMixingModelToTime!(p, b, duration, dat; dt=dt)

    if !success
        # Retry with a smaller dt before giving up
        success = RunBasalMixingModelToTime!(p, b, duration, dat; dt=dt*0.5)
    end

    if success
        # Final-state predictions live in slot 1 (set by RunBasalMixingModelToTime!)
        k81_age_pred := b.k81.dat[:, 1]
        dar40_pred   := b.dar40.dat[:, 1]
    else
        # Force a very low likelihood by returning predictions far from any plausible obs
        k81_age_pred := fill(1e8, length(k81_age_obs))
        dar40_pred   := fill(1e8, length(dar40_obs))
    end

    k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
    dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)

    return
end

"""
    basal_mixing_marginal(k81_age_obs, dar40_obs, bs, dat, dt, priors)

Marginal-time companion. `t_0` is *not* sampled; instead the forward model is
integrated 0 → 3000 kyr (positive elapsed) and the likelihood is
`log ∫ L(t | θ, data) dt`, treating elapsed time as a nuisance parameter under
a Uniform(t0, t1) prior. This is what we want when the data doesn't constrain
`t_0` past the equilibration time — see `docs/basal_mixing.tex` for the
derivation.

The `k81_age_obs` and `dar40_obs` arguments are kept for signature symmetry
with the other two models; they're not used directly here because the
log-likelihood is computed inside `RunBasalMixingModelMarginal!` against the
obs DataFrames in `dat`.
"""
@model function basal_mixing_marginal(k81_age_obs, dar40_obs, bs, dat, dt, priors)

    delta   = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old   = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40  = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40

    b = bs[Threads.threadid()]

    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40)
    log_Z = RunBasalMixingModelMarginal!(p, b, dat; dt=dt)

    if !isfinite(log_Z)
        # DomainError fallback: retry once at half dt.
        log_Z = RunBasalMixingModelMarginal!(p, b, dat; dt=dt*0.5)
    end

    # Add the marginal log-likelihood directly. No likelihood ~ statement is needed
    # because the observations have already been integrated against the predictions
    # inside RunBasalMixingModelMarginal!.
    Turing.@addlogprob! log_Z

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

## Model selection ##
# All three modes report `t_0` in *signed kyr BP* (negative = past, 0 = today).
#
# :profile  -> RunBasalMixingModel! integrates 0 -> t1=3000 kyr (positive
#              elapsed time) and the model reports the argmin-RMSE slice;
#              `t_0` is derived as -b.joint.time_min.
# :sampled  -> RunBasalMixingModelToTime! integrates the forward model once
#              over duration |t_0|; `t_0` is a proper sampled parameter with
#              the priors.t_0 prior (Uniform on signed kyr BP).
# :marginal -> RunBasalMixingModelMarginal! integrates 0 -> 3000 kyr and
#              computes log ∫ L(t | θ) dt with equilibration-tail shortcut.
#              `t_0` is a nuisance variable, integrated out — not sampled.
#              5-D MCMC (vs 6-D for :sampled). See docs/basal_mixing.tex.
model_kind = :sampled

if model_kind === :profile
    model = basal_mixing(k81.age, dar40.dar40, bs, (k81, dar40), 0.2, priors)
elseif model_kind === :sampled
    model = basal_mixing_sampled(k81.age, dar40.dar40, bs, (k81, dar40), 0.2, priors)
elseif model_kind === :marginal
    model = basal_mixing_marginal(k81.age, dar40.dar40, bs, (k81, dar40), 0.2, priors)
else
    error("Unknown model_kind: $model_kind")
end

## Sampler selection ##
# :emcee       -> Affine-invariant ensemble sampler (Goodman & Weare). Robust to
#                 correlations, no per-parameter tuning, no gradients required.
#                 In informal comparisons here, the only sampler that hits
#                 textbook ~25-40% acceptance and fills out the posterior. Slower
#                 per sample than the MH variants (serial over walkers), but
#                 wall-clock comparable for the same effective sample size.
#                 Default.
# :mh_ram      -> Robust Adaptive Metropolis (Vihola 2012): online-adapts the
#                 full proposal covariance toward target acceptance α. Needs a
#                 calibrated initial S (or many iterations) to be competitive;
#                 with default S=I its acceptance can stay near 0 for thousands
#                 of iterations.
# :mh_tuned    -> MH with per-parameter LinkedRW proposals. Tune `var_linked`
#                 (or pass `mh_variances=...`) toward ~25-40% acceptance.
# :mh_default  -> Turing's default MH (uses priors as proposals; explores poorly here).
sampler_choice = :emcee

n_samples = 10_000   # used for the MH variants; per-chain count
n_chains  = 4

# Where to persist the sampled chain (set to nothing to skip saving).
# Default name reflects model_kind so :sampled doesn't overwrite the :profile snapshot.
results_path = "results/emcee-chain-$(model_kind).jld2"

# Optional per-parameter linked-space variance overrides for :mh_tuned
# (target ~25-40% acceptance). Leave NamedTuple() to use var_linked for all.
mh_variances = NamedTuple()
# Example: mh_variances = (delta=0.05, m_clean=0.1, f_dirty=0.1, t_old=0.1, F_ar40=0.1)

if sampler_choice === :mh_default
    spl = MH()
    chain = sample(model, spl, MCMCThreads(), n_samples, n_chains)
elseif sampler_choice === :mh_tuned
    spl = tuned_mh(priors; var_linked=0.1, variances=mh_variances)
    chain = sample(model, spl, MCMCThreads(), n_samples, n_chains)
elseif sampler_choice === :mh_ram
    # α=0.234 is the asymptotically optimal RWMH acceptance rate for d>=5.
    # γ in (0.5, 1] controls the adaptation step; 2/3 is the default.
    # AutoFiniteDiff is used only for the initial-param-finding check (RAM
    # itself is gradient-free). ForwardDiff would propagate Duals through the
    # model's try/catch + mutating buffers and fail; FiniteDiff just calls the
    # log-density at real-valued points.
    spl = externalsampler(
        AMH.RobustAdaptiveMetropolis(; α=0.234, γ=2/3);
        adtype=Turing.ADTypes.AutoFiniteDiff(),
    )
    chain = sample(model, spl, MCMCThreads(), n_samples, n_chains)
elseif sampler_choice === :emcee
    # Emcee runs a single ensemble of walkers; n_walkers >= 2 * n_sampled_params.
    # Total draws = n_walkers * n_iter. Currently no MCMCThreads() support in
    # Turing's Emcee, so wall-clock is serial over walkers.
    n_walkers = 32
    n_iter    = 500
    spl = Emcee(n_walkers, 2.0)
    chain = sample(model, spl, n_iter)
else
    error("Unknown sampler_choice: $sampler_choice")
end

println("Sampler: $sampler_choice  |  acceptance ≈ ", round(acceptance_rate(chain), digits=3))

## Persist the chain for downstream plotting / reanalysis.
save_ensemble_results(results_path;
                 chain, k81, dar40, depth, setup, priors,
                 sampler_choice, model_kind)

## Quick summary (full plotting lives in plot_basalmixing_ensemble.jl):
##
##     include("plot_basalmixing_ensemble.jl")
##     plot_ensemble_from_workspace()         # uses chain/k81/dar40/depth/priors from Main
##     # or, separately, after restarting Julia:
##     plot_ensemble("results/emcee-chain.jld2")

df = DataFrame(chain)
df.logjoint = vec(chain[:logjoint])
best_idx = argmax(df.logjoint)
# t_0 only exists in :profile / :sampled chains; skip it for :marginal.
map_cols = filter(c -> string(c) in names(df), [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0])
@info "MAP (joint logp = $(round(maximum(df.logjoint); digits=2)))" df[best_idx, map_cols]
# Per-parameter quantile summary (FlexiChain doesn't implement StatsBase.describe in
# Turing 0.45; fall back to a quick DataFrame summary).
for p in map_cols
    v = df[!, p]
    qs = quantile(v, [0.025, 0.5, 0.975])
    @info "$(p): q025=$(round(qs[1]; sigdigits=4)) median=$(round(qs[2]; sigdigits=4)) q975=$(round(qs[3]; sigdigits=4)) std=$(round(std(v); sigdigits=4))"
end

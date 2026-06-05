## Probe NUTS run: 200 adapt + 100 sample, single chain, :combined.
##
## Goal: see whether the prior-boundary parking on `delta`/`f_dirty` from
## the smoke test (50 adapt only) clears once adaptation is given a real
## budget. If it does → kick off the full 4-chain production run. If it
## doesn't → widen those priors first.

cd(@__DIR__)
import Pkg; Pkg.activate(".")

log_line(s) = (println(s); flush(stdout))

log_line("[startup] loading…")
using Turing, DataFrames, CSV, JLD2
using LinearAlgebra, Random, Statistics
using DifferentialEquations, ForwardDiff
Random.seed!(42)

include("BasalMixingModel.jl")

depth, setup = generate_depths("highdirty"; step=0.25)
(k81, dar40) = load_basalmixing_data(depth=depth)
geom = BasalMixingGeometry(depth, 3040.0, 3053.44, k81.depth, dar40.depth)

priors = (
    delta   = Uniform(0.3, 3.0),                   # widened (was 0.3..2.0)
    m_clean = truncated(Normal(0.03, 0.002), lower=0.0),
    f_dirty = Uniform(3.0, 8.0),                   # widened (was 4.0..6.5)
    t_old   = truncated(Normal(250.0,25.0), lower=0.0),
    F_ar40  = Uniform(0.004, 0.007),
    σ_k81   = 30.0,
    σ_dar40 = 0.03,
    t_0     = Uniform(-3000.0, -700.0),
)

@model function basal_mixing_nuts_local(k81_age_obs, dar40_obs, geom, priors, likelihood::Symbol)
    delta     = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean   = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty   = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old     = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40    = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40
    σ_k81     = priors.σ_k81   isa Distribution ? σ_k81   ~ priors.σ_k81   : priors.σ_k81
    σ_dar40   = priors.σ_dar40 isa Distribution ? σ_dar40 ~ priors.σ_dar40 : priors.σ_dar40
    t_0       = priors.t_0     isa Distribution ? t_0     ~ priors.t_0     : priors.t_0
    duration = -t_0
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40, t_target=duration)
    # Loosened ODE tolerance: obs σ_k81=30 kyr, σ_dar40=0.03 ‰. reltol=1e-4
    # is still ~100× tighter than data precision, but should let Tsit5 take
    # larger steps than reltol=1e-6 default.
    k81_age_pred_v, dar40_pred_v, ok = run_basal_mixing_ode(p, geom;
                                                            reltol=1e-4, abstol=1e-7)
    if ok
        k81_age_pred := k81_age_pred_v
        dar40_pred   := dar40_pred_v
    else
        k81_age_pred := fill(typeof(k81_age_pred_v[1])(1e8), length(k81_age_obs))
        dar40_pred   := fill(typeof(dar40_pred_v[1])(1e8),   length(dar40_obs))
    end
    if likelihood === :combined
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    elseif likelihood === :kr81
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
    end
    return
end

model = basal_mixing_nuts_local(k81.age, dar40.dar40, geom, priors, :combined)

log_line("[NUTS] starting: 100 adapt + 50 sample, max_depth=8")
spl = NUTS(100, 0.65; max_depth=8, adtype=AutoForwardDiff())

# Instead of a callback (whose signature varies by Turing version), count
# log-density evals at the model-eval level. We can't see "iterations" from
# inside the @model, but we can see throughput: NUTS with avg tree depth d
# costs ~2^d log-density evals per draw, so eval rate × constant ≈ draw rate.
const EVAL_COUNTER = Ref(0)
const LAST_LOG_T   = Ref(time())
const LAST_LOG_E   = Ref(0)

@model function basal_mixing_nuts_traced(k81_age_obs, dar40_obs, geom, priors, likelihood::Symbol)
    delta     = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean   = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty   = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old     = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40    = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40
    σ_k81     = priors.σ_k81   isa Distribution ? σ_k81   ~ priors.σ_k81   : priors.σ_k81
    σ_dar40   = priors.σ_dar40 isa Distribution ? σ_dar40 ~ priors.σ_dar40 : priors.σ_dar40
    t_0       = priors.t_0     isa Distribution ? t_0     ~ priors.t_0     : priors.t_0
    duration = -t_0
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40, t_target=duration)
    k81_age_pred_v, dar40_pred_v, ok = run_basal_mixing_ode(p, geom; reltol=1e-4, abstol=1e-7)
    if ok
        k81_age_pred := k81_age_pred_v
        dar40_pred   := dar40_pred_v
    else
        k81_age_pred := fill(typeof(k81_age_pred_v[1])(1e8), length(k81_age_obs))
        dar40_pred   := fill(typeof(dar40_pred_v[1])(1e8),   length(dar40_obs))
    end
    if likelihood === :combined
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    elseif likelihood === :kr81
        k81_age_obs ~ MvNormal(k81_age_pred, σ_k81 * I)
    end

    # End-of-eval throughput log. Every 100 evals, dump (evals, s/eval) so we
    # can spot whether NUTS is making progress.
    EVAL_COUNTER[] += 1
    if EVAL_COUNTER[] % 100 == 0
        now = time()
        δe = EVAL_COUNTER[] - LAST_LOG_E[]
        δt = now - LAST_LOG_T[]
        log_line("[evals=$(EVAL_COUNTER[])] +$δe evals in $(round(δt;digits=1)) s " *
                 "($(round(1000*δt/δe;digits=1)) ms/eval)")
        LAST_LOG_T[] = now
        LAST_LOG_E[] = EVAL_COUNTER[]
    end

    return
end

# Use the traced model instead of the original local one.
model = basal_mixing_nuts_traced(k81.age, dar40.dar40, geom, priors, :combined)

t_chain = @elapsed begin
    chain = sample(model, spl, 50; progress=false)
end
log_line("[NUTS] done in $(round(t_chain;digits=1)) s ($(round(t_chain/50*1000;digits=0)) ms/sample); total evals=$(EVAL_COUNTER[])")

# Save for downstream plotting / inspection.
mkpath("results")
save_ensemble_results("results/nuts-probe-combined-widened.jld2";
                      chain, k81, dar40, depth, setup, priors,
                      sampler_choice=:nuts, likelihood=:combined)

# Summary
try; println("Mean accept:     ", round(mean(skipmissing(vec(chain[:acceptance_rate]))); digits=3)); catch; end
try; println("Mean tree depth: ", round(mean(skipmissing(vec(chain[:tree_depth]))); digits=2)); catch; end
println("logp range: ", round(minimum(vec(chain[:logjoint]));digits=2),
        " to ", round(maximum(vec(chain[:logjoint]));digits=2))

df = DataFrame(chain)
df.logjoint = vec(chain[:logjoint])
best_idx = argmax(df.logjoint)
println("\nMAP: logp=$(round(df.logjoint[best_idx]; digits=2))")
for p in [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0]
    string(p) in names(df) || continue
    v = df[!, p]
    println("  $(p): MAP=$(round(df[best_idx,p];sigdigits=4))  " *
            "median=$(round(median(v);sigdigits=4))  " *
            "q025=$(round(quantile(v,0.025);sigdigits=4))  " *
            "q975=$(round(quantile(v,0.975);sigdigits=4))")
end

# Boundary check
println("\nBoundary parking check (fraction of post-adapt samples within 1% of each prior edge):")
for p in [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0]
    string(p) in names(df) || continue
    prior = priors[p]
    if prior isa Distribution
        lo, hi = try (minimum(prior), maximum(prior)); catch; (-Inf, Inf); end
        if isfinite(lo) && isfinite(hi)
            span = hi - lo
            frac_lo = mean(df[!, p] .< lo + 0.01*span)
            frac_hi = mean(df[!, p] .> hi - 0.01*span)
            tag = (frac_lo + frac_hi) > 0.1 ? "  ⚠ pile-up" : ""
            println("  $(p): $(round(frac_lo*100;digits=1))% at lower, $(round(frac_hi*100;digits=1))% at upper$tag")
        end
    end
end

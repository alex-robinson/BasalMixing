## Smoke test for the NUTS path: confirm ForwardDiff works through the ODE,
## time a likelihood eval and a gradient eval, then run a tiny NUTS chain.

cd(@__DIR__)
import Pkg; Pkg.activate(".")

# Flushing print helper — macOS lacks `stdbuf`/`script -F`, so we make every
# log line flush stdout ourselves. Without this, redirected output stays
# block-buffered and a 6-minute Turing precompile looks indistinguishable
# from a hang.
log_line(s) = (println(s); flush(stdout))

log_line("[$(round(time();digits=1))] startup: loading Turing…")

using Turing
using DataFrames, CSV, JLD2
using LinearAlgebra
using Random
using Statistics
using DifferentialEquations
using ForwardDiff

log_line("[$(round(time();digits=1))] including BasalMixingModel.jl…")
include("BasalMixingModel.jl")

Random.seed!(42)

log_line("[$(round(time();digits=1))] loading data, building geometry…")
depth, setup = generate_depths("highdirty"; step=0.25)
(k81, dar40) = load_basalmixing_data(depth=depth)
geom = BasalMixingGeometry(depth, 3040.0, 3053.44, k81.depth, dar40.depth)

priors = (
    delta   = Uniform(0.3, 2.0),
    m_clean = truncated(Normal(0.03, 0.002), lower=0.0),
    f_dirty = Uniform(4.0, 6.5),
    t_old   = truncated(Normal(250.0,25.0), lower=0.0),
    F_ar40  = Uniform(0.004, 0.007),
    σ_k81   = 30.0,
    σ_dar40 = 0.03,
    t_0     = Uniform(-3000.0, -700.0),
)

# Re-include the @model from the runner — easier to keep one source of truth,
# but we just paste the model here to keep the smoke test self-contained.
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
    k81_age_pred_v, dar40_pred_v, ok = run_basal_mixing_ode(p, geom)
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
    elseif likelihood === :ar40
        dar40_obs   ~ MvNormal(dar40_pred,   σ_dar40 * I)
    end
    return
end

log_line("[$(round(time();digits=1))] instantiating model…")
model = basal_mixing_nuts_local(k81.age, dar40.dar40, geom, priors, :combined)

## --- Tiny NUTS chain (combines forward + gradient timing implicitly) ---
# max_depth=6 caps leapfrog steps per draw at 64 (vs default 1024) — the AD
# timing diagnostic shows ~60 ms/gradient, so depth 6 ≈ 3.8 s/draw worst case.
# Depth 7+ during adaptation was killing the previous run; cap fixes it.
log_line("[$(round(time();digits=1))] starting NUTS: 50 adapt + 50 sample, max_depth=6")
spl = NUTS(50, 0.65; max_depth=6, adtype=AutoForwardDiff())
t_chain = @elapsed begin
    chain = sample(model, spl, 50; progress=false)
end
log_line("[$(round(time();digits=1))] NUTS done")
println("NUTS 50 samples: ", round(t_chain; digits=1), " s (",
        round(t_chain/50*1000; digits=0), " ms/sample)")
flush(stdout)
try
    println("Mean accept:     ", round(mean(skipmissing(vec(chain[:acceptance_rate]))); digits=3))
catch; end
try
    println("Mean tree depth: ", round(mean(skipmissing(vec(chain[:tree_depth]))); digits=2))
catch; end
println("logp range: ", round(minimum(vec(chain[:logjoint]));digits=2),
        " to ", round(maximum(vec(chain[:logjoint]));digits=2))

df = DataFrame(chain)
println("\nParameter ranges (median, q025, q975):")
for p in [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :t_0]
    string(p) in names(df) || continue
    v = df[!, p]
    println("  $(p): ", round(median(v);sigdigits=4), "  (",
            round(quantile(v,0.025);sigdigits=4), ", ",
            round(quantile(v,0.975);sigdigits=4), ")")
end

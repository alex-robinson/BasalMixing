## Direct ForwardDiff-through-ODE timing diagnostic.
## Isolates the cost of (gradient eval) / (forward eval) for run_basal_mixing_ode.
## No Turing, no NUTS — just the raw AD path that NUTS depends on.

cd(@__DIR__)
import Pkg; Pkg.activate(".")

log_line(s) = (println(s); flush(stdout))

log_line("[startup] loading…")
using CSV, DataFrames, Statistics
using DifferentialEquations
using ForwardDiff
include("BasalMixingModel.jl")

depth, _ = generate_depths("highdirty"; step=0.25)
(k81, dar40) = load_basalmixing_data(depth=depth)
geom = BasalMixingGeometry(depth, 3040.0, 3053.44, k81.depth, dar40.depth)

log_line("[setup] grid N=$(geom.N), k81_obs=$(length(geom.k81_interp)), dar40_obs=$(length(geom.dar40_interp))")

# Test parameter point (near the chain MAP).
θ0 = [1.0, 0.03, 6.0, 250.0, 0.005, 2304.0]

# A scalar functional of the predictions, used as a dummy "loss" so ForwardDiff
# can produce a gradient. NUTS computes the log-likelihood; the cost of taking
# gradient(loss) is what we care about, not the loss value itself.
function loss(θ; reltol, abstol)
    p = (delta=θ[1], m_clean=θ[2], f_dirty=θ[3],
         t_old=θ[4], F_ar40=θ[5], t_target=θ[6])
    k81_pred, dar40_pred, ok = run_basal_mixing_ode(p, geom; reltol=reltol, abstol=abstol)
    return sum(k81_pred) + sum(dar40_pred)
end

function time_at_tol(reltol, abstol; n_warm=2, n_meas=5)
    log_line("\n[tolerance reltol=$reltol abstol=$abstol]")
    # Warm-up
    for _ in 1:n_warm; loss(θ0; reltol=reltol, abstol=abstol); end
    t_fwd = @elapsed for _ in 1:n_meas
        loss(θ0; reltol=reltol, abstol=abstol)
    end
    t_fwd_ms = 1000 * t_fwd / n_meas
    log_line("  forward (n=$n_meas): $(round(t_fwd_ms; digits=2)) ms/eval")

    f = θ -> loss(θ; reltol=reltol, abstol=abstol)
    cfg = ForwardDiff.GradientConfig(f, θ0, ForwardDiff.Chunk{6}())
    for _ in 1:n_warm; ForwardDiff.gradient(f, θ0, cfg); end
    t_grad = @elapsed for _ in 1:n_meas
        ForwardDiff.gradient(f, θ0, cfg)
    end
    t_grad_ms = 1000 * t_grad / n_meas
    ratio = t_grad_ms / t_fwd_ms
    log_line("  gradient (n=$n_meas): $(round(t_grad_ms; digits=2)) ms/eval ($(round(ratio; digits=2))× forward)")

    # Sanity: print the gradient itself
    g = ForwardDiff.gradient(f, θ0, cfg)
    log_line("  ∇loss at θ0: $(round.(g; sigdigits=4))")

    # NUTS budget extrapolation
    for tree_depth in [4, 5, 6, 7, 8]
        n_leap = 2^tree_depth
        per_draw_s = n_leap * t_grad_ms / 1000
        log_line("  NUTS @ tree_depth=$tree_depth ($n_leap leapfrogs/draw): " *
                 "$(round(per_draw_s; digits=2)) s/draw, " *
                 "$(round(per_draw_s*100/60; digits=1)) min for 100 draws")
    end
end

time_at_tol(1e-6, 1e-8)
time_at_tol(1e-4, 1e-7)
time_at_tol(1e-3, 1e-6)

## Validate the new ODE-based forward model against the existing
## hand-rolled Euler integrator. Same inputs → predictions should agree
## within Tsit5 tolerance.

cd(@__DIR__)
import Pkg; Pkg.activate(".")

using JLD2, DataFrames, Turing, CSV
using Statistics
using DifferentialEquations

include("BasalMixingModel.jl")

depth, setup = generate_depths("highdirty"; step=0.25)
(k81, dar40) = load_basalmixing_data(depth=depth)

geom = BasalMixingGeometry(depth, 3040.0, 3053.44, k81.depth, dar40.depth)

# A reasonable parameter point (close to the chain's combined MAP).
p_test = (
    delta    = 1.0,
    m_clean  = 0.03,
    f_dirty  = 6.0,
    t_old    = 250.0,
    F_ar40   = 0.005,
    t_target = 2304.0,   # |t_0| from the combined chain MAP
)

println("Test parameters: ", p_test)

# --- Euler path ---
b = BasalMixingModel(depth=depth, k81_obs_depths=k81.depth, dar40_obs_depths=dar40.depth)
p_euler = (delta=p_test.delta, m_clean=p_test.m_clean, f_dirty=p_test.f_dirty,
           t_old=p_test.t_old, F_ar40=p_test.F_ar40)
ok_euler = RunBasalMixingModelToTime!(p_euler, b, p_test.t_target, (k81, dar40); dt=0.1)
println("Euler success: $ok_euler")
k81_euler   = copy(b.k81.dat[:, 1])
dar40_euler = copy(b.dar40.dat[:, 1])

# --- ODE path ---
k81_ode, dar40_ode, ok = run_basal_mixing_ode(p_test, geom)
println("ODE success: $ok")

# --- Compare ---
println("\nk81 age (kyr):")
println("  Euler: ", round.(k81_euler;   digits=2))
println("  ODE:   ", round.(k81_ode;     digits=2))
println("  Δ:     ", round.(k81_ode .- k81_euler; digits=4))
println("  max|Δ| = ", maximum(abs.(k81_ode .- k81_euler)))

println("\ndar40 (per mil):")
println("  Euler: ", round.(dar40_euler; digits=5))
println("  ODE:   ", round.(dar40_ode;   digits=5))
println("  Δ:     ", round.(dar40_ode .- dar40_euler; digits=6))
println("  max|Δ| = ", maximum(abs.(dar40_ode .- dar40_euler)))

# --- Timing ---
println("\nTiming (warm):")
RunBasalMixingModelToTime!(p_euler, b, p_test.t_target, (k81, dar40); dt=0.2)  # warm
run_basal_mixing_ode(p_test, geom)
t_euler = @elapsed for _ in 1:5
    RunBasalMixingModelToTime!(p_euler, b, p_test.t_target, (k81, dar40); dt=0.2)
end
t_ode = @elapsed for _ in 1:5
    run_basal_mixing_ode(p_test, geom)
end
println("  Euler (5 calls, dt=0.1): ", round(t_euler*1000; digits=2), " ms (", round(t_euler*200; digits=2), " ms/call)")
println("  ODE   (5 calls, Tsit5): ", round(t_ode*1000;   digits=2), " ms (", round(t_ode*200;   digits=2), " ms/call)")

# --- Robustness sweep: 20 random draws from priors ---
println("\nRandom-parameter sweep (20 draws):")
using Random; Random.seed!(42)
import Distributions: Uniform, truncated, Normal
priors_test = (
    delta   = Uniform(0.3, 2.0),
    m_clean = truncated(Normal(0.03, 0.002), lower=0.0),
    f_dirty = Uniform(4.0, 6.5),
    t_old   = truncated(Normal(250.0,25.0), lower=0.0),
    F_ar40  = Uniform(0.004, 0.007),
    t_0     = Uniform(-3000.0, -700.0),
)
global max_dk81 = 0.0
global max_ddar40 = 0.0
global n_both_ok = 0
for trial in 1:20
    p = (
        delta    = rand(priors_test.delta),
        m_clean  = rand(priors_test.m_clean),
        f_dirty  = rand(priors_test.f_dirty),
        t_old    = rand(priors_test.t_old),
        F_ar40   = rand(priors_test.F_ar40),
        t_target = -rand(priors_test.t_0),
    )
    pe = (delta=p.delta, m_clean=p.m_clean, f_dirty=p.f_dirty, t_old=p.t_old, F_ar40=p.F_ar40)
    ok_e = RunBasalMixingModelToTime!(pe, b, p.t_target, (k81, dar40); dt=0.1)
    ke = copy(b.k81.dat[:, 1]);  de = copy(b.dar40.dat[:, 1])
    ko, do_, ok_o = run_basal_mixing_ode(p, geom)
    if ok_e && ok_o
        global n_both_ok += 1
        dk = maximum(abs.(ko .- ke));  dd = maximum(abs.(do_ .- de))
        global max_dk81 = max(max_dk81, dk);  global max_ddar40 = max(max_ddar40, dd)
    end
    println("  trial $(lpad(trial,2)): euler=$ok_e ode=$ok_o  t_target=$(round(p.t_target;digits=1))" *
            (ok_e && ok_o ? "  Δage_max=$(round(maximum(abs.(ko .- ke));digits=3))  Δdar40_max=$(round(maximum(abs.(do_ .- de));digits=5))" : ""))
end
println("\nBoth-ok count: $n_both_ok/20")
println("Max |Δ age|   over both-ok trials: $(round(max_dk81;   digits=3)) kyr")
println("Max |Δ dar40| over both-ok trials: $(round(max_ddar40; digits=5)) ‰")

## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise

using Turing
using LinearAlgebra
using Random

using Dates
using CairoMakie
using CSV
using DataFrames

Random.seed!(42)

# Include code for basal mixing model
include("BasalMixingModel.jl")

function plot_prior_line!(ax, prior::Distribution; color=:red, kwargs...)
    # Works for any Distributions.jl type
    lo, hi = quantile(prior, 0.001), quantile(prior, 0.999)
    x = range(lo, hi, length=300)
    lines!(ax, x, pdf.(prior, x); color=color, kwargs...)
end

# priors = (
#     L_ref       = Uniform(0.1, 1.5),    # 1 m
#     delta       = Uniform(0.1, 2.0),    # 1 m
#     m_clean     = truncated(Normal(0.03, 0.02), lower=0.0), # 0.03 m/kyr
#     m_dirty     = truncated(Normal(0.18, 0.10), lower=0.0), # 0.18 m/kyr  
#     σ           = Exponential(30.0),    # 30 kyr
# )
priors = (
    L_ref       = Uniform(0.1, 1.5),    # 1 m
    delta       = Uniform(0.1, 2.0),    # 1 m
    m_clean     = Uniform(0.0, 0.1),    # 0.03 m/kyr
    m_dirty     = Uniform(0.0,0.5),     # 0.18 m/kyr  
    σ           = Exponential(30.0),    # 30 kyr
)

@model function basal_mixing(age_obs, depth, dt, priors)
    L_ref       ~ priors.L_ref      
    #L_ref       = 1.0
    delta ~ priors.delta     
    #delta = 1.0
    m_clean     ~ priors.m_clean
    m_dirty     ~ priors.m_dirty 
    σ ~ priors.σ
    #σ = 30.0
    
    # Run the model
    p = (L_ref=L_ref, delta=delta, m_clean=m_clean, m_dirty=m_dirty)
    _, _, b2, success = RunBasalMixingModel(p; depth=depth, dt=dt, sampling=true)
    
    if !success
        # Try one more time with a smaller timestep
        _, _, b2, success = RunBasalMixingModel(p; depth=depth, dt=dt*0.5, sampling=true)
    end

    # Extract the best-fit age profile at the optimal time

    n_obs = length(age_obs)

    if success
        # Interpolate each depth's age time series at t_best
        #age_pred = [linterp(b2.time, b2.age_k81[i, :], t_best) for i in 1:n_obs]

        # Compute SSR at every time step across all observations
        n_times = length(b2.time)
        ssr = [sum((age_obs[i] - b2.age_k81[i, t])^2 for i in 1:n_obs) for t in 1:n_times]

        # Find the time index that minimises SSR
        t_best_idx = argmin(ssr)
        time_pred  = b2.time[t_best_idx]
        age_pred   = b2.age_k81[:, t_best_idx]
    else
        # Assign a very high age so that Likelihood is low
        time_pred = 0.0
        age_pred  = fill(1e8, n_obs)
    end

    # Likelihood: observed k81 ages ~ Normal(predicted, σ)
    age_obs ~ MvNormal(age_pred, σ * I)

    return (time_pred = time_pred, age_pred = age_pred)  # return whatever you want
end

## SCRIPT ##

# Load datasets for comparison
(k81, ar40) = load_basalmixing_data()
age_k81 = k81[!,"age"]

depth, setup = generate_depths("highdirty";step=0.25)

model = basal_mixing(age_k81, depth, 0.2, priors)

# Start with MH to verify it works, then switch to SMC for better exploration
chain = sample(model, MH(), MCMCThreads(), 5_000, 4)  # 4 chains in parallel

# Or SMC (no chain length needed — uses particle count)
#chain = @timed sample(model, SMC(1000), 1000)

## Analysis

begin
    df = DataFrame(chain)
    params = [:L_ref, :delta, :m_clean, :m_dirty, :t_best]
    labels = ["L_ref", "delta", "m_clean", "m_dirty", "t_best"]
    best_idx = argmax(df.logjoint)

    describe(chain)
end


begin
    p = (
        L_ref       = df.L_ref[best_idx],
        depth_scale = df.delta[best_idx],
        m_clean     = df.m_clean[best_idx],
        m_dirty     = df.m_dirty[best_idx],
    )

    b_cache = BasalMixingModel(depth=collect(depth))

    @btime b, b1, b2, rmse_k81 = RunBasalMixingModel(p, (k81, ar40); b=b_cache, depth=depth, dt=0.1)

    (time, rmse) = rmse_k81
    println("k81 (time, rmse): $time, $rmse")

    # Plot the results
    fig = plot_BasalMixingModelRun(b,b1,b2;k81=k81) #,ar40=ar40_data)
end

begin
    fig = Figure(size=(900, 700))

    ipar = 0
    for (i, (param, label)) in enumerate(zip(params, labels))
        if string(param) in names(df)
            ipar += 1
            row, col = divrem(ipar-1, 2)
            ax = Axis(fig[row+1, col+1], xlabel=label, ylabel="log p")
            ylims!(ax,(-50,0))
            scatter!(ax, df[!, param], df.logjoint, alpha=0.6, markersize=6, color=:steelblue)
            # Highlight best point in red
            scatter!(ax, [df[best_idx, param]], [df.logjoint[best_idx]], 
                    color=:red, markersize=12, label="MAP")
        end
    end

    display(fig)
end

begin
    fig = Figure(size=(900, 700))

    ipar = 0
    for (i, (param, label)) in enumerate(zip(params, labels))
        if string(param) in names(df)
            ipar += 1
            row, col = divrem(ipar-1, 2)
            ax = Axis(fig[row+1, col+1], xlabel=label, ylabel="density")
            hist!(ax, df[!, param], bins=20, normalization=:pdf, color=(:steelblue, 0.7))
            plot_prior_line!(ax, priors[param], label="Prior", linewidth=2, color=:grey50)
            vlines!(ax, [df[best_idx, param]], color=:red, linewidth=2, label="MAP")
        end
    end

    display(fig)
end

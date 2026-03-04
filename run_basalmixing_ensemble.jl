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

priors = (
    delta       = Uniform(0.1, 2.0),    # 1 m
    m_clean     = truncated(Normal(0.03, 0.02), lower=0.0), # 0.03 m/kyr
    m_dirty     = truncated(Normal(0.18, 0.10), lower=0.0), # 0.18 m/kyr  
    #σ           = Exponential(30.0),    # 30 kyr
    t_old = 250.0,
    σ = 30.0
)
# priors = (
#     delta       = Uniform(0.1, 2.0),    # 1 m
#     m_clean     = Uniform(0.0, 0.1),    # 0.03 m/kyr
#     m_dirty     = Uniform(0.0,0.5),     # 0.18 m/kyr  
#     #t_old       = Uniform(150.0,350.0), # 250 kyr
#     #σ           = Exponential(30.0),    # 30 kyr
#     t_old = 250.0,
#     σ = 30.0
# )

@model function basal_mixing(age_obs, b, dat, dt, priors)

    ## Set priors ##

    delta ~ priors.delta
    m_clean ~ priors.m_clean
    m_dirty ~ priors.m_dirty
    t_old = priors.t_old
    σ     = priors.σ

    ## Extract obs ##
    (k81, ar40) = dat
    #age_obs = k81.age # Must be passed in separately as a vector to be sampled

    # Run the model
    p = (delta=delta, m_clean=m_clean, m_dirty=m_dirty, t_old=t_old)
    success = RunBasalMixingModel!(p, b, dat; dt=dt, sampling=true)
    
    if !success
        # Try one more time with a smaller timestep
        success = RunBasalMixingModel!(p, b, dat; dt=dt*0.5, sampling=true)
    end

    if success
        # Extract the best-fit ages at the optimal time
        time_pred = b.b2.time[b.b2.kmin]
        age_pred  = b.b2.age_k81[:,b.b2.kmin]
    else
        # Assign a very high age so that Likelihood is very low
        time_pred = 0.0
        age_pred  = fill(1e8, length(age_obs))
    end

    # Likelihood: observed k81 ages ~ Normal(predicted, σ)
    age_obs ~ MvNormal(age_pred, σ * I)

    return (time_pred = time_pred, age_pred = age_pred)
end

## SCRIPT ##

# Load datasets for comparison
(k81, ar40) = load_basalmixing_data()

depth, setup = generate_depths("highdirty";step=0.25)
b = BasalMixingModel(depth=depth)

model = basal_mixing(k81.age, b, (k81, ar40), 0.2, priors)

# Start with MH to verify it works, then switch to SMC for better exploration
chain = sample(model, MH(), MCMCThreads(), 1000, 4)  # 4 chains in parallel

# Or SMC (no chain length needed — uses particle count)
#chain = @timed sample(model, SMC(1000), 1000)

## Analysis

begin
    df = DataFrame(chain)
    params = [:delta, :m_clean, :m_dirty, :t_old]
    labels = ["delta", "m_clean", "m_dirty", "t_old"]
    best_idx = argmax(df.logjoint)

    describe(chain)
end


begin
    p = (
        delta   = df.delta[best_idx],
        m_clean = df.m_clean[best_idx],
        m_dirty = df.m_dirty[best_idx],
        t_old   = 250.0,
    )

    b = BasalMixingModel(depth=depth)

    RunBasalMixingModel!(p, b, (k81, ar40); dt=0.1, sampling=false)

    # Plot the results
    fig = plot_BasalMixingModelRun(b; k81=k81) #,ar40=ar40_data)
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

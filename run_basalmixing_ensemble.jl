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
    f_dirty     = Uniform(4.0, 8.0),                        # 0.18 m/kyr / 0.03 m/kyr => f=6x
    t_old       = Normal(250.0,100.0),  # 250 kyr
    F_ar40      = Uniform(0.0,0.1), #Normal(0.075,0.01),   # 0.075 m^3 / kyr
    #σ_k81           = Exponential(30.0),    # 30 kyr
    #t_old = 250.0,
    σ_k81 = 30.0,
    σ_dar40 = 0.03,

    time_pred   = Uniform(0.0,3000.0),  # Anything is possible
)

@model function basal_mixing(k81_age_obs, dar40_obs, b, dat, dt, priors)

    ## Set priors ##
    delta   = priors.delta   isa Distribution ? delta   ~ priors.delta   : priors.delta
    m_clean = priors.m_clean isa Distribution ? m_clean ~ priors.m_clean : priors.m_clean
    f_dirty = priors.f_dirty isa Distribution ? f_dirty ~ priors.f_dirty : priors.f_dirty
    t_old   = priors.t_old   isa Distribution ? t_old   ~ priors.t_old   : priors.t_old
    F_ar40  = priors.F_ar40  isa Distribution ? F_ar40  ~ priors.F_ar40  : priors.F_ar40
    σ_k81   = priors.σ_k81   isa Distribution ? σ_k81   ~ priors.σ_k81   : priors.σ_k81
    σ_dar40 = priors.σ_dar40 isa Distribution ? σ_dar40 ~ priors.σ_dar40 : priors.σ_dar40

    ## Extract obs ##
    (k81, ar40) = dat
    #k81_age_obs = k81.age              # Must be passed in separately as a vector to be sampled
    #dar40_obs = ar40.dar40             # Must be passed in separately as a vector to be sampled

    # Run the model
    p = (delta=delta, m_clean=m_clean, f_dirty=f_dirty, t_old=t_old, F_ar40=F_ar40)
    success = RunBasalMixingModel!(p, b, dat; dt=dt, sampling=true)
    
    if !success
        # Try one more time with a smaller timestep
        success = RunBasalMixingModel!(p, b, dat; dt=dt*0.5, sampling=true)
    end

    if success
        # Extract the best-fit ages at the optimal time
        kmin = b.kmin
        time_pred := b.b2.time[kmin]
        k81_age_pred  := b.b2.age_k81[:,kmin]
        dar40_pred := b.b3.dar40[:,kmin]
    else
        # Assign a very high age so that Likelihood is very low
        time_pred := 0.0
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

## SCRIPT ##

depth, setup = generate_depths("highdirty";step=0.25)
b = BasalMixingModel(depth=depth)
(k81, ar40) = load_basalmixing_data(depth=b.depth)

model = basal_mixing(k81.age, ar40.dar40, b, (k81, ar40), 0.2, priors)

# Sample using Metropolis Hastings (MH)
chain = sample(model, MH(), MCMCThreads(), 1_000, 4)  # 4 chains in parallel

## Analysis

begin
    df = DataFrame(chain)
    params = [:delta, :m_clean, :f_dirty, :t_old, :F_ar40, :time_pred]
    labels = ["delta (m)", "m_clean (m/yr)", "f_dirty", "t_old (kyr)", "F_ar40 (m^3/kyr)", "time_pred (kyr)"]
    best_idx = argmax(df.logjoint)

    describe(chain)
end


begin
    p = (
        delta   = df.delta[best_idx],
        m_clean = df.m_clean[best_idx],
        f_dirty = df.f_dirty[best_idx],
        t_old   = df.t_old[best_idx],
        F_ar40  = df.F_ar40[best_idx],
    )

    b = BasalMixingModel(depth=depth)

    RunBasalMixingModel!(p, b, (k81, ar40); dt=0.1, sampling=false)

    # Plot the results
    fig = plot_BasalMixingModelRun(b; k81=k81,ar40=ar40)
    display(fig)
    mysave(plt_prefix()*"mixingmodel-ens-best.png",fig)
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
    mysave(plt_prefix()*"mixingmodel-ens-logp.png",fig)
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
    mysave(plt_prefix()*"mixingmodel-ens-hist.png",fig)
end

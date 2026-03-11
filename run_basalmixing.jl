## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise
using Dates
using CairoMakie
using CSV
using DataFrames

# Include code for basal mixing model
include("BasalMixingModel.jl")

begin
    
    #depth, setup = generate_depths("default")
    #depth, setup = generate_depths("high")
    depth, setup = generate_depths("highdirty";step=0.25)

    (k81_obs, dar40_obs) = load_basalmixing_data(depth=depth) # Load datasets for comparison

    p = (
        delta = 1.0,
        m_clean = 0.03,
        f_dirty = 6.0,
        t_old = 250.0,
        F_ar40 = 0.0075,
    )

    b = BasalMixingModel(depth=depth,k81_obs_depths=k81_obs.depth,dar40_obs_depths=dar40_obs.depth)
    
    RunBasalMixingModel!(p, b, (k81_obs, dar40_obs); dt=0.1, sampling=false)
    @btime RunBasalMixingModel!(p, b, (k81_obs, dar40_obs); dt=0.1, sampling=true)

    # Plot the results
    fig = plot_BasalMixingModelRun(b; k81_obs=k81_obs, dar40_obs=dar40_obs)

end

mysave(plt_prefix()*"mixingmodel_$setup.png",fig)
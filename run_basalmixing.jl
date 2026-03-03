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

    p = (
        L_ref = 1.0,
        depth_scale = 1.0,
        m_clean = 0.03,
        m_dirty = 0.18
    )

    b, b1, b2 = RunBasalMixingModel(p; depth=depth,t1=3000.0,dt=0.1)

    # Plot the results
    fig = plot_BasalMixingModelRun(b,b1,b2;k81=k81_data) #,ar40=ar40_data)
end

mysave(plt_prefix()*"mixingmodel_$setup.png",fig)
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

# Load datasets for comparison
(k81, ar40) = load_basalmixing_data()

begin
    #depth, setup = generate_depths("default")
    #depth, setup = generate_depths("high")
    depth, setup = generate_depths("highdirty";step=0.25)

    p = (
        delta = 1.0,
        m_clean = 0.03,
        m_dirty = 0.18,
        t_old = 250.0,
    )

    b = BasalMixingModel(depth=depth)

    RunBasalMixingModel!(p, b, (k81, ar40); dt=0.1, sampling=false)
    #@btime RunBasalMixingModel!(p, b, (k81, ar40); dt=0.1, sampling=true)

    # Plot the results
    fig = plot_BasalMixingModelRun(b; k81=k81) #,ar40=ar40_data)

end

mysave(plt_prefix()*"mixingmodel_$setup.png",fig)
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
        L_ref = 1.0,
        delta = 1.0,
        m_clean = 0.03,
        m_dirty = 0.18
    )

    b, b1, b2, rmse_k81 = RunBasalMixingModel(p, (k81, ar40); depth=depth, dt=0.1)

    (time, rmse) = rmse_k81
    println("k81 (time, rmse): $time, $rmse")

    # Plot the results
    fig = plot_BasalMixingModelRun(b,b1,b2;k81=k81) #,ar40=ar40_data)

end

mysave(plt_prefix()*"mixingmodel_$setup.png",fig)
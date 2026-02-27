## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise
using Dates
using CairoMakie
using CSV
using DataFrames


# Get data to compare with
ar40_data = CSV.read("data/Bender2010_ar40_data.txt",DataFrame;delim="|",ignorerepeated=true)
rename!(ar40_data, strip.(names(ar40_data)))

k81_data = CSV.read("data/k81_data.txt",DataFrame;delim=" ",ignorerepeated=true)
rename!(k81_data, strip.(names(k81_data)))
k81_data[!,:depth] = 0.5 .* (k81_data[!,"depth_top"] .+ k81_data[!,"depth_bottom"])

# Include code for basal mixing model
include("BasalMixingModel.jl")

begin
    #depth, setup = generate_depths("default")
    #depth, setup = generate_depths("high")
    depth, setup = generate_depths("highdirty";step=0.25)

    m_clean = 0.03
    m_dirty = m_clean*6
    depth_scale = 1.0

    #f_mixing_rate = make_mixing_rate_discrete(3040.0, m_clean, m_dirty, depth_scale)
    f_mixing_rate = make_mixing_rate_smooth(3040.0, m_clean, m_dirty, depth_scale)
    #f_mixing_rate = make_mixing_rate_continuous(3040.0, m_clean, m_dirty, depth_scale)
    #f_mixing_rate = make_mixing_rate_exponential(3040.0, m_clean, m_dirty, depth_scale)

    b, b1, b2 = RunBasalMixingModel(;depth=depth,f_mixing_rate=f_mixing_rate,t1=3000.0,dt=0.05)

    # Plot the results
    fig = plot_BasalMixingModelRun(b,b1,b2;k81=k81_data) #,ar40=ar40_data)
end

mysave(plt_prefix()*"mixingmodel_$setup.png",fig)
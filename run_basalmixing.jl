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

# Get data to compare with
ar40_data = CSV.read("data/Bender2010_ar40_data.txt",DataFrame;delim="|",ignorerepeated=true)
rename!(ar40_data, strip.(names(ar40_data)))

k81_data = CSV.read("data/k81_data.txt",DataFrame;delim=" ",ignorerepeated=true)
rename!(k81_data, strip.(names(k81_data)))
k81_data[!,:depth] = 0.5 .* (k81_data[!,"depth_top"] .+ k81_data[!,"depth_bottom"])

# Run the basal mixing model
# Return current state, summary1 and summary2

# Default depths
depth = 3035:1.0:3053
setup = "default"

# High resolution depths
depth = 3035:0.1:3053
setup = "high"

# Variable resolution depths
depths_clean = collect(3035:3040)
depths_dirty = range(3040.0,3053.0; step=0.2) #length=14)
depth = unique([depths_clean...,depths_dirty...])
setup = "high"

# Special depths
depth = collect(3035:1.0:3053) #.+ 0.5
setup = "special"


b, b1, b2 = RunBasalMixingModel(;depth=depth,t1=3000.0,dt=0.1)

# Plot the results

fig = plot_BasalMixingModelRun(b,b1,b2;k81=k81_data) #,ar40=ar40_data)

mysave(plt_prefix()*"mixingmodel_$setup.png",fig)
## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise
using CairoMakie
using CSV
using DataFrames

# Include code for basal mixing model
include("BasalMixingModel.jl")

# Test calculating concentration over time (half-life decay)
t = 0.0:3000.0
c = concentration.(1.0, t)

# Get data to compare with
ar40_data = CSV.read("data/Bender2010_ar40_data.txt",DataFrame;delim="|",ignorerepeated=true)
rename!(ar40_data, strip.(names(ar40_data)))

k81_data = CSV.read("data/k81_data.txt",DataFrame;delim=" ",ignorerepeated=true)
rename!(k81_data, strip.(names(k81_data)))
k81_data[!,:depth] = 0.5 .* (k81_data[!,"depth_top"] .+ k81_data[!,"depth_bottom"])

# Run the basal mixing model
# Return current state, summary1 and summary2

b, b1, b2 = RunBasalMixingModel(;t1=3000.0,dt=0.5)

# Plot the results

fig = plot_BasalMixingModelRun(b,b1,b2;k81=k81_data,ar40=ar40_data)

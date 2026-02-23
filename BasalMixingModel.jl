using CairoMakie

mutable struct BasalMixingModel
    n::Int
    layer::Vector{Int}
    depth::Vector{Float64}
    thickness::Vector{Float64}
    age_k81::Vector{Float64}
    c_k81::Vector{Float64}
    c_ar40::Vector{Float64}

    time::Float64
end

function BasalMixingModel(;
    depth::Union{Vector{Float64},UnitRange{Int64}} = collect(3040.0:3053.0),
    age_k81 = zeros(length(depth)),
    c_k81 = ones(length(depth)),
    c_ar40 = ones(length(depth)),
    time = 0.0
    )

    n = length(depth)
    layers = 1:n
    thickness = fill(diff(depth)[1],n)  # Uniform layer thickness
    thickness = [diff(depth)..., diff(depth)[end]]
    
    return BasalMixingModel(
        n,
        collect(layers),
        collect(depth),
        collect(thickness),
        age_k81,
        c_k81,
        c_ar40,
        time
    )
end

struct BasalMixingModelSummary1
    times::Vector{Float64}      # times of snapshots
    depth::Vector{Float64}      # depth axis
    age_k81::Array{Float64}     # [nt,nd]
    c_k81::Array{Float64}       # [nt,nd]
    c_ar40::Array{Float64}      # [nt,nd]
end

function BasalMixingModelSummary1(times::Vector{Float64},depth::Vector{Float64})
    nt = length(times)
    nd = length(depth)

    age_k81 = fill(0.0, nt, nd)
    c_k81 = fill(0.0, nt, nd)
    c_ar40 = fill(0.0, nt, nd)

    return BasalMixingModelSummary1(
        times, depth, age_k81, c_k81, c_ar40
    )
end

struct BasalMixingModelSummary2
    depths::Vector{Float64}     # depths of interest
    time::Vector{Float64}       # time axis
    age_k81::Array{Float64}     # [nd,nt]
    c_k81::Array{Float64}       # [nd,nt]
    c_ar40::Array{Float64}      # [nd,nt]
end

function BasalMixingModelSummary2(depths::Vector{Float64},time::Vector{Float64})
    nd = length(depths)
    nt = length(time)
    
    age_k81 = fill(0.0, nd, nt)
    c_k81 = fill(0.0, nd, nt)
    c_ar40 = fill(0.0, nd, nt)

    return BasalMixingModelSummary2(
        depths, time, age_k81, c_k81, c_ar40
    )
end

function concentration(R0::Float64, t::Float64; t_half::Float64 = 229.0)
    return R0 * 2.0^(-t / t_half)
end

function concentration_to_age(c::Float64, c0::Float64 = 1.0; t_half::Float64 = 229.0)
    return -t_half * log2(c / c0)
end

function decay_step(R::Float64, dt::Float64; t_half::Float64 = 229.0)
    return R * 2.0^(-dt / t_half)
end

function mixing_step(R0::Float64, R1::Float64, dz0::Float64, dz1::Float64, dt::Float64, mixing_rate::Float64)
    dV = mixing_rate * dt

    @assert dV < dz0 "Mixing volume ($dV) exceeds cell 0 thickness ($dz0) — reduce dt or mixing_rate"
    @assert dV < dz1 "Mixing volume ($dV) exceeds cell 1 thickness ($dz1) — reduce dt or mixing_rate"

    R0_new = (R0 * dz0 - dV * R0 + dV * R1) / dz0
    R1_new = (R1 * dz1 - dV * R1 + dV * R0) / dz1

    #@assert R0_new*dz0 + R1_new*dz1 == R0*dz0 + R1*dz1 "Concentration is not conserved!"
    return R0_new, R1_new
end

"""
    step_ar40(cc_ar40, flux, dt)

Step the total ⁴⁰Ar content forward in time.

Arguments:
- `cc_ar40`: current total ⁴⁰Ar content [cc m⁻²]
- `flux`: basal ⁴⁰Ar flux [cc m⁻² kyr⁻¹]
- `dt`: timestep [kyr]

Returns updated ⁴⁰Ar content [cc m⁻²]
"""
function step_ar40(cc_ar40::Float64, flux::Float64, dt::Float64)
    return cc_ar40 + flux * dt
end


"""
    cc_to_delta_ar40(ar40, ar40_ref)

Convert total ⁴⁰Ar content [cc m⁻²] to δ⁴⁰ArATM [‰].

The reference is the atmospheric ⁴⁰Ar content of the ice layer, assuming:
- 8% total air content by volume
- atmospheric ⁴⁰Ar fraction of 0.00934 (9340 ppm)

Arguments:
- `ar40`: total ⁴⁰Ar content [cc m⁻²]
- `ar40_ref`: total ⁴⁰Ar content [cc m⁻²] at present day

Returns δ⁴⁰ArATM [‰]
"""
function cc_to_delta_ar40(ar40::Float64, ar40_ref::Float64)
    return (ar40 / ar40_ref - 1.0) * 1000.0
end

# Atmospheric ⁴⁰Ar volume fraction as a function of age:
f_ar40_atm(t_kyr) = 9340e-6 * (1 - 0.066e-6 * t_kyr)  # returns volume fraction

function calc_ar40_ref(thickness::Float64;
    t_kyr = 0.0,                        # kyr, for which age do we want the reference?
    air_content::Float64 = 0.08,        # fraction
    f_ar40_atm::Float64 = 0.00934       # atmospheric ⁴⁰Ar volume fraction
    )

    # Reference ⁴⁰Ar content of the layer [cc m⁻²]
    # thickness [m] × 1e6 [cc m⁻³ per m] × air_content × f_ar40_atm
    ar40_ref = thickness * 1e6 * air_content * f_ar40_atm

    return ar40_ref
end

function RunBasalMixingModel(;t0=0.0,t1=1000.0,dt=1.0,mixing_rate_clean=0.03,mixing_rate_bottom=0.03*6,t_old=250.0)

    b = BasalMixingModel()

    # Get times to model
    time = t0:dt:t1

    # Get times and depths of interest
    times = collect(500.0:100.0:3000.0)
    depths = [3045.0, 3047.0, 3050.0]

    # Define summary objects
    b1 = BasalMixingModelSummary1(times,b.depth)
    b2 = BasalMixingModelSummary2(depths,collect(time))

    mixing_rate = fill(mixing_rate_bottom,b.n-1)
    mixing_rate[1] = mixing_rate_clean

    # Set initial values
    b.c_k81 .= 1.0
    #b.c_ar40 = 

    # Loop over time and advance model
    for (k, t) in enumerate(time)

        # Save previous clean ice concentration
        c_k81_clean = b.c_k81[1]

        # Update concentration from aging
        b.c_k81 .= decay_step.(b.c_k81,dt; t_half=229.0)

        # Update concentration from mixing
        for j in 1:b.n-1
            tmp0, tmp1 = mixing_step(
                b.c_k81[j], b.c_k81[j+1], 
                b.thickness[j], b.thickness[j+1], dt, 
                mixing_rate[j]
            )
            b.c_k81[j], b.c_k81[j+1] = tmp0, tmp1
        end

        # Restore clean ice age beyond t_old time
        if t > t_old
            b.c_k81[1] = c_k81_clean
        end

        # Get ages too
        b.age_k81 .= concentration_to_age.(b.c_k81,1.0)

        # Update to current time
        b.time = t

        ## Update summary objects

        # Store time slice of current variables
        if t in b1.times
            i = findall(t .== b1.times)[1]
            b1.age_k81[i,:] = b.age_k81
            b1.c_k81[i,:] = b.c_k81
        end

        # For each depth d of interest, store the value of the variables
        # at the current time
        for (i, d) in enumerate(b2.depths)
            j = findall( abs.(d .- b.depth) .< 0.4 )[1] # Index of current depth
            b2.age_k81[i,k] = b.age_k81[j]
            b2.c_k81[i,k] = b.c_k81[j]
        end

        
    end

    return b, b1, b2
end

function plot_BasalMixingModelRun(b,b1,b2;k81=nothing,ar40=nothing)

    col_data = "#BC401E"

    if !isnothing(ar40)
        fig = Figure(size=(900,600))
    else
        fig = Figure(size=(600,600))
    end

    ## PANEL 1: Depth versus closed-system age
    ax1 = Axis(fig[1,1], limits=((200,800),(-3053,-3035)), xlabel="⁸¹K closed system age (kyr)", ylabel="Depth (m)" )
    d = collect(-3052:2:-3036)
    ax1.yticks = (d,string.(abs.(d)))
    ax1.xticks = [200,400,600,800]

    # Plot time slices from model
    for (k, t) in enumerate(b1.times)
        lines!(ax1,b1.age_k81[k,:],-b1.depth,color=:grey50,linewidth=0.5)
        if t in [500.0,1000.0,1500.0]
            lines!(ax1,b1.age_k81[k,:],-b1.depth,color=:grey50,linewidth=1.5)
        end
    end

    # Plot data too
    errorbars!(ax1, k81.age, -k81.depth, k81.age_hi, k81.age_lo, color=col_data, direction=:x, whiskerwidth=8)
    scatter!(ax1, k81.age, -k81.depth, color=col_data, marker=:circle, markersize=12)
    
    ## PANEL 2 (optional): depth vs d40Ar_ATM concentration
    if !isnothing(ar40)
        ax3 = Axis(fig[1,end+1], limits=((-0.1,0.62),(-3053,-3035)), xlabel="δ⁴⁰ArATM (‰)", ylabel="Depth (m)" )
        d = collect(-3052:2:-3036)
        ax3.yticks = (d,string.(abs.(d)))
        ax3.xticks = 0.0:0.2:0.6

        # Plot data too
        errorbars!(ax3, ar40[!,"δ40/38atm"],-ar40[!,"depth"], ar40[!,"δ40/38atm_err"], ar40[!,"δ40/38atm_err"], color=col_data, direction=:x, whiskerwidth=8)
        scatter!(ax3, ar40[!,"δ40/38atm"],-ar40[!,"depth"], color=col_data, marker=:circle, markersize=12)
    
    end

    ## PANEL 2 or 3: Closed-system age versus time
    ax2 = Axis(fig[1,end+1], limits=((0,3000),(0,900)), xlabel="Time (kyr)", ylabel="⁸¹K closed system age (kyr)" )
    ax2.xticks = [0,1000,2000,3000]
    ax2.yticks = 0:100:900

    col = ["green","purple","teal"]
    for (i,d) in enumerate(b2.depths)
        lines!(ax2,b2.time,b2.age_k81[i,:],color=col[i],linewidth=2)
    end

    return fig
end

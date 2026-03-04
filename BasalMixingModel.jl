using CairoMakie
using DifferentialEquations
using LinearAlgebra

function plt_prefix(;path="plots")
    return joinpath(path,string(Dates.today())*"_")
end

function mysave(fout,fig;px_per_unit=2)
    println("Saving ",fout)
    save(fout,fig,px_per_unit=px_per_unit)
    return fout
end

function load_basalmixing_data()

    # Get data to compare with
    ar40 = CSV.read("data/Bender2010_ar40_data.txt",DataFrame;delim="|",ignorerepeated=true)
    rename!(ar40, strip.(names(ar40)))

    k81 = CSV.read("data/k81_data.txt",DataFrame;delim=" ",ignorerepeated=true)
    rename!(k81, strip.(names(k81)))
    k81[!,:depth] = 0.5 .* (k81[!,"depth_top"] .+ k81[!,"depth_bottom"])

    #return Dict(:k81=>k81, :ar40=>ar40)
    return (k81, ar40)
end

# Default depths
function generate_depths(setup="default";depth=nothing,step=0.2)
    if setup == "default"
        depth = 3035:1.0:3053
    elseif setup == "high"
        depth = 3035:0.1:3053
    elseif setup == "highdirty"
        depths_clean = collect(3035:3040)
        depths_dirty = range(3039.0,3053.0; step=step) #length=14)
        depth = unique(sort([depths_clean...,depths_dirty...]))
        dx = string(step)
        setup = setup*"-dx$dx"
    elseif setup == "highzoom"
        depths_clean = collect(3035.0:3053.0)
        depths_dirty = range(3040.0,3045.0; step=step) #length=14)
        depth = unique(sort([depths_clean...,depths_dirty...]))
        dx = string(step)
        setup = setup*"-dx$dx"
    else
        @assert !isnothing(depth)
        depth = depth
        setup = setup
    end

    return depth, setup
end

function cell_thickness(depth::AbstractVector{Float64}, depth_bedrock::Float64)

    # Midpoints between adjacent centers form the interior interfaces
    interfaces = (depth[1:end-1] .+ depth[2:end]) ./ 2

    # Top interface: mirror the first interior interface about the first cell center
    top = depth[1] - (interfaces[1] - depth[1])

    # Bottom interface is exactly depth_bedrock
    all_interfaces = [top; interfaces; depth_bedrock]

    return diff(all_interfaces)
end

struct BasalMixingModelSummary1
    times::Vector{Float64}      # times of snapshots
    times_set::Set{Float64}     # Set of times for faster indexing
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
        times, Set(times), depth, age_k81, c_k81, c_ar40
    )
end

mutable struct BasalMixingModelSummary2
    depths::Vector{Float64}     # depths of interest
    time::Vector{Float64}       # time axis
    time_set::Set{Float64}      # time axis set
    age_k81::Array{Float64}     # [nd,nt]
    c_k81::Array{Float64}       # [nd,nt]
    c_ar40::Array{Float64}      # [nd,nt]
    sse_k81::Vector{Float64}    # [nt]
    interp_idx::Vector{Int}     # [nd]
    kmin::Int                   # Index of minimum error
end

function BasalMixingModelSummary2(depths::Vector{Float64},time::Vector{Float64},depth::Vector{Float64})
    nd = length(depths)
    nt = length(time)
    
    age_k81 = fill(0.0, nd, nt)
    c_k81 = fill(0.0, nd, nt)
    c_ar40 = fill(0.0, nd, nt)
    sse_k81 = fill(0.0,nt)

    # Get depth interpolation indices to able to fill b2 with the right depth values
    interp_idx = [findlast(depth .<= d) for d in depths]

    return BasalMixingModelSummary2(
        depths, time, Set(time), age_k81, c_k81, c_ar40, sse_k81, interp_idx, 1
    )
end

mutable struct BasalMixingModel
    n::Int
    depth_lim::Float64
    depth_bedrock::Float64
    jj_clean::Vector{Int}
    layer::Vector{Int}
    depth::Vector{Float64}
    thickness::Vector{Float64}
    mixing_rate::Vector{Float64}
    age_k81::Vector{Float64}
    c_k81::Vector{Float64}
    c_ar40::Vector{Float64}

    time::Float64

    dRdt_mixing::Vector{Float64}
    dRdt_decay::Vector{Float64}

    b1::BasalMixingModelSummary1
    b2::BasalMixingModelSummary2

    rmse_k81::Float64
    time_k81::Float64
end

function BasalMixingModel(;
    depth::Union{Vector{Float64},UnitRange{Int64}} = collect(3035.0:3053.0),
    depth_lim = 3040.00,
    depth_bedrock = 3053.44,
    thickness = cell_thickness(depth,depth_bedrock),
    f_mixing_rate = nothing,
    m_clean = 0.03,
    m_dirty = m_clean*6,
    age_k81 = zeros(length(depth)),
    c_k81 = ones(length(depth)),
    c_ar40 = ones(length(depth)),
    time = 0.0
    )

    n = length(depth)
    layers = 1:n

    # Determine clean ice indices
    jj_clean = findall(depth .<= depth_lim)

    if isnothing(f_mixing_rate)
        mixing_rate = fill(0.0,n)
    else
        mixing_rate = f_mixing_rate(depth)
    end

    dRdt_mixing = zeros(length(depth))
    dRdt_decay  = zeros(length(depth))

    # Define summary objects
    b1 = BasalMixingModelSummary1(collect(500.0:100.0:3000.0), depth)
    b2 = BasalMixingModelSummary2([3044.8, 3047.4, 3049.84],collect(0.0:1.0:3000.0), depth)

    return BasalMixingModel(
        n,
        depth_lim,
        depth_bedrock,
        jj_clean,
        collect(layers),
        collect(depth),
        collect(thickness),
        collect(mixing_rate),
        age_k81,
        c_k81,
        c_ar40,
        time,
        dRdt_mixing,
        dRdt_decay,
        b1,
        b2,
        1e8,
        1e8
    )
end

function ResetBasalMixingModel!(b)

    b.age_k81 .= 0.0
    b.c_k81 .= 0.0
    b.c_ar40 .= 0.0
    b.dRdt_mixing .= 0.0
    b.dRdt_decay .= 0.0
    
    b.b1.age_k81 .= 0.0
    b.b1.c_k81 .= 0.0
    b.b1.c_ar40 .= 0.0
    
    b.b2.age_k81 .= 0.0
    b.b2.c_k81 .= 0.0
    b.b2.c_ar40 .= 0.0
    b.b2.sse_k81 .= 0.0

    return
end

function RunBasalMixingModel!(p, b, dat; t0=0.0,t1=3000.0,dt=1.0,sampling=false)

    # Extract model parameters
    (delta, m_clean, m_dirty, t_old) = p
    L_ref = 1.0     # [m] Use L_ref=1.0, since this just scales m_clean, can tune m_clean directly

    # Extract data for comparison
    (k81, ar40) = dat
    n_obs_k81 = length(k81.age)
    
    # Set the mixing rate
    mixing_rate_smooth!(b.mixing_rate, b.depth, b.depth_lim, m_clean, m_dirty, delta)

    # Get times to model
    time = t0:dt:t1

    ## Set initial values ##

    ResetBasalMixingModel!(b)

    b.c_k81 .= 1.0  # [c/m]
    #b.c_ar40 = 

    b.rmse_k81 = 1e8     # Initialize to a high value in case time loop exits prematurely
    b.time_k81 = 0.0 
    b.b2.kmin = 1

    k81_decay_constant = decay_constant(229.0)

    try
        # Loop over time and advance model
        for (k, t) in enumerate(time)

            # Get decay tendency
            @inline decay_tendency!(b.dRdt_decay, b.c_k81, dt; λ = k81_decay_constant)

            # Get mixing tendency
            @inline mixing_tendency!(b.dRdt_mixing, b.c_k81, b.mixing_rate, b.thickness; Lref=L_ref)
            for j in b.jj_clean; b.dRdt_mixing[j] = 0.0; end

            # Avoid aging in clean ice beyond t_old time too
            if t > t_old
                for j in b.jj_clean; b.dRdt_decay[j] = 0.0; end
            end

            # Update concentration and ages
            @. b.c_k81 = b.c_k81 + b.dRdt_decay * dt + b.dRdt_mixing * dt

            # Get ages too
            @. b.age_k81 = concentration_to_age(b.c_k81,1.0)

            # Update to current time
            b.time = t

            ## Update summary objects

            if !sampling
                if t in b.b1.times_set
                    # Store time slice of current variables
                    i = findfirst(==(t), b.b1.times)
                    b.b1.age_k81[i,:] = b.age_k81
                    b.b1.c_k81[i,:] = b.c_k81
                end
            end

            # For each depth d of interest, store the value of the variables
            # at the current time
            if t in b.b2.time_set

                #k1 = findfirst(==(t), b.b2.time)
                k1 = round(Int, (t - t0) / 1.0) + 1     # b.b2.time has dt=1.0

                for (i, d) in enumerate(b.b2.depths)
                    # b.b2.age_k81[i,k] = linterp(b.depth, b.age_k81, d)
                    # b.b2.c_k81[i,k] = linterp(b.depth, b.c_k81, d)
                    j = b.b2.interp_idx[i]
                    d_frac = (d - b.depth[j]) / (b.depth[j+1] - b.depth[j])
                    b.b2.age_k81[i,k1] = b.age_k81[j] + d_frac * (b.age_k81[j+1] - b.age_k81[j])
                    b.b2.c_k81[i,k1]   = b.c_k81[j]   + d_frac * (b.c_k81[j+1]   - b.c_k81[j])
                end
                
                # Get sum of squared errors over all observations for current time
                sse = 0.0
                for i in 1:n_obs_k81
                    sse += (b.b2.age_k81[i,k1] - k81.age[i])^2
                end
                b.b2.sse_k81[k1] = sse
                
                # Stop time loop early if sampling and relevant ages are too high already
                if sampling && minimum(@view b.b2.age_k81[:,k1]) >= 1000
                    break
                end
            end

        end
    
    catch e
        e isa DomainError || rethrow(e)  # only swallow DomainErrors, let everything else propagate
        return false            # false = integration failed
    end

    # Calculate metrics
    b.b2.kmin = kmin = argmin(b.b2.sse_k81)
    b.rmse_k81 = sqrt(b.b2.sse_k81[kmin]/n_obs_k81)
    b.time_k81 = b.b2.time[kmin]
    
    if !sampling
        rmse, time, ages = b.rmse_k81, b.time_k81, round.(b.b2.age_k81[:,b.b2.kmin])
        println("k81 (time, rmse, ages): $rmse, $time, $ages")
    end

    return true                 # true = integration succeeded
end

function mixing_rate_discrete!(m, depth, depth_lim, m_clean, m_dirty, delta)
    n = length(m)
    m .= 0.0
    for j in 1:n-1
        depth_interface = 0.5 * (depth[j] + depth[j+1])
        if depth_interface < depth_lim
            # Above transition zone: no mixing
            m[j] = 0.0
        elseif depth_interface < depth_lim + delta
            # Transition zone: clean ice mixing rate
            m[j] = m_clean
        else
            # Fully dirty ice
            m[j] = m_dirty
        end
    end

    m[end] = 0.0  # No mixing at the bottom of last layer

    return m
end

function mixing_rate_smooth!(m, depth, depth_lim, m_clean, m_dirty, delta; sharpness=50.0)
    n = length(m)
    m .= 0.0
    for j in 1:n-1
        d = 0.5 * (depth[j] + depth[j+1])
        w1 = 0.5 * (1 + tanh(sharpness * (d - depth_lim) / delta))         # 0→1 at depth_lim
        w2 = 0.5 * (1 + tanh(sharpness * (d - depth_lim - delta) / delta)) # 0→1 at depth_lim+delta
        m[j] = m_dirty * w1 * w2 + m_clean * w1 * (1 - w2)
    end
    m[end] = 0.0
    return m
end

function mixing_rate_continuous!(m, depth, depth_lim, m_clean, m_dirty, delta)
    n = length(m)
    m .= 0.0
    for j in 1:n-1

        # Get depth of lower cell boundary
        depth_now = 0.5*(depth[j] + depth[j+1])

        if depth_now < depth_lim
            # Clean ice
            m[j] = 0.0
        elseif depth_now < depth_lim + delta
            # Transition to dirty ice
            w = (depth_now - depth_lim)/delta
            m[j] = m_clean*(1-w) + m_dirty*w
        else
            # Fully dirty ice
            m[j] = m_dirty
        end

    end

    m[end] = 0.0        # No mixing at the bottom of last layer

    return m
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

decay_constant(t_half::Float64 = 229.0) = log(2) / t_half

function decay_tendency!(dRdt::Vector{Float64}, R::Vector{Float64}, dt::Float64; λ::Float64 = decay_constant())
    dRdt .= -λ .* R .* exp(-λ * dt)
    return
end

function mixing_tendency!(dRdt::Vector{Float64}, R::Vector{Float64}, Φ::Vector{Float64}, Δz::Vector{Float64}; Lref::Float64=1.0)
    # R: concentration at cell centers [R], length N
    # Φ: mixing rate at lower edge of each cell [m/kyr], length N
    #    - Φ[j] is at the boundary between cell j and cell j+1
    # Lref: length scale of diffusivity [m]
    # Δz: cell thickness [m], length N
    # Returns dR/dt [ [R] kyr⁻¹], length N

    N = length(dRdt)
    
    # Update flux contributions to each cell
    dRdt .= 0.0
    for j in 1:N-1
        D = Φ[j] * Lref                             # Diffusivity: mixing rate at lower boundary [m/kyr] * length scale [m] == [m^2/kyr]
        Δz_interface = 0.5 * (Δz[j] + Δz[j+1])      # Center-to-center distance [m]
        flux = D * (R[j+1] - R[j]) / Δz_interface   # [R/kyr]
        dRdt[j]   += flux / Δz[j]       # [R/kyr]
        dRdt[j+1] -= flux / Δz[j+1]     # [R/kyr]
    end

    return
end

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

function linterp(x, y, xi)
    i = findlast(x .<= xi)
    t = (xi - x[i]) / (x[i+1] - x[i])
    return y[i] + t * (y[i+1] - y[i])
end

### PLOTTING ###

function add_clean_dirty_boundary!(ax,x, y; with_label=true)

    hlines!(ax,y; color=:grey40,linewidth=1,linestyle=:solid)

    if with_label
        text!(ax, "clean ice",
            position = (x, y),
            align = (:right, :bottom),
            fontsize = 8, color = :grey40 )   # x=1.0 means right edge of the axis

        text!(ax, "dirty ice",
            position = (x, y),
            align = (:right, :top),
            fontsize = 8, color = :grey40 )
    end

    return
end
function plot_BasalMixingModelRun(b;k81=nothing,ar40=nothing)

    b1 = b.b1
    b2 = b.b2

    col_k81 = ["#487E3D","#8080F7","teal"]
    col_k81_transparent = [(c, 0.2) for c in col_k81]
    col_ar40 = "#BC401E"
    
    if !isnothing(ar40)
        fig = Figure(size=(1000,600))
    else
        fig = Figure(size=(700,600))
    end

    ## PANEL 0: mixing rate versus depth
    ax0 = Axis(fig[1,1], limits=((-0.05,0.25),(-3053,-3035)), xlabel="Mixing rate (m/yr)", ylabel="Depth (m)", ygridvisible = false )
    colsize!(fig.layout, 1, Auto(0.6))
    d = collect(-3052:2:-3036)
    ax0.yticks = (d,string.(abs.(d)))
    ax0.xticks = [0.0,0.1,0.2]

    add_clean_dirty_boundary!(ax0, 0.98, -b.depth_lim)
    hlines!(ax0,-b.depth;color=(:orange,0.5),linewidth=1.5,linestyle=:dash)
    
    jj = findall(b.depth .>= b.depth_lim)
    scatter!(ax0,b.mixing_rate[jj],-b.depth[jj];color=:black,markersize=5)

    ## PANEL 1: Depth versus closed-system age
    ax1 = Axis(fig[1,end+1], limits=((200,800),(-3053,-3035)), xlabel="⁸¹K closed system age (kyr)", ylabel="Depth (m)", ygridvisible = false )
    d = collect(-3052:2:-3036)
    ax1.yticks = (d,string.(abs.(d)))
    ax1.xticks = [200,400,600,800]

    add_clean_dirty_boundary!(ax1, 0.28, -b.depth_lim,with_label=false)
    hlines!(ax1,-b.depth;color=(:orange,0.5),linewidth=1.5,linestyle=:dash)
    
    # Plot time slices from model
    for (k, t) in enumerate(b1.times)
        lines!(ax1,b1.age_k81[k,:],-b1.depth,color=:grey50,linewidth=0.5)
        if t in [500.0,1000.0,1500.0]
            lines!(ax1,b1.age_k81[k,:],-b1.depth,color=:grey50,linewidth=1.5)
        end
    end

    # Plot data too
    errorbars!(ax1, k81.age, -k81.depth, k81.age_hi, k81.age_lo, color=col_k81, direction=:x, whiskerwidth=8)
    scatter!(ax1, k81.age, -k81.depth, color=col_k81, marker=:circle, markersize=12)
    
    ## PANEL 2 (optional): depth vs d40Ar_ATM concentration
    if !isnothing(ar40)
        ax3 = Axis(fig[1,end+1], limits=((-0.1,0.62),(-3053,-3035)), xlabel="δ⁴⁰ArATM (‰)", ylabel="Depth (m)" )
        d = collect(-3052:2:-3036)
        ax3.yticks = (d,string.(abs.(d)))
        ax3.xticks = 0.0:0.2:0.6

        # Plot data too
        errorbars!(ax3, ar40[!,"δ40/38atm"],-ar40[!,"depth"], ar40[!,"δ40/38atm_err"], ar40[!,"δ40/38atm_err"], color=col_ar40, direction=:x, whiskerwidth=8)
        scatter!(ax3, ar40[!,"δ40/38atm"],-ar40[!,"depth"], color=col_ar40, marker=:circle, markersize=12)
    
    end

    ## PANEL 2 or 3: Closed-system age versus time
    ax2 = Axis(fig[1,end+1], limits=((0,3000),(0,900)), xlabel="Time (kyr)", ylabel="⁸¹K closed system age (kyr)" )
    ax2.xticks = [0,1000,2000,3000]
    ax2.yticks = 0:100:900

    for (i,d) in enumerate(b2.depths)
        lines!(ax2,b2.time,b2.age_k81[i,:],color=col_k81[i],linewidth=2)
    end

    # Plot closed-system age from data
    hlines!(ax2,k81.age;color=col_k81,linewidth=2,linestyle=:dash)
    hspan!(ax2, k81.age .- k81.age_lo, k81.age .+ k81.age_hi; color=col_k81_transparent)

    return fig
end
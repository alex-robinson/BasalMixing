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

function load_basalmixing_data(;depth=[0.0,5000.0])

    # Get data to compare with
    ar40 = CSV.read("data/Bender2010_ar40_data.txt",DataFrame;delim="|",ignorerepeated=true)
    rename!(ar40, strip.(names(ar40)))
    ar40[!,:dar40] = ar40.var"δ40/38atm"
    ar40[!,:dar40_err] = ar40.var"δ40/38atm_err"
    
    # Limit to depth range of interest
    idx = findall(ar40.depth .>= minimum(depth) .&& ar40.depth .<= maximum(depth))
    ar40 = ar40[idx,:]
    
    k81 = CSV.read("data/k81_data.txt",DataFrame;delim=" ",ignorerepeated=true)
    rename!(k81, strip.(names(k81)))
    k81[!,:depth] = 0.5 .* (k81[!,"depth_top"] .+ k81[!,"depth_bottom"])

    # Limit to depth range of interest
    idx = findall(k81.depth .>= minimum(depth) .&& k81.depth .<= maximum(depth))
    k81 = k81[idx,:]

    # Get errors too
    n_obs_k81 = length(k81.age)
    k81_sigma = (sum(k81.age_hi) + sum(k81.age_lo)) / (2*n_obs_k81)
    n_obs_dar40 = length(ar40.dar40)
    dar40_sigma = sum(ar40.dar40_err) / n_obs_dar40

    # Single values, but store in DataFrames
    k81[!,:age_sigma] = fill(k81_sigma,n_obs_k81)
    ar40[!,:dar40_sigma] = fill(dar40_sigma,n_obs_dar40)

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

mutable struct BasalMixingModelState
    time::Vector{Float64}
    age_k81::Array{Float64}
    c_k81::Array{Float64}
    c_ar40::Array{Float64}
    dar40::Array{Float64}
    
    k::Int
    nt::Int
end

function BasalMixingModelState(n; time=[0.0])
    @assert length(time) == 1
    age_k81 = zeros(n)
    c_k81 = ones(n)
    c_ar40 = ones(n)
    dar40 = zeros(n)
    k = 1
    nt = length(time)
    return BasalMixingModelState(time, age_k81, c_k81, c_ar40, dar40,k,nt)
end

function BasalMixingModelState(n,nt; time=zeros(nt))
    @assert length(time) == nt
    age_k81 = zeros(n,nt)
    c_k81 = ones(n,nt)
    c_ar40 = ones(n,nt)
    dar40 = zeros(n,nt)
    k = 1
    return BasalMixingModelState(time, age_k81, c_k81, c_ar40, dar40,k,nt)
end

function reset!(state::BasalMixingModelState)
    state.age_k81 .= 0.0
    state.c_k81 .= 0.0
    state.c_ar40 .= 0.0
    state.dar40 .= 0.0
    state.k = 1
    return
end

mutable struct BasalMixingModelPredictions
    depth::Vector{Float64}          # [nd] depths of interest for this variable
    time::Vector{Float64}           # [nt] time axis
    dat::Array{Float64}             # [nd,nt] predicted value of variable
    rmse::Vector{Float64}           # [nt]
    interp_idx::Vector{Int}         # [nd]
    kmin::Int                       # Time index of minimum error
    rmse_min::Float64               # minimum rmse at t=t[kmin]
    time_min::Float64               # t=t[kmin]

    k::Int
    nd::Int
    nt::Int
end

function BasalMixingModelPredictions(depth::Vector{Float64}, time::Vector{Float64}; depth_ref=nothing)
    nd, nt = length(depth), length(time)
    dat = zeros(nd,nt)
    rmse = fill(1e8,nt)
    if isnothing(depth_ref)
        interp_idx = zeros(nd)
    else
        # Get depth interpolation indices to able to easily match reference depths to predicted depths of interest
        interp_idx = [findlast(depth_ref .<= d) for d in depth]
    end
    kmin = 1
    rmse_min = 1e8
    time_min = 0.0
    k = 1
    return BasalMixingModelPredictions(depth,time,dat,rmse,interp_idx, kmin, rmse_min, time_min, k, nd, nt)
end

function BasalMixingModelPredictions(time::Vector{Float64})
    depth = [0.0]
    nd, nt = 1, length(time)
    dat = zeros(nt)
    rmse = fill(1e8,nt)
    interp_idx = [0]
    kmin = 1
    rmse_min = 1e8
    time_min = 0.0
    k = 1
    return BasalMixingModelPredictions(depth,time,dat,rmse,interp_idx, kmin, rmse_min, time_min, k, nd, nt)
end

function reset!(pred::BasalMixingModelPredictions)
    pred.dat .= 0.0
    pred.rmse .= 0.0
    pred.kmin = 1
    pred.rmse_min = 1e8
    pred.time_min = 0.0
    pred.k = 1
    return
end

function store!(pred::BasalMixingModelPredictions,k,depth,dat)

    @assert k <= pred.nt

    for (i, d) in enumerate(pred.depth)
        j = pred.interp_idx[i]
        d_frac = (d - depth[j]) / (depth[j+1] - depth[j])
        pred.dat[i,k] = dat[j] + d_frac * (dat[j+1] - dat[j])
    end
    
    return
end

mutable struct BasalMixingModel
    n::Int
    depth_lim::Float64
    depth_bedrock::Float64

    layer::Vector{Int}
    depth::Vector{Float64}
    thickness::Vector{Float64}
    mixing_rate::Vector{Float64}
    dRdt_mixing::Vector{Float64}
    dRdt_decay::Vector{Float64}
    idx_clean::Vector{Int}
    
    state::BasalMixingModelState
    states::BasalMixingModelState

    k81::BasalMixingModelPredictions
    dar40::BasalMixingModelPredictions
    joint::BasalMixingModelPredictions
end

function BasalMixingModel(;
    depth = collect(3035.0:3053.0),
    depth_lim = 3040.00,
    depth_bedrock = 3053.44,
    time = 0.0,
    time_states = 500.0:100.0:3000.0,
    time_pred = 0.0:1.0:3000.0,
    k81_obs_depths = [3044.8, 3047.4, 3049.84],
    dar40_obs_depths = [3036.5, 3038.0, 3042.4, 3044.6, 3045.29, 3047.0, 3048.8, 3049.2, 3051.04, 3052.39],
    )

    n = length(depth)
    layers = 1:n
    thickness = cell_thickness(depth,depth_bedrock)
    mixing_rate = fill(0.0,n)
    dRdt_mixing = zeros(n)
    dRdt_decay  = zeros(n)
    idx_clean = findall(depth .<= depth_lim) # Determine clean ice indices
    
    # Get state objects
    state = BasalMixingModelState(n;time=[time])

    nt = length(time_states)
    states = BasalMixingModelState(n, nt; time=time_states)

    # Get prediction objects for comparing with observations

    k81 = BasalMixingModelPredictions(k81_obs_depths, collect(time_pred); depth_ref=depth)
    dar40 = BasalMixingModelPredictions(dar40_obs_depths, collect(time_pred); depth_ref=depth)
    joint = BasalMixingModelPredictions(collect(time_pred))  # No depth for joint comparison

    return BasalMixingModel(
        n,
        depth_lim,
        depth_bedrock,
        collect(layers),
        collect(depth),
        collect(thickness),
        collect(mixing_rate),
        dRdt_mixing,
        dRdt_decay,
        idx_clean,
        state,
        states,
        k81,
        dar40,
        joint
    )
end

function ResetBasalMixingModel!(b)

    reset!(b.state)
    reset!(b.states)
    reset!(b.k81)
    reset!(b.dar40)
    reset!(b.joint)

    return
end

function RunBasalMixingModel!(p, b, dat; t0=0.0,t1=3000.0,dt=1.0,sampling=false)

    # Extract model parameters
    (delta, m_clean, f_dirty, t_old, F_ar40) = p
    L_ref = 1.0     # [m] Use L_ref=1.0, since this just scales m_clean, can tune m_clean directly

    k81_decay_constant = decay_constant(229.0)

    # Extract data for comparison
    (k81_obs_df, dar40_obs_df) = dat

    # Check consistency with our BasalMixingModel info
    @assert b.k81.depth == k81_obs_df.depth
    @assert b.dar40.depth == dar40_obs_df.depth
    
    k81_obs   = k81_obs_df.age
    k81_sigma = k81_obs_df.age_sigma[1]
    k81_var   = k81_sigma^2
    dar40_obs   = dar40_obs_df.dar40
    dar40_sigma = dar40_obs_df.dar40_sigma[1]
    dar40_var   = dar40_sigma^2

    n_obs_k81 = length(k81_obs)
    n_obs_dar40 = length(dar40_obs)

    # Set the mixing rate
    mixing_rate_smooth!(b.mixing_rate, b.depth, b.depth_lim, m_clean, m_clean * f_dirty, delta)

    # Generate all times to simulate
    time = t0:dt:t1

    ## Set initial values ##

    ResetBasalMixingModel!(b)

    b.state.c_k81 .= 1.0  # [c/m]

    Ar40_00 = calc_ar40_with_aging(0.0, 0.0)    # Modern concentration [cc/m³]
    b.state.c_ar40 .= Ar40_00                         # store uniform modern value initially
    
    @. b.state.dar40 = calc_delta_ar40(b.state.c_ar40, Ar40_00, t_old)

    # Loop over time and advance model
    try
        for t in time

            # Get decay tendency
            @inline decay_tendency!(b.dRdt_decay, b.state.c_k81, dt; λ = k81_decay_constant)

            # Get mixing tendency
            @inline mixing_tendency!(b.dRdt_mixing, b.state.c_k81, b.mixing_rate, b.thickness; Lref=L_ref)
            for j in b.idx_clean; b.dRdt_mixing[j] = 0.0; end

            # Avoid aging in clean ice beyond t_old time too
            if t > t_old
                for j in b.idx_clean; b.dRdt_decay[j] = 0.0; end
            end

            # Update concentration and ages
            @. b.state.c_k81 = b.state.c_k81 + b.dRdt_decay * dt + b.dRdt_mixing * dt

            # Get ages too
            @. b.state.age_k81 = concentration_to_age(b.state.c_k81,1.0)
            
            ## AR40 ##
            
            # Ar40 mixing tendency, overwrites b.dRdt_mixing with Ar40 values
            b.dRdt_mixing .= 0.0
            @inline mixing_tendency!(b.dRdt_mixing, b.state.c_ar40, b.mixing_rate, b.thickness; Lref=L_ref)
            for j in b.idx_clean; b.dRdt_mixing[j] = 0.0; end   # no mixing in clean ice

            # Bottom source flux into the deepest box only
            b.dRdt_mixing[end] += F_ar40 / b.thickness[end]

            # Advance Ar40
            @. b.state.c_ar40 = b.state.c_ar40 + b.dRdt_mixing * dt

            # Update delta Ar40
            @. b.state.dar40 = calc_delta_ar40(b.state.c_ar40, Ar40_00, t_old)

            ##########

            # Update to current time
            b.state.time .= t

            ## Update summary objects

            if !sampling && t == b.states.time[b.states.k]
                # Store time slice of current variables
                k = b.states.k
                b.states.age_k81[:,k] = b.state.age_k81
                b.states.c_k81[:,k] = b.state.c_k81
                b.states.c_ar40[:,k] = b.state.c_ar40
                b.states.dar40[:,k] = b.state.dar40
                b.states.k += 1 # Advance to next time index
            end

            # Ar40 ##

            if t == b.dar40.time[b.dar40.k]
                k = b.dar40.k

                store!(b.dar40,k,b.depth,b.state.dar40)
                
                # Get sum of squared errors over all observations for current time
                sse = 0.0
                for i in 1:n_obs_dar40
                    sse += (b.dar40.dat[i,k] - dar40_obs[i])^2 / dar40_var
                end
                b.dar40.rmse[k] = sqrt(sse/n_obs_dar40)

                b.dar40.k += 1
            end

            # K81 ##
            
            if t == b.k81.time[b.k81.k]
                k = b.k81.k

                store!(b.k81,k,b.depth,b.state.age_k81)
                
                # Get sum of squared errors over all observations for current time
                sse = 0.0
                for i in 1:n_obs_k81
                    sse += (b.k81.dat[i,k] - k81_obs[i])^2 / k81_var
                end
                b.k81.rmse[k] = sqrt(sse/n_obs_k81)
                
                # Stop time loop early if sampling and relevant ages are too high already
                if sampling && minimum(@view b.k81.dat[:,k]) >= 1000
                    break
                end

                b.k81.k += 1
            end

        end
    
    catch e
        e isa DomainError || rethrow(e)  # only swallow DomainErrors, let everything else propagate
        return false            # false = integration failed
    end

    # Calculate summary metrics
    b.k81.kmin = kmin = argmin(b.k81.rmse)
    b.k81.rmse_min = b.k81.rmse[kmin]
    b.k81.time_min = b.k81.time[kmin]

    b.dar40.kmin = kmin = argmin(b.dar40.rmse)
    b.dar40.rmse_min = b.dar40.rmse[kmin]
    b.dar40.time_min = b.dar40.time[kmin]

    b.joint.rmse .= b.k81.rmse .+ b.dar40.rmse
    b.joint.kmin = kmin = argmin(b.joint.rmse)
    b.joint.rmse_min = b.joint.rmse[kmin]
    b.joint.time_min = b.joint.time[kmin]

    if !sampling
        kmin = b.joint.kmin
        rmse_min, rmses, time_min, ages = round(b.joint.rmse_min,digits=3), round.([b.k81.rmse_min, b.dar40.rmse_min],digits=3), b.joint.time_min, round.(b.k81.dat[:,kmin])
        println("k81&dAr40 (rmse_min, rmses, time_min, ages): $rmse_min, $rmses, $time_min, $ages")
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

function calc_delta_ar40(ar40::Float64, ar40_ref::Float64, t_kyr)
    return (ar40 / ar40_ref - 1.0) * 1000.0 - (0.066/1000.0) * max(t_kyr, 0.0)
end

function calc_ar40_with_aging(t_kyr::Float64, t_old::Float64;
    air_content::Float64 = 0.08,        # fraction
    f_ar40_atm_pd::Float64 = 0.00934    # atmospheric ⁴⁰Ar volume fraction today
    )
    # Assumes δ40ar=0, that surface concentration is in equilibrium with the air at that time.

    # Reference ⁴⁰Ar content of the layer [cc m⁻³]
    # volume [100^3 cc m⁻³] × air_content × f_ar40_atm
    ar40 = 100^3 * (air_content * f_ar40_atm_pd) * (1 - 0.066e-6 * max(t_kyr-t_old, 0.0) )

    #%% calculate amount of argon (in ccs) starting in ice
    #TAC=0.08; % average TAC in basal ice
    #Ar40cc_mod=TAC*.00934*100^3; % ccs of 40 argon per cubic meter for modern ice

    return ar40
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
function plot_BasalMixingModelRun(b;k81_obs=nothing,dar40_obs=nothing)

    states = b.states
    k81 = b.k81
    dar40 = b.dar40
    joint = b.joint

    col_k81 = ["#487E3D","#8080F7","teal"]
    col_k81_transparent = [(c, 0.2) for c in col_k81]
    col_dar40 = "#BC401E"
    
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
    for (k, t) in enumerate(states.time)
        lines!(ax1,states.age_k81[:,k],-b.depth,color=:grey50,linewidth=0.5)
        if t in [500.0,1000.0,1500.0]
            lines!(ax1,states.age_k81[:,k],-b.depth,color=:grey50,linewidth=1.5)
        end
    end
    k = argmin(abs.(states.time .- joint.time_min))
    lines!(ax1,states.age_k81[:,k],-b.depth,color=:black,linewidth=2)

    # Plot data too
    errorbars!(ax1, k81_obs.age, -k81_obs.depth, k81_obs.age_hi, k81_obs.age_lo, color=col_k81, direction=:x, whiskerwidth=8)
    scatter!(ax1, k81_obs.age, -k81_obs.depth, color=col_k81, marker=:circle, markersize=12)
    
    ## PANEL 2 (optional): depth vs d40Ar_atm concentration
    if !isnothing(ar40)
        ax3 = Axis(fig[1,end+1], limits=((-0.1,0.62),(-3053,-3035)), xlabel="δ⁴⁰ArATM (‰)", ylabel="Depth (m)" )
        d = collect(-3052:2:-3036)
        ax3.yticks = (d,string.(abs.(d)))
        ax3.xticks = 0.0:0.2:0.6

        # Plot time slices from model
        for (k, t) in enumerate(states.time)
            lines!(ax3,states.dar40[:,k],-b.depth,color=:grey50,linewidth=0.5)
            if t in [500.0,1000.0,1500.0]
                lines!(ax3,states.dar40[:,k],-b.depth,color=:grey50,linewidth=1.5)
            end
        end
        k = argmin(abs.(states.time .- joint.time_min))
        lines!(ax3,states.dar40[:,k],-b.depth,color=:black,linewidth=2.5)

        # Plot data too
        errorbars!(ax3, dar40_obs[!,:dar40],-dar40_obs[!,"depth"], dar40_obs[!,:dar40_err], dar40_obs[!,:dar40_err], color=col_dar40, direction=:x, whiskerwidth=8)
        scatter!(ax3, dar40_obs[!,:dar40],-dar40_obs[!,"depth"], color=col_dar40, marker=:circle, markersize=12)
    end

    ## PANEL 2 or 3: Closed-system age versus time
    ax2 = Axis(fig[1,end+1], limits=((0,3000),(0,900)), xlabel="Time (kyr)", ylabel="⁸¹K closed system age (kyr)" )
    ax2.xticks = [0,1000,2000,3000]
    ax2.yticks = 0:100:900

    for (j,d) in enumerate(k81.depth)
        lines!(ax2,k81.time,k81.dat[j,:],color=col_k81[j],linewidth=2)
    end

    # Plot closed-system age from data
    hlines!(ax2,k81_obs.age;color=col_k81,linewidth=2,linestyle=:dash)
    hspan!(ax2, k81_obs.age .- k81_obs.age_lo, k81_obs.age .+ k81_obs.age_hi; color=col_k81_transparent)

    return fig
end
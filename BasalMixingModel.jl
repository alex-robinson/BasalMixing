using CairoMakie
using DifferentialEquations
using JLD2
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
    rmse = fill(Inf, nt)
    if isnothing(depth_ref)
        interp_idx = zeros(nd)
    else
        # Get depth interpolation indices to able to easily match reference depths to predicted depths of interest
        interp_idx = [findlast(depth_ref .<= d) for d in depth]
    end
    kmin = 1
    rmse_min = Inf
    time_min = 0.0
    k = 1
    return BasalMixingModelPredictions(depth,time,dat,rmse,interp_idx, kmin, rmse_min, time_min, k, nd, nt)
end

function BasalMixingModelPredictions(time::Vector{Float64})
    depth = [0.0]
    nd, nt = 1, length(time)
    dat = zeros(nt)
    rmse = fill(Inf, nt)
    interp_idx = [0]
    kmin = 1
    rmse_min = Inf
    time_min = 0.0
    k = 1
    return BasalMixingModelPredictions(depth,time,dat,rmse,interp_idx, kmin, rmse_min, time_min, k, nd, nt)
end

function reset!(pred::BasalMixingModelPredictions)
    pred.dat .= 0.0
    # Use Inf (not 0) as the "not yet computed" sentinel so argmin(pred.rmse)
    # later never picks a slot that the integration didn't actually fill —
    # important when t1 stops short of the prediction grid (e.g. the :sampled
    # path which only integrates 0 → |t_0|). Matches the constructor.
    pred.rmse .= Inf
    pred.kmin = 1
    pred.rmse_min = Inf
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

"""
    RunBasalMixingModel!(p, b, dat; t0=0.0, t1=3000.0, dt=1.0)

Trajectory-tracing forward run used by the plotting code. Integrates the model
on `t0:dt:t1` and records the full trajectory into `b.states`, `b.k81.dat`,
and `b.dar40.dat`. Also computes per-slice RMSE against the observations
(handy as a diagnostic).

Used only for visualisation — inference goes through
`RunBasalMixingModelToTime!`, which evaluates a single endpoint.

Returns `true` on success, `false` on a `DomainError` (concentration went
non-positive, etc.); other exceptions propagate.
"""
function RunBasalMixingModel!(p, b, dat; t0=0.0,t1=3000.0,dt=1.0)

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

    try
        for t in time

            @inline decay_tendency!(b.dRdt_decay, b.state.c_k81, dt; λ = k81_decay_constant)
            @inline mixing_tendency!(b.dRdt_mixing, b.state.c_k81, b.mixing_rate, b.thickness; Lref=L_ref, idx_clean=b.idx_clean)

            if t > t_old
                for j in b.idx_clean; b.dRdt_decay[j] = 0.0; end
            end

            @. b.state.c_k81 = b.state.c_k81 + b.dRdt_decay * dt + b.dRdt_mixing * dt
            @. b.state.age_k81 = concentration_to_age(b.state.c_k81,1.0)

            @inline mixing_tendency!(b.dRdt_mixing, b.state.c_ar40, b.mixing_rate, b.thickness; Lref=L_ref, idx_clean=b.idx_clean)
            b.dRdt_mixing[end] += F_ar40 / b.thickness[end]
            @. b.state.c_ar40 = b.state.c_ar40 + b.dRdt_mixing * dt
            @. b.state.dar40 = calc_delta_ar40(b.state.c_ar40, Ar40_00, t_old)

            b.state.time .= t

            if t == b.states.time[b.states.k]
                k = b.states.k
                b.states.age_k81[:,k] = b.state.age_k81
                b.states.c_k81[:,k] = b.state.c_k81
                b.states.c_ar40[:,k] = b.state.c_ar40
                b.states.dar40[:,k] = b.state.dar40
                b.states.k += 1
            end

            if t == b.dar40.time[b.dar40.k]
                k = b.dar40.k
                store!(b.dar40,k,b.depth,b.state.dar40)
                sse = 0.0
                for i in 1:n_obs_dar40
                    sse += (b.dar40.dat[i,k] - dar40_obs[i])^2 / dar40_var
                end
                b.dar40.rmse[k] = sqrt(sse/n_obs_dar40)
                b.dar40.k += 1
            end

            if t == b.k81.time[b.k81.k]
                k = b.k81.k
                store!(b.k81,k,b.depth,b.state.age_k81)
                sse = 0.0
                for i in 1:n_obs_k81
                    sse += (b.k81.dat[i,k] - k81_obs[i])^2 / k81_var
                end
                b.k81.rmse[k] = sqrt(sse/n_obs_k81)
                b.k81.k += 1
            end

        end

    catch e
        e isa DomainError || rethrow(e)
        return false
    end

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

    return true
end

"""
    RunBasalMixingModelToTime!(p, b, t_target, dat; t0=0.0, dt=0.2)

Sampled-time variant of `RunBasalMixingModel!`. Integrates the forward model
from `t0` up to `t_target` (a sampled parameter, not a grid) and stores the
final-state predictions at the observation depths in slot 1 of `b.k81.dat`
and `b.dar40.dat`. No RMSE / argmin scan — the only "time of comparison"
is `t_target`.

Returns `true` on success, `false` on a `DomainError` (concentration went
non-positive, etc.); other exceptions propagate as in `RunBasalMixingModel!`.

`b.joint.time_min`, `b.k81.time_min`, `b.dar40.time_min`, and the `kmin`
fields are set to `t_target` / `1` respectively so the existing plotting code
(which reads `b.joint.kmin` etc.) still works.
"""
function RunBasalMixingModelToTime!(p, b, t_target, dat; t0::Float64=0.0, dt::Float64=0.2)

    (delta, m_clean, f_dirty, t_old, F_ar40) = p
    L_ref = 1.0
    k81_decay_constant = decay_constant(229.0)

    (k81_obs_df, dar40_obs_df) = dat
    @assert b.k81.depth == k81_obs_df.depth
    @assert b.dar40.depth == dar40_obs_df.depth

    # Set the mixing rate from the parameters
    mixing_rate_smooth!(b.mixing_rate, b.depth, b.depth_lim, m_clean, m_clean * f_dirty, delta)

    ## Initial values ##
    ResetBasalMixingModel!(b)
    b.state.c_k81 .= 1.0
    Ar40_00 = calc_ar40_with_aging(0.0, 0.0)
    b.state.c_ar40 .= Ar40_00
    @. b.state.dar40 = calc_delta_ar40(b.state.c_ar40, Ar40_00, t_old)

    # Integrate t0 -> t_target. The grid may not land exactly on t_target; the
    # discretization error is at most dt and is absorbed in the likelihood's σ.
    time = t0:dt:t_target

    try
        for t in time
            # k81 decay
            @inline decay_tendency!(b.dRdt_decay, b.state.c_k81, dt; λ = k81_decay_constant)
            # k81 mixing
            @inline mixing_tendency!(b.dRdt_mixing, b.state.c_k81, b.mixing_rate, b.thickness;
                                     Lref=L_ref, idx_clean=b.idx_clean)
            # No aging in clean ice past t_old
            if t > t_old
                for j in b.idx_clean; b.dRdt_decay[j] = 0.0; end
            end
            @. b.state.c_k81 = b.state.c_k81 + b.dRdt_decay * dt + b.dRdt_mixing * dt
            @. b.state.age_k81 = concentration_to_age(b.state.c_k81, 1.0)

            # Ar40 mixing (reuses b.dRdt_mixing buffer)
            @inline mixing_tendency!(b.dRdt_mixing, b.state.c_ar40, b.mixing_rate, b.thickness;
                                     Lref=L_ref, idx_clean=b.idx_clean)
            b.dRdt_mixing[end] += F_ar40 / b.thickness[end]
            @. b.state.c_ar40 = b.state.c_ar40 + b.dRdt_mixing * dt
            @. b.state.dar40  = calc_delta_ar40(b.state.c_ar40, Ar40_00, t_old)

            b.state.time .= t
        end
    catch e
        e isa DomainError || rethrow(e)
        return false
    end

    # Drop the final-state predictions into slot 1 of b.k81 / b.dar40.
    store!(b.k81,   1, b.depth, b.state.age_k81)
    store!(b.dar40, 1, b.depth, b.state.dar40)

    # Make the existing plot code (which reads b.joint.kmin, b.joint.time_min,
    # b.k81.dat[:, b.joint.kmin], ...) still produce the right slice.
    b.k81.kmin     = 1
    b.dar40.kmin   = 1
    b.joint.kmin   = 1
    b.k81.time_min   = float(t_target)
    b.dar40.time_min = float(t_target)
    b.joint.time_min = float(t_target)

    return true
end

### ODE-based forward model (NUTS / ForwardDiff path) ##############
##
## The mutable-buffer code above is the Emcee path: hand-rolled Euler in
## positive elapsed time, all buffers `Vector{Float64}`. It can't take a
## `Dual`-valued `t_target` (`0:dt:t_target` breaks) and the `Vector{Float64}`
## buffers can't hold ForwardDiff `Dual`s anyway.
##
## The functions below are the AD-friendly path: allocating, `Vector{T}`-typed,
## and call `OrdinaryDiffEq.Tsit5` so `tspan = (0, t_target)` is differentiable
## w.r.t. `t_target`. Geometry (depth, thickness, idx_clean_mask) stays
## `Float64`; only the state and the (delta, m_clean, f_dirty, t_old, F_ar40)
## dependent quantities are typed by `T`.

"""
    compute_mixing_rate_smooth(depth, depth_lim, m_clean, m_dirty, delta; sharpness=50.0)

Allocating, AD-friendly twin of `mixing_rate_smooth!`. Returns a fresh
`Vector{T}` where `T = promote_type(typeof(m_clean), typeof(m_dirty), typeof(delta))`.
Matches `mixing_rate_smooth!` exactly for `Float64` inputs.
"""
function compute_mixing_rate_smooth(depth::AbstractVector{<:Real},
                                    depth_lim::Real,
                                    m_clean::Real, m_dirty::Real, delta::Real;
                                    sharpness::Real=50.0)
    T = promote_type(typeof(m_clean), typeof(m_dirty), typeof(delta))
    n = length(depth)
    m = zeros(T, n)
    for j in 1:n-1
        d = 0.5 * (depth[j] + depth[j+1])
        w1 = 0.5 * (1 + tanh(sharpness * (d - depth_lim) / delta))
        w2 = 0.5 * (1 + tanh(sharpness * (d - depth_lim - delta) / delta))
        m[j] = m_dirty * w1 * w2 + m_clean * w1 * (1 - w2)
    end
    # m[end] = 0 already from `zeros`
    return m
end

"""
    basal_mixing_rhs!(du, u, p, t)

ODE right-hand side for the basal-mixing forward model. State vector `u`
packs `c_k81` (entries 1:N) and `c_ar40` (entries N+1:2N). Parameters `p` is
a `NamedTuple` carrying both the sampled scalars (`t_old`, `F_ar40`) and the
precomputed geometry (`thickness`, `mixing_rate`, `idx_clean_mask`,
`k81_decay_constant`, `N`). `mixing_rate` is precomputed once per likelihood
eval by `compute_mixing_rate_smooth` and threaded in via `p`.

Mirrors the Euler integrator's semantics exactly:
- k81 decays at rate λ, zeroed for clean-ice cells once `t > t_old`.
- Mixing flux is centred-difference between adjacent cells; clean-ice cells
  receive no mixing flux (their endpoint contributions are skipped).
- ⁴⁰Ar receives a bottom-source flux `F_ar40 / thickness[end]` in the
  deepest cell.
"""
function basal_mixing_rhs!(du, u, p, t)
    N = p.N
    c_k81  = @view u[1:N]
    c_ar40 = @view u[N+1:2N]
    du_k81  = @view du[1:N]
    du_ar40 = @view du[N+1:2N]

    λ          = p.k81_decay_constant
    t_old      = p.t_old
    F_ar40     = p.F_ar40
    thickness  = p.thickness
    mixing_rate = p.mixing_rate
    clean_mask  = p.idx_clean_mask
    Tz = eltype(du)
    zero_T = zero(Tz)

    # k81 decay (zeroed for clean ice past t_old). ar40 starts at zero.
    decay_kills_clean = (t > t_old)
    @inbounds for j in 1:N
        if decay_kills_clean && clean_mask[j]
            du_k81[j] = zero_T
        else
            du_k81[j] = -λ * c_k81[j]
        end
        du_ar40[j] = zero_T
    end

    # Mixing flux between adjacent cells; skip endpoints that are clean-ice.
    @inbounds for j in 1:N-1
        Δz_interface = 0.5 * (thickness[j] + thickness[j+1])
        D = mixing_rate[j]                  # already has Lref=1 absorbed
        flux_k81  = D * (c_k81[j+1]  - c_k81[j])  / Δz_interface
        flux_ar40 = D * (c_ar40[j+1] - c_ar40[j]) / Δz_interface
        if !clean_mask[j]
            du_k81[j]   += flux_k81  / thickness[j]
            du_ar40[j]  += flux_ar40 / thickness[j]
        end
        if !clean_mask[j+1]
            du_k81[j+1]  -= flux_k81  / thickness[j+1]
            du_ar40[j+1] -= flux_ar40 / thickness[j+1]
        end
    end

    # ar40 bottom source flux into the deepest cell only.
    du_ar40[N] += F_ar40 / thickness[N]

    return nothing
end

"""
    BasalMixingGeometry(depth, depth_lim, depth_bedrock, k81_obs_depths, dar40_obs_depths)
        -> NamedTuple

Precompute everything that's `Float64`-only and depends solely on the depth
grid (not on sampled parameters). Pass the returned NamedTuple to
`run_basal_mixing_ode` to avoid recomputing it for every likelihood eval.

Fields: `N`, `depth`, `depth_lim`, `thickness`, `idx_clean_mask`,
`k81_interp` (`(j_idx, frac)` tuples for each obs depth), `dar40_interp`.
"""
function BasalMixingGeometry(depth::Vector{Float64}, depth_lim::Float64,
                             depth_bedrock::Float64,
                             k81_obs_depths::Vector{Float64},
                             dar40_obs_depths::Vector{Float64})
    N = length(depth)
    thickness = cell_thickness(depth, depth_bedrock)
    idx_clean_mask = depth .<= depth_lim   # BitVector

    function _interp(obs_depths)
        out = Tuple{Int,Float64}[]
        for d in obs_depths
            j = findlast(depth .<= d)
            (isnothing(j) || j == N) && error("obs depth $d outside grid")
            frac = (d - depth[j]) / (depth[j+1] - depth[j])
            push!(out, (j, frac))
        end
        return out
    end

    k81_interp   = _interp(k81_obs_depths)
    dar40_interp = _interp(dar40_obs_depths)

    return (; N, depth, depth_lim, thickness, idx_clean_mask, k81_interp, dar40_interp)
end

# Helper: linear interpolation of a field `f` at the precomputed obs locations.
@inline function _interp_at_obs(f::AbstractVector, interp::Vector{Tuple{Int,Float64}})
    T = eltype(f)
    out = Vector{T}(undef, length(interp))
    @inbounds for (i, (j, frac)) in enumerate(interp)
        out[i] = f[j] + frac * (f[j+1] - f[j])
    end
    return out
end

"""
    run_basal_mixing_ode(params, geom; solver=Tsit5(), reltol=1e-6, abstol=1e-8)
        -> (k81_age_pred, dar40_pred, success)

ODE-based forward model used by the NUTS @model. Integrates the coupled
(k81, ar40) system on positive elapsed time `0 → params.t_target` using
`OrdinaryDiffEq`, then linearly interpolates the predicted ⁸¹Kr closed-system
age and δ⁴⁰Ar at the observation depths.

`params` is a NamedTuple with fields `delta, m_clean, f_dirty, t_old, F_ar40,
t_target`. All are `Real` (possibly `Dual`). `geom` is what
`BasalMixingGeometry` returned.

`success` is `true` if the ODE integration finished cleanly. On failure
(e.g. solver max-iters hit), returns predictions filled with `1e8` so the
likelihood is driven to ~-Inf without throwing.
"""
function run_basal_mixing_ode(params, geom;
                              solver=Tsit5(),
                              reltol::Float64=1e-6,
                              abstol::Float64=1e-8)
    delta    = params.delta
    m_clean  = params.m_clean
    f_dirty  = params.f_dirty
    t_old    = params.t_old
    F_ar40   = params.F_ar40
    t_target = params.t_target

    N         = geom.N
    depth     = geom.depth
    depth_lim = geom.depth_lim
    thickness = geom.thickness
    clean_mask = geom.idx_clean_mask

    T = promote_type(typeof(delta), typeof(m_clean), typeof(f_dirty),
                     typeof(t_old), typeof(F_ar40), typeof(t_target))

    # Precompute mixing_rate (T-valued).
    mixing_rate = compute_mixing_rate_smooth(depth, depth_lim,
                                             m_clean, m_clean*f_dirty, delta)

    # Modern ⁴⁰Ar reference for the initial condition (Float64).
    Ar40_00 = calc_ar40_with_aging(0.0, 0.0)

    # Initial state: c_k81 = 1, c_ar40 = Ar40_00.
    u0 = Vector{T}(undef, 2N)
    @inbounds for j in 1:N
        u0[j]     = one(T)
        u0[N+j]   = T(Ar40_00)
    end

    p_ode = (; N,
             k81_decay_constant = decay_constant(229.0),
             t_old, F_ar40,
             thickness, mixing_rate, idx_clean_mask=clean_mask)

    tspan = (zero(T), t_target)
    prob = ODEProblem{true}(basal_mixing_rhs!, u0, tspan, p_ode)

    sol = solve(prob, solver;
                reltol=reltol, abstol=abstol,
                save_everystep=false, save_start=false,
                saveat=[t_target])

    if sol.retcode !== ReturnCode.Success
        return (fill(T(1e8), length(geom.k81_interp)),
                fill(T(1e8), length(geom.dar40_interp)),
                false)
    end

    u_final = sol.u[end]
    c_k81_final  = @view u_final[1:N]
    c_ar40_final = @view u_final[N+1:2N]

    # Concentration → age, Ar40 → δ⁴⁰Ar.
    # `concentration_to_age` is `log2`, which DomainErrors on c ≤ 0. With
    # reasonable params and a well-resolved Tsit5 solve this shouldn't happen;
    # if it does we want a finite (huge) prediction rather than an exception.
    age_k81_final = similar(c_k81_final)
    @inbounds for j in 1:N
        c = c_k81_final[j]
        age_k81_final[j] = c > 0 ? -229.0 * log2(c) : T(1e8)
    end
    dar40_final = similar(c_ar40_final)
    @inbounds for j in 1:N
        dar40_final[j] = calc_delta_ar40(c_ar40_final[j], Ar40_00, t_old)
    end

    k81_pred   = _interp_at_obs(age_k81_final, geom.k81_interp)
    dar40_pred = _interp_at_obs(dar40_final,   geom.dar40_interp)

    return k81_pred, dar40_pred, true
end

####################################################################

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

function mixing_tendency!(dRdt::Vector{Float64}, R::Vector{Float64}, Φ::Vector{Float64}, Δz::Vector{Float64}; Lref::Float64=1.0, idx_clean=1)
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

    dRdt[idx_clean] .= 0.0

    return
end

function calc_delta_ar40(ar40::Real, ar40_ref::Real, t_kyr::Real)
    return (ar40 / ar40_ref - 1.0) * 1000.0 - (0.066/1000.0) * max(t_kyr, zero(t_kyr))
end

function calc_ar40_with_aging(t_kyr::Real, t_old::Real;
    air_content::Real = 0.08,        # fraction
    f_ar40_atm_pd::Real = 0.00934    # atmospheric ⁴⁰Ar volume fraction today
    )
    # Assumes δ40ar=0, that surface concentration is in equilibrium with the air at that time.

    # Reference ⁴⁰Ar content of the layer [cc m⁻³]
    # volume [100^3 cc m⁻³] × air_content × f_ar40_atm
    ar40 = 100^3 * (air_content * f_ar40_atm_pd) * (1 - 0.066e-6 * max(t_kyr-t_old, zero(t_kyr-t_old)) )

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
"""
    plot_BasalMixingModelRun(b; k81_obs, dar40_obs, t_max=nothing, overlay=nothing)

Render the standard four-panel figure for a single forward run `b`.

`t_max` truncates the time-series and snapshot panels so that runs which only
integrated 0 → t_max (e.g. an MAP re-run with duration |t_0|) don't show
uninitialised buffer tails as a misleading zero trajectory.

`overlay`, if given, is a NamedTuple `(b=b2, t_max=t_max2, label="kr81")`
whose best-fit age and δ⁴⁰Ar are drawn on top of the primary in a distinct
dashed style, so the secondary chain (e.g. kr81-only) can be compared
visually to the primary (combined).
"""
function plot_BasalMixingModelRun(b;k81_obs=nothing,dar40_obs=nothing,t_max=nothing,overlay=nothing)

    states = b.states
    k81 = b.k81
    dar40 = b.dar40
    joint = b.joint

    t_max_val = t_max === nothing ? Float64(k81.time[end]) : Float64(t_max)
    k_last_pred  = searchsortedlast(k81.time, t_max_val)
    k_last_state = searchsortedlast(states.time, t_max_val)

    col_k81 = ["#487E3D","#8080F7","teal"]
    col_k81_transparent = [(c, 0.2) for c in col_k81]
    col_dar40 = "#BC401E"
    col_overlay = :black

    if !isnothing(dar40_obs)
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

    if overlay !== nothing
        b2 = overlay.b
        jj2 = findall(b2.depth .>= b2.depth_lim)
        lines!(ax0, b2.mixing_rate[jj2], -b2.depth[jj2];
               color=col_overlay, linewidth=1.5, linestyle=:dash)
    end

    ## PANEL 1: Depth versus closed-system age
    ax1 = Axis(fig[1,end+1], limits=((200,800),(-3053,-3035)), xlabel="⁸¹K closed system age (kyr)", ylabel="Depth (m)", ygridvisible = false )
    d = collect(-3052:2:-3036)
    ax1.yticks = (d,string.(abs.(d)))
    ax1.xticks = [200,400,600,800]

    add_clean_dirty_boundary!(ax1, 0.28, -b.depth_lim,with_label=false)
    hlines!(ax1,-b.depth;color=(:orange,0.5),linewidth=1.5,linestyle=:dash)

    for k in 1:k_last_state
        t = states.time[k]
        lines!(ax1,states.age_k81[:,k],-b.depth,color=:grey50,linewidth=0.5)
        if t in [500.0,1000.0,1500.0]
            lines!(ax1,states.age_k81[:,k],-b.depth,color=:grey50,linewidth=1.5)
        end
    end
    k = argmin(abs.(states.time[1:k_last_state] .- joint.time_min))
    lines!(ax1,states.age_k81[:,k],-b.depth,color=:black,linewidth=2)

    if overlay !== nothing
        b2 = overlay.b
        t2 = overlay.t_max === nothing ? Float64(b2.k81.time[end]) : Float64(overlay.t_max)
        k_last_state2 = searchsortedlast(b2.states.time, t2)
        k2 = argmin(abs.(b2.states.time[1:k_last_state2] .- b2.joint.time_min))
        lines!(ax1, b2.states.age_k81[:,k2], -b2.depth;
               color=col_overlay, linewidth=2, linestyle=:dash)
    end

    errorbars!(ax1, k81_obs.age, -k81_obs.depth, k81_obs.age_hi, k81_obs.age_lo, color=col_k81, direction=:x, whiskerwidth=8)
    scatter!(ax1, k81_obs.age, -k81_obs.depth, color=col_k81, marker=:circle, markersize=12)

    ## PANEL 2 (optional): depth vs d40Ar_atm concentration
    if !isnothing(dar40_obs)
        ax3 = Axis(fig[1,end+1], limits=((-0.1,0.62),(-3053,-3035)), xlabel="δ⁴⁰ArATM (‰)", ylabel="Depth (m)" )
        d = collect(-3052:2:-3036)
        ax3.yticks = (d,string.(abs.(d)))
        ax3.xticks = 0.0:0.2:0.6

        for k in 1:k_last_state
            t = states.time[k]
            lines!(ax3,states.dar40[:,k],-b.depth,color=:grey50,linewidth=0.5)
            if t in [500.0,1000.0,1500.0]
                lines!(ax3,states.dar40[:,k],-b.depth,color=:grey50,linewidth=1.5)
            end
        end
        k = argmin(abs.(states.time[1:k_last_state] .- joint.time_min))
        lines!(ax3,states.dar40[:,k],-b.depth,color=:black,linewidth=2.5)

        if overlay !== nothing
            b2 = overlay.b
            t2 = overlay.t_max === nothing ? Float64(b2.k81.time[end]) : Float64(overlay.t_max)
            k_last_state2 = searchsortedlast(b2.states.time, t2)
            k2 = argmin(abs.(b2.states.time[1:k_last_state2] .- b2.joint.time_min))
            lines!(ax3, b2.states.dar40[:,k2], -b2.depth;
                   color=col_overlay, linewidth=2, linestyle=:dash)
        end

        errorbars!(ax3, dar40_obs[!,:dar40],-dar40_obs[!,"depth"], dar40_obs[!,:dar40_err], dar40_obs[!,:dar40_err], color=col_dar40, direction=:x, whiskerwidth=8)
        scatter!(ax3, dar40_obs[!,:dar40],-dar40_obs[!,"depth"], color=col_dar40, marker=:circle, markersize=12)
    end

    ## PANEL 2 or 3: Closed-system age versus time
    # Display in signed kyr (negative = past, 0 = present day). The integrator
    # runs in positive elapsed time τ ∈ [0, t_max_val]; the curve shown here is
    # x(τ) = τ - t_max_val, so the trajectory starts at t_0 = -t_max_val and
    # ends at t_end = 0.
    ax2 = Axis(fig[1,end+1], limits=((-3000,0),(0,900)), xlabel="Time (kyr)", ylabel="⁸¹K closed system age (kyr)" )
    ax2.xticks = -3000:1000:0
    ax2.yticks = 0:100:900

    for (j,d) in enumerate(k81.depth)
        lines!(ax2,k81.time[1:k_last_pred] .- t_max_val,k81.dat[j,1:k_last_pred],color=col_k81[j],linewidth=2)
    end

    if overlay !== nothing
        b2 = overlay.b
        t2 = overlay.t_max === nothing ? Float64(b2.k81.time[end]) : Float64(overlay.t_max)
        k_last_pred2 = searchsortedlast(b2.k81.time, t2)
        for (j,d) in enumerate(b2.k81.depth)
            lines!(ax2, b2.k81.time[1:k_last_pred2] .- t2, b2.k81.dat[j,1:k_last_pred2];
                   color=col_overlay, linewidth=1.5, linestyle=:dash)
        end
    end

    hlines!(ax2,k81_obs.age;color=col_k81,linewidth=2,linestyle=:dash)
    hspan!(ax2, k81_obs.age .- k81_obs.age_lo, k81_obs.age .+ k81_obs.age_hi; color=col_k81_transparent)

    return fig
end

"""
    save_ensemble_results(path; chain, k81, dar40, depth, setup, priors,
                                sampler_choice, likelihood=:combined)

Persist a sampled ensemble to a JLD2 file. If `path === nothing`, no-op. Pass
all other arguments as keywords. Companion to `load_ensemble_results`.

`likelihood` records which observation channels the chain was conditioned on
(`:combined`, `:kr81`, or `:ar40`). Defaults to `:combined`.
"""
function save_ensemble_results(path; chain, k81, dar40, depth, setup, priors,
                               sampler_choice, likelihood::Symbol=:combined)
    if path === nothing
        println("Chain not saved - no path given.")
        return
    end
    mkpath(dirname(path))
    JLD2.jldsave(path;
                 chain, k81, dar40, depth, setup, priors,
                 sampler_choice, likelihood)
    println("Saved chain to ", path)
    return
end


"""
    load_ensemble_results(path::String) -> NamedTuple

Read a JLD2 snapshot written by `run_basalmixing_ensemble.jl`. Returns a
NamedTuple with `(chain, k81, dar40, depth, setup, priors, sampler_choice,
likelihood)`. Old snapshots that pre-date the `likelihood` field default to
`:combined`.
"""
function load_ensemble_results(path::String)
    isfile(path) || error("ensemble results file not found: $path")
    JLD2.jldopen(path, "r") do f
        chain = f["chain"]
        k81 = f["k81"]
        dar40 = f["dar40"]
        depth = f["depth"]
        setup = haskey(f, "setup") ? f["setup"] : ""
        priors = f["priors"]
        sampler_choice = haskey(f, "sampler_choice") ? f["sampler_choice"] : :unknown
        likelihood = haskey(f, "likelihood") ? f["likelihood"] : :combined
        return (; chain, k81, dar40, depth, setup, priors, sampler_choice, likelihood)
    end
end
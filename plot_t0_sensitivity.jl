## Sensitivity of the inception age t_0 to the other parameters.
##
## t_0 is the headline number but is non-identifiable on its own: the ⁸¹Kr ages
## saturate, so the likelihood is flat in t_0 above a lower bound, and what
## residual constraint exists comes from δ⁴⁰Ar — which trades off with the
## bottom-source flux F_ar40 (accumulated ⁴⁰Ar ≈ F_ar40 × duration). F_ar40 is a
## per-area ⁴⁰Ar flux in cc(STP) m⁻² kyr⁻¹: it enters the ODE as
## `du_ar40[end] += F_ar40 / thickness[end]`, where du_ar40 is the ⁴⁰Ar
## concentration rate in cc m⁻³ kyr⁻¹ and thickness is in m. See `f_ar40_to_mol`
## (BasalMixingModel.jl) to convert to mol m⁻² yr⁻¹ for context. These two
## figures make that explicit:
##
##   1. t0-vs-params : posterior scatter of t_0 vs each parameter (which
##                     parameters does t_0 actually trade off with?).
##   2. t0-vs-far40  : joint-loglik map over (F_ar40, t_0) with the other params
##                     pinned at their posterior medians, posterior samples
##                     overlaid (how sharply would knowing F_ar40 date the onset?).
##
## Usage:  julia --project=. plot_t0_sensitivity.jl [combined_chain.jld2]

cd(@__DIR__); import Pkg; Pkg.activate(".")
using JLD2, DataFrames, CairoMakie, Statistics, Dates, Printf
using Turing   # loads FlexiChains so JLD2 can reconstruct the chain type
include("BasalMixingModel.jl")

const CHAIN_PATH = length(ARGS) >= 1 ? ARGS[1] : "output/combined/chain.jld2"
const OUTDIR     = "plots"

log_line(s) = (println(s); flush(stdout))
log_line("[load] $CHAIN_PATH")
res  = load_ensemble_results(CHAIN_PATH)
k81, dar40, depth = res.k81, res.dar40, res.depth
geom = BasalMixingGeometry(depth, 3040.0, 3053.44, k81.depth, dar40.depth)
lik  = get(res, :likelihood, :combined)
lik_label(l) = l === :combined ? "joint ⁸¹Kr–⁴⁰Ar sampling" :
               l === :kr81    ? "⁸¹Kr-only sampling" :
               l === :ar40    ? "⁴⁰Ar-only sampling" : string(l)

df = DataFrame(res.chain)
df.loglikelihood = vec(res.chain[:loglikelihood])

params = [:delta, :m_clean, :f_dirty, :t_old, :F_ar40]
plabel = Dict(:delta=>"δ transition width (m)", :m_clean=>"m_clean (m/kyr)",
              :f_dirty=>"f_dirty", :t_old=>"t_old (kyr)",
              :F_ar40=>"F_ar40 (cc m⁻² kyr⁻¹)")
med = Dict(p => quantile(df[!, p], 0.5) for p in params)
t0_med = quantile(df.t_0, 0.5)

# ── shared likelihood machinery (used by both figures) ──────────────────────
σk = Float64.(k81.age_sigma);     ok_obs = Float64.(k81.age)
σa = Float64.(dar40.dar40_sigma); da_obs = Float64.(dar40.dar40)
thr68 = -1.15      # ΔlogL, 68% joint region (Δχ²=2.30, 2 dof)
thr95 = -2.996     # ΔlogL, 95% joint region (Δχ²=5.99, 2 dof)

# Joint loglik with one parameter set to `pval` and the rest fixed at their
# posterior medians; t_0 supplied separately (t_target = -t0).
function cond_ll(pname::Symbol, pval, t0)
    v = (delta   = pname===:delta   ? pval : med[:delta],
         m_clean = pname===:m_clean ? pval : med[:m_clean],
         f_dirty = pname===:f_dirty ? pval : med[:f_dirty],
         t_old   = pname===:t_old   ? pval : med[:t_old],
         F_ar40  = pname===:F_ar40  ? pval : med[:F_ar40])
    p = (v..., t_target=-t0)
    age, d40, ok = run_basal_mixing_ode(p, geom; reltol=1e-4, abstol=1e-7)
    ok || return -Inf
    return -0.5*sum(((age .- ok_obs)./σk).^2) - 0.5*sum(((d40 .- da_obs)./σa).^2)
end

# 2-D conditional loglik grid over (param, t_0), normalised to its own max.
function cond_grid(pname, pg, tg)
    LL = Array{Float64}(undef, length(tg), length(pg))
    Threads.@threads for j in eachindex(pg)
        for i in eachindex(tg)
            LL[i, j] = cond_ll(pname, pg[j], tg[i])
        end
    end
    return LL .- maximum(filter(isfinite, LL))
end

# ───────────────────────── Figure 1: t_0 vs each parameter ──────────────────
log_line("[fig1] t_0 vs params: scatter + 95% band + 95% conditional contour")
fig1 = Figure(size=(1000, 620))
Label(fig1[0, 1:3], "Posterior: inception age t₀ vs each parameter  ($(lik_label(lik)))";
      fontsize=16, font=:bold)
crange = (quantile(df.loglikelihood, 0.05), maximum(df.loglikelihood))
tg1 = collect(range(-3000.0, -300.0; length=90))
const R_CUTOFF = 0.5                              # contour only where t₀ truly trades off
rho = Dict(p => cor(df[!, p], df.t_0) for p in params)
params_sorted = sort(params; by = p -> rho[p])    # panels ordered by increasing r
sc1 = nothing
for (i, p) in enumerate(params_sorted)
    global sc1
    row, col = divrem(i-1, 3)
    ax = Axis(fig1[row+1, col+1]; xlabel=plabel[p], ylabel="t₀ (kyr)")
    # 95% posterior CI of this parameter (orange band).
    plo, phi = quantile(df[!, p], [0.025, 0.975])
    vspan!(ax, plo, phi; color=(:orange, 0.12))
    sc1 = scatter!(ax, df[!, p], df.t_0; color=df.loglikelihood,
                   colormap=:viridis, colorrange=crange, markersize=6, alpha=0.5)
    # 95% conditional contour — only where t₀ genuinely trades off (|r| > cutoff);
    # elsewhere it's a compact, uninformative blob, so skip it.
    if abs(rho[p]) > R_CUTOFF
        pg = collect(range(extrema(df[!, p])...; length=36))
        dLLp = cond_grid(p, pg, tg1)
        contour!(ax, pg, tg1, permutedims(dLLp); levels=[thr95], color=:gray25, linewidth=1.2)
    end
    hlines!(ax, [t0_med]; color=(:firebrick, 0.7), linestyle=:dash, linewidth=1.5)
    text!(ax, 0.04, 0.96; text=@sprintf("r = %+.2f", rho[p]), space=:relative,
          align=(:left, :top), fontsize=13, font=:bold,
          color = abs(rho[p]) > R_CUTOFF ? :firebrick : :gray30)
end
Colorbar(fig1[1:2, 4], sc1; label="loglikelihood")
Label(fig1[3, 1:3],
      "Orange band: 95% posterior CI of the parameter. Dashed red: posterior median t₀.\n" *
      "Contour: 95% joint region with the other parameters fixed at their posterior " *
      "medians (when |r| > 0.5). ";
      fontsize=11, color=:gray30, halign=:left)

# ──────────────────── Figure 2: (F_ar40, t_0) joint-loglik map ──────────────
log_line("[fig2] joint-loglik map over (F_ar40, t_0) — others at posterior median")
F_grid  = collect(range(0.0025, 0.0095; length=44))
t0_grid = collect(range(-3000.0, -300.0; length=150))
dLL = cond_grid(:F_ar40, F_grid, t0_grid)        # normalised to its own max
# Conditional ridge: best t_0 for each assumed F_ar40.
ridge_t0 = [t0_grid[argmax(view(dLL, :, j))] for j in eachindex(F_grid)]

fig2 = Figure(size=(900, 640))
ax = Axis(fig2[1, 1]; xlabel="F_ar40  (bottom-source ⁴⁰Ar flux, cc m⁻² kyr⁻¹)",
          ylabel="t₀  (inception age, kyr)",
          title="Inception age t₀ vs ⁴⁰Ar flux  ($(lik_label(lik)))")
# Show only the 95% joint region (Δχ²=5.99, 2 dof → ΔlogL=-3.0); white outside.
dLL_masked = map(x -> x >= thr95 ? x : NaN, dLL)
hm = heatmap!(ax, F_grid, t0_grid, permutedims(dLL_masked);
              colormap=:viridis, colorrange=(thr95, 0))
contour!(ax, F_grid, t0_grid, permutedims(dLL); levels=[thr68],  # 68% (2 dof)
         color=:white, linewidth=1.2)
# Posterior samples: small solid semi-transparent dots. Points landing on the
# white (masked) area are samples that fell outside the 95% region — i.e. drawn
# at low likelihood — which is informative to see directly.
scatter!(ax, df.F_ar40, df.t_0; color=(:black, 0.35), markersize=4)
# Posterior 95% band on F_ar40 and median t_0.
Flo, Fhi = quantile(df.F_ar40, [0.025, 0.975])
vspan!(ax, Flo, Fhi; color=(:orange, 0.12))
hlines!(ax, [t0_med]; color=:orange, linestyle=:dash, linewidth=1.5,
        label=@sprintf("posterior median t₀ = %d kyr", round(Int, t0_med)))
axislegend(ax; position=:rt, framevisible=true, backgroundcolor=(:white, 0.8))
Colorbar(fig2[1, 2], hm; label="Δ loglikelihood (vs best)")
Label(fig2[2, 1],
      "Coloured = 95% joint region; white contour = 68% (others fixed at posterior median). " *
      "Orange band: posterior 95% on F_ar40.";
      fontsize=11, color=:gray30, halign=:left, tellwidth=false)

prefix = joinpath(OUTDIR, string(Dates.today()) * "_")
mkpath(OUTDIR)
save(prefix * "t0-vs-params.png", fig1; px_per_unit=2)
save(prefix * "t0-vs-far40.png", fig2; px_per_unit=2)
log_line("[saved] " * prefix * "t0-vs-params.png")
log_line("[saved] " * prefix * "t0-vs-far40.png")

# Console summary: how tightly would a given F_ar40 precision pin t_0?
log_line("\nConditional best-fit t₀ along the ridge:")
for (Fv, tv) in zip(F_grid[1:6:end], ridge_t0[1:6:end])
    log_line(@sprintf("  F_ar40 = %.4f  ->  t₀ ≈ %d kyr", Fv, round(Int, tv)))
end
log_line("[done]")

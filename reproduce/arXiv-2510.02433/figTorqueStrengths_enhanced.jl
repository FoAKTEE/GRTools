# arXiv:2510.02433 fig:TorqueStrengths + tab:FitParameters WITH the relativistic amplitude
# enhancement — the piece the fast driver `figTorqueStrengths.jl` deliberately leaves out.
#
# The enhancement is E^{ILR/OLR}(m, r₀) = |N_Teukolsky / N_Newtonian|, i.e. the ratio of the FULL
# relativistic Hirata overlap (`AccretionDisc.lindblad_torque`, the Teukolsky ℓ-sum through the
# GrTools Perturbation MST engine) to the Newtonian-amplitude hybrid kernel
# (`AccretionDisc.lindblad_strength`) that the total-torque assembly uses. It is built once on a
# grid at a = 0 and applied to every disc through the `enhancement` hook of
# `total_lindblad_torque` (the Shakura–Sunyaev disc is left UN-enhanced: its Newtonian amplitude is
# already the right one, and applying a relativistic enhancement to a Newtonian disc would
# double-count — that is also why the SS slope is already exact in the fast driver).
#
# TIER heavy — expected wall time at -t4 on an idle machine:
#   * enhancement grid: 9 m-values × 6 radii × 2 resonances = 108 `lindblad_torque` calls at
#     lmax = max(16, m+8), ≈ 3–8 s each  →  ~8–12 min. The ℓ-truncation MATTERS: at lmax = 12 the
#     ratio is 52 % low at m = 10 and 32 % low at m = 8 (measured); it is converged by lmax = 16
#     for m ≤ 10 (E(28)/E(12) = 1.000 / 1.089 / 1.324 / 1.518 at m = 2 / 5 / 8 / 10, r₀ = 30M);
#   * five 50-point torque sweeps  →  ~5 min (measured 4m43s in the fast driver);
#   total ≈ 10–12 min, roughly double under the mission's concurrent load.
# MEASURED RESULT (this driver, converged grid): E < 1 everywhere and FALLING with m
# (0.86 at m=2, 0.67 at m=5, 0.48 at m=8, 0.32 at m=10 at r₀=30M), so the self-consistent
# relativistic amplitude correction FLATTENS the total-torque slope instead of steepening it:
# n_r moves from 3.19/3.31/3.49/3.70 (un-enhanced) further AWAY from the published
# 3.57/3.66/3.77/3.92. The legacy "enhancement grid" was > 1 (1.9–9.6) because it divided by the
# PURE-Newtonian N (fixed ς_m^∓), double-counting the relativistic resonance-geometry shift that
# the hybrid kernel already carries. The published slope steepening therefore is NOT explained by
# the Teukolsky amplitude — an open obligation, recorded rather than tuned.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2510.02433/figTorqueStrengths_enhanced.jl --regen --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, NovikovThorneDisc, nt_isco,
                             lindblad_torque, lindblad_strength, total_lindblad_torque, powerlaw_fit
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const MS_E  = 2:10
const R0S_E = [12.0, 16.0, 22.0, 30.0, 45.0, 70.0]
const SPINS = (0.0, 0.3, 0.6, 0.9)
const NPTS  = 50
const PUB_A  = Dict(0.0 => 2.33, 0.3 => 1.69, 0.6 => 1.16, 0.9 => 0.721, :SS => 0.176)
const PUB_NR = Dict(0.0 => 3.57, 0.3 => 3.66, 0.6 => 3.77, 0.9 => 3.92, :SS => 4.39)

d0 = RelativisticDisc(a = 0.0)
EI = zeros(length(MS_E), length(R0S_E))
EO = zeros(length(MS_E), length(R0S_E))
println("enhancement grid E(m, r₀) = |N_Teukolsky / N_Newtonian| at a = 0 (lmax = max(16, m+8))")
for (i, m) in enumerate(MS_E), (j, r0) in enumerate(R0S_E)
    t0 = time()
    for (mat, kind) in ((EI, :ILR), (EO, :OLR))
        rel = lindblad_torque(d0, r0, m; kind = kind, lmax = max(16, m + 8))
        newt = lindblad_strength(d0, r0, m; kind = kind)
        mat[i, j] = (isfinite(rel) && isfinite(newt) && newt != 0) ? abs(rel/newt) : 1.0
    end
    @printf("  m=%2d r₀=%5.1f  E_ILR=%7.3f  E_OLR=%7.3f  [%.0f s]\n", m, r0, EI[i, j], EO[i, j], time() - t0)
end
@repro_result "E_grid_m" collect(MS_E)
@repro_result "E_grid_r0" R0S_E
@repro_result "E_ILR" EI
@repro_result "E_OLR" EO

"nearest-m, log-linear-in-r interpolation of the enhancement grid (the legacy scheme)"
function Einterp(mat, m, r)
    i = clamp(searchsortedfirst(collect(MS_E), clamp(m, first(MS_E), last(MS_E))), 1, length(MS_E))
    rc = clamp(r, first(R0S_E), last(R0S_E))
    j = clamp(searchsortedlast(R0S_E, rc), 1, length(R0S_E) - 1)
    t = (log(rc) - log(R0S_E[j])) / (log(R0S_E[j+1]) - log(R0S_E[j]))
    max((1 - t)*mat[i, j] + t*mat[i, j+1], 0.0)
end
enh(m, r, kind) = m <= 1 && kind === :ILR ? 1.0 : Einterp(kind === :ILR ? EI : EO, m, r)

curves = Dict{Any,Any}()
for a in SPINS
    nt = NovikovThorneDisc(a = a)
    g = exp10.(range(log10(nt_isco(nt) + 0.6), log10(100.0); length = NPTS))
    t0 = time()
    T = [total_lindblad_torque(nt, rs; enhancement = enh) for rs in g]
    curves[a] = (g = g, T = T)
    @repro_result "rs_grid_a$(a)" g
    @repro_result "TLR_enh_a$(a)" T
    @printf("a=%.1f enhanced sweep done in %.1f s\n", a, time() - t0)
end
let nt = NovikovThorneDisc(a = 0.0)
    g = exp10.(range(log10(6.5), log10(100.0); length = NPTS))
    T = [total_lindblad_torque(nt, rs; profile = :ss) for rs in g]      # SS: NOT enhanced (see header)
    curves[:SS] = (g = g, T = T)
    @repro_result "rs_grid_SS" g
    @repro_result "TLR_SS" T
end

fits = Dict{Any,Any}()
for k in (SPINS..., :SS)
    d = curves[k]
    msk = (d.g .>= 15.0) .& (d.T .> 0)
    fits[k] = powerlaw_fit(d.g[msk], d.T[msk])
end
κU = 2.33e-9 / fits[0.0].A
@repro_result "kappa_U" κU
println("\n=== tab:FitParameters \"All\" rows, WITH enhancement ===")
@printf("%-5s | %10s %10s | %8s %8s\n", "a/M", "A (this)", "A (pub)", "n (this)", "n (pub)")
for k in (SPINS..., :SS)
    A = fits[k].A*κU/1e-9; n = fits[k].n
    @repro_result "fit_A_$(k)" A
    @repro_result "fit_n_$(k)" n
    @printf("%-5s | %10.3f %10.3f | %8.2f %8.2f\n", string(k), A, PUB_A[k], n, PUB_NR[k])
end
@repro_result "slope_gap_to_published" [fits[k].n - PUB_NR[k] for k in (SPINS..., :SS)]

emri_defaults!()
COL = Dict(0.0 => RGB(0.30, 0.55, 0.95), 0.3 => RGB(0.15, 0.55, 0.20),
           0.6 => RGB(0.55, 0.20, 0.55), 0.9 => RGB(0.95, 0.65, 0.10), :SS => RGB(0.85, 0.15, 0.15))
plt = plot(; xlabel = "\$r_s/M\$", ylabel = "\$|T_\\mathrm{LR}|/q^2\$", xscale = :log10, yscale = :log10,
           legend = :bottomright, xlims = (5.5, 105), ylims = (1e-13, 1e-5), size = (720, 560))
for k in (SPINS..., :SS)
    d = curves[k]; pos = d.T .> 0; Tp = abs.(d.T) .* κU
    lab = k === :SS ? "Shakura - Sunyaev" : "\$a/M=$(k)\$"
    plot!(plt, d.g[pos], Tp[pos]; color = COL[k], lw = 2.2, ls = :solid, label = lab)
    plot!(plt, d.g[.!pos], Tp[.!pos]; color = COL[k], lw = 2.2, ls = :dot, label = "")
    plot!(plt, d.g, fits[k].A*κU .* (d.g ./ 10) .^ fits[k].n; color = COL[k], ls = :dash, lw = 1.3, label = "")
end
png = joinpath(@__DIR__, "expected", "figTorqueStrengths_enhanced.png")
savefig(plt, png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

exit(ReproHarness.conclude!(opts; paper = "2510.02433", id = "figTorqueStrengths_enhanced",
                            script = @__FILE__, rtol = 1e-8))

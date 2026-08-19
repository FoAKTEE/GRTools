# arXiv:2510.02433 fig:TorqueStrengths (left, "all r") + tab:FitParameters (the "All" rows) —
# the total signed Lindblad torque |T_LR|/q² vs r_s for Novikov–Thorne discs at a/M = 0, 0.3, 0.6,
# 0.9 and for the Shakura–Sunyaev comparison disc, and the single power-law fits T = A (r/10M)^{n_r}
# on the positive (OLR-dominated) branch.
#
# THIN driver: `AccretionDisc.total_lindblad_torque` (Σ_m Σ_{ILR,OLR} 2π r_LR Σ(r_LR) N^{(m)}
# × pressure_cutoff, with the signed hybrid strength `lindblad_strength` and the Novikov–Thorne
# Σ, h) and `AccretionDisc.powerlaw_fit`.
#
# SCOPE / FINDING — this driver runs the kernel WITHOUT the relativistic amplitude enhancement
# E(m, r) = |N_Teukolsky/N_Newtonian|: `lindblad_strength` uses Hirata's Newtonian-limit overlap at
# the relativistic resonance geometry. The published slopes (n_r = 3.57, 3.66, 3.77, 3.92, 4.39) are
# steeper than the un-enhanced ones registered here; the legacy analysis closed roughly half of that
# gap with an a=0 enhancement grid (see figTorqueStrengths_enhanced.jl, heavy tier). The amplitude
# A is not predicted at all in absolute units — the paper's 10⁻⁹/q² axis is fixed by ONE global unit
# factor κ_U, here anchored (as in the legacy analysis) so that the a=0 all-r amplitude equals the
# published 2.33e-9. Only the SHAPE (slopes, spin ordering, sign reversal) is a prediction.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2510.02433/figTorqueStrengths.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: NovikovThorneDisc, nt_isco, total_lindblad_torque, powerlaw_fit
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const SPINS = (0.0, 0.3, 0.6, 0.9)
const NPTS = 50
const PUB_A  = Dict(0.0 => 2.33, 0.3 => 1.69, 0.6 => 1.16, 0.9 => 0.721, :SS => 0.176)   # 1e-9
const PUB_NR = Dict(0.0 => 3.57, 0.3 => 3.66, 0.6 => 3.77, 0.9 => 3.92, :SS => 4.39)

"last sign change of T_LR below 30M (the ILR→OLR reversal radius)"
function reversal_radius(g, T)
    ir = findlast(i -> T[i]*T[i+1] < 0, 1:length(T)-1)
    ir === nothing ? NaN : 0.5*(g[ir] + g[ir+1])
end

curves = Dict{Any,Any}()
for a in SPINS
    nt = NovikovThorneDisc(a = a)
    g = exp10.(range(log10(nt_isco(nt) + 0.6), log10(100.0); length = NPTS))
    t0 = time()
    T = [total_lindblad_torque(nt, rs) for rs in g]
    curves[a] = (g = g, T = T, rrev = reversal_radius(g, T))
    @printf("a=%.1f done in %5.1f s (reversal at r_s = %.3f M)\n", a, time() - t0, curves[a].rrev)
    @repro_result "rs_grid_a$(a)" g
    @repro_result "TLR_a$(a)" T
    @repro_result "reversal_radius_a$(a)" curves[a].rrev
end
let nt = NovikovThorneDisc(a = 0.0)
    g = exp10.(range(log10(6.5), log10(100.0); length = NPTS))
    T = [total_lindblad_torque(nt, rs; profile = :ss) for rs in g]
    curves[:SS] = (g = g, T = T, rrev = reversal_radius(g, T))
    @printf("SS   done (reversal at r_s = %.3f M)\n", curves[:SS].rrev)
    @repro_result "rs_grid_SS" g
    @repro_result "TLR_SS" T
    @repro_result "reversal_radius_SS" curves[:SS].rrev
end

# power-law fits on the positive branch above the reversal (or 15M, whichever is larger)
fits = Dict{Any,Any}()
for k in (SPINS..., :SS)
    d = curves[k]
    rlo = max(isnan(d.rrev) ? 0.0 : d.rrev + 0.3, 15.0)
    msk = (d.g .>= rlo) .& (d.T .> 0)
    fits[k] = powerlaw_fit(d.g[msk], d.T[msk])
end
κU = 2.33e-9 / fits[0.0].A          # ONE global unit factor, anchored to the published a=0 amplitude
@repro_result "kappa_U" κU
@repro_result "raw_A_a0" fits[0.0].A

println("\n=== tab:FitParameters, \"All\" rows:  T = A (r/10M)^{n_r}   [A in 1e-9 after κ_U] ===")
@printf("%-5s | %10s %10s | %8s %8s\n", "a/M", "A (this)", "A (pub)", "n (this)", "n (pub)")
for k in (SPINS..., :SS)
    A = fits[k].A*κU/1e-9; n = fits[k].n
    @repro_result "fit_A_$(k)" A
    @repro_result "fit_n_$(k)" n
    @printf("%-5s | %10.3f %10.3f | %8.2f %8.2f\n", string(k), A, PUB_A[k], n, PUB_NR[k])
end
@repro_result "slope_gap_to_published" [fits[k].n - PUB_NR[k] for k in (SPINS..., :SS)]
@printf("slope gap n(this) − n(published): %s\n", string(round.([fits[k].n - PUB_NR[k] for k in (SPINS..., :SS)], digits = 3)))

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
png = joinpath(@__DIR__, "expected", "figTorqueStrengths.png")
savefig(plt, png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

exit(ReproHarness.conclude!(opts; paper = "2510.02433", id = "figTorqueStrengths",
                            script = @__FILE__, rtol = 1e-8))

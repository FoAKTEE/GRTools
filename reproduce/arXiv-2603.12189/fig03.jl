# arXiv:2603.12189 fig:fluxcomparen (fig. 3, main.tex:834): the EXACT linear-in-spin flux
# corrections |δF^{E,∞}|, |δF^{E,H}|, |δF^{Lz,∞}|, |δF^{Lz,H}| of the (l,m,k) = (2,2,0) modes
# versus the radial harmonic n = 1..11, against the SHIFTED-GEODESIC (SG) approximation, for
# a = 0.9 M, p = 12, e = 0.2, x = √3/2.
#
# THIN driver, both curves from Perturbation.SpinFluxes:
#   EXACT = TeukolskySpinModeCorrection(2, 2, n, 0, oc).FluxesCorrection   (full numeric trajectory
#           with the oscillatory δχ / δr libration — the exact O(σ) correction)
#   SG    = central finite difference in the spin parameter of the PEX-GAUGE full flux,
#           TeukolskySpinModeAnalytical(...; gauge = "pex"), step s = 2e-6 — the secular shift only.
# The CTilde-locked correction entry point is a DIFFERENT gauge here and lands ~7.7× off with the
# opposite sign, so it is not the SG curve.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2603.12189/fig03.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: KerrSpinOrbitCorrection, TeukolskySpinModeCorrection,
                            TeukolskySpinModeAnalytical
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, P, E, X = 0.9, 12.0, 0.2, sqrt(3) / 2
const NS = collect(1:11)
const SFD = 2e-6
const CHANNELS = ("Einf", "EH", "Lzinf", "LzH")

oc = KerrSpinOrbitCorrection(A, P, E, X, 8, 8)

EX = [zeros(length(NS)) for _ in 1:4]
SG = [zeros(length(NS)) for _ in 1:4]
Threads.@threads for i in eachindex(NS)
    n = NS[i]
    ex = TeukolskySpinModeCorrection(2, 2, n, 0, oc)
    EX[1][i] = abs(ex.FluxesCorrection.Energy.I)
    EX[2][i] = abs(ex.FluxesCorrection.Energy.H)
    EX[3][i] = abs(ex.FluxesCorrection.AngularMomentum.I)
    EX[4][i] = abs(ex.FluxesCorrection.AngularMomentum.H)
    fp = TeukolskySpinModeAnalytical(2, 2, n, 0, A, P, E, X,  SFD; gauge = "pex")
    fm = TeukolskySpinModeAnalytical(2, 2, n, 0, A, P, E, X, -SFD; gauge = "pex")
    SG[1][i] = abs((fp.Fluxes.Energy.I          - fm.Fluxes.Energy.I)          / (2 * SFD))
    SG[2][i] = abs((fp.Fluxes.Energy.H          - fm.Fluxes.Energy.H)          / (2 * SFD))
    SG[3][i] = abs((fp.Fluxes.AngularMomentum.I - fm.Fluxes.AngularMomentum.I) / (2 * SFD))
    SG[4][i] = abs((fp.Fluxes.AngularMomentum.H - fm.Fluxes.AngularMomentum.H) / (2 * SFD))
end

for (j, nm) in enumerate(CHANNELS)
    println("  ", rpad(nm, 6), " exact(n=1) = ", EX[j][1], "  SG(n=1) = ", SG[j][1],
            "   exact(n=11) = ", EX[j][end], "  SG(n=11) = ", SG[j][end])
    flush(stdout)
end

ylabs = (raw"$|\delta\mathcal{F}^{E,\infty}_{lmnk}|$", raw"$|\delta\mathcal{F}^{E,H}_{lmnk}|$",
         raw"$|\delta\mathcal{F}^{L_z,\infty}_{lmnk}|$", raw"$|\delta\mathcal{F}^{L_z,H}_{lmnk}|$")
panels = Any[]
for j in 1:4
    pl = plot(xlabel = "n", ylabel = ylabs[j], xscale = :log10, yscale = :log10,
              legend = (j == 1 ? :bottomleft : false),
              xticks = ([1, 2, 5, 10], ["1", "2", "5", "10"]))
    plot!(pl, NS, EX[j], marker = (:circle, 4), line = (:dash, 2),    color = MMA_RED,  label = "Exact")
    plot!(pl, NS, SG[j], marker = (:rect, 4),   line = (:dashdot, 2), color = MMA_GOLD, label = "Shifted geodesic")
    push!(panels, pl)
end
savefig(plot(panels..., layout = (2, 2), size = (1100, 750)), joinpath(@__DIR__, "fig03.png"))

@repro_result "n" NS
for (j, nm) in enumerate(CHANNELS)
    @repro_result "exact_$(nm)" EX[j]
    @repro_result "sg_$(nm)"    SG[j]
end
# the published message: SG tracks the exact correction at low n and deviates at high n
@repro_result "rel_SG_vs_exact_n1"  [abs(1 - SG[j][1] / EX[j][1]) for j in 1:4]
@repro_result "rel_SG_vs_exact_n11" [abs(1 - SG[j][end] / EX[j][end]) for j in 1:4]
@repro_result "SG_closer_at_low_n"  all(j -> abs(1 - SG[j][1] / EX[j][1]) <
                                            abs(1 - SG[j][end] / EX[j][end]), 1:4)

exit(ReproHarness.conclude!(opts; paper = "2603.12189", id = "fig03", script = @__FILE__,
                            rtol = 1e-8))

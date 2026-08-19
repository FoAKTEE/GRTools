# arXiv:2603.13482 fig:Clmnk_k (main.tex:671), the "Analytical" (green) curves of both panels:
#   FE1_220k     |F^{E,1}_{220k}| — the linear-in-spin energy flux vs k
#   Cplus1_220k  |C^{+,1}_{220k}| — the linear-in-spin infinity amplitude vs k
# for (l,m,n) = (2,2,0), k = −10..10, at a = 0.9 M, p = 12, e = 0.2, x = 0.5.
#
# The published panels ALSO overlay Drummond–Hughes numeric-trajectory curves at k^s_max = 8 and 16;
# that solver is external and un-vendored, so those overlays are [FUTURE] and are NOT drawn here.
#
# THIN driver: the analytical curve is the central finite difference in the spin parameter of the
# PEX-gauge `Perturbation.TeukolskySpinModeAnalytical` (the published "Analytical" convention; the
# CTilde-locked entry point is a different gauge and lands ~8–25× off). Step s_fd = 2e-6.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2603.13482/fig_FE1_Cplus1_220k.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: TeukolskySpinModeAnalytical
using Plots
include(joinpath(@__DIR__, "..", "_style", "mma_style.jl"))
gr(); mma_defaults!()

const A, P, E, X = 0.9, 12.0, 0.2, 0.5
const KS = -10:10
const SFD = 2e-6

Flux = fill(NaN, length(KS))
Camp = fill(NaN, length(KS))
Threads.@threads for i in eachindex(KS)
    k = KS[i]
    try
        fp = TeukolskySpinModeAnalytical(2, 2, 0, k, A, P, E, X,  SFD; gauge = "pex", Nmax = 40)
        fm = TeukolskySpinModeAnalytical(2, 2, 0, k, A, P, E, X, -SFD; gauge = "pex", Nmax = 40)
        Flux[i] = abs((fp.Fluxes.Energy.I + fp.Fluxes.Energy.H -
                       fm.Fluxes.Energy.I - fm.Fluxes.Energy.H) / (2 * SFD))
        Camp[i] = abs((fp.Amplitudes.I - fm.Amplitudes.I) / (2 * SFD))
    catch                                   # continued-fraction limit at large |k|
    end
end

fin = isfinite.(Flux)
println("  FE1    peak k=", collect(KS)[argmax(replace(Flux, NaN => -Inf))],
        " = ", maximum(Flux[fin]), "   @k=+10 ", Flux[end], "   @k=−10 ", Flux[1])
println("  Cplus1 peak k=", collect(KS)[argmax(replace(Camp, NaN => -Inf))],
        " = ", maximum(Camp[isfinite.(Camp)]), "   @k=+10 ", Camp[end], "   @k=−10 ", Camp[1])

pf = plot(xlabel = "k", ylabel = raw"$|\mathcal{F}^{E,1}_{220k}|$", yscale = :log10,
          legend = :topright, xticks = -10:5:10)
plot!(pf, collect(KS), Flux, marker = :diamond, ms = 5, color = MMA_GREEN, label = "Analytical")
savefig(pf, joinpath(@__DIR__, "FE1_220k.png"))

pc = plot(xlabel = "k", ylabel = raw"$|C^{+,1}_{220k}|$", yscale = :log10,
          legend = :topright, xticks = -10:5:10)
plot!(pc, collect(KS), Camp, marker = :diamond, ms = 5, color = MMA_GREEN, label = "Analytical")
savefig(pc, joinpath(@__DIR__, "Cplus1_220k.png"))

@repro_result "k"            collect(KS)
@repro_result "absFE1"       Flux
@repro_result "absCplus1"    Camp
@repro_result "FE1_at_k0"    Flux[findfirst(==(0), KS)]
@repro_result "Cplus1_at_k0" Camp[findfirst(==(0), KS)]
@repro_result "FE1_peak_k"   collect(KS)[argmax(replace(Flux, NaN => -Inf))]
# the published message: both curves peak near k = 0 and fall by many decades toward |k| = 10
@repro_result "FE1_decades_k0_to_k10"    log10(Flux[findfirst(==(0), KS)] / Flux[end])
@repro_result "Cplus1_decades_k0_to_k10" log10(Camp[findfirst(==(0), KS)] / Camp[end])

exit(ReproHarness.conclude!(opts; paper = "2603.13482", id = "fig_FE1_Cplus1_220k",
                            script = @__FILE__, rtol = 1e-8))

# arXiv:1010.0759 fig:kerrm1 (LEFT panel) — Kerr m=1 inner-Lindblad resonance LOCATIONS r₁(r₀) for
# seven spins a/M ∈ {−0.9,…,0.9} (prograde perturber; a<0 = retrograde).
#
# THIN driver: `AccretionDisc.lindblad_locations` on the Kerr equatorial geodesic frequencies
# (identical to `Spacetime.kerr_frequencies(a,r,0,1).Omegaph/Omegar`), `Spacetime.kerr_isco` for the
# inner edge. Published anchors (l.1100): r₀ = 250 ⇒ r₁ = 39.8 (a = +0.9), 46.4 (a = −0.9); the a = 0
# cross-check is tab:ilrS r₁ = 43.25.
#
# The RIGHT (strength) panel needs the Kerr relativistic overlap at every spin; the two published
# anchor spins are registered here through `lindblad_torque`, the full panel is the heavy driver
# figKerrm1Strengths.jl.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0759/figKerrm1.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_locations
using GrTools.Spacetime: kerr_isco
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
using Plots

const SPINS = [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9]
r0s = exp10.(range(log10(10.0), log10(250.0); length = 34))
curves = Dict(a => [first(lindblad_locations(RelativisticDisc(a = a), r0, 1)) for r0 in r0s] for a in SPINS)

@repro_result "r0_grid" r0s
for a in SPINS
    @repro_result "r1_m1_ILR_a$(a)" curves[a]
end
@repro_result "isco_vs_spin" [kerr_isco(a) for a in SPINS]

# published anchors
anch = (aplus = first(lindblad_locations(RelativisticDisc(a = 0.9), 250.0, 1)),
        aminus = first(lindblad_locations(RelativisticDisc(a = -0.9), 250.0, 1)),
        azero = first(lindblad_locations(RelativisticDisc(a = 0.0), 250.0, 1)))
@repro_result "r1_at_r0_250_a_plus_0.9" anch.aplus
@repro_result "r1_at_r0_250_a_minus_0.9" anch.aminus
@repro_result "r1_at_r0_250_a_0" anch.azero
@printf("anchors at r₀=250: a=+0.9 r₁=%.3f (paper 39.8, rel %.3f) | a=−0.9 r₁=%.3f (paper 46.4, rel %.3f) | a=0 r₁=%.3f (tab:ilrS 43.25, rel %.4f)\n",
        anch.aplus, abs(anch.aplus-39.8)/39.8, anch.aminus, abs(anch.aminus-46.4)/46.4,
        anch.azero, abs(anch.azero-43.25)/43.25)

emri_defaults!()
plt = plot(xscale = :log10, yscale = :log10, xlabel = "r_0 (secondary)",
           ylabel = "r_1 (m=1 ILR location)", title = "Kerr m=1 ILR locations", legend = :topleft)
for a in SPINS
    plot!(plt, r0s, curves[a]; label = "a=$a")
end
png = joinpath(@__DIR__, "expected", "figKerrm1.png")
savefig(plt, png)
@printf("wrote %s (%d bytes)\n", relpath(png, @__DIR__), filesize(png))

exit(ReproHarness.conclude!(opts; paper = "1010.0759", id = "figKerrm1",
                            script = @__FILE__, rtol = 1e-9))

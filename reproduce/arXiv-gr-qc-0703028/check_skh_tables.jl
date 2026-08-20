# Sundararajan, Khanna & Hughes [arXiv:gr-qc/0703028] — the PUBLISHED frequency-domain circular-orbit
# energy-flux rows that anchor the Kerr sourced time-domain sector of `Perturbation.TimeDomain`
# (the DG-Teukolsky s = −2 Kerr flux campaign of emri-disk `results/dg-teukolsky-spin/`, whose
# `data-flux/dgs2_skh_anchor.out` recorded exactly these numbers).
#
# WHAT IS REPRODUCED — the |m|-mode energy fluxes to infinity of a circular equatorial orbit about
# Kerr, summed over ℓ and doubled for ±m:
#   Table I  (column FD, |m| = 2): (a = 0.9, r₀ = 10M) Ė = 4.456e-5 ; (a = 0.5, r₀ = 10M) 4.812e-5
#   Table II (column FD, |m| = 3): (a = 0.9, r₀ = 10M) Ė = 5.044e-6
# from the frequency-domain point-particle Teukolsky modes (`TeukolskyPointParticleModeCircular`).
# The published values carry four digits; the driver stores the full computed sums, the published
# values and the relative differences, so the comparison is visible rather than asserted.
#
# WHY THIS PAPER DIRECTORY: the DG-Teukolsky "F1.7-analogue" and "T1.3-analogue" artefacts of the
# legacy campaign are NOT reproductions of a published object — the DG papers publish no s = −2
# figure or table (recorded finding of that campaign); their reference is this SKH table plus the
# in-package frequency-domain oracle. The published object is therefore the SKH table, and it lives
# here. The time-domain side of the same comparison is registered by `check_kerr_td_flux.jl`.
#
#   julia --project=reproduce -t4 reproduce/arXiv-gr-qc-0703028/check_skh_tables.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Perturbation: TeukolskyPointParticleModeCircular

"Ė^∞ of the |m| pair, summed over the ℓ band (the ±m modes contribute equally)."
mode_sum(m, a, r0, ls) = 2 * sum(TeukolskyPointParticleModeCircular(-2, l, m, a, r0).EnergyFlux.I
                                 for l in ls)

# published rows: (label, m, a, r0, ℓ band, published Ė)
const ROWS = (("tabI_a0p9_r10_m2", 2, 0.9, 10.0, 2:8, 4.456e-5),
              ("tabI_a0p5_r10_m2", 2, 0.5, 10.0, 2:8, 4.812e-5),
              ("tabII_a0p9_r10_m3", 3, 0.9, 10.0, 3:9, 5.044e-6))

for (name, m, a, r0, ls, published) in ROWS
    t0 = time()
    s = mode_sum(m, a, r0, ls)
    @repro_result "$(name)_computed" s
    @repro_result "$(name)_published" published
    @repro_result "$(name)_reldiff" s / published - 1
    println("[repro] $(name): computed $(s)  published $(published)  rel $(s/published - 1)  ",
            "($(round(time()-t0, digits=1)) s)"); flush(stdout)
end

# per-ℓ breakdown of the headline row (the ℓ-tail that bounds the band-truncation error)
@repro_result "tabI_a0p9_r10_m2_per_l" [2 * TeukolskyPointParticleModeCircular(-2, l, 2, 0.9, 10.0).EnergyFlux.I
                                        for l in 2:10]
# horizon flux of the same row (negative: the hole gains energy for a prograde orbit at r₀ = 10M)
@repro_result "tabI_a0p9_r10_m2_horizon" 2 * sum(TeukolskyPointParticleModeCircular(-2, l, 2, 0.9, 10.0).EnergyFlux.H
                                                 for l in 2:8)

exit(ReproHarness.conclude!(opts; paper = "gr-qc-0703028", id = "skh_tables",
                            script = @__FILE__, rtol = 1e-8))

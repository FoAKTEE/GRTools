# arXiv:1301.1591 Tab. m_neq_0_rates_sm2 (ID1, m = 2) analogue — the non-axisymmetric rows of the
# s = −2, a = 0.8 mode-scattering decay-rate matrix, seeds ℓ′ = 2,3.
#
# MEASURED, stated precisely: on the committed P4 window (700, 1300) NO entry of this block reaches
# its printed rate — ℓ′ = 3 → ℓ = 2 at scri fits −4.665 (raw) / −4.674 (envelope) against the
# printed 6, and ℓ′ = 2 → ℓ = 3 at scri fits −3.741 / −4.615 against 7. What the legacy campaign
# reported as "≈ −6.0, dev 0.01" for the ℓ′ = 3 → ℓ = 2 scri entry is the MINIMUM over the
# sliding-window scan (−6.010), not the primary fit; the golden stores raw, envelope, adaptive,
# plateau, min and max side by side so the distinction is visible. Everything here is junk-floor
# limited exactly as in the m = 0 driver, and the paper itself marks its deep entries × even at
# quadruple precision.
#
# tier HEAVY: two tf = 3000 runs, ≈ 16 min at -t4 (legacy 454 s each).
#
#   julia --project=reproduce -t4 reproduce/arXiv-1301.1591/scatter_m2_ID1.jl --check --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

include(joinpath(@__DIR__, "common_scatter.jl"))

for lp in (2, 3)
    t0 = time()
    τs, Ψ, maxu, dt = run_scatter(lp, 2, :ID1)
    println("[repro] lp=$(lp) m=2 ID1 wall=$(round(time()-t0, digits=1)) s maxu_end=$(maxu[end])")
    flush(stdout)
    register_scatter!(lp, 2, :ID1, τs, Ψ, maxu, dt, T_ID1_M2)
end

exit(ReproHarness.conclude!(opts; paper = "1301.1591", id = "scatter_m2_ID1",
                            script = @__FILE__, rtol = 1e-8))

# arXiv:1301.1591 Tab. s2_LPI (ID0, m = 0) analogue — the s = −2, a = 0.8 mode-scattering
# decay-rate matrix, seeds ℓ′ = 2..5 (the DG artefact `T-S1.1_ScatteringMatrix_s-2_a0.8`, m = 0 rows).
#
# MEASURED VERDICT, stated as the legacy verifier table records it and not more: at a = 0.8 the junk
# floor sits high enough that only ONE of the sixteen registered entries falls inside the committed
# ±1.0 gate — ℓ′ = 2 → ℓ = 2 at scri, raw fit −5.147 against the printed 6 (dev 0.853). The rest
# deviate by 1.0 to 7.0, monotonically worse the deeper the printed rate, and every entry with a
# printed rate μ ≥ 8 is simply below the floor: the paper itself marks its μ ≥ 10 entries × even at
# QUADRUPLE precision. The golden stores the measured raw / envelope / adaptive / plateau / min /
# max fits AND the printed laws side by side, so the mismatch is a recorded measurement, never a
# claim of agreement.
#
# tier HEAVY: four tf = 3000 runs with a 5-mode coupled window, ≈ 32 min at -t4 (legacy 455 s each).
#
#   julia --project=reproduce -t4 reproduce/arXiv-1301.1591/scatter_m0_ID0.jl --check --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

include(joinpath(@__DIR__, "common_scatter.jl"))

for lp in (2, 3, 4, 5)
    t0 = time()
    τs, Ψ, maxu, dt = run_scatter(lp, 0, :ID0)
    println("[repro] lp=$(lp) m=0 ID0 wall=$(round(time()-t0, digits=1)) s maxu_end=$(maxu[end])")
    flush(stdout)
    register_scatter!(lp, 0, :ID0, τs, Ψ, maxu, dt, T_ID0_M0)
end

exit(ReproHarness.conclude!(opts; paper = "1301.1591", id = "scatter_m0_ID0",
                            script = @__FILE__, rtol = 1e-8))

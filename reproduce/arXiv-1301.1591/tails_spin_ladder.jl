# arXiv:1301.1591 — Campaign-A spin ladder of the vacuum s = −2 Kerr tail (the DG figure-set
# analogue `F-S1.4n_Tails_spin_ladder`): the SAME committed configuration as `tails_a0p9.jl`
# swept over a ∈ {0.99, 0.995, 0.9995} at K = 24, which documents where the Float64 RT-chart
# sector stops resolving the paper's empirical laws Eqs. (28),(29).
#
# MEASURED BOUNDARY of the legacy campaign, reproduced here rather than asserted: the finite-radius
# probes (R₊, R_mid) track the ℓ = 2 law −7 up to a = 0.995; at a = 0.9995 the horizon probe is
# junk-flooded (fit ≈ 0, i.e. a flat record) while scri still measures ≈ −6. No knob is changed
# between spins — the ceiling is a property of the sector at Float64, not of a tuning.
#
# tier HEAVY: three tf = 3000 runs, ≈ 12 min at -t4 (legacy wall times 166 s / 242 s / 280 s plus
# compilation).
#
#   julia --project=reproduce -t4 reproduce/arXiv-1301.1591/tails_spin_ladder.jl --check --tier heavy

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

include(joinpath(@__DIR__, "common_tails.jl"))

const SPINS = (0.99, 0.995, 0.9995)
const K = 24

for a in SPINS
    t0 = time()
    τs, Ψ, maxu, dt, excise = run_tail(a, K)
    tag = spin_tag(a)
    println("[repro] a=$(a) K=$(K) excise=$(round(excise, digits=5)) dt=$(dt) wall=$(round(time()-t0, digits=1)) s")
    @repro_result "$(tag)_excise" excise
    register_tail_fits!(tag, τs, Ψ, maxu, dt; series_stride = 20, series_ells = (2,))
    flush(stdout)
end

exit(ReproHarness.conclude!(opts; paper = "1301.1591", id = "tails_spin_ladder",
                            script = @__FILE__, rtol = 1e-8))

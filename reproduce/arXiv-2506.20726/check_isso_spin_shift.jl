# arXiv:2506.20726, printed result T1: the secondary-spin shift of the ISSO at the Schwarzschild
# limit is δr_ISSO/s∥ → −√(8/3) ≈ −1.63299. Reproduced as the a → 0 limit of the CTilde-gauge δp of
# `Spacetime.SpinningSecondary.spin_linear_orbital_params` evaluated at the ISSO radius
# `Spacetime.kerr_isso` (which is 6 at a = 0 for every inclination).
#
#   julia --project=reproduce -t4 reproduce/arXiv-2506.20726/check_isso_spin_shift.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using GrTools
using GrTools.Spacetime: kerr_isso, kerr_isco
using GrTools.Spacetime.SpinningSecondary: spin_linear_orbital_params

dr_isso(a, x) = spin_linear_orbital_params(a, kerr_isso(a, x), 0.0, x)[1]

const AS = [1e-3, 1e-4, 1e-5, 1e-6]
shifts_eq   = [dr_isso(a, 1.0) for a in AS]     # equatorial approach to Schwarzschild
shifts_incl = [dr_isso(a, 0.5) for a in AS]     # inclined approach — same limit
target = -sqrt(8 / 3)

println("  −√(8/3) = ", target)
for (a, s, si) in zip(AS, shifts_eq, shifts_incl)
    println("  a=", a, "  δr_ISSO/s∥ (x=1) = ", s, "   (x=0.5) = ", si)
end

@repro_result "spins"                 AS
@repro_result "dr_isso_equatorial"    shifts_eq
@repro_result "dr_isso_inclined_x0.5" shifts_incl
@repro_result "schwarzschild_limit"   target
@repro_result "isso_at_a0_is_six"     isapprox(kerr_isso(1e-12, 0.5), 6.0; rtol = 1e-9)
@repro_result "isco_isso_agree_equatorial" isapprox(kerr_isso(0.9, 1.0), kerr_isco(0.9); rtol = 1e-10)
@repro_result "rel_to_minus_sqrt8over3_at_a1e-6" abs(shifts_eq[end] / target - 1)
@repro_result "limit_reached_to_1e-5" abs(shifts_eq[end] / target - 1) < 1e-5 &&
                                      abs(shifts_incl[end] / target - 1) < 1e-5

exit(ReproHarness.conclude!(opts; paper = "2506.20726", id = "check_isso_spin_shift",
                            script = @__FILE__, rtol = 1e-9))

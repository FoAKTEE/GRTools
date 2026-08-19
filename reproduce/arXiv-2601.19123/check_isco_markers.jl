# arXiv:2601.19123 F7/F9 ISCO markers — the Bardeen–Press–Teukolsky (1972) equatorial circular ISCO
# radius r_ISCO(a)/M used as the vertical marker of fig:densitiesspindep and
# fig:AdvAngSpinVariation, for the marker spins a = 0.01, 0.6, 0.99, plus the textbook cross-checks.
#
# THIN driver: `Spacetime.kerr_isco` (the CANONICAL ISCO of GrTools; `AccretionDisc.risco_disc` and
# the SpiralWaves `r_isco` both delegate to it).
#
#   julia --project=reproduce -t4 reproduce/arXiv-2601.19123/check_isco_markers.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.Spacetime: kerr_isco

println("=== marker spins (prograde) ===")
for a in (0.01, 0.6, 0.99)
    @repro_result "isco_a$(a)" kerr_isco(a)
    @printf("  a=%.2f : r_isco = %.6f M\n", a, kerr_isco(a))
end
println("\n=== textbook cross-check (Bardeen 1972, prograde) ===")
gold = [(0.0, 6.0), (0.5, 4.23300), (0.9, 2.32088), (0.99, 1.45450), (0.998, 1.23708)]
errs = Float64[]
for (a, g) in gold
    r = kerr_isco(a); push!(errs, abs(r - g)/g)
    @printf("  a=%.3f : r_isco=%.6f (textbook %.5f) relerr=%.1e\n", a, r, g, errs[end])
end
@repro_result "textbook_max_relerr" maximum(errs)
@repro_result "isco_extremal_limit" kerr_isco(0.9999999)
println("\n=== retrograde cross-check ===")
retro = [(0.0, 6.0), (-0.9, 8.71735), (-0.99, 8.97175)]
rerrs = Float64[]
for (a, g) in retro
    r = kerr_isco(a); push!(rerrs, abs(r - g)/g)
    @repro_result "isco_retro_a$(a)" r
    @printf("  a=%+.2f : r_isco=%.6f (textbook %.5f) relerr=%.1e\n", a, r, g, rerrs[end])
end
@repro_result "retrograde_max_relerr" maximum(rerrs)

exit(ReproHarness.conclude!(opts; paper = "2601.19123", id = "check_isco_markers",
                            script = @__FILE__, rtol = 1e-10))

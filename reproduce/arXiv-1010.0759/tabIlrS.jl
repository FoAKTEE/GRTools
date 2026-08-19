# arXiv:1010.0759 tab:ilrS — Schwarzschild Lindblad resonances of a circular perturber at r₀:
# the resonance RADII r₁ (m=1 ILR, m=2 ILR, m=2 OLR, m=1 OLR) for the 16 published r₀ rows, and the
# ILR resonant STRENGTHS N (m = 1, 2) from the full relativistic Teukolsky overlap.
#
# THIN driver: `AccretionDisc.lindblad_locations` solves m(Ω−Ω₀) = ±κ on the Kerr equatorial
# geodesic frequencies, `AccretionDisc.lindblad_torque` assembles the Hirata GW overlap
# (eq:P-ILR + eq:SP + eq:Tnorm) through the GrTools Perturbation MST engine.
#
# SCOPE — the OLR strengths are NOT registered. The OLR overlap eq:P-OLR is a near-cancellation of
# two large terms (Z^out_pert(Z^out*_test − c₁₃*Z^down*_test)); in Float64 the legacy run reproduces
# the published OLR column only to factors 1.5–260 (results/emri-disk/figures/1010.0759/
# repro_tab_ilrS.out), a documented BigFloat-precision gap, so registering those numbers would
# register noise. The ILR column and all four location columns are the reproduced content.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0759/tabIlrS.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, lindblad_locations, lindblad_torque

# published tab:ilrS: r₀ => (r₁, N) for m1 ILR, m2 ILR, m2 OLR, m1 OLR
const PUB = [
 (8.0,   6.48,-7.11e-1,  6.97,-4.69e1,   9.55, 2.28e1,  11.33, 7.05e-2),
 (9.0,   6.80,-2.32e-1,  7.57,-2.55e1,  10.91, 1.79e1,  12.99, 8.76e-2),
 (10.0,  7.13,-1.10e-1,  8.19,-1.74e1,  12.26, 1.51e1,  14.62, 1.01e-1),
 (12.0,  7.78,-4.08e-2,  9.44,-1.09e1,  14.92, 1.30e1,  17.85, 1.19e-1),
 (14.0,  8.40,-2.10e-2, 10.70,-8.27e0,  17.57, 1.18e1,  21.07, 1.31e-1),
 (16.0,  9.01,-1.28e-2, 11.97,-6.85e0,  20.21, 1.11e1,  24.26, 1.40e-1),
 (18.0,  9.58,-8.57e-3, 13.23,-5.96e0,  22.84, 1.07e1,  27.46, 1.46e-1),
 (20.0, 10.14,-6.12e-3, 14.50,-5.37e0,  25.47, 1.04e1,  30.65, 1.51e-1),
 (30.0, 12.67,-1.89e-3, 20.82,-4.00e0,  38.60, 9.60e0,  46.56, 1.66e-1),
 (40.0, 14.89,-8.69e-4, 27.13,-3.49e0,  51.72, 9.30e0,  62.45, 1.73e-1),
 (50.0, 16.91,-4.86e-4, 33.44,-3.22e0,  64.83, 9.13e0,  78.33, 1.77e-1),
 (75.0, 21.36,-1.72e-4, 49.21,-2.89e0,  97.59, 8.94e0, 118.03, 1.82e-1),
 (100.0,25.25,-8.33e-5, 64.97,-2.75e0, 130.36, 8.85e0, 157.72, 1.85e-1),
 (150.0,32.00,-3.01e-5, 96.47,-2.61e0, 195.88, 8.77e0, 237.10, 1.87e-1),
 (200.0,37.91,-1.46e-5,127.98,-2.54e0, 261.40, 8.73e0, 316.47, 1.88e-1),
 (250.0,43.25,-8.35e-6,159.48,-2.50e0, 326.92, 8.70e0, 395.84, 1.89e-1)]

const D0 = RelativisticDisc(a = 0.0)
r0s = [row[1] for row in PUB]
rI1 = Float64[]; rI2 = Float64[]; rO2 = Float64[]; rO1 = Float64[]
for r0 in r0s
    a1, b1 = lindblad_locations(D0, r0, 1)
    a2, b2 = lindblad_locations(D0, r0, 2)
    push!(rI1, a1); push!(rO1, b1); push!(rI2, a2); push!(rO2, b2)
end
println("=== r₁ LOCATIONS (computed vs published, % error) ===")
@printf("%6s | %17s | %17s | %17s | %17s\n", "r0", "m1 ILR", "m2 ILR", "m2 OLR", "m1 OLR")
locmax = 0.0
for (i, row) in enumerate(PUB)
    e = (abs(rI1[i]-row[2])/row[2], abs(rI2[i]-row[4])/row[4], abs(rO2[i]-row[6])/row[6], abs(rO1[i]-row[8])/row[8])
    global locmax = max(locmax, maximum(e))
    @printf("%6.0f | %8.3f(%+5.2f%%) | %8.3f(%+5.2f%%) | %8.3f(%+5.2f%%) | %8.3f(%+5.2f%%)\n",
            row[1], rI1[i], 100*(rI1[i]-row[2])/row[2], rI2[i], 100*(rI2[i]-row[4])/row[4],
            rO2[i], 100*(rO2[i]-row[6])/row[6], rO1[i], 100*(rO1[i]-row[8])/row[8])
end
@repro_result "r0" r0s
@repro_result "r1_m1_ILR" rI1
@repro_result "r1_m2_ILR" rI2
@repro_result "r1_m2_OLR" rO2
@repro_result "r1_m1_OLR" rO1
@repro_result "location_max_rel_err_vs_published" locmax
@printf("max |Δr₁|/r₁ vs the published 3-digit table = %.3e (the table is rounded to 2 decimals)\n\n", locmax)

println("=== ILR STRENGTHS N (computed vs published, ratio) ===")
N1 = Float64[]; N2 = Float64[]
@printf("%6s | %13s %10s %6s | %13s %10s %6s\n", "r0", "m1 ILR", "pub", "ratio", "m2 ILR", "pub", "ratio")
for (i, row) in enumerate(PUB)
    n1 = lindblad_torque(D0, row[1], 1; kind = :ILR)
    n2 = lindblad_torque(D0, row[1], 2; kind = :ILR)
    push!(N1, n1); push!(N2, n2)
    @printf("%6.0f | %13.4e %10.3e %6.3f | %13.4e %10.3e %6.3f\n",
            row[1], n1, row[3], n1/row[3], n2, row[5], n2/row[5])
end
@repro_result "N_m1_ILR" N1
@repro_result "N_m2_ILR" N2
rat1 = [N1[i]/PUB[i][3] for i in eachindex(N1)]
rat2 = [N2[i]/PUB[i][5] for i in eachindex(N2)]
@repro_result "N_m1_ILR_ratio_to_published" rat1
@repro_result "N_m2_ILR_ratio_to_published" rat2
@printf("\nratio to published: m1 ILR ∈ [%.3f, %.3f]; m2 ILR ∈ [%.3f, %.3f]\n",
        minimum(rat1), maximum(rat1), minimum(rat2), maximum(rat2))

exit(ReproHarness.conclude!(opts; paper = "1010.0759", id = "tabIlrS",
                            script = @__FILE__, rtol = 1e-8))

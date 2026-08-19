# arXiv:1010.0758 (Hirata I) RT8 — the paper's ONLY quantitative EXTERNAL validation: in the
# point-mass Keplerian limit r/M ≫ 1 the relativistic Lindblad torque density collapses onto
# Goldreich & Tremaine (1980, ApJ 241, 425) Eq. (18),
#
#   dT^L_{m,m}/dr = sgn(r−a) κ² r Σ /(2³A⁴) · (GM_s)²/(a−r)⁴ · [2(Ω/κ)K₀(κ/2|A|) + K₁(κ/2|A|)]²,
#
# whose PURE COEFFICIENT is κ²Ω²/(2³A⁴) = Ω⁴/(2³(¾Ω)⁴) = 32/81 = 0.395061728…  with the Keplerian
# Oort constant A = (r/2)dΩ/dr = −(3/4)Ω and κ = Ω. The derivation certificate is
# derivations/RT8_GT80_newtonian_check.tex (the Bessel moment K₀(2/3)=0.6967699885980613,
# K₁(2/3)=1.0658787418850444 is literature-grounded there, mpmath 25 dps, and is NOT recomputed here
# — no GrTools package owns modified Bessel functions).
#
# THIN driver: the frequencies come from `AccretionDisc.Modes` (disc_OmegaK, kappar2 — the exact Kerr
# equatorial geodesic values, whose a=0, r→∞ limit is Keplerian) and the resonance radii from
# `lindblad_locations`; the script only forms ratios and compares with the published constants.
#
#   julia --project=reproduce -t4 reproduce/arXiv-1010.0758/check_RT8_GT80.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf
using GrTools.AccretionDisc: RelativisticDisc, disc_OmegaK, kappar2, lindblad_locations

const D0 = RelativisticDisc(a = 0.0)

# Oort A(r) = (r/2) dΩ/dr, central difference with a relative step
oortA(r; h = 1e-6) = (r/2) * (disc_OmegaK(D0, r*(1+h)) - disc_OmegaK(D0, r*(1-h))) / (2*h*r)

println("Keplerian limit of the GT80 prefactor  κ²Ω²/(2³A⁴)  →  32/81 = 0.3950617284")
coef = Float64[]
for r in (1.0e3, 1.0e4, 1.0e5, 1.0e6)
    Ω = disc_OmegaK(D0, r); κ2 = kappar2(D0, r); A = oortA(r)
    c = κ2*Ω^2/(8*A^4)
    push!(coef, c)
    @printf("  r=%9.3g  Ω=%.10e  κ/Ω=%.12f  A/Ω=%+.12f  κ²Ω²/(2³A⁴)=%.12f  rel=%.2e\n",
            r, Ω, sqrt(κ2)/Ω, A/Ω, c, abs(c - 32/81)/(32/81))
end
@repro_result "GT80_coefficient_r1e3" coef[1]
@repro_result "GT80_coefficient_r1e6" coef[4]
@repro_result "GT80_coefficient_limit_relerr" abs(coef[4] - 32/81)/(32/81)
@repro_result "kappa_over_Omega_r1e6" sqrt(kappar2(D0, 1.0e6))/disc_OmegaK(D0, 1.0e6)
@repro_result "oortA_over_Omega_r1e6" oortA(1.0e6)/disc_OmegaK(D0, 1.0e6)

# RT4/RT5 tie: the Keplerian Lindblad radii ς_m^∓ = ((m∓1)/m)^{2/3} (the resonance geometry the
# closed-form S^{(m)} of RT5 is evaluated at) emerge from the same relativistic root solve.
println("\nKeplerian Lindblad ratios ς = r₁/r₀ → ((m∓1)/m)^{2/3}")
sI = Float64[]; sO = Float64[]
for m in 2:5
    rI, _ = lindblad_locations(D0, 1.0e6, m)
    _, rO = lindblad_locations(D0, 1.0e6, m)
    push!(sI, rI/1.0e6); push!(sO, rO/1.0e6)
    @printf("  m=%d  ς⁻=%.12f (%.12f)  ς⁺=%.12f (%.12f)\n",
            m, sI[end], ((m-1)/m)^(2/3), sO[end], ((m+1)/m)^(2/3))
end
@repro_result "sigma_ILR_m2to5" sI
@repro_result "sigma_OLR_m2to5" sO
@repro_result "sigma_ILR_max_relerr" maximum(abs(sI[i] - (i/(i+1))^(2/3))/(i/(i+1))^(2/3) for i in 1:4)
@repro_result "sigma_OLR_max_relerr" maximum(abs(sO[i] - ((i+2)/(i+1))^(2/3))/((i+2)/(i+1))^(2/3) for i in 1:4)

# the literature-grounded Bessel moment of GT80 Eq. (18) (certificate: derivations/RT8_*.tex)
@repro_result "GT80_bessel_K0_2over3" 0.6967699885980613
@repro_result "GT80_bessel_K1_2over3" 1.0658787418850444
@repro_result "GT80_bessel_moment" (2*0.6967699885980613 + 1.0658787418850444)^2
@printf("\nGT80 Bessel moment [2K₀(2/3)+K₁(2/3)]² = %.12f   (literature-grounded, mpmath 25 dps)\n",
        (2*0.6967699885980613 + 1.0658787418850444)^2)

exit(ReproHarness.conclude!(opts; paper = "1010.0758", id = "check_RT8_GT80",
                            script = @__FILE__, rtol = 1e-9))

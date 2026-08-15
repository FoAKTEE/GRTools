# Kerr bound geodesics: constants of motion, fundamental frequencies, ISCO/separatrix.
#
# Run from the repository root:
#     julia --project=. examples/01_kerr_geodesics.jl
#
# Conventions: G = c = M = 1; orbits parameterized by (a, p, e, x) — spin, semi-latus
# rectum, eccentricity, cos(inclination).

using GrTools                 # umbrella: every suite package is GrTools.<Name>
using GrTools.Spacetime
using Printf

a, p, e, x = 0.9, 8.0, 0.3, 0.8    # inclined eccentric orbit around a fast-spinning hole

com = constants_of_motion(a, p, e, x)
@printf "Kerr geodesic (a=%.2f, p=%.1f, e=%.2f, x=%.2f)\n" a p e x
@printf "  energy          E  = %.12f\n" com.E
@printf "  ang. momentum   Lz = %.12f\n" com.Lz
@printf "  Carter constant Q  = %.12f\n" com.Q

freq = kerr_frequencies(a, p, e, x)   # Boyer–Lindquist-time frequencies
@printf "  Omega_r  = %.12f\n" freq.Omegar
@printf "  Omega_th = %.12f\n" freq.Omegath
@printf "  Omega_ph = %.12f\n" freq.Omegaph

# Circular equatorial limit: Omega_phi -> the Kepler-like 1/(p^{3/2} + a) formula.
pc = 10.0
fc = kerr_frequencies(a, pc, 0.0, 1.0)
@printf "\nCircular equatorial p=%.1f:  Omega_ph = %.12f  (1/(p^1.5+a) = %.12f)\n" pc fc.Omegaph 1 / (pc^1.5 + a)

# ISCO and the separatrix between bound and plunging orbits.
@printf "\nISCO radius:  a=0 -> %.6f,  a=%.2f (prograde) -> %.6f\n" kerr_isco(0.0) a kerr_isco(a)
@printf "Separatrix p_s(a=%.2f, e=%.2f, x=%.2f) = %.8f\n" a e x separatrix(a, e, x)

# Type-generic: BigFloat inputs inherit the working precision end-to-end.
setprecision(BigFloat, 256) do
    comb = constants_of_motion(big"0.9", big"8.0", big"0.3", big"0.8")
    Δ = Float64(abs(comb.E - com.E))
    @printf "\nBigFloat-256 energy agrees with Float64 to |dE| = %.2e\n" Δ
    println("  E (BigFloat-256) = ", comb.E)
end

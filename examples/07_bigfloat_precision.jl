# End-to-end BigFloat precision inheritance across the suite.
#
# Run from the repository root:
#     julia --project=. examples/07_bigfloat_precision.jl
#
# Every numerical path in GrTools is type-generic in the element type: BigFloat
# inputs inherit the active setprecision(BigFloat, bits) end-to-end, with no
# silent Float64 waypoints. This example measures it: the same quantity computed
# at 256 and at 512 bits must agree to ~the 256-bit floor, far beyond Float64.

using Spacetime
using Harmonics
using Wigner
using Printf

# Self-convergence: |x(256 bits) - x(512 bits)| ~ eps(256 bits) ~ 1e-77.
compute(f) = (setprecision(() -> f(), BigFloat, 256), setprecision(() -> f(), BigFloat, 512))

diff_digits(lo, hi) = -log10(Float64(max(abs(BigFloat(lo) - hi), big"1e-120")))

# 1. Kerr geodesic energy (root-finding + algebra chain through Spacetime).
Elo, Ehi = compute(() -> constants_of_motion(big"0.9", big"8.0", big"0.3", big"0.8").E)
@printf "geodesic energy      : 256 vs 512 bits agree to %5.1f digits (Float64 has ~16)\n" diff_digits(Elo, Ehi)

# 2. Spheroidal eigenvalue (Leaver continued fraction through Harmonics).
λlo, λhi = compute(() -> SpinWeightedSpheroidalEigenvalue(-2, 2, 2, big"0.35"))
@printf "spheroidal eigenvalue: 256 vs 512 bits agree to %5.1f digits\n" diff_digits(λlo, λhi)

# 3. Wigner 3j (BigInt-exact combinatorics, downcast only at the boundary).
wlo, whi = compute(() -> wigner3j(BigFloat, 8, 6, 4, 2, -1, -1))
@printf "wigner3j(8,6,4;...)  : 256 vs 512 bits agree to %5.1f digits\n" diff_digits(wlo, whi)

# The Float64 entry points are the same code path at T = Float64.
E64 = constants_of_motion(0.9, 8.0, 0.3, 0.8).E
@printf "\nFloat64 energy differs from BigFloat-512 by %.2e (Float64 roundoff)\n" Float64(abs(E64 - Ehi))

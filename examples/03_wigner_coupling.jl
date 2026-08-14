# Angular-momentum coupling: 3j/6j symbols, Clebsch-Gordan, harmonic coupling.
#
# Run from the repository root:
#     julia --project=. examples/03_wigner_coupling.jl
#
# All primitives take a leading element type; the combinatorics are BigInt-exact
# and downcast only at the boundary, so BigFloat requests are exact to precision.

using Wigner
using Printf

# A classic closed form: (1 1 2; 1 -1 0) = 1/sqrt(30).
w = wigner3j(Float64, 1, 1, 2, 1, -1, 0)
@printf "wigner3j(1,1,2; 1,-1,0) = %.15f   (1/sqrt(30) = %.15f)\n" w 1 / sqrt(30)

# 3j orthogonality: sum_{m1,m2} (2j3+1) (j1 j2 j3; m1 m2 m3)^2 = 1.
j1, j2, j3, m3 = 2, 3, 4, 1
acc = sum((2j3 + 1) * wigner3j(Float64, j1, j2, j3, m1, m2, -m3)^2
          for m1 in -j1:j1, m2 in -j2:j2 if m1 + m2 == m3)
@printf "3j orthogonality sum (j=%d,%d,%d, m3=%d) = %.15f (exact 1)\n" j1 j2 j3 m3 acc

# Clebsch-Gordan from the same Racah core: clebsch_gordan(T, j1, m1, j2, m2, j3, m3).
cg = clebsch_gordan(Float64, 1, 1, 1, 1, 2, 2)
@printf "<1 1; 1 1 | 2 2> = %.15f (exact 1)\n" cg

# 6j symbol.
w6 = wigner6j(Float64, 1, 2, 3, 2, 1, 2)
println("wigner6j(1,2,3; 2,1,2) = ", w6)

# 3-harmonic coupling, coupling(T, lc,mc,sc, la,ma,sa, lb,mb,sb) — the canonical home
# of the second-order Ricci source coupling. Y_10 x Y_10 -> monopole: exactly 1/sqrt(4 pi).
c = coupling(Float64, 0, 0, 0, 1, 0, 0, 1, 0, 0)
@printf "coupling(Y_10 x Y_10 -> Y_00) = %.15f   (1/sqrt(4pi) = %.15f)\n" c 1 / sqrt(4π)
c2 = coupling(Float64, 4, 4, -4, 2, 2, -2, 2, 2, -2)
println("coupling((2,2,s=-2) x (2,2,s=-2) -> (4,4,s=-4)) = ", c2)

# BigFloat: the same symbol at 256-bit working precision.
setprecision(BigFloat, 256) do
    wb = wigner3j(BigFloat, 1, 1, 2, 1, -1, 0)
    @printf "\nBigFloat-256 3j: |w - 1/sqrt(big(30))| = %.2e\n" Float64(abs(wb - 1 / sqrt(big(30))))
end

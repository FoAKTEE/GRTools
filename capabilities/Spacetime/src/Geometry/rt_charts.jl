"""
    Spacetime.RTChart

Pure chart maps of the Rácz–Tóth (RT) horizon-penetrating hyperboloidal coordinate system
[Harms, Bernuzzi & Brügmann, arXiv:1301.1591] — the coordinate surface of the generic-spin
time-domain Teukolsky formulation. Re-homed here from `Perturbation.TimeDomain.rtchart` (F5 Part B,
MINIMAL move): ONLY the four PURE maps `rt_r_of_R`, `rt_R_of_r`, `rt_Rplus`, `rt_F_of_R` (no spin
weight `s`, no mode list `ℓs`, no master-equation coefficients) live here. The RT master-equation
COEFFICIENT fields (`rt_C*`, `rt_A*`, `rt_C0diag`, …) stay in `Perturbation.TimeDomain` — they are
operator-APPLICATION geometry (coefficients-vs-application principle), not pure chart maps.

Nested + namespaced (the NPGHP / SpinningSecondary collision-avoidance precedent): the RT symbols
are re-imported and re-exported by `Perturbation.TimeDomain` as the SAME objects (`===`), so the
locked time-domain physics keeps its `rt_*` call sites bit-untouched (Move-A Lorenz shim pattern).

CHART [arXiv:1301.1591, Eqs. (3)–(6); classic RT parameters κ = S = 1]:
  ingoing-Kerr:  τ = t − r + ∫dr (r²+a²)/Δ,   φ̃ = φ + a∫dr/Δ
  RT:            τ = T − 4M log(1−R²) + (1+R²)/(1−R²),   r = 2R/(1−R²)
scri⁺ at R = 1, horizon at R₊ [Eq. (7)].
"""
module RTChart

"""
    rt_r_of_R(R) -> r

Boyer–Lindquist radius from the RT compactified radius, `r = 2R/(1−R²)`
[arXiv:1301.1591, Eq. (6)]; r(1⁻) = ∞ (scri⁺ at R = 1).
"""
rt_r_of_R(R) = 2R / (1 - R^2)

"RT compactified radius of a BL radius (inverse of `rt_r_of_R`): R = (√(1+r²) − 1)/r."
rt_R_of_r(r) = (sqrt(1 + r^2) - 1) / r

"""
    rt_Rplus(M, a) -> R₊

RT coordinate location of the outer horizon r₊ = M + √(M²−a²)
[arXiv:1301.1591, Eq. (7)].
"""
rt_Rplus(M, a) =
    (2 * sqrt(2M * sqrt(M^2 - a^2) - a^2 + 2M^2 + 1) - 2) / (2 * (sqrt(M^2 - a^2) + M))

"""
    rt_F_of_R(R; M=1) -> t − T   (Schwarzschild, a=0)

Composition of the chart maps [arXiv:1301.1591, Eqs. (3),(5),(6)] at a=0: with
τ = t − r + r_* (Eq. (3), ∫dr r²/Δ = r_* = r + 2M log(r/2M − 1)) and
τ = T − 4M log(1−R²) + (1+R²)/(1−R²) (Eq. (5)),

    t = T + F(R),   F = −4M log(1−R²) + (1+R²)/(1−R²) − r_*(r(R)) + r(R).

Used to push the Bardeen–Press (t,r_*) operator [arXiv:1702.04204, Eqs. (67)–(70)] into the RT
chart, and by the sourced-sector monochromatic jump amplitudes.
"""
function rt_F_of_R(R; M=1.0)
    r = rt_r_of_R(R)
    rstar = r + 2M * log(r / (2M) - 1)
    return -4M * log(1 - R^2) + (1 + R^2) / (1 - R^2) - rstar + r
end

end # module RTChart

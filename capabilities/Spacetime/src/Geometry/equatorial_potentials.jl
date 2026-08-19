# Geometry — Page–Thorne-form potentials of the Kerr EQUATORIAL line element (arXiv:2608.17003
# main.tex:148 eq (4); appendix.tex:7-21 eqs A1a–e), M = 1 like every Spacetime neighbor:
#
#     ds² = −e^{2ν} dt² + e^{2ψ} (dφ − ω̃ dt)² + e^{2μ} dr² + dz²,          (θ = π/2, z = r cosθ → 0)
#
# so that e^{2ψ} = g_φφ, ω̃ = −g_tφ/g_φφ (frame-dragging angular velocity), e^{2ν} = −g_tt +
# g_tφ²/g_φφ (lapse², = −1/g^{tt}), e^{2μ} = g_rr. Pure arithmetic over `delta_kerr`, so BigFloat
# precision is inherited; the circular dt/dτ REUSES `kerr_circular_ut_uphi` (no second derivation).
# Additive (EGD::spacetime-kerr-equatorial-potentials): no existing Spacetime function is modified.

"""
    kerr_equatorial_potentials(a, r) -> (ν, ψ, μ, ω̃)

Page–Thorne potentials of the Kerr equatorial line element (M = 1, χ = a, arXiv:2608.17003
eq (4) + appendix A1a–e):

    ds² = −e^{2ν} dt² + e^{2ψ} (dφ − ω̃ dt)² + e^{2μ} dr² + dz²,

    ω̃     = 2a / D,            D ≡ 2a² + 2r² + rΔ = r³ + a² r + 2a²,
    e^{2ψ} = 2a²/r + 2r + Δ = D/r            (= g_φφ),
    e^{2ν} = 1 − 2/r + 2aω̃/r = rΔ/D          (= lapse² = −1/g^{tt}; the second form is an exact
                                               rearrangement, manifestly positive for Δ > 0),
    e^{2μ} = r²/Δ                            (= g_rr),   Δ = r² − 2r + a² = `delta_kerr`.

Valid on the exterior r > r₊ (Δ > 0), INCLUDING the equatorial ergoregion r₊ < r < 2 where
g_tt > 0 but the lapse stays positive. Inside the horizon (Δ < 0) `log` of a negative argument
throws `DomainError`. `a > 0` prograde, `a < 0` retrograde (ω̃ flips sign). Schwarzschild
limit a = 0: ν = ½log(1−2/r), ψ = log r, μ = −½log(1−2/r), ω̃ = 0. Type-generic: BigFloat and
Float32 inherit, Int inputs promote.
"""
function kerr_equatorial_potentials(a, r)
    T = float(promote_type(typeof(a), typeof(r)))
    aa, rr = T(a), T(r)
    Δ = delta_kerr(aa, rr)
    D = 2 * aa^2 + 2 * rr^2 + rr * Δ           # r³ + a² r + 2a² > 0 on the exterior
    ω̃ = 2 * aa / D
    (ν = log(rr * Δ / D) / 2, ψ = log(D / rr) / 2, μ = log(rr^2 / Δ) / 2, ω̃ = ω̃)
end

"""
    kerr_equatorial_potentials_dr(a, r) -> (dν, dψ, dμ, dω̃)

Analytic radial derivatives of [`kerr_equatorial_potentials`](@ref) (M = 1), from the
logarithmic forms e^{2ν} = rΔ/D, e^{2ψ} = D/r, e^{2μ} = r²/Δ, ω̃ = 2a/D with Δ' = 2(r−1),
D' = 3r² + a²:

    dν/dr = ½ (1/r + Δ'/Δ − D'/D) = [(r² + a²)² − 4a² r] / (r Δ D),
    dψ/dr = ½ (D'/D − 1/r)        = (r³ − a²) / (r D),
    dμ/dr = 1/r − (r−1)/Δ         = (a² − r) / (r Δ),
    dω̃/dr = −ω̃ D'/D              = −2a (3r² + a²) / D².

The right-hand forms are the EXACT single-fraction rearrangements (numerators expanded
symbolically; a = 0 ⇒ 1/(r(r−2)), 1/r, −1/(r(r−2)), 0) and are the ones evaluated: the
logarithmic sums cancel catastrophically at large r (the three O(1/r) terms of dν sum to
O(1/r²), losing ~r·eps relative accuracy — 1e-13 at r = 10³ in Float64), whereas the
single-fraction forms are cancellation-free (relative error a few eps at every r > r₊).

Same domain (r > r₊) and type-genericity as the potentials (no `log`, so it is finite but
sign-indefinite across the horizon).
"""
function kerr_equatorial_potentials_dr(a, r)
    T = float(promote_type(typeof(a), typeof(r)))
    aa, rr = T(a), T(r)
    a2 = aa^2
    Δ = delta_kerr(aa, rr)
    D = 2 * a2 + 2 * rr^2 + rr * Δ              # r³ + a² r + 2a²
    dD = 3 * rr^2 + a2
    (dν = ((rr^2 + a2)^2 - 4 * a2 * rr) / (rr * Δ * D),
     dψ = (rr^3 - a2) / (rr * D),
     dμ = (a2 - rr) / (rr * Δ),
     dω̃ = -2 * aa * dD / D^2)
end

"""
    kerr_circular_dtdtau(a, r) -> dt/dτ

`u^t = dt/dτ` of the prograde (`a > 0`) / retrograde (`a < 0`) circular equatorial Kerr geodesic
at radius `r` (M = 1), arXiv:2608.17003 main.tex:673 with a⋆ = r/M, χ = a:

    dt/dτ = (a⋆^{3/2} + χ) / (a⋆^{3/4} √(√a⋆ (a⋆ − 3) + 2χ)) ≡ (a + r^{3/2}) / √(r³ − 3r² + 2a r^{3/2}),

which is EXACTLY [`kerr_circular_ut_uphi`](@ref)`(a, r).ut` (the denominator identity
a⋆^{3/4}√(√a⋆(a⋆−3)+2χ) = √(r³−3r²+2a r^{3/2})); this is a named accessor for the paper's
worldline four-velocity `v^μ = (dt/dτ, 0, Ω dt/dτ, 0)`, NOT a second implementation.
"""
kerr_circular_dtdtau(a, r) = kerr_circular_ut_uphi(a, r).ut

# Geodesics — circular equatorial Kerr four-velocity closed forms + the Kerr azimuthal chart
# offset φ̄ − φ (arXiv:2606.04998 main.tex:192-196, 223-226; the chart follows Dolan 2011b).
# M = 1 like every Spacetime neighbor; prograde/retrograde selected by the SIGN of a (the
# azimuthal motion is always +φ, the spin is mirrored — main.tex:196-198). Pure arithmetic over
# the existing Geometry primitives `horizons`/`delta_kerr`, so BigFloat precision is inherited.
# Additive W1-C (elliptic-SSF): no existing Spacetime function is modified.

"""
    kerr_circular_ut_uphi(a, r0) -> (ut, uphi)

Boyer–Lindquist four-velocity components `u^t`, `u^φ` of a circular equatorial Kerr geodesic at
radius `r0` (M = 1); `a > 0` prograde, `a < 0` retrograde — the motion is always in +φ
(`uphi > 0`), the spin axis is mirrored. Closed forms (arXiv:2606.04998 main.tex:192-196):

    u^t = (a + r0^{3/2}) / √(r0³ − 3 r0² + 2 a r0^{3/2}),   u^φ = 1 / √(⋯),

so the orbital frequency is `Ω = uphi/ut = 1/(a + r0^{3/2})` (= `kerr_frequencies` circ-eq
`Omegaph`). Valid for `r0` outside the photon orbit, where the radicand is positive (a real-sqrt
`DomainError` signals otherwise). Type-generic: BigFloat is inherited, never down-cast.
"""
function kerr_circular_ut_uphi(a, r0)
    aa, rr = promote(a, r0)
    v = sqrt(rr) * rr                                  # r0^{3/2}
    den = sqrt(rr^2 * (rr - 3) + 2 * aa * v)           # √(r0³ − 3r0² + 2a·r0^{3/2})
    (ut = (aa + v) / den, uphi = one(den) / den)
end

"""
    kerr_varphi_offset(a, r) -> φ̄ − φ

Offset of Kerr's original azimuthal coordinate `φ̄` from Boyer–Lindquist `φ` (M = 1), the
horizon-regular chart used by the m-mode ansatz (arXiv:2606.04998 main.tex:223-226, following
Dolan 2011b):

    φ̄ − φ = a/(r₊−r₋) · ln((r−r₊)/(r−r₋)),

on the exterior `r > r₊` with `|a| < 1` (at a = 0 it vanishes identically; at the extremal
limit r₊ = r₋ the closed form degenerates). Diverges like −sign(a)·∞ at the horizon.
Type-generic via `horizons`.
"""
function kerr_varphi_offset(a, r)
    aa, rr = promote(a, r)
    rp, rm = horizons(aa)
    aa / (rp - rm) * log((rr - rp) / (rr - rm))
end

"""
    kerr_varphi_offset_dr(a, r) -> a/Δ

Radial derivative `d(φ̄−φ)/dr = a/Δ(a,r)` of `kerr_varphi_offset` (M = 1) — exact, since
Δ = (r−r₊)(r−r₋). Well-defined for all `r` off the horizons, including the extremal limit.
"""
kerr_varphi_offset_dr(a, r) = a / delta_kerr(a, r)

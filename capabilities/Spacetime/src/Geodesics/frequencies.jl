# Geodesics — radial/polar roots + Fujita–Hikida Mino/BL/proper orbital frequencies.
# Method: Fujita & Hikida, CQG 26 (2009) 135002 (arXiv:0906.1420). Faithful port of
# KerrGeodesics`OrbitalFrequencies, here built on NumericsCore's Mathematica-convention
# ellipK/E/Pi (the same complex Carlson layer the legacy elliptic.jl path used) and the
# E/Lz/Q closed forms in constants.jl. Elliptic *parameters* are mathematically real for
# bound orbits; coerce with `real` before the Legendre call (reproduces Mathematica's branch).

@inline _erc(z) = real(z)   # elliptic-parameter real coercion (args are real for bound orbits)

"""
    radial_roots(a, p, e, x[, E, Q]) -> (r1, r2, r3, r4)

The four roots of the Kerr radial potential (M = 1): r1,r2 the real turning points
p/(1∓e), r3,r4 (possibly `Complex`) from the quartic. Replaces `KerrGeoRadialRoots`.
"""
function radial_roots(a, p, e, x, En, Q)
    r1 = p / (1 - e)
    r2 = p / (1 + e)
    b  = 2 / (1 - En^2) - r1 - r2
    disc = b^2 - (4 * a^2 * (1 - e) * (1 + e) * Q) / ((1 - En^2) * p^2)
    s  = sqrt(complex(disc))
    r3 = (b + s) / 2
    r4 = (2 * a^2 * (1 - e) * (1 + e) * Q) / ((1 - En^2) * p^2 * (b + s))
    (r1, r2, r3, r4)
end
radial_roots(a, p, e, x) =
    radial_roots(a, p, e, x, geodesic_energy(a, p, e, x), carter_constant(a, p, e, x))

"""
    polar_roots(a, p, e, x) -> (zp, zm)

Polar turning-point parameters (zm = sin of the min polar angle, zp the second root of the
θ-potential). Replaces `KerrGeoPolarRoots`.
"""
function polar_roots(a, p, e, x)
    En = geodesic_energy(a, p, e, x); L = geodesic_angular_momentum(a, p, e, x); Q = carter_constant(a, p, e, x)
    zm = sqrt(1 - x^2)
    zp = x == 0 ? sqrt(Q) : sqrt(a^2*(1 - En^2) + L^2/(1 - zm^2))
    (zp, zm)
end

# ---- Mino-time fundamental frequencies (generic Kerr) ----------------------------------------
function _mino_freq_r(a, En, r1, r2, r3, r4)
    kr = _erc((r1 - r2)/(r1 - r3) * (r3 - r4)/(r2 - r4))
    (pi*sqrt((1 - En^2)*(r1 - r3)*(r2 - r4))) / (2*ellipK(kr))
end

function _mino_freq_theta(a, En, zp, zm)
    ktheta = _erc(a^2*(1 - En^2)*(zm/zp)^2)
    (pi*zp) / (2*ellipK(ktheta))
end

function _mino_freq_phi_r(a, En, L, r1, r2, r3, r4)
    rin  = 1 - sqrt(1 - a^2); rout = 1 + sqrt(1 - a^2)
    kr   = _erc((r1 - r2)/(r1 - r3) * (r3 - r4)/(r2 - r4))
    hout = _erc((r1 - r2)/(r1 - r3) * (r3 - rout)/(r2 - rout))
    hin  = _erc((r1 - r2)/(r1 - r3) * (r3 - rin)/(r2 - rin))
    Kk = ellipK(kr)
    a/(2*sqrt(1 - a^2)) * (
        (2*En*rout - a*L)/(r3 - rout) * (1 - (r2 - r3)/(r2 - rout) * ellipPi(hout, kr)/Kk) -
        (2*En*rin  - a*L)/(r3 - rin ) * (1 - (r2 - r3)/(r2 - rin ) * ellipPi(hin , kr)/Kk))
end

function _mino_freq_phi_theta(a, En, L, zp, zm)
    ktheta = _erc(a^2*(1 - En^2)*(zm/zp)^2)
    L*ellipPi(_erc(zm^2), ktheta)/ellipK(ktheta)
end

function _mino_freq_t_r(a, En, L, r1, r2, r3, r4)
    rin  = 1 - sqrt(1 - a^2); rout = 1 + sqrt(1 - a^2)
    kr   = _erc((r1 - r2)/(r1 - r3) * (r3 - r4)/(r2 - r4))
    hout = _erc((r1 - r2)/(r1 - r3) * (r3 - rout)/(r2 - rout))
    hin  = _erc((r1 - r2)/(r1 - r3) * (r3 - rin)/(r2 - rin))
    hr   = _erc((r1 - r2)/(r1 - r3))
    Kk = ellipK(kr); Ek = ellipE(kr); Phr = ellipPi(hr, kr)
    (a^2 + 4)*En + En*(
        (r3*(r1 + r2 + r3) - r1*r2 + (r1 + r2 + r3 + r4)*(r2 - r3)*Phr/Kk +
         (r1 - r3)*(r2 - r4)*Ek/Kk)/2 +
        2*(r3 + (r2 - r3)*Phr/Kk) +
        1/sqrt(1 - a^2) * (
            ((4 - a*L/En)*rout - 2*a^2)/(r3 - rout) * (1 - (r2 - r3)/(r2 - rout)*ellipPi(hout, kr)/Kk) -
            ((4 - a*L/En)*rin  - 2*a^2)/(r3 - rin ) * (1 - (r2 - r3)/(r2 - rin )*ellipPi(hin , kr)/Kk)))
end

function _mino_freq_t_theta(a, En, Q, zp, zm, x)
    x^2 == 1 && return complex(-a^2*En)
    ktheta = _erc(a^2*(1 - En^2)*(zm/zp)^2)
    (En*Q)/((1 - En^2)*zm^2) * (1 - ellipE(ktheta)/ellipK(ktheta)) - a^2*En
end

# Schwarzschild (a=0) closed forms.
function _mino_freqs_schwarz(p, e, x)
    sgnx = x < 0 ? -one(p) : one(p)
    if e == 0
        return (Ur = sqrt((-6 + p)*p/(-3 + p)),
                Uth = p*sqrt(1/(-3 + p)),
                Uph = p*sgnx/sqrt(-3 + p),
                Ut  = sqrt(p^5/(-3 + p)))
    end
    m  = (4*e)/(-6 + 2*e + p)
    Km = ellipK(m); Em = ellipE(m)
    Pi1 = ellipPi((2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), m)
    Pi2 = ellipPi((16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), m)
    Ur = (sqrt(-((p*(-6 + 2*e + p))/(3 + e^2 - p)))*pi)/(2*Km)
    Uth = p/sqrt(-3 - e^2 + p)
    Uph = (p*sgnx)/sqrt(-3 - e^2 + p)
    inner = -(((-4 + p)*p^2*(-6 + 2*e + p)*Em)/(-1 + e^2)) +
            (p^2*(28 + 4*e^2 - 12*p + p^2)*Km)/(-1 + e^2) -
            (2*(6 + 2*e - p)*(3 + e^2 - p)*p^2*Pi1)/((-1 + e)*(1 + e)^2) +
            (4*(-4 + p)*p*(2*(1 + e)*Km + (-6 - 2*e + p)*Pi1))/(1 + e) +
            2*(-4 + p)^2*((-4 + p)*Km - ((6 + 2*e - p)*p*Pi2)/(2 + 2*e - p))
    Ut = sqrt((-4*e^2 + (-2 + p)^2)/(p*(-3 - e^2 + p)))/2 * (8 + inner/((-4 + p)^2*Km))
    (Ur = Ur, Uth = complex(Uth), Uph = complex(Uph), Ut = Ut)
end

# --- the circular-equatorial Υ_r radicand and its exact ISCO zero -----------------------------
#
# On the stable branch the sign factor `disc/|disc|` in `_mino_freqs_circeq` is −1, because
#     disc = a² − 4a√p − (p−4)p = −(p + 2√p − a)(p − 2√p + a) < 0
# for every p above the marginally bound orbit r_mb = 2 − a + 2√(1−a) (both factors are then
# positive). With that sign the Υ_r radicand collapses EXACTLY — no approximation — to
#
#     Υ_r² = p · P_isco(p) / (2a√p + (p−3)p),      P_isco(p) = p² − 6p + 8a√p − 3a² ,
#
# and P_isco is precisely the Bardeen–Press–Teukolsky (1972, ApJ 178 347) polynomial whose root
# is `kerr_isco(a)` — it is r² times the ISCO condition 0 = 1 − 6M/r + 8a√(M/r³) − 3a²/r² of
# arXiv:2606.04998 main.tex:207-209 (d²V_eff/dr² = 0). So Υ_r has a SIMPLE ZERO at
# the ISCO, is positive outside it, and is negative only for the unstable circular orbits inside
# it. Two consequences, handled by `_circeq_isco_radicand` and nowhere else:
#
#  1. At p = kerr_isco(a) the computed radicand is pure cancellation error — measured −2.2e−16
#     (a = 0.9) and −1.7e−14 (a = 0.998), which made `sqrt` raise a `DomainError`, and +1.3e−15
#     (a = 0.5), which returned a spurious Ω_r = 1.2e−9. There Υ_r is set to exactly 0, its
#     physical value. The trigger is the exact identity `p == kerr_isco(a)`, NOT a tolerance:
#     that is what makes it precision-independent (at 256 bits the residual is still ~1e−76,
#     i.e. nonzero, and no fixed tolerance would cover both regimes).
#  2. A NEGATIVE radicand within the first-order rounding bound of that zero carries no
#     significant digit either, so it is snapped to 0 as well — this catches a caller whose `p`
#     is the ISCO only to a few ulp. The bound is K·eps(T)·Σ|terms|/|denominator| over the terms
#     actually summed, with K = 16; the measured residual at `kerr_isco(a)` uses at most 0.19 of
#     it for |a| ≥ 0.5, and the window stays inside |δp/p| ≤ 6.2e−13 even at a = 0.998. It is a
#     roundoff guard only — near a → 0 `kerr_isco`'s own cancellation error exceeds this bound,
#     which is precisely why (1) and not (2) is the mechanism that makes the ISCO exact.
#
# A radicand below −tol is a genuine inside-the-ISCO (unstable) orbit and keeps raising the
# original `DomainError`; this is deliberately NOT a blanket `max(x, 0)` clamp. Every p that
# previously returned a finite Υ_r is bit-for-bit unchanged: the radicand expression itself is
# untouched, and the only value that moves is the ISCO point where 0 is exact.
const _CIRCEQ_ISCO_ROUNDOFF_K = 16

function _circeq_isco_radicand(rad, a, p)
    T = float(promote_type(typeof(a), typeof(p)))
    p == kerr_isco(T(a)) && return zero(rad)     # Υ_r ≡ 0 at the ISCO (exact, any precision)
    rad < 0 || return rad                        # positive radicand: never touched
    sq = sqrt(p)
    scale = (p*(2*a^2 + 6*a*sq + abs((-5 + p)*p) + (a - sq)^2)) / abs(2*a*sq + (-3 + p)*p)
    rad ≥ -_CIRCEQ_ISCO_ROUNDOFF_K*eps(one(rad))*scale ? zero(rad) : rad
end

# Circular-equatorial Kerr (a≠0, e=0, x=1) closed forms.
function _mino_freqs_circeq(a, p)
    den = sqrt(2*a + (-3 + p)*sqrt(p))
    disc = a^2 - 4*a*sqrt(p) - (-4 + p)*p
    rad = (p*(-2*a^2 + 6*a*sqrt(p) + (-5 + p)*p + ((a - sqrt(p))^2*disc)/abs(disc))) /
          (2*a*sqrt(p) + (-3 + p)*p)
    Ur = sqrt(_circeq_isco_radicand(rad, a, p))
    Uth = abs((p^(1//4)*sqrt(3*a^2 - 4*a*sqrt(p) + p^2))/den)
    Uph = p^(5//4)/den
    Ut  = (p^(5//4)*(a + p^(3//2)))/den
    (Ur = complex(Ur), Uth = complex(Uth), Uph = complex(Uph), Ut = complex(Ut))
end

function _mino_frequencies(a, p, e, x)
    if a < 0   # retrograde-spin symmetry: flip a→-a, x→-x, negate Υφ
        f = _mino_frequencies(-a, p, e, -x)
        return (Ur = f.Ur, Uth = f.Uth, Uph = -f.Uph, Ut = f.Ut)
    end
    if a == 0
        f = _mino_freqs_schwarz(p, e, x)
        return (Ur = real(f.Ur), Uth = abs(f.Uth), Uph = real(f.Uph), Ut = real(f.Ut))
    end
    if e == 0 && x == 1
        f = _mino_freqs_circeq(a, p)
        return (Ur = real(f.Ur), Uth = abs(f.Uth), Uph = real(f.Uph), Ut = real(f.Ut))
    end
    En = geodesic_energy(a, p, e, x); L = geodesic_angular_momentum(a, p, e, x); Q = carter_constant(a, p, e, x)
    r1, r2, r3, r4 = radial_roots(a, p, e, x, En, Q)
    zp, zm = polar_roots(a, p, e, x)
    Ur  = _mino_freq_r(a, En, r1, r2, r3, r4)
    Uth = _mino_freq_theta(a, En, zp, zm)
    Uph = _mino_freq_phi_r(a, En, L, r1, r2, r3, r4) + _mino_freq_phi_theta(a, En, L, zp, zm)
    Ut  = _mino_freq_t_r(a, En, L, r1, r2, r3, r4) + _mino_freq_t_theta(a, En, Q, zp, zm, x)
    (Ur = real(Ur), Uth = abs(Uth), Uph = real(Uph), Ut = real(Ut))
end

"""
    kerr_frequencies(a, p, e, x; Time="BoyerLindquist") -> NamedTuple

Fundamental orbital frequencies of a generic bound Kerr geodesic (M = 1). `Time="Mino"` →
`(Ur, Uth, Uph, Ut)`; `"BoyerLindquist"` → `(Omegar, Omegath, Omegaph)` = Mino/Ut; `"Proper"`
→ Mino/⟨dτ/dλ⟩. Replaces `KerrGeoFrequencies`.

ISCO limit: for a prograde circular equatorial orbit (`a > 0, e = 0, x = 1`) the radial
frequency vanishes at the innermost stable circular orbit, and `kerr_frequencies(a,
kerr_isco(a), 0, 1)` returns `Υ_r = Ω_r = 0` EXACTLY in every precision — the closed-form Υ_r
radicand is `p·(p² − 6p + 8a√p − 3a²)/(2a√p + (p−3)p)`, whose numerator is the Bardeen ISCO
polynomial itself. Radii inside the ISCO are unstable circular orbits with a genuinely negative
radicand and still raise a `DomainError`. Retrograde circular equatorial orbits (`a < 0`) route
through the generic branch, where the ISCO leaves a `|Ω_r| ≲ 1e-9` numerical residual instead.
"""
function kerr_frequencies(a, p, e, x; Time = "BoyerLindquist")
    f = _mino_frequencies(a, p, e, x)
    if Time == "Mino"
        return (Ur = f.Ur, Uth = f.Uth, Uph = f.Uph, Ut = f.Ut)
    elseif Time == "BoyerLindquist"
        return (Omegar = f.Ur/f.Ut, Omegath = f.Uth/f.Ut, Omegaph = f.Uph/f.Ut)
    elseif Time == "Proper"
        P = proper_frequency_factor(a, p, e, x)
        return (Omegar = f.Ur/P, Omegath = f.Uth/P, Omegaph = f.Uph/P)
    else
        throw(ArgumentError("kerr_frequencies: unknown Time=$(repr(Time)) (use \"Mino\", \"BoyerLindquist\", or \"Proper\")"))
    end
end

"""
    proper_frequency_factor(a, p, e, x) -> ⟨dτ/dλ⟩

Proper-time/Mino-time conversion factor (M = 1). Replaces `KerrGeoProperFrequencyFactor`.
"""
function proper_frequency_factor(a, p, e, x)
    e >= 1 && return Inf
    if a == 0
        e == 0 && return p^2
        m  = (4*e)/(-6 + 2*e + p)
        Km = ellipK(m); Em = ellipE(m)
        Pi1 = ellipPi((2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), m)
        val = (p^2*((1 + e)*(28 + 4*e^2 + (-12 + p)*p) -
               ((1 + e)*(-4 + p)*(-6 + 2*e + p)*Em + 2*(6 + 2*e - p)*(3 + e^2 - p)*Pi1)/Km)) /
              (2*(-1 + e)*(1 + e)^2*(-4 + p)^2)
        return real(val)
    end
    r1, r2, r3, r4 = radial_roots(a, p, e, x)
    zp, zm = polar_roots(a, p, e, x)
    En = geodesic_energy(a, p, e, x)
    kr = _erc((r1 - r2)/(r1 - r3) * (r3 - r4)/(r2 - r4))
    ktheta = _erc(a^2*(1 - En^2)*(zm/zp)^2)
    hr = _erc((r1 - r2)/(r1 - r3))
    val = (-((2*zp^2)/(-1 + En^2)) + r1*(-r2 + r3) + r3*(r2 + r3))/2 +
          ((r1 - r3)*(r2 - r4)*ellipE(kr))/(2*ellipK(kr)) +
          (zp^2*ellipE(ktheta))/((-1 + En^2)*ellipK(ktheta)) +
          ((r2 - r3)*(r1 + r2 + r3 + r4)*ellipPi(hr, kr))/(2*ellipK(kr))
    real(val)
end

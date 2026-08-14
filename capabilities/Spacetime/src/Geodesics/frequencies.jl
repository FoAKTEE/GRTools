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

# Circular-equatorial Kerr (a≠0, e=0, x=1) closed forms.
function _mino_freqs_circeq(a, p)
    den = sqrt(2*a + (-3 + p)*sqrt(p))
    disc = a^2 - 4*a*sqrt(p) - (-4 + p)*p
    Ur = sqrt((p*(-2*a^2 + 6*a*sqrt(p) + (-5 + p)*p + ((a - sqrt(p))^2*disc)/abs(disc))) /
              (2*a*sqrt(p) + (-3 + p)*p))
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

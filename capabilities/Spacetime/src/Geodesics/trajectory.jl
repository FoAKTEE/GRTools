# Geodesics — closed-form bound-orbit worldlines: Mino-phase parametrised r(qr)/θ(qz) plus the
# oscillating t, φ, τ deviations from their secular drift. Faithful consolidation of
# KerrGeodesics`kerrgeo_trajectory.jl (KerrGeoTrajectory* + KerrGeoProperTime), rebuilt on
# NumericsCore's Mathematica-convention elliptic/Jacobi layer and the in-package constants/roots.
# `_erc` (real-coercion of mathematically-real elliptic parameters) is shared from frequencies.jl.
# M = 1, bound orbits. Type-generic {T}: arithmetic only ⇒ BigFloat inherited.

"""
    geodesic_trajectory_equatorial(a, p, e) -> Function

Radial worldline of a bound equatorial Kerr geodesic (M = 1). Returns
`qr -> (; r, Ur, Δtr, Δφr, Δτr)`: BL radius `r(qr)`, signed Mino velocity `Ur = dr/dλ`, and the oscillating
deviations of t, φ, τ (the SpinningTrajectory "TrajectoryDeltas" convention), as functions of the
radial Mino phase `qr ∈ [0,2π)`. Replaces `KerrGeoTrajectoryEquatorial`.
"""
function geodesic_trajectory_equatorial(a, p, e)
    radial = geodesic_trajectory_radial(a, p, e, 1)
    # Match the signed-spin symmetry used by kerr_frequencies.
    En, _, Q = constants_of_motion(abs(a), p, e, a < 0 ? -1 : 1)
    r1, r2, r3, r4 = real.(radial_roots(a, p, e, 1, En, Q))
    kr = (r1 - r2) * (r3 - r4) / ((r1 - r3) * (r2 - r4))
    hr = (r1 - r2) / (r1 - r3)
    Kk = real(ellipK(kr)); Ek = real(ellipE(kr)); Pihr = real(ellipPi(hr, kr))
    sden = sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
    return function (qr)
        ψr = jacobiAM(Kk / pi * qr, kr)
        sn, cn = sincos(ψr)
        ΔE = Ek * qr / pi - ellipEinc(ψr, kr)
        ΔΠhr = Pihr * qr / pi - ellipPiinc(hr, ψr, kr)
        # cn carries the inward sign on the second half-period. The same continued
        # amplitude is used by the shared radial t/φ trajectory.
        Δτr = -1 / sden * ((r2 - r3) * 2 / (1 - En^2) * ΔΠhr
            + (r2 - r4) * ((r1 - r3) * ΔE + (r1 - r2) * sn / (1 - hr * sn^2) *
                           cn * sqrt(1 - kr * sn^2)))
        return (; radial(qr)..., Δτr)
    end
end

"""
    geodesic_trajectory_radial(a, p, e, x) -> Function

Radial worldline of a generic (inclined) bound Kerr geodesic (M = 1). Returns
`qr -> (; r, Ur, Δtr, Δφr)`: BL radius and the oscillating radial deviations of t, φ (the
KerrGeoOrbit "Phases" convention). Replaces `KerrGeoTrajectoryRadial`.
"""
function geodesic_trajectory_radial(a, p, e, x)
    if a < 0 && x^2 == 1
        # Equatorial signed-spin chart: consistent with the frequency routine.
        radial = geodesic_trajectory_radial(-a, p, e, -x)
        return function (qr)
            g = radial(qr)
            return (r=g.r, Ur=g.Ur, Δtr=g.Δtr, Δφr=-g.Δφr)
        end
    end
    En, L, Q = constants_of_motion(a, p, e, x)
    r1c, r2c, r3c, r4c = radial_roots(a, p, e, x, En, Q)
    r1 = real(r1c); r2 = real(r2c); r3 = real(r3c); r4 = real(r4c)
    rp = 1 + sqrt(1 - a^2); rm = 1 - sqrt(1 - a^2)
    kr = _erc((r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4))
    hr = (r1 - r2) / (r1 - r3)
    hp = ((r1 - r2) * (r3 - rp)) / ((r1 - r3) * (r2 - rp))
    hm = ((r1 - r2) * (r3 - rm)) / ((r1 - r3) * (r2 - rm))
    Kk = real(ellipK(kr)); Ek = real(ellipE(kr))
    Πhr = real(ellipPi(hr, kr)); Πhp = real(ellipPi(hp, kr)); Πhm = real(ellipPi(hm, kr))
    sden = sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
    return function (qr)
        u = Kk / pi * qr
        ψr = jacobiAM(u, kr); sn = sin(ψr)
        ΔE  = Ek  * qr / pi - ellipEinc(ψr, kr)
        ΔΠhr = Πhr * qr / pi - ellipPiinc(hr, ψr, kr)
        ΔΠhm = Πhm * qr / pi - ellipPiinc(hm, ψr, kr)
        ΔΠhp = Πhp * qr / pi - ellipPiinc(hp, ψr, kr)
        r = (r3 * (r1 - r2) * sn^2 - r2 * (r1 - r3)) / ((r1 - r2) * sn^2 - (r1 - r3))
        # Differentiate r(qr) analytically instead of subtracting large terms in R(r).
        # This is signed on both halves and vanishes at the turning points.
        Ur = sden * hr * (r2 - r3) * sn * cos(ψr) * sqrt(1 - kr * sn^2) /
             (1 - hr * sn^2)^2
        bnd = hr * (sin(ψr) * cos(ψr) * sqrt(1 - kr * sin(ψr)^2)) / (1 - hr * sin(ψr)^2)
        Δtr = -En / sden * (
            4 * (r2 - r3) * ΔΠhr
            - 4 * (r2 - r3) / (rp - rm) * (
              -(1 / ((-rm + r2) * (-rm + r3))) * (-2a^2 + rm * (4 - a * L / En)) * ΔΠhm
              +(1 / ((-rp + r2) * (-rp + r3))) * (-2a^2 + rp * (4 - a * L / En)) * ΔΠhp)
            + (r2 - r3) * (r1 + r2 + r3 + r4) * ΔΠhr
            + (r1 - r3) * (r2 - r4) * (ΔE + bnd))
        Δφr = (2a * En * (
              -(1 / ((-rm + r2) * (-rm + r3))) * (2rm - a * L / En) * (r2 - r3) * ΔΠhm
              +(1 / ((-rp + r2) * (-rp + r3))) * (2rp - a * L / En) * (r2 - r3) * ΔΠhp)) /
              ((-rm + rp) * sden)
        return (r = r, Ur = Ur, Δtr = Δtr, Δφr = Δφr)
    end
end

"""
    geodesic_trajectory_polar(a, p, e, x) -> Function

Polar worldline of a bound Kerr geodesic (M = 1). Returns `qz -> (; z, θ, Uz, Δtz, Δφz)`:
`z = cosθ` via Jacobi sn, the BL polar angle, and the oscillating polar deviations of t, φ.
Equatorial (x² = 1) ⇒ z ≡ 0, θ ≡ π/2, Δtz = Δφz = 0. Replaces `KerrGeoTrajectoryPolar`.
"""
function geodesic_trajectory_polar(a, p, e, x)
    En, L, Q = constants_of_motion(a, p, e, x)
    Kc = Q + (L - a * En)^2                       # Carter K = Q + (Lz − aE)²
    zp, zm = polar_roots(a, p, e, x)
    kθ = _erc(a^2 * (1 - En^2) * (zm / zp)^2)
    Kkθ = real(ellipK(kθ)); Ekθ = real(ellipE(kθ)); Πkθ = real(ellipPi(_erc(zm^2), kθ))
    return function (qz)
        u = (2 / pi) * Kkθ * (qz + pi / 2)        # sn argument: qz=0 ⇒ u=K(kθ) ⇒ z=zm (turning point)
        f = 2 * (qz + pi / 2) / pi                 # number of K(kθ)'s traversed (complete-integral weight)
        z = zm * jacobiSN(u, kθ)
        ψz = jacobiAM(u, kθ)
        Zz = -((1 - z^2) * a * En - L)^2 + (1 - z^2) * (Kc - a^2 * z^2)   # polar potential Z(z)
        Uz = (mod(qz, 2pi) < pi ? -1 : 1) * sqrt(max(Zz, zero(Zz)))       # signed dz/dλ
        Δtz = (En * zp / (1 - En^2)) * (Ekθ * f - ellipEinc(ψz, kθ))
        Δφz = (-L / zp) * (Πkθ * f - ellipPiinc(_erc(zm^2), ψz, kθ))
        return (z = z, θ = acos(z), Uz = Uz, Δtz = Δtz, Δφz = Δφz)
    end
end

"""
    geodesic_trajectory(a, p, e, x) -> NamedTuple of functions

Full generic worldline of a bound Kerr geodesic (M = 1), the closed-form Mino-phase
"Phases" representation. Returns `(; t, r, θ, φ)` with `t(qt,qr,qz)=qt+Δtr(qr)+Δtz(qz)`,
`r(qr)`, `θ(qz)`, `φ(qφ,qr,qz)=qφ+Δφr(qr)+Δφz(qz)`. Assembles the radial + polar halves.
Replaces `KerrGeoTrajectory`.
"""
function geodesic_trajectory(a, p, e, x)
    rad = geodesic_trajectory_radial(a, p, e, x)
    pol = geodesic_trajectory_polar(a, p, e, x)
    return (t = (qt, qr, qz) -> qt + rad(qr).Δtr + pol(qz).Δtz,
            r = qr -> rad(qr).r,
            θ = qz -> pol(qz).θ,
            φ = (qφ, qr, qz) -> qφ + rad(qr).Δφr + pol(qz).Δφz)
end

"""
    geodesic_proper_time(a, p, e, x) -> NamedTuple

Proper-time structure of a generic (inclined) bound Kerr geodesic (M = 1). Returns
`(; ProperTimeFrequency=Υτ, Υτz, Δτr, Δτz, ProperTime)` with
`τ(qτ,qr,qz)=qτ+Δτr(qr)+Δτz(qz)` (dτ/dλ = r²+a²z² separates). Inclined-only: degenerate at
x²=1 (z1=0); equatorial proper time is `geodesic_trajectory_equatorial`'s Δτr. Replaces
`KerrGeoProperTime` (crosswalk row 350).
"""
function geodesic_proper_time(a, p, e, x)
    En, L, Q = constants_of_motion(a, p, e, x)
    K = Q + (L - a * En)^2
    r1c, r2c, r3c, r4c = radial_roots(a, p, e, x, En, Q)
    r1 = real(r1c); r2 = real(r2c); r3 = real(r3c); r4 = real(r4c)
    z1 = sqrt(1 - x^2)
    z2 = sqrt(Q / (a^2 * (1 - En^2) * z1^2))                  # K − (Lz−aE)² = Q
    krsq = _erc((r1 - r2) * (r3 - r4) / ((r1 - r3) * (r2 - r4)))
    kzsq = _erc(z1^2 / z2^2)
    Kkr = real(ellipK(krsq)); Kkz = real(ellipK(kzsq))
    hr = (r1 - r2) / (r1 - r3)
    Υτr = 0.5 * ((r1 * (-r2 + r3) + r3 * (r2 + r3))
          + (r2 - r3) * (r1 + r2 + r3 + r4) * real(ellipPi(hr, krsq)) / Kkr
          + (r1 - r3) * (r2 - r4) * real(ellipE(krsq)) / Kkr)
    Υτz = a^2 * z2^2 * (1 - real(ellipE(kzsq)) / Kkz)
    Υτ = Υτr + Υτz
    sdr = sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
    τhatz(ξz) = -z2 * a / sqrt(1 - En^2) * ellipEinc(ξz, kzsq)
    τhatr(ξr) = ((r2 - r3) * (r1 + r2 + r3 + r4) * ellipPiinc(hr, ξr, krsq)
                 + (r1 - r3) * (r2 - r4) * ellipEinc(ξr, krsq)
                 - (r1 - r2) * (r2 - r4) * sin(ξr) * cos(ξr) * sqrt(1 - krsq * sin(ξr)^2) /
                   (1 - hr * sin(ξr)^2)) / sdr
    τhatrπ = τhatr(pi); τhatzπ = τhatz(pi)
    Δτr(qr) = τhatr(jacobiAM(Kkr * qr / pi, krsq)) - τhatrπ / (2pi) * qr
    Δτz(qz) = τhatz(jacobiAM(2 * Kkz * (qz + pi / 2) / pi, kzsq)) - τhatzπ / pi * (qz + pi / 2)
    return (ProperTimeFrequency = Υτ, Υτz = Υτz, Δτr = Δτr, Δτz = Δτz,
            ProperTime = (qτ, qr, qz) -> qτ + Δτr(qr) + Δτz(qz))
end

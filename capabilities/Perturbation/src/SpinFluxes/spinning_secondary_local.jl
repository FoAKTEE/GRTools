# LOCAL spinning-secondary kinematics for SpinFluxes — the closed-form spinning worldline
# (KerrSpinOrbit, equatorial + INCLINED) and the generic linear-in-spin orbital-parameter shift
# pex1CTilde. Faithful 1:1 of SpinningOrbit/src/{linear_parts,spin_orbit,spin_deviations}.jl
# (SpinningTrajectory.wl KerrSpinOrbit L320-398, pex1CTildeFun, L272-275), rebuilt on Spacetime's
# geodesic layer + NumericsCore's elliptic layer.
#
# FIREWALL NOTE (design §4.3/§7.7): Spacetime.SpinningSecondary carries only the EQUATORIAL slice
# (spin_orbit throws for x²≠1); the inclined/spherical/numeric-correction kinematics the
# TeukolskySpinMode* entrypoints need are therefore implemented LOCALLY here — Perturbation gains
# NO new [deps] (pure arithmetic + the existing Spacetime/NumericsCore closure). These stay
# unexported `_`-conceptual internals of Perturbation.SpinFluxes. M = 1, bound orbits.

# shared linear-in-spin constants-of-motion shift (En,Jz,K) ← (EnT,JzT,KT), SpinningTrajectory.wl L339-341
function _kso_spin_shifted_constants(EnT, JzT, KT, a, spar)
    sK = sqrt(KT)
    En = EnT - spar * (1 - EnT^2) / (2 * sK)
    Jz = JzT - spar * (a - EnT * JzT / 2) / sK
    K  = KT  - spar * (3a * (JzT - a * EnT) - EnT * KT) / sK
    return (En, Jz, K, sK)
end

"""
    _kso_pex1_ctilde(a, p, e, x) -> (Δp, Δe, Δx)

Linear-in-spin correction (coefficient of `s∥`) to the orbital parameters `(p, e, x)` in the
CTilde gauge — the GENERIC `pex1CTildeFun` (SpinningTrajectory.wl), equatorial + inclined.
Equatorial (x²=1): z1=0 ⟹ intΣr=1/r², no elliptic integrals, Δx=0. Inclined: the `az1IntΣ` term
needs a complex-characteristic `EllipticPi` (NumericsCore elliptic layer). Local port of
`SpinningOrbit.linear_orbital_params_CTilde` (Spacetime ships only the equatorial branch).
"""
function _kso_pex1_ctilde(a, p, e, x)
    En, Lz, Q = constants_of_motion(a, p, e, x)
    K = Q + (Lz - a * En)^2
    sK = sqrt(K); r1 = p / (1 - e); r2 = p / (1 + e)
    if x^2 == 1                                              # equatorial: z1=0 ⟹ intΣr=1/r², Δx=0
        dr1 = -((r1^2 + a^2) * En - a * Lz) / (r1 * sK)
        dr2 = -((r2^2 + a^2) * En - a * Lz) / (r2 * sK)
        return ((2 * r2^2 * dr1 + 2 * r1^2 * dr2) / (r1 + r2)^2,
                (2 * r2 * dr1 - 2 * r1 * dr2) / (r1 + r2)^2, 0.0)
    end
    # inclined: full pex1CTildeFun with elliptic integrals (NumericsCore elliptic layer)
    _, _, r3, r4 = radial_roots(a, p, e, x, En, Q)           # r3,r4 (Complex carrier; real for bound)
    z1 = sqrt(1 - x^2)
    az2 = sqrt(a^2 * (1 - z1^2) + (Lz^2 + Q) / (1 - En^2))
    kzsq = a^2 * z1^2 / az2^2
    krsq = real((r1 - r2) * (r3 - r4) / ((r1 - r3) * (r2 - r4)))    # mathematically real for bound orbits
    intSr1 = ellipPi(-a^2 * z1^2 / r1^2, kzsq) / ellipK(kzsq) / r1^2
    intSr2 = ellipPi(-a^2 * z1^2 / r2^2, kzsq) / ellipK(kzsq) / r2^2
    nC = ((r1 - r2) * (r3 + im * a * z1)) / ((r1 - r3) * (r2 + im * a * z1))   # complex characteristic
    prefac = (im * (-r2 + r3)) / (2 * (r2 + im * a * z1) * (r3 + im * a * z1))
    az1IntS = real(a * z1 / (r3^2 + a^2 * z1^2)) +
              2 * real(prefac * ellipPi(nC, krsq)) / real(ellipK(krsq))
    dr1 = -r1 * ((r1^2 + a^2) * En - a * Lz) * intSr1 / sK
    dr2 = -r2 * ((r2^2 + a^2) * En - a * Lz) * intSr2 / sK
    dz1 = -(Lz - a * En * x^2) * az1IntS / sK
    dp = (2 * r2^2 * dr1) / (r1 + r2)^2 + (2 * r1^2 * dr2) / (r1 + r2)^2
    de = (2 * r2 * dr1) / (r1 + r2)^2 - (2 * r1 * dr2) / (r1 + r2)^2
    dx = -(z1 * dz1) / x
    return (real(dp), real(de), real(dx))
end

# Linear-in-spin worldline deviations (spin-curvature coupling) δt3/δr3/δz3/δφ3 — faithful 1:1 of
# SpinningTrajectory.wl L272-275 (local port of SpinningOrbit/src/spin_deviations.jl).
_kso_δt3(r, z, Ur, Uz, K, a) =
    (r * (a^2 + r^2) * Ur - a^2 * (a^2 + r * (-2 + r)) * z * Uz) /
    (sqrt(K) * (a^2 + r * (-2 + r)) * (r^2 + a^2 * z^2))
_kso_δr3(r, z, En, Jz, K, a) =
    (r * (a^2 * En - a * Jz + En * r^2)) / (sqrt(K) * (r^2 + a^2 * z^2))
_kso_δz3(r, z, En, Jz, K, a) =
    (a * z * (Jz - a * En * (1 - z^2))) / (sqrt(K) * (r^2 + a^2 * z^2))
_kso_δφ3(r, z, Ur, Uz, K, a) =
    (a * (r * (1 - z^2) * Ur - (a^2 + r * (-2 + r)) * z * Uz)) /
    (sqrt(K) * (a^2 + r * (-2 + r)) * (1 - z^2) * (r^2 + a^2 * z^2))

"""
    KerrSpinOrbit(a, p, e, x, spar; gauge="pex") -> NamedTuple

Trajectory + scalars of a spinning secondary on a Kerr orbit, linear in the secondary spin `spar`,
EQUATORIAL + INCLINED (M = 1). Local 1:1 port of `SpinningOrbit.KerrSpinOrbit`
(SpinningTrajectory.wl L320-398) on the Spacetime geodesic layer — the inclined branch feeds
`TeukolskySpinModeAnalytical`. `gauge ∈ {"pex","CTilde"}`. Unexported (internal kinematics).
"""
function KerrSpinOrbit(a, p, e, x, spar; gauge::AbstractString="pex")
    if gauge == "pex"
        δpex = _kso_pex1_ctilde(a, p, e, x)
        pT = p - spar * δpex[1]; eT = e - spar * δpex[2]; xT = x - spar * δpex[3]
    elseif gauge == "CTilde"
        pT = p; eT = e; xT = x
    else
        throw(ArgumentError("KerrSpinOrbit: unknown gauge $(repr(gauge)) (use \"pex\" or \"CTilde\")"))
    end
    if x^2 == 1
        xT = sign(xT)                                        # equatorial stays exactly ±1
        EnT = geodesic_energy(a, pT, eT, xT); JzT = geodesic_angular_momentum(a, pT, eT, xT)
        KT  = carter_constant(a, pT, eT, xT) + (JzT - a * EnT)^2
        En, Jz, K, sK = _kso_spin_shifted_constants(EnT, JzT, KT, a, spar)
        fr  = kerr_frequencies(a, pT, eT, xT; Time="Mino")
        Uτ  = proper_frequency_factor(a, pT, eT, xT)
        Υr = fr.Ur; Υz = fr.Uth; Υφ = fr.Uph
        Υt = fr.Ut - 3 * spar / (2 * sK) * Uτ                # the only nonzero Mino freq correction
        geo = geodesic_trajectory_equatorial(a, pT, eT)
        traj = function (qr)
            g = geo(qr)
            return (r = g.r, Δtr = g.Δtr - 3 * spar / (2 * sK) * g.Δτr, Δφr = g.Δφr)
        end
        return (a = float(a), p = float(p), e = float(e), x = float(x), spar = float(spar),
                Energy = En, AngularMomentum = Jz, CarterConstantK = K,
                EnTilde = EnT, JzTilde = JzT, KTilde = KT, pTilde = pT, eTilde = eT, xTilde = float(xT),
                Frequencies = (Ur = Υr, Uz = Υz, Uph = Υφ, Ut = Υt),
                BLFrequencies = (Υr / Υt, Υz / Υt, Υφ / Υt),
                ProperTimeFrequency = Uτ, trajectory = traj)
    end
    return _kerr_spin_orbit_inclined(a, p, e, x, spar, pT, eT, xT)
end

"""
    _kerr_spin_orbit_inclined(a, p, e, x, spar, pT, eT, xT) -> NamedTuple

Inclined (x²<1) closed-form spinning orbit — SpinningTrajectory.wl `KerrSpinOrbit` L320-398.
Adds the polar sector to the equatorial path: polar Mino frequencies Υtz/Υφz, the spin-corrected
polar/radial t-deltas, the geodesic four-velocities Ur(qr)/Uz(qz), and the spin-curvature
deviation field δx(qr,qz)={δt3,δr3,δz3,δφ3}.
"""
function _kerr_spin_orbit_inclined(a, p, e, x, spar, pT, eT, xT)
    EnT = geodesic_energy(a, pT, eT, xT); JzT = geodesic_angular_momentum(a, pT, eT, xT)
    QT  = carter_constant(a, pT, eT, xT)
    KT  = QT + (JzT - a * EnT)^2
    En, Jz, K, sK = _kso_spin_shifted_constants(EnT, JzT, KT, a, spar)
    fr  = kerr_frequencies(a, pT, eT, xT; Time="Mino")       # Υr, Υθ, Υφ, Υt (geodesic, tilde)
    pt  = geodesic_proper_time(a, pT, eT, xT)                # Υτ, Υτz, Δτr, Δτz, τ
    Uτ  = pt.ProperTimeFrequency
    zp, zm = polar_roots(a, pT, eT, xT)
    Υr = fr.Ur; Υz = fr.Uth; Υφ = fr.Uph
    Υt  = fr.Ut - 3 * spar / (2 * sK) * Uτ
    Υtz = Spacetime._mino_freq_t_theta(a, EnT, QT, zp, zm, xT) + a * JzT -
          3 * spar / (2 * sK) * pt.Υτz
    Υφz = Spacetime._mino_freq_phi_theta(a, EnT, JzT, zp, zm) - a * EnT
    rad = geodesic_trajectory_radial(a, pT, eT, xT)          # r, Ur, Δtr, Δφr
    pol = geodesic_trajectory_polar(a, pT, eT, xT)           # z, θ, Uz, Δtz, Δφz
    c3 = 3 * spar / (2 * sK)
    Δtr(qr) = rad(qr).Δtr - c3 * pt.Δτr(qr)
    Δtz(qz) = pol(qz).Δtz - c3 * pt.Δτz(qz)
    Δφr(qr) = rad(qr).Δφr
    Δφz(qz) = pol(qz).Δφz
    deltax = function (qr, qz)
        g = rad(qr); h = pol(qz); r = g.r; z = h.z; Ur = g.Ur; Uz = h.Uz
        return (δt3 = _kso_δt3(r, z, Ur, Uz, KT, a), δr3 = _kso_δr3(r, z, EnT, JzT, KT, a),
                δz3 = _kso_δz3(r, z, EnT, JzT, KT, a), δφ3 = _kso_δφ3(r, z, Ur, Uz, KT, a))
    end
    return (a = float(a), p = float(p), e = float(e), x = float(x), spar = float(spar),
            Energy = En, AngularMomentum = Jz, CarterConstantK = K,
            EnTilde = EnT, JzTilde = JzT, KTilde = KT, pTilde = pT, eTilde = eT, xTilde = float(xT),
            Frequencies = (Ur = Υr, Uz = Υz, Uph = Υφ, Ut = Υt),
            Upsilontz = Υtz, Upsilonphiz = Υφz,
            BLFrequencies = (Υr / Υt, Υz / Υt, Υφ / Υt), ProperTimeFrequency = Uτ,
            TrajectoryDeltas = (Δtr = Δtr, Δtz = Δtz, Δφr = Δφr, Δφz = Δφz),
            trajectory = (qr -> rad(qr).r, qz -> pol(qz).z),
            FourVelocities = (qr -> rad(qr).Ur, qz -> pol(qz).Uz),
            deltax = deltax)
end

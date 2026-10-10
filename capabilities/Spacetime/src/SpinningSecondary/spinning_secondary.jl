"""
    Spacetime.SpinningSecondary

Linear-in-spin (Mathisson–Papapetrou–Dixon) secondary-orbit corrections to a Kerr geodesic
(M = 1), per design §4.1: linear-in-spin constants/frequencies, the closed-form spinning
worldline in the EQUATORIAL case, AND the ~60 spin-curvature force terms with the numeric
orbit-correction machinery they drive — the 4 worldline deviations `_δt3/_δr3/_δz3/_δφ3`, the
47 force/coefficient terms + 2-torus DFT `Correction` engine (`_numeric_correction`), and the
11 `_oc_*` four-velocity/metric helpers + orbit-correction assembly (`spin_orbit_correction`,
`spin_orbit_numeric`), all for GENERIC (inclined) bound orbits. Consolidates
`SpinningOrbit/src/{linear_parts,spin_orbit,spin_deviations,numeric_correction,orbit_correction}.jl`
(`KerrSpinLinearParts`, `KerrSpinOrbit`, `KerrSpinOrbitEquatorial`, `KerrSpinOrbitCorrection`,
`KerrSpinOrbitNum`). Namespaced as `Spacetime.SpinningSecondary.*` to avoid the root-package
`SpinningOrbit` module-name collision (design §4.1).

NOT in this module ([FUTURE], scoped out of the S5-1 force-term slice — audit 2026-07: the legacy
oracle NOW implements these too, so this is a scope boundary, not an oracle-absent stub): the
INCLINED closed-form `spin_orbit` branch (x² ≠ 1, `_kerr_spin_orbit_inclined`) and the SPHERICAL
(e = 0, θ-only) engine (`KerrSpinOrbit[Correction]Spherical` / `_orbit_derivatives_spherical`).
NOTE the numeric `Correction` engine here IS the generic inclined force machinery (evaluated at
x² < 1); only the closed-form-inclined/spherical assemblies stay out.

Type-generic `{T}` for all closed-form/algebraic parts (BigFloat inherited); the DFT/least-squares
engine is Float64-tier (LAPACK `pinv`, faithful to the oracle) and serial (no threaded LAPACK
region ⇒ no §8 BLAS-guard required).
"""
module SpinningSecondary

using LinearAlgebra: pinv, norm     # M-matrix least-squares solve in the numeric correction engine

using NumericsCore: ellipK, ellipE, ellipPi, RD    # action-Jacobian elliptic layer (Mathematica m=k²)

using ..Spacetime: constants_of_motion, geodesic_energy, geodesic_angular_momentum, carter_constant,
                   kerr_frequencies, proper_frequency_factor, geodesic_trajectory_equatorial,
                   geodesic_trajectory_radial, geodesic_trajectory_polar, radial_roots, orbital_params_jacobian

export spin_linear_constants, spin_linear_frequencies, spin_mino_time_correction,
       spin_linear_orbital_params,
       spin_linear_actions, spin_actions_pex, spin_frequencies_action_gauge,
       spin_frequencies_constants_gauge, spin_frequencies_turning_point_gauge,
       spin_action_jacobian, spin_frequency_jacobian, spin_orbital_params_jacobian,
       spin_orbit, spin_orbit_equatorial,
       spin_orbit_correction, spin_orbit_numeric, shifted_geodesic_correction

include("linear_parts.jl")
include("action_parts.jl")        # linear-in-spin action parts + the J / pex gauge transforms
include("spin_deviations.jl")     # spin-curvature worldline deviations δt3/δr3/δz3/δφ3
include("numeric_correction.jl")  # 47 force terms + Darwin-angle maps + 2-torus DFT + M-matrix solve
include("orbit_correction.jl")    # 11 _oc_* helpers + spin_orbit_correction/_numeric assembly

# shared linear-in-spin constants-of-motion shift (En,Jz,K) ← (EnT,JzT,KT) (SpinningTrajectory.wl L339-341)
function _spin_shifted_constants(EnT, JzT, KT, a, spar)
    sK = sqrt(KT)
    En = EnT - spar * (1 - EnT^2) / (2 * sK)
    Jz = JzT - spar * (a - EnT * JzT / 2) / sK
    K  = KT  - spar * (3a * (JzT - a * EnT) - EnT * KT) / sK
    return (En, Jz, K, sK)
end

"""
    spin_orbit(a, p, e, x, spar; gauge="pex") -> NamedTuple

Trajectory + scalars of a spinning secondary on a Kerr orbit, linear in the secondary spin `spar`,
EQUATORIAL (x² = 1, M = 1). Returns `(; a,p,e,x,spar, Energy,AngularMomentum,CarterConstantK,
EnTilde,JzTilde,KTilde,pTilde,eTilde,xTilde, Frequencies=(Ur,Uz,Uph,Ut), BLFrequencies,
ProperTimeFrequency, trajectory)`, where `trajectory(qr) = (; r, Δtr, Δφr)` with the analytic spin
shift `Δtr = Δtr_geo − 3 spar/(2√K̃)·Δτr` folded in. `gauge ∈ {"pex","CTilde"}`. Inclined (x²≠1) is
outside the equatorial-first slice. Replaces `SpinningOrbit.KerrSpinOrbit` (equatorial branch).
"""
function spin_orbit(a, p, e, x, spar; gauge::AbstractString = "pex")
    x^2 == 1 || throw(ArgumentError(
        "spin_orbit: only the equatorial case x²=1 is in this slice (got x=$x); " *
        "inclined/spherical spinning orbits are NotYetPorted."))
    if gauge == "pex"
        δp, δe, _ = _spin_orbital_params_eq(a, p, e)
        pT = p - spar * δp; eT = e - spar * δe; xT = sign(x)
    elseif gauge == "CTilde"
        pT = p; eT = e; xT = sign(x)
    else
        throw(ArgumentError("spin_orbit: unknown gauge $(repr(gauge)) (use \"pex\" or \"CTilde\")"))
    end
    EnT = geodesic_energy(a, pT, eT, xT); JzT = geodesic_angular_momentum(a, pT, eT, xT)
    KT  = carter_constant(a, pT, eT, xT) + (JzT - a * EnT)^2
    En, Jz, K, sK = _spin_shifted_constants(EnT, JzT, KT, a, spar)
    fr = kerr_frequencies(a, pT, eT, xT; Time = "Mino")
    Uτ = proper_frequency_factor(a, pT, eT, xT)
    Υr = fr.Ur; Υz = fr.Uth; Υφ = fr.Uph
    Υt = fr.Ut - 3 * spar / (2 * sK) * Uτ                       # the only nonzero Mino-freq correction
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

"""
    spin_orbit_equatorial(a, p, e) -> NamedTuple

Geodesic-limit (spin = 0) bound equatorial Kerr orbit (M = 1) in the spinning-orbit assembly
convention. Returns `(; a,p,e,Inclination=1, Energy,AngularMomentum,CarterConstant,
Frequencies=(Ur,Uth,Uph,Ut), ProperTimeFrequency, trajectory)` with `trajectory(qr) =
(; r,Ur,Δtr,Δφr,Δτr)`. Replaces `SpinningOrbit.KerrSpinOrbitEquatorial`.
"""
function spin_orbit_equatorial(a, p, e)
    fr = kerr_frequencies(a, p, e, 1; Time = "Mino")
    return (a = float(a), p = float(p), e = float(e), Inclination = 1,
            Energy = geodesic_energy(a, p, e, 1),
            AngularMomentum = geodesic_angular_momentum(a, p, e, 1),
            CarterConstant = carter_constant(a, p, e, 1),
            Frequencies = (Ur = fr.Ur, Uth = fr.Uth, Uph = fr.Uph, Ut = fr.Ut),
            ProperTimeFrequency = proper_frequency_factor(a, p, e, 1),
            trajectory = geodesic_trajectory_equatorial(a, p, e))
end

end # module SpinningSecondary

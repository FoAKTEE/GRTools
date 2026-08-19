# Linear-in-spin (MPD) corrections in the "CTilde" (fixed shifted-constants) gauge — the
# foundational closed-form parts of a spinning secondary. Consolidates SpinningOrbit/src/
# linear_parts.jl (C1CTildeFun, Ω1CTildeFun). Paper: arXiv:2603.13482 (fixed-Ċ gauge). M = 1.
# K = Q + (Lz − aE)² is the squared Carter–Killing constant. Pure arithmetic + the in-package
# Geodesics frequency layer ⇒ BigFloat inherited.

"`_ctilde_K(a,p,e,x)` → `(En, Lz, Q, K)` with `K = Q + (Lz − aEn)²` (M = 1)."
function _ctilde_K(a, p, e, x)
    En, Lz, Q = constants_of_motion(a, p, e, x)
    return En, Lz, Q, Q + (Lz - a * En)^2
end

"""
    spin_linear_constants(a, p, e, x) -> (δE, δLz, δK)

Linear-in-spin correction (coefficient of the aligned spin `s∥`) to the constants of motion in
the CTilde gauge: `−{(1−E²)/2, a − Lz E/2, 3a(Lz − aE) − K E} / √K` (SpinningTrajectory.wl
`C1CTildeFun`). Replaces `SpinningOrbit.linear_constants_CTilde` (`KerrSpinLinearParts.C1CTilde`).
"""
function spin_linear_constants(a, p, e, x)
    En, Lz, _, K = _ctilde_K(a, p, e, x)
    sK = sqrt(K)
    return (-(1 - En^2) / 2 / sK,
            -(a - Lz * En / 2) / sK,
            -(3a * (Lz - a * En) - K * En) / sK)
end

"""
    spin_linear_frequencies(a, p, e, x) -> (δΩr, δΩθ, δΩφ)

Linear-in-spin correction (coefficient of `s∥`) to the Boyer–Lindquist fundamental frequencies
in the CTilde gauge: `(Ωr, Ωθ, Ωφ) · 3 Υτ / (2 Υt √K)` (SpinningTrajectory.wl `Ω1CTildeFun`),
Υτ = proper-frequency factor, Υt = Mino-time `Υ_t`. Replaces `SpinningOrbit.linear_frequencies_CTilde`
(`KerrSpinLinearParts.Omega1CTilde`).
"""
function spin_linear_frequencies(a, p, e, x)
    _, _, _, K = _ctilde_K(a, p, e, x)
    Uptau = proper_frequency_factor(a, p, e, x)                     # Υτ
    Upt   = kerr_frequencies(a, p, e, x; Time = "Mino").Ut          # Υt
    Ω     = kerr_frequencies(a, p, e, x)                            # (Omegar, Omegath, Omegaph)
    f = 3 * Uptau / (2 * Upt * sqrt(K))
    return (Ω.Omegar * f, Ω.Omegath * f, Ω.Omegaph * f)
end

"""
    spin_mino_time_correction(a, p, e, x) -> δΥt

Linear-in-spin correction (coefficient of `s∥`) to the Mino-time frequency factor `Υt` in the
CTilde gauge: `Υt¹ = −3 Υτ / (2 √K)` — the only nonzero Mino-frequency correction (Υr¹=Υθ¹=Υφ¹=0);
the BL corrections satisfy `Ω_k¹ = −Ω_k · Υt¹/Υt`. Replaces `SpinningOrbit.linear_mino_time_CTilde`.
"""
function spin_mino_time_correction(a, p, e, x)
    _, _, _, K = _ctilde_K(a, p, e, x)
    return -3 * proper_frequency_factor(a, p, e, x) / (2 * sqrt(K))
end

"""
    spin_linear_orbital_params(a, p, e, x) -> (δp, δe, δx)

Linear-in-spin correction (coefficient of the aligned spin `s∥`) to the orbital parameters
`(p, e, x)` in the CTilde gauge — the GENERIC `pex1CTildeFun` (SpinningTrajectory.wl), equatorial
AND inclined. Equatorial (`x² = 1`): `z1 = 0` ⟹ `intΣr = 1/r²`, no elliptic integrals, `δx = 0`.
Inclined: the `az1IntΣ` term needs a complex-characteristic `EllipticPi` (NumericsCore elliptic
layer). Feeds the "pex" (fixed turning points) gauge of [`spin_orbit`](@ref) and the action
transforms in [`spin_actions_pex`](@ref). Replaces `SpinningOrbit.linear_orbital_params_CTilde`
(`KerrSpinLinearParts.pex1CTilde`).
"""
function spin_linear_orbital_params(a, p, e, x)
    En, Lz, Q, K = _ctilde_K(a, p, e, x)
    sK = sqrt(K); r1 = p / (1 - e); r2 = p / (1 + e)
    if x^2 == 1                                    # equatorial: z1=0 ⟹ intΣr = 1/r², δx = 0
        dr1 = -((r1^2 + a^2) * En - a * Lz) / (r1 * sK)
        dr2 = -((r2^2 + a^2) * En - a * Lz) / (r2 * sK)
        return ((2 * r2^2 * dr1 + 2 * r1^2 * dr2) / (r1 + r2)^2,
                (2 * r2 * dr1 - 2 * r1 * dr2) / (r1 + r2)^2,
                zero(sK))
    end
    _, _, r3, r4 = radial_roots(a, p, e, x, En, Q)          # r3,r4 (Complex carrier; real for bound)
    z1  = sqrt(1 - x^2)
    az2 = sqrt(a^2 * (1 - z1^2) + (Lz^2 + Q) / (1 - En^2))
    kzsq = a^2 * z1^2 / az2^2
    krsq = real((r1 - r2) * (r3 - r4) / ((r1 - r3) * (r2 - r4)))   # real for bound orbits
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

# equatorial-only alias kept for the closed-form `spin_orbit` call site (x² = 1 branch).
_spin_orbital_params_eq(a, p, e) = spin_linear_orbital_params(a, p, e, one(float(p)))

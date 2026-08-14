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

# Equatorial (x²=1) linear-in-spin orbital-parameter correction (δp, δe, 0): z1=0 ⟹ intΣr = 1/r²,
# no elliptic integrals, δx=0 (SpinningTrajectory.wl pex1CTildeFun, equatorial branch). Feeds the
# "pex" gauge of spin_orbit. The inclined elliptic branch is out of the equatorial-first slice.
function _spin_orbital_params_eq(a, p, e)
    En, Lz, _, K = _ctilde_K(a, p, e, 1)
    sK = sqrt(K); r1 = p / (1 - e); r2 = p / (1 + e)
    dr1 = -((r1^2 + a^2) * En - a * Lz) / (r1 * sK)
    dr2 = -((r2^2 + a^2) * En - a * Lz) / (r2 * sK)
    return ((2 * r2^2 * dr1 + 2 * r1^2 * dr2) / (r1 + r2)^2,
            (2 * r2 * dr1 - 2 * r1 * dr2) / (r1 + r2)^2,
            zero(sK))
end

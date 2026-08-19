# Modes/PressurelessLindblad.jl — the ANALYTIC pressureless relativistic Lindblad torque of
# Hegade K. R., Gammie & Yunes 2025, "Paper II" (arXiv:2510.03564, Kerr; PRD 112, 124068), whose
# Schwarzschild special case is "Paper I" (arXiv:2509.20457).  This is the pressureless comparison
# BASELINE that arXiv:2608.17003 (the SpiralWaves linear theory) tests in its Fig. 9: 2608 main.tex
# :984 states verbatim "The exact analytical result in Kerr spacetime can be obtained from
# Eqs. (18 b) and Eq. (20) of [HegadeKR:2025rpr]", and the per-mode / summed forms are 2608 (77),
# (78) with the cutoffs (79) m_max ≈ 1/h0 and (80) m_45 ≈ 0.35/h0.
#
# NOT a duplicate of `Modes.lindblad_torque` (Forcing.jl): that one is the Hirata 1010.0758/0759
# Teukolsky GW-overlap torque of a single (m, ILR/OLR) resonance evaluated through the `Perturbation`
# MST engine.  This file is the CLOSED-FORM secular Hamiltonian-perturbation-theory result of the
# disturbing function summed over resonances, with a torque cutoff — a different kernel, a different
# evidence chain, one exported name each.
#
# Symbol map (their notation → ours): their mode number j ≡ 2608's m; their j_max ≡ 2608's m_max;
# their p' (dimensionless orbital semi-latus rectum, circular ⇒ orbital radius in M) ≡ 2608's a•/M.
# M = 1 geometric units, χ = a/M (χ > 0 prograde, χ < 0 retrograde), q = m•/M the mass ratio.
#
# REUSE (no second derivation of a frequency): the printed auxiliary polynomials are EXACTLY the
# existing Kerr primitives — with x ≡ (p')^{1/4},
#     (χ + x⁶)              = 1/Ω          `Modes.disc_OmegaK`
#     (−3χ² + x⁸ − 6x⁴ + 8χx²) = p² κ²/Ω²  `Modes.kappar2` / `Modes.disc_OmegaK`
#     (χ² + x⁸ − 2x⁴)       = Δ            `Spacetime.delta_kerr` (via `Modes.delta_r`)
#     (2χ + x⁶ − 3x²)       = 1/(p^{3/2} Ω² (u^t)²)   `Spacetime.kerr_circular_dtdtau`
# so this file DEFINES them through those primitives and the test file asserts the printed-polynomial
# identity.  Every algebraic primitive is type-generic `{T}` (BigFloat-inherited); the modified
# Bessel kernel K_{0,1} is INJECTABLE because `SpecialFunctions.besselk` is Float64-only
# [HOLE: H-HKR-BIGK — arbitrary-precision K_ν is not yet a GrTools primitive].

using SpecialFunctions: besselk
using ..Modes: RelativisticDisc, disc_OmegaK, kappar2, delta_r
import NumericsCore: bisect

# ---- surface-density profile -----------------------------------------------------------------

"""
    SurfaceDensityProfile{T}

Disc surface-density profile Σ(p') entering the analytic Lindblad torque, where `p'` is the
Boyer–Lindquist radius in units of M.  A concrete profile must implement
[`surface_density`](@ref) and [`surface_density_dp`](@ref); [`surface_density_slope`](@ref)
(dlogΣ/dlogp') then follows.
"""
abstract type SurfaceDensityProfile{T<:Real} end

"""
    PowerLawSigma(Σ_ref, p_ref, slope)
    PowerLawSigma(; Σ_ref = 1, p_ref = 10, slope = 0)

Power-law surface density `Σ(p') = Σ_ref (p'/p_ref)^{-slope}` — arXiv:2510.03564 Eq. (22a)
(`Σ_0 (r/10M)^{-Σ_p}`, so `p_ref = 10`, `slope = Σ_p`) and, in the arXiv:2608.17003
normalisation, the height-integrated energy density (12a) `ε_0 (r/a•)^{-𝔭}` (`p_ref = a•/M`,
`slope = 𝔭`, `Σ_ref = ε_0`).  Type-generic; BigFloat-inherited.
"""
struct PowerLawSigma{T<:Real} <: SurfaceDensityProfile{T}
    Σ_ref::T
    p_ref::T
    slope::T
end

function PowerLawSigma(Σ_ref, p_ref, slope)
    T = float(promote_type(typeof(Σ_ref), typeof(p_ref), typeof(slope)))
    PowerLawSigma{T}(T(Σ_ref), T(p_ref), T(slope))
end
PowerLawSigma(; Σ_ref = 1, p_ref = 10, slope = 0) = PowerLawSigma(Σ_ref, p_ref, slope)

"Σ(p') of the profile (Eq. (22a) for [`PowerLawSigma`](@ref))."
surface_density(σ::PowerLawSigma, p) = σ.Σ_ref * (p / σ.p_ref)^(-σ.slope)

"dΣ/dp' of the profile — the gradient term of Eq. (18b)."
surface_density_dp(σ::PowerLawSigma, p) = -σ.slope * surface_density(σ, p) / p

"""
    surface_density_slope(σ, p) -> dlogΣ/dlogp'

Logarithmic surface-density slope `p' Σ'(p')/Σ(p')` (`-slope` for [`PowerLawSigma`](@ref)); the
combination that appears in the Newtonian torque of arXiv:2510.03564 Eq. (26).
"""
surface_density_slope(σ::SurfaceDensityProfile, p) =
    p * surface_density_dp(σ, p) / surface_density(σ, p)

# ---- shear-frame auxiliaries, Eq. (19) ---------------------------------------------------------

# Thin read-only Kerr view for the EXISTING Modes geodesic generics (they depend on the spin only).
_disc_view(χ::T) where {T<:Real} = RelativisticDisc{T}(χ, one(T) / 10, one(T), one(T))

"""
    hkr_epicyclic_poly(χ, p) -> −3χ² + x⁸ − 6x⁴ + 8χx²

The epicyclic auxiliary of arXiv:2510.03564 Eq. (18)/(19) (`x ≡ (p')^{1/4}`).  Identically
`p'² κ²/Ω²` for the Kerr circular geodesic, so it is BUILT from `Modes.kappar2` /
`Modes.disc_OmegaK` rather than re-derived; positive strictly outside the ISCO, zero at it.
"""
function hkr_epicyclic_poly(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    d = _disc_view(T(χ)); pp = T(p)
    pp^2 * kappar2(d, pp) / disc_OmegaK(d, pp)^2
end

"""
    hkr_binding_poly(χ, p) -> 2χ + x⁶ − 3x²

The binding-energy auxiliary of arXiv:2510.03564 Eq. (19) — identically
`1/(p'^{3/2} Ω² (u^t)²)` for the Kerr circular geodesic (`u^t = ` `Spacetime.kerr_circular_dtdtau`);
positive strictly outside the photon ring, and outside the ISCO in particular.
"""
function hkr_binding_poly(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    2 * T(χ) + T(p)^(3//2) - 3 * sqrt(T(p))
end

"""
    hkr_delta_poly(χ, p) -> χ² + x⁸ − 2x⁴

The horizon auxiliary of arXiv:2510.03564 Eq. (19) — identically the Kerr Δ = p'² − 2p' + χ²
(`Modes.delta_r` / `Spacetime.delta_kerr`).
"""
function hkr_delta_poly(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    delta_r(_disc_view(T(χ)), T(p))
end

"""
    hkr_shear_argument(χ, p) -> X

Eq. (19b) of arXiv:2510.03564:

    X = 2 x √(2χ + x⁶ − 3x²) √(−3χ² + x⁸ − 6x⁴ + 8χx²) / [3 (χ² + x⁸ − 2x⁴)],   x = (p')^{1/4}

the argument of the modified Bessel functions.  The paper notes `X = κ/(2|A|)` with `|A|` the local
shear rate in the Fermi frame (Gammie 2004); the Newtonian limit is `X → 2/3`.  In Schwarzschild it
reduces to `2√((p'−6)(p'−3))/(3(p'−2))` (arXiv:2509.20457).
"""
function hkr_shear_argument(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    x = T(p)^(1//4)
    2 * x * sqrt(hkr_binding_poly(χ, p)) * sqrt(hkr_epicyclic_poly(χ, p)) /
        (3 * hkr_delta_poly(χ, p))
end

"""
    hkr_b0(χ, p) -> b₀

Eq. (19c) of arXiv:2510.03564: `b₀ = √(−3χ² + x⁸ − 6x⁴ + 8χx²) / [2x √(2χ + x⁶ − 3x²)]`, the
`K₁` weight of the universal bracket `K₀(X) + b₀K₁(X)`.  Newtonian limit `b₀ → 1/2`;
Schwarzschild `√(p'−6)/(2√(p'−3))`.
"""
function hkr_b0(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    x = T(p)^(1//4)
    sqrt(hkr_epicyclic_poly(χ, p)) / (2 * x * sqrt(hkr_binding_poly(χ, p)))
end

# arXiv:2510.03564 Eq. (C2d) `b_{2,d}` — transcribed VERBATIM from src/appendix.tex:262-271.  The
# same polynomial is the brace of `b_{1,n}` (C2a) up to the prefactor −x¹⁶(2χ+x⁶−3x²)³; the two are
# transcribed independently below and the test file asserts they agree (transcription self-check).
function _b2_denominator(χ::T, x::T) where {T<:Real}
    72χ^7 - 8χ^5 * (51x^4 + 64) * x^4 + 24χ * (288 - 5x^4 * (3 * (x^4 - 8) * x^4 + 68)) * x^12 -
        8χ^3 * (9x^8 + 264x^4 - 220) * x^8 + 3χ^6 * (57x^4 - 49) * x^2 +
        3χ^4 * (59x^8 + 103x^4 + 390) * x^6 +
        χ^2 * (113x^12 - 1485x^8 + 6936x^4 - 6900) * x^10 +
        (-21x^16 + 331x^12 - 1626x^8 + 3348x^4 - 2376) * x^14
end

# arXiv:2510.03564 Eq. (C2a) — the brace of `b_{1,n}`, transcribed VERBATIM from appendix.tex:240-248.
function _b1_numerator_brace(χ::T, x::T) where {T<:Real}
    72χ^7 - 21x^30 + 331x^26 - 360χ * x^24 + (113χ^2 - 1626) * x^22 + 2880χ * x^20 -
        27 * (55χ^2 - 124) * x^18 - 24χ * (3χ^2 + 340) * x^16 +
        3 * (59χ^4 + 2312χ^2 - 792) * x^14 - 192χ * (11χ^2 - 36) * x^12 +
        (309χ^4 - 6900χ^2) * x^10 - 8χ^3 * (51χ^2 - 220) * x^8 +
        9χ^4 * (19χ^2 + 130) * x^6 - 512χ^5 * x^4 - 147χ^6 * x^2
end

# arXiv:2510.03564 Eq. (C2c) `b_{2,n}` — appendix.tex:251-257.
function _b2_numerator(χ::T, x::T) where {T<:Real}
    2 * (63χ^7 + χ^5 * (63x^4 - 236) * x^4 - χ^3 * (15x^8 + 96x^4 + 244) * x^8 +
         45χ^6 * (x^4 - 4) * x^2 + 3χ * (-69x^12 + 548x^8 - 1380x^4 + 1056) * x^12 +
         15χ^4 * (6x^8 - 38x^4 + 75) * x^6 +
         χ^2 * (65x^12 - 864x^8 + 2946x^4 - 2520) * x^10 +
         (-12x^16 + 190x^12 - 927x^8 + 1800x^4 - 1188) * x^14)
end

"""
    hkr_b1(χ, p) -> b₁

Eq. (C1a) of arXiv:2510.03564, `b₁ = b_{1,n}/b_{1,d}` with `b_{1,n}` (C2a) and
`b_{1,d} = (χ+x⁶)³ (−3χ²+x⁸−6x⁴+8χx²)² (χ²+x⁸−2x⁴)⁴` (C2b) — the overall amplitude of the
`Σ(p')` (non-gradient) term of Eq. (18b) and hence the scale of the torque-reversal criterion
[`hkr_torque_reversal_criterion`](@ref).  Newtonian limit `b₁ → 21/√p'`.
"""
function hkr_b1(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    χT = T(χ); x = T(p)^(1//4)
    num = -x^16 * _b1_numerator_brace(χT, x) * hkr_binding_poly(χ, p)^3
    den = (χT + x^6)^3 * hkr_epicyclic_poly(χ, p)^2 * hkr_delta_poly(χ, p)^4
    num / den
end

"""
    hkr_b2(χ, p) -> b₂

Eq. (C1b) of arXiv:2510.03564,
`b₂ = (b_{2,n}/b_{2,d}) [x√(2χ+x⁶−3x²)/√(−3χ²+x⁸−6x⁴+8χx²)]^{-1}` with (C2c)/(C2d) — the `K₁`
weight of the SECOND bracket `K₀(X) + b₂K₁(X)` of Eq. (18b).  Newtonian limit `b₂ → 8/7`.
"""
function hkr_b2(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    χT = T(χ); x = T(p)^(1//4)
    ratio = _b2_numerator(χT, x) / _b2_denominator(χT, x)
    ratio * sqrt(hkr_epicyclic_poly(χ, p)) / (x * sqrt(hkr_binding_poly(χ, p)))
end

"The universal bracket K₀(X) + b₀K₁(X) of Eqs. (18a)/(18b)."
function _bracket0(χ, p, bessel)
    X = hkr_shear_argument(χ, p)
    bessel(0, X) + hkr_b0(χ, p) * bessel(1, X)
end

"The second bracket K₀(X) + b₂K₁(X) of Eq. (18b)."
function _bracket2(χ, p, bessel)
    X = hkr_shear_argument(χ, p)
    bessel(0, X) + hkr_b2(χ, p) * bessel(1, X)
end

# ---- the two Lindblad power coefficients, Eqs. (18a)/(18b) --------------------------------------

"""
    hkr_lindblad_power_quadratic(χ, p, σ; k = 1, q = 1, bessel = besselk) -> ⟨Ė'⟩^{(0)}_{Lr,k}

Eq. (18a) of arXiv:2510.03564 — the `O(j²)` coefficient of the secular power exchanged at the
`k = ±1` Lindblad resonances of a circular equatorial EMRI (mass ratio `q`) in a pressureless
Kerr disc of surface density `σ` (a [`SurfaceDensityProfile`](@ref)), at orbital radius `p'`
(in M).  `k = +1` is the outer, `k = −1` the inner family; the expression is ODD in `k`, so with
equal inner/outer cutoffs the `O(j²)` contributions cancel in the total (Eq. (24)).

M = 1 geometric units.  `bessel(ν, x)` is the modified-Bessel-K kernel (default
`SpecialFunctions.besselk`, Float64 only).
"""
function hkr_lindblad_power_quadratic(χ::Real, p::Real, σ::SurfaceDensityProfile;
                                      k::Integer = 1, q = 1, bessel = besselk)
    _check_outside_isco(χ, p)
    T = float(promote_type(typeof(χ), typeof(p)))
    x = T(p)^(1//4)
    A = hkr_epicyclic_poly(χ, p)
    B = hkr_binding_poly(χ, p)
    C = hkr_delta_poly(χ, p)
    Σ = surface_density(σ, p)
    -16 * q^2 * k * x^26 * Σ * B^4 / (3 * (T(χ) + x^6)^4 * A^(3//2)) *
        _bracket0(χ, p, bessel)^2 / C^2
end

"""
    hkr_lindblad_power_linear(χ, p, σ; q = 1, bessel = besselk) -> ⟨Ė'⟩^{(1)}_{Lr,±}

Eq. (18b) of arXiv:2510.03564 — the `O(j)` coefficient of the secular power exchanged at a
Lindblad resonance (it is EVEN in `k`, i.e. the same for the inner and outer families).  Two
pieces: a gradient term `∝ dΣ/dp'` (negative definite outside the ISCO) and a term `∝ Σ(p')`
whose sign is set by [`hkr_torque_reversal_criterion`](@ref).

This is the single object arXiv:2608.17003 Fig. 9 rests on: with the equal-cutoff assumption the
total analytic Lindblad torque is `m_max²/Ω(p') × ⟨Ė'⟩^{(1)}` (Eq. (24) there, (78) in 2608).
"""
function hkr_lindblad_power_linear(χ::Real, p::Real, σ::SurfaceDensityProfile;
                                   q = 1, bessel = besselk)
    _check_outside_isco(χ, p)
    T = float(promote_type(typeof(χ), typeof(p)))
    x = T(p)^(1//4)
    A = hkr_epicyclic_poly(χ, p)
    B = hkr_binding_poly(χ, p)
    C = hkr_delta_poly(χ, p)
    br0 = _bracket0(χ, p, bessel)
    gradient = -32 * q^2 * x^20 * B^4 * surface_density_dp(σ, p) /
               (9 * (T(χ) + x^6)^3 * A) * br0^2 / C^2
    density = -(8 * hkr_b1(χ, p) * q^2 / 27) * surface_density(σ, p) *
              br0 * _bracket2(χ, p, bessel)
    gradient + density
end

"""
    hkr_lindblad_torque_rate(χ, p, σ; order = :linear, k = 1, q = 1, bessel = besselk)

Eq. (20) of arXiv:2510.03564: for a circular orbit `⟨L̇'_z⟩^{(0,1)}_{Lr,k} = ⟨Ė'⟩^{(0,1)}_{Lr,k}/Ω(p')`
with `Ω` the Kerr circular geodesic frequency (`Modes.disc_OmegaK`).  `order = :linear` returns the
`O(j)` coefficient (Eq. (18b)), `order = :quadratic` the `O(j²)` one (Eq. (18a)).
"""
function hkr_lindblad_torque_rate(χ::Real, p::Real, σ::SurfaceDensityProfile;
                                           order::Symbol = :linear, k::Integer = 1, q = 1,
                                           bessel = besselk)
    T = float(promote_type(typeof(χ), typeof(p)))
    Ω = disc_OmegaK(_disc_view(T(χ)), T(p))
    Ė = if order === :linear
        hkr_lindblad_power_linear(χ, p, σ; q = q, bessel = bessel)
    elseif order === :quadratic
        hkr_lindblad_power_quadratic(χ, p, σ; k = k, q = q, bessel = bessel)
    else
        throw(ArgumentError("hkr_lindblad_torque_rate: order must be :linear or :quadratic, got $(order)"))
    end
    Ė / Ω
end

# ---- the arXiv:2608.17003 per-mode and summed analytic torque, Eqs. (77)-(80) -------------------

"""
    hkr_lindblad_torque_mode(χ, p, m, σ; q = 1, bessel = besselk) -> T^{analytical}_{L,m}

Per-mode analytic Lindblad torque of arXiv:2608.17003 Eq. (77), `T^{analytical}_{L,m} ∝ m`:
summing the inner (`k = −1`) and outer (`k = +1`) families of arXiv:2510.03564 Eq. (17b) at fixed
mode number `m` (their `j`), the `O(m²)` pieces cancel because Eq. (18a) is odd in `k`, leaving

    T^{analytical}_{L,m} = 2 m ⟨Ė'⟩^{(1)}_{Lr,+} / Ω(p')  .

Exactly linear in `m` — the paper's Eq. (77) — with the whole radial/spin dependence in
[`hkr_lindblad_power_linear`](@ref).
"""
function hkr_lindblad_torque_mode(χ::Real, p::Real, m::Real,
                                           σ::SurfaceDensityProfile; q = 1, bessel = besselk)
    2 * m * hkr_lindblad_torque_rate(χ, p, σ; order = :linear, q = q, bessel = bessel)
end

"""
    hkr_lindblad_torque(χ, p, σ; m_max, asymptotic = true, q = 1, bessel = besselk)

Total analytic Lindblad torque of arXiv:2608.17003 Eq. (78) = arXiv:2510.03564 Eq. (24),

    T_L^{analytical} = Σ_{m=1}^{m_max} T^{analytical}_{L,m} ≈ m_max² ⟨Ė'⟩^{(1)}_{Lr,+} / Ω(p')  ,

where the last step uses `Σ_{m=1}^{m_max} m = m_max²/2 + O(m_max)` (2608 main.tex:983, "the
discussion near Eq. (24)").  `asymptotic = true` (the paper's choice, and what Fig. 9 plots) uses
`m_max²`; `asymptotic = false` uses the exact partial sum `m_max(m_max+1)`, i.e. the literal
`Σ_{m=1}^{m_max}` of Eq. (78) with integer `m_max`.

`m_max` is the torque cutoff parameter — see [`hkr_mode_cutoff`](@ref) for the
`≈ 1/h0` prescription of Eq. (79) and the tuned `m_45 ≈ 0.35/h0` of Eq. (80).
"""
function hkr_lindblad_torque(χ::Real, p::Real, σ::SurfaceDensityProfile;
                                      m_max::Real, asymptotic::Bool = true, q = 1,
                                      bessel = besselk)
    m_max ≥ 0 || throw(DomainError(m_max, "hkr_lindblad_torque: need m_max ≥ 0"))
    rate = hkr_lindblad_torque_rate(χ, p, σ; order = :linear, q = q, bessel = bessel)
    weight = asymptotic ? m_max^2 : m_max * (m_max + 1)
    weight * rate
end

"""
    hkr_mode_cutoff(h0; coefficient = 1) -> m_max

Torque cutoff parameter of arXiv:2608.17003: `m_max ≈ coefficient/h0`.  `coefficient = 1` is
Eq. (79), the HegadeKR/Duque prescription `m_max ≈ 1/h0` (the ratio of orbital to sound speed at
the EMRI radius, arXiv:2510.03564 Eq. (23) evaluated where `h = h0`); `coefficient =`
[`HKR_M45_COEFFICIENT`](@ref) is the value `m_45 ∼ 0.35/h0` of Eq. (80), tuned in 2608 so that
the analytic and numerical Lindblad torques match most closely for `(𝔭,𝔮) = (0,0)`.
"""
function hkr_mode_cutoff(h0::Real; coefficient = 1)
    h0 > 0 || throw(DomainError(h0, "hkr_mode_cutoff: need h0 > 0"))
    coefficient / h0
end

"""
    HKR_M45_COEFFICIENT

`0.35`, the tuned cutoff coefficient `m_45 ∼ 0.35/h0` of arXiv:2608.17003 Eq. (80) (main.tex:999):
`m_K` is defined there as the smallest azimuthal mode number at which the NUMERICAL cumulative
torque reaches `K%` of its converged value, and `K = 45` is the choice that best matches the
analytic model for `(𝔭, 𝔮) = (0, 0)`, `h0 = 0.05`.
"""
const HKR_M45_COEFFICIENT = 0.35

"""
    hkr_mode_cutoff_powerlaw(h0, p; Σ_h = 1, p_ref = 10) -> m_max(p')

The radius-dependent cutoff ansatz of arXiv:2510.03564 Eq. (23):
`j_max = (rΩ/c_s)|_{r=p'M} ≈ 1/h|_{r=p'M} = (1/h0)(p'/p_ref)^{Σ_h}` for the aspect-ratio power law
(22b) `h = h0 (r/10M)^{-Σ_h}`.  Reduces to [`hkr_mode_cutoff`](@ref) at `p' = p_ref`
(or `Σ_h = 0`), which is the flat `m_max ≈ 1/h0` used by arXiv:2608.17003 Eq. (79).
"""
function hkr_mode_cutoff_powerlaw(h0::Real, p::Real; Σ_h = 1, p_ref = 10)
    h0 > 0 || throw(DomainError(h0, "hkr_mode_cutoff_powerlaw: need h0 > 0"))
    (p / p_ref)^Σ_h / h0
end

# ---- torque reversal, Eq. (25) -----------------------------------------------------------------

"""
    hkr_torque_reversal_criterion(χ, p; bessel = besselk) -> 𝒢₁(p')

Eq. (25) of arXiv:2510.03564: `𝒢₁ = −b₁[K₀(X) + b₀K₁(X)][K₀(X) + b₂K₁(X)]`, the sign of the
`Σ(p')` term of Eq. (18b) — i.e. the sign of the total pressureless Lindblad torque whenever the
surface-density gradient is not large.  `𝒢₁ < 0` in the Newtonian limit (net negative torque);
in Schwarzschild `𝒢₁ < 0` only for `p' > 8.92`, so the torque REVERSES sign inside that radius.
"""
hkr_torque_reversal_criterion(χ::Real, p::Real; bessel = besselk) =
    -hkr_b1(χ, p) * _bracket0(χ, p, bessel) * _bracket2(χ, p, bessel)

"""
    hkr_torque_reversal_radius(χ; p_max = 40, bessel = besselk, atol = nothing) -> p'_rev

Root of [`hkr_torque_reversal_criterion`](@ref) `𝒢₁(p') = 0` outside the ISCO — the pressureless
torque-reversal location of arXiv:2510.03564 Sec. IV.A (their Fig. 3/4).  `p'_rev ≈ 8.92` in
Schwarzschild (main.tex:662).  Bracketed bisection (`NumericsCore.bisect`) between `r_ISCO(χ)` and
`p_max`; throws if `𝒢₁` does not change sign there.
"""
function hkr_torque_reversal_radius(χ::Real; p_max = 40, bessel = besselk, atol = nothing)
    T = float(promote_type(typeof(χ), typeof(p_max)))
    lo = _isco(T(χ)) * (1 + eps(T)^(1//2))
    hi = T(p_max)
    f(p) = hkr_torque_reversal_criterion(T(χ), p; bessel = bessel)
    bisect(f, lo, hi; atol = atol === nothing ? eps(T) * hi : T(atol))
end

# ---- domain ------------------------------------------------------------------------------------

_isco(χ) = risco_disc(_disc_view(χ))

function _check_outside_isco(χ::Real, p::Real)
    T = float(promote_type(typeof(χ), typeof(p)))
    ri = _isco(T(χ))
    T(p) > ri || throw(DomainError(p,
        "pressureless Lindblad torque: the analytic model of arXiv:2510.03564 needs p' > r_ISCO(χ)/M = $(ri) (κ → 0 at the ISCO makes X → 0 and K₀(X) → ∞)"))
    nothing
end

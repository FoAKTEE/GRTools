# SpiralWaves/background.jl — the stationary, axisymmetric, purely rotating thin-disc equilibrium of
# arXiv:2608.17003 Sec. II.A (main.tex:121-205) on the Page–Thorne equatorial line element (4):
#
#     ds² = −e^{2ν}dt² + e^{2ψ}(dφ − ω̃ dt)² + e^{2μ}dr² + dz²,      (ν, ψ, μ, ω̃)(r) from Spacetime
#
# Power-law equation of state (12a-c), the EMRI orbital frequency (13), the adiabatic index (22), and
# the relativistic hydrostatic balance (11) solved EXACTLY (quadratic in v) for the rotation profile
# v(r), Ω(r) and the four-velocity (8)-(10). M = 1 geometric units, χ = a/M, a• = `aEMRI`; every
# primitive is type-generic `{T}` so BigFloat precision is inherited. Kerr geometry is reached ONLY via
# `Spacetime.kerr_equatorial_potentials(_dr)`, `kerr_isco`; the circular-geodesic frequency and the
# radial epicyclic frequency are the EXISTING `Modes.disc_OmegaK` / `Modes.kappar2` generics (new
# methods on `DiscBackground`, one name = one concept — no second derivation). Additive
# (EGD::sw-background).

using Spacetime: kerr_equatorial_potentials, kerr_equatorial_potentials_dr, kerr_isco
import ..Modes: disc_OmegaK, kappar2, csound, eden0, RelativisticDisc

"""
    DiscBackground{T}(; χ, aEMRI, h0, 𝔭=0, 𝔮=0, b=0.6, εp=1e-5, σε=aEMRI/10, ε0=1)

Background disc + EMRI parameters of arXiv:2608.17003 (M = 1): `χ` = SMBH dimensionless spin
(|χ| < 1; χ < 0 ⇒ the EMRI/disc counter-rotate with the hole, the fluid ALWAYS rotates in +φ),
`aEMRI` = a• the EMRI semi-major axis (> r_ISCO(χ)), `h0` = fiducial aspect ratio, (`𝔭`, `𝔮`) = the
energy-density / sound-speed power-law slopes of (12a,c), `b` = Plummer softening parameter (59)
(nominal 0.6, M:884), `εp` = Landau pole displacement (63) (nominal 1e-5), `σε` = pole-profile width
(nominal 0.1 a•, M:735), `ε0` = ε(a•) the fiducial height-integrated energy density. The mass ratio
q = m•/M is NOT stored: every linear result scales as q² and is quoted in units of (76).

All fields share one type `T` (promoted from the keyword arguments; `DiscBackground{BigFloat}(; …)`
forces it), so BigFloat precision is inherited by every profile.
"""
struct DiscBackground{T<:Real}
    χ::T
    aEMRI::T
    h0::T
    𝔭::T
    𝔮::T
    b::T
    εp::T
    σε::T
    ε0::T
    function DiscBackground{T}(χ, aEMRI, h0, 𝔭, 𝔮, b, εp, σε, ε0) where {T<:Real}
        abs(χ) < 1 || throw(DomainError(χ, "DiscBackground: need |χ| < 1"))
        h0 > 0 || throw(DomainError(h0, "DiscBackground: need h0 > 0"))
        b ≥ 0 || throw(DomainError(b, "DiscBackground: need b ≥ 0"))
        εp ≥ 0 || throw(DomainError(εp, "DiscBackground: need εp ≥ 0"))
        σε > 0 || throw(DomainError(σε, "DiscBackground: need σε > 0"))
        ε0 > 0 || throw(DomainError(ε0, "DiscBackground: need ε0 > 0"))
        aEMRI > kerr_isco(χ) || throw(DomainError(aEMRI, "DiscBackground: need aEMRI > r_ISCO(χ) = $(kerr_isco(χ))"))
        new{T}(χ, aEMRI, h0, 𝔭, 𝔮, b, εp, σε, ε0)
    end
end

function DiscBackground{T}(; χ, aEMRI, h0, 𝔭 = 0, 𝔮 = 0, b = 0.6, εp = 1e-5, σε = nothing,
                           ε0 = 1) where {T<:Real}
    σεv = σε === nothing ? T(aEMRI) / 10 : T(σε)
    DiscBackground{T}(T(χ), T(aEMRI), T(h0), T(𝔭), T(𝔮), T(b), T(εp), σεv, T(ε0))
end

function DiscBackground(; χ, aEMRI, h0, 𝔭 = 0, 𝔮 = 0, b = 0.6, εp = 1e-5, σε = nothing, ε0 = 1)
    T = float(promote_type(typeof(χ), typeof(aEMRI), typeof(h0), typeof(𝔭), typeof(𝔮),
                           typeof(b), typeof(εp), typeof(ε0), σε === nothing ? Bool : typeof(σε)))
    DiscBackground{T}(; χ, aEMRI, h0, 𝔭, 𝔮, b, εp, σε, ε0)
end

# Thin read-only view of the spin for the EXISTING `Modes` Kerr-geodesic generics (disc_OmegaK,
# kappar2 depend on `a` only; h0/r_I/eden_I are irrelevant there and carry no 2601 semantics here).
_modes_view(bg::DiscBackground{T}) where {T} =
    RelativisticDisc{T}(bg.χ, bg.h0, bg.aEMRI, bg.ε0)

"""
    disc_OmegaK(bg::DiscBackground, r) -> Ω_K

Circular equatorial Kerr geodesic frequency Ω_K = 1/(χ + r^{3/2}) (M = 1) — the SAME concept and
formula as `Modes.disc_OmegaK(::RelativisticDisc, r)` (delegated, not re-derived).
"""
disc_OmegaK(bg::DiscBackground, r) = disc_OmegaK(_modes_view(bg), r)

"""
    kappar2(bg::DiscBackground, r) -> κ_r²

Kerr circular-geodesic radial epicyclic frequency squared (BL coordinate time, Torok–Stuchlik closed
form) — delegated to `Modes.kappar2`. See [`kappa2`](@ref).
"""
kappar2(bg::DiscBackground, r) = kappar2(_modes_view(bg), r)

"""
    kappa2(bg, r) -> κ²

Epicyclic frequency squared κ² of the background used by the Newtonian-limit master matrices (50b-d)
and the WKB wavenumber heuristic of the adaptive grid (App. C). IMPLEMENTATION CHOICE
[PRELIMINARY]: the Kerr circular-GEODESIC κ_r² (`Modes.kappar2`, a = 0: (r−6)/r⁴ → 1/r³); the
O(h0²) pressure correction 2rΩΩ' − 2rΩ_KΩ_K' to the epicyclic frequency of the pressure-supported
flow is neglected here (the paper never defines κ beyond "the epicyclic frequency"; the relativistic
master matrices of `master.jl` carry the exact background Ω' and do not consume this).
"""
kappa2(bg::DiscBackground, r) = kappar2(bg, r)

"Kerr equatorial ISCO radius r_ISCO(χ) (Bardeen 1972) — the fiducial inner edge (C1); `Spacetime.kerr_isco`."
r_isco(bg::DiscBackground) = kerr_isco(bg.χ)

"""
    Omega_emri(bg) -> Ω•(a•)

EMRI orbital frequency (13): Ω• = 1/[(a•)^{3/2} + χ] (M = 1) — the circular Kerr geodesic frequency
at r = a• (`disc_OmegaK(bg, aEMRI)`; also the forcing frequency ω = mΩ• of (47)).
"""
Omega_emri(bg::DiscBackground) = disc_OmegaK(bg, bg.aEMRI)

"Height-integrated energy density (12a): ε(r) = ε0 (r/a•)^{−𝔭}."
eps_profile(bg::DiscBackground, r) = bg.ε0 * (r / bg.aEMRI)^(-bg.𝔭)

"""
    eden0(bg::DiscBackground, r)

Background energy density — the same generic as `Modes.eden0`, here the power law (12a)
(`eps_profile`).
"""
eden0(bg::DiscBackground, r) = eps_profile(bg, r)

"""
    csound(bg::DiscBackground, r) -> c_s

Sound speed (12c): c_s(r) = h0 a• Ω•(a•) (r/a•)^{−𝔮/2} — the same generic name as `Modes.csound`
(there the constant-aspect-ratio thermodynamics of arXiv:2601.19123); one name, one concept.
"""
csound(bg::DiscBackground, r) = bg.h0 * bg.aEMRI * Omega_emri(bg) * (r / bg.aEMRI)^(-bg.𝔮 / 2)

"Height-integrated pressure (12b): p = c_s² ε."
pressure(bg::DiscBackground, r) = csound(bg, r)^2 * eps_profile(bg, r)

"""
    gamma1(bg, r) -> Γ1

Adiabatic index parameterization (22): Γ1 = (ε + p) c_s²/p, which with (12b) is exactly 1 + c_s².
"""
gamma1(bg::DiscBackground, r) = (eps_profile(bg, r) + pressure(bg, r)) / pressure(bg, r) * csound(bg, r)^2

"""
    dlogp_dr_over_enthalpy(bg, r) -> p'/(ε + p)

The pressure-gradient source of the hydrostatic balance (11) from the power laws (12a-c):
p ∝ r^{−(𝔭+𝔮)} ⇒ p' = −(𝔭+𝔮) p/r and ε + p = ε(1 + c_s²), so
p'/(ε+p) = −(𝔭+𝔮) c_s² / [r (1 + c_s²)]. Exactly zero for 𝔭 = 𝔮 = 0 (pressureless-gradient disc).
"""
function dlogp_dr_over_enthalpy(bg::DiscBackground, r)
    cs2 = csound(bg, r)^2
    -(bg.𝔭 + bg.𝔮) * cs2 / (r * (1 + cs2))
end

"""
    rotation_profile(bg, r) -> (Ω, v, ut, ut_low, uφ)

Equilibrium rotation of the disc at radius `r` from the relativistic hydrostatic balance (11),

    p'/(ε+p) = −(ν' − v²ψ' + e^{ψ−ν} v ω̃')/(1 − v²),       v ≡ (Ω − ω̃) e^{ψ−ν}   (9),

with (ν, ψ, ω̃)(r) and their r-derivatives from `Spacetime.kerr_equatorial_potentials(_dr)` and
Π ≡ p'/(ε+p) from [`dlogp_dr_over_enthalpy`](@ref). Multiplying (11) by (1−v²) gives the QUADRATIC

    A v² + B v + C = 0,   A = Π + ψ',   B = −e^{ψ−ν} ω̃',   C = −(Π + ν').

BRANCH / UNIQUENESS: on the exterior ψ' > 0 and ν' > 0; for a thin disc |Π| = O(h0²/r) ≪ ψ', ν', so
A > 0 and C < 0, hence AC < 0, the discriminant B² − 4AC > B² ≥ 0 and the two real roots have
OPPOSITE signs. The physical root is the unique POSITIVE one (the fluid rotates in +φ faster than
the frame dragging: Ω > ω̃ for every χ, including χ < 0 where ω̃ < 0), which is continuous with the
Π → 0 geodesic root v_K = [e^{ψ−ν}ω̃' + √(e^{2(ψ−ν)}ω̃'² + 4ψ'ν')]/(2ψ') ⇔ Ω = 1/(χ + r^{3/2})
(`disc_OmegaK`; the negative root is the retrograde geodesic −1/(r^{3/2} − χ)). It is evaluated in
the cancellation-free form (−B + √D)/(2A) for B ≤ 0 and 2C/(−B − √D) for B > 0. A `DomainError` is
thrown if A ≤ 0 or Π + ν' ≤ 0 (pressure gradient too strong for a purely rotating equilibrium) or
v ≥ 1. Then (8)-(10):

    Ω = ω̃ + v e^{ν−ψ},   u^t = e^{−ν}/√(1−v²),   u_t = −e^{ν}(1 + e^{ψ−ν} ω̃ v)/√(1−v²),
    u_φ = e^{ψ} v/√(1−v²)   (u^φ = Ω u^t; u^t u_t + u^φ u_φ = −1).

Returned as a NamedTuple `(Ω, v, ut, ut_low, uφ)` (`ut` = u^t, `ut_low` = u_t, `uφ` = u_φ), all of
the promoted type `T` (BigFloat inherits).
"""
function rotation_profile(bg::DiscBackground{T}, r) where {T}
    rr = convert(promote_type(T, typeof(r)), r)
    pot = kerr_equatorial_potentials(bg.χ, rr)
    d = kerr_equatorial_potentials_dr(bg.χ, rr)
    Π = dlogp_dr_over_enthalpy(bg, rr)
    e_ψν = exp(pot.ψ - pot.ν)
    A = Π + d.dψ
    B = -e_ψν * d.dω̃
    C = -(Π + d.dν)
    (A > 0 && C < 0) || throw(DomainError(rr, "rotation_profile: no purely-rotating equilibrium (A=$A, C=$C) at r=$rr"))
    sq = sqrt(B^2 - 4 * A * C)
    v = B ≤ 0 ? (-B + sq) / (2 * A) : 2 * C / (-B - sq)
    v < 1 || throw(DomainError(rr, "rotation_profile: v = $v ≥ 1 at r = $rr"))
    Ω = pot.ω̃ + v / e_ψν
    γ = 1 / sqrt(1 - v^2)
    ut = exp(-pot.ν) * γ
    ut_low = -exp(pot.ν) * γ * (1 + e_ψν * pot.ω̃ * v)
    uφ = exp(pot.ψ) * v * γ
    (Ω = Ω, v = v, ut = ut, ut_low = ut_low, uφ = uφ)
end

"Rotational angular velocity Ω(r) of the equilibrium disc (`rotation_profile(bg, r).Ω`)."
Omega(bg::DiscBackground, r) = rotation_profile(bg, r).Ω

"""
    omega_d(bg, m, r) -> ω_d

Doppler-shifted forcing frequency (51) in the fluid frame, ω_d ≡ ω − mΩ(r) = m[Ω•(a•) − Ω(r)] with
the EMRI forcing frequency ω = mΩ• (47). Vanishes at the corotation radius r_c (61).
"""
omega_d(bg::DiscBackground, m::Integer, r) = m * (Omega_emri(bg) - Omega(bg, r))

"Plummer softening length (59): ε_soft = a• h0 b."
softening_length(bg::DiscBackground) = bg.aEMRI * bg.h0 * bg.b

"""
    F_J0(bg; q=1) -> F_{J,0}

Characteristic angular-momentum flux (76a): F_{J,0} = ε0 a•⁴ Ω•(a•)² q²/h0³ (Goldreich–Tremaine
units), `q` = m•/M the mass ratio.
"""
F_J0(bg::DiscBackground; q = 1) = bg.ε0 * bg.aEMRI^4 * Omega_emri(bg)^2 * q^2 / bg.h0^3

"Characteristic torque density (76b): F_{J,1} = F_{J,0}/a• = ε0 a•³ Ω•² q²/h0³."
F_J1(bg::DiscBackground; q = 1) = F_J0(bg; q = q) / bg.aEMRI

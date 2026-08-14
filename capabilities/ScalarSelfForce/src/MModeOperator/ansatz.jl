# MModeOperator/ansatz.jl — the m-mode ansatz layer (ESSF node `mode-ansatz-decomposition`,
# W3-C; arXiv:2606.04998 Sec. II.B, main.tex:210-244 + eq:toTilde main.tex:286-289). This file
# is the SINGLE locked home of the mode-decomposition conventions used by every downstream
# module (t-sliced operator, puncture, effective source, extraction):
#
#   Φ(t, r, θ, φ̄) = Σ_m Ψ̃_m(r, z) · sin^{|m|}θ / r · e^{+i m (φ̄ − Ω t)}     (eq:mode_sum)
#   Ψ̃_m = r e^{+i m Ω t} / (2π sin^{|m|}θ) · ∫_{-π}^{π} Φ e^{−i m φ̄} dφ̄     (eq:modeIntegral)
#   Ψ_m  = Ψ̃_m · e^{−i m Ω h(r_*)}                                          (eq:toTilde)
#
# with z = cosθ (θ = acos z ∈ [0, π], sinθ = +√(1−z²) ≥ 0) and the horizon-regular azimuth
# φ̄ = φ + (a/(r₊−r₋))·ln((r−r₊)/(r−r₋)) (main.tex:223-226, Spacetime.kerr_varphi_offset).
# The module NEVER integrates: azimuthal quadrature belongs to the caller; this file only
# applies/strips the locked ansatz factors. See the `mode_term` docstring for the full
# sign/phase convention statement, the conjugate-symmetry derivation, and the spectre
# doxygen sign trap this convention guards against.

"""
    sintheta_of_z(z) -> sinθ

The locked z-chart of the m-mode ansatz (main.tex:237-241): `z = cosθ` with
`θ = acos(z) ∈ [0, π]`, so `sinθ = +√(1−z²) ≥ 0` — this is the ONLY branch used
package-wide (`dz = −sinθ dθ`). Computed as `√((1−z)(1+z))` for accuracy near the
poles. `|z| > 1` raises the natural `DomainError`.
"""
sintheta_of_z(z) = sqrt((1 - z) * (1 + z))

"""
    sinm_theta(m, z) -> sin^{|m|}θ

The polar regularity factor of the ansatz, `sin^{|m|}θ = (1−z²)^{|m|/2}`, factored out so
each mode approaches a nonzero constant towards the poles (main.tex:227, following
Macedo 2024). `m` may be any integer (only `|m|` enters). At the poles `z = ±1` the
factor is `0` for `|m| ≥ 1` and exactly `1` for `m = 0`.
"""
sinm_theta(m::Integer, z) = sintheta_of_z(z)^abs(m)

"""
    sinm_theta_dz(m, z) -> d/dz sin^{|m|}θ

First z-derivative of [`sinm_theta`](@ref), closed form `−|m| z (1−z²)^{(|m|−2)/2}`.
Exact structural limits: `≡ 0` for `m = 0`; `−2z` for `|m| = 2`; diverges at the poles
for `|m| = 1` with the correct signs (`−Inf` at `z = 1`, `+Inf` at `z = −1`).
"""
function sinm_theta_dz(m::Integer, z)
    am = abs(m)
    s = sintheta_of_z(z)
    am == 0 && return zero(s)
    -am * z * s^(am - 2)
end

"""
    sinm_theta_d2z(m, z) -> d²/dz² sin^{|m|}θ

Second z-derivative of [`sinm_theta`](@ref), closed form
`|m| ((|m|−1) z² − 1) (1−z²)^{(|m|−4)/2}`. The `|m| = 2` case is returned as the exact
constant `−2` (the generic form is an indeterminate `0·Inf` at the poles); `m = 0` is
identically `0`; odd `|m| ≤ 3` diverges at the poles (correct behavior of the factor).
"""
function sinm_theta_d2z(m::Integer, z)
    am = abs(m)
    s = sintheta_of_z(z)
    am == 0 && return zero(s)
    am == 2 && return -2 * one(s)
    am * ((am - 1) * z^2 - 1) * s^(am - 4)
end

"""
    kerr_varphi(a, r, phi) -> φ̄

Kerr's original (horizon-regular) azimuthal coordinate `φ̄` from Boyer–Lindquist `φ`
(main.tex:223-226, following Dolan 2011b), composed from the Spacetime closed form:
`φ̄ = φ + kerr_varphi_offset(a, r)`. At `a = 0` the offset is a hard zero (`φ̄ ≡ φ`);
the offset diverges like `−sign(a)·∞` at the horizon. All ansatz phases in this package
are phases in `φ̄`, never in `φ`.
"""
kerr_varphi(a, r, phi) = phi + kerr_varphi_offset(a, r)

"""
    mode_term(m, psim; Omega, ts, r, z, varphi)

One summand of the locked m-mode ansatz (eq:mode_sum, main.tex:228-231):

    mode_term = psim · sin^{|m|}θ / r · e^{+i m (φ̄ − Ω·ts)}

where `varphi` is the horizon-regular azimuth `φ̄` ([`kerr_varphi`](@ref)) and `ts` is the
slicing time matching the amplitude: pass `ts = t` with a t-sliced `Ψ̃_m`, or `ts = s`
(`s = t − h(r_*)`) with a hyperboloidally-sliced `Ψ_m = Ψ̃_m e^{−imΩh}`
([`psi_from_tilde`](@ref)); both give the SAME physical Φ (eq:hyperboloidal_mode_sum,
main.tex:280-289).

LOCKED SIGN/PHASE CONVENTION (the single normative statement for this package):

  * the ansatz phase is `e^{+i m (φ̄ − Ω t)}`; the forward integral (eq:modeIntegral,
    [`mode_amplitude`](@ref)) analyzes with `e^{−i m φ̄}` and restores `e^{+i m Ω t}`;
  * conjugate symmetry — for REAL Φ, `Ψ̃_{−m} = conj(Ψ̃_m)`. Derivation: conjugating
    eq:modeIntegral gives `conj(Ψ̃_m) = r e^{−imΩt}/(2π sin^{|m|}θ) ∫ Φ e^{+imφ̄} dφ̄`
    (Φ real), which is eq:modeIntegral evaluated at `−m` since `|−m| = |m|`; hence
    `Ψ̃_{−m} = conj(Ψ̃_m)`. The hyperboloidal phase preserves it:
    `Ψ_{−m} = Ψ̃_{−m} e^{+imΩh} = conj(Ψ̃_m e^{−imΩh}) = conj(Ψ_m)`;
  * a real field reflection-symmetric about the particle azimuth `φ̄_p` therefore has
    `Ψ_m = (real) · e^{−i m φ̄_p}` — the phase structure measured on the baked
    EffectiveSource fixtures (test/fixtures/effsource headers, W2-D: `Ψ^P_m` pairs are
    `(v cos m φ_p, −v sin m φ_p)`, v real);
  * SPECTRE DOXYGEN SIGN TRAP (audit-spectre-ssf §9): spectre's CircularOrbit.hpp
    doc-strings write `Ψ^P_m ∝ e^{+im(φ−φ_p)}` but `S^eff_m ∝ e^{−im(φ−φ_p)}` — the two
    signs are inconsistent and one is a typo; the spectre CODE applies one common
    rotation to both (validated by its unit test). Trust the fixture-measured
    `e^{−i m φ_p}` structure above, never the doxygen sign.
"""
mode_term(m::Integer, psim; Omega, ts, r, z, varphi) =
    psim * sinm_theta(m, z) / r * cis(m * (varphi - Omega * ts))

"""
    mode_amplitude(m, cm; Omega, ts, r, z)

Forward mode integral (eq:modeIntegral, main.tex:233-236) given the azimuthal Fourier
coefficient `cm = (1/2π) ∫_{-π}^{π} Φ e^{−i m φ̄} dφ̄` at fixed `(ts, r, z)` — the
quadrature producing `cm` belongs to the caller; this strips the locked ansatz factors:

    Ψ_m = cm · r · e^{+i m Ω·ts} / sin^{|m|}θ

`ts` is the slicing time of the samples (`t` → t-sliced `Ψ̃_m`; `s` → hyperboloidal
`Ψ_m`). Valid at `|z| < 1` for `m ≠ 0` (the sin factor vanishes at the poles); sign/phase
conventions and conjugate symmetry: see [`mode_term`](@ref).
"""
mode_amplitude(m::Integer, cm; Omega, ts, r, z) =
    cm * r * cis(m * Omega * ts) / sinm_theta(m, z)

"""
    reconstruct_field(ms, psims; Omega, ts, r, z, varphi)

Truncated mode-sum reconstruction (eq:mode_sum, main.tex:228-231):
`Φ ≈ Σ_k mode_term(ms[k], psims[k]; ...)` over the supplied mode list. For a real field
supply conjugate pairs `(−m, conj(Ψ_m))` explicitly (see the [`mode_term`](@ref)
conjugate-symmetry convention; the ×2-real-part shortcut is the extraction layer's
business, not this one's).
"""
reconstruct_field(ms, psims; Omega, ts, r, z, varphi) =
    sum(mode_term(ms[k], psims[k]; Omega, ts, r, z, varphi) for k in eachindex(ms, psims))

"""
    psi_from_tilde(m, tildepsi, Omega, h)

Hyperboloidal slicing phase (eq:toTilde, main.tex:286-289): `Ψ_m = Ψ̃_m e^{−i m Ω h}`
with `h = h(r_*)` the height function (`s = t − h`). Inverse: [`tilde_from_psi`](@ref).
"""
psi_from_tilde(m::Integer, tildepsi, Omega, h) = tildepsi * cis(-m * Omega * h)

"""
    tilde_from_psi(m, psim, Omega, h)

Inverse of [`psi_from_tilde`](@ref): `Ψ̃_m = Ψ_m e^{+i m Ω h}` (eq:toTilde).
"""
tilde_from_psi(m::Integer, psim, Omega, h) = psim * cis(m * Omega * h)

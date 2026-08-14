# SphericalHarmonics — the γ=0 spin-weighted spherical-Y front door + the UNIFIED m<0 parity
# surface (design §4.2). A thin convention layer over NumericsCore's already-promoted angular
# seed (`spin_weighted_Y`/`Ylm_real`/`Ylm_eq`): this file carries NO Legendre/half-angle copy.
#
# It collapses the per-package negative-m forks — FEW `ylm`, KerrBackwards `ylm_few`/`ylms_few`,
# DiscModes real-Y, Punctures `(-1)^m·conj` — into ONE parity-flagged primitive
# `spin_weighted_spherical_Y(s,l,m,θ,φ; convention)` + the sign helper `angular_parity`. The
# legacy names are kept as re-export shims so the oracle packages stay green.

"""
    angular_parity(l, m; convention=:none) -> Int

Multiplicative parity SIGN applied to the bare seed `{}_sY_{lm}` when forming the negative-`m`
harmonic under a packaging `convention`. `m ≥ 0` is always a no-op (`+1`).

- `:none`   — no flip (the BHPT / bare-seed convention).
- `:few`    — `(-1)^l` for `m<0` (FEW `GetYlms` / KerrBackwards `ylm_few` parity operator).
- `:reality`— `(-1)^m` for `m<0`, the `m`-part of the general spin-weighted reality relation
  `{}_sY_{l,-|m|} = (-1)^{s+|m|} conj({}_{-s}Y_{l,|m|})`. This helper supplies ONLY the `m`-part
  `(-1)^m`; the `s`-part `(-1)^s` and the conjugation/`s→−s` flip are applied by the front door
  [`spin_weighted_spherical_Y`](@ref), so the full spin-weighted phase is `(-1)^{s+m}`.

`(-1)^l` genuinely needs `l`, so the helper takes `(l, m)` (the design's loose `angular_parity(m)`
name realized with the required degree).
"""
function angular_parity(l::Integer, m::Integer; convention::Symbol=:none)
    m >= 0 && return 1
    convention === :none && return 1
    convention === :few && return isodd(l) ? -1 : 1
    convention === :reality && return isodd(m) ? -1 : 1
    throw(ArgumentError("unknown parity convention $(convention); use :none, :few or :reality"))
end

"""
    spin_weighted_spherical_Y(s, l, m, θ, φ; convention=:none) -> Complex

Unified γ=0 spin-weighted spherical-harmonic front door — `{}_sY_{lm}(θ,φ)` from the canonical
NumericsCore seed, with the negative-`m` packaging `convention` folded in (the single tested
parity flag replacing the FEW/KerrBackwards/DiscModes/Punctures forks).

- `:none`   — bare seed `{}_sY_{lm}` (default; bit-identical to the BHPT
  `SpinWeightedSphericalHarmonicY`).
- `:few`    — `m<0 ⇒ (-1)^l {}_sY_{lm}` (FEW/KerrBackwards parity).
- `:reality`— GENERAL spin-weighted negative-`m` reality relation (scalar / DiscModes /
  Punctures use): `m<0` is built as `{}_sY_{l,-|m|} = (-1)^{s+|m|} conj({}_{-s}Y_{l,|m|})` — the
  `s→−s`, `+|m|` seed, conjugated, with the full `(-1)^{s+m}` phase. Value-equivalent to (indeed
  bit-identical with) the bare seed. For `s=0` it reduces to the scalar
  `{}_0Y_{l,-|m|} = (-1)^{|m|} conj({}_0Y_{l,|m|})`; `m≥0` is the bare seed.

Type-generic: BigFloat colatitude/precision is inherited from the core seed.
"""
function spin_weighted_spherical_Y(s::Integer, l::Integer, m::Integer, θ::Real, φ::Real;
                                   convention::Symbol=:none)
    if convention === :reality
        m >= 0 && return spin_weighted_Y(s, l, m, θ, φ)
        am = abs(m)
        # {}_sY_{l,-|m|} = (-1)^{s+m} conj({}_{-s}Y_{l,|m|}); angular_parity supplies (-1)^m,
        # the (-1)^s front-door factor completes the phase (s=0 ⇒ the scalar (-1)^m conj case).
        s_phase = isodd(s) ? -1 : 1
        return s_phase * angular_parity(l, m; convention=:reality) *
               conj(spin_weighted_Y(-s, l, am, θ, φ))
    end
    return angular_parity(l, m; convention=convention) * spin_weighted_Y(s, l, m, θ, φ)
end

"""
    SpinWeightedSphericalHarmonicY(s, l, m, θ, φ; parity=:none) -> Complex

BHPT-named spin-weighted spherical harmonic `{}_sY_{lm}(θ,φ)` (the γ=0 seed). Default
`parity=:none` is bit-identical to the legacy `SpinWeightedSpheroidalHarmonics`
implementation; `parity=:few`/`:reality` route through the unified
[`spin_weighted_spherical_Y`](@ref) front door. The high-`l` BigFloat path is inherited from the
NumericsCore seed.
"""
SpinWeightedSphericalHarmonicY(s::Integer, l::Integer, m::Integer, θ::Real, φ::Real;
                               parity::Symbol=:none) =
    spin_weighted_spherical_Y(s, l, m, θ, φ; convention=parity)

# --- legacy convention shims (kept so the FEW/KerrBackwards oracle tests stay green) ---

"""
    ylm(s, l, m, θ, φ) -> Complex
    ylm(l, m, θ, φ)    -> Complex   # s = -2 (gravitational, FEW default)

FEW `few.utils.ylm` harmonic — the bare seed (`convention=:none`). The 2-argument form fixes
`s=-2`, the gravitational spin weight used throughout the FEW waveform pipeline.
"""
ylm(s::Integer, l::Integer, m::Integer, θ::Real, φ::Real) =
    spin_weighted_spherical_Y(s, l, m, θ, φ; convention=:none)
ylm(l::Integer, m::Integer, θ::Real, φ::Real) = ylm(-2, l, m, θ, φ)

"""
    ylm_few(l, m, θ, φ) -> Complex

KerrBackwards `(-2)` spin-weighted harmonic with the FEW parity convention
(`convention=:few`): `m<0 ⇒ (-1)^l {}_{-2}Y_{lm}`.
"""
ylm_few(l::Integer, m::Integer, θ::Real, φ::Real) =
    spin_weighted_spherical_Y(-2, l, m, θ, φ; convention=:few)

"""
    ylms_few(l_arr, m_arr, θ, φ) -> Vector{Complex}

Vectorised [`ylm_few`](@ref) over matching `(l,m)` arrays at one viewing angle (mirrors
`few.utils.ylm.GetYlms.__call__`).
"""
function ylms_few(l_arr::AbstractVector{<:Integer}, m_arr::AbstractVector{<:Integer},
                  θ::Real, φ::Real)
    length(l_arr) == length(m_arr) ||
        throw(ArgumentError("l_arr and m_arr must have equal length"))
    return [ylm_few(l_arr[i], m_arr[i], θ, φ) for i in eachindex(l_arr)]
end

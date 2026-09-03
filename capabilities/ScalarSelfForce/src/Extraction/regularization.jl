# Extraction/regularization.jl — Heffernan m-mode regularization parameters of the RADIAL
# scalar self-force for circular equatorial Kerr geodesics (ESSF node
# `heffernan-regularization`, GLOBAL_DAG row 45; DESIGN §2.2 `Extraction/regularization.jl`).
#
# TRANSCRIBED FROM (read-only mirror `ref-paper/arxiv-1211.6446/src/KerrRegularization.tex`
# = A. Heffernan, A. Ottewill and B. Wardell, "High-order expansions of the Detweiler-Whiting
# singular field in Kerr spacetime", Phys. Rev. D 89, 024030 (2014), arXiv:1211.6446 —
# hereafter [H14]; the CRP arXiv:2606.04998 cites it as \ccite{Heffernan_2014} at
# main.tex:522 and main.tex:723):
#   :1675-1680  the explicit m-dependence  F^m_{a[2j]}  (orders 2, 4, 6, 8)  -> mmode_reg_factor
#   :1689-1699  F_{r[2]}, circular orbit                                     -> reg_param_Fr(2, …)
#   :1700-1731  F_{r[4]}, circular orbit                                     -> reg_param_Fr(4, …)
#   :1610-1671  the derivation (folded m ≥ 0 modes at :1583-1589)
#   :1784-1798  the a = 0.6M, r0 = 10M worked example (m^-6 residual after subtraction)
# The analytic tail completion (`mmode_reg_factor_tail`) is DERIVED here from the same printed
# m-dependence by telescoping; it agrees with the independently published closed form of
# Thornburg & Wardell (arXiv:1610.09319, `ksf.tex:2225-2240`,
# S^{(Φ)}_{α,m+1} = 1/((α−1)(m − (α−3)/2)_{α−1})) and is gated exactly in
# test/test_regularization.jl.
#
# NAME COLLISION WARNING (a real trap in [H14]): the symbols `F_{r[2]}` / `F_{r[4]}` are used
# TWICE in that paper for DIFFERENT objects — the ℓ-MODE parameters of §IV
# (:756-798, built from complete elliptic integrals 𝓔, 𝓚) and the m-MODE parameters of §V
# (:1689-1731, purely algebraic in (a, r0, M)). THIS FILE OWNS ONLY THE §V m-MODE ONES; no
# elliptic integral is needed anywhere below.
#
# SCOPE: parameters + tail ONLY. Extracting F_r^m from a solved mode is node
# `mode-self-force-extraction` (row 44) and performing Σ_m is `mode-sum-total-force` (row 46).

const _HEFF_TEX = "ref-paper/arxiv-1211.6446/src/KerrRegularization.tex"

"""
    HEFFERNAN_MMODE_ORDERS

The regularization-parameter orders `2j` whose explicit m-dependence [H14] prints
(`KerrRegularization.tex:1675-1680`): `(2, 4, 6, 8)`. [`mmode_reg_factor`](@ref) and
[`mmode_reg_factor_tail`](@ref) accept exactly these.
"""
const HEFFERNAN_MMODE_ORDERS = (2, 4, 6, 8)

"""
    HEFFERNAN_CIRCULAR_FR_ORDERS

The orders for which [H14] prints the CIRCULAR-orbit Kerr coefficient `F_{r[2j]}` in closed
form (`KerrRegularization.tex:1689-1731`): `(2, 4)`. [`reg_param_Fr`](@ref) accepts exactly
these; the printed higher orders `F_{r[6]}`, `F_{r[8]}` are numerical fits in [H14]
(`:1795-1797`), not closed forms, and are deliberately not shipped.
"""
const HEFFERNAN_CIRCULAR_FR_ORDERS = (2, 4)

# ---------------------------------------------------------------------- m-dependence
"""
    mmode_reg_factor_over_pi(order, m) -> Rational{BigInt}

`π · F^m_{r[order]} / F_{r[order]}` — the EXACT rational m-dependence of the m-mode
regularization parameters, transcribed verbatim from [H14]
(`$(_HEFF_TEX):1675-1680`):

    F^m_{a[2]} = −4   F_{a[2]} / (π (2m−1)(2m+1))
    F^m_{a[4]} = +24  F_{a[4]} / (π (2m−3)(2m−1)(2m+1)(2m+3))
    F^m_{a[6]} = −480 F_{a[6]} / (π (2m−5)(2m−3)(2m−1)(2m+1)(2m+3)(2m+5))
    F^m_{a[8]} = 20160 F_{a[8]} / (π (2m−7)(2m−5)(2m−3)(2m−1)(2m+1)(2m+3)(2m+5)(2m+7))

`order ∈ `[`HEFFERNAN_MMODE_ORDERS`](@ref) and `m ≥ 0`; the denominator never vanishes at
integer `m`. Returning the π-stripped rational keeps mode-sum bookkeeping EXACT — divide by
`π` once at the end, or use [`mmode_reg_factor`](@ref).

LOCKED CONVENTION (the m-mode counterpart of the normative mode-decomposition statement in
[`mode_term`](@ref); [H14]:1583-1589 and CRP main.tex:228-231): these are FOLDED `m ≥ 0`
parameters — the ×2 real-part factor of `Φ = Φ^{(0)} + 2 Σ_{m≥1} Re[Φ^{(m)} e^{imφ_0}]` is
already inside them, so a mode sum adds `F^m_{r[2j]}` once per `m ≥ 1`. Retrograde orbits are
encoded by `a < 0` (motion always `+φ`), never by `m < 0`.

`m = 0` SPECIAL CASE (operational rule): the closed form is a single analytic function of
`m` with no `m = 0` branch, so the ×2 real-part folding is INSIDE the parameter at every `m`,
`m = 0` included. Node 46 (`mode-sum-total-force`) must therefore use
`mmode_reg_Fr(order, orbit, 0)/2` wherever its `m = 0` force term is the UNFOLDED
`q ∇(Ψ^ℛ_0 / r)` of Osburn–Nishimura Eq. (3.10) (first branch); for `m ≥ 1` use the value as
returned. Equivalently the doubly-infinite sum vanishes,
`Σ_{m∈ℤ} F^{|m|}/2 = F^0/2 + Σ_{m≥1} F^m = 0`, i.e.
`mmode_reg_factor_tail_over_pi(order, 1) == −mmode_reg_factor_over_pi(order, 0)//2` (gated).
The subtraction is in any case a LARGE-m model: [H14] uses it for the tail, never at small `m`.
"""
function mmode_reg_factor_over_pi(order::Integer, m::Integer)
    m >= 0 || throw(ArgumentError(
        "m-mode regularization parameters are folded (m ≥ 0; retrograde is a < 0), got m = $m"))
    mm = big(m)
    order == 2 && return Rational{BigInt}(-4, (2mm - 1) * (2mm + 1))
    order == 4 && return Rational{BigInt}(24, (2mm - 3) * (2mm - 1) * (2mm + 1) * (2mm + 3))
    order == 6 && return Rational{BigInt}(-480,
        (2mm - 5) * (2mm - 3) * (2mm - 1) * (2mm + 1) * (2mm + 3) * (2mm + 5))
    order == 8 && return Rational{BigInt}(20160,
        (2mm - 7) * (2mm - 5) * (2mm - 3) * (2mm - 1) *
        (2mm + 1) * (2mm + 3) * (2mm + 5) * (2mm + 7))
    throw(ArgumentError("order must be one of $HEFFERNAN_MMODE_ORDERS ($(_HEFF_TEX):1675-1680), got $order"))
end

"""
    mmode_reg_factor([T=Float64,] order, m) -> T

`F^m_{r[order]} / F_{r[order]}` in floating-point type `T` — [`mmode_reg_factor_over_pi`](@ref)
divided by `π`, so `T = BigFloat` inherits the ambient `precision`. Multiply by
[`reg_param_Fr`](@ref) (or call [`mmode_reg_Fr`](@ref)) for the parameter itself.
Asymptotically `F^m_{r[2]} ∼ −m^{-2}/π` and `F^m_{r[4]} ∼ 3 m^{-4}/(2π)` times the
coefficient — the `m^{-4}` tail the CRP subtracts (main.tex:520-522, 723).
"""
mmode_reg_factor(::Type{T}, order::Integer, m::Integer) where {T<:AbstractFloat} =
    T(mmode_reg_factor_over_pi(order, m)) / T(π)
mmode_reg_factor(order::Integer, m::Integer) = mmode_reg_factor(Float64, order, m)

# ------------------------------------------------------------------- analytic tail sum
"""
    mmode_reg_factor_tail_over_pi(order, mstart) -> Rational{BigInt}

`π · Σ_{m=mstart}^∞ F^m_{r[order]} / F_{r[order]}` — the EXACT closed form of the analytic
large-m tail, the completion `mode-sum-total-force` (row 46) adds to its truncated numerical
sum `Σ_{m=0}^{m_max}` with `mstart = m_max + 1`.

DERIVATION (from the printed m-dependence alone, [H14]:1675-1680): writing
`D_j(m) = Π_{k=1}^{j} (2m−(2k−1))(2m+(2k−1))` and
`h_j(m) = 1 / [Π_{k=1}^{j}(2m−(2k−1)) · Π_{k=1}^{j−1}(2m+(2k−1))]`, one has
`h_j(m) − h_j(m+1) = (4j−2)/D_j(m)`, so the sum telescopes to
`Σ_{m≥mstart} 1/D_j(m) = h_j(mstart)/(4j−2)`. Hence

    order 2: −2   / (2M−1)
    order 4: +4   / ((2M−3)(2M−1)(2M+1))
    order 6: −48  / ((2M−5)(2M−3)(2M−1)(2M+1)(2M+3))
    order 8: +1440/ ((2M−7)(2M−5)(2M−3)(2M−1)(2M+1)(2M+3)(2M+5))

with `M = mstart`. Equivalent to the published `S^{(Φ)}_{α,m+1} = 1/((α−1)(m−(α−3)/2)_{α−1})`
of Thornburg & Wardell (arXiv:1610.09319, `ksf.tex:2225-2240`), and gated against
brute-force partial sums in exact rational arithmetic. Valid for any `mstart ≥ 0`; see the
`m = 0` folding note in [`mmode_reg_factor_over_pi`](@ref).
"""
function mmode_reg_factor_tail_over_pi(order::Integer, mstart::Integer)
    mstart >= 0 || throw(ArgumentError("mstart must be ≥ 0 (folded m ≥ 0 modes), got $mstart"))
    mm = big(mstart)
    order == 2 && return Rational{BigInt}(-2, 2mm - 1)
    order == 4 && return Rational{BigInt}(4, (2mm - 3) * (2mm - 1) * (2mm + 1))
    order == 6 && return Rational{BigInt}(-48,
        (2mm - 5) * (2mm - 3) * (2mm - 1) * (2mm + 1) * (2mm + 3))
    order == 8 && return Rational{BigInt}(1440,
        (2mm - 7) * (2mm - 5) * (2mm - 3) * (2mm - 1) * (2mm + 1) * (2mm + 3) * (2mm + 5))
    throw(ArgumentError("order must be one of $HEFFERNAN_MMODE_ORDERS ($(_HEFF_TEX):1675-1680), got $order"))
end

"""
    mmode_reg_factor_tail([T=Float64,] order, mstart) -> T

`Σ_{m=mstart}^∞ F^m_{r[order]} / F_{r[order]}` in floating-point type `T` —
[`mmode_reg_factor_tail_over_pi`](@ref) divided by `π`. See
[`mmode_reg_tail_Fr`](@ref) for the orbit-level tail.
"""
mmode_reg_factor_tail(::Type{T}, order::Integer, mstart::Integer) where {T<:AbstractFloat} =
    T(mmode_reg_factor_tail_over_pi(order, mstart)) / T(π)
mmode_reg_factor_tail(order::Integer, mstart::Integer) =
    mmode_reg_factor_tail(Float64, order, mstart)

# ------------------------------------------------------- circular-Kerr coefficients
"""
    ScalarSelfForce.reg_param_Fr_bracket(order, a, r0, M, s)

The square-bracket polynomial of the printed circular-orbit coefficients — `[H14]:1691-1698`
for `order == 2` and `[H14]:1702-1730` for `order == 4` — with `s = √(M r0)` supplied by the
caller. Written in the paper's FACTORED display form and containing only `+ - * ^`, so it is
EXACT on `Rational{BigInt}` inputs (that is how the transcription re-diff in
`test/test_regularization.jl` compares it against an independently expanded monomial table).
Homogeneous of degree 9 / 17 in `(a, M, r0, s)`. Internal: the assembled parameter is
[`reg_param_Fr`](@ref).
"""
function reg_param_Fr_bracket(order::Integer, a, r0, M, s)
    if order == 2
        return (24 * a^7 * M^2
                - 24 * a^6 * M * s * (M - 2 * r0)
                - 4 * a^5 * M * r0 * (23 * M^2 + M * r0 - 6 * r0^2)
                + 2 * a^4 * M * r0 * s * (45 * M^2 - 112 * M * r0 + 31 * r0^2)
                + 2 * a^3 * M * r0^2 * (45 * M^3 + 45 * M^2 * r0 - 73 * M * r0^2 + 19 * r0^3)
                - 3 * a^2 * r0^2 * s * (29 * M^4 - 88 * M^3 * r0 + 38 * M^2 * r0^2
                                        - 4 * M * r0^3 + r0^4)
                - 6 * a * M * r0^4 * (29 * M^3 - 43 * M^2 * r0 + 21 * M * r0^2 - 3 * r0^3)
                - 3 * r0^5 * s * (29 * M^3 - 25 * M^2 * r0 + 3 * M * r0^2 + r0^3))
    elseif order == 4
        return (-23040 * a^14 * M^2 * s
                + 11520 * a^13 * M^2 * r0 * (M - 8 * r0)
                + 384 * a^12 * M * r0 * s * (461 * M^2 - 81 * M * r0 - 360 * r0^2)
                - 192 * a^11 * M * r0^2 * (307 * M^3 - 3780 * M^2 * r0 + 1233 * M * r0^2
                                           + 480 * r0^3)
                - 64 * a^10 * r0^2 * s * (8549 * M^4 - 3593 * M^3 * r0 - 16212 * M^2 * r0^2
                                          + 6336 * M * r0^3 + 360 * r0^4)
                + 32 * a^9 * M * r0^3 * (2835 * M^4 - 69401 * M^3 * r0 + 46565 * M^2 * r0^2
                                         + 13779 * M * r0^3 - 8748 * r0^4)
                + 192 * a^8 * r0^3 * s * (4470 * M^5 - 3621 * M^4 * r0 - 15645 * M^3 * r0^2
                                          + 12662 * M^2 * r0^3 - 1529 * M * r0^4 - 342 * r0^5)
                + 16 * a^7 * r0^4 * (-1479 * M^6 + 210966 * M^5 * r0 - 224760 * M^4 * r0^2
                                     - 49213 * M^3 * r0^3 + 93619 * M^2 * r0^4
                                     - 19953 * M * r0^5 + 180 * r0^6)
                - 16 * a^6 * r0^4 * s * (43101 * M^6 - 61443 * M^5 * r0 - 271980 * M^4 * r0^2
                                         + 343776 * M^3 * r0^3 - 100489 * M^2 * r0^4
                                         - 7763 * M * r0^5 + 4158 * r0^6)
                + 12 * a^5 * r0^5 * (-2367 * M^7 - 221220 * M^6 * r0 + 337457 * M^5 * r0^2
                                     + 71894 * M^4 * r0^3 - 262111 * M^3 * r0^4
                                     + 111498 * M^2 * r0^5 - 14459 * M * r0^6 + 588 * r0^7)
                + 3 * a^4 * r0^5 * s * (76125 * M^7 - 176307 * M^6 * r0 - 1157559 * M^5 * r0^2
                                        + 1949709 * M^4 * r0^3 - 855873 * M^3 * r0^4
                                        - 26505 * M^2 * r0^5 + 76235 * M * r0^6 - 8065 * r0^7)
                + 12 * a^3 * r0^7 * (76125 * M^7 - 152637 * M^6 * r0 - 93174 * M^5 * r0^2
                                     + 281414 * M^4 * r0^3 - 166063 * M^3 * r0^4
                                     + 36555 * M^2 * r0^5 - 3480 * M * r0^6 + 460 * r0^7)
                + 18 * a^2 * r0^8 * s * (76125 * M^6 - 145182 * M^5 * r0 + 70771 * M^4 * r0^2
                                         + 16696 * M^3 * r0^3 - 19905 * M^2 * r0^4
                                         + 3190 * M * r0^5 + 225 * r0^6)
                + 36 * a * r0^10 * (25375 * M^6 - 47369 * M^5 * r0 + 31856 * M^4 * r0^2
                                    - 8692 * M^3 * r0^3 + 705 * M^2 * r0^4 - 75 * M * r0^5
                                    + 40 * r0^6)
                + 9 * r0^11 * s * (25375 * M^5 - 47015 * M^4 * r0 + 29014 * M^3 * r0^2
                                   - 4814 * M^2 * r0^3 - 1365 * M * r0^4 + 405 * r0^5))
    end
    throw(ArgumentError("order must be one of $HEFFERNAN_CIRCULAR_FR_ORDERS " *
                        "($(_HEFF_TEX):1689-1731), got $order"))
end

"""
    reg_param_Fr(order, a, r0; M = 1) -> T
    reg_param_Fr(order, orbit::CircularOrbit) -> T

The m-INDEPENDENT circular-equatorial-Kerr coefficient `F_{r[order]}` of the radial
scalar-self-force m-mode regularization parameters, transcribed from [H14]
(`$(_HEFF_TEX):1689-1699` for `order == 2`, `:1700-1731` for `order == 4`):

    F_{r[2]} = M · B₂(a, r0, M) /
               (24 r0⁴ [a M + r0 √(M r0)] [2a√(M r0) + r0(r0−3M)]^{3/2} [a² + r0(r0−2M)]^{3/2})
    F_{r[4]} = M² · B₄(a, r0, M) /
               (1440 r0⁹ [a M + r0 √(M r0)] [2a√(M r0) + r0(r0−3M)]^{7/2} [a² + r0(r0−2M)]^{3/2})

with `B₂`, `B₄` the printed square brackets ([`ScalarSelfForce.reg_param_Fr_bracket`](@ref)).
Everything here is ALGEBRAIC in `(a, r0, M)` — the m-mode parameters need no elliptic
integrals (unlike the ℓ-mode ones of [H14] §IV, see this file's name-collision warning).

The full parameter is `F^m_{r[order]} = F_{r[order]} · `[`mmode_reg_factor`](@ref)`(order, m)`
= [`mmode_reg_Fr`](@ref). Signs and normalization follow [H14] exactly, in the folded `m ≥ 0`
×2 real-part convention locked in [`mmode_reg_factor_over_pi`](@ref) (and in
[`mode_term`](@ref) for the field itself). Retrograde orbits are `a < 0`; the coefficients are
smooth through `a = 0` and are NOT symmetric under `a → −a`.

Units: `q²/M²`, with the scalar charge `q = 1` package-wide (the same normalization as the
Warburton–Barack reference tables, which print `(M²/q²)F_r`).

VALIDATION ASYMMETRY (read before trusting either order): `order == 4` is cross-checked
against the INDEPENDENT Wardell–Thornburg notebook oracle (`Fmr[4]`, arXiv:1610.09319
`anc/m-ModeKerrRegularizationParameters.nb`, baked as `test/fixtures/heffernan_mmode_reg.jl`).
`order == 2` has NO independent source — the notebook ships no `[2]` cell and no published
Schwarzschild m-mode statement exists — so it rests on the double transcription (factored
display vs an independently hand-expanded monomial table, exact in `Rational{BigInt}`) plus
the `a → 0` reduction below [HOLE]. The CRP uses only `[4]` (main.tex:723), its fourth-order
puncture having already removed the `[2]` parameter through the effective source
([H14]:1786-1789).

`T = float(promote_type(typeof(a), typeof(r0), typeof(M)))`, so `BigFloat` precision is
inherited. `order ∈ `[`HEFFERNAN_CIRCULAR_FR_ORDERS`](@ref) (otherwise `ArgumentError`);
radii at or inside the circular photon orbit — where the printed factor
`[2a√(M r0) + r0(r0−3M)]` is non-positive and no circular geodesic exists — raise a
`DomainError`.

SCHWARZSCHILD LIMIT (`a → 0`), obtained by cancelling the surviving bracket group against the
prefactor and gated in `test/test_regularization.jl`:

    F_{r[2]}|_{a=0} = −M (r0³ + 3M r0² − 25M² r0 + 29M³) / (8 r0³ (r0−3M)^{3/2} (r0−2M)^{3/2})
    F_{r[4]}|_{a=0} =  M² (405 r0⁵ − 1365M r0⁴ − 4814M² r0³ + 29014M³ r0²
                           − 47015M⁴ r0 + 25375M⁵) / (160 r0⁴ (r0−3M)^{7/2} (r0−2M)^{3/2})
"""
function reg_param_Fr(order::Integer, a::Real, r0::Real;
                      M::Real = one(float(promote_type(typeof(a), typeof(r0)))))
    order in HEFFERNAN_CIRCULAR_FR_ORDERS || throw(ArgumentError(
        "circular-orbit F_r coefficients are printed only for orders " *
        "$HEFFERNAN_CIRCULAR_FR_ORDERS ($(_HEFF_TEX):1689-1731), got $order"))
    T = float(promote_type(typeof(a), typeof(r0), typeof(M)))
    aT, rT, MT = T(a), T(r0), T(M)
    (rT > 0 && MT > 0) || throw(DomainError((a, r0, M), "need r0 > 0 and M > 0"))
    s = sqrt(MT * rT)
    u = 2 * aT * s + rT * (rT - 3 * MT)      # [2a√(Mr0) + r0(r0−3M)] — circular-geodesic factor
    d = aT^2 + rT * (rT - 2 * MT)            # [a² + r0(r0−2M)] = Δ(r0)
    u > 0 || throw(DomainError(r0,
        "[2a√(Mr0) + r0(r0−3M)] ≤ 0: no circular geodesic at or inside the photon orbit " *
        "(a = $a, M = $M)"))
    d > 0 || throw(DomainError(r0, "Δ(r0) = a² + r0(r0−2M) ≤ 0: r0 is at or inside r₊ (a = $a)"))
    b = reg_param_Fr_bracket(order, aT, rT, MT, s)
    return order == 2 ?
        MT * b / (24 * rT^4 * (aT * MT + rT * s) * (u * sqrt(u)) * (d * sqrt(d))) :
        MT^2 * b / (1440 * rT^9 * (aT * MT + rT * s) * (u^3 * sqrt(u)) * (d * sqrt(d)))
end

reg_param_Fr(order::Integer, orbit::CircularOrbit) =
    reg_param_Fr(order, orbit.a, orbit.r0; M = orbit.M)

# --------------------------------------------------------- assembled m-mode parameters
"""
    mmode_reg_Fr(order, a, r0, m; M = 1) -> T
    mmode_reg_Fr(order, orbit::CircularOrbit, m) -> T

The m-mode regularization parameter itself,
`F^m_{r[order]} = `[`reg_param_Fr`](@ref)`(order, …) · `[`mmode_reg_factor`](@ref)`(order, m)`
([H14] `$(_HEFF_TEX):1675-1731`), in the folded `m ≥ 0` ×2 real-part convention.

This is the PER-MODE SUBTRACTION TERM of the mode sum: with the CRP's fourth-order puncture
the `[2]` parameter is already removed by the effective source ([H14]:1786-1789), so
`mmode_reg_Fr(4, orbit, m)` is the leading unsubtracted large-m content of `F_r^m`
(CRP main.tex:520-522, 723) and `|F_r^m − mmode_reg_Fr(4, orbit, m)|` should fall off as
`m^{-6}`. The matching infinite completion is [`mmode_reg_tail_Fr`](@ref).
"""
mmode_reg_Fr(order::Integer, orbit::CircularOrbit{T}, m::Integer) where {T} =
    reg_param_Fr(order, orbit) * mmode_reg_factor(T, order, m)

function mmode_reg_Fr(order::Integer, a::Real, r0::Real, m::Integer;
                      M::Real = one(float(promote_type(typeof(a), typeof(r0)))))
    F = reg_param_Fr(order, a, r0; M = M)
    return F * mmode_reg_factor(typeof(F), order, m)
end

"""
    mmode_reg_tail_Fr(order, a, r0, mstart; M = 1) -> T
    mmode_reg_tail_Fr(order, orbit::CircularOrbit, mstart) -> T

The analytic tail completion `Σ_{m=mstart}^∞ F^m_{r[order]}` in closed form
= [`reg_param_Fr`](@ref)`(order, …) · `[`mmode_reg_factor_tail`](@ref)`(order, mstart)`.

`mode-sum-total-force` (row 46) uses it as `F_r ≈ Σ_{m=0}^{m_max} F_r^m + mmode_reg_tail_Fr(4,
orbit, m_max + 1)`, replacing the `m^{-4}`-convergent remainder of the numerical sum by its
exact value (CRP main.tex:723). Because the doubly-infinite sum of the parameter vanishes,
this equals `−[F^0/2 + Σ_{m=1}^{m_max} F^m]` — see the `m = 0` note in
[`mmode_reg_factor_over_pi`](@ref).
"""
mmode_reg_tail_Fr(order::Integer, orbit::CircularOrbit{T}, mstart::Integer) where {T} =
    reg_param_Fr(order, orbit) * mmode_reg_factor_tail(T, order, mstart)

function mmode_reg_tail_Fr(order::Integer, a::Real, r0::Real, mstart::Integer;
                           M::Real = one(float(promote_type(typeof(a), typeof(r0)))))
    F = reg_param_Fr(order, a, r0; M = M)
    return F * mmode_reg_factor_tail(typeof(F), order, mstart)
end

"""
    mmode_reg_term_Fr(order, orbit::CircularOrbit, m) -> T
    mmode_reg_term_Fr(order, a, r0, m; M = 1) -> T

THE ONE HELPER for the folded/unfolded bookkeeping: the term to subtract from an
Osburn–Nishimura Eq. (3.10) radial force mode `F_r^m`, namely
[`mmode_reg_Fr`](@ref)`(order, …, m)` for `m ≥ 1` and HALF of it at `m = 0`.

WHY THE HALVING (the operational rule stated in [`mmode_reg_factor_over_pi`](@ref), made
executable here so no caller re-derives it): the printed m-dependence is a single analytic
function of `m` with no `m = 0` branch, so the ×2 real-part folding sits INSIDE
`F^m_{r[2j]}` at every `m`, `m = 0` included — while Eq. (3.10)'s `m = 0` branch,
`q ∇_α(Ψ^ℛ_0/r)`, is UNFOLDED. Subtracting the unhalved parameter from an unfolded force mode
would double-count the `m = 0` term.

The rule is equivalent to the vanishing of the doubly-infinite parameter sum,
`F^0/2 + Σ_{m≥1} F^m = 0`, i.e. `mmode_reg_term_Fr(order, orbit, 0) ==
−`[`mmode_reg_tail_Fr`](@ref)`(order, orbit, 1)`, which is how it is gated in
`test/test_regularization.jl`. [`regularized_mode_force`](@ref) (row 44) and the mode sum
(row 46) both go through this function and never apply the factor themselves.
"""
mmode_reg_term_Fr(order::Integer, orbit::CircularOrbit, m::Integer) =
    m == 0 ? mmode_reg_Fr(order, orbit, 0) / 2 : mmode_reg_Fr(order, orbit, m)

function mmode_reg_term_Fr(order::Integer, a::Real, r0::Real, m::Integer;
                           M::Real = one(float(promote_type(typeof(a), typeof(r0)))))
    F = mmode_reg_Fr(order, a, r0, m; M = M)
    return m == 0 ? F / 2 : F
end

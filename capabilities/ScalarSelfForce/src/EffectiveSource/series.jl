# EffectiveSource/series.jl — the NEAR-WORLDLINE POWER SERIES of the m-mode effective source
# and the deterministic route SWITCH (ESSF node `effective-source-series`, W8-B;
# arXiv:2606.04998 main.tex:488-504). This file owns THREE things and nothing else:
#
#   1. the assembly of the paper's eleven series coefficients from the mechanically
#      transcribed `_series_coeffs` (series_coeffs.jl) — in particular the exact
#      half-integer harmonic number `H_{−m−1/2}` that completes S^[2]_0 and S^[2]_2;
#   2. the series evaluation of `S^eff_m(r, z)` in EXACTLY the convention the sealed
#      [`effective_source_m`](@ref) produces (the C's real series scalar carries no rotation,
#      so the full `e^{−im(φ_p + χ(r))}` chart factor of `puncture_psi_chart` applies here);
#   3. THE SWITCH — see the [`effective_source_switched_bundle`](@ref) docstring, the single
#      normative statement of the criterion.
#
# WHY A SERIES AT ALL (main.tex:488-492): the direct assembly cancels four orders of
# singularity, so its roundoff amplification grows like λ^{-2} and Float64 becomes unusable
# below λ ~ 1e-2. In the series the singular terms are cancelled ANALYTICALLY, so its error
# is the truncation O(λ³ log λ) — which vanishes towards the worldline instead of exploding.
# The two routes therefore have opposite λ-trends and cross; the crossover moves inward like
# `eps(T)^{1/5}` (measured log₁₀λ* ≈ −2 at Float64, −5.5 at BigFloat(106), −14 at
# BigFloat(256)), which is why the switch is precision-aware rather than a fixed radius.

"""
    _series_harmonic_mhalf(T, m) -> T

`H_{−m−1/2} = γ + ψ(1/2 − m)`, the harmonic number the C's `HarmonicNumber(-0.5 - m)`
evaluates through GSL's digamma (kerr-circular-comoving.c:638-641).

At half-integer argument the recurrence `ψ(x) = ψ(x+1) − 1/x` terminates on `ψ(1/2) = −γ −
2log2`, so this is an EXACT finite reduction, not a re-implementation of a special function:

    H_{−m−1/2} = 2·Σ_{k=1}^{|m|} 1/(2k−1) − 2 log 2 .

It is used in preference to `BigSpecial.digamma_big` deliberately: that primitive is a fixed
7-term Bernoulli asymptotic series with a ≈1e-16 accuracy floor (see its docstring) and would
cap the BigFloat precision inheritance of every coefficient. The two agree to 1e-13 at
Float64, gated in the node tests.
"""
function _series_harmonic_mhalf(::Type{T}, m::Integer) where {T<:AbstractFloat}
    s = zero(T)
    for k in 1:abs(m)
        s += one(T) / (2k - 1)
    end
    return 2 * s - 2 * log(T(2))
end

"""
    effective_source_series_coefficients(p::Puncture, m::Integer) -> NamedTuple

The ELEVEN coefficients of the near-worldline series of `S^eff_m` (arXiv:2606.04998
main.tex:493-504), for the orbit `p.orbit` and mode `|m|`:

    s(Δr, Δθ) = S0 + S1_1·Δr
              + [ S2_2·Δr² + S2_0·ϱ² + (S2_2L·Δr² + S2_0L·ϱ²)·log ϱ
                  + S2_4·Δr⁴/ϱ² + S2_6·Δr⁶/ϱ⁴ + S2_8·Δr⁸/ϱ⁶
                  + S2_10·Δr¹⁰/ϱ⁸ + S2_12·Δr¹²/ϱ¹⁰ ]  + O(λ³ log λ) ,

with `ϱ² = orbit.rho2_dr2·Δr² + orbit.rho2_dtheta2·Δθ²` (eq:s2). `s` is the C's REAL
`effsource_calc_m_series` scalar — `□_m` applied to the bare azimuthal integral — and
[`effective_source_series_m`](@ref) maps it to `S^eff_m` through the convention chain.

TRANSCRIPTION. The algebra lives in the GENERATED `series_coeffs.jl` (`_series_coeffs`,
mechanically transcribed from `effsource_calc_m_series` with the λ-structure divided out);
this function only adds the harmonic/log regrouping. The C writes that block as
`(X0 ϱ² + X2 Δr²)·(H_{−m−1/2} + log(ϱ/z_c))`, so

    S2_0L = X0 ,  S2_2L = X2 ,  S2_0 = S2_0rat + X0·(H_{−m−1/2} − log z_c) ,

and likewise for S2_2 — i.e. the `−log z_c` remainder is absorbed into the non-log
coefficients, exactly as the paper's form requires.

m-DEPENDENCE (main.tex:503-504, asserted exactly in `Rational{BigInt}` by the node tests and
by the generator on the transcription itself): `S0` and `S1_1` carry `m` ONLY through
`1/((2m−1)(2m+1))`; `S2_0` and `S2_2` are a degree-1 polynomial in `m²` over that same
factor PLUS a `H_{−m−1/2}` term; the other seven are m-independent. Everything reads `|m|`
(the package's `Ψ_{−m} = conj(Ψ_m)` convention), and every coefficient is exactly linear in
`z_c`. Type-generic; a BigFloat orbit gives BigFloat coefficients.
"""
function effective_source_series_coefficients(p::Puncture, m::Integer)
    o = p.orbit
    T = typeof(o.r0)
    zc = sqrt(o.zc2)
    c = _series_coeffs(o.M, o.a, o.r0, o.Omega, zc, abs(m))
    HL = _series_harmonic_mhalf(T, m) - log(zc)
    return (S0 = c.S0, S1_1 = c.S1_1,
            S2_2 = c.S2_2rat + c.S2_2L * HL, S2_0 = c.S2_0rat + c.S2_0L * HL,
            S2_2L = c.S2_2L, S2_0L = c.S2_0L,
            S2_4 = c.S2_4, S2_6 = c.S2_6, S2_8 = c.S2_8,
            S2_10 = c.S2_10, S2_12 = c.S2_12)
end

# The C series scalar → `S^eff_m` chart factor. The C's `src` is real and carries NO azimuthal
# rotation (its `cos/sin(m φ_p)` pair is stripped by the caller), so the FULL horizon-regular
# field-point rotation of `puncture_psi_chart` applies: `e^{−im(φ_p + χ(r))}`, χ =
# kerr_varphi_offset. The remaining factor is the `−Σ/(r²+a²)` mode-projection scale of
# eq:elliptic_eqn times the ansatz stripping `r/(2π sin^{|m|}θ)` — spectre
# CircularOrbit.cpp:322-332, the same composition the direct route reaches through
# `mode_amplitude` (both are gated against the baked C fixtures).
function _series_chart_factor(o::CircularOrbit, m::Integer, r, z, phi_p)
    Sigma = r^2 + o.a^2 * z^2
    return cis(-m * kerr_varphi(o.a, r, phi_p)) * r * Sigma /
           (2 * oftype(r, pi) * (r^2 + o.a^2) * sinm_theta(m, z))
end

# ONE evaluation of the series: the complex source, the relative size `q` of the retained
# λ² bracket (the switch's raw material) and ϱ².
function _series_eval(p::Puncture, m::Integer, r::Real, z::Real, phi_p::Real)
    o = p.orbit
    abs(z) < 1 || throw(DomainError(z,
        "effective_source_series_m: |z| < 1 required (z = cosθ; the sin^{|m|}θ ansatz factor \
         of the mode-projection chain vanishes at the poles)"))
    dr = r - o.r0
    dtheta = -asin(z)                                    # Δθ = θ − π/2 (locked branch)
    rho2 = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dtheta^2
    A = effective_source_series_coefficients(p, m)
    dr2 = dr^2
    # the O(λ²) bracket; ϱ = 0 is the worldline circle, where the C guards its λ-dependent
    # orders to zero (`(dr==0&&dth==0)? 0 :`) and the series reduces to S^[0]
    br = iszero(rho2) ? zero(rho2) :
         (A.S2_2 * dr2 + A.S2_0 * rho2 + (A.S2_2L * dr2 + A.S2_0L * rho2) * log(rho2) / 2 +
          A.S2_4 * dr2^2 / rho2 + A.S2_6 * dr2^3 / rho2^2 + A.S2_8 * dr2^4 / rho2^3 +
          A.S2_10 * dr2^5 / rho2^4 + A.S2_12 * dr2^6 / rho2^5)
    s = A.S0 + A.S1_1 * dr + br
    q = abs(br) / abs(s)
    return (source = s * _series_chart_factor(o, m, r, z, phi_p),
            q = isfinite(q) ? q : oftype(rho2, Inf), rho2 = rho2)
end

"""
    effective_source_series_m(p::Puncture, m::Integer, r, z; phi_p = 0) -> Complex

The NEAR-WORLDLINE SERIES approximation to the m-mode effective source `S^eff_m(r, z)`
(arXiv:2606.04998 main.tex:493-504) — the same quantity the sealed [`effective_source_m`](@ref)
computes directly, in the same convention, but with the four cancelling orders of singularity
removed ANALYTICALLY instead of numerically.

Accurate through order `λ²`; the remainder is `O(λ³ log ϱ)` (the fitted exponent is
`3.00 ± 0.05` against a BigFloat(256) direct evaluation over `λ ∈ [1e-7, 1e-4]` at every
config and `m` of the frozen fixture grid). Use it INSTEAD of the direct route wherever the
direct route's `cancellation` diagnostic says its roundoff exceeds this truncation — the
[`effective_source_switched_bundle`](@ref) criterion, which callers should normally use
rather than choosing by hand.

DOMAIN. Unlike the direct route this is FINITE ON THE WORLDLINE: at `ϱ = 0` it returns
exactly `S^[0]` times the chart factor (`S^eff_m` is finite there — it is the puncture FIELD
that diverges logarithmically, not its effective source). `|z| < 1` is still required and
raises `DomainError`: the `sin^{|m|}θ` ansatz stripping of the mode-projection chain is
singular at the poles. `phi_p` is the particle's BOYER–LINDQUIST azimuth, exactly as for
[`effective_source_m`](@ref); the horizon-regular field-point rotation is applied internally
and must never be applied twice. Type-generic and allocation-free on the Float64 path.
"""
effective_source_series_m(p::Puncture, m::Integer, r::Real, z::Real; phi_p::Real = 0) =
    _series_eval(p, m, r, z, phi_p).source

# The measured roundoff amplification of the DIRECT assembly, in units of its own
# `eps·cancellation` diagnostic. `cancellation = term_scale/|source|` counts ONE ulp per
# assembled term; the seven terms, their own inputs and the λ^{-2}-amplified inputs of
# `puncture_d2phi_m` push the realised loss above that. MEASURED over the node's 3240-point
# (r, z) calibration sweep, restricted to the 2148 points where the direct value still has
# ≥ 3 correct digits (so that its `cancellation` is itself meaningful): median 91, 90th
# percentile 437, worst 2310 — and 307 on the sealed direct node's own near-puncture box.
# The constant enters the crossover only as `ULPS^{1/5}` and is NOT sensitive: 1000, 2310
# and 3000 all give the SAME worst switched relative error on that sweep (6.95e-4), with
# worst achieved/optimal factors 42.1, 36.3 and 36.3. 1000 is kept as the working value
# because it also gives the smallest near-crossover route disagreement (1.2e-3 vs 4.0e-2).
const _ESS_DIRECT_ULPS = 1000

"""
    effective_source_switched_bundle(p::Puncture, m::Integer, r, z; phi_p = 0)
        -> (source, route, err_series, err_direct)

THE COMBINED EVALUATOR of `S^eff_m(r, z)` and the single normative statement of the near-worldline
ROUTE SWITCH (arXiv:2606.04998 main.tex:488-492). `route` is `:series` or `:direct`; `source`
is that route's value; `err_series` and `err_direct` are the two PREDICTED RELATIVE errors the
decision compares.

THE CRITERION — pick the route with the smaller predicted error:

    err_direct = _ESS_DIRECT_ULPS · eps(T) · cancellation            (roundoff)
    err_series = q · max(q, (1 + |m|)·ϱ/r₀)                          (truncation)
    route      = err_series ≤ err_direct ? :series : :direct

with `cancellation` the sealed [`effective_source_bundle`](@ref) diagnostic (`term_scale/|S|`,
the λ^{-2} amplification of the four cancelled orders) and `q = |bracket|/|s|` the relative
size of the LAST RETAINED order of the series. `q` estimates `(λ/L)²` for the local expansion
scale `L`, so `q·(1+|m|)λ/r₀` is the next (λ³) order; the `max` with `q` itself is the
non-asymptotic guard — where the retained order is not yet small, the expansion is not
converging and the series must not be trusted whatever λ says.

The rule is precision-aware BY CONSTRUCTION: `err_series ~ λ³` and `err_direct ~ eps·λ^{-2}`
cross at `λ* ~ eps(T)^{1/5}`, so raising the precision moves the switch INWARD — measured
log₁₀λ* = −1.923 at Float64, −5.201 at BigFloat(106) and −14.320 at BigFloat(256) (a = 0.5,
r₀ = 10, m = 2); it also moves inward with `m` and with decreasing r₀.

MEASURED QUALITY (node evidence `calibration_switch.log`: a 3240-point sweep, 6 fixture
configs × m ∈ {0,1,2,5,10,20} × 9 λ ∈ [1e-5, 1e-1] × 5 directions × both signs of Δr, against
a BigFloat(256) truth). The switched evaluator's worst relative error is 6.9542e-4, against
5.8608e-4 for an ORACLE that always picks the better route (ratio 1.187) — while on the same
points the direct route ALONE reaches 6.0385e+2 and the series ALONE 4.0122e+0. The worst
pointwise achieved/optimal ratio is 42.1, and for all 1080 points with λ ≤ 3.2e-4 the rule
picks the series with ratio exactly 1.

CONTINUITY. The two routes are not identical, so `source` jumps by |series − direct| at the
crossover. That jump is BY CONSTRUCTION of the order of the local accuracy limit itself
(both predicted errors are equal there): measured worst |series − direct|/|S| = 3.63e-4 at
Float64 over the node test's grid and 1.2e-3 over the full calibration sweep, both at m = 20
and λ ≈ 1e-2, where NEITHER route does better than ~1e-4. It is a jump in the DG fixed
source at a measure-zero set of quadrature points inside the worldtube; it is NOT gated below
the DG discretization error, because at Float64 nothing about this source is more accurate
than that jump — a caller needing continuity below it must raise the working precision, which
moves the crossover (and the jump) down like `eps^{1/5}`.

DOMAIN. `|z| < 1` (`DomainError`, as for both routes). On the worldline circle `ϱ = 0` the
direct route is undefined and the series is exact in the limit, so the route is `:series`
with `err_direct = Inf` — this evaluator, unlike [`effective_source_m`](@ref), is finite
there. Type-generic (`eps(T)` makes the decision precision-aware); allocation-free on the
Float64 path.
"""
function effective_source_switched_bundle(p::Puncture, m::Integer, r::Real, z::Real;
                                          phi_p::Real = 0)
    f = _series_eval(p, m, r, z, phi_p)
    T = typeof(f.rho2)
    err_series = f.q * max(f.q, (1 + abs(m)) * sqrt(f.rho2) / p.orbit.r0)
    if iszero(f.rho2)
        return (source = f.source, route = :series, err_series = err_series,
                err_direct = T(Inf))
    end
    b = effective_source_bundle(p, m, r, z; phi_p = phi_p)
    err_direct = _ESS_DIRECT_ULPS * eps(T) * b.cancellation
    if err_series <= err_direct
        return (source = f.source, route = :series, err_series = err_series,
                err_direct = err_direct)
    end
    return (source = b.source, route = :direct, err_series = err_series,
            err_direct = err_direct)
end

"""
    effective_source_m_switched(p::Puncture, m::Integer, r, z; phi_p = 0) -> Complex

`S^eff_m(r, z)` through the near-worldline route SWITCH — THE evaluator the DG fixed-source
assembly should call (`+S^eff_m` inside the worldtube, zero outside; the sign contract of
node `effective-source-direct`). Equal to
`effective_source_switched_bundle(p, m, r, z; phi_p).source`; that docstring states the
criterion, its measured quality and the crossover behaviour.
"""
effective_source_m_switched(p::Puncture, m::Integer, r::Real, z::Real; phi_p::Real = 0) =
    effective_source_switched_bundle(p, m, r, z; phi_p = phi_p).source

# Self-contained confluent hypergeometric ₁F₁(a;b;z) (Kummer's function M) — the special-function
# primitive the Teukolsky-MST / Leaver radial solutions need WITHOUT live-depping the heavy
# Teukolsky/FEW tree. Genuinely new (the Teukolsky hypergeom layer only carries ₂F₁ / U); the
# incremental-Pochhammer term recurrence mirrors the structural pattern of that file's `_f21_maclaurin`
# so it is exact at ANY precision (no factorials, no Γ/ψ). Type-generic so BigFloat is inherited.

"""
    _kummer_series(a, b, z; maxterms=1_000_000) -> T

Raw Maclaurin partial sum of `₁F₁(a;b;z) = Σ_{n≥0} (a)_n/((b)_n n!) zⁿ` via the incremental term
ratio `tₙ₊₁ = tₙ·(a+n)/((b+n)(n+1))·z`, summed until `|tₙ| ≤ ε·|s|` (working-precision `ε`). No
Γ/factorials, so it is exact at any precision; the `(a)ₙ` Pochhammer makes it terminate naturally when
`a` is a non-positive integer. NOTE: for `Re z < 0` the terms alternate and cancel — call `kummer_1F1`
(which applies Kummer's transform) rather than this raw kernel there.
"""
function _kummer_series(a, b, z; maxterms::Integer=1_000_000)
    T = promote_type(typeof(a), typeof(b), typeof(z))
    aT = T(a); bT = T(b); zT = T(z)
    tol = eps(real(T))
    term = one(T); s = term; n = 0
    while n < maxterms
        term *= (aT + n) / ((bT + n) * (n + 1)) * zT
        s += term
        abs(term) <= tol * abs(s) && break
        n += 1
    end
    return s
end

# ---------------------------------------------------------------------------------------------------
# Large-|z| asymptotics: the Γ-free ₂F₀ kernel + the Tricomi-U (recessive) bridge + the dominant M term
# ---------------------------------------------------------------------------------------------------

_is_nonpos_int(x::Real) = isinteger(x) && x <= 0
_is_nonpos_int(x::Complex) = imag(x) == 0 && _is_nonpos_int(real(x))

"""
    _hyper2F0_asymp(α, β, w) -> T

Optimally-truncated divergent hypergeometric `₂F₀(α,β;;w) = Σ_{n≥0} (α)_n(β)_n/n! wⁿ` via the
incremental-Pochhammer term ratio `tₙ₊₁ = tₙ·(α+n)(β+n)/(n+1)·w` (no Γ/factorials). When `α` or `β`
is a non-positive integer the series TERMINATES (the `(α)_n`/`(β)_n` Pochhammer hits zero) and the
finite sum is EXACT. Otherwise the series is asymptotic (diverges for every `w≠0`): it is summed to
the SMALLEST term and truncated there (superasymptotic), giving relative accuracy `~e^{−1/|w|}`.
Shared kernel of both the `M` and `U` large-argument asymptotics. Type-generic (BigFloat-stable).
"""
function _hyper2F0_asymp(α, β, w)
    T = promote_type(typeof(α), typeof(β), typeof(w))
    αT = T(α); βT = T(β); wT = T(w)
    terminating = _is_nonpos_int(α) || _is_nonpos_int(β)
    tol = eps(real(T))
    term = one(T); s = term; prevabs = abs(term)
    nmax = terminating ? 1_000_000 : 100_000
    n = 0
    while n < nmax
        term *= (αT + n) * (βT + n) / (n + 1) * wT
        na = abs(term)
        if terminating
            s += term
            na == 0 && break
        else
            na > prevabs && break          # past the smallest term ⇒ discard the divergent tail
            s += term
            na <= tol * abs(s) && break
        end
        prevabs = na
        n += 1
    end
    return s
end

"""
    kummer_U(a, b, z) -> T

Tricomi confluent hypergeometric `U(a,b,z)` (the recessive solution of Kummer's equation), via the
Γ-FREE leading large-`z` asymptotic `U(a,b,z) ~ z^{−a}·₂F₀(a, a−b+1;; −1/z)` (DLMF 13.7.3). Exact
(finite `₂F₀`) for the degenerate families `U(a,a+1,z)=z^{−a}`, `U(0,b,z)=1`, and the truncating
`U(−n,b,z)`; asymptotically accurate to `~e^{−|z|}` otherwise (valid for large `|z|`, `Re z ≥ 0`).
The DLMF 13.2.40 transform `U(a,b,z)=z^{1−b}U(a−b+1,2−b,z)` is an equivalent route (round-trip-verified)
covering the small-|z| / integer-`b` branch. Type-generic; no Γ, so BigFloat is inherited.
"""
function kummer_U(a, b, z)
    T = promote_type(typeof(a), typeof(b), typeof(z))
    aT = T(a); bT = T(b); zT = T(z)
    return zT^(-aT) * _hyper2F0_asymp(aT, aT - bT + one(T), -one(T) / zT)
end

"""
    _kummer_asymp(a, b, z) -> T

Dominant large-`z` asymptotic of `M(a,b,z) = ₁F₁(a;b;z)` (DLMF 13.7.2),
`M(a,b,z) ~ Γ(b)/Γ(a) eᶻ z^{a−b} ₂F₀(b−a,1−a;;1/z)`, reusing [`_hyper2F0_asymp`](@ref). The `eᶻ`
term dominates for `Re z → +∞` (the recessive `z^{−a}`/`U` term is exponentially subdominant and
dropped). The `Γ(b)/Γ(a)` prefactor is Γ-FREE only for the ratio-1 anchor families it is routed for:
`b=a` (⇒ `₂F₀(0,·)=1`, `z^0` ⇒ `eᶻ` exactly) and `a=1,b=2` (⇒ `₂F₀(1,0)=1`, `z^{−1}` ⇒ `eᶻ/z`).
Callers must gate on [`_is_gammaratio_anchor`](@ref) so the sink stays Γ-free.
"""
function _kummer_asymp(a, b, z)
    T = promote_type(typeof(a), typeof(b), typeof(z))
    aT = T(a); bT = T(b); zT = T(z)
    return exp(zT) * zT^(aT - bT) * _hyper2F0_asymp(bT - aT, one(T) - aT, one(T) / zT)
end

"""`_is_gammaratio_anchor(a,b)` — the `M` asymptotic families with `Γ(b)/Γ(a)=1` (Γ-free): `b=a` or `a=1,b=2`."""
_is_gammaratio_anchor(a, b) = (a == b) || (a == oneunit(a) && b == 2 * oneunit(b))

# Re z threshold beyond which the dropped recessive term (~e^{−Re z}) is below working precision.
_asymp_z_threshold(::Type{T}) where {T} = -log(eps(real(T))) + 2

"""
    kummer_1F1(a, b, z) -> T

Confluent hypergeometric `₁F₁(a;b;z)` (Kummer's `M`), type-generic and BigFloat-stable, valid in ALL
`z` regimes:
  - `Re z < 0`: Kummer's transform `₁F₁(a;b;z) = eᶻ ₁F₁(b−a;b;−z)` moves to the monotone-term series
    (the direct Maclaurin cancels catastrophically there);
  - large `|z|`, `Re z ≥ 0`, ratio-1 anchor family (`b=a` or `a=1,b=2`): the Γ-free dominant asymptotic
    [`_kummer_asymp`](@ref) `~ eᶻ z^{a−b} ₂F₀(b−a,1−a;;1/z)`, where the recessive `z^{−a}` (`U`-function)
    term is exponentially subdominant — routed only past `_asymp_z_threshold` so that dropped term is
    below working precision;
  - otherwise (`Re z ≥ 0`): the direct Maclaurin series (cancellation-free for `Re z ≥ 0`).
`b` must not be a non-positive integer (pole of `M`). The Tricomi partner is [`kummer_U`](@ref).
"""
function kummer_1F1(a, b, z)
    T = promote_type(typeof(a), typeof(b), typeof(z))
    aT = T(a); bT = T(b); zT = T(z)
    if real(zT) < 0
        return exp(zT) * _kummer_series(bT - aT, bT, -zT)
    end
    if real(zT) > _asymp_z_threshold(T) && _is_gammaratio_anchor(aT, bT)
        return _kummer_asymp(aT, bT, zT)
    end
    return _kummer_series(aT, bT, zT)
end

# series.jl — truncated Laurent series over ℚ, and a direct evaluator for Symbolics expressions.
#
# Why this exists: the Newtonian limit of the derived entries is a LIMIT s → 0, and doing it by
# symbolic substitution of truncated Taylor polynomials followed by `expand`/`simplify_fractions`
# blows the expressions up (measured: >10 GB resident, no result in 12 min — recorded as a finding of
# the node).  Evaluating the SAME expression tree in the ring ℚ((s)) of truncated Laurent series is
# exact, allocation-cheap, and yields the Laurent order and the leading coefficients directly.
#
# `eval_series` walks the SymbolicUtils tree with only +, *, /, ^(integer) — the operations the whole
# derivation is built from — so it doubles as an assertion that the entries really are rational
# functions of the alphabet (any transcendental would raise instead of being approximated).

"Number of retained coefficients."
const SPREC = 16
"`n0` of a series that is EXACTLY zero (known to vanish at every order)."
const SZERO_N0 = 1_000_000

"""
    LSeries(n0, c)

Σ_i c[i] s^{n0+i-1} + O(s^{n0+length(c)}).  `c[1] ≠ 0` after normalisation; an empty `c` means "known
to vanish below order `n0`, unknown above" (with `n0 = SZERO_N0` for an exact zero).
"""
struct LSeries <: Number
    n0::Int
    c::Vector{RatB}
end

_lim(a::LSeries) = a.n0 + length(a.c)

function _mk(n0::Int, c::Vector{RatB})
    i = findfirst(!iszero, c)
    i === nothing && return LSeries(n0 + length(c), RatB[])
    LSeries(n0 + i - 1, c[i:min(length(c), i - 1 + SPREC)])
end

series_zero() = LSeries(SZERO_N0, RatB[])
_pad(q) = RatB[i == 1 ? RatB(q) : zero(RatB) for i in 1:SPREC]
series_const(q) = iszero(q) ? series_zero() : LSeries(0, _pad(q))
series_mono(q, k::Integer) = iszero(q) ? series_zero() : LSeries(Int(k), _pad(q))

"Laurent order (lowest order with a nonzero coefficient); `nothing` if the series is known-zero."
lorder(a::LSeries) = isempty(a.c) ? nothing : a.n0

Base.zero(::Type{LSeries}) = series_zero()
Base.one(::Type{LSeries}) = series_const(one(RatB))

function Base.:+(a::LSeries, b::LSeries)
    (isempty(a.c) && isempty(b.c)) && return LSeries(min(a.n0, b.n0), RatB[])
    n0 = min(a.n0, b.n0)
    lim = min(_lim(a), _lim(b))
    lim <= n0 && return LSeries(lim, RatB[])
    c = zeros(RatB, lim - n0)
    for i in eachindex(a.c)
        k = a.n0 + i - n0
        1 <= k <= length(c) && (c[k] += a.c[i])
    end
    for i in eachindex(b.c)
        k = b.n0 + i - n0
        1 <= k <= length(c) && (c[k] += b.c[i])
    end
    _mk(n0, c)
end

Base.:-(a::LSeries) = LSeries(a.n0, -a.c)
Base.:-(a::LSeries, b::LSeries) = a + (-b)

function Base.:*(a::LSeries, b::LSeries)
    (isempty(a.c) || isempty(b.c)) && return LSeries(min(a.n0 + _lim(b), b.n0 + _lim(a)), RatB[])
    n0 = a.n0 + b.n0
    L = min(length(a.c), length(b.c))
    c = zeros(RatB, L)
    for i in eachindex(a.c), j in eachindex(b.c)
        k = i + j - 1
        k <= L && (c[k] += a.c[i] * b.c[j])
    end
    _mk(n0, c)
end

function Base.inv(a::LSeries)
    isempty(a.c) && error("inv(LSeries): the series is known only as zero")
    L = length(a.c)
    b = zeros(RatB, L)
    b[1] = inv(a.c[1])
    for k in 2:L
        acc = zero(RatB)
        for j in 2:k
            acc += a.c[j] * b[k-j+1]
        end
        b[k] = -b[1] * acc
    end
    _mk(-a.n0, b)
end

Base.:/(a::LSeries, b::LSeries) = a * inv(b)

function Base.:^(a::LSeries, n::Integer)
    n == 0 && return one(LSeries)
    n < 0 && return inv(a)^(-n)
    r = a
    for _ in 2:n
        r = r * a
    end
    r
end

for op in (:+, :-, :*, :/)
    @eval Base.$op(a::LSeries, x::Union{Integer,Rational}) = $op(a, series_const(RatB(x)))
    @eval Base.$op(x::Union{Integer,Rational}, a::LSeries) = $op(series_const(RatB(x)), a)
end

"exp of a series of strictly positive order (Σ aᵏ/k!, exact in ℚ)."
function sexp(a::LSeries)
    isempty(a.c) && return one(LSeries)
    a.n0 >= 1 || error("sexp: needs order ≥ 1")
    r = one(LSeries)
    term = one(LSeries)
    for k in 1:SPREC
        term = term * a / k
        isempty(term.c) && break
        r = r + term
        term.n0 > SPREC && break
    end
    r
end

"(1 − x)^(−1/2) for a series `x` of order ≥ 1: Σ_k binomial(2k,k)/4ᵏ xᵏ, exact in ℚ."
function sinvsqrt1m(x::LSeries)
    isempty(x.c) && return one(LSeries)
    x.n0 >= 1 || error("sinvsqrt1m: needs order ≥ 1")
    r = one(LSeries)
    xp = one(LSeries)
    for k in 1:SPREC
        xp = xp * x
        isempty(xp.c) && break
        r = r + xp * (RatB(binomial(big(2k), big(k))) // RatB(big(4)^k))
        xp.n0 > SPREC && break
    end
    r
end

# ---------------------------------------------------------------------------------------------
# direct evaluation of a Symbolics expression tree in the series ring
# ---------------------------------------------------------------------------------------------

"""
    eval_series(e, env) -> LSeries  (or `(re, im)` for a complex expression)

Evaluate the expression tree of `e` with each symbol bound to the `LSeries` in `env` (keyed by
`Symbolics.value(sym)`).  Only `+`, `*`, `/` and integer powers are accepted.
"""
eval_series(e::Num, env) = _evs(Symbolics.value(e), env)
eval_series(e::Complex{Num}, env) = (eval_series(real(e), env), eval_series(imag(e), env))

function _evs(x, env)
    x isa LSeries && return x
    x isa Number && return series_const(RatB(x))
    haskey(env, x) && return env[x]
    # SymbolicUtils wraps numeric literals in a `Const` leaf (`-1`, `2//3`, …) that is neither a
    # Julia `Number` nor a `Sym`; unwrap it with `Symbolics.value`.
    SymbolicUtils.isconst(x) && return series_const(RatB(Symbolics.value(x)))
    if SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        as = SymbolicUtils.arguments(x)
        op === (+) && return reduce(+, (_evs(a, env) for a in as))
        op === (*) && return reduce(*, (_evs(a, env) for a in as))
        op === (/) && return _evs(as[1], env) / _evs(as[2], env)
        if op === (^)
            b = _evs(as[1], env)
            n = as[2] isa Number ? as[2] : Symbolics.value(as[2])
            (n isa Number && isinteger(n)) || error("eval_series: non-integer exponent $(as[2])")
            return b^Int(n)
        end
        if op === (-)
            return length(as) == 1 ? -_evs(as[1], env) : _evs(as[1], env) - _evs(as[2], env)
        end
        error("eval_series: unsupported operation $(op)")
    end
    error("eval_series: unbound symbol $(x)")
end

"Laurent order of a real or complex expression under `env` (`nothing` ⇒ known-zero)."
function expr_order(e, env)
    r = eval_series(e, env)
    r isa LSeries && return lorder(r)
    a, b = r
    oa, ob = lorder(a), lorder(b)
    oa === nothing && return ob
    ob === nothing && return oa
    min(oa, ob)
end

"Known-limit (exclusive) of the expansion of `e` under `env` — how far a comparison is meaningful."
function expr_limit(e, env)
    r = eval_series(e, env)
    r isa LSeries && return _lim(r)
    a, b = r
    min(_lim(a), _lim(b))
end

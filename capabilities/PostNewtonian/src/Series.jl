"""
    PostNewtonian.Series

The `PNSeries` Num-algebra engine: truncated Laurent/Taylor series in a PN parameter whose
coefficients are `Symbolics.Num` — the Julia mirror of Mathematica `SeriesData`. Provides arithmetic
(`+ - * / inv ^`), `Normal` materialization, `polyToSeries`, `series_take/collect/terms/of`,
change-of-parameter, power-counting, and `expand_log`/`expand_gamma`.

`SeriesData[x, x0, {a_0..a_k}, nmin, nmax, den]` represents
`Sum_{i=0}^{k} a_i (x-x0)^((nmin+i)/den) + O[(x-x0)^(nmax/den)]`, with `k = (nmax-nmin)/den - 1`
stored coefficients. This is the foundation every other PN submodule builds on; `Symbolics` is its
only algebra dep (`SpecialFunctions.gamma` is the rewrite target of `expand_gamma`).
"""
module Series

using Symbolics
using Symbolics: Num, value, expand, simplify, substitute
import SpecialFunctions

export PNSeries, Normal, svar,
    series_min_order, series_max_order, series_length,
    series_take, series_collect, series_of, series_terms,
    polyToSeries, change_series_parameter, power_counting,
    expand_log, expand_gamma

# ----------------------------------------------------------------------------
# Type
# ----------------------------------------------------------------------------

"""
    PNSeries(var, x0, coeffs, nmin, nmax, den)

Truncated Laurent series in `var` about `x0`: mirror of MMA `SeriesData`. `coeffs[i]` (1-based) is
the coefficient of `(var-x0)^((nmin+i-1)/den)`; the series carries `O[(var-x0)^(nmax/den)]`.
Invariant `length(coeffs) == nmax-nmin` (the consecutive num-exponents `nmin..nmax-1`).
"""
struct PNSeries
    var::Symbol
    x0::Num
    coeffs::Vector{Num}
    nmin::Int
    nmax::Int
    den::Int
    function PNSeries(var::Symbol, x0, coeffs::AbstractVector, nmin::Integer, nmax::Integer, den::Integer)
        den > 0 || throw(ArgumentError("den must be positive, got $den"))
        c = Num[Num(x) for x in coeffs]
        # MMA invariant: coeffs sit at CONSECUTIVE num-exponents nmin..nmax-1 (step 1 in num-units;
        # real order = num/den). Stored count == nmax-nmin, the O[] term is at real order nmax/den.
        nslots = nmax - nmin
        nslots >= 0 || throw(ArgumentError("nmax ($nmax) < nmin ($nmin)"))
        if length(c) > nslots
            c = c[1:nslots]            # MMA never stores coeffs at/after the O[] order
        elseif length(c) < nslots
            c = vcat(c, fill(Num(0), nslots - length(c)))
        end
        new(var, Num(x0), c, Int(nmin), Int(nmax), Int(den))
    end
end

# convenience: build with x0 defaulting to 0
PNSeries(var::Symbol, coeffs::AbstractVector, nmin::Integer, nmax::Integer, den::Integer) =
    PNSeries(var, 0, coeffs, nmin, nmax, den)

"""
    svar(s::PNSeries) -> Num

The expansion variable of `s` as a `Symbolics.Num` (`Symbolics.variable(s.var)`), so callers can
form `(var - x0)`-type expressions from the stored `Symbol`.
"""
svar(s::PNSeries) = Symbolics.variable(s.var)

# real-valued leading/trailing orders
"""
    series_min_order(s::PNSeries) -> Rational

Leading (lowest stored) real order of `s`, `nmin // den`. Mirror of MMA `SeriesMinOrder`. The
scalar method `series_min_order(::Int) = 0` matches MMA's convention when `SeriesTerms` is applied
to a plain scalar.
"""
series_min_order(s::PNSeries) = s.nmin // s.den

"""
    series_max_order(s::PNSeries) -> Rational

Real order of the `O[]` truncation term, `nmax // den` — the first order NOT stored. Mirror of MMA
`SeriesMaxOrder`. Distinct from `SelfForce.pn_max_order`, which is `(nmax - 1) // den`, the last
RETAINED order.
"""
series_max_order(s::PNSeries) = s.nmax // s.den

"""
    series_length(s::PNSeries) -> Rational

Order-span of the stored coefficient window, `series_max_order(s) - series_min_order(s)`
(= `length(s.coeffs) // s.den` by the `PNSeries` invariant `length(coeffs) == nmax - nmin`).
"""
series_length(s::PNSeries)    = series_max_order(s) - series_min_order(s)
series_min_order(::Int) = 0   # MMA: SeriesMinOrder[1]=0 (used by SeriesTerms on a scalar)

# ----------------------------------------------------------------------------
# Normal: drop O[] and return the Symbolics polynomial/Laurent expression
# ----------------------------------------------------------------------------

"""
    Normal(s) -> Num

The polynomial (Laurent) expression `Sum_i coeffs[i] (var-x0)^((nmin+i-1)/den)`, dropping the O[]
term. Mirror of MMA `Normal[SeriesData[...]]`. Scalars pass through. Extended by
`SelfForce.ResummedSeriesData` over the same method table.
"""
function Normal(s::PNSeries)
    x = svar(s)
    base = x - s.x0
    acc = Num(0)
    for (i, c) in enumerate(s.coeffs)
        iszero(c) && continue
        e = (s.nmin + i - 1) // s.den
        acc += c * base^e
    end
    acc
end
Normal(x) = x  # scalars/Num pass through

# ----------------------------------------------------------------------------
# Common-denominator helpers for arithmetic
# ----------------------------------------------------------------------------

# Re-express s on denominator `den`: return (step-1 coeff grid in new-den units, newmin, newmax, den).
function _on_den(s::PNSeries, den::Int)
    f = div(den, s.den)                 # den = s.den * f, f>=1
    newmin = s.nmin * f
    newmax = s.nmax * f
    n = length(s.coeffs)
    grid = Vector{Num}(undef, n == 0 ? 0 : (n - 1) * f + 1)
    if n > 0
        fill!(grid, Num(0))
        for i in 1:n
            grid[(i - 1) * f + 1] = s.coeffs[i]
        end
    end
    (grid, newmin, newmax, den)
end

_lcm(a::Int, b::Int) = lcm(a, b)

# ----------------------------------------------------------------------------
# Addition / subtraction (result truncation = min of the two O[] orders)
# ----------------------------------------------------------------------------

function Base.:+(a::PNSeries, b::PNSeries)
    a.var == b.var || throw(ArgumentError("cannot add series in different variables ($(a.var) vs $(b.var))"))
    isequal(a.x0, b.x0) || throw(ArgumentError("cannot add series about different points"))
    den = _lcm(a.den, b.den)
    (ca, amin, amax, _) = _on_den(a, den)
    (cb, bmin, bmax, _) = _on_den(b, den)
    nmin = min(amin, bmin)
    nmax = min(amax, bmax)             # truncate at the smaller O[] order
    nslots = nmax - nmin
    out = fill(Num(0), max(nslots, 0))
    for (coeffs, mn) in ((ca, amin), (cb, bmin))
        for (i, c) in enumerate(coeffs)
            e = mn + (i - 1)
            e >= nmax && break
            e < nmin && continue
            out[e - nmin + 1] += c
        end
    end
    PNSeries(a.var, a.x0, out, nmin, nmax, den)
end

Base.:-(s::PNSeries) = PNSeries(s.var, s.x0, [-c for c in s.coeffs], s.nmin, s.nmax, s.den)
Base.:-(a::PNSeries, b::PNSeries) = a + (-b)

# ----------------------------------------------------------------------------
# Scalar / Num arithmetic (a scalar == c + O[var]^Infinity: does NOT lower the truncation)
# ----------------------------------------------------------------------------

const ScalarLike = Union{Number, Num}

function Base.:+(s::PNSeries, c::ScalarLike)
    cN = Num(c)
    iszero(cN) && return s
    den = s.den
    nmin = min(s.nmin, 0)
    nmax = s.nmax
    nmax > nmin || throw(ArgumentError("scalar add would have empty range"))
    nslots = nmax - nmin
    out = fill(Num(0), nslots)
    for (i, cc) in enumerate(s.coeffs)
        e = s.nmin + (i - 1)
        out[e - nmin + 1] += cc
    end
    if nmin <= 0 < nmax
        out[0 - nmin + 1] += cN
    end
    PNSeries(s.var, s.x0, out, nmin, nmax, den)
end
Base.:+(c::ScalarLike, s::PNSeries) = s + c
Base.:-(s::PNSeries, c::ScalarLike) = s + (-Num(c))
Base.:-(c::ScalarLike, s::PNSeries) = (-s) + c

Base.:*(s::PNSeries, c::ScalarLike) =
    PNSeries(s.var, s.x0, [Num(c) * cc for cc in s.coeffs], s.nmin, s.nmax, s.den)
Base.:*(c::ScalarLike, s::PNSeries) = s * c

Base.:/(s::PNSeries, c::ScalarLike) =
    PNSeries(s.var, s.x0, [cc / Num(c) for cc in s.coeffs], s.nmin, s.nmax, s.den)

# ----------------------------------------------------------------------------
# Multiplication (Cauchy product with MMA truncation)
# ----------------------------------------------------------------------------

function Base.:*(a::PNSeries, b::PNSeries)
    a.var == b.var || throw(ArgumentError("cannot multiply series in different variables"))
    isequal(a.x0, b.x0) || throw(ArgumentError("cannot multiply series about different points"))
    den = _lcm(a.den, b.den)
    (ca, amin, amax, _) = _on_den(a, den)
    (cb, bmin, bmax, _) = _on_den(b, den)
    nmin = amin + bmin
    nmax = min(amax + bmin, bmax + amin)   # MMA O[] propagation
    nslots = nmax - nmin
    nslots <= 0 && return PNSeries(a.var, a.x0, Num[], nmin, nmin, den)
    out = fill(Num(0), nslots)
    for (i, x) in enumerate(ca)
        iszero(x) && continue
        ei = amin + (i - 1)
        for (j, y) in enumerate(cb)
            iszero(y) && continue
            e = ei + bmin + (j - 1)
            e >= nmax && continue
            e < nmin && continue
            out[e - nmin + 1] += x * y
        end
    end
    PNSeries(a.var, a.x0, out, nmin, nmax, den)
end

# ----------------------------------------------------------------------------
# Inverse and division
# ----------------------------------------------------------------------------

function Base.inv(s::PNSeries)
    li = findfirst(c -> !iszero(c), s.coeffs)
    li === nothing && throw(ArgumentError("cannot invert a series with all-zero coefficients (pure O[])"))
    den = s.den
    lead = s.nmin + (li - 1)
    # Leading zero slots do not provide extra information.  If
    #   s = x^lead (b0 + ... + O[x]^(nmax-lead)),
    # then inv(s) ends at O[x]^(nmax-2lead), with nmax-lead reciprocal slots.
    span = s.nmax - lead
    b = s.coeffs[li:end]
    nrm = length(b)
    inmin = -lead
    inmax = -lead + span
    m = span
    c = fill(Num(0), max(m, 0))
    m == 0 && return PNSeries(s.var, s.x0, c, inmin, inmax, den)
    b0 = b[1]
    c[1] = 1 / b0
    for k in 2:m
        acc = Num(0)
        for j in 1:(k - 1)
            (k - j + 1) <= nrm || continue
            acc += c[j] * b[k - j + 1]
        end
        c[k] = -acc / b0
    end
    PNSeries(s.var, s.x0, c, inmin, inmax, den)
end

Base.:/(a::PNSeries, b::PNSeries) = a * inv(b)
Base.:/(c::ScalarLike, s::PNSeries) = Num(c) * inv(s)

# ----------------------------------------------------------------------------
# Integer powers
# ----------------------------------------------------------------------------

function Base.:^(s::PNSeries, n::Integer)
    n == 0 && return _one_series(s)
    if n > 0
        r = s
        for _ in 2:n
            r = r * s
        end
        return r
    else
        return inv(s^(-n))
    end
end

# the series "1 + O[var]^span" matching s's truncation footprint (used by s^0)
function _one_series(s::PNSeries)
    nmin = 0
    nmax = s.nmax - s.nmin
    nmax <= 0 && (nmax = s.den)
    nslots = nmax - nmin
    out = fill(Num(0), nslots)
    nslots >= 1 && (out[1] = Num(1))
    PNSeries(s.var, s.x0, out, nmin, nmax, s.den)
end

# ----------------------------------------------------------------------------
# polyToSeries
# ----------------------------------------------------------------------------

"""
    polyToSeries(poly, var; x0=0) -> PNSeries

A polynomial in `var` becomes a series with O[] one past its top degree:
`result = poly + O[var]^(Exponent[poly,var]+1)`. Leading zero coefficients raise `nmin`
(MMA SeriesData auto-detects the leading order). `polyToSeries(0)` returns the scalar `0`.
"""
function polyToSeries(poly, var::Symbol; x0 = 0)
    p = Num(poly)
    iszero(p) && return p
    x = Symbolics.variable(var)
    maxpow = _degree_in(p, x)
    nmax = maxpow + 1
    coeffs = Num[_coeff_at(p, x, k) for k in 0:maxpow]
    nmin = 0
    while !isempty(coeffs) && isequal(coeffs[1], Num(0))
        popfirst!(coeffs)
        nmin += 1
    end
    PNSeries(var, x0, coeffs, nmin, nmax, 1)
end

# ----------------------------------------------------------------------------
# SeriesTake
# ----------------------------------------------------------------------------

"""
    series_take(s, order=1) -> PNSeries

Keep `order` terms from the leading order, i.e. set the O[] order to `nmin + order*den`.
`order=0` gives pure `O[var]^minorder` (all-zero coeffs).
"""
function series_take(s::PNSeries, order::Integer = 1)
    if order == 0
        return PNSeries(s.var, s.x0, Num[], s.nmin, s.nmin, s.den)
    end
    newmax = min(s.nmax, s.nmin + order * s.den)
    nslots = newmax - s.nmin
    PNSeries(s.var, s.x0, s.coeffs[1:min(nslots, length(s.coeffs))], s.nmin, newmax, s.den)
end

# ----------------------------------------------------------------------------
# series_of / SeriesTerms
# ----------------------------------------------------------------------------

"""
    series_of(expr, var; x0=0, order) -> PNSeries

Taylor/Laurent-expand `expr` in `var` about `x0`, keeping terms up to real order `order` (the O[] is
at `var^(order+1)`). Integer powers only (den=1).
"""
function series_of(expr, var::Symbol; x0 = 0, order::Integer)
    iszero(Num(x0)) || throw(ArgumentError("series_of currently supports x0=0 only (Symbolics.taylor about a nonzero point is unstable in this version)"))
    x = Symbolics.variable(var)
    e = Num(expr)
    lead = _leading_order(e, x, x0)
    shifted = lead == 0 ? e : simplify(e * x^(-lead))
    hi = order - lead
    hi < 0 && return PNSeries(var, x0, Num[], lead, lead, 1)
    te = expand(Symbolics.taylor(shifted, x, 0:hi))
    coeffs = Num[_coeff_at(te, x, k) for k in 0:hi]
    PNSeries(var, x0, coeffs, lead, order + 1, 1)
end

# leading (integer) order of e in (x-x0): smallest p with lim (x-x0)^(-p) e finite & nonzero.
function _leading_order(e, x, x0; search = -64:64)
    es = simplify(e)
    iszero(es) && throw(ArgumentError("leading order is undefined for the zero expression"))
    for p in search
        g = es * (x - x0)^(-p)
        val = try
            Symbolics.substitute(g, Dict(x => x0))
        catch
            continue
        end
        if !_is_zero_num(val) && _is_finite_num(val)
            return p
        end
    end
    throw(ArgumentError(
        "could not determine an integer leading order in search range $(first(search)):$(last(search))"))
end

_is_zero_num(v) = (v isa Number && iszero(v)) || (v isa Num && iszero(v))
function _is_finite_num(v)
    vv = v isa Num ? Symbolics.value(v) : v
    vv isa Number ? isfinite(vv) : true
end

"""
    series_terms(expr, var, x0, nterms) -> PNSeries

Expand to `nterms` terms counting from the leading order: `min_order = leading order`, the result
goes to real order `min_order + nterms - 1`. Mirror of MMA `SeriesTerms[expr, {x, x0, nterms}]`.
"""
function series_terms(expr, var::Symbol, x0, nterms::Integer)
    x = Symbolics.variable(var)
    minord = _leading_order(Num(expr), x, x0)
    series_of(expr, var; x0 = x0, order = minord + nterms - 1)
end

# ----------------------------------------------------------------------------
# SeriesCollect
# ----------------------------------------------------------------------------

"""
    series_collect(s, var, func=identity) -> PNSeries

Apply `Collect[#, var, func]` to each coefficient: collect by powers of `var`, applying `func` to each
grouped sub-coefficient (`func=Symbolics.expand` mirrors `SeriesCollect[s, v, Expand]`).
"""
function series_collect(s::PNSeries, var, func = identity)
    PNSeries(s.var, s.x0, [_collect(c, var, func) for c in s.coeffs], s.nmin, s.nmax, s.den)
end

function _collect(expr, var, func)
    x = var isa Symbol ? Symbolics.variable(var) : Num(var)
    e = expand(Num(expr))
    maxp = _degree_in(e, x)
    if maxp <= 0
        return func(e)
    end
    acc = Num(0)
    for k in maxp:-1:1
        ck = _coeff_at(e, x, k)
        if !iszero(ck)
            acc += func(ck) * x^k
        end
    end
    c0 = _coeff_at(e, x, 0)
    acc += func(c0)
    acc
end

# ----------------------------------------------------------------------------
# ChangeSeriesParameter
# ----------------------------------------------------------------------------

"""
    change_series_parameter(s, newvar) -> PNSeries

The four MMA `ChangeSeriesParameter` forms:
  * `newvar::Symbol`            — rename the expansion variable.
  * `(:rescale, fac, sym)`      — `fac*sym`: keep symbol, scale coeff i by `fac^exp_i`.
  * `(:power, sym, pow)`        — `sym^pow`: insert `pow_num-1` zero coeffs, scale orders.
  * `(:facpower, fac, sym, pow)`— `fac*sym^pow`: both, scaling THEN zero-insertion.

`pow` may be a rational `p//q`.
"""
change_series_parameter(s::PNSeries, newvar::Symbol) =
    PNSeries(newvar, s.x0, copy(s.coeffs), s.nmin, s.nmax, s.den)

function change_series_parameter(s::PNSeries, ::Val{:rescale}, fac, sym::Symbol)
    f = Num(fac)
    nc = Vector{Num}(undef, length(s.coeffs))
    for (i, c) in enumerate(s.coeffs)
        exp_i = (s.nmin + (i - 1)) // s.den
        nc[i] = c * f^exp_i
    end
    PNSeries(sym, s.x0, nc, s.nmin, s.nmax, s.den)
end

function change_series_parameter(s::PNSeries, ::Val{:power}, sym::Symbol, pow::Union{Integer,Rational})
    pr = Rational(pow)
    powNum = numerator(pr); powDen = denominator(pr)
    nc = _interleave_zeros(s.coeffs, powNum)
    PNSeries(sym, s.x0, nc, s.nmin * powNum, s.nmax * powNum, s.den * powDen)
end

function change_series_parameter(s::PNSeries, ::Val{:facpower}, fac, sym::Symbol, pow::Union{Integer,Rational})
    f = Num(fac)
    pr = Rational(pow)
    powNum = numerator(pr); powDen = denominator(pr)
    scaled = Vector{Num}(undef, length(s.coeffs))
    for (i, c) in enumerate(s.coeffs)
        exp_i = (s.nmin + (i - 1)) // s.den
        scaled[i] = c * f^exp_i
    end
    nc = _interleave_zeros(scaled, powNum)
    PNSeries(sym, s.x0, nc, s.nmin * powNum, s.nmax * powNum, s.den * powDen)
end

# place each coeff then (powNum-1) zeros after it
function _interleave_zeros(coeffs::Vector{Num}, powNum::Integer)
    powNum == 1 && return copy(coeffs)
    out = Num[]
    for c in coeffs
        push!(out, c)
        for _ in 1:(powNum - 1)
            push!(out, Num(0))
        end
    end
    out
end

# ----------------------------------------------------------------------------
# PowerCounting
# ----------------------------------------------------------------------------

"""
    power_counting(s, newvar) -> PNSeries

Replace the expansion variable by `newvar` while keeping the ORIGINAL `sym^exp_i` as a constant factor
folded into each coefficient. Mirror of MMA `PowerCounting`.
"""
function power_counting(s::PNSeries, newvar::Symbol)
    x = svar(s)
    nc = Vector{Num}(undef, length(s.coeffs))
    for (i, c) in enumerate(s.coeffs)
        exp_i = (s.nmin + (i - 1)) // s.den
        nc[i] = c * x^exp_i
    end
    PNSeries(newvar, s.x0, nc, s.nmin, s.nmax, s.den)
end

# ----------------------------------------------------------------------------
# Polynomial coefficient helpers (degree, coeff-at) over Symbolics
# ----------------------------------------------------------------------------

# Highest integer power of x appearing in expanded e (>=0). Returns 0 if x absent.
# `Symbolics.degree` is the robust polynomial-degree primitive: under Symbolics 7 a power exponent is
# a symbolic constant literal (not a Julia `Integer`), so a hand-rolled `args[2] isa Integer` tree
# walk silently misses every `x^k` and collapses the degree — `degree` unwraps it correctly.
function _degree_in(e, x)
    ex = expand(Num(e))
    iszero(ex) && return 0
    d = try
        Int(Symbolics.degree(ex, x))
    catch
        0
    end
    max(d, 0)
end

# coefficient of x^k in e (k>=0). The x^0 term is e|_{x=0}.
function _coeff_at(e, x, k)
    en = Num(e)
    k == 0 && return Symbolics.substitute(en, Dict(x => 0))
    Num(Symbolics.coeff(en, x^k))
end

# ----------------------------------------------------------------------------
# ExpandLog and ExpandGamma
# ----------------------------------------------------------------------------

"""
    expand_log(expr) -> Num

Rewrite every `log` of a product/power into a sum/multiple of logs (PowerExpand semantics, assuming
positivity): `log(a*b)->log a+log b`, `log(a^p)->p log a`.
"""
function expand_log(expr)
    e = Num(expr)
    _maptree(e) do t
        if Symbolics.iscall(t) && Symbolics.operation(t) === log
            return _powerexpand_log(Symbolics.arguments(t)[1])
        end
        return nothing
    end
end

function _powerexpand_log(arg)
    a = value(arg)
    if a isa Number
        return Num(log(a))
    end
    if Symbolics.iscall(a)
        op = Symbolics.operation(a)
        args = Symbolics.arguments(a)
        if op === (*)
            return sum(_powerexpand_log(x) for x in args)
        elseif op === (^)
            base, p = args[1], args[2]
            return Num(p) * _powerexpand_log(base)
        elseif op === (/)
            return _powerexpand_log(args[1]) - _powerexpand_log(args[2])
        end
    end
    return Num(log(a))
end

"""
    expand_gamma(expr; shift=0) -> Num

Rewrite `Gamma[Integer n + x]` via the rising/falling factorial so the residual is `Gamma[x+shift]`.
For `shift=0`, `Gamma[n+x] -> (prod_{i=0}^{n-1}(x+i)) Gamma[x]` for integer n (negative n giving
inverse factors).
"""
function expand_gamma(expr; shift::Integer = 0)
    e = Num(expr)
    _maptree(e) do t
        if Symbolics.iscall(t) && Symbolics.operation(t) === SpecialFunctions.gamma
            return _expandgamma_one(Symbolics.arguments(t)[1], shift)
        end
        return nothing
    end
end

function _expandgamma_one(arg, j)
    a = value(arg)
    n, x = _split_integer_plus(a)
    num = Num(1)
    for i in 0:(n - j - 1)
        num *= (Num(x) + j + i)
    end
    den = Num(1)
    for i in (n - j):(-1)
        den *= (Num(x) + j + i)
    end
    num / den * Num(SpecialFunctions.gamma(Num(x) + j))
end

# Recover a concrete Julia integer from a (possibly symbolic-wrapped) constant, else `nothing`.
# Under Symbolics 7 the integer summands of `n + x` are `BasicSymbolic` constants, not Julia
# `Integer`s, so a bare `isa Integer` test never fires; `Symbolics.value` unwraps the literal.
function _const_int(t)
    v = t isa Number ? t : Symbolics.value(t)
    v isa Integer && return Int(v)
    try
        return Int(v)
    catch
        return nothing
    end
end

# arg = Integer n + rest  ->  (n, rest);  else (0, arg).
function _split_integer_plus(a)
    av = value(a)
    if av isa Number
        return (0, av)   # a pure number arg is taken wholesale as the residual x (no split)
    end
    if Symbolics.iscall(av) && Symbolics.operation(av) === (+)
        args = Symbolics.arguments(av)
        intpart = 0
        rest = Num(0)
        nfound = 0
        for t in args
            ti = _const_int(t)
            if ti !== nothing
                intpart += ti
                nfound += 1
            else
                rest += Num(t)
            end
        end
        if nfound > 0
            return (intpart, rest)
        end
    end
    return (0, Num(av))
end

# ----------------------------------------------------------------------------
# Generic post-order tree rewrite over Symbolics expressions.
# ----------------------------------------------------------------------------
_maptree(f, e::Num) = Num(_maptree_val(f, value(e)))

function _maptree_val(f, t)
    if !Symbolics.iscall(t)
        r = f(t)
        return r === nothing ? t : value(r)
    end
    op = Symbolics.operation(t)
    args = Symbolics.arguments(t)
    newargs = [_maptree_val(f, a) for a in args]
    rebuilt = op(map(Num, newargs)...)
    r = f(value(rebuilt))
    return r === nothing ? value(rebuilt) : value(r)
end

end # module Series

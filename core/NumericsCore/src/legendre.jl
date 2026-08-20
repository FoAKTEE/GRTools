# Normalized associated-Legendre + spin-weighted spherical-harmonic angular SEED.
#
# The canonical {T}-generic angular primitive that Harmonics / MetricReconstruction / DiscModes /
# SelfForce build on. Consolidates (AHA) the duplicate Condon-Shortley upward-recurrence Legendre
# kernels that lived in H1Lorenz/src/math_funcs.jl and DiscModes/src/metric_recon.jl into ONE
# `_legendre_Plm` primitive, plus the SpinWeightedSpheroidalHarmonics half-angle term-list machinery
# for the spin-weighted seed. Everything is type-generic so the source coefficients inherit BigFloat
# precision; the spin-weighted seed additionally exact-integers the half-angle binomials in BigFloat at
# high l (Int64 overflow / Float64 cancellation regime) and at any BigFloat element type.

"""
    factorial_int(::Type{T}, n::Integer) -> T

Exact integer factorial `n!` accumulated in `T` (overflow-safe at high `l` where an `Int64`
product would wrap); `factorial_int(T, 0) = factorial_int(T, 1) = one(T)`. The single home of
the integer-product factorial (BigSpecial.factorial_big and the Lorenz-field static series forward
to it).
"""
factorial_int(::Type{T}, n::Integer) where {T} = n <= 1 ? one(T) : prod(T(2):T(n))

"""
    _legendre_Plm(l, m, x) -> T

Associated Legendre `P_l^m(x)` for `m≥0`, WITH the Condon-Shortley `(-1)^m` phase, via the stable
upward recurrence. Type-generic in `x::T`; `m>l ⇒ 0`. This is the single shared kernel (the former
`H1Lorenz._legendre_Plm` and `DiscModes.assoc_legendre` collapsed into one).
"""
function _legendre_Plm(l::Integer, m::Integer, x::T) where {T}
    m = abs(m)
    m > l && return zero(T)
    pmm = one(T)
    if m > 0
        s = sqrt(1 - x^2); ff = one(T)
        for _ in 1:m
            pmm *= -ff * s; ff += 2
        end
    end
    l == m && return pmm
    pmmp1 = x * (2m + 1) * pmm
    l == m + 1 && return pmmp1
    pll = zero(T)
    for ll in m+2:l
        pll = ((2ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm = pmmp1; pmmp1 = pll
    end
    pll
end

"`assoc_legendre(l, m, x)` — DiscModes name for the consolidated `_legendre_Plm` primitive."
assoc_legendre(l::Integer, m::Integer, x::T) where {T} = _legendre_Plm(l, m, x)

"Normalized GSL `sphPlm` prefactor `N_l^|m| = sqrt((2l+1)(l-|m|)!/(4π(l+|m|)!))`, computed as the
reciprocal falling factorial in `T` so it is overflow-free for `l+|m| ≥ 21` (Int64 factorial caps at 20!)."
function _norm_factor(::Type{T}, l::Integer, am::Integer) where {T}
    rat = one(T)
    for k in (l - am + 1):(l + am)
        rat /= k
    end
    sqrt((2l + 1) / (4 * T(π)) * rat)
end

"""
    assoc_legendre_P(l, m, θ; normalized=false) -> T

Associated Legendre at the colatitude `θ`. `normalized=false` is the GSL `Plm` convention
`P_l^{|m|}(cosθ)` (Condon-Shortley phase, `(-1)^|m|` for `m<0`); `normalized=true` is the GSL `sphPlm`
convention `N_l^|m| P_l^|m|(cosθ)` (same `m<0` sign). Includes the `θ=π/2` parity shortcut
(`l+m` odd ⇒ exactly 0). Consolidates `H1Lorenz.associated_legendre_P` +
`normalized_associated_legendre_P`.
"""
function assoc_legendre_P(l::Integer, m::Integer, θ::T; normalized::Bool=false) where {T}
    abs(m) > l && return zero(T)
    if θ == T(π) / 2 && isodd(l + m)
        return zero(T)
    end
    am = abs(m)
    base = _legendre_Plm(l, am, cos(θ))
    val = normalized ? _norm_factor(T, l, am) * base : base
    m >= 0 ? val : (iseven(am) ? val : -val)
end

"""
    legendre_P_thetaderiv(l, m, θ; normalized=false) -> T

The `(l-1)` coupling closed form `-(-l cosθ P_l^m + (l+m) P_{l-1}^m)/(1-cos²θ)` (math_funcs.c
L198-204), preserved byte-for-byte from `H1Lorenz.associated_legendre_P_theta_deriv` (+ its
`normalized=true` `sphPlm`-prefactor variant). NOTE this legacy quantity equals `-dP_l^m/d(cosθ)`,
i.e. `(1/sinθ)·dP_l^m/dθ` — the genuine colatitude derivative is `sinθ` times this value (the legacy
name is a misnomer kept for downstream-call compatibility). Consolidates both legacy θ-deriv forks.
"""
function legendre_P_thetaderiv(l::Integer, m::Integer, θ::T; normalized::Bool=false) where {T}
    x = cos(θ)
    d = -(-l * x * assoc_legendre_P(l, m, θ) + (l + m) * assoc_legendre_P(l - 1, m, θ)) / (1 - x^2)
    normalized ? _norm_factor(T, l, abs(m)) * d : d
end

"""
    Ylm_real(l, m, θ) -> T

Real scalar spherical-harmonic angular part `Y_lm(θ, φ=0) = N_l^|m| P_l^|m|(cosθ)` (the `e^{imφ}`
phase is applied separately by the caller). Uses `|m|` with no `m<0` sign flip (DiscModes convention);
the normalization is the overflow-free reciprocal falling factorial.
"""
Ylm_real(l::Integer, m::Integer, θ::T) where {T} =
    _norm_factor(T, l, abs(m)) * _legendre_Plm(l, abs(m), cos(θ))

"""
    Ylm_eq(l, m, θ; h=T(1e-5)) -> (Y, ∂_θY, ∂²_θY)

`Y_lm(θ,0)` and its first two θ-derivatives via central finite differences (type-generic step `h`).
"""
function Ylm_eq(l::Integer, m::Integer, θ::T; h::T=T(1e-5)) where {T}
    Y(t) = Ylm_real(l, m, t)
    (Y(θ), (Y(θ + h) - Y(θ - h)) / (2h), (Y(θ + h) - 2Y(θ) + Y(θ - h)) / h^2)
end

"""
    pl_deriv(T, l, m, sign) -> T

`m`-th `x`-derivative of the Legendre polynomial `P_l(x)` at `x = ±1`:
`Derivative[0,m][P_l](±1) = (1/(2^m m!)) ∏_{k=-m+1}^{m}(l+k)`, with the `x=-1` branch carrying the
extra `(-1)^(l+m)`. (Punctures `_pl_deriv`.)
"""
function pl_deriv(::Type{T}, l::Integer, m::Integer, sign::Integer) where {T}
    prd = one(T)
    for k in (-m + 1):m
        prd *= T(l + k)
    end
    val = prd / (T(2)^m * T(factorial(m)))
    sign == 1 ? val : T((-1)^(l + m)) * val
end

# --- spin-weighted spherical-harmonic seed: half-angle (coeff, a, b) term-list machinery (SWSH §2.1) ---

"Above this l the Float64 half-angle binomial sum loses too many digits to cancellation; switch to the
exact-integer BigFloat path. (Also used whenever the element type is itself BigFloat.)"
const _LHIGH_SWSHY = 20

"half-angle power with the Mathematica `0^0 -> 1` guard."
@inline _pow0(base, e::Integer) = e == 0 ? one(base) : base^e

"prefactor + Float64 half-angle term list `(coeff, a, b)` of `{}_sY_{lm}(θ,0)` (fast path, `l<_LHIGH_SWSHY`)."
function _y_terms(s::Integer, l::Integer, m::Integer)
    num = factorial(big(l + m)) * factorial(big(l - m)) * (2l + 1)
    den = factorial(big(l + s)) * factorial(big(l - s))
    pref = (isodd(m) ? -1.0 : 1.0) * sqrt(Float64(num // den) / (4 * pi))
    terms = Tuple{Float64,Int,Int}[]
    for r in max(m - s, 0):min(l - s, l + m)
        c = binomial(l - s, r) * binomial(l + s, r + s - m) * (isodd(l - r - s) ? -1 : 1)
        push!(terms, (float(c), 2l - 2r - s + m, 2r + s - m))
    end
    return pref, terms
end

"d/dθ of the Float64 term list: `(c,a,b) -> (c a/2, a-1, b+1) + (-c b/2, a+1, b-1)` (`du=v/2, dv=-u/2`)."
function _dtheta_terms(terms::Vector{Tuple{Float64,Int,Int}})
    out = Tuple{Float64,Int,Int}[]
    for (c, a, b) in terms
        a != 0 && push!(out, (c * a / 2, a - 1, b + 1))
        b != 0 && push!(out, (-c * b / 2, a + 1, b - 1))
    end
    return out
end

"evaluate `pref * Σ coeff sin(θ/2)^a cos(θ/2)^b` (Float64 fast path)."
function _eval_terms(pref::Real, terms::Vector{Tuple{Float64,Int,Int}}, θ::Real)
    u = sin(θ / 2); v = cos(θ / 2)
    return pref * sum(c * _pow0(u, a) * _pow0(v, b) for (c, a, b) in terms; init=0.0)
end

"""
    _y_eval_big(T, s, l, m, θ, d1) -> T

`d^{d1}/dθ^{d1} {}_sY_{lm}(θ,0)` (real) with exact-integer half-angle binomials (`BigInt` →
`Rational{BigInt}` under differentiation), accumulated in `BigFloat` and returned as `T`. No Int64
overflow nor Float64 cancellation at high `l`; when `T=BigFloat` the full BigFloat precision is kept.
"""
function _y_eval_big(::Type{T}, s::Integer, l::Integer, m::Integer, θ::Real, d1::Integer) where {T}
    num = factorial(big(l + m)) * factorial(big(l - m)) * (2l + 1)
    den = factorial(big(l + s)) * factorial(big(l - s))
    pref = (isodd(m) ? -one(BigFloat) : one(BigFloat)) * sqrt(BigFloat(num // den) / (4 * BigFloat(pi)))
    terms = Tuple{Rational{BigInt},Int,Int}[]
    for r in max(m - s, 0):min(l - s, l + m)
        c = binomial(big(l - s), big(r)) * binomial(big(l + s), big(r + s - m)) * (isodd(l - r - s) ? -1 : 1)
        push!(terms, (Rational{BigInt}(c), 2l - 2r - s + m, 2r + s - m))
    end
    for _ in 1:d1
        out = Tuple{Rational{BigInt},Int,Int}[]
        for (c, a, b) in terms
            a != 0 && push!(out, (c * a // 2, a - 1, b + 1))
            b != 0 && push!(out, (-c * b // 2, a + 1, b - 1))
        end
        terms = out
    end
    u = sin(BigFloat(θ) / 2); v = cos(BigFloat(θ) / 2)
    acc = sum(BigFloat(c) * _pow0(u, a) * _pow0(v, b) for (c, a, b) in terms; init=zero(BigFloat))
    return T(pref * acc)
end

"""
    _y_dtheta(s, l, m, θ, d1) -> T

`d^{d1}/dθ^{d1} {}_sY_{lm}(θ,0)` (real, type-generic). Uses the fast Float64 term list for
`l < _LHIGH_SWSHY` with a Float64 element type, else the exact-integer BigFloat path.
"""
function _y_dtheta(s::Integer, l::Integer, m::Integer, θ::T, d1::Integer) where {T}
    if l >= _LHIGH_SWSHY || T <: BigFloat
        return _y_eval_big(T, s, l, m, θ, d1)
    end
    pref, terms = _y_terms(s, l, m)
    for _ in 1:d1
        terms = _dtheta_terms(terms)
    end
    return _eval_terms(pref, terms, θ)
end

"""
    spin_weighted_Y(s, l, m, θ, φ) -> Complex

Spin-weighted spherical harmonic `{}_sY_{lm}(θ,φ)` (the γ=0 seed). Half-angle binomial sum; the fast
Float64 term list is used for `l < _LHIGH_SWSHY` with a Float64 colatitude, otherwise the exact-integer
BigFloat path (and full precision when the element type is BigFloat). Requires `l ≥ |s|` and `|m| ≤ l`.
"""
function spin_weighted_Y(s::Integer, l::Integer, m::Integer, θ::T, φ::Real) where {T}
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    phase = cis(m * float(φ))
    if l >= _LHIGH_SWSHY || T <: BigFloat
        return complex(_y_eval_big(T, s, l, m, θ, 0)) * phase
    end
    pref, terms = _y_terms(s, l, m)
    return complex(_eval_terms(pref, terms, θ)) * phase
end

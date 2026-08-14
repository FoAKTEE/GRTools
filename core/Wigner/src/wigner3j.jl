# Wigner 3j and 6j symbols — hand-rolled Racah closed forms, BigInt-exact.
#
# The canonical sink replacing the WignerSymbols PACKAGE (design §5.3 / §7.4: package
# ELIMINATED from every runtime path). Type-generic {T}: the combinatorics run in
# BigInt/Rational{BigInt} arithmetic and narrow to T only at the boundary, so BigFloat
# precision (set via setprecision) is inherited by the √-prefactors.

# Exact factorial in BigInt (loop bounds guarantee a nonnegative argument at every call site).
_fbig(n::Integer) = factorial(big(n))

# Triangle coefficient Δ(a,b,c) = (a+b-c)!(a-b+c)!(-a+b+c)! / (a+b+c+1)!, as a Rational{BigInt}.
_delta(a::Integer, b::Integer, c::Integer) =
    _fbig(a + b - c) * _fbig(a - b + c) * _fbig(-a + b + c) // _fbig(a + b + c + 1)

# Triangle inequality |a-b| ≤ c ≤ a+b (integer arguments).
_triangle(a::Integer, b::Integer, c::Integer) = (c >= abs(a - b)) && (c <= a + b)

# Narrow a Rational{BigInt} square-root to T at the active precision: do the division and
# sqrt in BigFloat (honours setprecision(BigFloat,…)) then convert, so T=BigFloat is exact
# to the working precision and T=Float64 is full double accuracy.
_sqrt_rat(::Type{T}, r::Rational{BigInt}) where {T} =
    T(sqrt(BigFloat(numerator(r)) / BigFloat(denominator(r))))

"""
    wigner3j(::Type{T}, j1, j2, j3, m1, m2, m3) -> T

Wigner 3j symbol `(j1 j2 j3; m1 m2 m3)` via the Racah single-sum closed form. Returns
`zero(T)` on any selection-rule violation (triangle inequality, `m1+m2+m3 ≠ 0`, or
`|mᵢ| > jᵢ`) to match GSL / the `coupling` caller. Integer arguments.
"""
function wigner3j(::Type{T}, j1::Integer, j2::Integer, j3::Integer,
                  m1::Integer, m2::Integer, m3::Integer) where {T}
    (m1 + m2 + m3 != 0) && return zero(T)
    (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3) && return zero(T)
    _triangle(j1, j2, j3) || return zero(T)

    pref2 = _delta(j1, j2, j3) *
            (_fbig(j1 + m1) * _fbig(j1 - m1) * _fbig(j2 + m2) * _fbig(j2 - m2) *
             _fbig(j3 + m3) * _fbig(j3 - m3))

    kmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    kmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    ssum = big(0) // 1
    for k in kmin:kmax
        ssum += (iseven(k) ? 1 : -1) //
                (_fbig(k) * _fbig(j1 + j2 - j3 - k) * _fbig(j1 - m1 - k) *
                 _fbig(j2 + m2 - k) * _fbig(j3 - j2 + m1 + k) * _fbig(j3 - j1 - m2 + k))
    end

    sgn = iseven(j1 - j2 - m3) ? one(T) : -one(T)
    return sgn * _sqrt_rat(T, pref2) * T(ssum)
end

"""
    wigner6j(::Type{T}, j1, j2, j3, j4, j5, j6) -> T

Wigner 6j symbol `{j1 j2 j3; j4 j5 j6}` via the Racah closed form: the four triangle Δ
coefficients times the single Racah W-sum. Returns `zero(T)` unless all four triads
`(j1 j2 j3)`, `(j1 j5 j6)`, `(j4 j2 j6)`, `(j4 j5 j3)` satisfy the triangle inequality.
Integer arguments.
"""
function wigner6j(::Type{T}, j1::Integer, j2::Integer, j3::Integer,
                  j4::Integer, j5::Integer, j6::Integer) where {T}
    (_triangle(j1, j2, j3) && _triangle(j1, j5, j6) &&
     _triangle(j4, j2, j6) && _triangle(j4, j5, j3)) || return zero(T)

    Δprod = _delta(j1, j2, j3) * _delta(j1, j5, j6) *
            _delta(j4, j2, j6) * _delta(j4, j5, j3)

    tmin = max(j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3)
    tmax = min(j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4)
    ssum = big(0) // 1
    for t in tmin:tmax
        ssum += (iseven(t) ? 1 : -1) * _fbig(t + 1) //
                (_fbig(t - j1 - j2 - j3) * _fbig(t - j1 - j5 - j6) *
                 _fbig(t - j4 - j2 - j6) * _fbig(t - j4 - j5 - j3) *
                 _fbig(j1 + j2 + j4 + j5 - t) * _fbig(j2 + j3 + j5 + j6 - t) *
                 _fbig(j3 + j1 + j6 + j4 - t))
    end

    return _sqrt_rat(T, Δprod) * T(ssum)
end

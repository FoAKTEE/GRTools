# Terminating-2F1 polar-angle α-integration: the per-(p,l,m) angular integrals the first-order
# puncture mode integrands contract against. Port of the "Integration over α" subsection of
# Punctures/alpha_integration.jl (Wardell et al. arXiv:2508.00087), REDESIGNED:
#
#   • the local generalised-Pochhammer fork `_poch` is DROPPED onto the one canonical
#     `BigSpecial.pochhammer` (AHA — the SAME rising/falling-product decision Teukolsky + Punctures
#     both made; one leaf primitive).
#   • the local terminating Gauss `_f21_term` is DROPPED onto `BigSpecial.f21_terminating` (the same
#     a=-l-non-positive-integer finite sum; re-exported here as `_f21_term` for the internal callers
#     and the collapse cross-check).
#
# Everything reduces to closed forms: Legendre values at ±1, a terminating Gauss 2F1 (a = -l a
# non-positive integer), and the regulator polynomial obtained from an incomplete Beta with integer
# parameters (norder=4, mpmax=10). All evaluators are numeric in a real δ and inherit the float
# precision of δ (pass BigFloat for extended precision). The upstream `reg` discards its second
# argument and always uses mpmax, so αintmp2/αintmp4 coincide with the even-m αint; that identity
# is reproduced rather than re-derived.

import BigSpecial

const NORDER = 4
const MPMAX = 10

# Generalised Pochhammer (a)_n — collapsed onto BigSpecial.pochhammer (handles negative n via the
# reciprocal-falling product, identical to the legacy `_poch`).
const _poch = BigSpecial.pochhammer

# Terminating Gauss 2F1(-l, 1+l; c; z) — collapsed onto BigSpecial.f21_terminating.
const _f21_term = BigSpecial.f21_terminating

# Binomial C(n,k) for small non-negative integers, returned in type T (exact).
_binom(::Type{T}, n::Integer, k::Integer) where {T} = T(binomial(n, k))

# Coefficient C_{m,n}(δ) common to dint/rest:
#   (-1)^(m/2+n) Σ_{i=0}^{n} C(m/2,i) C(m/2,n-i) |δ|^(m-2i) (2+δ²)^(m/2-n+i)  (m even).
function _dint_coeff(m::Integer, n::Integer, δ::T) where {T}
    h = m ÷ 2
    d2 = δ * δ
    s = zero(T)
    for i in 0:n
        s += _binom(T, h, i) * _binom(T, h, n - i) *
             abs(δ)^(m - 2i) * (T(2) + d2)^(h - n + i)
    end
    return T((-1)^(h + n)) * s
end

# m-th x-derivative of the Legendre polynomial P_l(x) evaluated at x = ±1.
function _pl_deriv(::Type{T}, l::Integer, m::Integer, sign::Integer) where {T}
    prod = one(T)
    for k in (-m + 1):m
        prod *= T(l + k)
    end
    val = prod / (T(2)^m * T(factorial(m)))
    return sign == 1 ? val : T((-1)^(l + m)) * val
end

# r-th x-derivative of R(x,δ)^(p/2) at x = ±1 (R has dR/dx = -1).
function _dR(x::Integer, p::Integer, r::Integer, δ::T) where {T}
    pre = T((-1)^r) * _poch(one(T) + T(p) / 2 - r, r)
    if x == 1
        return pre * abs(δ)^(p - 2r)
    else
        return pre * (T(2) + δ * δ)^(T(p) / 2 - r)
    end
end

# dint[x,p,m,r,δ] = Σ_{n=0}^{m} C_{m,n}(δ) dR[x,2n+p,r,δ]  (x = ±1, numeric).
function _dint_pm(x::Integer, p::Integer, m::Integer, r::Integer, δ::T) where {T}
    s = zero(T)
    for n in 0:m
        s += _dint_coeff(m, n, δ) * _dR(x, 2n + p, r, δ)
    end
    return s
end

# Boundary term Σ_{r=0}^{m-1} (-1)^r [ dint(1,…) Pl^(m-1-r)(l,1) − dint(-1,…) Pl^(m-1-r)(l,-1) ].
function _bound(p::Integer, l::Integer, m::Integer, δ::T) where {T}
    s = zero(T)
    for r in 0:(m - 1)
        sgn = T((-1)^r)
        s += sgn * (_dint_pm(1, p, m, r, δ) * _pl_deriv(T, l, m - 1 - r, 1) -
                    _dint_pm(-1, p, m, r, δ) * _pl_deriv(T, l, m - 1 - r, -1))
    end
    return s
end

"""
    exactint(p, l, δ)

Closed form of the radial-exact angular integral, the `R^n → exactint` building block of `rest`.
`p` is an (odd) integer; `δ` is real. Inherits the float precision of `δ`.
"""
function exactint(p::Integer, l::Integer, δ::T) where {T<:Real}
    d2 = δ * δ
    z = -d2 / 2
    half = one(T) / 2
    pre1 = T((-1)^((p + 1) ÷ 2)) * _poch(half, (p + 1) ÷ 2)^2 /
           _poch(T(l) - T(p) / 2, 2 + p)
    t1 = pre1 * (T(2) + d2)^(T(p) / 2 + 1) * _f21_term(l, -T(p) / 2, z)
    t2 = -T(2) / (2 + p) * δ^(1 + p) * abs(δ) * _f21_term(l, 2 + T(p) / 2, z)
    return t1 + t2
end
exactint(p::Integer, l::Integer, δ::Integer) = exactint(p, l, float(δ))

# rest[p,l,m,δ] = Σ_{n=0}^{m} C_{m,n}(δ) (-1)^m (1+(2n+p)/2-m)_m exactint(2n+p-2m,l,δ).
function _rest(p::Integer, l::Integer, m::Integer, δ::T) where {T}
    s = zero(T)
    for n in 0:m
        q = 2n + p
        s += _dint_coeff(m, n, δ) * T((-1)^m) *
             _poch(one(T) + T(q) / 2 - m, m) * exactint(q - 2m, l, δ)
    end
    return s
end

"""
    modeint(p, l, m, δ)

Angular mode integral ∫ ρ^p P_l^m(cos α) … dα. For `m ≥ 0` it is `bound + rest`; for `m < 0` it
uses the reflection `modeint[p,l,−m] = (−1)^m (l+m)!/(l−m)! modeint[p,l,m]`.
"""
function modeint(p::Integer, l::Integer, m::Integer, δ::T) where {T<:Real}
    if m ≥ 0
        return _bound(p, l, m, δ) + _rest(p, l, m, δ)
    else
        ratio = one(T)
        for j in (l + m + 1):(l - m)
            ratio *= T(j)
        end
        return T((-1)^m) / ratio * modeint(p, l, -m, δ)
    end
end
modeint(p::Integer, l::Integer, m::Integer, δ::Integer) = modeint(p, l, m, float(δ))

# Coefficients of the regulator polynomial reg(y) in powers of y, for the fixed (norder, mpmax).
_reg_ycoeffs(::Type{T}) where {T} = _reg_ycoeffs_n(T, NORDER)

function _reg_ycoeffs_n(::Type{T}, n::Integer) where {T}
    n == 0 && return T[one(T)]     # reg(0,·,y) = 1 identically
    a = n ÷ 2
    b = MPMAX ÷ 2
    pref = T(n) / 2 * _binom(T, (MPMAX + n - 2) ÷ 2, n ÷ 2)
    deg = a + (b - 1)
    c = zeros(T, deg + 1)          # c[q+1] = coeff of y^q
    c[1] = one(T)
    for j in 0:(b - 1)
        q = a + j
        c[q + 1] -= pref * _binom(T, b - 1, j) * T((-1)^j) / T(q)
    end
    return c
end

# modeintreg: substitute y = (R - δ²)/2 into reg(y), expand in R, and map R^k → modeint[p+2k,l,m,δ].
function _modeintreg(p::Integer, l::Integer, m::Integer, δ::T) where {T}
    yc = _reg_ycoeffs(T)
    d2 = δ * δ
    s = zero(T)
    for q in 0:(length(yc) - 1)
        cq = yc[q + 1]
        cq == 0 && continue
        for k in 0:q
            coeff = cq * _binom(T, q, k) * (-d2)^(q - k) / T(2)^q
            s += coeff * modeint(p + 2k, l, m, δ)
        end
    end
    return s
end

"""
    alphaint(p, l, m, δ)

Angular integral `αint[p,l,m,δ]` (scalar sector). Even `m` uses the regulated mode integral
`modeintreg`; odd `m` uses the lowering recurrence
`αint[p,l,m] = 1/(2l+1)((l−m+1)(l−m+2)αint[p,l+1,m−1] − (l+m−1)(l+m)αint[p,l−1,m−1])`.
"""
function alphaint(p::Integer, l::Integer, m::Integer, δ::T) where {T<:Real}
    if iseven(m)
        return _modeintreg(p, l, m, δ)
    else
        return one(T) / (2l + 1) *
               ((l - m + 1) * (l - m + 2) * alphaint(p, l + 1, m - 1, δ) -
                (l + m - 1) * (l + m) * alphaint(p, l - 1, m - 1, δ))
    end
end
alphaint(p::Integer, l::Integer, m::Integer, δ::Integer) = alphaint(p, l, m, float(δ))

# αintmp2/αintmp4 coincide with the even-m αint (upstream reg ignores its arg).
_alphaint_mp(p::Integer, l::Integer, m::Integer, δ::T) where {T} = _modeintreg(p, l, m, δ)

"""
    alphaintv1(p, l, m, δ)

Vector-sector integral of type 1 (odd m):
`½((l−m+1)(l+m) αintmp2[p,l,m−1] − αint[p,l,m+1])`.
"""
alphaintv1(p::Integer, l::Integer, m::Integer, δ::Real) =
    let δf = float(δ)
        (1 // 2) * ((l - m + 1) * (l + m) * _alphaint_mp(p, l, m - 1, δf) -
                    alphaint(p, l, m + 1, δf))
    end

"""
    alphaintv2(p, l, m, δ)

Vector-sector integral of type 2 (odd m):
`−½(αint[p,l−1,m+1] + (l+m−1)(l+m) αintmp2[p,l−1,m−1])`.
"""
alphaintv2(p::Integer, l::Integer, m::Integer, δ::Real) =
    let δf = float(δ)
        -(1 // 2) * (alphaint(p, l - 1, m + 1, δf) +
                     (l + m - 1) * (l + m) * _alphaint_mp(p, l - 1, m - 1, δf))
    end

"""
    alphaintv3(p, l, m, δ)

Vector-sector integral of type 3:
`1/(2l+1)((l+1)(l+m) αint[p,l−1,m] − l(l−m+1) αint[p,l+1,m])`.
"""
alphaintv3(p::Integer, l::Integer, m::Integer, δ::Real) =
    let δf = float(δ)
        (1 // (2l + 1)) * ((l + 1) * (l + m) * alphaint(p, l - 1, m, δf) -
                           l * (l - m + 1) * alphaint(p, l + 1, m, δf))
    end

"""
    alphaintv4(p, l, m, δ)

Vector-sector integral of type 4: `m αint[p,l,m]`.
"""
alphaintv4(p::Integer, l::Integer, m::Integer, δ::Real) = m * alphaint(p, l, m, float(δ))

"""
    alphaintt1(p, l, m, δ)

Tensor-sector integral of type 1 (even m):
`½(1+l−m)(2+l−m)(−1+l+m)(l+m) αintmp4[p,l,m−2] + m² αintmp2[p,l,m] + ½ αint[p,l,m+2]`.
"""
alphaintt1(p::Integer, l::Integer, m::Integer, δ::Real) =
    let δf = float(δ)
        (1 // 2) * (1 + l - m) * (2 + l - m) * (-1 + l + m) * (l + m) *
        _alphaint_mp(p, l, m - 2, δf) +
        m^2 * _alphaint_mp(p, l, m, δf) +
        (1 // 2) * alphaint(p, l, m + 2, δf)
    end

"""
    alphaintt2(p, l, m, δ)

Tensor-sector integral of type 2:
`¼(−l⁴+2l³(−3+m)+l²(−11+6m)+m(−6+11m−6m²+m³)−2l(3−3m²+m³)) αintmp4[p,l+1,m−2]
 − ½(1+l−m)m αintmp2[p,l+1,m] + ¼ αint[p,l+1,m+2]`.
"""
alphaintt2(p::Integer, l::Integer, m::Integer, δ::Real) =
    let δf = float(δ)
        coeff = (1 // 4) * (-l^4 + 2 * l^3 * (-3 + m) + l^2 * (-11 + 6m) +
                            m * (-6 + 11m - 6 * m^2 + m^3) -
                            2 * l * (3 - 3 * m^2 + m^3))
        coeff * _alphaint_mp(p, l + 1, m - 2, δf) -
        (1 // 2) * (1 + l - m) * m * _alphaint_mp(p, l + 1, m, δf) +
        (1 // 4) * alphaint(p, l + 1, m + 2, δf)
    end

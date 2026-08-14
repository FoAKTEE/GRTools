# Reduced Wigner-d matrix d^j_{m1,m2}(β) (general β) + the β=π/2 special case, and the fixed-angle
# rotation D^j_{m1,m2}(π,π/2,π/2). Promotes Punctures/src/mode_mixing.jl `_wignerd_halfpi` /
# `_wignerD_fixed` to the canonical Wigner home. The general-β form is the AHA consolidation:
# the β=π/2 collapse (trig product → 2^{-j}) is recovered as the wignerd_halfpi special case.
# Type-generic {T} (BigFloat-inherited); the factorials run in BigInt and narrow to T.

"""
    wignerd(::Type{T}, j, m1, m2, β) -> T

Reduced Wigner-d matrix element `d^j_{m1,m2}(β)` via Wigner's formula (the cos(β/2)/sin(β/2)
trig-product single sum). Returns `zero(T)` for `|m1| > j` or `|m2| > j`. The angle `β` is
converted to `T`, so BigFloat precision is inherited. `wignerd_halfpi` is the `β=π/2`
special case.
"""
function wignerd(::Type{T}, j::Integer, m1::Integer, m2::Integer, β) where {T}
    (abs(m1) > j || abs(m2) > j) && return zero(T)
    pref = sqrt(T(_fbig(j + m1)) * T(_fbig(j - m1)) * T(_fbig(j + m2)) * T(_fbig(j - m2)))
    βT = T(β)
    c = cos(βT / 2)
    s = sin(βT / 2)
    smin = max(0, m1 - m2)
    smax = min(j + m1, j - m2)
    acc = zero(T)
    for k in smin:smax
        denom = T(_fbig(j + m1 - k)) * T(_fbig(k)) * T(_fbig(m2 - m1 + k)) * T(_fbig(j - m2 - k))
        acc += T((-1)^(m2 - m1 + k)) / denom * c^(2j + m1 - m2 - 2k) * s^(m2 - m1 + 2k)
    end
    return pref * acc
end

"""
    wignerd_halfpi(::Type{T}, j, m1, m2) -> T

Reduced Wigner-d at `β=π/2`, where `cos(β/2)=sin(β/2)=1/√2` collapses the trig product to the
`s`-independent prefactor `2^{-j}`. Promoted from Punctures `_wignerd_halfpi`; exact rational
sum × √, narrowed to `T`.
"""
function wignerd_halfpi(::Type{T}, j::Integer, m1::Integer, m2::Integer) where {T}
    (abs(m1) > j || abs(m2) > j) && return zero(T)
    pref = sqrt(T(_fbig(j + m1)) * T(_fbig(j - m1)) * T(_fbig(j + m2)) * T(_fbig(j - m2)))
    smin = max(0, m1 - m2)
    smax = min(j + m1, j - m2)
    acc = zero(T)
    for s in smin:smax
        denom = T(_fbig(j + m1 - s)) * T(_fbig(s)) * T(_fbig(m2 - m1 + s)) * T(_fbig(j - m2 - s))
        acc += T((-1)^(m2 - m1 + s)) / denom
    end
    return pref * acc / T(2)^j
end

"""
    wignerD_fixed(::Type{T}, j, m1, m2) -> Complex{T}

Wigner-D at the fixed Euler angles `(π, π/2, π/2)` in Mathematica's convention, which factorises
as `WignerD[{j,m1,m2}, π, π/2, π/2] = (-1)^{m1} i^{m2} d^j_{m1,m2}(π/2)`. Promoted from Punctures
`_wignerD_fixed`; grounded against the Mathematica WignerD oracle fixture.
"""
function wignerD_fixed(::Type{T}, j::Integer, m1::Integer, m2::Integer) where {T}
    d = wignerd_halfpi(T, j, m1, m2)
    phase = Complex{T}((-1)^m1) * (Complex{T}(0, 1))^m2
    return phase * d
end

# Clebsch-Gordan coefficient ⟨j1 m1; j2 m2 | j3 m3⟩ — Racah formula, BigInt-exact.
#
# Promotes KerrLorenzCirc/src/angular.jl `cg` (the spheroidal→spherical mixing CG) to the
# canonical Wigner home (DAG §6 crosswalk), generalised off hard-Float64 to type-generic {T}.
# Reuses the shared _fbig / _sqrt_rat combinatorial helpers from wigner3j.jl (DRY).

"""
    clebsch_gordan(::Type{T}, j1, m1, j2, m2, j3, m3) -> T

Clebsch-Gordan coefficient `⟨j1 m1; j2 m2 | j3 m3⟩` via the Racah single-sum closed form,
exact in `Rational{BigInt}` and narrowed to `T` (so BigFloat precision is inherited). Returns
`zero(T)` on `m1+m2 ≠ m3`, triangle violation, or `|mᵢ| > jᵢ`. Integer arguments.
"""
function clebsch_gordan(::Type{T}, j1::Integer, m1::Integer, j2::Integer, m2::Integer,
                        j3::Integer, m3::Integer) where {T}
    (m1 + m2 != m3) && return zero(T)
    (j3 < abs(j1 - j2) || j3 > j1 + j2) && return zero(T)
    (abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3) && return zero(T)

    pref2 = (2j3 + 1) *
            (_fbig(j1 + j2 - j3) * _fbig(j1 - j2 + j3) * _fbig(-j1 + j2 + j3) //
             _fbig(j1 + j2 + j3 + 1)) *
            (_fbig(j1 + m1) * _fbig(j1 - m1) * _fbig(j2 + m2) * _fbig(j2 - m2) *
             _fbig(j3 + m3) * _fbig(j3 - m3))

    ssum = big(0) // 1
    for k in max(0, j2 - j3 - m1, j1 - j3 + m2):min(j1 + j2 - j3, j1 - m1, j2 + m2)
        ssum += (iseven(k) ? 1 : -1) //
                (_fbig(k) * _fbig(j1 + j2 - j3 - k) * _fbig(j1 - m1 - k) *
                 _fbig(j2 + m2 - k) * _fbig(j3 - j2 + m1 + k) * _fbig(j3 - j1 - m2 + k))
    end

    return _sqrt_rat(T, pref2) * T(ssum)
end

"""
    clebsch_gordan(j1, m1, j2, m2, j3, m3) -> Float64

Float64 convenience method matching the legacy `cg` signature (for facade re-export).
"""
clebsch_gordan(j1::Integer, m1::Integer, j2::Integer, m2::Integer, j3::Integer, m3::Integer) =
    clebsch_gordan(Float64, j1, m1, j2, m2, j3, m3)

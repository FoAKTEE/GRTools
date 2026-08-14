# 3-harmonic spin-weighted angular coupling — the canonical home of SecondOrderRicci `coupling`
# (ref-code/SecondOrderRicci/src/Coupling.cc, originally GSL gsl_sf_coupling_3j). Re-implemented
# over the in-package wigner3j so SecondOrderRicci drops its WignerSymbols [deps] (design §7.4).

"""
    coupling(::Type{T}, lc, mc, sc, la, ma, sa, lb, mb, sb) -> T

Angular 3-harmonic coupling coefficient: couples two spin-weighted harmonics `(la,ma,sa)` and
`(lb,mb,sb)` onto a third `(lc,mc,sc)`,

    (-1)^(mc+sc) √((2la+1)(2lb+1)(2lc+1)/4π) · 3j(lc la lb; sc -sa -sb) · 3j(lc la lb; -mc ma mb).

`T` sets the working precision (feeds both the √-normalisation and the two 3j evaluations).
Returns `zero(T)` on any selection-rule violation (triangle, m-sum ≠ 0, s-sum ≠ 0, or `|m|>j` /
`|s|>j` in a slot — the last guarded explicitly to match GSL, which returns 0 there).
"""
function coupling(::Type{T}, lc::Integer, mc::Integer, sc::Integer,
                  la::Integer, ma::Integer, sa::Integer,
                  lb::Integer, mb::Integer, sb::Integer) where {T}
    (abs(sc) <= lc && abs(sa) <= la && abs(sb) <= lb &&
     abs(mc) <= lc && abs(ma) <= la && abs(mb) <= lb) || return zero(T)
    sgn = isodd(mc + sc) ? -one(T) : one(T)
    norm = sqrt(T(2la + 1) * T(2lb + 1) * T(2lc + 1) / (4 * T(pi)))
    w_s = wigner3j(T, lc, la, lb, sc, -sa, -sb)
    w_m = wigner3j(T, lc, la, lb, -mc, ma, mb)
    return sgn * norm * w_s * w_m
end

"""
    coupling(lc, mc, sc, la, ma, sa, lb, mb, sb) -> Float64

Float64 convenience method matching the legacy `coupling` signature (for facade re-export).
"""
coupling(lc::Integer, mc::Integer, sc::Integer,
         la::Integer, ma::Integer, sa::Integer,
         lb::Integer, mb::Integer, sb::Integer) =
    coupling(Float64, lc, mc, sc, la, ma, sa, lb, mb, sb)

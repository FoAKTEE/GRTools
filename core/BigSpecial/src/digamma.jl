# Hand-rolled complex digamma ψ and a PolyGamma wrapper. SpecialFunctions.digamma has no
# Complex{BigFloat} method, so the WaSABI(Hybrid) forcing terms and the Punctures regsum tail
# need this. Promotes WaSABI `_digamma`/`PolyGamma`; the Bernoulli coefficients are written as
# exact ratios in the element type so a BigFloat argument keeps precision (no Float64 literal cap
# — the accuracy floor is the fixed shift-to-Re≥12 + 7-term asymptotic truncation, ≈1e-16).

"""
    digamma_big(z::Complex) -> typeof(z)

Digamma `ψ(z)` for complex `z` (incl. `Complex{BigFloat}`, which upstream lacks). Recurrence
`ψ(z) = ψ(z+1) − 1/z` shifts `Re z` up to `≥ 12`, then the Bernoulli asymptotic series
`ψ(z) = log z − 1/(2z) − Σ_{k≥1} B_{2k}/(2k z^{2k})` (through `B_14`). Type-generic in the
element float type.
"""
function digamma_big(z::Complex)
    T = real(typeof(z))
    r = zero(z)
    while real(z) < 12
        r -= one(z) / z
        z += one(z)
    end
    inv = one(z) / z
    inv2 = inv * inv
    r += log(z) - inv / 2
    t = inv2                                   # 1/z^{2k}, k=1..7  (coeffs = B_{2k}/(2k))
    r -= t / 12;                  t *= inv2     # B_2/2  = 1/12
    r += t / 120;                 t *= inv2     # B_4/4  = -1/120 (sign in the −Σ)
    r -= t / 252;                 t *= inv2     # B_6/6  = 1/252
    r += t / 240;                 t *= inv2     # B_8/8
    r -= t / 132;                 t *= inv2     # B_10/10
    r += (T(691) / 32760) * t;    t *= inv2     # B_12/12
    r -= t / 12                                 # B_14/14
    return r
end

"""
    PolyGamma(n::Integer, z) -> Complex

Polygamma `ψ^(n)(z)`. Only `n = 0` (digamma) is implemented — the order WaSABI-C's closed forms
need; any other order is an `ArgumentError` (preserved from the legacy contract).
"""
function PolyGamma(n::Integer, z::Number)
    n == 0 || throw(ArgumentError("PolyGamma order $n not implemented (n=0 only)"))
    return digamma_big(complex(z))
end

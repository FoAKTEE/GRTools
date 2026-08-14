# Pochhammer (rising-factorial) primitive. The leaf of the BigSpecial stratum: it has no
# special-function dependency and is referenced by the same product identity that biggamma /
# hyp2f1 build on. Consolidates Teukolsky `_pochhammer` and Punctures `_poch` (the SAME
# decision in two legacy packages) into one canonical type-generic form.

"""
    pochhammer(a, n::Integer) -> typeof(a)

Rising Pochhammer symbol `(a)_n`. For `n > 0` the explicit rising product
`∏_{k=0}^{n-1}(a+k)`; for `n = 0` the unit `one(a)`; for `n < 0` the reciprocal-falling
`(a)_{-n} = 1 / ∏_{k=1}^{n}(a-k)`. Type-generic in `a`, so `Complex{BigFloat}` precision is
inherited (no factorials / Γ — exact at any precision).
"""
function pochhammer(a, n::Integer)
    n == 0 && return one(a)
    p = one(a)
    if n > 0
        for k in 0:(n - 1)
            p *= (a + k)
        end
        return p
    else
        for k in 1:(-n)
            p *= (a - k)
        end
        return one(a) / p
    end
end

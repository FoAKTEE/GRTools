# Self-contained Gauss ₂F₁ — the Complex{BigFloat} fallback layer that frees the Perturbation /
# SWSH callers from upstream's precision gaps. Consolidates Teukolsky `_f21_maclaurin`/`_f21_zinf`
# + Punctures `_f21_term`. AHA: Teukolsky's `_safe_2f1` |z|>6 & near-integer branch-SELECTION
# heuristic is a per-caller decision and stays LOCAL to Teukolsky; the canonical `safe_2f1` here
# is the generic robust dispatcher (HypGeo pass-through, Pfaff+Maclaurin BigFloat fallback).

"""
    f21_maclaurin(a, b, c, z; maxterms=2_000_000) -> T

Gauss ₂F₁ Maclaurin series `Σ_{n≥0} (a)_n (b)_n / ((c)_n n!) zⁿ` (|z|<1) via the incremental
Pochhammer ratio `tₙ₊₁ = tₙ·(a+n)(b+n)/((c+n)(n+1))·z`. No factorials / Γ / ψ, so it is exact at
ANY precision; summed to working tolerance. Terminates naturally when `a` or `b` is a
non-positive integer. The working type promotes all of `a,b,c,z` (a ComplexF64 `z` with BigFloat
`a,b,c` still runs to BigFloat tolerance).
"""
function f21_maclaurin(a, b, c, z; maxterms::Integer=2_000_000)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
    aT = T(a); bT = T(b); cT = T(c); zT = T(z)
    tol = eps(real(T))
    term = one(T); s = term; n = 0
    while n < maxterms
        term *= (aT + n) * (bT + n) / ((cT + n) * (n + 1)) * zT
        s += term
        abs(term) <= tol * abs(s) && break
        n += 1
    end
    return s
end

"""
    f21_zinf(a, b, c, z; maxterms=200000) -> T

₂F₁ via the z→∞ connection (DLMF 15.8.2): two series in `1/z` (|1/z|<1 for |z|>1), FAST when
`|z| ≫ 1`. Non-degenerate form — needs `a−b ∉ ℤ` (the caller guards near-integer `a−b`). Γ via
`gamma_big`, so the real-negative `b−a` of real-ν arguments is handled.
"""
function f21_zinf(a, b, c, z; maxterms::Integer=200000)
    T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
    aT = T(a); bT = T(b); cT = T(c); zT = T(z); tol = eps(real(T))
    ser(p, q, r, x) = begin       # ₂F₁(p,q,r;x) Maclaurin (|x|=|1/z|<1)
        t = one(x); s = t; n = 0
        while n < maxterms
            t *= (p + n) * (q + n) / ((r + n) * (n + 1)) * x
            s += t; abs(t) <= tol * abs(s) && break; n += 1
        end
        s
    end
    iz = one(zT) / zT
    gamma_big(cT) * gamma_big(bT - aT) / (gamma_big(bT) * gamma_big(cT - aT)) * (-zT)^(-aT) *
        ser(aT, aT - cT + 1, aT - bT + 1, iz) +
    gamma_big(cT) * gamma_big(aT - bT) / (gamma_big(aT) * gamma_big(cT - bT)) * (-zT)^(-bT) *
        ser(bT, bT - cT + 1, bT - aT + 1, iz)
end

"""
    safe_2f1(a, b, c, z) -> Complex

Robust Gauss ₂F₁: byte-identical to `HypergeometricFunctions._₂F₁` for `Float64`/`ComplexF64`
(it just calls it). When upstream throws on a `Complex{BigFloat}` argument (the analytic-
continuation branches reach a `loggamma`/`digamma` of a `Complex{BigFloat}`), it falls back to
the Pfaff transform `z→w=z/(z−1)` (|w|<1 for Re z<1/2) + the no-Γ/ψ Maclaurin series, which is
robust for ALL `a,b,c` incl. integer `a−b`. The fast large-|z| `f21_zinf` route and the
selection heuristic between them are a per-caller decision and live in the caller, not here.
"""
function safe_2f1(a, b, c, z)
    try
        return HypergeometricFunctions._₂F₁(a, b, c, z)
    catch err
        (err isa DomainError || err isa MethodError) || rethrow(err)
        T = promote_type(typeof(a), typeof(b), typeof(c), typeof(z))
        aT = T(a); bT = T(b); cT = T(c); zT = T(z)
        w = zT / (zT - one(zT))
        return (one(zT) - zT)^(-aT) * f21_maclaurin(aT, cT - bT, cT, w)
    end
end

"""
    f21_terminating(l::Integer, c, z) -> T

Terminating Gauss ₂F₁`(−l, 1+l; c; z)` (`a=−l` a non-positive integer ⇒ finite `l+1`-term sum).
The closed-form angular-integral building block Punctures uses; type-generic in `c,z`.
"""
function f21_terminating(l::Integer, c, z)
    T = promote_type(typeof(c), typeof(z))
    cT = T(c); zT = T(z)
    a = -l; b = l + 1
    term = one(T); s = one(T)
    for k in 1:l
        term *= T(a + k - 1) * T(b + k - 1) / ((cT + (k - 1)) * T(k)) * zT
        s += term
    end
    return s
end

# Shared continued-fraction primitives (dependency-free, type-generic `{T}`, BigFloat-capable).
#
# Consolidates the private Teukolsky copy (`continued_fraction.jl/_cf_finite`,`_cf_convergent`) used by
# the MST renormalized-angular-momentum / Leaver paths into one reusable family, and adds the
# GRTOOLS_DESIGN §5.1 "shared Lentz" modified-Lentz iteration. All three evaluate the Gauss K-form
# (NO leading term): `a(n₀)/(b(n₀) + a(n₀+1)/(b(n₀+1) + ⋯))`.

"""
    cf_finite(f, g, n, nmax) -> Number

Finite Gauss continued fraction `ContinuedFractionK[f_i, g_i, {i, n, n+nmax}]` =
`f(n) / (g(n) + f(n+1)/(g(n+1) + ⋯ + f(m)/g(m)))`, `m = n + nmax`. `f`,`g` are functions of the
integer index. K-form: NO leading term, so `nmax = 0` collapses to `f(n)/g(n)`. (Teukolsky
`_cf_finite`.)
"""
function cf_finite(f, g, n::Integer, nmax::Integer)
    m = n + nmax
    acc = g(m)
    for i in (m - 1):-1:n
        acc = g(i) + f(i + 1) / acc
    end
    return f(n) / acc
end

"""
    cf_convergent(a, b, n0; rtol=1e-14, maxiter=100000) -> Number

Convergent Gauss continued fraction `a(n0)/(b(n0) + a(n0+1)/(b(n0+1) + ⋯))` via the fundamental
A/B recurrence (`A_k = b_k A_{k-1} + a_k A_{k-2}`, `B_k = b_k B_{k-1} + a_k B_{k-2}`; seeds
`A_{n0-2}=1, B_{n0-2}=0, A_{n0-1}=0, B_{n0-1}=1`; convergent `A_k/B_k`) iterated to relative
stability. Periodic rescaling (ratio-invariant) guards over/underflow. (Teukolsky `_cf_convergent`.)
"""
function cf_convergent(a, b, n0::Integer; rtol::Real=1e-14, maxiter::Integer=100_000)
    T = typeof(b(n0) / oneunit(b(n0)) + a(n0))
    Akm2, Bkm2 = one(T), zero(T)
    Akm1, Bkm1 = zero(T), one(T)
    ratio = zero(T)
    for k in n0:(n0 + maxiter)
        bk = b(k); ak = a(k)
        A = bk * Akm1 + ak * Akm2
        B = bk * Bkm1 + ak * Bkm2
        if !iszero(B)
            r = A / B
            k > n0 && abs(r - ratio) <= rtol * abs(r) && return r
            ratio = r
        end
        Akm2, Bkm2, Akm1, Bkm1 = Akm1, Bkm1, A, B
        sc = max(abs(Akm1), abs(Bkm1))                 # rescale (ratio-invariant) to avoid overflow
        if sc > 1e100 || (sc < 1e-100 && sc > 0)
            Akm2 /= sc; Bkm2 /= sc; Akm1 /= sc; Bkm1 /= sc
        end
    end
    error("cf_convergent did not converge in $maxiter iterations")
end

# Small nonzero placeholder used by modified Lentz to step over an exact zero partial denominator.
@inline _cf_tiny(::Type{T}) where {T} = eps(real(float(T)))^2

"""
    cf_lentz(a, b, n0; rtol=1e-14, maxiter=100000) -> Number

Modified-Lentz evaluation of the SAME Gauss K-form `a(n0)/(b(n0) + a(n0+1)/(b(n0+1) + ⋯))`
(GRTOOLS_DESIGN §5.1 shared-Lentz primitive). Iterates the running `C_j = b_j + a_j/C_{j-1}`,
`D_j = 1/(b_j + a_j D_{j-1})` and accumulates `f ← f·C_j·D_j`, stepping over exact zeros with a
tiny placeholder; stops when `|Δ_j − 1| ≤ rtol`. Numerically equals `cf_convergent` (different
iteration, same limit). Type-generic; BigFloat-capable.
"""
function cf_lentz(a, b, n0::Integer; rtol::Real=1e-14, maxiter::Integer=100_000)
    T = typeof(b(n0) / oneunit(b(n0)) + a(n0))
    tiny = T(_cf_tiny(T))
    f = tiny                                           # leading b_0 = 0 → tiny
    C = f
    D = zero(T)
    for j in 0:maxiter
        aj = a(n0 + j); bj = b(n0 + j)
        D = bj + aj * D
        iszero(D) && (D = tiny)
        D = one(T) / D
        C = bj + aj / C
        iszero(C) && (C = tiny)
        Δ = C * D
        f *= Δ
        abs(Δ - one(T)) <= rtol && return f
    end
    error("cf_lentz did not converge in $maxiter iterations")
end

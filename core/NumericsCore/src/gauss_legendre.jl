# Type-generic Gauss-Legendre quadrature rule (interior roots of P_n). NumericsCore previously had
# only the Gauss-LOBATTO rule (`_lgl_nodes_weights`, dg.jl, endpoints ±1 included); the compactified
# angular collocation ported from TeukEvolution.jl/src/Sphere.jl needs the plain Gauss-Legendre rule
# with NO endpoints. Mirrors the self-contained `_lgl_nodes_weights` pattern (Newton on the Legendre
# recurrence) and REUSES the same `_legendre_PdP` helper (dg.jl) — no FastGaussQuadrature dependency,
# no duplicated recurrence. Type-generic {T}, so BigFloat precision is inherited (the upstream
# FastGaussQuadrature.gausslegendre returned Float64 unconditionally).

"""
    gl_nodes_weights(::Type{T}, n) -> (x, w)

The `n`-point Gauss-Legendre quadrature rule on `[-1, 1]` in precision `T`, exact for polynomials of
degree `≤ 2n-1`. Nodes `x` are the `n` interior roots of the Legendre polynomial `P_n` (NO
endpoints — contrast [`_lgl_nodes_weights`](@ref) which pins `±1`), found by Newton iteration on
`P_n` (`δ = P_n(x)/P'_n(x)`, both from `_legendre_PdP`) starting from the asymptotic seed
`x_k ≈ -cos(π(4k-1)/(4n+2))`; weights `w_i = 2 / ((1 - x_i²) · P'_n(x_i)²)`. Nodes ascend `-1 → 1`,
matching `FastGaussQuadrature.gausslegendre(n)`. Self-contained (LinearAlgebra-only), so
`BigFloat`(256+) inputs give genuine BigFloat nodes/weights (converges to ~`eps(T)` residual).
"""
function gl_nodes_weights(::Type{T}, n::Int) where {T<:AbstractFloat}
    n >= 1 || throw(ArgumentError("number of nodes n must be ≥ 1"))
    x = Vector{T}(undef, n)
    w = Vector{T}(undef, n)
    @inbounds for k in 1:n
        # asymptotic seed, ascending order (-cos of the Chebyshev-like angle)
        xk = -cospi(T(4k - 1) / T(4n + 2))
        for _ in 1:100                       # quadratic convergence: a handful of iterations
            P, dP = _legendre_PdP(n, xk)
            δ = P / dP                       # Newton on P_n
            xk -= δ
            abs(δ) <= 4 * eps(T) * max(one(T), abs(xk)) && break
        end
        _, dP = _legendre_PdP(n, xk)
        x[k] = xk
        w[k] = 2 / ((1 - xk^2) * dP^2)
    end
    return x, w
end

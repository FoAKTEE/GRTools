# chebyshev.jl — Chebyshev–Gauss–Lobatto collocation. {T}-generic, LinearAlgebra-only.
# D2 is built by the Weideman–Reddy ℓ=2 recursion, NOT D1*D1, so D2 ≈ D1² is a genuine
# spectral-consistency gate. CGL barycentric weights: c_j = 1/w_j = (−1)^j·(2 endpoint / 1 interior).

"""
    cheb_lobatto(T, n; lo=zero(T), hi=one(T)) -> Vector{T}

Chebyshev–Gauss–Lobatto nodes on `[lo, hi]`:
`x_j = ½(lo+hi) + ½(hi−lo)·cos(jπ/n)`, `j = 0…n` (`n+1` nodes, DESCENDING — `x[1]=hi`,
`x[end]=lo`, both endpoints represented EXACTLY). The affine image of the reference `[−1,1]` CGL
grid; `[0, R₊]` is the hyperboloidal radial domain (scri `R=0` at `x[end]`, horizon `R=R₊` at
`x[1]`). {T}-generic (BigFloat-clean).
"""
function cheb_lobatto(::Type{T}, n::Integer; lo=zero(T), hi=one(T)) where {T<:Real}
    n >= 1 || throw(ArgumentError("cheb_lobatto needs n ≥ 1 (got $n)"))
    lo_ = T(lo); hi_ = T(hi)
    mid = (lo_ + hi_) / 2; half = (hi_ - lo_) / 2
    return T[mid + half * cos(T(j) * T(pi) / T(n)) for j in 0:n]
end

"""
    cheb_collocation(T, n; lo=zero(T), hi=one(T)) -> (x, D1, D2)

CGL nodes `x` (`cheb_lobatto`) plus the first- and second-derivative spectral differentiation
matrices on `[lo, hi]`. `D1` is the barycentric (Trefethen) first-derivative matrix; `D2` is
constructed by the Weideman–Reddy ℓ=2 recursion from `D1`'s off-diagonals and the negative-sum
diagonal — an INDEPENDENT build whose agreement with `D1*D1` (to the differentiation-matrix
conditioning floor `~κ·eps`) is a spectral-consistency diagnostic, not a definitional identity.
`D1*fvals` / `D2*fvals` are the spectral first/second derivatives, exact for polynomials of degree
≤ n (`D2` to degree ≤ n as well). {T}-generic: pass `T=BigFloat` for the BF-256 tier.
"""
function cheb_collocation(::Type{T}, n::Integer; lo=zero(T), hi=one(T)) where {T<:Real}
    x = cheb_lobatto(T, n; lo=lo, hi=hi)
    # c_j = 1/w_j (barycentric): (−1)^j × (2 at endpoints, 1 interior).
    c = T[(j == 0 || j == n) ? T(2) : one(T) for j in 0:n]
    for j in 0:n
        c[j+1] *= isodd(j) ? -one(T) : one(T)
    end
    D1 = zeros(T, n + 1, n + 1)
    @inbounds for i in 0:n, j in 0:n
        i != j && (D1[i+1, j+1] = (c[i+1] / c[j+1]) / (x[i+1] - x[j+1]))
    end
    for i in 0:n
        D1[i+1, i+1] = -sum(@view D1[i+1, :])          # negative-sum trick (row-sum ≡ 0)
    end
    D2 = zeros(T, n + 1, n + 1)
    @inbounds for i in 0:n, j in 0:n
        i != j && (D2[i+1, j+1] = 2 * D1[i+1, j+1] * (D1[i+1, i+1] - one(T) / (x[i+1] - x[j+1])))
    end
    for i in 0:n
        D2[i+1, i+1] = -sum(@view D2[i+1, :])
    end
    return x, D1, D2
end

"""
    cheb_solve(A, b) -> x

Dense LU solve of `A x = b` via a generic `lu` factorization — the collocation / shift-invert inner
solve for the hyperboloidal QNM path. Kept LinearAlgebra-only (no `\\`-triggered BLAS assumptions on
generic eltypes): works uniformly for `Float64`, `ComplexF64`, and `Complex{BigFloat}` (generic LU),
so the RQI polish inherits BigFloat precision. `b` may be a vector or a matrix.
"""
cheb_solve(A::AbstractMatrix, b) = lu(A) \ b

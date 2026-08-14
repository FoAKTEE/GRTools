# refelement.jl — 1D Legendre-Gauss reference element over NumericsCore.gl_nodes_weights:
# nodes/weights, barycentric Lagrange basis, differentiation matrix D (eq:diffmat), face
# interpolation to xi = +-1, tensor-product 2D indexing helpers.
# DAG node: mesh-basis-quadrature (elliptic-ssf GLOBAL_DAG row 8, TOPOLOGICAL_SCHEDULE W2-A).
#
# The CRP/spectre discretization collocates on Legendre-GAUSS points (interior roots of P_N,
# main.tex:556) — NOT Lobatto. NumericsCore supplies the type-generic quadrature rule
# (`gl_nodes_weights`) and the barycentric weight/eval primitives (`bary_weights`); everything
# here is exact-polynomial machinery derived from them, so BigFloat precision is inherited.
#
# NOTE on the name: NumericsCore also exports a `RefElement` (the 1D LGL *evolution* element in
# dg.jl). This module defines its own binding, which shadows the implicit import inside
# EllipticCore — per DESIGN.md §4.1 the two are distinct types with distinct node families;
# cross-package consumers qualify. Confirmed shadow-clean (no warning) on Julia 1.12.

export RefElement, refelement, interpolation_matrix, tensor_lin, tensor_ij, tensor_line

"""
    RefElement{RT<:AbstractFloat}

Immutable 1D Legendre-**Gauss** reference element on `ξ ∈ [-1, 1]` with `N = length(nodes)`
collocation points (the interior roots of the Legendre polynomial `P_N` — the element carries
**no boundary points**, contrast Lobatto). Fields:

- `nodes::Vector{RT}` — Gauss nodes, ascending, strictly inside `(-1, 1)`
- `weights::Vector{RT}` — Gauss quadrature weights, exact for polynomial degree `≤ 2N-1`
- `barywts::Vector{RT}` — barycentric weights of the node set (basis for all interpolation here)
- `D::Matrix{RT}` — nodal differentiation matrix, exact for degree `≤ N-1` (barycentric form
  with the negative-sum-trick diagonal, Berrut & Trefethen 2004 §9)
- `face_lo::Vector{RT}`, `face_hi::Vector{RT}` — face-interpolation vectors `ℓ_j(∓1)` /
  `ℓ_j(+1)` (DESIGN.md §4.1 `Iface_lo/hi`), exact for degree `≤ N-1`; consumed later by the
  DG lifting operators.

Generalized summation-by-parts (the precise identity for Gauss nodes, which exclude ±1): with
`W = Diagonal(weights)` and the face-lift matrix `B = face_hi*face_hi' - face_lo*face_lo'`,

    W*D + D'*W == B          (exact matrix identity)

because nodal vectors correspond bijectively to polynomials of degree `≤ N-1`, the integrand
`u'v + uv'` has degree `≤ 2N-3` (quadrature-exact), and boundary evaluation is interpolated.
The **pure** antisymmetry `W*D == -(W*D)'` is FALSE — its defect is exactly the rank-≤2 `B`.

Construct via [`refelement`](@ref) (cached); geometry is always a real type `RT = real(T)`.
"""
struct RefElement{RT<:AbstractFloat}
    nodes::Vector{RT}
    weights::Vector{RT}
    barywts::Vector{RT}
    D::Matrix{RT}
    face_lo::Vector{RT}
    face_hi::Vector{RT}
end

# Lagrange cardinal-function row [ℓ_1(xq) … ℓ_N(xq)] in second barycentric form; an exact node
# collision returns the exact unit row (the d == 0 branch of NumericsCore.bary_eval, matrixified).
function _lagrange_row(x::AbstractVector{RT}, bw::AbstractVector{RT}, xq::RT) where {RT}
    N = length(x)
    row = Vector{RT}(undef, N)
    @inbounds for j in 1:N
        if xq == x[j]
            fill!(row, zero(RT))
            row[j] = one(RT)
            return row
        end
    end
    s = zero(RT)
    @inbounds for j in 1:N
        t = bw[j] / (xq - x[j])
        row[j] = t
        s += t
    end
    row ./= s
    return row
end

function _build_refelement(::Type{RT}, N::Int) where {RT<:AbstractFloat}
    N >= 1 || throw(ArgumentError("RefElement needs N ≥ 1 nodes, got $N"))
    x, w = gl_nodes_weights(RT, N)
    bw = bary_weights(x)
    # barycentric differentiation matrix: D[i,j] = (bw[j]/bw[i])/(x[i]-x[j]) for i ≠ j,
    # diagonal by the negative-sum trick (rows annihilate constants exactly).
    D = Matrix{RT}(undef, N, N)
    @inbounds for i in 1:N
        acc = zero(RT)
        for j in 1:N
            i == j && continue
            Dij = (bw[j] / bw[i]) / (x[i] - x[j])
            D[i, j] = Dij
            acc += Dij
        end
        D[i, i] = -acc
    end
    return RefElement{RT}(x, w, bw, D,
                          _lagrange_row(x, bw, -one(RT)),
                          _lagrange_row(x, bw, one(RT)))
end

# ---- global element cache -----------------------------------------------------------------
# Keyed on (RT, N, bits) where RT = real(float(T)) — the ComplexF64 solver path and its Float64
# geometry share one entry — and bits = precision(BigFloat) for RT === BigFloat (a 256-bit
# element must never serve a 512-bit caller; fixed-width floats key bits = 0). Lock-protected:
# the DG matrix/apply layers above are threaded.
const _REFELEMENT_CACHE = Dict{Tuple{DataType,Int,Int},RefElement}()
const _REFELEMENT_LOCK = ReentrantLock()
const _refelement_hits = Ref(0)
const _refelement_misses = Ref(0)

"""
    refelement(T, N) -> RefElement{real(float(T))}

The `N`-point Legendre-Gauss [`RefElement`](@ref) for scalar type `T`, from the global cache
(`===`-identical on repeated calls). `T` may be complex (`ComplexF64` is the CRP solver type):
geometry lives in `RT = real(float(T))`. Cache keys include `precision(BigFloat)` when
`RT === BigFloat`, so elements built under different `setprecision` regimes never alias.
"""
function refelement(::Type{T}, N::Integer) where {T<:Number}
    RT = real(float(T))
    key = (RT, Int(N), RT === BigFloat ? precision(BigFloat) : 0)
    lock(_REFELEMENT_LOCK) do
        cached = get(_REFELEMENT_CACHE, key, nothing)
        if cached !== nothing
            _refelement_hits[] += 1
            return cached::RefElement{RT}
        end
        _refelement_misses[] += 1
        elem = _build_refelement(RT, Int(N))
        _REFELEMENT_CACHE[key] = elem
        return elem
    end
end

"""
    refelement_cache_stats() -> (hits = …, misses = …)

Counters of [`refelement`](@ref) cache hits/misses since the last
[`refelement_cache_reset!`](@ref) (unexported; test/diagnostic seam).
"""
refelement_cache_stats() = (hits = _refelement_hits[], misses = _refelement_misses[])

"""
    refelement_cache_reset!() -> nothing

Empty the [`refelement`](@ref) cache and zero its hit/miss counters (unexported; tests only).
"""
function refelement_cache_reset!()
    lock(_REFELEMENT_LOCK) do
        empty!(_REFELEMENT_CACHE)
        _refelement_hits[] = 0
        _refelement_misses[] = 0
    end
    return nothing
end

"""
    interpolation_matrix(elem::RefElement{RT}, pts) -> Matrix{RT}

Barycentric interpolation matrix `P` of size `(length(pts), N)` from the element's Gauss nodes
to the arbitrary points `pts` (converted to `RT`): `(P * f_nodes)[q] = ∑_j ℓ_j(pts[q]) f_j`,
exact for polynomials of degree `≤ N-1`. A query point that coincides exactly with a node
yields that node's exact unit row.
"""
function interpolation_matrix(elem::RefElement{RT}, pts::AbstractVector{<:Real}) where {RT}
    P = Matrix{RT}(undef, length(pts), length(elem.nodes))
    for (q, p) in enumerate(pts)
        P[q, :] = _lagrange_row(elem.nodes, elem.barywts, RT(p))
    end
    return P
end

# ---- tensor-product 2D index helpers --------------------------------------------------------
# Element-local ordering is Julia-column-major over (i, j) = (ξ-index, η-index): i fastest.
# Pure index arithmetic (no bounds checks — these are inner-loop primitives; callers own ranges).

"""
    tensor_lin(i, j, Nx) -> k

Linear index of tensor-product node `(i, j)` on an `Nx × Ny` element, column-major with `i`
(the ξ direction) fastest: `k = i + (j-1)*Nx`. Agrees with `LinearIndices((Nx, Ny))[i, j]`.
"""
tensor_lin(i::Integer, j::Integer, Nx::Integer) = Int(i) + (Int(j) - 1) * Int(Nx)

"""
    tensor_ij(k, Nx) -> (i, j)

Inverse of [`tensor_lin`](@ref): the `(i, j)` pair of linear index `k` on an element with
`Nx` nodes in the ξ direction.
"""
tensor_ij(k::Integer, Nx::Integer) = ((Int(k) - 1) % Int(Nx) + 1, (Int(k) - 1) ÷ Int(Nx) + 1)

"""
    tensor_line(dim, t, Nx, Ny) -> AbstractRange{Int}

Linear indices of the 1D node line running along dimension `dim` of an `Nx × Ny` element at
transverse index `t`:

- `dim == 1`: the ξ-line at fixed `j = t` — the contiguous range `(t-1)*Nx+1 : (t-1)*Nx+Nx`
  (`== reshape`'s column `A[:, t]`);
- `dim == 2`: the η-line at fixed `i = t` — the strided range `t : Nx : t+(Ny-1)*Nx`
  (`== A[t, :]`).

Ranges (not vectors), so `view`s over them are allocation-free in the DG kernels.
"""
function tensor_line(dim::Integer, t::Integer, Nx::Integer, Ny::Integer)
    if dim == 1
        lo = (Int(t) - 1) * Int(Nx)
        return (lo + 1):(lo + Int(Nx))
    elseif dim == 2
        return Int(t):Int(Nx):(Int(t) + (Int(Ny) - 1) * Int(Nx))
    end
    throw(ArgumentError("tensor_line: dim must be 1 or 2, got $dim"))
end

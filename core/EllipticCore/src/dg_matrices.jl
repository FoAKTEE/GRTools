# dg_matrices.jl — per-element DG matrix layer of CRP eq:dgres_strong_matrix (main.tex:647-664):
# diagonal mass M, stiffness/derivative MD_i, face lifting ML, derivative-transpose lift MLD^T,
# 1D differentiation matrix D (eq:diffmat) and volume→face interpolation I, with the Jacobian /
# surface-Jacobian / inverse-Jacobian factors of the tensor-product block maps (mesh2d.jl).
# DAG node: dg-operator-matrices (elliptic-ssf GLOBAL_DAG row 10, TOPOLOGICAL_SCHEDULE W3-A).
#
# Matrix-free discipline (main.tex:668): only the 1D constituents (D, weights, face-interpolation
# vectors) are cached globally — via the wave-2 `refelement` cache keyed (RT, N, precision) —
# and the small dense 2D matrices are formed per element on demand. No global sparse assembly
# lives here; that is the column-probe utility of the matrix-free-operator node.
#
# CRP blocks are tensor products of two 1D coordinate maps (mesh2d.jl Block), so the Jacobian
# matrix is DIAGONAL: ∂x_i/∂ξ^î = δ_i^î j_î, J = j_1 j_2, invjac^î_i = δ^î_i / j_î, and the
# surface Jacobian of a face with normal in dimension î is the transverse Jacobian Ĵ = j_ĵ,
# ĵ ≠ î. Under this structure the main.tex:649-652 formulas specialize per physical/logical
# dimension d (grid point p = (p₁, p₂), i-fastest tensor_lin ordering):
#
#   M[p,p]        = j₁(p₁) j₂(p₂) w_{p₁} w_{p₂}                          (main.tex:649)
#   MD_d[p,q]     = w_{p₁} w_{p₂} j_t(p_t) D_d[p_d, q_d] δ_{p_t q_t}     (main.tex:650, t ≠ d;
#                    the column-side J·invjac^d_d = j_t collapses onto p_t through the δ)
#   I⁽ᵈ'ˢ⁾[p,q̄]   = ℓ_{p_d}(s) δ_{p_t q̄}                                 (main.tex:661-664)
#   ML⁽ᵈ'ˢ⁾[p,q̄]  = ℓ_{p_d}(s) δ_{p_t q̄} j_t(q̄) w_{q̄}                    (main.tex:651)
#   MLDT⁽ᵈ'ˢ⁾[p,q̄]= ℓ′_{p_d}(s) δ_{p_t q̄} j_t(q̄) w_{q̄} / j_d(face)       (main.tex:652)
#
# with s = ±1 the face side, j_d(face) the composed normal Jacobian AT the face (the invjac
# factor of MLD^T is evaluated at the face point q̄, not at volume nodes), and
# ℓ′_k(s) = (Dᵀ · ℓ(s))_k — exact, since ℓ′_k has degree ≤ N−2 and Gauss interpolation is
# exact through degree N−1. Only the face-NORMAL derivative direction of MLD^T is built:
# the flux matrix 𝔄^{ij} of eq:dgres_strong is diagonal (main.tex:575) and the lifted term
# carries n_j, so for tensor-product maps every other component is multiplied by zero.
#
# 2D Jacobian-weighted generalized SBP (the load-bearing structural identity, extending the
# wave-2 1D identity  W·D + Dᵀ·W == face_hi·face_hiᵀ − face_lo·face_loᵀ): for each d,
#
#   MD_d + MD_dᵀ == ML⁽ᵈ'⁺⁾·(I⁽ᵈ'⁺⁾)ᵀ − ML⁽ᵈ'⁻⁾·(I⁽ᵈ'⁻⁾)ᵀ        (exact up to roundoff)
#
# because both sides factor per transverse index t into [W D + Dᵀ W − B]·j_t(t) w_t(t): the
# volume Jacobian and the inverse-Jacobian combine into exactly the surface Jacobian. Equally,
# MD_d == M · Dphys_d with Dphys_d = diag(1/j_d)(D_d ⊗ δ). Gated in test_dg_matrices.jl on
# affine AND compaction-mapped elements, Float64 (1e-13) and BigFloat-256 (1e-70).

export ElementMatrices, FaceMatrices, element_matrices, face_matrices

"""
    FaceMatrices{RT<:AbstractFloat}

The face-coupling matrices of one element face (CRP main.tex:651-652, 661-664), for the face
on side `side ∈ (-1, +1)` in logical dimension `dim ∈ (1, 2)`. With `n = Nx·Ny` volume points
and `nf` face points (the transverse dimension's collocation nodes — Legendre-Gauss faces at
`ξ = ±1` are OFF the volume grid, hence genuinely interpolatory):

- `I::Matrix{RT}` (`n × nf`) — volume→face interpolation: face values are `transpose(I) * u`;
  its transpose-action `I * g` scatters face data into the volume (the lifting skeleton).
- `ML::Matrix{RT}` (`n × nf`) — face lifting: interpolation × surface Jacobian × face weights,
  `ML = I · Diagonal(Ĵ(q̄) w_q̄)`. Applies the `(ψ_p, ·)_∂Ω` face inner product of
  eq:dgres_strong_matrix to face data.
- `MLDT::Matrix{RT}` (`n × nf`) — the derivative-transpose lift `MLD^T` applying
  `(∂_n ψ_p, ·)_∂Ω`: normal-derivative face row × surface Jacobian × face weights ×
  the inverse normal Jacobian evaluated AT the face (`njac`).
- `njac::RT` — the composed normal-dimension Jacobian at the face.

Sign conventions (outward normal `n_i`) are applied by the residual assembly (dg-strong-
residuals), not here: all three matrices are unsigned quadrature/interpolation objects.
"""
struct FaceMatrices{RT<:AbstractFloat}
    dim::Int
    side::Int
    njac::RT
    I::Matrix{RT}
    ML::Matrix{RT}
    MLDT::Matrix{RT}
end

"""
    ElementMatrices{RT<:AbstractFloat}

Per-element DG matrices of CRP eq:dgres_strong_matrix (main.tex:647-664) on the tensor-product
Legendre-Gauss grid of one mesh element (`N = (Nx, Ny)` points, `tensor_lin` i-fastest
ordering), geometry from a [`Block`](@ref) map composed with the element's dyadic
[`logical_interval`](@ref):

- `N::NTuple{2,Int}` — collocation points per dimension.
- `x::NTuple{2,Vector{RT}}` — physical coordinates of the per-dimension collocation nodes.
- `jac::NTuple{2,Vector{RT}}` — composed per-dimension Jacobians `j_d = dx_d/dξ_d` at the
  nodes (block-map Jacobian × half the dyadic segment width), strictly positive.
- `M::Diagonal{RT,Vector{RT}}` — mass matrix `M_pp = j₁ j₂ w_{p₁} w_{p₂}` (main.tex:649;
  exactly diagonal for Gauss collocation).
- `MD::NTuple{2,Matrix{RT}}` — stiffness/derivative matrices `MD_d` (main.tex:650), equal to
  `M · Dphys_d` with `Dphys_d` the chain-rule physical differentiation matrix.
- `faces::NTuple{2,NTuple{2,FaceMatrices{RT}}}` — `faces[d][s]` with `s = 1` the lower
  (`ξ_d = -1`) and `s = 2` the upper (`ξ_d = +1`) face; prefer [`face_matrices`](@ref).

Real-typed (`RT = real(float(T))` promoted with the block's geometry type): the complex solver
path multiplies complex nodal vectors by these real matrices. Satisfies, per dimension `d`,
the 2D Jacobian-weighted generalized-SBP identity
`MD_d + MD_dᵀ == ML⁽ᵈ'⁺⁾(I⁽ᵈ'⁺⁾)ᵀ − ML⁽ᵈ'⁻⁾(I⁽ᵈ'⁻⁾)ᵀ` (file header; gated to 1e-13/1e-70).
"""
struct ElementMatrices{RT<:AbstractFloat}
    N::NTuple{2,Int}
    x::NTuple{2,Vector{RT}}
    jac::NTuple{2,Vector{RT}}
    M::Diagonal{RT,Vector{RT}}
    MD::NTuple{2,Matrix{RT}}
    faces::NTuple{2,NTuple{2,FaceMatrices{RT}}}
end

# composed 1D geometry of the element along dimension `dim`: block map ∘ dyadic segment map.
# Returns (x nodes, jacobian at nodes, jacobian at the lo face, jacobian at the hi face).
function _dim_geometry(::Type{RT}, el::RefElement{RT}, m::CoordinateMap1D,
                       id::ElementId, dim::Int) where {RT}
    lo, hi = logical_interval(RT, id, dim)
    half = (hi - lo) / 2                                   # dyadic, exact
    ζ = [(lo * (1 - ξ) + hi * (1 + ξ)) / 2 for ξ in el.nodes]
    x = RT[map_forward(m, ζk) for ζk in ζ]
    j = RT[map_jacobian(m, ζk) * half for ζk in ζ]
    return x, j, RT(map_jacobian(m, lo) * half), RT(map_jacobian(m, hi) * half)
end

# one face (normal dimension `dim`, side ±1): interpolation, lifting, derivative-transpose
# lift from the normal element's cached face vector and the transverse Jacobian/weights.
function _face_matrices(dim::Int, side::Int, eln::RefElement{RT}, njac::RT,
                        jt::Vector{RT}, wt::Vector{RT}, Nx::Int, Ny::Int) where {RT}
    fv = side == -1 ? eln.face_lo : eln.face_hi            # ℓ_k(±1)
    dv = transpose(eln.D) * fv                             # ℓ′_k(±1), exact (degree ≤ N−2)
    n, nn, nf = Nx * Ny, length(fv), length(jt)
    If = zeros(RT, n, nf)
    ML = zeros(RT, n, nf)
    MLDT = zeros(RT, n, nf)
    for t in 1:nf
        s = jt[t] * wt[t]                                  # surface Jacobian × face weight
        for k in 1:nn
            p = dim == 1 ? tensor_lin(k, t, Nx) : tensor_lin(t, k, Nx)
            If[p, t] = fv[k]
            ML[p, t] = fv[k] * s
            MLDT[p, t] = dv[k] * s / njac
        end
    end
    return FaceMatrices{RT}(dim, side, njac, If, ML, MLDT)
end

"""
    element_matrices(T, blk::Block, id::ElementId, N::NTuple{2,Integer}) -> ElementMatrices

Build the per-element DG matrices (CRP main.tex:647-664) for element `id` of block `blk` with
`N = (Nx, Ny)` Legendre-Gauss points per dimension. `T` is the solver scalar type: matrices
are real, in `RT = promote_type(real(float(T)), block geometry type)` — a `ComplexF64` request
shares the `Float64` constituents, and BigFloat precision is inherited end to end.

The 1D constituents (differentiation matrix `D`, quadrature weights, face-interpolation
vectors) come from the global [`refelement`](@ref) cache — exactly one lookup per dimension —
and the 2D matrices are assembled per element (matrix-free discipline, main.tex:668).
Throws `ArgumentError` if `id` does not belong to `blk`.
"""
function element_matrices(::Type{T}, blk::Block{RTB}, id::ElementId,
                          N::NTuple{2,Integer}) where {T<:Number,RTB}
    id.block == blk.id ||
        throw(ArgumentError("element $(id) does not belong to block $(blk.id)"))
    RT = promote_type(real(float(T)), RTB)
    Nx, Ny = Int(N[1]), Int(N[2])
    el1 = refelement(RT, Nx)                               # cached 1D constituents
    el2 = refelement(RT, Ny)
    x1, j1, j1lo, j1hi = _dim_geometry(RT, el1, blockmap(blk, 1), id, 1)
    x2, j2, j2lo, j2hi = _dim_geometry(RT, el2, blockmap(blk, 2), id, 2)
    w1, w2 = el1.weights, el2.weights
    n = Nx * Ny

    m = Vector{RT}(undef, n)                               # M, main.tex:649
    for j in 1:Ny, i in 1:Nx
        m[tensor_lin(i, j, Nx)] = j1[i] * j2[j] * w1[i] * w2[j]
    end
    MD1 = zeros(RT, n, n)                                  # MD_d, main.tex:650 (tensor form)
    for t in 1:Ny, q in 1:Nx, p in 1:Nx
        MD1[tensor_lin(p, t, Nx), tensor_lin(q, t, Nx)] = w1[p] * w2[t] * j2[t] * el1.D[p, q]
    end
    MD2 = zeros(RT, n, n)
    for q in 1:Ny, p in 1:Ny, t in 1:Nx
        MD2[tensor_lin(t, p, Nx), tensor_lin(t, q, Nx)] = w1[t] * w2[p] * j1[t] * el2.D[p, q]
    end
    faces = ((_face_matrices(1, -1, el1, j1lo, j2, w2, Nx, Ny),
              _face_matrices(1, +1, el1, j1hi, j2, w2, Nx, Ny)),
             (_face_matrices(2, -1, el2, j2lo, j1, w1, Nx, Ny),
              _face_matrices(2, +1, el2, j2hi, j1, w1, Nx, Ny)))
    return ElementMatrices{RT}((Nx, Ny), (x1, x2), (j1, j2), Diagonal(m), (MD1, MD2), faces)
end

"""
    face_matrices(em::ElementMatrices, dim, side) -> FaceMatrices

The [`FaceMatrices`](@ref) of the face with normal in logical dimension `dim ∈ (1, 2)` on
side `side ∈ (-1, +1)` (mesh2d `neighbor_id` direction convention).
"""
function face_matrices(em::ElementMatrices, dim::Integer, side::Integer)
    (dim == 1 || dim == 2) ||
        throw(ArgumentError("face_matrices: dim must be 1 or 2, got $dim"))
    side == -1 || side == 1 ||
        throw(ArgumentError("face_matrices: side must be -1 or +1, got $side"))
    return em.faces[dim][side == -1 ? 1 : 2]
end

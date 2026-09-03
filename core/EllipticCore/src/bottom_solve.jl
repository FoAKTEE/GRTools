# bottom_solve.jl — the multigrid BOTTOM SOLVER of CRP main.tex:674: "At the coarsest
# level, the multigrid solver explicitly constructs the matrix A^bottom column-by-column
# and inverts it directly using the \texttt{Eigen} sparse LU decomposition routines
# (bottom solver). This helps provide an accurate large-scale baseline solution. The LU
# decomposition of A^bottom is cached and reused in all multigrid V-cycles across AMR
# iterations, since the bottom grid remains the same even when AMR adds finer levels."
# (fig:precond caption, main.tex:632: "The coarsest level is solved directly"; the yaml
# twin anc/ScalarSelfForce.yaml:190 `UseBottomSolver: True`.)
# DAG node: multigrid-bottom-solver (elliptic-ssf GLOBAL_DAG row 20, TOPOLOGICAL_SCHEDULE
# W8-A, CRITICAL PATH). Built TDD — gates calibrated pre-red in
# evidence/multigrid-bottom-solver/calibrate_bottom{1,2,3}.jl.
#
# ---- what this module is ----------------------------------------------------------------
# It fills the SEALED callable seam of `MultigridV(hier; bottom = (z, b) -> …)` — nothing in
# multigrid.jl changes, and its `dense_bottom_solver` default stays the reference the gates
# measure against. The route is the paper's: the sealed column probes [`build_matrix`](@ref)
# → `SparseArrays.sparse` → `lu` (UMFPACK through the stdlib, no new dependency) → a
# closure that applies the cached factorization.
#
# ---- the cache key (the whole design, in one place) --------------------------------------
# main.tex:674 caches the factorization "across AMR iterations". A cache that survives an
# AMR step cannot be keyed on object identity — every step builds a NEW mesh, a NEW operator
# and a NEW hierarchy — so it is keyed on what actually determines the matrix:
#
#   (1) the operator's scalar type `S`   — a Float64 and a ComplexF64 operator on the same
#       mesh are different matrices (and different factorization types);
#   (2) the coarsest MESH: per block `(id, map type, map parameters)`, the element id
#       vector, the per-element orders, and the block connections — everything the layout
#       and the geometry depend on;
#   (3) the penalty parameter `C`;
#   (4) the caller's `epoch` — the explicit invalidator, exactly as for
#       [`SubdomainCache`](@ref);
#   (5) an OPERATOR FINGERPRINT: `hash(A·v)` for one deterministic pseudo-random `v`, i.e.
#       ONE extra [`apply!`](@ref) of the n that assembling the matrix costs anyway
#       (the whole cache-HIT path measures 0.262 ms against a 17.8 ms build at the CRP
#       coarsest shape: 1.5 %).
#
# (5) is what makes the cache SELF-VALIDATING, and it is why this node's falsifier gate can
# be stated numerically: coefficients, hooks and boundary-condition TYPES are caller-owned
# closures that no structural key can see, so with (1)-(4) alone a coefficient change at an
# unchanged epoch silently reuses a stale factorization (measured: relative residual 3.5e-2
# against the new matrix). The fingerprint catches every such change that moves the matrix —
# calibrated: β 1e-10 relative, 𝔄¹¹ + 1e-12, C 1.5 → 2.0, and one boundary face switched
# Dirichlet → Neumann all rebuild, while an identically rediscretized operator hits. It is
# reliable BECAUSE the sealed rediscretization is bitwise reproducible (ESSF::multigrid-
# hierarchy M4: ‖A_c(hierarchy) − A_c(factory)‖_∞ = 0 exactly), which is re-measured here
# across 8 AMR rebuilds. It is a fingerprint, not a proof: two different matrices with the
# same image of one vector would collide (measure zero), and `epoch` remains the caller's
# explicit override.
#
# ---- scalar types ------------------------------------------------------------------------
# UMFPACK takes EXACTLY `Float64` and `ComplexF64`. `lu` of a sparse `BigFloat` matrix
# throws, and — the reason the dispatch is on the exact pair rather than on `AbstractFloat` —
# `lu` of a sparse `Float32`/`ComplexF32` matrix SILENTLY PROMOTES to a `Float64`/`ComplexF64`
# factorization, changing the scalar type under the caller's feet. Everything that is not
# exactly `{Float64, ComplexF64}` therefore goes to the generic DENSE LU (the coarsest grid
# is the irreducible block structure: small by construction), which keeps `Complex{BigFloat}`
# in `Complex{BigFloat}` at the inherited precision (measured 4.9e-77 at 256 bits).
#
# ---- cost (measured, evidence/measure_landing_w8a.jl) -------------------------------------
# At the CRP PRODUCTION coarsest shape (10 blocks × one 4×4 element = 160 dofs, 8.88 % dense)
# sparse and dense are within a factor ≈ 2 of each other in both directions, and both are
# ≈ 2 % of the column-probe assembly that precedes them (0.40 ms against 23.6 ms): the
# sparse route buys nothing measurable at that size. It is the right default anyway because
# it is the only one that scales — at 640 dofs (40 blocks) it factors 12.2× and solves 6.7×
# faster than dense, at 2560 dofs 172.9× and 38.6× (and its factorization holds 0.99 MB
# against the dense 100.0 MB) — which is why the paper's general-operator solver uses a
# sparse LU.
# The dense intermediate of `build_matrix` (O(n²) memory) is the remaining cost ceiling of
# the assembly; a directly column-wise sparse assembly is [FUTURE], not this node.

export BottomSolverCache, sparse_bottom_solver

"""
    BottomSolverCache()

Cache of coarsest-level ("bottom") LU factorizations, owned by the CALLER so that it
survives the hierarchy rebuild of an AMR step — the CRP main.tex:674 mechanism, "the LU
decomposition of ``A^\\mathrm{bottom}`` is cached and reused in all multigrid V-cycles
across AMR iterations, since the bottom grid remains the same even when AMR adds finer
levels". `entries` maps the key described below to a factorization; `builds` counts
factorizations constructed and `hits` counts reuses (the two counters are the executable
form of that claim: a mock 8-step AMR run that refines only the finest level measures
`builds == 1`, `hits == 7`).

The key of [`sparse_bottom_solver`](@ref) is

1. the operator's scalar type,
2. the mesh: per block `(id, map type, map parameters)`, the element id vector, the
   per-element orders, and the [`BlockConnection`](@ref)s,
3. the penalty parameter `C`,
4. the caller's `epoch`,
5. an operator FINGERPRINT `hash(A·v)` for one deterministic `v` — one extra
   [`apply!`](@ref), the whole hit path costing 1.5 % of a build (measured). It makes
   the cache self-validating against the coefficient/hook/boundary-condition closures
   no structural key can see: β perturbed by 1e-10 relative, `𝔄¹¹` by 1e-12, `C`
   changed, or a single boundary face switched Dirichlet → Neumann all rebuild, while
   an identically rediscretized operator hits.

The fingerprint is not a proof (two matrices with the same image of one vector would
collide — measure zero), so `epoch` stays the explicit invalidator, exactly as for
[`SubdomainCache`](@ref): bump it whenever you want a rebuild regardless. Sharing one
cache across operators of different scalar types is safe (the type is in the key).
"""
mutable struct BottomSolverCache
    entries::Dict{Any,Any}
    builds::Int
    hits::Int
end
BottomSolverCache() = BottomSolverCache(Dict{Any,Any}(), 0, 0)

# ---- the key ------------------------------------------------------------------------------

# a coordinate map as plain, structurally comparable data: its concrete type (which carries
# the geometry type) plus its parameters. Generic over any future CoordinateMap1D.
_bs_mapkey(m::CoordinateMap1D) =
    (typeof(m), ntuple(i -> getfield(m, i), fieldcount(typeof(m))))

# the mesh's identity: blocks (id + maps), the element layout, the orders, the seams.
# Vectors are COPIED into the key so a later mesh cannot mutate a live dictionary key.
function _bs_meshkey(mesh::Mesh2D)
    blocks = [(b.id, _bs_mapkey(blockmap(b, 1)), _bs_mapkey(blockmap(b, 2)))
              for b in mesh.blocks]
    conns = [(c.block, c.dim, c.neighbor) for c in mesh.connections]
    return (blocks, copy(mesh.elements), copy(mesh.orders), conns)
end

# one deterministic pseudo-random probe vector (Knuth-MMIX LCG, fixed seed — the same
# vector for a given (type, length) in every process and every run)
function _bs_probe_vector(::Type{S}, n::Int) where {S<:Number}
    RT = real(S)
    s = 0x9e3779b97f4a7c15
    nxt = () -> (s = 0x5851f42d4c957f2d * s + 0x14057b7ef767814f;
                 RT(Int64(s >> 11) / 2.0^53 * 2 - 1))
    return S[S <: Complex ? S(nxt(), nxt()) : S(nxt()) for _ in 1:n]
end

# the operator fingerprint: hash of A·v through the LINEARIZED matrix-free path. Serial by
# construction — `apply!`'s threaded path is bitwise identical to its serial one (sealed),
# and this must not depend on the thread count of the process that fills the cache.
function _bs_probe(op::MeshOperator)
    S = _scalar(op)
    v = _bs_probe_vector(S, op.mesh.ndofs)
    y = similar(v)
    apply!(y, op, v; threaded = false)
    return hash(y)
end

_bs_key(op::MeshOperator, epoch::Integer) =
    (_scalar(op), _bs_meshkey(op.mesh), op.C, Int(epoch), _bs_probe(op))

# ---- factorization ------------------------------------------------------------------------

# UMFPACK takes EXACTLY these two (see the file header: the small floats are silently
# promoted, which would change the caller's scalar type)
_bs_factor(A::SparseMatrixCSC{T}) where {T<:Union{Float64,ComplexF64}} = lu(A)
# everything else: the generic dense LU, in the caller's own scalar type and precision
_bs_factor(A::SparseMatrixCSC) = lu(Matrix(A))

# function barrier: the closure captures a CONCRETELY typed factorization even though the
# cache stores `Any`, so the apply path is allocation-free (measured 0 bytes, constant)
function _bs_apply(F, n::Int)
    return function (z::AbstractVector, b::AbstractVector)
        length(z) == n && length(b) == n ||
            throw(ArgumentError("bottom solver: expected length-$n vectors, got " *
                                "$(length(z)) and $(length(b))"))
        # `ldiv!(z, F, b)` is the allocation-free route and leaves `b` intact; the aliased
        # call falls back to the in-place form (UMFPACK copies internally: 2632 B)
        z === b ? ldiv!(F, z) : ldiv!(z, F, b)
        return z
    end
end

"""
    sparse_bottom_solver(op::MeshOperator; cache = BottomSolverCache(), epoch = 0)
        -> bottom!(z, b)

The CRP bottom solver (main.tex:674): the coarsest-level matrix assembled by the sealed
column probes ([`build_matrix`](@ref)), converted to `SparseArrays.sparse`, factored with
the stdlib sparse `lu` (UMFPACK — through `SparseArrays`, no new dependency), and returned
as the callable `bottom!(z, b)` that writes the solution of `op·z = b` into `z`. Relative
residual ‖A z − b‖/‖b‖ ≤ 6.3e-16 at ComplexF64 and Float64 on the calibration hierarchies
(gated < 1e-12), agreeing with the dense [`dense_bottom_solver`](@ref) to 1.1e-15.

This is the production filling of [`MultigridV`](@ref)'s `bottom` seam — no API change:

    cache = BottomSolverCache()                 # caller-owned: it survives AMR steps
    hier  = MultigridHierarchy(op, build_op; caches = subdomain_cache)
    mg    = MultigridV(hier; bottom = sparse_bottom_solver(hier.levels[end].op; cache))
    gmres(op, b; Minv! = mg)

`b` must already carry any [`operator_on_zero`](@ref) correction: like the smoother and the
V-cycle, this level is LINEARIZED-only. The solve is exact, so GMRES ∘ MultigridV keeps the
sealed iteration counts (7 at both 3 and 4 levels on the multigrid-hierarchy model). One
apply allocates NOTHING (measured 0 bytes, constant over repeats) and never touches the
cache; `bottom!(z, z)` is allowed and is bitwise the un-aliased call, at the cost of one
internal copy.

CACHING (the point of the node). The factorization is stored in `cache` under the key
documented on [`BottomSolverCache`](@ref) — scalar type, mesh identity, penalty, `epoch`,
and an operator fingerprint `hash(A·v)`. It is therefore reused across every V-cycle AND
across hierarchy rebuilds whose coarsest grid and coarsest operator are unchanged, which is
exactly the CRP claim that the bottom LU survives AMR iterations: a mock 8-step AMR run that
refines the finest level only measures `builds == 1, hits == 7`. Any change of the coarsest
matrix — an `epoch` bump, a coefficient/penalty/boundary-condition change, an order change,
or a coarsening that lands on a different bottom grid — rebuilds. The default `cache` is a
fresh one, so a caller who does not keep it gets no reuse (the [`SchwarzSmoother`](@ref)
convention).

SCALAR TYPES. UMFPACK handles exactly `Float64` and `ComplexF64`; every other type
(`BigFloat`, `Complex{BigFloat}`, the 32-bit floats) falls back to the generic DENSE LU in
the caller's own type and precision — the coarsest grid is small by construction. Measured
at `Complex{BigFloat}` with 256 bits: relative residual 4.9e-77.

COST (measured, not gated). At the CRP production coarsest shape (10 blocks × one 4×4
element = 160 dofs, 8.9 % dense) the sparse and dense factorizations are within a factor ≈ 2
of each other in both directions and each costs ≈ 2 % of the `build_matrix` assembly that
precedes it — at THAT size the sparse route buys nothing measurable, and it is the default
because it is the one that scales: 12.2× (factor) / 6.7× (solve) faster than dense at 640
dofs, 172.9× / 38.6× at 2560 dofs, where the sparse factorization also holds 0.99 MB against
the dense 100.0 MB.
"""
function sparse_bottom_solver(op::MeshOperator;
                              cache::BottomSolverCache = BottomSolverCache(),
                              epoch::Integer = 0)
    key = _bs_key(op, epoch)
    F = get(cache.entries, key, nothing)
    if F === nothing
        F = _bs_factor(sparse(build_matrix(op)))
        cache.entries[key] = F
        cache.builds += 1
    else
        cache.hits += 1
    end
    return _bs_apply(F, op.mesh.ndofs)
end

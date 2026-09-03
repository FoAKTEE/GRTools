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
# The KEY IS EXACTLY (1)-(4), and `epoch` is AUTHORITATIVE — the same contract as
# [`SubdomainCache`](@ref), deliberately: a solver stack carries ONE epoch for both caches
# and bumps it on ANY change of coefficients, hooks, boundary conditions or penalty.
#
# Stored WITH each entry, and never part of the key, is an OPERATOR FINGERPRINT `hash(A·v)`
# for one deterministic pseudo-random `v` — ONE extra [`apply!`](@ref) out of the n that
# assembling the matrix costs anyway (the whole cache-HIT path measures 0.262 ms against a
# 17.8 ms build at the CRP coarsest shape: 1.5 %). On a key hit whose fingerprint differs,
# the solver THROWS: the caller broke the epoch contract, and a loud ArgumentError naming
# the rule is the only safe answer — silently rebuilding would make `epoch` decorative and
# leave the SubdomainCache (which cannot check anything) holding stale factorizations for
# the same operator change.
#
# THE FINGERPRINT IS NOT A PROOF, and this module does not claim it self-validates. A
# single deterministic probe cannot certify matrix identity: for a face whose probe traces
# are (q1, q2), the legitimate linear-homogeneous hook adding q2·ψ_int[1] − q1·ψ_int[2] to
# the received data is bitwise zero ON THE PROBE (IEEE multiplication commutes) and nonzero
# in general — measured in evidence/multigrid-bottom-solver/calibrate_bottom4.jl: A changes
# by 1.02 relative, A·v does not move a bit, the cache hits, and the reused factorization
# then solves the wrong system at a relative residual of 1.83. That construction is pinned
# as a documented-limitation test. What the fingerprint DOES catch (calibrated: β perturbed
# by 1e-10 relative, 𝔄¹¹ by 1e-12, a boundary face switched Dirichlet → Neumann) it catches
# only because the sealed rediscretization is bitwise reproducible (ESSF::multigrid-
# hierarchy M4: ‖A_c(hierarchy) − A_c(factory)‖_∞ = 0 exactly), re-measured here across 8
# AMR rebuilds. It is a smoke alarm, not a lock.
#
# What it must NOT fire on is equally pinned: inhomogeneous boundary DATA. The linearized
# matrix is bitwise identical when Dirichlet data move (only `operator_on_zero` changes), so
# the entry is REUSED — the fingerprint rides the LINEARIZED path exactly like the matrix.
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
levels". `builds` counts factorizations constructed and `hits` counts reuses (the two
counters are the executable form of that claim: a mock 8-step AMR run that refines only
the finest level measures `builds == 1`, `hits == 7`).

`entries` maps the STRUCTURAL key

1. the operator's scalar type,
2. the mesh: per block `(id, map type, map parameters)`, the element id vector, the
   per-element orders, and the [`BlockConnection`](@ref)s,
3. the penalty parameter `C`,
4. the caller's `epoch`,

to `(fingerprint, factorization)`. **`epoch` is AUTHORITATIVE**, exactly as for
[`SubdomainCache`](@ref), and deliberately so: a solver stack carries ONE epoch for both
caches and bumps it on ANY change of coefficients, hooks, boundary conditions or penalty —
the things a structural key cannot see, because they are caller-supplied closures.

The fingerprint `hash(A·v)` stored beside each factorization (one extra [`apply!`](@ref);
a whole cache hit costs 1.5 % of a build) is a CONSISTENCY ASSERTION, not part of the key:
on a key hit whose fingerprint differs, [`sparse_bottom_solver`](@ref) THROWS an
`ArgumentError` naming the epoch rule instead of silently rebuilding or silently returning
a stale factorization. It is **not a proof and does not make the cache self-validating** —
a single deterministic probe cannot certify matrix identity, and the module header carries
the measured counterexample (a legitimate linear-homogeneous hook that changes `A` by 1.02
relative while leaving `A·v` bitwise fixed). Inhomogeneous boundary DATA correctly does NOT
fire it: the linearized matrix is bitwise unchanged and the entry is reused.

Lifetime, threads, singularity:

- entries are NEVER evicted; the caller owns the cache's lifetime and may drop it whenever
  the problem changes (the [`SubdomainCache`](@ref) convention). An hp-AMR step that
  p-raises the coarsest grid strands one factorization per step, which is why a long AMR
  run should hold the cache for one solve campaign rather than forever.
- NOT thread-safe: `entries`, `builds` and `hits` are mutated without a lock. Use one cache
  per thread, or serialize the construction calls (the V-cycle apply path never touches it).
- a singular coarse matrix raises `SingularException` from `lu`, propagated unchanged; no
  entry is stored and the counters do not move.
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
# vector for a given (type, length) in every process and every run). Two methods over a
# plain local `UInt64` state: a captured-and-mutated closure variable would be boxed, and
# the resulting type instability cost one whole extra `apply!` (measured 8784 B / 28 µs at
# n = 32, against 656 B / 0.4 µs here).
_bs_lcg(s::UInt64) = 0x5851f42d4c957f2d * s + 0x14057b7ef767814f
_bs_draw(s::UInt64) = (s, Int64(s >> 11) / 2.0^53 * 2 - 1)

function _bs_probe_vector(::Type{S}, n::Int) where {S<:Real}
    v = Vector{S}(undef, n)
    s = 0x9e3779b97f4a7c15
    for i in 1:n
        s = _bs_lcg(s)
        v[i] = S(_bs_draw(s)[2])
    end
    return v
end

function _bs_probe_vector(::Type{S}, n::Int) where {S<:Complex}
    RT = real(S)
    v = Vector{S}(undef, n)
    s = 0x9e3779b97f4a7c15
    for i in 1:n
        s = _bs_lcg(s)
        re = RT(_bs_draw(s)[2])
        s = _bs_lcg(s)
        v[i] = S(re, RT(_bs_draw(s)[2]))
    end
    return v
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

# The cache key is STRUCTURAL ONLY — the fingerprint is NOT part of it (H1): `epoch` is
# authoritative, exactly as for `SubdomainCache`, and the fingerprint is stored WITH the
# entry as a consistency assertion.
_bs_key(op::MeshOperator, epoch::Integer) =
    (_scalar(op), _bs_meshkey(op.mesh), op.C, Int(epoch))

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

CACHING (the point of the node). The factorization is stored in `cache` under the
STRUCTURAL key documented on [`BottomSolverCache`](@ref) — scalar type, mesh identity,
penalty, `epoch` — and is therefore reused across every V-cycle AND across hierarchy
rebuilds whose coarsest grid is unchanged, which is exactly the CRP claim that the bottom LU
survives AMR iterations: a mock 8-step AMR run that refines the finest level only measures
`builds == 1, hits == 7`. A change of the key — an `epoch` bump, a penalty change, an order
change, or a coarsening that lands on a different bottom grid — rebuilds cleanly.

**`epoch` is AUTHORITATIVE and is the caller's job**, exactly as for [`SubdomainCache`](@ref)
(carry ONE epoch for both caches): bump it on ANY change of coefficients, hooks, boundary
conditions or penalty. If you do not, and the change is one the stored fingerprint `hash(A·v)`
can see, this function THROWS an `ArgumentError` naming the rule rather than returning a
stale factorization — but the fingerprint is an assertion, NOT a proof, and cannot see every
change (the module header carries the measured counterexample). Inhomogeneous boundary DATA
is not a matrix change and correctly hits. The default `cache` is a fresh one, so a caller who
does not keep it gets no reuse (the [`SchwarzSmoother`](@ref) convention); entries are never
evicted, the cache is not thread-safe, and a singular coarse matrix propagates
`SingularException` from `lu`.

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
    fp = _bs_probe(op)
    entry = get(cache.entries, key, nothing)
    if entry === nothing
        F = _bs_factor(sparse(build_matrix(op)))
        cache.entries[key] = (fp, F)                     # the fingerprint rides along
        cache.builds += 1
        return _bs_apply(F, op.mesh.ndofs)
    end
    fp0, F = entry
    fp0 == fp || throw(ArgumentError(
        "sparse_bottom_solver: the coarsest operator changed but the cache key did not — " *
        "its fingerprint A·v differs from the cached one. `epoch` is AUTHORITATIVE: bump " *
        "the shared epoch on ANY change of coefficients, hooks, boundary conditions or " *
        "penalty (the same rule as SubdomainCache — use ONE epoch for both caches). This " *
        "check is an assertion, not a proof: it cannot detect every such change, so never " *
        "rely on it in place of the epoch."))
    cache.hits += 1
    return _bs_apply(F, op.mesh.ndofs)
end

"""
    sparse_bottom_solver(hier::MultigridHierarchy; cache = BottomSolverCache(), epoch = 0)
        -> bottom!(z, b)

The RECOMMENDED entry point: the bottom solver of `hier`'s coarsest level, i.e.
`sparse_bottom_solver(hier.levels[end].op; cache, epoch)` — it removes the one way the
operator form can be misused (handing it the FINEST operator, which builds an n × n dense
matrix by n column probes and then refuses every coarse-level right-hand side at call
time). Both forms share a cache and its entries verbatim.

    mg = MultigridV(hier; bottom = sparse_bottom_solver(hier; cache = bcache, epoch = e))
"""
sparse_bottom_solver(hier::MultigridHierarchy;
                     cache::BottomSolverCache = BottomSolverCache(),
                     epoch::Integer = 0) =
    sparse_bottom_solver(hier.levels[end].op; cache = cache, epoch = epoch)

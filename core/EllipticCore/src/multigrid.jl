# multigrid.jl — geometric h-coarsening multigrid V-cycle preconditioner of CRP
# main.tex:670-676: "Each GMRES iteration is preconditioned with one multigrid V-cycle,
# where the multigrid hierarchy is built by h-coarsening the domain up to the initial
# level (the irreducible block structure)"; an element-centered additive Schwarz smoother
# (2 passes, main.tex:676) runs on every finer level and a direct solve on the coarsest.
# The yaml twin (anc/ScalarSelfForce.yaml:186-194, mirrored by the ScalarSelfForce
# `MultigridConfig`): Iterations 1, PreSmoothing True, UseBottomSolver True,
# PostSmoothingAtBottom False, ComplexShift 0.0 (the multigrid-level shift is a documented
# no-op — the complex shift is SUBDOMAIN-level, schwarz_shift.jl).
# DAG node: multigrid-hierarchy (elliptic-ssf GLOBAL_DAG row 19, TOPOLOGICAL_SCHEDULE
# W7-A, CRITICAL PATH). Built TDD — gates calibrated pre-red in
# evidence/multigrid-hierarchy/calibrate_multigrid{,2,3,4}.jl.
#
# ---- coarsening (spectre LinearSolver::multigrid::coarsen, Hierarchy.hpp) --------------
# spectre decrements the refinement level uniformly over the domain, doing nothing for
# blocks already fully coarsened. That is defined on a UNIFORM per-block refinement spec;
# an AMR mesh needs a rule that also preserves the tiling and the 2:1 balance. The rule
# implemented here is the natural generalization: CAP the per-dimension refinement level
# at (max level − 1) and replace every deeper element by its ancestor AT the cap. Because
# `level ↦ min(level, cap)` is MONOTONE and every element maps to an ancestor:
#   * targets are equal or disjoint ⇒ the coarse elements tile exactly what the fine ones did;
#   * |min(a,cap) − min(b,cap)| ≤ |a − b| ⇒ within-block 2:1 balance AND the cross-block
#     transverse compatibility are preserved automatically (no repair loop);
#   * the per-dimension parent gap is always 0 or 1, so the transfer half-intervals are
#     exactly the mortar `MORTAR_LOWER`/`MORTAR_UPPER` cases;
#   * on a uniform mesh it reproduces spectre's decrement exactly, and iterating it
#     terminates at ONE element per block — the irreducible block structure (main.tex:672).
# Polynomial orders p-coarsen by the MIN rule ("we only p-coarsen the grid in the sense
# that we choose the smaller of the two polynomial degrees when combining elements",
# spectre Multigrid.hpp), which also guarantees the mortar precondition N_fine ≥ N_coarse.
#
# ---- transfers (spectre Spectral::projection_matrix_{parent_to_child,child_to_parent}) --
# Per fine element and per dimension the transfer operators ARE the sealed mortar
# projections of mortars.jl, with the coarse (parent) element playing the "face" and the
# fine (child) element the "mortar":
#   P = projection_prolong(RT, N_coarse, N_fine, msize)   — nodal Lagrange interpolation,
#       exact on the coarse element's polynomial space (spectre's parent_to_child);
#   R = projection_restrict(RT, N_coarse, N_fine, msize)  — the mass-conservative L2
#       adjoint jac·W_c⁻¹·Pᵀ·W_f (spectre's child_to_parent with operand_is_massive=false).
# The 2D operators are the tensor products (P₂ ⊗ P₁) and (R₂ ⊗ R₁), applied as
# U_f = P₁·U_c·P₂ᵀ on the i-fastest `tensor_lin` element blocks.
#
# MASSIVE vs NODAL restriction (BINDING; the Galerkin seam). This package's operator is
# the MASSIVE Gauss formulation (dg_operator.jl `R(Ψ) = M·S`, the CRP production config
# Massive: True), so a residual vector r represents the linear FUNCTIONAL v ↦ vᵀr. Its
# coarse-level image is therefore v_c ↦ (P v_c)ᵀ r, i.e.
#
#     b_coarse = Pᵀ · r_fine          (`mg_restrict_residual!`)
#
# — spectre's `operand_is_massive` branch, which is literally `blaze::trans(prolongation)`
# (Projection.cpp:107-118), and the restriction the V-cycle uses. It is EXACT for any
# block map. The nodal L2 restriction R is the transfer for FIELDS (and the object of the
# `mg_transfer_check` trio); using it on a massive residual double-counts the mass and
# degrades the contraction from 0.171 to 0.676 (calibrated falsifier, round 2D).
#
# SCOPE of the mass-weighted identities (measured, round 1): with the diagonal Legendre-
# GAUSS mass, `M_c·R == Pᵀ·M_f` and `Pᵀ == M_c ∘ R ∘ M_f⁻¹` hold to 5e-16 on AFFINE block
# maps but FAIL at O(1) (1.08 / 1.69 relative) on the CRP eq:compaction 1/r map: the Gauss
# rule integrates the polynomial part exactly, not a non-constant Jacobian. Hence
# `mg_transfer_check` reports adjointness in the LOGICAL (dyadic-measure) mass, which is
# exact for every map and equals the physical L2 inner product on affine blocks — and the
# V-cycle restricts with Pᵀ, which needs no mass at all. Likewise "prolongation is exact on
# polynomials" is a statement in the element's LOGICAL coordinate (a physical-coordinate
# polynomial is not in the coarse space when the map is curved: measured 1.2e+03).
#
# ---- V-cycle ---------------------------------------------------------------------------
# Descend: pre-smooth (2 Schwarz passes), restrict the residual with Pᵀ, recurse from a
# zero coarse guess. Bottom: the `bottom` seam (default: a dense LU of build_matrix,
# factored ONCE at construction) — `PostSmoothingAtBottom: False` means the bottom result
# is prolongated straight up. Ascend: add the prolongated correction, post-smooth. Started
# from z = 0 the whole map is a FIXED linear operator, so ordinary right-preconditioned
# GMRES stays valid (no flexible GMRES needed) — `MultigridV` is itself the `Minv!`
# callable of `gmres`, exactly like `schwarz_preconditioner`.
#
# Measured (evidence/multigrid-hierarchy/calibration.log): contraction 0.1715 / 0.1718 /
# 0.1719 at 3 / 4 / 5 levels (LEVEL-INDEPENDENT, spread 4e-4) on a two-block complex model;
# GMRES ∘ MultigridV needs 7 iterations at both 3 and 4 levels where GMRES ∘ Schwarz needs
# 18 and 32 and unpreconditioned GMRES 213 / >400. The 2 Schwarz passes are load-bearing
# (1 pass: 0.306; 3 passes: 0.105).
#
# Threading: this layer is SERIAL — the only threaded work is inside the sealed `apply!`
# and `smooth!`, whose threaded paths are gated BITWISE against serial, so a V-cycle is
# bitwise independent of the thread count (calibrated: identical iterate hash at -t1 and
# -t4). BigFloat is never threaded (owner directive), inherited through those defaults.

export MultigridHierarchy, MultigridV, coarsen_mesh, mg_prolong!,
       mg_restrict_residual!, mg_restrict_field!, mg_transfer_check,
       dense_bottom_solver, vcycle!

# ---- h-coarsening ---------------------------------------------------------------------

"""
    coarsen_mesh(mesh::Mesh2D) -> (coarse::Mesh2D, parent::Vector{Int})  |  nothing

One geometric h-coarsening step: merge every COMPLETE SIBLING GROUP by one refinement
level and restore the two-to-one balance with the sealed [`balance_two_to_one`](@ref).
Returns the coarse [`Mesh2D`](@ref) — same blocks, same [`BlockConnection`](@ref)s, orders
p-coarsened by the MIN rule over the merged children (spectre `Multigrid.hpp`: "we choose
the smaller of the two polynomial degrees when combining elements") — together with
`parent[e]`, the coarse-element index each fine element `e` belongs to, or `nothing` when
no further progress is possible (one element per block: the irreducible block structure of
CRP main.tex:672).

A group is the 2, 2 or 4 elements that share an address except for the last index bit of
the coarsened dimension(s), at IDENTICAL levels. Groups are tried largest-mask first (both
dimensions, then dimension 1, then dimension 2) and are CONSUMED as a whole: letting each
element pick its own mask independently produces OVERLAPPING targets (with A, B, C present
and D refined, A picks {1} → A∪B while C picks {2} → A∪C). Elements with no complete group
are carried over unchanged, so a coarse element is always an ancestor-or-equal of the fine
elements it covers, at a per-dimension level gap of 0 or 1 — exactly the `MORTAR_FULL` /
`MORTAR_LOWER` / `MORTAR_UPPER` cases the mortar transfers implement.

Merging enlarges elements, so previously non-adjacent elements can become face-adjacent
with a level gap of 2; the sealed ANISOTROPIC [`balance_two_to_one`](@ref) repairs that by
splitting only the offending dimension of the offending element. Because the fine mesh is
itself a balanced refinement of the merged set, the repair never splits past the fine
level, and every coarse element stays an ancestor-or-equal of fine elements (gated: a
300-mesh census over AMR-reachable anisotropic meshes builds 300/300 valid, balanced,
nested, exactly-tiling hierarchies). A step that cannot make progress returns `nothing`.

On a uniform mesh this reproduces spectre's uniform decrement exactly (64 → 16 → 4 → 1,
`Hierarchy.cpp:20-46`). The earlier per-dimension level CAP is NOT equivalent: it is
monotone only for pairs that are already face-adjacent, and on the CRP's own anisotropic
AMR class (`anc/ScalarSelfForce.yaml:145` `Isotropy: Anisotropic`) it produced nested,
non-tiling ancestors — 244 of 300 census meshes failed, and an 18-element hand
counterexample threw "invalid face coverage".
"""
function coarsen_mesh(mesh::Mesh2D{RT}) where {RT}
    isempty(mesh.elements) && return nothing
    merged = _mg_merge_targets(mesh)
    cand = ElementId[]
    seen = Set{ElementId}()
    for el in mesh.elements
        t = merged[el]
        if !(t in seen)
            push!(seen, t)
            push!(cand, t)
        end
    end
    ids = balance_two_to_one(cand; connections = mesh.connections)
    length(ids) == length(mesh.elements) && Set(ids) == Set(mesh.elements) && return nothing
    lookup = Dict{ElementId,Int}(id => k for (k, id) in enumerate(ids))
    parent = Vector{Int}(undef, length(mesh.elements))
    for (e, el) in enumerate(mesh.elements)
        anc = _amr_ancestor(lookup, el)
        anc === nothing &&
            throw(ArgumentError("coarsen_mesh: fine element $el has no ancestor in the " *
                                "coarsened set — this is a bug in the coarsening rule"))
        for d in 1:2
            Δ = el.segments[d].level - anc.segments[d].level
            Δ <= 1 ||
                throw(ArgumentError("coarsen_mesh: level gap $Δ between $el and $anc " *
                                    "(the transfers require 0 or 1)"))
        end
        parent[e] = lookup[anc]
    end
    ords = Dict{ElementId,NTuple{2,Int}}()
    for (e, k) in enumerate(parent)
        o = mesh.orders[e]
        prev = get(ords, ids[k], (typemax(Int), typemax(Int)))
        ords[ids[k]] = (min(prev[1], o[1]), min(prev[2], o[2]))
    end
    cmesh = Mesh2D(mesh.blocks, ids, ords; connections = mesh.connections)
    return cmesh, parent
end

_mg_up(s::SegmentId) = SegmentId(s.level - 1, s.index >> 1)

# the members of `el`'s sibling group under `mask` (identical levels, low bits varied)
function _mg_group(el::ElementId, mask::NTuple{2,Bool})
    s1, s2 = el.segments
    b1 = mask[1] ? (s1.index >> 1) << 1 : s1.index
    b2 = mask[2] ? (s2.index >> 1) << 1 : s2.index
    out = ElementId[]
    for j in (mask[2] ? (0:1) : (0:0)), i in (mask[1] ? (0:1) : (0:0))
        push!(out, ElementId(el.block, SegmentId(s1.level, b1 + i),
                             SegmentId(s2.level, b2 + j)))
    end
    return out
end

_mg_group_parent(el::ElementId, mask::NTuple{2,Bool}) =
    ElementId(el.block, mask[1] ? _mg_up(el.segments[1]) : el.segments[1],
              mask[2] ? _mg_up(el.segments[2]) : el.segments[2])

# fine element => its merge target. Groups are consumed as a whole, and the identity is
# assigned only AFTER every merge has been tried, so an element is never consumed as a
# singleton while it is still needed to complete another element's group.
function _mg_merge_targets(mesh::Mesh2D)
    present = Set(mesh.elements)
    consumed = Set{ElementId}()
    target = Dict{ElementId,ElementId}()
    for el in mesh.elements
        el in consumed && continue
        s1, s2 = el.segments
        for mask in ((true, true), (true, false), (false, true))
            (mask[1] && s1.level == 0) && continue
            (mask[2] && s2.level == 0) && continue
            grp = _mg_group(el, mask)
            all(g -> (g in present) && !(g in consumed), grp) || continue
            p = _mg_group_parent(el, mask)
            for g in grp
                push!(consumed, g)
                target[g] = p
            end
            break
        end
    end
    for el in mesh.elements                       # identity for whatever did not merge
        haskey(target, el) || (target[el] = el)
    end
    return target
end

# ---- per-step transfer matrices --------------------------------------------------------

# per fine element, per dimension: (prolongation coarse→fine, nodal L2 restriction fine→coarse)
function _mg_transfer_matrices(::Type{RT}, mesh::Mesh2D, cmesh::Mesh2D,
                               parent::Vector{Int}) where {RT}
    n = length(mesh.elements)
    P = Vector{NTuple{2,Matrix{RT}}}(undef, n)
    R = Vector{NTuple{2,Matrix{RT}}}(undef, n)
    for e in 1:n
        el = mesh.elements[e]
        c = cmesh.elements[parent[e]]
        nf = mesh.orders[e]
        nc = cmesh.orders[parent[e]]
        ops = ntuple(2) do d
            ms = _mg_msize(el.segments[d], c.segments[d])
            (projection_prolong(RT, nc[d], nf[d], ms),
             projection_restrict(RT, nc[d], nf[d], ms))
        end
        P[e] = (ops[1][1], ops[2][1])
        R[e] = (ops[1][2], ops[2][2])
    end
    return P, R
end

# which half of the parent the child covers (one coarsening STEP gives a gap of 0 or 1)
function _mg_msize(sf::SegmentId, sc::SegmentId)
    Δ = sf.level - sc.level
    Δ == 0 && return MORTAR_FULL
    Δ == 1 || throw(ArgumentError("multigrid transfer: parent level gap $Δ (expected 0 or 1)"))
    return isodd(sf.index) ? MORTAR_UPPER : MORTAR_LOWER
end

# the child's logical sub-interval inside the parent's [-1, 1], for any level gap (a
# hierarchy LEVEL may compose several coarsening steps — see the `shrink` policy);
# diagnostics only ([`mg_transfer_check`](@ref)).
function _mg_child_interval(::Type{RT}, sf::SegmentId, sc::SegmentId) where {RT}
    Δ = sf.level - sc.level
    Δ == 0 && return (-one(RT), one(RT))
    w = ldexp(one(RT), 1 - Δ)                     # width 2/2^Δ, exact
    lo = -one(RT) + (sf.index - (sc.index << Δ)) * w
    return (lo, lo + w)
end

# ---- hierarchy -------------------------------------------------------------------------

"""
    MultigridLevel

One level of a [`MultigridHierarchy`](@ref). Public fields:

- `mesh::Mesh2D` — this level's mesh (level 1 is the caller's finest mesh).
- `op::MeshOperator` — the operator on this level: the SAME caller coefficient closures
  rediscretized here (built by the hierarchy's `build_op` factory; level 1 is the
  caller's own operator object).
- `smoother::SchwarzSmoother` — this level's smoother, with its OWN
  [`SubdomainCache`](@ref) (one LU per subdomain, built once at hierarchy construction).
- `parent::Vector{Int}` — `parent[e]` is the index in the NEXT (coarser) level's mesh of
  the element fine element `e` merges into; empty on the coarsest level.
- `P`, `R` — per fine element, the per-dimension prolongation / L2-restriction matrices
  ([`mg_prolong!`](@ref), [`mg_restrict_residual!`](@ref),
  [`mg_restrict_field!`](@ref)); empty on the coarsest level. A level may compose several
  coarsening steps — see the `shrink` policy of [`MultigridHierarchy`](@ref).
"""
struct MultigridLevel{OP,SM,RT<:AbstractFloat}
    mesh::Mesh2D{RT}
    op::OP
    smoother::SM
    parent::Vector{Int}
    P::Vector{NTuple{2,Matrix{RT}}}
    R::Vector{NTuple{2,Matrix{RT}}}
end

"""
    MultigridHierarchy(op::MeshOperator, build_op; maxlevels = nothing, overlap = 2,
                       shift = 0, epoch = 0) -> MultigridHierarchy

The geometric h-coarsening multigrid hierarchy of `op` (CRP main.tex:672): level 1 is
`op` on its own mesh, and each further level is [`coarsen_mesh`](@ref) applied to the one
above, down to the irreducible block structure (one element per block) or to `maxlevels`
levels. `hier.levels` is the vector of [`MultigridLevel`](@ref)s, finest first.

`maxlevels` counts the TOTAL number of levels (`≥ 1`, the finest included). The yaml
`LinearSolver.Multigrid.InitialCoarseLevels: L` counts the COARSER levels only (spectre
`Multigrid/Tags.hpp:30-37`, `ElementsAllocator.hpp:128-136`: `num_levels = min(…, value + 1)`;
`ScalarSelfForce` stores it verbatim as `MultigridConfig.initial_coarse_levels`), so the
yaml twin is `maxlevels = L + 1`.

`build_op(mesh::Mesh2D) -> MeshOperator` is the caller's operator FACTORY: it must
rediscretize the same problem on the mesh it is handed — the same coefficient closures,
penalty `C` and hook/BC policy, with the boundary conditions rebuilt for that mesh's own
[`boundary_faces`](@ref) (boundary-condition tables are keyed by element id, so they
cannot be reused verbatim across levels). Physics stays entirely caller-side; this package
never inspects the coefficients.

Every level gets a [`SchwarzSmoother`](@ref) with the given `overlap`, subdomain `shift`
and coefficient `epoch` — exactly one LU per subdomain, built here and reused by every
V-cycle. `caches` chooses where those factorizations live and is what makes main.tex:676
("the matrix inverses can remain cached across AMR iterations") true ACROSS hierarchy
rebuilds:

- `nothing` (default) — a fresh [`SubdomainCache`](@ref) per level: nothing survives a rebuild.
- a single `SubdomainCache` — shared by every level and reusable across AMR steps. This is
  the recommended form: the cache key is `(element id, subdomain shape, epoch, shift)`, so
  entries from different levels and different AMR iterations can never collide, and a level
  whose neighbourhood is unchanged is all hits. Measured on the 4-level model: rebuilding on
  identical meshes costs 0 new factorizations and 170 hits; after refining the FINEST level
  only, 512 new factorizations are built and all 170 coarser subdomains are reused.
- a `Vector{SubdomainCache}` — one cache per level, extended in place as levels are built
  (note that adding an AMR level shifts the level indices, so entries then miss rather than
  collide; a single shared cache is preferred).

The `epoch` stays CALLER-managed for every form: bump it whenever coefficients, hooks,
boundary conditions or the penalty change — the cache cannot see them.

`shrink ∈ (0, 1]` is the LEVEL policy: a coarsening step that fails to shrink the problem
to `shrink · ndofs` is not worth a smoother, so further steps are taken and their transfer
matrices COMPOSED (`P = P₁·P₂`, `R = R₂·R₁`) before the level is materialized. The default
`1` never skips (one coarsening step = one level). Measured tradeoff on the graded puncture
cascade (uniform base + corner refinement rounds): `shrink = 1` gives 10 levels and a total
cost `Σ ndofs / ndofs_fine` of 2.84 at contraction 0.1712; `shrink = 0.5` gives 5 levels and
1.72 at the same 0.1712. On uniform meshes the two agree exactly. On randomly ANISOTROPIC
meshes skipping costs convergence (0.168 → 0.222 and 0.242 → 0.387), which is why the
conservative `1` is the default. The `epoch` discipline is the caller's, as for the
smoother itself: bump it whenever coefficients, hooks, boundary conditions or the penalty
change. `shift` is the SUBDOMAIN complex shift of [`helmholtz_shift`](@ref)'s literal seam;
the yaml's multigrid-level `ComplexShift` is a documented no-op (adopted decision 8).

Consumed by [`MultigridV`](@ref) (the V-cycle preconditioner), [`mg_prolong!`](@ref),
[`mg_restrict_residual!`](@ref), [`mg_restrict_field!`](@ref) and
[`mg_transfer_check`](@ref).
"""
struct MultigridHierarchy{L}
    levels::Vector{L}
end

function MultigridHierarchy(op::MeshOperator{S,RT}, build_op;
                            maxlevels::Union{Nothing,Integer} = nothing,
                            overlap::Integer = 2, shift::Real = 0,
                            epoch::Integer = 0, shrink::Real = 1,
                            caches = nothing) where {S,RT}
    maxlevels === nothing || maxlevels >= 1 ||
        throw(ArgumentError("MultigridHierarchy: maxlevels must be ≥ 1, got $maxlevels"))
    0 < shrink <= 1 ||
        throw(ArgumentError("MultigridHierarchy: shrink must be in (0, 1], got $shrink"))
    built = Any[]
    mesh = op.mesh
    cur = op
    while true
        co = (maxlevels !== nothing && length(built) + 1 >= Int(maxlevels)) ? nothing :
             _mg_coarsen_level(RT, mesh, shrink)
        sm = SchwarzSmoother(cur; overlap = overlap, shift = shift, epoch = epoch,
                             cache = _mg_cache(caches, length(built) + 1))
        if co === nothing
            push!(built, MultigridLevel(mesh, cur, sm, Int[],
                                        NTuple{2,Matrix{RT}}[], NTuple{2,Matrix{RT}}[]))
            break
        end
        cmesh, parent, P, R = co
        push!(built, MultigridLevel(mesh, cur, sm, parent, P, R))
        mesh = cmesh
        cur = build_op(cmesh)::MeshOperator
    end
    levels = identity.(built)                 # narrow to the concrete level type
    return MultigridHierarchy{eltype(levels)}(levels)
end

# the cache for level ℓ: a fresh one per level (default), one shared cache, or the ℓ-th
# entry of a caller-owned vector (extended in place as levels are built)
_mg_cache(::Nothing, ::Int) = SubdomainCache()
_mg_cache(c::SubdomainCache, ::Int) = c
function _mg_cache(c::AbstractVector, ℓ::Int)
    while length(c) < ℓ
        push!(c, SubdomainCache())
    end
    return c[ℓ]::SubdomainCache
end

# ONE HIERARCHY LEVEL = as many coarsening STEPS as it takes to shrink the problem by the
# factor `shrink`, with the per-dimension transfer matrices COMPOSED across the skipped
# steps (P = P₁·P₂, R = R₂·R₁, so R∘P = R₂(R₁P₁)P₂ = I still). `shrink = 1` (the default)
# never skips: every coarsening step is a level.
function _mg_coarsen_level(::Type{RT}, mesh::Mesh2D, shrink::Real) where {RT}
    co = coarsen_mesh(mesh)
    co === nothing && return nothing
    cur, par = co
    P, R = _mg_transfer_matrices(RT, mesh, cur, par)
    while cur.ndofs > shrink * mesh.ndofs
        co2 = coarsen_mesh(cur)
        co2 === nothing && break
        nxt, par2 = co2
        Pn, Rn = _mg_transfer_matrices(RT, cur, nxt, par2)
        P = [(P[e][1] * Pn[par[e]][1], P[e][2] * Pn[par[e]][2]) for e in eachindex(P)]
        R = [(Rn[par[e]][1] * R[e][1], Rn[par[e]][2] * R[e][2]) for e in eachindex(R)]
        par = [par2[par[e]] for e in eachindex(par)]
        cur = nxt
    end
    return cur, par, P, R
end

# ---- transfers -------------------------------------------------------------------------

_mg_scratch(::Type{T}, a::Int, b::Int) where {T} = Matrix{T}(undef, a, b)

"""
    mg_prolong!(u_f, hier::MultigridHierarchy, ℓ, u_c) -> u_f

Prolongate the coarse-level (`ℓ+1`) vector `u_c` onto level `ℓ`, OVERWRITING `u_f`: per
fine element, the nodal Lagrange interpolation `U_f = P₁·U_c·P₂ᵀ` of its parent's
polynomial onto its own Gauss grid (spectre `projection_matrix_parent_to_child`; the
sealed [`projection_prolong`](@ref) with the parent as face and the child as mortar).
EXACT on the coarse element's polynomial space — in the element's LOGICAL coordinate
(file header). This is the V-cycle's correction transfer.
"""
function mg_prolong!(uf::AbstractVector, hier::MultigridHierarchy, ℓ::Integer,
                     uc::AbstractVector)
    lev = hier.levels[ℓ]
    isempty(lev.parent) &&
        throw(ArgumentError("mg_prolong!: level $ℓ is the coarsest — nothing to prolong from"))
    cm = hier.levels[ℓ+1].mesh
    length(uf) == lev.mesh.ndofs ||
        throw(ArgumentError("mg_prolong!: u_f has length $(length(uf)), level $ℓ has " *
                            "$(lev.mesh.ndofs) dofs"))
    length(uc) == cm.ndofs ||
        throw(ArgumentError("mg_prolong!: u_c has length $(length(uc)), level $(ℓ+1) has " *
                            "$(cm.ndofs) dofs"))
    T = eltype(uf)
    fill!(uf, zero(T))
    tmp = _mg_scratch(T, maximum(o[1] for o in lev.mesh.orders),
                      maximum(o[2] for o in cm.orders))
    @inbounds for e in eachindex(lev.mesh.elements)
        k = lev.parent[e]
        ncx, ncy = cm.orders[k]
        nfx, nfy = lev.mesh.orders[e]
        Uc = reshape(view(uc, _element_range(cm, k)), ncx, ncy)
        Uf = reshape(view(uf, _element_range(lev.mesh, e)), nfx, nfy)
        P1, P2 = lev.P[e]
        W = view(tmp, 1:nfx, 1:ncy)
        mul!(W, P1, Uc)
        mul!(Uf, W, transpose(P2))
    end
    return uf
end

"""
    mg_restrict_residual!(b_c, hier::MultigridHierarchy, ℓ, r_f) -> b_c

Restrict the level-`ℓ` MASSIVE residual `r_f` onto the coarse level `ℓ+1`, OVERWRITING
`b_c` and accumulating over the children of each coarse element: the GALERKIN transfer
`b_c = Pᵀ r_f`, applied per element as `B_c += P₁ᵀ·R_f·P₂`.

This is the V-cycle's residual transfer and the ONLY restriction on its path. This
package's operator is the MASSIVE Gauss formulation (`dg_operator.jl`, `R(Ψ) = M·S`), so a
residual vector represents the linear FUNCTIONAL `v ↦ vᵀr_f`; its image on the coarse space
is `v_c ↦ (P v_c)ᵀ r_f`, i.e. exactly `Pᵀ r_f`. That is spectre's `operand_is_massive`
branch of `Spectral::projection_matrix_child_to_parent` (`Projection.cpp:107-118`, literally
the transpose of the prolongation), and it is EXACT for any block map — no mass matrix
enters. Restricting a massive residual with the NODAL [`mg_restrict_field!`](@ref) instead
double-counts the mass and degrades the V-cycle contraction from 0.171 to 0.676 (gated
falsifier).
"""
mg_restrict_residual!(vc::AbstractVector, hier::MultigridHierarchy, ℓ::Integer,
                      vf::AbstractVector) = _mg_restrict_impl!(vc, hier, ℓ, vf, true)

"""
    mg_restrict_field!(v_c, hier::MultigridHierarchy, ℓ, v_f) -> v_c

Restrict the level-`ℓ` NODAL FIELD `v_f` onto the coarse level `ℓ+1`, OVERWRITING `v_c`
and accumulating over the children of each coarse element: the mass-conservative L2
projection `V_c += R₁·V_f·R₂ᵀ` built from the sealed [`projection_restrict`](@ref)
(spectre's `projection_matrix_child_to_parent` with `operand_is_massive = false`).

This is the transfer for FIELDS, not for residuals: it satisfies `R∘P = I` and the
adjointness `⟨R v, u⟩_c = ⟨v, P u⟩_f` in the LOGICAL (dyadic-measure) Gauss mass
`2^{-level}·w` — which is the physical L2 mass of an AFFINE block up to its constant
Jacobian and, unlike the physical mass, makes the identity exact for CURVED maps too (the
physical-mass form fails at O(1) on the CRP `eq:compaction` 1/r block: 1.08 relative,
pinned by the test suite). Both properties are reported by [`mg_transfer_check`](@ref).
Use [`mg_restrict_residual!`](@ref) for the massive residual of a V-cycle.
"""
mg_restrict_field!(vc::AbstractVector, hier::MultigridHierarchy, ℓ::Integer,
                   vf::AbstractVector) = _mg_restrict_impl!(vc, hier, ℓ, vf, false)

function _mg_restrict_impl!(vc::AbstractVector, hier::MultigridHierarchy, ℓ::Integer,
                            vf::AbstractVector, massive::Bool)
    lev = hier.levels[ℓ]
    isempty(lev.parent) &&
        throw(ArgumentError("mg_restrict: level $ℓ is the coarsest — nothing to restrict to"))
    cm = hier.levels[ℓ+1].mesh
    length(vf) == lev.mesh.ndofs ||
        throw(ArgumentError("mg_restrict: v_f has length $(length(vf)), level $ℓ has " *
                            "$(lev.mesh.ndofs) dofs"))
    length(vc) == cm.ndofs ||
        throw(ArgumentError("mg_restrict: v_c has length $(length(vc)), level $(ℓ+1) has " *
                            "$(cm.ndofs) dofs"))
    T = eltype(vc)
    fill!(vc, zero(T))
    tmp = _mg_scratch(T, maximum(o[1] for o in cm.orders),
                      maximum(o[2] for o in lev.mesh.orders))
    @inbounds for e in eachindex(lev.mesh.elements)
        k = lev.parent[e]
        ncx, ncy = cm.orders[k]
        nfx, nfy = lev.mesh.orders[e]
        Vf = reshape(view(vf, _element_range(lev.mesh, e)), nfx, nfy)
        Vc = reshape(view(vc, _element_range(cm, k)), ncx, ncy)
        W = view(tmp, 1:ncx, 1:nfy)
        if massive
            P1, P2 = lev.P[e]
            mul!(W, transpose(P1), Vf)
            mul!(Vc, W, P2, true, true)
        else
            R1, R2 = lev.R[e]
            mul!(W, R1, Vf)
            mul!(Vc, W, transpose(R2), true, true)
        end
    end
    return vc
end

# deterministic diagnostic generator (no Random dependency — the dep budget is
# {LinearAlgebra, SparseArrays, NumericsCore}); Knuth-MMIX LCG in [-1, 1]
mutable struct _MGDiagRNG
    s::UInt64
end
function _mg_rand!(r::_MGDiagRNG, ::Type{RT}) where {RT}
    r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f
    return RT(Int64(r.s >> 11) / 2.0^53 * 2 - 1)
end

"""
    mg_transfer_check(hier::MultigridHierarchy) -> (; prolong_exact, restrict_identity,
                                                      adjointness)

The three transfer diagnostics of the multigrid hierarchy, as the WORST value over all
level pairs (deterministic — an internal LCG seeds the probe vectors), in the hierarchy's
real type (`BigFloat` in, `BigFloat` out):

- `prolong_exact` — max nodal error of prolonging a separable polynomial of the coarse
  element's degree, sampled on the coarse grid and compared against its exact values on
  the fine grid. The polynomial lives in the element's LOGICAL coordinate: that is the
  coarse polynomial SPACE for any block map (a physical-coordinate polynomial is not in
  it when the map is curved).
- `restrict_identity` — max‖(R∘P − I) u_c‖∞ on random coarse data: the children of each
  coarse element reconstruct the coarse space exactly.
- `adjointness` — the relative defect of `⟨R v, u⟩_c = ⟨v, P u⟩_f` in the mass-weighted
  inner product, with the LOGICAL (dyadic-measure) Gauss mass `2^{-level}·w`. That mass
  is the physical L2 mass of an AFFINE block up to its constant Jacobian, and unlike the
  physical mass it makes the identity exact for CURVED maps too (file header scope note).

All three are ≤ 1.6e-15 at Float64 and ≤ 1.1e-76 at BigFloat(256) on the calibration
meshes (affine multi-block, 2:1 h + p nonconforming, and the 1/r compactified block);
they are the falsifiers for a wrong transfer sign or a dropped half-interval Jacobian.
"""
function mg_transfer_check(hier::MultigridHierarchy)
    RT = _mg_realtype(hier)
    pe = zero(RT)
    ri = zero(RT)
    ad = zero(RT)
    rng = _MGDiagRNG(0x9e3779b97f4a7c15)
    for ℓ in 1:(length(hier.levels)-1)
        lev = hier.levels[ℓ]
        fm = lev.mesh
        cm = hier.levels[ℓ+1].mesh
        # (a) prolongation exactness in the coarse element's logical coordinate
        nmax = maximum(max(o[1], o[2]) for o in cm.orders)     # per-dimension coeff blocks
        cf = [_mg_rand!(rng, RT) for _ in 1:2nmax]
        f(ξ, η, nx, ny) = sum(cf[i+1] * ξ^i for i in 0:nx-1) *
                          sum(cf[nmax+i+1] * η^i for i in 0:ny-1)
        uc = zeros(RT, cm.ndofs)
        for k in eachindex(cm.elements)
            nx, ny = cm.orders[k]
            ξ = refelement(RT, nx).nodes
            η = refelement(RT, ny).nodes
            rgn = _element_range(cm, k)
            for j in 1:ny, i in 1:nx
                uc[rgn[tensor_lin(i, j, nx)]] = f(ξ[i], η[j], nx, ny)
            end
        end
        uf = zeros(RT, fm.ndofs)
        mg_prolong!(uf, hier, ℓ, uc)
        for e in eachindex(fm.elements)
            el = fm.elements[e]
            c = cm.elements[lev.parent[e]]
            nx, ny = fm.orders[e]
            ncx, ncy = cm.orders[lev.parent[e]]
            lo1, hi1 = _mg_child_interval(RT, el.segments[1], c.segments[1])
            lo2, hi2 = _mg_child_interval(RT, el.segments[2], c.segments[2])
            ξ = refelement(RT, nx).nodes
            η = refelement(RT, ny).nodes
            rgn = _element_range(fm, e)
            for j in 1:ny, i in 1:nx
                ξp = (lo1 * (1 - ξ[i]) + hi1 * (1 + ξ[i])) / 2
                ηp = (lo2 * (1 - η[j]) + hi2 * (1 + η[j])) / 2
                pe = max(pe, abs(uf[rgn[tensor_lin(i, j, nx)]] - f(ξp, ηp, ncx, ncy)))
            end
        end
        # (b) R∘P = I on random coarse data
        ucr = [_mg_rand!(rng, RT) for _ in 1:cm.ndofs]
        mg_prolong!(uf, hier, ℓ, ucr)
        back = zeros(RT, cm.ndofs)
        mg_restrict_field!(back, hier, ℓ, uf)
        ri = max(ri, maximum(abs.(back .- ucr)))
        # (c) adjointness in the logical (dyadic-measure) mass
        Wf = _mg_logical_mass(RT, fm)
        Wc = _mg_logical_mass(RT, cm)
        v = [_mg_rand!(rng, RT) for _ in 1:fm.ndofs]
        u = [_mg_rand!(rng, RT) for _ in 1:cm.ndofs]
        Rv = zeros(RT, cm.ndofs)
        mg_restrict_field!(Rv, hier, ℓ, v)
        Pu = zeros(RT, fm.ndofs)
        mg_prolong!(Pu, hier, ℓ, u)
        rhs = dot(Pu, Wf .* v)
        ad = max(ad, abs(dot(u, Wc .* Rv) - rhs) / abs(rhs))
    end
    return (prolong_exact = pe, restrict_identity = ri, adjointness = ad)
end

_mg_realtype(hier::MultigridHierarchy) = typeof(hier.levels[1].mesh).parameters[1]

# the LOGICAL (dyadic-measure) Gauss mass of every dof: 2^{-level_d}·w per dimension
function _mg_logical_mass(::Type{RT}, mesh::Mesh2D) where {RT}
    out = Vector{RT}(undef, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        nx, ny = mesh.orders[e]
        wx = refelement(RT, nx).weights .* ldexp(one(RT), -el.segments[1].level)
        wy = refelement(RT, ny).weights .* ldexp(one(RT), -el.segments[2].level)
        rgn = _element_range(mesh, e)
        for j in 1:ny, i in 1:nx
            out[rgn[tensor_lin(i, j, nx)]] = wx[i] * wy[j]
        end
    end
    return out
end

# ---- bottom solver seam -----------------------------------------------------------------

"""
    dense_bottom_solver(op::MeshOperator) -> bottom!(z, b)

The default coarsest-level solver: `lu(build_matrix(op))` factored ONCE here, returned as
the callable `bottom!(z, b)` that writes the exact solution of `op·z = b` into `z`
(generic — a `Complex{BigFloat}` operator factors and solves in `Complex{BigFloat}`).
Relative residual ≤ 1.3e-15 on the calibration hierarchies.

This is the API SEAM of [`MultigridV`](@ref)'s `bottom` keyword: any callable
`(z, b) -> z` solving the coarsest system may be substituted — in particular the cached
sparse UMFPACK factorization of the multigrid-bottom-solver node (CRP main.tex:674
"\\texttt{Eigen} sparse LU ... cached and reused in all multigrid V-cycles across AMR
iterations"). `b` must already carry any [`operator_on_zero`](@ref) correction: like the
smoother, this level is LINEARIZED-only.
"""
function dense_bottom_solver(op::MeshOperator)
    F = lu(build_matrix(op))
    n = size(F, 1)
    return function (z::AbstractVector, b::AbstractVector)
        length(z) == n && length(b) == n ||
            throw(ArgumentError("bottom solver: expected length-$n vectors, got " *
                                "$(length(z)) and $(length(b))"))
        copyto!(z, b)
        ldiv!(F, z)
        return z
    end
end

# ---- the V-cycle preconditioner ---------------------------------------------------------

"""
    MultigridV(hier::MultigridHierarchy; passes = 2, iterations = 1, pre_smoothing = true,
               post_smoothing_at_bottom = false,
               bottom = dense_bottom_solver(hier.levels[end].op)) -> MultigridV

The multigrid V-cycle preconditioner over `hier` (CRP main.tex:670-676). Defaults are the
CRP production configuration (anc/ScalarSelfForce.yaml:186-197): ONE V-cycle per
application (`iterations`, yaml `Multigrid.Iterations: 1`), `pre_smoothing` on
(`PreSmoothing: True`), a direct bottom solve (`UseBottomSolver: True` — pass
`bottom = nothing` for the `False` twin, which smooths at the bottom instead), NO
post-smoothing at the bottom (`PostSmoothingAtBottom: False`), and `passes = 2` Schwarz
smoothing passes per level (yaml `SchwarzSmoother.Iterations: 2`). The 2 passes are
load-bearing: 1 pass measures a contraction of 0.306, 2 passes 0.171, 3 passes 0.105.

`MultigridV` IS the right-preconditioner callable: `mg(z, v)` writes one full V-cycle
started from `z = 0` into `z`, so it plugs straight into the `Minv!` hook of
[`gmres`](@ref) exactly like [`schwarz_preconditioner`](@ref):

    hier = MultigridHierarchy(op, mesh -> MeshOperator(T, mesh, A, β, γ; C, bcs = bcs(mesh)))
    gmres(op, b; Minv! = MultigridV(hier))

Started from zero the V-cycle is a FIXED linear operator, so ordinary (non-flexible)
right-preconditioned GMRES stays valid. Measured: 7 GMRES iterations at both 3 and 4
levels, against 18 and 32 for the Schwarz smoother alone.

LINEARIZED-only, inherited from [`smooth!`](@ref): the source `b` must already carry the
[`operator_on_zero`](@ref) correction of any fixed boundary data — the coarse levels solve
for a CORRECTION (zero initial guess, restricted-residual source), so the V-cycle is a
strictly LINEAR map and no fixed data may enter it (gated: `M(2v) == 2M(v)` bitwise on an
operator whose fine AND coarse levels carry nonzero Dirichlet data).

An apply allocates a CONSTANT number of bytes after warm-up but not zero: `MultigridV`
preallocates its own per-level buffers here, and 100 % of the remaining allocation is the
per-call workspace of the sealed [`smooth!`](@ref)/[`apply!`](@ref). Preallocating those is
the solver-stack node's scope. The per-level work
buffers are preallocated here, so one `MultigridV` is NOT safe to apply from two threads
at once (build one per thread if you need that); the V-cycle itself is serial and its
result is bitwise independent of the thread count (file header).

Apply with [`vcycle!`](@ref) to smooth an existing iterate in place.
"""
struct MultigridV{H<:MultigridHierarchy,B,T<:Number}
    hier::H
    passes::Int
    iterations::Int
    pre_smoothing::Bool
    post_smoothing_at_bottom::Bool
    bottom::B
    x::Vector{Vector{T}}
    b::Vector{Vector{T}}
    r::Vector{Vector{T}}
end

function MultigridV(hier::MultigridHierarchy; passes::Integer = 2,
                    iterations::Integer = 1, pre_smoothing::Bool = true,
                    post_smoothing_at_bottom::Bool = false,
                    bottom = dense_bottom_solver(hier.levels[end].op))
    passes >= 1 ||
        throw(ArgumentError("MultigridV: passes must be ≥ 1, got $passes"))
    iterations >= 1 ||
        throw(ArgumentError("MultigridV: iterations must be ≥ 1, got $iterations"))
    T = _scalar(hier.levels[1].op)
    x = [zeros(T, l.mesh.ndofs) for l in hier.levels]
    b = [zeros(T, l.mesh.ndofs) for l in hier.levels]
    r = [zeros(T, l.mesh.ndofs) for l in hier.levels]
    return MultigridV{typeof(hier),typeof(bottom),T}(hier, Int(passes), Int(iterations),
                                                     pre_smoothing,
                                                     post_smoothing_at_bottom, bottom,
                                                     x, b, r)
end

# BigFloat is never threaded (owner directive); mirrors the schwarz/matrixfree default
_mg_default_threaded(::MultigridV{H,B,T}, u, b) where {H,B,T} =
    Threads.nthreads() > 1 && isbitstype(promote_type(T, eltype(u), eltype(b)))

"""
    vcycle!(u, mg::MultigridV, b; threaded = auto) -> u

Run `mg.iterations` multigrid V-cycles of `mg` on `A·u = b` (the LINEARIZED operator of
the hierarchy's finest level), updating `u` IN PLACE from its current value — the
stationary multigrid iteration, whose asymptotic error-contraction factor is 0.172 on the
calibration model, independent of the number of levels.

One V-cycle descends the hierarchy — pre-smooth with `mg.passes` Schwarz passes, restrict
the residual with the Galerkin transfer `Pᵀ` ([`mg_restrict_residual!`](@ref)), recurse from a zero
coarse guess — solves the coarsest level with the `bottom` seam, then ascends, adding the
prolongated correction ([`mg_prolong!`](@ref)) and post-smoothing on every non-bottom
level. `threaded` is forwarded to the sealed [`apply!`](@ref)/[`smooth!`](@ref) and
defaults to `Threads.nthreads() > 1` for bits scalar types (never for `BigFloat`); the
result is BITWISE independent of it.

`b` must carry any [`operator_on_zero`](@ref) correction (LINEARIZED-only, as for
[`smooth!`](@ref)). Called with `u = 0` this is the preconditioner apply — which is what
the callable form `mg(z, v)` does.
"""
function vcycle!(u::AbstractVector, mg::MultigridV, b::AbstractVector;
                 threaded::Bool = _mg_default_threaded(mg, u, b))
    n = mg.hier.levels[1].mesh.ndofs
    length(u) == n ||
        throw(ArgumentError("vcycle!: u has length $(length(u)), the finest mesh has $n dofs"))
    length(b) == n ||
        throw(ArgumentError("vcycle!: b has length $(length(b)), the finest mesh has $n dofs"))
    copyto!(mg.x[1], u)
    copyto!(mg.b[1], b)
    for _ in 1:mg.iterations
        _mg_vcycle!(mg, 1, threaded)
    end
    copyto!(u, mg.x[1])
    return u
end

"""
    (mg::MultigridV)(z, v) -> z

The right-preconditioner adapter for the `Minv!` hook of [`gmres`](@ref): one V-cycle
(or `mg.iterations` of them) started from `z = 0`, i.e. `z = M⁻¹ v` for the fixed linear
multigrid operator `M⁻¹`. Same contract as [`schwarz_preconditioner`](@ref)'s closure.
"""
function (mg::MultigridV)(z::AbstractVector, v::AbstractVector)
    fill!(z, zero(eltype(z)))
    return vcycle!(z, mg, v)
end

function _mg_vcycle!(mg::MultigridV{H,B,T}, ℓ::Int, threaded::Bool) where {H,B,T}
    lev = mg.hier.levels[ℓ]
    x = mg.x[ℓ]
    b = mg.b[ℓ]
    if ℓ == length(mg.hier.levels)
        if mg.bottom === nothing
            smooth!(x, lev.smoother, b; passes = mg.passes, threaded = threaded)
        else
            mg.bottom(x, b)
        end
        mg.post_smoothing_at_bottom &&
            smooth!(x, lev.smoother, b; passes = mg.passes, threaded = threaded)
        return nothing
    end
    mg.pre_smoothing &&
        smooth!(x, lev.smoother, b; passes = mg.passes, threaded = threaded)
    r = mg.r[ℓ]
    apply!(r, lev.op, x; threaded = threaded)
    @. r = b - r
    mg_restrict_residual!(mg.b[ℓ+1], mg.hier, ℓ, r)
    fill!(mg.x[ℓ+1], zero(T))
    _mg_vcycle!(mg, ℓ + 1, threaded)
    mg_prolong!(r, mg.hier, ℓ, mg.x[ℓ+1])
    @. x = x + r
    smooth!(x, lev.smoother, b; passes = mg.passes, threaded = threaded)
    return nothing
end

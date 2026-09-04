# amr.jl — h-split (halve element in both/either dims) + p-increment mesh
# transformations, 2:1 enforcement by splitting offenders, Lagrange interpolation of
# the current solution onto the new mesh; criteria are caller-supplied predicates.
# DAG nodes: amr-h-refinement, amr-p-refinement. Built TDD — no code lands before
# their failing tests.
#
# The CRP hp-AMR strategy (arXiv:2606.04998 main.tex:552-560, 714): each AMR step splits
# every element containing the puncture in half in BOTH dimensions (h) while raising the
# polynomial order of every OTHER element by one (p, from the initial N = 4), and restores
# the two-to-one balance at element interfaces by splitting the offending elements as well
# ("neighboring elements can differ by at most one h refinement level in any dimension",
# main.tex:554). The solution of the previous AMR level is carried over to the new mesh so
# each solve starts from a good initial guess (main.tex:714, fig:linsolv).
#
# Scope of this module (binding): mesh TRANSFORMATIONS and solution CARRY-OVER only.
#   * the refinement CRITERION is caller-supplied (`flags`, see [`h_refine`](@ref)) — the
#     CRP "element contains the puncture" predicate lives in ScalarSelfForce;
#   * the AMR LOOP (solve → flag → refine → carry over → re-solve) is the `amr-driver`
#     node, which composes [`h_split`](@ref), [`balance_two_to_one`](@ref),
#     [`balanced_flags`](@ref), [`h_refine`](@ref), [`p_refine`](@ref), [`hp_refine`](@ref)
#     and [`transfer_solution`](@ref);
#   * coarsening / element JOIN is NOT implemented — the CRP strategy never coarsens, and
#     `transfer_solution` REFUSES a target element with no ancestor in the source mesh
#     rather than silently producing a projection.
#
# Numerical contracts (test_amr.jl gates):
#   * every transformation is a PURE function: the source `Mesh2D` and the source solution
#     vector are never mutated, and the result shares no mutable state with the source.
#   * transformations return a mesh built through the PUBLIC `Mesh2D` constructor, so the
#     full invariant set (within-block 2:1 balance, face tiling, cross-block transverse
#     balance, layout/offsets/index consistency) is re-validated on every AMR step.
#   * the carry-over of one dyadic descendant is the parent's Lagrange polynomial evaluated
#     at the child's Gauss nodes: with the child's node ξ and the descendant's dyadic
#     address (Δ levels below the ancestor, sub-index j), the parent-frame point is
#     (ξ + c)/2^Δ with the ODD integer c = 2j + 1 − 2^Δ (c = 0 for Δ = 0). The SCALING is
#     exact (an `ldexp` by a power of two); the point itself carries at most ONE rounding,
#     from the integer shift `ξ + c` (none at Δ = 0, where c = 0). At Δ = 1 the form reduces
#     BITWISE to the sealed mortar half-interval prolongation
#     `projection_prolong(T, N, N, MORTAR_LOWER/UPPER)` (gated) — the same interpolation
#     semantics, not a second implementation.
#   * consequently the transfer is EXACT (round-off) on every polynomial the source element
#     represents (degree ≤ N_source − 1), and on smooth data the error is the source
#     element's Gauss interpolation error, ≤ H^N 2^N N!/(2N)! max|f^(N)| per dimension.
#   * geometry type and precision are inherited: a `Mesh2D{BigFloat}` carries a BigFloat
#     transfer at the ambient `setprecision`; a ComplexF64 solution over Float64 geometry
#     stays ComplexF64.
#
# 2:1 CLOSURE POLICY (decided from the CRP's own solver input, not by taste). The paper's
# ancillary spectre input file `anc/ScalarSelfForce.yaml:145` sets `Amr: Policies: Isotropy:
# Anisotropic`, and spectre defines that enum as "each dimension can be refined
# independently" (ParallelAlgorithms/Amr/Policies/Isotropy.hpp:24) — the policy is applied in
# `amr::enforce_policies` (EnforcePolicies.cpp:24-26), where the ONLY thing `Isotropic` does
# is equalize the per-dimension flags to their maximum. This module is therefore ANISOTROPIC
# by construction: `_amr_ask!` marks ONLY the offending dimension, and an element is quartered
# only when it genuinely violates the balance in both dimensions.
#
# Whether the requested dimensions are applied together or one at a time is a SEPARATE
# question, and the reference implementation settles it: `amr::update_amr_decision`
# (Domain/Amr/UpdateAmrDecision.cpp) raises per-dimension flags to a FIXPOINT over the
# PRE-SPLIT neighbour decisions, after which the domain is adjusted once — so an element whose
# decision ends up Split in both dimensions is quartered in a single step. The sweep loop
# below is that algorithm. A one-split-at-a-time closure over the same violation predicate is
# equally `Mesh2D`-legal and yields FEWER elements (measured 0-16.7 % fewer, worst 1.167×, on
# the off-centre two-block cascade of test_amr.jl gate (b2): step 4 raw 19 → 28 batched vs 24
# minimal), because splitting one dimension first can shrink an offender out of the second
# dimension's overlap. That saving is NOT taken here: it would deviate from the solver the
# CRP actually ran. The inflation is gated against a structurally independent minimal-closure
# oracle in test_amr.jl gate (b3).
#
# COST (measured, single thread, uniform meshes): `balance_two_to_one` 0.010 / 0.106 / 1.615 s
# and `h_refine` 0.37 / 0.85 / 11.1 s at n = 1024 / 4096 / 16384 elements — O(n²) from the
# pairwise adjacency scan (and `Mesh2D`'s own O(n²) face topology inside `h_refine`). The CRP
# runs a few hundred elements per m-mode (main.tex:552-560, 8 AMR iterations), so this is
# irrelevant at the paper's sizes; a driver pushing past ~10⁴ elements wants a spatial index.
#
# CALLER-MANAGED CACHE INVALIDATION (binding, for the `amr-driver` node): every AMR step
# produces a NEW `Mesh2D`, so every object built ON a mesh is stale — the `MeshOperator`,
# its element matrices/coefficients/mortars, and any `SchwarzSmoother` over it must be
# REBUILT against the new mesh. The Schwarz subdomain/LU cache is keyed by
# (element id, subdomain shape, epoch, shift) and its `epoch` is CALLER-managed (wave-6
# contract): a refinement that changes an element's order or its neighbours' shapes must
# bump the epoch, or a stale factorization is reused for an id whose subdomain changed.
# Element ids that survive a step unchanged (same id, same order, same neighbour shapes)
# may legitimately keep their cached factorization — that reuse is the CRP caching claim of
# main.tex:716 and is the driver's decision, not this module's.

export h_split, balance_two_to_one, balanced_flags, h_refine, p_refine, hp_refine,
       transfer_solution

# ---- refinement flags ---------------------------------------------------------------------
# `flags` is the caller-supplied criterion. Accepted forms, all reduced to a per-dimension
# `NTuple{2,Bool}` by `_amr_dims`: a Bool (both dimensions), an NTuple{2,Bool} (per
# dimension), and as containers a Set/Vector of flagged ElementIds, a Dict ElementId => flag
# (absent ⇒ false), or any callable predicate `f(::ElementId)`.

_amr_dims(b::Bool) = (b, b)
_amr_dims(t::NTuple{2,Bool}) = t
_amr_dims(x) = throw(ArgumentError("AMR flags must be a Bool or an NTuple{2,Bool} per " *
                                   "element, got $(typeof(x))"))

_amr_flag(flags::AbstractDict, id::ElementId) = _amr_dims(get(flags, id, false))
_amr_flag(flags::AbstractSet, id::ElementId) = _amr_dims(id in flags)
_amr_flag(flags::AbstractVector{ElementId}, id::ElementId) = _amr_dims(id in flags)
# a Bool MASK aligned with `mesh.elements` is a plausible mistake and would otherwise die with
# an opaque MethodError on `flags(id)`; name the accepted forms instead.
_amr_flag(::AbstractVector{Bool}, ::ElementId) =
    throw(ArgumentError("AMR flags: a Vector{Bool} mask is not accepted — element order is " *
                        "not part of the criterion contract. Pass a callable predicate " *
                        "`f(::ElementId)`, an AbstractDict `ElementId => flag`, or an " *
                        "AbstractSet/AbstractVector of flagged `ElementId`s"))
_amr_flag(flags, id::ElementId) = _amr_dims(flags(id))

# children of one element along the flagged dimensions, appended to `out` in the canonical
# order (dimension-1 index fastest — the Mesh2D layout convention)
function _amr_children!(out::Vector{ElementId}, id::ElementId, dims::NTuple{2,Bool})
    if dims[1] && dims[2]
        for jy in 0:1, jx in 0:1
            push!(out, child_id(child_id(id, 1, jx), 2, jy))
        end
    elseif dims[1]
        for jx in 0:1
            push!(out, child_id(id, 1, jx))
        end
    elseif dims[2]
        for jy in 0:1
            push!(out, child_id(id, 2, jy))
        end
    else
        push!(out, id)
    end
    return out
end

"""
    h_split(ids::AbstractVector{ElementId}, flags) -> Vector{ElementId}

Halve every flagged element of the leaf set `ids` along its flagged dimensions
(CRP main.tex:552) — pure refinement-tree address arithmetic over
[`child_id`](@ref), with NO 2:1 enforcement and no mesh. The result generally VIOLATES the
two-to-one balance; [`balance_two_to_one`](@ref) restores it and [`h_refine`](@ref)
composes both into a mesh transformation.

`flags` is the caller-supplied refinement criterion (the CRP "element contains the
puncture" predicate lives in the application package). Accepted forms, per element:

- a callable predicate `flags(id)`,
- an `AbstractDict` `id => flag` (an absent key is `false`),
- an `AbstractSet{ElementId}` or `AbstractVector{ElementId}` of flagged ids,

each yielding either a `Bool` (halve in BOTH dimensions — the CRP strategy) or an
`NTuple{2,Bool}` (halve in the listed dimensions only).

Order is deterministic and locality-preserving: each element is replaced in place by its
children with the dimension-1 index fastest, so the global nodal layout of the refined mesh
follows the layout of the source mesh.

STALE IDS ARE A SILENT NO-OP. A container form is consulted by element identity, and every
refinement step REPLACES the ids it splits, so a `Set`/`Dict`/`Vector` built against an
earlier mesh simply matches nothing and the step becomes an expensive identity. A driver must
re-evaluate its criterion against the CURRENT `mesh.elements` on every step and never cache a
flag container across steps (a callable predicate is immune, which is why it is the
recommended form).
"""
function h_split(ids::AbstractVector{ElementId}, flags)
    out = ElementId[]
    sizehint!(out, length(ids))
    for id in ids
        _amr_children!(out, id, _amr_flag(flags, id))
    end
    return out
end

# ---- 2:1 balance enforcement -----------------------------------------------------------------

# record "element `id` must be halved along dimension `d`" for the current sweep
function _amr_ask!(req::Dict{ElementId,NTuple{2,Bool}}, id::ElementId, d::Int)
    cur = get(req, id, (false, false))
    req[id] = d == 1 ? (true, cur[2]) : (cur[1], true)
    return nothing
end

"""
    balance_two_to_one(ids; connections = BlockConnection[], maxsweeps = 128)
        -> Vector{ElementId}

Restore the CRP two-to-one balance on the leaf set `ids` by splitting the offending
elements (main.tex:554: "in cases where the two-to-one balance would be violated, the
offending elements are split in half as well to restore the two-to-one balance").

Two constraints are enforced, matching exactly what [`Mesh2D`](@ref) validates:

- WITHIN a block, every face-adjacent pair (touching along one dimension and overlapping
  with positive measure in the other — the [`two_to_one_balanced`](@ref) adjacency) must
  differ by at most one refinement level in EACH dimension; the coarser element of a
  violating pair is halved along the offending dimension.
- ACROSS a [`BlockConnection`](@ref), the elements that meet at the seam must differ by at
  most one level in the TRANSVERSE dimension (the mortar representability condition
  validated by `Mesh2D`); levels along the connected dimension are not comparable across
  blocks (different block widths), so they are left alone.

The enforcement is a FIXPOINT: each sweep collects all violations, applies every split at
once, and repeats until a sweep finds none — so the result is independent of the pair
visitation order and the function is IDEMPOTENT (a second call returns the identical
vector, gated). It terminates because a split never raises an element above the level of
the neighbour that forced it, so the maximum level is invariant; `maxsweeps` is a defensive
cap (measured: ≤ 3 sweeps over 150 random multi-block layouts) and an exceeded cap is an
`ErrorException`, never a silent unbalanced return.

ANISOTROPIC, BATCHED, AND NOT MINIMAL — the CRP's own policy, see the module header for the
sources. Only the OFFENDING dimension of a violating pair is marked (`Isotropy: Anisotropic`,
`anc/ScalarSelfForce.yaml:145`), so an element is quartered only when it genuinely violates in
both dimensions; and the dimensions marked in one sweep are applied TOGETHER, which is what
spectre's `amr::update_amr_decision` fixpoint over pre-split neighbour decisions plus a single
domain adjustment does. A one-split-at-a-time closure over the same violation predicate is
equally valid and gives 0-16.7 % fewer elements (gated against an independent oracle in
test_amr.jl gate (b3)); it is deliberately NOT used, because it is not the algorithm the
paper's solver ran.

Element order is preserved: offenders are replaced in place by their children.
"""
function balance_two_to_one(ids::AbstractVector{ElementId};
                            connections::AbstractVector{BlockConnection} = BlockConnection[],
                            maxsweeps::Integer = 128)
    cur = collect(ids)
    for _ in 1:maxsweeps
        req = Dict{ElementId,NTuple{2,Bool}}()
        n = length(cur)
        @inbounds for i in 1:n-1
            a = cur[i]
            for j in i+1:n
                b = cur[j]
                a.block == b.block || continue
                t1, o1 = _touch_overlap(a.segments[1], b.segments[1])
                t2, o2 = _touch_overlap(a.segments[2], b.segments[2])
                (t1 && o2) || (t2 && o1) || continue         # not face-adjacent
                for d in 1:2
                    δ = a.segments[d].level - b.segments[d].level
                    δ >= 2 && _amr_ask!(req, b, d)
                    δ <= -2 && _amr_ask!(req, a, d)
                end
            end
        end
        for c in connections                                  # cross-block: transverse only
            d = c.dim
            t = 3 - d
            for a in cur
                a.block == c.block || continue
                sa = a.segments[d]
                sa.index == (Int64(1) << sa.level) - 1 || continue    # at the hi seam face
                for b in cur
                    b.block == c.neighbor || continue
                    b.segments[d].index == 0 || continue              # at the lo seam face
                    _mf_overlap(a.segments[t], b.segments[t]) || continue
                    δ = a.segments[t].level - b.segments[t].level
                    δ >= 2 && _amr_ask!(req, b, t)
                    δ <= -2 && _amr_ask!(req, a, t)
                end
            end
        end
        isempty(req) && return cur
        nxt = ElementId[]
        sizehint!(nxt, length(cur) + 3 * length(req))
        for id in cur
            _amr_children!(nxt, id, get(req, id, (false, false)))
        end
        cur = nxt
    end
    throw(ErrorException("balance_two_to_one: the 2:1 fixpoint was not reached in " *
                         "$maxsweeps sweeps — this is a bug (the maximum refinement " *
                         "level is invariant under enforcement, so the sweep count is " *
                         "bounded by it)"))
end

"""
    balanced_flags(mesh::Mesh2D, flags; maxsweeps = 128) -> Dict{ElementId,NTuple{2,Bool}}

The caller's refinement `flags` CLOSED under the two-to-one rule: the per-dimension decision
each element of `mesh` actually ends up with once [`balance_two_to_one`](@ref) has forced the
offenders to split as well. Absent keys mean "no split in either dimension".

# WHY THE CLOSURE IS A SEPARATE OBJECT (normative)

Because a balance-induced split is a SPLIT, not a p refinement, and the two are mutually
exclusive per dimension. spectre keeps ONE `amr::Flag` per element per dimension:
`amr::update_amr_decision` REPLACES that flag by `Flag::Split` when a neighbour's decision
would break the 2:1 rule (`Domain/Amr/UpdateAmrDecision.cpp:82-87`), and the mesh projector
grows the extents only where the FINAL flag is still `Flag::IncreaseResolution`
(`ParallelAlgorithms/Amr/Projectors/Mesh.cpp:30-34`). An element the balancer splits therefore
KEEPS its polynomial order in the split dimension, and [`hp_refine`](@ref) needs this closure —
not the caller's raw flags — to hold the right orders.

The closure is ANISOTROPIC, matching `Amr.Policies.Isotropy: Anisotropic` and
`EnforceTwoToOneBalanceInNormalDirection: True` (`anc/ScalarSelfForce.yaml:145-146`): only the
OFFENDING dimension is raised, so an element quartered by the balancer is one that genuinely
violates the rule in both dimensions. [`balance_two_to_one`](@ref) marks exactly that
dimension (`_amr_ask!`), so this function inherits the policy rather than restating it.

SCOPE — THIS IS FAITHFUL FOR THE CRP POLICY, NOT A GENERAL TRANSCRIPTION of
`amr::update_amr_decision`. It implements the refinement-only, normal-direction-enabled case
the paper runs: there is no `Flag::Join` branch (the CRP never coarsens), no
difference-of-three or potential-sibling reasoning, and no
`EnforceTwoToOneBalanceInNormalDirection = false` path. A configuration asking for any of
those is refused upstream by the criterion's own policy check rather than silently
approximated here.

# HOW IT IS COMPUTED

A fixpoint over the SPLIT LEAF SET, which is the same algorithm [`balance_two_to_one`](@ref)
already runs: apply the current decisions with [`h_split`](@ref), balance the result, and read
each leaf's level rise against its ancestor in `mesh.elements` back as a per-dimension
decision. Adding decisions only adds splits, so the iteration is monotone and terminates;
`maxsweeps` is a defensive cap. The result is a FIXPOINT — `balanced_flags(mesh, closed)`
returns `closed` (gated) — and it is idempotent under composition with `h_refine`, whose own
balancing pass then has nothing left to do except any SECOND split of an element the balancer
must halve twice in one step (the decision container holds one bit per dimension; the extra
split is still applied, and its children still inherit the parent's held order).

`flags` accepts every form [`h_split`](@ref) documents.

COST: the 2:1 fixpoint is run TWICE per step — once here to obtain the closure and once inside
[`h_refine`](@ref), which balances the split leaf set again. The second pass is a no-op except
for the rare element the balancer must halve twice in one step, so the redundancy is correct
and cheap at the CRP's few-hundred-element meshes; a driver pushing past ~10⁴ elements should
hand `h_refine` the closure's own leaf set instead.
"""
function balanced_flags(mesh::Mesh2D, flags; maxsweeps::Integer = 128)
    ids = mesh.elements
    dec = Dict{ElementId,NTuple{2,Bool}}()
    for id in ids
        d = _amr_flag(flags, id)
        (d[1] || d[2]) && (dec[id] = d)
    end
    lookup = Dict{ElementId,Bool}(id => true for id in ids)
    for _ in 1:maxsweeps
        leaves = balance_two_to_one(h_split(ids, dec); connections = mesh.connections)
        changed = false
        for leaf in leaves
            anc = _amr_ancestor(lookup, leaf)
            anc === nothing && continue
            cur = get(dec, anc, (false, false))
            nxt = (cur[1] || leaf.segments[1].level > anc.segments[1].level,
                   cur[2] || leaf.segments[2].level > anc.segments[2].level)
            if nxt != cur
                dec[anc] = nxt
                changed = true
            end
        end
        changed || return dec
    end
    throw(ErrorException("balanced_flags: the decision fixpoint was not reached in " *
                         "$maxsweeps sweeps — this is a bug (decisions only ever gain " *
                         "splits, so the sweep count is bounded by the element count)"))
end

# ---- ancestry ---------------------------------------------------------------------------------

# the unique ancestor-or-equal of `id` present in `lookup` (a Dict keyed by ElementId), or
# `nothing`. Candidates are the (level+1)² addresses obtained by reducing each dimension's
# level; the identity is tried first.
function _amr_ancestor(lookup, id::ElementId)
    s1 = id.segments[1]
    s2 = id.segments[2]
    for a in 0:s1.level, b in 0:s2.level
        cand = ElementId(id.block, SegmentId(s1.level - a, s1.index >> a),
                                   SegmentId(s2.level - b, s2.index >> b))
        haskey(lookup, cand) && return cand
    end
    return nothing
end

_amr_geometry(::Mesh2D{RT}) where {RT} = RT

# per-element orders of a mesh as a Dict (the `Mesh2D` constructor's `orders` form)
_amr_orderdict(mesh::Mesh2D) =
    Dict{ElementId,NTuple{2,Int}}(mesh.elements[e] => mesh.orders[e]
                                  for e in eachindex(mesh.elements))

# orders of a refined leaf set, inherited from each element's ancestor in `ord0`
function _amr_inherit(ord0::AbstractDict{ElementId,NTuple{2,Int}},
                      ids::AbstractVector{ElementId})
    ords = Dict{ElementId,NTuple{2,Int}}()
    for id in ids
        anc = _amr_ancestor(ord0, id)
        anc === nothing &&
            throw(ArgumentError("h_refine: refined element $id has no ancestor in the " *
                                "source mesh — the flag set produced an element outside " *
                                "the source leaf set"))
        ords[id] = ord0[anc]
    end
    return ords
end

"""
    h_refine(mesh::Mesh2D, flags) -> Mesh2D

One CRP h-refinement step (main.tex:552-554) as a PURE mesh transformation: halve every
flagged element ([`h_split`](@ref), see there for the accepted `flags` forms), restore the
two-to-one balance by splitting the offenders ([`balance_two_to_one`](@ref), with the
mesh's own [`BlockConnection`](@ref)s), and rebuild the mesh through the public
[`Mesh2D`](@ref) constructor — so every mesh invariant is re-validated on every AMR step.

Children INHERIT their parent's polynomial order (`h` refinement changes the element size,
not the order; the CRP raises the order only in the un-flagged elements, which is
[`p_refine`](@ref), and the combined step is [`hp_refine`](@ref) — which closes the flag set
with [`balanced_flags`](@ref) FIRST, so that an element the balancer splits keeps its order in
the split dimension instead of inheriting a p increment it never earned). Blocks and connections
are carried over unchanged. `mesh` is not mutated and shares no mutable state with the
result.

Carry the previous level's solution over with [`transfer_solution`](@ref); rebuild every
operator and bump the caller-managed Schwarz cache epoch (module header).
"""
function h_refine(mesh::Mesh2D, flags)
    ord0 = _amr_orderdict(mesh)
    ids = balance_two_to_one(h_split(mesh.elements, flags); connections = mesh.connections)
    return Mesh2D(mesh.blocks, ids, _amr_inherit(ord0, ids); connections = mesh.connections)
end


# ---- p refinement ---------------------------------------------------------------------------

"""
    p_refine(mesh::Mesh2D; hold = Set{ElementId}(), increment = 1) -> Mesh2D

One CRP p-refinement step (main.tex:558-560) as a PURE mesh transformation: "increase the
polynomial order in all elements *that do not contain the puncture* by one in both
dimensions in each refinement step", starting from the CRP's initial `4 × 4` grid. The
element set is UNCHANGED; only the per-element orders move, so after `k` steps with a fixed
flag set the un-flagged elements carry `N = 4 + k` and the flagged ones still `N = 4`.

`hold` names the elements whose order is HELD — the same criterion that
[`h_refine`](@ref) SPLITS (the CRP puncture elements are h-refined *instead of* p-refined),
so one criterion drives both halves of the step and [`hp_refine`](@ref) takes exactly that
one set. It is a KEYWORD, and named for what it does, precisely because the polarity is the
opposite of `h_refine`'s positional `flags` ("split these"): `p_refine(mesh; hold = puncture)`
cannot be misread the way `p_refine(mesh, puncture)` could. The accepted forms are the ones
documented on [`h_split`](@ref); a per-dimension `NTuple{2,Bool}` holds only the flagged
dimension's order, and the default holds nothing (a uniform p step).

The result is built through the public [`Mesh2D`](@ref) constructor, so p-nonconforming
interfaces are validated there and handled by the mortar layer of [`MeshOperator`](@ref);
`mesh` is not mutated. Carry the solution over with [`transfer_solution`](@ref) (exact on
every polynomial the source element represents) and rebuild every operator, bumping the
caller-managed Schwarz cache epoch (module header) — a p step changes every un-flagged
element's subdomain matrix.
"""
function p_refine(mesh::Mesh2D; hold = Set{ElementId}(), increment::Integer = 1)
    increment >= 0 ||
        throw(ArgumentError("p_refine: increment must be ≥ 0, got $increment"))
    ords = Dict{ElementId,NTuple{2,Int}}()
    for (e, id) in enumerate(mesh.elements)
        ords[id] = _amr_praise(mesh.orders[e], _amr_flag(hold, id), Int(increment))
    end
    return Mesh2D(mesh.blocks, mesh.elements, ords; connections = mesh.connections)
end

# per-dimension order increment: a flagged (held) dimension keeps its order
_amr_praise(o::NTuple{2,Int}, hold::NTuple{2,Bool}, inc::Int) =
    (hold[1] ? o[1] : o[1] + inc, hold[2] ? o[2] : o[2] + inc)

"""
    hp_refine(mesh::Mesh2D, flags; increment = 1) -> Mesh2D

One full CRP hp-AMR step (main.tex:552-560) from a SINGLE caller-supplied criterion: the
flagged elements are halved ([`h_refine`](@ref)) and keep their polynomial order, every
other element keeps its extent and gains `increment` orders ([`p_refine`](@ref)), and the
two-to-one balance is restored by splitting the offenders — whose children inherit their
parent's already-incremented order.

LITERALLY the composition `h_refine(p_refine(mesh; hold = closed, increment), closed)` with
`closed = `[`balanced_flags`](@ref)`(mesh, flags)` — ONE criterion, used as "hold the order"
by the p half and "split" by the h half. It is also equal to the other order once the closed
flags are lifted along ancestry (gated): the two halves act on disjoint element sets, so they
commute. This is the step the `amr-driver` node iterates; the criterion behind `flags` stays
the caller's.

THE CLOSURE IS NOT OPTIONAL. An element the 2:1 rule forces to split must NOT also gain a
polynomial order in the split dimension: spectre replaces the per-dimension flag by
`Flag::Split` (`Domain/Amr/UpdateAmrDecision.cpp:82-87`) and grows extents only where the FINAL
flag is `IncreaseResolution` (`ParallelAlgorithms/Amr/Projectors/Mesh.cpp:30-34`). This follows
the reference IMPLEMENTATION where the paper's prose does not settle the case: main.tex:554's
"increase the polynomial order in all elements that do not contain the puncture" would, read
literally, include the elements the balance rule splits. Passing the
caller's RAW flags to the p half instead — as this function did before — p-increments the
offenders and their children inherit that order, which moves the order map, the dof count, the
Schwarz cache keys and the iteration counts away from the reference. The decision stays
per-dimension throughout, so an anisotropic balance split holds the order in the split
dimension and still increments the other one.
"""
function hp_refine(mesh::Mesh2D, flags; increment::Integer = 1)
    closed = balanced_flags(mesh, flags)
    return h_refine(p_refine(mesh; hold = closed, increment = increment), closed)
end

# ---- solution carry-over ------------------------------------------------------------------------

# 1D nodal transfer matrix (N_new × N_old) from an ancestor element to a descendant that
# sits Δ levels below it at sub-index j: the ancestor's Lagrange basis evaluated at the
# descendant's Gauss nodes expressed in the ancestor's reference frame,
#   ξ_ancestor = (ξ_child + c) / 2^Δ,   c = 2j + 1 − 2^Δ   (c = 0 when Δ = 0).
# `ldexp` makes the SCALING exact (a power of two); the point carries at most one rounding
# from `ξ + c` (none at Δ = 0). At Δ = 1 this is BITWISE the sealed mortar prolongation.
function _amr_xfer1d(::Type{RT}, Nold::Int, Nnew::Int, Δ::Int, j::Int,
                     cache::Dict{NTuple{4,Int},Matrix{RT}}) where {RT<:AbstractFloat}
    return get!(cache, (Nold, Nnew, Δ, j)) do
        c = 2 * j + 1 - (Int64(1) << Δ)
        cRT = RT(c)
        pts = [ldexp(ξ + cRT, -Δ) for ξ in refelement(RT, Nnew).nodes]
        return interpolation_matrix(refelement(RT, Nold), pts)
    end
end

"""
    transfer_solution(source::Mesh2D, target::Mesh2D, u::AbstractVector) -> Vector

Carry the nodal solution `u` (in `source`'s global layout) over to `target`'s layout by
Lagrange interpolation (CRP main.tex:714 — the previous AMR level's solution is the initial
guess for the next). Every `target` element must be the SAME element as, or a dyadic
descendant of, exactly one `source` element; its values are the source element's
tensor-product Lagrange polynomial evaluated at the target element's Gauss nodes,

    v = P_ξ · U · P_ηᵀ,    P_d = ℓ^source( (ξ_target + c_d)/2^{Δ_d} ),

built from the sealed [`interpolation_matrix`](@ref) of the source element's
[`RefElement`](@ref). Consequences (gated):

- EXACT (to round-off) whenever the source element represents the field exactly, i.e. for
  every polynomial of degree ≤ `N_source − 1` per dimension — h refinement, p increase, and
  the two combined;
- on smooth data the error is the SOURCE element's Gauss interpolation error, bounded by
  `H^N 2^N N!/(2N)! max|f^(N)|` per dimension (`H` the source element's half-width) and
  decaying by `2^-N` per refinement level;
- the identity map (`source === target` in content) is BITWISE the input;
- one-level, same-order h refinement is BITWISE
  `projection_prolong(T, N, N, MORTAR_LOWER/UPPER)` — the same interpolation as the mortar
  layer, not a second implementation.

The element type follows `promote_type(eltype(u), RT)` with `RT` the source mesh's geometry
type, so a ComplexF64 solution over Float64 geometry stays ComplexF64 and a
`Mesh2D{BigFloat}` transfers in BigFloat at the ambient precision.

Coarsening / element JOIN is out of scope (the CRP strategy never coarsens): a `target`
element with no ancestor in `source` is an `ArgumentError`, never a silent projection.
"""
function transfer_solution(source::Mesh2D, target::Mesh2D, u::AbstractVector)
    length(u) == source.ndofs ||
        throw(ArgumentError("transfer_solution: u has length $(length(u)), the source " *
                            "mesh has $(source.ndofs) dofs"))
    RT = _amr_geometry(source)
    S = promote_type(eltype(u), RT)
    v = Vector{S}(undef, target.ndofs)
    cache = Dict{NTuple{4,Int},Matrix{RT}}()
    for (e, id) in enumerate(target.elements)
        anc = _amr_ancestor(source.index, id)
        anc === nothing &&
            throw(ArgumentError("transfer_solution: target element $id has no ancestor " *
                                "in the source mesh — coarsening/element join is not " *
                                "supported (the CRP hp-AMR strategy never coarsens)"))
        eo = source.index[anc]
        No = source.orders[eo]
        Nn = target.orders[e]
        all(Nn .>= No) ||
            throw(ArgumentError("transfer_solution: target element $id has order $Nn " *
                                "below its source element's $No — an order DECREASE " *
                                "discards the modes above the target degree, so nodal " *
                                "interpolation is NOT the L2 projection of the source " *
                                "field (deviation is O(1) on data that uses the dropped " *
                                "modes). p-coarsening needs `projection_restrict`, not " *
                                "this function"))
        P1, P2 = ntuple(2) do d
            so = anc.segments[d]
            sn = id.segments[d]
            Δ = sn.level - so.level
            _amr_xfer1d(RT, No[d], Nn[d], Δ, sn.index - (so.index << Δ), cache)
        end
        U = reshape(view(u, source.offsets[eo]+1:source.offsets[eo]+prod(No)), No[1], No[2])
        copyto!(view(v, target.offsets[e]+1:target.offsets[e]+prod(Nn)),
                vec(P1 * U * transpose(P2)))
    end
    return v
end

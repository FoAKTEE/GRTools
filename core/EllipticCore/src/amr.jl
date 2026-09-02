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
#     [`h_refine`](@ref), [`p_refine`](@ref), [`hp_refine`](@ref) and
#     [`transfer_solution`](@ref);
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
#     (ξ + c)/2^Δ with the ODD integer c = 2j + 1 − 2^Δ (c = 0 for Δ = 0). This form is
#     exact in binary FP: the scaling is an `ldexp`, and at Δ = 1 it reduces BITWISE to the
#     sealed mortar half-interval prolongation `projection_prolong(T, N, N, MORTAR_LOWER/
#     UPPER)` (gated) — the same interpolation semantics, not a second implementation.
#   * consequently the transfer is EXACT (round-off) on every polynomial the source element
#     represents (degree ≤ N_source − 1), and on smooth data the error is the source
#     element's Gauss interpolation error, ≤ H^N 2^N N!/(2N)! max|f^(N)| per dimension.
#   * geometry type and precision are inherited: a `Mesh2D{BigFloat}` carries a BigFloat
#     transfer at the ambient `setprecision`; a ComplexF64 solution over Float64 geometry
#     stays ComplexF64.
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

export h_split, balance_two_to_one, h_refine, p_refine, hp_refine, transfer_solution

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
[`p_refine`](@ref), and the combined step is [`hp_refine`](@ref)). Blocks and connections
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
    p_refine(mesh::Mesh2D, flags; increment = 1) -> Mesh2D

One CRP p-refinement step (main.tex:558-560) as a PURE mesh transformation: "increase the
polynomial order in all elements *that do not contain the puncture* by one in both
dimensions in each refinement step", starting from the CRP's initial `4 × 4` grid. The
element set is UNCHANGED; only the per-element orders move, so after `k` steps with a fixed
flag set the un-flagged elements carry `N = 4 + k` and the flagged ones still `N = 4`.

`flags` marks the elements whose order is HELD — the same criterion that
[`h_refine`](@ref) splits (the CRP puncture elements are h-refined *instead of* p-refined),
so one flag set drives both halves of the step and [`hp_refine`](@ref) takes exactly that
one set. The accepted flag forms are the ones documented on [`h_split`](@ref); a
per-dimension `NTuple{2,Bool}` holds only the flagged dimension's order.

The result is built through the public [`Mesh2D`](@ref) constructor, so p-nonconforming
interfaces are validated there and handled by the mortar layer of [`MeshOperator`](@ref);
`mesh` is not mutated. Carry the solution over with [`transfer_solution`](@ref) (exact on
every polynomial the source element represents) and rebuild every operator, bumping the
caller-managed Schwarz cache epoch (module header) — a p step changes every un-flagged
element's subdomain matrix.
"""
function p_refine(mesh::Mesh2D, flags; increment::Integer = 1)
    increment >= 0 ||
        throw(ArgumentError("p_refine: increment must be ≥ 0, got $increment"))
    ords = Dict{ElementId,NTuple{2,Int}}()
    for (e, id) in enumerate(mesh.elements)
        ords[id] = _amr_praise(mesh.orders[e], _amr_flag(flags, id), Int(increment))
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

Exactly the composition `h_refine(p_refine(mesh, flags; increment), flags)` (gated), and
equally the other order once the flags are lifted along ancestry — the two halves of the
step act on disjoint element sets, so they commute. This is the step the `amr-driver` node
iterates; the criterion behind `flags` stays the caller's.
"""
function hp_refine(mesh::Mesh2D, flags; increment::Integer = 1)
    increment >= 0 ||
        throw(ArgumentError("hp_refine: increment must be ≥ 0, got $increment"))
    ord0 = Dict{ElementId,NTuple{2,Int}}()
    for (e, id) in enumerate(mesh.elements)
        ord0[id] = _amr_praise(mesh.orders[e], _amr_flag(flags, id), Int(increment))
    end
    ids = balance_two_to_one(h_split(mesh.elements, flags); connections = mesh.connections)
    return Mesh2D(mesh.blocks, ids, _amr_inherit(ord0, ids); connections = mesh.connections)
end

# ---- solution carry-over ------------------------------------------------------------------------

# 1D nodal transfer matrix (N_new × N_old) from an ancestor element to a descendant that
# sits Δ levels below it at sub-index j: the ancestor's Lagrange basis evaluated at the
# descendant's Gauss nodes expressed in the ancestor's reference frame,
#   ξ_ancestor = (ξ_child + c) / 2^Δ,   c = 2j + 1 − 2^Δ   (c = 0 when Δ = 0).
# `ldexp` makes the scaling exact; at Δ = 1 this is BITWISE the sealed mortar prolongation.
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

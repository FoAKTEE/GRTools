# matrixfree.jl — the OWNED global composition of the elliptic hp-DG operator (CRP
# main.tex:668 matrix-free discipline): Mesh2D-level connectivity containers and the
# two-pass sweep that waves 4-5 kept test-side, composing dg_operator (strong residual,
# conventions N1-N7 inherited VERBATIM) + dg_penalty (eq:numflux, receiver-σ boundary
# ghosts) + mortars (nonconforming h/p coupling) + boundary_hooks (ModifyBoundaryData +
# ghost BCs) into a single matrix-free A(u) over a global nodal vector layout.
# DAG node: matrix-free-operator (elliptic-ssf GLOBAL_DAG row 15, TOPOLOGICAL_SCHEDULE
# W6-A, CRITICAL PATH). Built TDD — gates calibrated pre-red in scratchpad
# calibrate_matrixfree.jl against an independent quadrature-level dense assembly.
#
# The two-pass sweep (per apply):
#   pass 1 — every element builds its ElementState (derivatives, fluxes, own-face traces).
#   pass 2 — every element gathers its exterior data per face: the SENDER's face state is
#     projected to each interface mortar (project_to_mortar with the sender's own P —
#     "both neighbors project their interface data to the mortar", main.tex:678-686), the
#     receiver applies its ModifyBoundaryData hook to the RAW received mortar data (with
#     its own trace projected to the same mortar as the hook's `int`), packs the mortar
#     states in receiver-mortar order, and strong_residual! lifts the numerical flux.
#     Block-boundary faces synthesize their exterior from the ghost-BC seam instead.
# Both passes are element-parallel with NO cross-element reductions, so the threaded
# apply is BITWISE identical to serial (gated). BigFloat is never threaded (binding
# owner directive: BigFloat scales via processes, not threads).
#
# EVERY interior interface goes through the mortar path — on conforming faces the mortar
# operators are exact identities and the flux is BITWISE PenaltyFlux (the W5-A gate), so
# a single code path serves conforming, p-jump, and 2:1 h-nonconforming interfaces while
# reproducing the wave-4 test-side sweep bitwise (gated in test_matrixfree.jl).
#
# W5-A pairwise-(σ, Am) contract (BINDING, inherited): both sides of an interface mortar
# hold IDENTICAL (σ, Am) — Am is evaluated ONCE per interface mortar (coefficient callable
# at the mortar's physical points, shared by both sides' Mortar objects) and σ comes from
# the mirrored Mortar constructor arguments (bitwise-equal by the max/min symmetry of
# penalty_sigma). The receiver-side-σ shortcut breaks SIPG symmetry at O(1) (W5-A
# falsifier); this composition makes it unconstructible through the public API.
#
# W5-B linearized-seam contract (BINDING, inherited): every MATRIX application — apply!,
# the GMRES operator callable, build_matrix column probes, direct_solve, and the later
# Schwarz/multigrid consumers — routes received data through modify_boundary_data_
# linearized and boundary ghosts through ghost_state_linearized. Fixed (trial-independent)
# data act only in apply_full! and reach the RHS via OPERATOR-ON-ZERO:
# b -= operator_on_zero(op) (spectre DgOperator.hpp:1081-1130 pattern).
#
# Cross-block closure (the mesh2d-blocks caveat, owned here): blocks are joined by
# BlockConnection (hi face of one block to the lo face of its neighbor in the same
# dimension) under two EXACT construction requirements — identical transverse coordinate
# maps (the transverse dyadic segment arithmetic then carries across the seam unchanged)
# and exact physical endpoint continuity along the connected dimension (map endpoints are
# exact by the mesh2d contract, so valid domains satisfy this bitwise). Cross-block
# balance is enforced on the TRANSVERSE dimension only (level difference ≤ 1, the mortar
# representability condition); normal-dimension levels are NOT comparable across blocks
# (different block widths — the pairwise σ handles the size contrast). One connection per
# block face, axis-aligned full-edge contact; richer layouts (partial-edge contact,
# multiple neighbor blocks per face) are the domain-blocks-map node's scope.
# Hooks bind per FACE and are applied per mortar (uniform on a split face) — per-mortar
# data granularity (distinct fixed data on the two halves of a 2:1 face) is the
# documented extension point for the dg-jump-conditions node.

export BlockConnection, Mesh2D, boundary_faces,
       MeshOperator, apply!, apply_full!, operator_on_zero, build_matrix, direct_solve

"""
    BlockConnection(block, dim, neighbor)

Axis-aligned inter-block face identification: the `dim`-upper (`ξ_dim = +1`) face of block
`block` abuts the `dim`-lower face of block `neighbor` along its full edge. Requirements
enforced at [`Mesh2D`](@ref) construction: both blocks share the IDENTICAL transverse
coordinate map (same type, equal parameters — the transverse dyadic segment arithmetic is
then valid across the seam), and the physical endpoints match EXACTLY
(`map_forward(blockmap(a, dim), +1) == map_forward(blockmap(b, dim), -1)`; endpoints are
exact by the mesh2d contract). One connection per block face.
"""
struct BlockConnection
    block::Int
    dim::Int
    neighbor::Int
    function BlockConnection(block::Integer, dim::Integer, neighbor::Integer)
        (dim == 1 || dim == 2) ||
            throw(ArgumentError("BlockConnection: dim must be 1 or 2, got $dim"))
        block == neighbor &&
            throw(ArgumentError("BlockConnection: self-connections are not supported"))
        return new(Int(block), Int(dim), Int(neighbor))
    end
end

# One interior element face as seen from the RECEIVER: per mortar (transverse order,
# lower first) the neighbor's element index, this side's / the neighbor side's mortar
# coverage, the mortar's transverse dyadic segment (block-local logical coordinate,
# shared bitwise across a connection by the identical-transverse-map requirement), and
# the index of this interface in the NEIGHBOR face's mortar list.
struct MeshInterface
    nbs::Vector{Int}
    msize::Vector{MortarSize}
    nbmsize::Vector{MortarSize}
    span::Vector{SegmentId}
    nbslot::Vector{Int}
end

"""
    Mesh2D(blocks, elements, orders; connections = BlockConnection[])
    Mesh2D(block::Block, elements, orders; ...)

The multi-block 2D mesh connectivity container (deferred from mesh2d-blocks to its first
consumer): `blocks` (one or a vector of [`Block`](@ref)s, all with the same geometry
type), `elements` (a vector of leaf [`ElementId`](@ref)s — the order DEFINES the global
nodal vector layout), and `orders` — per-element collocation points `(Nx, Ny)`, either
one tuple for all elements or an `AbstractDict` `ElementId => (Nx, Ny)`.

Public fields (the documented layout contract):
- `elements::Vector{ElementId}`, `index::Dict{ElementId,Int}` — layout order and lookup.
- `orders::Vector{NTuple{2,Int}}`, `offsets::Vector{Int}`, `ndofs::Int` — element `e`
  (by position) owns global dofs `offsets[e]+1 : offsets[e]+prod(orders[e])`, i-fastest
  `tensor_lin` ordering within the element.
- `connections` — the cross-block seams ([`BlockConnection`](@ref)).

Construction validates and precomputes the FULL face topology: within-block 2:1 balance
(`two_to_one_balanced`, both dimensions), tiling (every interior face covered by exactly
one full or exactly two half neighbor faces; holes and >2:1 transverse jumps —
including across connections — are `ArgumentError`s), and the boundary-face list
([`boundary_faces`](@ref)). Face keys throughout are `(ElementId, d, sidx)` with
`d ∈ (1, 2)` and `sidx` 1 ↔ side −1, 2 ↔ side +1 (the wave-4 face indexing).
"""
struct Mesh2D{RT<:AbstractFloat}
    blocks::Vector{Block}
    blockindex::Dict{Int,Int}
    elements::Vector{ElementId}
    index::Dict{ElementId,Int}
    orders::Vector{NTuple{2,Int}}
    offsets::Vector{Int}
    ndofs::Int
    connections::Vector{BlockConnection}
    faces::Vector{NTuple{2,NTuple{2,Union{MeshInterface,Nothing}}}}
end

# segment endpoints in integer units of 2^-L (exact; the mesh2d SegmentId level cap)
_mf_ep(s::SegmentId, L::Int) =
    (Int64(s.index) << (L - s.level), Int64(s.index + 1) << (L - s.level))
function _mf_touchat(sa::SegmentId, sb::SegmentId, side::Int)
    L = max(sa.level, sb.level)
    alo, ahi = _mf_ep(sa, L)
    blo, bhi = _mf_ep(sb, L)
    return side == 1 ? ahi == blo : bhi == alo
end
function _mf_overlap(sa::SegmentId, sb::SegmentId)
    L = max(sa.level, sb.level)
    alo, ahi = _mf_ep(sa, L)
    blo, bhi = _mf_ep(sb, L)
    return max(alo, blo) < min(ahi, bhi)
end

# identical-map check for the cross-block transverse requirement: same concrete type and
# equal parameters (BigFloat maps compare by value, not object identity)
_mf_samemap(a::CoordinateMap1D, b::CoordinateMap1D) =
    typeof(a) == typeof(b) &&
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))

_mf_orders(o::Tuple, elements) = [(Int(o[1]), Int(o[2])) for _ in elements]
_mf_orders(o::AbstractDict, elements) =
    [(Int(o[el][1]), Int(o[el][2])) for el in elements]

function Mesh2D(blocks::AbstractVector, elements::AbstractVector{ElementId}, orders;
                connections::AbstractVector{BlockConnection} = BlockConnection[])
    isempty(blocks) && throw(ArgumentError("Mesh2D: needs at least one block"))
    RT = typeof(blocks[1]).parameters[1]
    all(typeof(b).parameters[1] == RT for b in blocks) ||
        throw(ArgumentError("Mesh2D: all blocks must share one geometry type"))
    blockindex = Dict{Int,Int}()
    for (k, b) in enumerate(blocks)
        haskey(blockindex, b.id) &&
            throw(ArgumentError("Mesh2D: duplicate block id $(b.id)"))
        blockindex[b.id] = k
    end
    index = Dict{ElementId,Int}()
    for (e, el) in enumerate(elements)
        haskey(index, el) && throw(ArgumentError("Mesh2D: duplicate element $el"))
        haskey(blockindex, el.block) ||
            throw(ArgumentError("Mesh2D: element $el references unknown block $(el.block)"))
        index[el] = e
    end
    ords = _mf_orders(orders, elements)
    offsets = zeros(Int, length(elements))
    for e in 2:length(elements)
        offsets[e] = offsets[e-1] + prod(ords[e-1])
    end
    ndofs = isempty(elements) ? 0 : offsets[end] + prod(ords[end])
    two_to_one_balanced(collect(elements)) ||
        throw(ArgumentError("Mesh2D: element set violates within-block 2:1 balance " *
                            "(CRP main.tex:554)"))
    conn = Dict{Tuple{Int,Int,Int},Int}()
    for c in connections
        (haskey(blockindex, c.block) && haskey(blockindex, c.neighbor)) ||
            throw(ArgumentError("Mesh2D: connection references unknown block"))
        (haskey(conn, (c.block, c.dim, +1)) || haskey(conn, (c.neighbor, c.dim, -1))) &&
            throw(ArgumentError("Mesh2D: duplicate connection on a block face"))
        a = blocks[blockindex[c.block]]
        b = blocks[blockindex[c.neighbor]]
        _mf_samemap(blockmap(a, 3 - c.dim), blockmap(b, 3 - c.dim)) ||
            throw(ArgumentError("Mesh2D: connected blocks $(c.block)→$(c.neighbor) must " *
                                "share an IDENTICAL transverse map (dim $(3 - c.dim))"))
        map_forward(blockmap(a, c.dim), one(RT)) == map_forward(blockmap(b, c.dim), -one(RT)) ||
            throw(ArgumentError("Mesh2D: connected blocks $(c.block)→$(c.neighbor) must " *
                                "meet EXACTLY along dim $(c.dim) (endpoints are exact " *
                                "by the mesh2d contract)"))
        conn[(c.block, c.dim, +1)] = c.neighbor
        conn[(c.neighbor, c.dim, -1)] = c.block
    end
    faces = [_mf_face_topology(elements, conn, e, d, sidx)
             for e in eachindex(elements), d in 1:2, sidx in 1:2]
    packed = [ntuple(d -> ntuple(sidx -> faces[e, d, sidx], 2), 2)
              for e in eachindex(elements)]
    return Mesh2D{RT}(Block[b for b in blocks], blockindex, collect(elements), index,
                      ords, offsets, ndofs, collect(connections), packed)
end
Mesh2D(block::Block, elements::AbstractVector{ElementId}, orders; kwargs...) =
    Mesh2D([block], elements, orders; kwargs...)

# classify one face of one element: nothing (boundary) or the MeshInterface
function _mf_face_topology(elements, conn, e, d, sidx)
    el = elements[e]
    s = sidx == 1 ? -1 : +1
    t = 3 - d
    cands = Int[]
    for (f, ot) in enumerate(elements)
        (f != e && ot.block == el.block) || continue
        _mf_touchat(el.segments[d], ot.segments[d], s) || continue
        _mf_overlap(el.segments[t], ot.segments[t]) || continue
        push!(cands, f)
    end
    if isempty(cands)
        seg = el.segments[d]
        atedge = s == -1 ? seg.index == 0 : seg.index == (Int64(1) << seg.level) - 1
        if atedge && haskey(conn, (el.block, d, s))
            nbblk = conn[(el.block, d, s)]
            for (f, ot) in enumerate(elements)
                ot.block == nbblk || continue
                os = ot.segments[d]
                facing = s == 1 ? os.index == 0 :
                                  os.index == (Int64(1) << os.level) - 1
                facing || continue
                _mf_overlap(el.segments[t], ot.segments[t]) || continue
                push!(cands, f)
            end
            isempty(cands) &&
                throw(ArgumentError("Mesh2D: tiling hole across the block connection at " *
                                    "$el face ($d, $s)"))
        elseif atedge
            return nothing                            # domain-boundary face
        else
            throw(ArgumentError("Mesh2D: tiling hole at $el face ($d, $s)"))
        end
    end
    myt = el.segments[t]
    rel = map(cands) do f
        nbt = elements[f].segments[t]
        if nbt.level == myt.level && nbt.index == myt.index
            (:equal, 0)
        elseif nbt.level == myt.level + 1 && (nbt.index >> 1) == myt.index
            (:child, nbt.index & 1)
        elseif nbt.level == myt.level - 1 && (myt.index >> 1) == nbt.index
            (:parent, myt.index & 1)
        else
            throw(ArgumentError("Mesh2D: transverse 2:1 violation at $el face ($d, $s)"))
        end
    end
    if length(cands) == 1 && rel[1][1] == :equal
        return MeshInterface([cands[1]], [MORTAR_FULL], [MORTAR_FULL], [myt], [1])
    elseif length(cands) == 1 && rel[1][1] == :parent
        which = rel[1][2]
        return MeshInterface([cands[1]], [MORTAR_FULL],
                             [which == 0 ? MORTAR_LOWER : MORTAR_UPPER], [myt], [which + 1])
    elseif length(cands) == 2 && rel[1][1] == :child && rel[2][1] == :child &&
           rel[1][2] != rel[2][2]
        ord = sortperm([r[2] for r in rel])
        return MeshInterface(cands[ord], [MORTAR_LOWER, MORTAR_UPPER],
                             [MORTAR_FULL, MORTAR_FULL],
                             [SegmentId(myt.level + 1, 2 * myt.index),
                              SegmentId(myt.level + 1, 2 * myt.index + 1)], [1, 1])
    end
    throw(ArgumentError("Mesh2D: invalid face coverage at $el face ($d, $s): " *
                        "$(length(cands)) neighbors"))
end

_element_range(mesh::Mesh2D, e::Int) =
    mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])
_element_block(mesh::Mesh2D, e::Int) = mesh.blocks[mesh.blockindex[mesh.elements[e].block]]

"""
    boundary_faces(mesh::Mesh2D) -> Vector{Tuple{ElementId,Int,Int}}

The domain-boundary faces of the mesh as `(element, d, sidx)` keys (sidx 1 ↔ side −1,
2 ↔ side +1), in element-layout order. Exactly these faces must receive a
[`GhostBoundary`](@ref) in the [`MeshOperator`](@ref) `bcs` table.
"""
function boundary_faces(mesh::Mesh2D)
    out = Tuple{ElementId,Int,Int}[]
    for e in eachindex(mesh.elements), d in 1:2, sidx in 1:2
        mesh.faces[e][d][sidx] === nothing && push!(out, (mesh.elements[e], d, sidx))
    end
    return out
end

"""
    MeshOperator(T, mesh::Mesh2D, A, β, γ; C = 1.5, hooks = Dict(), bcs) -> MeshOperator

The global matrix-free elliptic DG operator over `mesh`: the owned composition of the
wave-3/4/5 layers into a single `A(u)` callable on the global nodal layout. `T` is the
requested solver scalar type (as in [`element_matrices`](@ref)); the operator scalar `S`
promotes `T` with the evaluated coefficient types (complex β ⇒ complex `S`). `A = (𝔄¹¹,
𝔄²²)`, `β`, `γ = (γ¹, γ²)` are the coefficient callables of the physical coordinates
(evaluated per element via [`element_coefficients`](@ref) — matrix-free discipline, CRP
main.tex:668).

Interface coupling: EVERY interior face goes through a [`MortarFlux`](@ref) built here
with the pairwise `(σ, Am)` constructed IDENTICALLY on both sides (W5-A contract): `Am` =
`𝔄` evaluated ONCE per interface mortar at the mortar's physical points and shared, `σ`
from mirrored [`Mortar`](@ref) constructor arguments. On conforming faces this is bitwise
[`PenaltyFlux`](@ref)`(C)` (gated). Boundary faces use `PenaltyFlux(C)` with the ghost
exterior from `bcs[(el, d, sidx)]::GhostBoundary` (every boundary face REQUIRED, nothing
else allowed). `hooks[(el, d, sidx)]::BoundaryDataHook` post-process received interior
interface data (per mortar, uniform on split faces; default [`NoHook`](@ref)).

Applications: [`apply!`](@ref) (LINEARIZED — the matrix; also available as the callable
`op(y, u)` for [`gmres`](@ref)), [`apply_full!`](@ref) (affine), [`operator_on_zero`](@ref)
(the fixed-data RHS lift), [`build_matrix`](@ref), [`direct_solve`](@ref).
"""
struct MeshOperator{S<:Number,RT<:AbstractFloat,TC<:Number,MT<:Mesh2D}
    mesh::MT
    C::RT
    ems::Vector{ElementMatrices{RT}}
    coeffs::Vector{ElementCoefficients{TC}}
    flux::Vector{PerFaceFlux}
    sendP::Vector{NTuple{2,NTuple{2,Vector{Matrix{RT}}}}}
    hooks::Vector{NTuple{2,NTuple{2,BoundaryDataHook}}}
    bcs::Vector{NTuple{2,NTuple{2,Union{Nothing,GhostBoundary}}}}
end

function MeshOperator(::Type{T}, mesh::Mesh2D, A::Tuple{Any,Any}, β, γ::Tuple{Any,Any};
                      C::Real = 1.5,
                      hooks::AbstractDict = Dict{Tuple{ElementId,Int,Int},BoundaryDataHook}(),
                      bcs::AbstractDict) where {T<:Number}
    nel = length(mesh.elements)
    # seam validation: bcs cover EXACTLY the boundary faces; hooks bind interior faces only
    bset = Set(boundary_faces(mesh))
    for k in keys(bcs)
        k in bset || throw(ArgumentError("MeshOperator: bc key $k is not a boundary face"))
    end
    for f in bset
        haskey(bcs, f) || throw(ArgumentError("MeshOperator: missing boundary condition " *
                                              "for boundary face $f"))
    end
    for k in keys(hooks)
        el, d, sidx = k
        haskey(mesh.index, el) && mesh.faces[mesh.index[el]][d][sidx] !== nothing ||
            throw(ArgumentError("MeshOperator: hook key $k is not an interior face"))
    end

    ems = [element_matrices(T, _element_block(mesh, e), mesh.elements[e], mesh.orders[e])
           for e in 1:nel]
    coeffs = [element_coefficients(ems[e], _element_block(mesh, e), mesh.elements[e],
                                   A, β, γ) for e in 1:nel]
    RT = typeof(ems[1]).parameters[1]

    # interface mortars, built ONCE per interface from its hi side: shared Am (evaluated
    # at the mortar's physical points), pairwise σ via mirrored constructor args
    nm(e, d, sidx) = mesh.faces[e][d][sidx] === nothing ? 0 :
                     length(mesh.faces[e][d][sidx].nbs)
    morts = [ntuple(d -> ntuple(sidx -> Vector{Any}(undef, nm(e, d, sidx)), 2), 2)
             for e in 1:nel]
    sends = [ntuple(d -> ntuple(sidx -> Vector{Matrix{RT}}(undef, nm(e, d, sidx)), 2), 2)
             for e in 1:nel]
    for e in 1:nel, d in 1:2
        iface = mesh.faces[e][d][2]
        iface === nothing && continue
        el = mesh.elements[e]
        blk = _element_block(mesh, e)
        tmap = blockmap(blk, 3 - d)
        _, hi_d = logical_interval(RT, el.segments[d])
        xf = RT(map_forward(blockmap(blk, d), hi_d))
        fm_int = face_matrices(ems[e], d, +1)
        for k in 1:length(iface.nbs)
            nb = iface.nbs[k]
            fm_ext = face_matrices(ems[nb], d, -1)
            Nm = max(size(fm_int.I, 2), size(fm_ext.I, 2))
            lo, hi = logical_interval(RT, iface.span[k])
            ξ = refelement(RT, Nm).nodes
            xt = [map_forward(tmap, (lo * (1 - x) + hi * (1 + x)) / 2) for x in ξ]
            Am = d == 1 ? [A[1](xf, tt) for tt in xt] : [A[2](tt, xf) for tt in xt]
            m_hi = Mortar(fm_int, fm_ext, iface.msize[k], Am; C = C)
            m_lo = Mortar(fm_ext, fm_int, iface.nbmsize[k], Am; C = C)
            morts[e][d][2][k] = m_hi
            sends[e][d][2][k] = m_lo.P
            slot = iface.nbslot[k]
            morts[nb][d][1][slot] = m_lo
            sends[nb][d][1][slot] = m_hi.P
        end
    end
    pen = PenaltyFlux(RT(C))
    flux = [PerFaceFlux(ntuple(d -> ntuple(sidx ->
                mesh.faces[e][d][sidx] === nothing ? pen :
                    MortarFlux([m for m in morts[e][d][sidx]]), 2), 2)) for e in 1:nel]
    hk = [ntuple(d -> ntuple(sidx ->
              get(hooks, (mesh.elements[e], d, sidx), NoHook())::BoundaryDataHook, 2), 2)
          for e in 1:nel]
    bc = [ntuple(d -> ntuple(sidx -> mesh.faces[e][d][sidx] === nothing ?
              (bcs[(mesh.elements[e], d, sidx)]::GhostBoundary) : nothing, 2), 2)
          for e in 1:nel]
    TC = typeof(coeffs[1]).parameters[1]
    S = promote_type(float(T), TC)
    return MeshOperator{S,RT,TC,typeof(mesh)}(mesh, RT(C), ems, coeffs, flux,
                                              [s for s in sends], hk, bc)
end

_scalar(::MeshOperator{S}) where {S} = S

# BigFloat is never threaded (binding owner directive: processes, not threads)
_default_threaded(op::MeshOperator{S}, u) where {S} =
    Threads.nthreads() > 1 && isbitstype(promote_type(S, eltype(u)))

# hook application on one mortar's data: NoHook skips the interior projection entirely
_hooked(::NoHook, P, int::FaceState, raw::FaceState, lin::Bool) = raw
function _hooked(h::BoundaryDataHook, P, int::FaceState, raw::FaceState, lin::Bool)
    intm = project_to_mortar(P, int)
    return lin ? modify_boundary_data_linearized(h, intm, raw) :
                 modify_boundary_data(h, intm, raw)
end

# exterior FaceState of one face: ghost (boundary) or hook-processed packed mortar data
function _exterior(op::MeshOperator, sts, e::Int, d::Int, sidx::Int, lin::Bool)
    iface = op.mesh.faces[e][d][sidx]
    int = sts[e].faces[d][sidx]
    if iface === nothing
        bc = op.bcs[e][d][sidx]
        return lin ? ghost_state_linearized(bc, int) : ghost_state(bc, int)
    end
    mf = op.flux[e].faces[d][sidx]::MortarFlux
    hook = op.hooks[e][d][sidx]
    parts = [_hooked(hook, mf.mortars[k].P, int,
                     project_to_mortar(op.sendP[e][d][sidx][k],
                                       sts[iface.nbs[k]].faces[d][3-sidx]), lin)
             for k in 1:length(iface.nbs)]
    length(parts) == 1 && return parts[1]
    return FaceState(vcat([p.psi for p in parts]...), vcat([p.nflux for p in parts]...))
end

_state(op::MeshOperator, u, e::Int) =
    element_state(op.ems[e], op.coeffs[e], view(u, _element_range(op.mesh, e)))

function _residual!(y, op::MeshOperator, sts, e::Int, lin::Bool)
    ext = ntuple(d -> ntuple(sidx -> _exterior(op, sts, e, d, sidx, lin), 2), 2)
    strong_residual!(view(y, _element_range(op.mesh, e)), op.ems[e], op.coeffs[e],
                     sts[e], ext, op.flux[e])
    return nothing
end

function _apply!(y::AbstractVector, op::MeshOperator{S,RT,TC}, u::AbstractVector,
                 lin::Bool, threaded::Bool) where {S,RT,TC}
    n = op.mesh.ndofs
    length(u) == n ||
        throw(ArgumentError("apply: u has length $(length(u)), the mesh has $n dofs"))
    length(y) == n ||
        throw(ArgumentError("apply: y has length $(length(y)), the mesh has $n dofs"))
    nel = length(op.mesh.elements)
    Tp = promote_type(eltype(u), TC, RT)
    sts = Vector{ElementState{Tp}}(undef, nel)
    if threaded
        Threads.@threads for e in 1:nel               # pass 1: build + send face data
            sts[e] = _state(op, u, e)
        end
        Threads.@threads for e in 1:nel               # pass 2: receive + flux + lift
            _residual!(y, op, sts, e, lin)
        end
    else
        for e in 1:nel
            sts[e] = _state(op, u, e)
        end
        for e in 1:nel
            _residual!(y, op, sts, e, lin)
        end
    end
    return y
end

"""
    apply!(y, op::MeshOperator, u; threaded = auto) -> y

The LINEARIZED global operator application `y = A·u` — the matrix every solver component
sees (W5-B contract): received interface data route through
`modify_boundary_data_linearized`, boundary ghosts through `ghost_state_linearized`, so
`A` is linear-homogeneous with fixed data excluded. Two-pass element sweep; `threaded`
defaults to `Threads.nthreads() > 1` for bits types (never for BigFloat — owner
directive) and is BITWISE identical to serial. `op(y, u)` forwards here, making the
operator directly usable as the [`gmres`](@ref) `applyA!`.
"""
apply!(y::AbstractVector, op::MeshOperator, u::AbstractVector;
       threaded::Bool = _default_threaded(op, u)) = _apply!(y, op, u, true, threaded)

"""
    apply_full!(y, op::MeshOperator, u; threaded = auto) -> y

The FULL affine operator application: hooks act with their fixed data
(`modify_boundary_data`) and ghosts carry their boundary data (`ghost_state`), so
`apply_full!(y, op, u) == apply!(y, op, u) + operator_on_zero(op)` (gated at 1e-13).
Used for the operator-on-zero RHS lift and affine-residual checks — matrix applications
use [`apply!`](@ref).
"""
apply_full!(y::AbstractVector, op::MeshOperator, u::AbstractVector;
            threaded::Bool = _default_threaded(op, u)) = _apply!(y, op, u, false, threaded)

(op::MeshOperator)(y::AbstractVector, u::AbstractVector) = apply!(y, op, u)

"""
    operator_on_zero(op::MeshOperator) -> Vector

The fixed-data source lift `A_full(0)`: applying the FULL operator to the zero trial
field yields exactly the lifted fixed-hook shifts and inhomogeneous boundary data (the
spectre `impose_inhomogeneous_boundary_conditions_on_source` pattern, DgOperator.hpp:
1081-1130). Solve the inhomogeneous problem as `A_lin·u = b − operator_on_zero(op)` with
the linearized [`apply!`](@ref) matrix; equals the W5-B per-face closed forms (gated at
1e-13) and is exactly zero for purely linearized hook sets.
"""
function operator_on_zero(op::MeshOperator)
    S = _scalar(op)
    n = op.mesh.ndofs
    return apply_full!(zeros(S, n), op, zeros(S, n))
end

"""
    build_matrix(op::MeshOperator) -> Matrix

Dense global matrix of the LINEARIZED operator by unit-vector column probes through the
matrix-free [`apply!`](@ref) path (small meshes: tests, fig-matrix, the multigrid bottom
solver, [`direct_solve`](@ref)). Linearized hooks are IN the matrix; fixed data are not
(falsifiers gated both ways in test_matrixfree.jl).
"""
function build_matrix(op::MeshOperator)
    S = _scalar(op)
    n = op.mesh.ndofs
    A = zeros(S, n, n)
    u = zeros(S, n)
    y = zeros(S, n)
    for j in 1:n
        u[j] = 1
        apply!(y, op, u)
        A[:, j] = y
        u[j] = 0
    end
    return A
end

"""
    direct_solve(op::MeshOperator, b) -> Vector

Direct-solve bypass for small problems: `build_matrix(op) \\ b` (dense LU; generic — a
BigFloat operator solves in BigFloat). `b` must already carry the operator-on-zero
correction for fixed data (`b − operator_on_zero(op)`). The spectre
`EnableDirectSolve`-style twin; the cached-LU production path is the
multigrid-bottom-solver node.
"""
direct_solve(op::MeshOperator, b::AbstractVector) = build_matrix(op) \ Vector(b)

# mesh2d.jl — 2D multiblock mesh layer: Block coordinate maps (linear + inverse-radial
# eq:compaction family), refinement-tree ElementId (block id + per-dimension level/index),
# parent/child/neighbor arithmetic, dyadic logical intervals, 2:1-balance predicate.
# DAG node: mesh2d-blocks (elliptic-ssf GLOBAL_DAG row 9, TOPOLOGICAL_SCHEDULE W2-B).
#
# The CRP domain (main.tex:541-554) splits the (r, z) rectangle into irreducible blocks, each
# carrying a coordinate map to logical coordinates ξ ∈ [-1,1]: linear in the v/t regions, and
# radially ∝ 1/r in the u wavezone (eq:compaction, main.tex:546-549). Blocks are h-refined by
# repeated halving with a two-to-one balance at element interfaces (main.tex:552-554).
#
# Numerical contracts (test_mesh2d.jl gates):
#   * map endpoints are EXACT: the linear map is evaluated in the endpoint-exact form
#     (a(1-ξ) + b(1+ξ))/2 (at ξ = ∓1 one term is exactly 0 and the other exactly 2a or 2b —
#     scaling by 2 is exact in binary FP); the inverse-radial map snaps ξ == ±1 to r_u / R,
#     since no FP rearrangement of eq:compaction guarantees fl(2 r_u R / (2R)) == r_u.
#   * the inverse-radial denominator is computed as r_u(1+ξ) + R(1-ξ) — a sum of two
#     nonnegative terms (no cancellation near ξ = +1, where the naive r_u + R + (r_u - R)ξ
#     loses ~eps·R absolute → ~eps·R/(2r_u) relative accuracy at the steep outer end).
#   * geometry is real (RT), inherited from the map's construction type; a BigFloat logical
#     coordinate promotes the result (DESIGN §4 precision discipline).
#   * refinement-tree arithmetic is pure Int64: segment (level, index) with
#     0 ≤ index < 2^level, level ≤ 62 so dyadic interval endpoints stay exact under shifts.

export CoordinateMap1D, LinearMap1D, InverseRadialMap1D,
       map_forward, map_inverse, map_jacobian,
       Block, blockmap,
       SegmentId, ElementId, child_id, parent_id, neighbor_id,
       logical_interval, two_to_one_balanced

# ---- coordinate maps -------------------------------------------------------------------------

"""
    CoordinateMap1D{RT<:AbstractFloat}

Abstract 1D coordinate map from the logical coordinate `ξ ∈ [-1, 1]` to a physical coordinate,
with geometry stored in the real type `RT`. Concrete maps ([`LinearMap1D`](@ref),
[`InverseRadialMap1D`](@ref)) implement [`map_forward`](@ref), [`map_inverse`](@ref), and the
analytic [`map_jacobian`](@ref). All maps are orientation-preserving (positive Jacobian on
`[-1, 1]`), enforced at construction. A 2D block map is the tensor product of two 1D maps
(one per dimension), matching the CRP domain where slicing regions differ only radially.
"""
abstract type CoordinateMap1D{RT<:AbstractFloat} end

"""
    LinearMap1D(a, b)

Affine map `x(ξ) = (a(1-ξ) + b(1+ξ))/2` with **exact** endpoints `x(-1) == a`, `x(+1) == b`
(this form makes endpoint evaluation exact in binary FP: one term vanishes exactly and
`2a/2 == a`). Requires `a < b`. Used for the radial v/t regions and the angular dimension
of every CRP block (main.tex:545).
"""
struct LinearMap1D{RT<:AbstractFloat} <: CoordinateMap1D{RT}
    a::RT
    b::RT
    function LinearMap1D{RT}(a::RT, b::RT) where {RT<:AbstractFloat}
        a < b || throw(ArgumentError("LinearMap1D needs a < b (orientation-preserving), got a=$a, b=$b"))
        return new{RT}(a, b)
    end
end
function LinearMap1D(a::Real, b::Real)
    af, bf = promote(float(a), float(b))
    return LinearMap1D{typeof(af)}(af, bf)
end

"""
    InverseRadialMap1D(r_u, R)

The CRP wavezone compactification map (eq:compaction, main.tex:546-549)

    r(ξ) = 2 r_u R / (r_u + R + (r_u - R) ξ),

i.e. `1/r` is the affine interpolant between `1/r_u` (at `ξ = -1`) and `1/R` (at `ξ = +1`),
resolving the large `u`-domain wavezone with radially increasing grid spacing. Requires
`0 < r_u < R`. Endpoints are **exact**: `r(-1) == r_u`, `r(+1) == R` (snapped — no FP
rearrangement of the formula guarantees them). The denominator is evaluated in the
cancellation-free form `r_u(1+ξ) + R(1-ξ)`.
"""
struct InverseRadialMap1D{RT<:AbstractFloat} <: CoordinateMap1D{RT}
    r_u::RT
    R::RT
    function InverseRadialMap1D{RT}(r_u::RT, R::RT) where {RT<:AbstractFloat}
        0 < r_u < R || throw(ArgumentError("InverseRadialMap1D needs 0 < r_u < R, got r_u=$r_u, R=$R"))
        return new{RT}(r_u, R)
    end
end
function InverseRadialMap1D(r_u::Real, R::Real)
    rf, Rf = promote(float(r_u), float(R))
    return InverseRadialMap1D{typeof(rf)}(rf, Rf)
end

"""
    map_forward(m, ξ) -> x

Physical coordinate of logical `ξ ∈ [-1, 1]`. Endpoints `ξ == ±1` are exact by contract.
The result type promotes the map's geometry type with `typeof(ξ)` (a `BigFloat` argument on a
`Float64` map yields `BigFloat`).
"""
map_forward(m::LinearMap1D, ξ::Real) = (m.a * (1 - ξ) + m.b * (1 + ξ)) / 2

function map_forward(m::InverseRadialMap1D, ξ::Real)
    den = m.r_u * (1 + ξ) + m.R * (1 - ξ)     # = r_u + R + (r_u - R)ξ, cancellation-free
    r = 2 * m.r_u * m.R / den
    ξ == -1 && return oftype(r, m.r_u)        # eq:compaction endpoints are EXACT by contract
    ξ == 1 && return oftype(r, m.R)
    return r
end

"""
    map_inverse(m, x) -> ξ

Logical coordinate of physical `x`; analytic inverse of [`map_forward`](@ref). The round trip
`map_forward(m, map_inverse(m, x))` reproduces `x` to ≤ 1e-13 relative at Float64 across both
map families (gated in test_mesh2d.jl; the worst case is the steep outer end of the
inverse-radial map, calibrated at 5.1e-14 for the paper's (r_u, R) = (25.3, 1e4) contrast).
"""
map_inverse(m::LinearMap1D, x::Real) = (2x - (m.a + m.b)) / (m.b - m.a)
map_inverse(m::InverseRadialMap1D, r::Real) = (m.r_u + m.R - 2 * m.r_u * m.R / r) / (m.R - m.r_u)

"""
    map_jacobian(m, ξ) -> dx/dξ

Analytic Jacobian of [`map_forward`](@ref) at `ξ`; strictly positive on `[-1, 1]` for every
constructible map. Linear: `(b-a)/2` (constant). Inverse-radial: `2 r_u R (R - r_u) / den²`
with `den = r_u(1+ξ) + R(1-ξ)`, equivalently `(R - r_u) r(ξ)² / (2 r_u R)`.
"""
map_jacobian(m::LinearMap1D, ξ::Real) = ((m.b - m.a) / 2) * one(ξ)

function map_jacobian(m::InverseRadialMap1D, ξ::Real)
    den = m.r_u * (1 + ξ) + m.R * (1 - ξ)
    return 2 * m.r_u * m.R * (m.R - m.r_u) / den^2
end

# ---- blocks ------------------------------------------------------------------------------------

"""
    Block(id, xmap, ymap)

Irreducible domain block (CRP main.tex:543): an integer `id` plus one [`CoordinateMap1D`](@ref)
per dimension mapping the logical square `[-1,1]²` to the block's physical rectangle. Both maps
must share the same geometry type `RT` (enforced by the type system — mixing `Float64` and
`BigFloat` maps is a `MethodError`). Access per-dimension maps with [`blockmap`](@ref).
"""
struct Block{RT<:AbstractFloat,MX<:CoordinateMap1D{RT},MY<:CoordinateMap1D{RT}}
    id::Int
    maps::Tuple{MX,MY}
end
Block(id::Integer, xmap::CoordinateMap1D{RT}, ymap::CoordinateMap1D{RT}) where {RT} =
    Block{RT,typeof(xmap),typeof(ymap)}(Int(id), (xmap, ymap))

"""
    blockmap(blk::Block, dim) -> CoordinateMap1D

The block's coordinate map along dimension `dim ∈ (1, 2)`.
"""
blockmap(blk::Block, dim::Integer) = blk.maps[dim]

# ---- refinement-tree addresses -------------------------------------------------------------------

"""
    SegmentId(level, index)

Address of a dyadic segment of the logical interval `[-1, 1]` in one dimension: at refinement
`level ≥ 0` the interval is split into `2^level` equal segments and `index ∈ 0:2^level-1`
counts from the lower end. `level ≤ 62` keeps all dyadic-endpoint arithmetic exact in `Int64`
(the 2:1 predicate shifts indices to a common level).
"""
struct SegmentId
    level::Int
    index::Int
    function SegmentId(level::Integer, index::Integer)
        0 <= level <= 62 || throw(ArgumentError("SegmentId level must be in 0:62, got $level"))
        0 <= index < (Int64(1) << level) ||
            throw(ArgumentError("SegmentId index must be in 0:2^$level-1, got $index"))
        return new(Int(level), Int(index))
    end
end

"""
    ElementId(block, sx::SegmentId, sy::SegmentId)

Refinement-tree address of a mesh element: the owning [`Block`](@ref) id plus one
[`SegmentId`](@ref) per dimension (per-dimension level/index — CRP elements are produced by
repeatedly halving the logical square along either dimension, main.tex:552). Plain immutable
data: equality/hash are structural, so ids serve directly as `Dict`/`Set` keys.
"""
struct ElementId
    block::Int
    segments::NTuple{2,SegmentId}
end
ElementId(block::Integer, sx::SegmentId, sy::SegmentId) = ElementId(Int(block), (sx, sy))

# replace one dimension's segment, keeping block and the other dimension untouched
_with_segment(id::ElementId, dim::Integer, s::SegmentId) =
    ElementId(id.block, Base.setindex(id.segments, s, dim))

"""
    child_id(id::ElementId, dim, which) -> ElementId

The `which ∈ (0, 1)` half (lower/upper) of element `id` after halving along dimension `dim`:
level `+1`, index `2·index + which` in that dimension; the other dimension and the block are
untouched. Inverse of [`parent_id`](@ref): `parent_id(child_id(id, dim, w), dim) == id`.
"""
function child_id(id::ElementId, dim::Integer, which::Integer)
    which == 0 || which == 1 || throw(ArgumentError("child_id: which must be 0 or 1, got $which"))
    s = id.segments[dim]
    return _with_segment(id, dim, SegmentId(s.level + 1, 2 * s.index + which))
end

"""
    parent_id(id::ElementId, dim) -> ElementId

The element that `id` was split from along dimension `dim` (level `-1`, index halved).
Throws `ArgumentError` at the block root (`level == 0`) — coarsening beyond the root is a
logic error in any refinement driver.
"""
function parent_id(id::ElementId, dim::Integer)
    s = id.segments[dim]
    s.level >= 1 || throw(ArgumentError("parent_id: dimension $dim is at root level 0"))
    return _with_segment(id, dim, SegmentId(s.level - 1, s.index >> 1))
end

"""
    neighbor_id(id::ElementId, dim, dir) -> Union{ElementId, Nothing}

The same-level face-neighbor address of `id` along dimension `dim` in direction
`dir ∈ (-1, +1)`, or `nothing` when the face lies on the block boundary. Symmetric:
`neighbor_id(neighbor_id(id, dim, dir), dim, -dir) == id` whenever the first call is not
`nothing`. (An actual mesh neighbor may live at a different level — that is the 2:1
predicate's concern, not the address arithmetic's.)
"""
function neighbor_id(id::ElementId, dim::Integer, dir::Integer)
    dir == 1 || dir == -1 || throw(ArgumentError("neighbor_id: dir must be +1 or -1, got $dir"))
    s = id.segments[dim]
    j = s.index + dir
    0 <= j < (Int64(1) << s.level) || return nothing
    return _with_segment(id, dim, SegmentId(s.level, j))
end

"""
    logical_interval(RT, s::SegmentId) -> (lo, hi)
    logical_interval(RT, id::ElementId, dim) -> (lo, hi)

The dyadic subinterval of `[-1, 1]` addressed by segment `s` (or by `id`'s segment along
`dim`), in the real type `RT`. Endpoints are integer multiples of `2^(1-level)`, hence EXACT
in any binary type whose mantissa holds the index (all practical levels): children bisect the
parent exactly and same-level neighbors share their face point bitwise.
"""
function logical_interval(::Type{RT}, s::SegmentId) where {RT<:AbstractFloat}
    w = ldexp(one(RT), 1 - s.level)           # segment width 2^(1-level), exact
    lo = s.index * w - 1                       # integer × power-of-two − 1: exact dyadic
    return (lo, lo + w)
end
logical_interval(::Type{RT}, id::ElementId, dim::Integer) where {RT<:AbstractFloat} =
    logical_interval(RT, id.segments[dim])

# ---- 2:1 balance ---------------------------------------------------------------------------------

# segment endpoints in integer units of 2^-L (grid of the common refinement level L ≥ level);
# (index+1) << (L - level) ≤ 2^L ≤ 2^62: exact in Int64 by the SegmentId level cap
function _endpoints_at(s::SegmentId, L::Int)
    sh = L - s.level
    return (Int64(s.index) << sh, Int64(s.index + 1) << sh)
end

# (touch, overlap) of two same-dimension segments on their common dyadic integer grid
function _touch_overlap(sa::SegmentId, sb::SegmentId)
    L = max(sa.level, sb.level)
    alo, ahi = _endpoints_at(sa, L)
    blo, bhi = _endpoints_at(sb, L)
    return (ahi == blo || bhi == alo, max(alo, blo) < min(ahi, bhi))
end

"""
    two_to_one_balanced(ids::AbstractVector{ElementId}) -> Bool

CRP two-to-one balance predicate (main.tex:554): `true` iff every pair of face-adjacent
same-block elements differs by at most one refinement level in **each** dimension
("neighboring elements can differ by at most one h refinement level in any dimension").
Face adjacency means the elements touch along one dimension AND overlap with positive measure
in the other — corner-only contact is not an interface. All interval arithmetic is exact
(integer dyadic shifts to the pair's common level).

Scope: `ids` is a leaf set of per-block refinement trees (elements tile, never overlap —
not validated here); adjacency is same-block only. Balance across block boundaries needs the
block-connectivity/mortar layer (`dg-mortars`), which owns that check.
"""
function two_to_one_balanced(ids::AbstractVector{ElementId})
    n = length(ids)
    for i in 1:n-1
        a = ids[i]
        for j in i+1:n
            b = ids[j]
            a.block == b.block || continue
            t1, o1 = _touch_overlap(a.segments[1], b.segments[1])
            t2, o2 = _touch_overlap(a.segments[2], b.segments[2])
            (t1 && o2) || (t2 && o1) || continue                  # not face-adjacent
            abs(a.segments[1].level - b.segments[1].level) <= 1 || return false
            abs(a.segments[2].level - b.segments[2].level) <= 1 || return false
        end
    end
    return true
end

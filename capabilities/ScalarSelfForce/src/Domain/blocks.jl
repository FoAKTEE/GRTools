# Domain/blocks.jl — domain-blocks-map (ESSF node, schedule W7-D, GLOBAL_DAG row 39): the
# PURE DATA mapping from the sealed yaml-twin [`RunConfig`](@ref) (Domain/config.jl) to an
# EllipticCore [`Mesh2D`](@ref) — the CRP computational domain of arXiv:2606.04998
# (main.tex:541-554, fig:domain main.tex:529).
#
# The domain is the (r, z) rectangle from the inner (horizon-penetrating) radius to the outer
# radius R and from the equator to the pole, cut into IRREDUCIBLE BLOCKS on a rectangular
# lattice — 5 radial columns × 2 angular rows in the equatorial-symmetry (r, z²) chart of the
# production yaml, or 5 × 3 in the full (r, z) chart (main.tex:541-543). The radial columns
# are, in order, the horizon-penetrating wavezone (v slicing), a transition region, the
# worldtube region, another transition region, and the outer wavezone (u slicing); the angular
# rows are the worldtube row and the polar row(s). Each block carries a coordinate map to the
# logical square: LINEAR everywhere except the outermost radial column, which carries the 1/r
# wavezone compactification eq:compaction (main.tex:546-549) so the huge u domain is resolved
# with radially increasing spacing.
#
# ZERO new physics and ZERO duplicated coefficient algebra. This file decides only WHICH
# sealed object each block gets:
#   * the vtu slicing region of a block comes from `vtu_domain(vtu_slicing(orbit, config), r)`
#     (vtu.jl) evaluated at the block's radial midpoint — never from a re-derived H;
#   * a block's DG coefficients are the sealed `fluxform_callables(orbit, config, region)`
#     bundle (fluxform.jl / z2.jl — which also picks the (r, z) or (r, z²) chart from
#     `ImposeEquatorialSymmetry`);
#   * every number in the lattice (bounds, refinement levels, grid points) is read from the
#     sealed `RunConfig`, which already reproduces the yaml's Jinja2 derivations exactly.
#
# BLOCK NUMBERING (binding, the yaml's own): the lattice index (i, j) — i the radial column,
# j the angular row, both 0-based — maps to the flat block id `i + nr·j`, i.e. the RADIAL
# index runs fastest. This is spectre's `AlignedLattice` ordering, pinned by the production
# yaml's `NullSlicingBlocks: [Block0, Block5, Block4, Block9]` (yaml:39-45): the v-slicing
# blocks are the radial column 0 of both rows and the u-slicing blocks are the radial column 4
# of both rows, which is only consistent with i-fastest numbering.
#
# ANGULAR ROWS WITHOUT EQUATORIAL SYMMETRY: the yaml lattice is specified in z², rows
# {0, dz², 1}. The full-z lattice covers the SAME physical domain, so its rows are the SIGNED
# SQUARE ROOTS {−1, −√dz², +√dz², 1} — exactly the angular `BlockBounds` of spectre's
# non-equatorial input file (`tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml:59`,
# [−1, −0.5, 0.5, 1] for dz² = 0.25). The equator is then an interior POINT of the middle
# block row — z = 0 is not a block bound in the full chart — and that middle row is the
# worldtube row.
#
# WHAT LATER NODES BIND TO (the contract this node establishes):
#   * `boundary_face_classes` / `boundary_class` — every domain-boundary face tagged
#     `:horizon` (r = inner bound), `:outer` (r = R), `:equator` (z² = 0, equatorial chart
#     only) or `:pole` (z² = 1, or z = ±1). Node `dg-boundary-conditions` binds a
#     `GhostBoundary` to exactly these faces (the `MeshOperator` `bcs` table must cover them
#     all): Bayliss-Turkel on `:outer` for m ≠ 0 (`None` at m = 0, yaml:115-125), and nothing
#     physical on `:horizon`/`:equator`/`:pole` — the flux coefficients vanish there
#     (F^r(r₊) = 0, F^z(±1) = 0, main.tex:385), which is why the paper needs no boundary
#     condition on those faces. For `:horizon` that is a PRECONDITION, not a consequence: it
#     holds only because the inner radial bound is placed AT r₊ (main.tex:529), which
#     `block_domain` now ENFORCES — at r_in = 2 with a = 0.5 (r₊ = 1.866) the sealed
#     coefficient measures A^rr = 5.88e-2, and inside the horizon (a = 0.9, r_in = 1.3) it is
#     NEGATIVE (−4.0e-2), so those faces would silently need data node 41 does not supply.
#   * `worldtube_blocks` — the regularized-region block (the red rectangle of fig:domain), the
#     one containing the particle at (r₀, equator); node `worldtube-formulation` imposes the
#     Ψ^R = Ψ − Ψ^P jump on its interfaces with the rest of the domain.
#   * `block_coefficients` / `domain_coefficients` — the per-block and the whole-domain
#     coefficient bundles the solver driver feeds to `MeshOperator`.

export BlockDomain, block_domain, block_lattice_bounds, block_lattice_index, block_region,
       worldtube_blocks, boundary_face_classes, boundary_class, boundary_class_rule,
       block_coefficients, domain_coefficients

"""
    BlockDomain{RT}

The instantiated CRP block lattice (arXiv:2606.04998 main.tex:541-554): an EllipticCore
[`Mesh2D`](@ref) plus the per-block metadata the worldtube, boundary-condition, and solver
nodes bind to. Build it with [`block_domain`](@ref).

Public fields:
- `mesh::Mesh2D{RT}` — the blocks, their initial elements (layout order: block-major with the
  radial index fastest, then `ix` fastest within a block), the block connections, and the
  precomputed face topology.
- `nr::Int`, `nz::Int` — lattice size (5 × 2 with equatorial symmetry, 5 × 3 without).
- `chart::Symbol` — `:z2` (angular coordinate `w = z²` on `[0, 1]`, equatorial symmetry) or
  `:z` (angular coordinate `z = cosθ` on `[−1, 1]`). Matches the chart the sealed
  [`fluxform_callables`](@ref)`(orbit, config, region)` selects.
- `bounds_r::Vector{RT}`, `bounds_z::Vector{RT}` — the `nr+1` radial and `nz+1` angular block
  bounds (the radial ones are the config's own `BlockBounds[r]`).
- `block_ids::Matrix{Int}` — `block_ids[i+1, j+1]` is the flat id `i + nr·j` of lattice cell
  `(i, j)` (yaml block numbering, radial index fastest).
- `regions::Dict{Int,Symbol}` — block id → vtu slicing region `:v` / `:t` / `:u`.
- `refinements::Dict{Int,NTuple{2,Int}}` — block id → initial h-refinement levels `(nx, ny)`,
  so the block holds `2^nx × 2^ny` elements.
- `worldtube::Vector{Int}` — the regularized-region block id(s).
- `boundary_classes::Dict{Tuple{ElementId,Int,Int},Symbol}` — every domain-boundary face key
  of `mesh` tagged `:horizon` / `:outer` / `:equator` / `:pole`.
"""
struct BlockDomain{RT<:AbstractFloat}
    mesh::Mesh2D{RT}
    nr::Int
    nz::Int
    chart::Symbol
    bounds_r::Vector{RT}
    bounds_z::Vector{RT}
    block_ids::Matrix{Int}
    regions::Dict{Int,Symbol}
    refinements::Dict{Int,NTuple{2,Int}}
    worldtube::Vector{Int}
    boundary_classes::Dict{Tuple{ElementId,Int,Int},Symbol}
end

"""
    block_lattice_bounds(config::RunConfig) -> (bounds_r, bounds_z, chart)

The block lattice of `config` as pure data: the radial bounds (the config's own
`DomainCreator.AlignedLattice.BlockBounds[r]`, yaml:55-61 — `{inner, r₀−1.5w, r₀−w, r₀+w,
r₀+3w, R}` at the production point), the angular bounds, and the angular `chart`.

With `ImposeEquatorialSymmetry` (yaml:36, the production setting) the angular chart is
`w = z²` and the rows are the config's `BlockBounds[z²]` verbatim (`{0, dz², 1}`, yaml:63).
Without it the chart is `z = cosθ` over the full range and the rows are the SIGNED SQUARE
ROOTS of the z² rows — `{−1, −√dz², +√dz², 1}`, the same physical domain, matching the
angular `BlockBounds` of spectre's non-equatorial input file — so the lattice grows from
`nr × nz` to `nr × (2nz − 1)` rows and the equator becomes an interior interface.
"""
function block_lattice_bounds(config::RunConfig{T}) where {T}
    r = config.domain.block_bounds_r
    z2 = config.domain.block_bounds_z2
    config.orbit.impose_equatorial_symmetry && return (r, z2, :z2)
    (length(z2) >= 2 && iszero(z2[1])) || throw(ArgumentError(
        "block_lattice_bounds: the z² lattice must start at the equator z² = 0; got $z2"))
    zpos = T[sqrt(w) for w in z2[2:end]]              # {√dz², …, 1}
    return (r, vcat(-reverse(zpos), zpos), :z)
end

# per-radial-column initial refinement from DomainCreator.AlignedLattice.RefinedLevels
# (yaml:78-110). Both reference lattices — the production 5×2 table and spectre's
# non-equatorial 5×3 twin — specify the refinement as a function of the RADIAL COLUMN ALONE
# (the same (nx, ny) in every angular row), which is what lets the same config drive either
# number of rows; a table that is not row-uniform is refused rather than guessed at.
# Columns no patch covers keep `InitialLevels` (yaml:75). Both corner indices are bounds-
# checked against the lattice, as spectre does (`AlignedLattice.cpp:123-131` PARSE_ERRORs when
# `upper_corner_index[i] >= block_bounds_[i].size()` in EITHER dimension). DELIBERATE
# DIVERGENCE from spectre: where overlapping patches would silently leave spectre's last write
# standing, a conflicting overlap is an `ArgumentError` here — an ambiguous refinement table
# must not be resolved by input ordering.
function _column_refinements(config::RunConfig, nr::Int, nz::Int)
    L = config.domain.initial_levels
    col = Dict{Int,NTuple{2,Int}}()
    for p in config.domain.refined_levels
        lo, hi = p.lower_corner_index, p.upper_corner_index
        (hi[1] > lo[1] && hi[2] > lo[2]) || throw(ArgumentError(
            "block_domain: RefinedLevels patch $lo → $hi is empty"))
        hi[1] <= nr || throw(ArgumentError(
            "block_domain: RefinedLevels patch $lo → $hi exceeds the $nr radial columns"))
        hi[2] <= nz || throw(ArgumentError(
            "block_domain: RefinedLevels patch $lo → $hi exceeds the $nz angular rows"))
        for i in lo[1]:hi[1]-1
            r = (Int(p.refinement[1]), Int(p.refinement[2]))
            get(col, i, r) == r || throw(ArgumentError(
                "block_domain: RefinedLevels is not a function of the radial column alone " *
                "(column $i carries both $(col[i]) and $r) — the lattice cannot be " *
                "instantiated for a different number of angular rows"))
            col[i] = r
        end
    end
    return [get(col, i, (Int(L[1]), Int(L[2]))) for i in 0:nr-1]
end

# the vtu region of every radial column, from the SEALED classifier at the column midpoint.
# Requires each slicing transition to BE a block bound (main.tex:323 — "they coincide with
# radial block bounds"): otherwise a block would straddle two slicing slabs and no single
# per-block boost H would be correct.
function _column_regions(v::VtuSlicing, config::RunConfig, bounds_r::AbstractVector)
    t = config.orbit.hyperboloidal_slicing_transitions
    for (name, x) in (("v→t", t[1]), ("t→u", t[3]))
        any(==(x), bounds_r) || throw(ArgumentError(
            "block_domain: the $name slicing transition $x is not a radial block bound " *
            "$(bounds_r) — a block would straddle two slicing regions (main.tex:323)"))
    end
    regs = [vtu_domain(v, (bounds_r[i] + bounds_r[i+1]) / 2) for i in 1:length(bounds_r)-1]
    rank(s) = s === :v ? 1 : s === :t ? 2 : 3
    issorted(regs; by = rank) || throw(ArgumentError(
        "block_domain: the radial slicing regions $regs are not contiguous v…t…u"))
    return regs
end

"""
    boundary_class_rule(chart::Symbol, d::Integer, sidx::Integer) -> Symbol

THE boundary-face classification rule of the CRP lattice, as a pure function of the angular
`chart` (`:z2` or `:z`) and the face key's dimension/side — the single normative statement
that [`block_domain`](@ref) uses to build [`boundary_face_classes`](@ref) and that
`Solve/bc_glue.jl` re-applies to the COARSER meshes a multigrid hierarchy hands its
`build_op` (whose element ids do not exist in `dom.boundary_classes`):

    (d = 1, sidx = 1) → :horizon      (r = r₊, the inner radial bound)
    (d = 1, sidx = 2) → :outer        (r = R)
    (d = 2, sidx = 1) → :equator in the `:z2` chart, :pole in the `:z` chart
    (d = 2, sidx = 2) → :pole

It is valid only for faces ON the lattice perimeter; the callers establish that separately
(`block_domain` from [`boundary_faces`](@ref) of a fully connected rectangle,
[`boundary_conditions`](@ref) from an explicit perimeter check). Factored out so the two
call sites cannot drift — that they agree on `dom.mesh` AND on a coarsened mesh is gated in
test_bc_glue.jl. `d`/`sidx` outside `(1, 2)` throw `ArgumentError`.
"""
function boundary_class_rule(chart::Symbol, d::Integer, sidx::Integer)
    (d == 1 || d == 2) ||
        throw(ArgumentError("boundary_class_rule: dim must be 1 or 2, got $d"))
    (sidx == 1 || sidx == 2) ||
        throw(ArgumentError("boundary_class_rule: sidx must be 1 or 2, got $sidx"))
    return d == 1 ? (sidx == 1 ? :horizon : :outer) :
                    (sidx == 1 ? (chart === :z2 ? :equator : :pole) : :pole)
end

"""
    block_domain(orbit::CircularOrbit, config::RunConfig; orders = InitialGridPoints)
        -> BlockDomain

Instantiate the CRP block lattice of `config` as an EllipticCore [`Mesh2D`](@ref)
(main.tex:541-554). `orbit` is the physical point the config describes — it must agree with
`config` on `(M, a, r₀)` (checked by the sealed [`vtu_slicing`](@ref)), because the slicing
regions and hence the per-block coefficients belong to that orbit.

What is built, all of it read from the sealed config:

- **Blocks.** One per lattice cell `(i, j)`, id `i + nr·j` (radial index fastest — the yaml's
  own numbering, see the file header). The angular map is always
  [`LinearMap1D`](@ref)`(z_j, z_{j+1})`; the radial map is `LinearMap1D(r_i, r_{i+1})` except
  in the outermost column, which is the wavezone [`InverseRadialMap1D`](@ref)`(r_u, R)` of
  eq:compaction (main.tex:546-549) — matching `Distributions[r] = [Linear ×(nr−1), Inverse]`
  (yaml:66), which is verified rather than assumed. Map endpoints are exact, so neighbouring
  blocks meet bitwise and the `Mesh2D` connection check passes by construction.
- **Elements.** The initial h refinement of column `i` is `RefinedLevels`' `(nx, ny)` (the
  source column carries `L+1` in both dimensions, yaml:86-88), giving `2^nx × 2^ny` elements
  per block, laid out block-major with `ix` fastest. All elements carry `InitialGridPoints`
  points (yaml:74) unless `orders` overrides it; `RefinedGridPoints` (per-element p
  refinement) must be empty — hp-AMR is the `amr-driver` node's business.
- **Connections.** One [`BlockConnection`](@ref) per interior lattice face, both directions
  of the lattice: `(nr−1)·nz + nr·(nz−1)` in total.
- **Metadata.** Per-block vtu region (from [`vtu_domain`](@ref) at the column midpoint), the
  worldtube block (the cell containing the particle at `(r₀, equator)`), and the class of
  every domain-boundary face (`:horizon`, `:outer`, `:equator`, `:pole`).

Type-generic: the geometry type is the config's `T` (a `RunConfig{BigFloat}` yields a
`BlockDomain{BigFloat}` whose bounds and maps carry the ambient precision).
"""
function block_domain(orbit::CircularOrbit, config::RunConfig{T};
                      orders::NTuple{2,Integer} = config.domain.initial_grid_points) where {T}
    bounds_r, bounds_z, chart = block_lattice_bounds(config)
    nr, nz = length(bounds_r) - 1, length(bounds_z) - 1
    (nr >= 1 && nz >= 1) ||
        throw(ArgumentError("block_domain: the lattice needs at least one block per dimension"))
    isempty(config.domain.refined_grid_points) || throw(ArgumentError(
        "block_domain: RefinedGridPoints $(config.domain.refined_grid_points) is not empty — " *
        "per-element p refinement belongs to the AMR driver, not the initial lattice"))
    isempty(config.domain.blocks_to_exclude) || throw(ArgumentError(
        "block_domain: BlocksToExclude $(config.domain.blocks_to_exclude) is not supported — " *
        "the CRP lattice is a full rectangle"))
    # Distributions[r] (yaml:66): Linear in the v/t columns, Inverse in the wavezone
    dr = config.domain.distributions_r
    length(dr) == nr && all(==(:Linear), dr[1:nr-1]) && dr[nr] === :Inverse ||
        throw(ArgumentError("block_domain: Distributions[r] must be " *
                            "[Linear ×$(nr-1), Inverse] (yaml:66); got $dr"))
    all(==(:Linear), config.domain.distributions_z2) ||
        throw(ArgumentError("block_domain: the angular Distributions must all be Linear " *
                            "(yaml:68); got $(config.domain.distributions_z2)"))

    v = vtu_slicing(orbit, config)                      # validates orbit ≡ config point
    _require_horizon_inner_bound(orbit, bounds_r)
    regs = _column_regions(v, config, bounds_r)
    colref = _column_refinements(config, nr, nz)

    block_ids = [i + nr * j for i in 0:nr-1, j in 0:nz-1]
    blocks = Block[]
    regions = Dict{Int,Symbol}()
    refinements = Dict{Int,NTuple{2,Int}}()
    for j in 0:nz-1, i in 0:nr-1
        xmap = i == nr - 1 ? InverseRadialMap1D(bounds_r[i+1], bounds_r[i+2]) :
                             LinearMap1D(bounds_r[i+1], bounds_r[i+2])
        ymap = LinearMap1D(bounds_z[j+1], bounds_z[j+2])
        bid = block_ids[i+1, j+1]
        push!(blocks, Block(bid, xmap, ymap))
        regions[bid] = regs[i+1]
        refinements[bid] = colref[i+1]
    end

    elements = ElementId[]
    for j in 0:nz-1, i in 0:nr-1
        nx, ny = colref[i+1]
        bid = block_ids[i+1, j+1]
        for iy in 0:(1<<ny)-1, ix in 0:(1<<nx)-1
            push!(elements, ElementId(bid, SegmentId(nx, ix), SegmentId(ny, iy)))
        end
    end

    connections = BlockConnection[]
    for j in 0:nz-1, i in 0:nr-2
        push!(connections, BlockConnection(block_ids[i+1, j+1], 1, block_ids[i+2, j+1]))
    end
    for j in 0:nz-2, i in 0:nr-1
        push!(connections, BlockConnection(block_ids[i+1, j+1], 2, block_ids[i+1, j+2]))
    end

    mesh = Mesh2D(blocks, elements, (Int(orders[1]), Int(orders[2]));
                  connections = connections)

    # the regularized region: the lattice cell holding the particle at (r₀, equator)
    icol = findall(i -> bounds_r[i] < orbit.r0 < bounds_r[i+1], 1:nr)
    length(icol) == 1 || throw(ArgumentError(
        "block_domain: the orbital radius r₀ = $(orbit.r0) is not strictly inside exactly " *
        "one radial block of $(bounds_r)"))
    jrow = findall(j -> bounds_z[j] <= 0 <= bounds_z[j+1], 1:nz)
    length(jrow) == 1 || throw(ArgumentError(
        "block_domain: the equator is not inside exactly one angular block of $(bounds_z)"))
    worldtube = [block_ids[only(icol), only(jrow)]]

    # boundary-face classification. The lattice is a full rectangle with EVERY interior block
    # face connected, so a face of `mesh` is a domain boundary exactly when it lies on the
    # lattice perimeter, and its class follows from the face's dimension and side.
    bclasses = Dict{Tuple{ElementId,Int,Int},Symbol}()
    for f in boundary_faces(mesh)
        _, d, sidx = f
        bclasses[f] = boundary_class_rule(chart, d, sidx)
    end

    RT = typeof(mesh).parameters[1]
    return BlockDomain{RT}(mesh, nr, nz, chart, collect(RT, bounds_r), collect(RT, bounds_z),
                           block_ids, regions, refinements, worldtube, bclasses)
end

# The `:horizon` precondition. The CRP domain starts AT the outer horizon (main.tex:529,
# fig:domain), where the SEALED radial flux coefficient A^rr = Δ/(r²+a²) is exactly zero
# because Δ is carried as (r − r₊)(r − r₋) off the orbit's cached horizons — so F^r = A^rr ∂_rΨ
# carries no data through the inner radial faces and they need no boundary condition
# (main.tex:385). That is a property of WHERE the bound is placed, not of the face: at
# r_in = 2 with a = 0.5 (r₊ = 1.8660254037844386) the same coefficient measures 5.88e-2, at
# a = −0.7 it measures 1.09e-1, and INSIDE the horizon (a = 0.9, r_in = 1.3) it is −4.0e-2 —
# the operator is not even elliptic there. Enforcing the placement here hands node
# `dg-boundary-conditions` the FACT that `:horizon` faces need no data, instead of a bare tag.
# Tolerance: 8 ulp relative (the bound must BE the horizon, not merely near it); the measured
# A^rr is reported in the error so a caller sees the size of what it would have to supply.
function _require_horizon_inner_bound(orbit::CircularOrbit, bounds_r::AbstractVector)
    rin = bounds_r[1]
    RT = promote_type(typeof(rin), typeof(orbit.rp))
    rtol = 8 * eps(one(RT))
    isapprox(rin, orbit.rp; rtol = rtol) || throw(ArgumentError(
        "block_domain: the inner radial block bound must be the outer horizon " *
        "r₊ = $(orbit.rp) (main.tex:529; enforced to $(rtol) relative), got $rin — at which " *
        "the sealed radial flux coefficient is A^rr = $(flux_coeff_r(orbit, rin)) ≠ 0, so " *
        "the :horizon faces would carry flux and need boundary data this domain does not " *
        "supply (main.tex:385 holds only AT r₊)"))
    return nothing
end

"""
    block_lattice_index(dom::BlockDomain, id::Integer) -> (i, j)

The 0-based lattice cell `(radial column, angular row)` of block `id` — the inverse of the
yaml block numbering `id = i + nr·j` (radial index fastest). Throws `ArgumentError` for an id
outside the lattice.
"""
function block_lattice_index(dom::BlockDomain, id::Integer)
    0 <= id < dom.nr * dom.nz || throw(ArgumentError(
        "block_lattice_index: block id $id is outside the $(dom.nr)×$(dom.nz) lattice"))
    return (Int(id) % dom.nr, Int(id) ÷ dom.nr)
end

"""
    block_region(dom::BlockDomain, id::Integer) -> Symbol

The vtu slicing region of block `id`: `:v` (horizon-penetrating wavezone, boost H = −1), `:t`
(the Boyer–Lindquist transition/worldtube region, H = 0) or `:u` (outer wavezone, H = +1),
as classified by the sealed [`vtu_domain`](@ref) at the block's radial midpoint. The v and u
blocks are exactly the yaml's `NullSlicingBlocks` (yaml:39-45).
"""
block_region(dom::BlockDomain, id::Integer) = dom.regions[Int(id)]

"""
    worldtube_blocks(dom::BlockDomain) -> Vector{Int}

The regularized-region ("worldtube") block ids — the lattice cell containing the particle at
`(r₀, equator)`, the red rectangle of fig:domain (main.tex:529) inside which the solver
carries `Ψ^R_m = Ψ_m − Ψ^P_m` (eq:jump_worldtube). Node `worldtube-formulation` imposes the
jump on the interfaces between these blocks and the rest of the domain.
"""
worldtube_blocks(dom::BlockDomain) = dom.worldtube

"""
    boundary_face_classes(dom::BlockDomain) -> Dict{Tuple{ElementId,Int,Int},Symbol}

Every domain-boundary face of `dom.mesh` (the keys are exactly
[`boundary_faces`](@ref)`(dom.mesh)`) tagged with the physical boundary it lies on:

- `:horizon` — the inner radial boundary `r = BlockBounds[r][1]`, which
  [`block_domain`](@ref) REQUIRES to be the outer horizon `r₊` (main.tex:529). No boundary
  data are needed there — but only because of that placement: `A^rr = Δ/(r²+a²)` is exactly
  zero at `r₊`, so `F^r` carries nothing through the face (main.tex:385). Away from `r₊` the
  coefficient is nonzero (5.88e-2 at `r_in = 2`, `a = 0.5`) and negative inside the horizon,
  so the precondition is enforced at construction rather than assumed here.
- `:outer` — the outer radial boundary `r = R`, where the order-2 Bayliss–Turkel/Sommerfeld
  radiation condition lives for `m ≠ 0` (`None` at `m = 0`, yaml:115-125).
- `:equator` — `z² = 0`, present only in the equatorial-symmetry chart; the angular flux
  coefficient vanishes there in the even sector, so no boundary data are needed.
- `:pole` — `z² = 1` (or `z = ±1` in the full chart); `F^z(±1) = 0` exactly (main.tex:385).

Node `dg-boundary-conditions` consumes this table to build the `MeshOperator` `bcs` entry of
every boundary face.
"""
boundary_face_classes(dom::BlockDomain) = dom.boundary_classes

"""
    boundary_class(dom::BlockDomain, face::Tuple{ElementId,Int,Int}) -> Symbol

The class of one domain-boundary `face = (element, dim, sidx)` (`sidx` 1 ↔ side −1, 2 ↔ side
+1, the EllipticCore face key). See [`boundary_face_classes`](@ref); a face that is not a
domain boundary raises `KeyError`.
"""
boundary_class(dom::BlockDomain, face::Tuple{ElementId,Int,Int}) = dom.boundary_classes[face]

"""
    block_coefficients(dom::BlockDomain, orbit, config, id::Integer)
        -> (A = (A¹¹, A²²), beta, gamma = (γ¹, γ²))

The sealed first-order flux-form coefficient bundle of block `id`: exactly
[`fluxform_callables`](@ref)`(orbit, config, block_region(dom, id))` — the per-region assembly
with the block's constant vtu boost `H ∈ {−1, 0, +1}` baked in, in the `(r, z)` or `(r, z²)`
chart the config selects, in the shape EllipticCore's `element_coefficients` consumes. No new
algebra: this function only chooses the region.
"""
block_coefficients(dom::BlockDomain, orbit::CircularOrbit, config::RunConfig, id::Integer) =
    fluxform_callables(orbit, config, block_region(dom, id))

"""
    domain_coefficients(dom::BlockDomain, orbit, config)
        -> (A = (A¹¹, A²²), beta, gamma = (γ¹, γ²))

ONE coefficient bundle for the WHOLE domain — what a single [`MeshOperator`](@ref) over
`dom.mesh` needs, since it takes one set of coefficient callables for all elements. The
bundle dispatches on the radius through the sealed [`vtu_domain`](@ref) classifier, and it is
BITWISE equal to the per-block bundles ([`block_coefficients`](@ref)) at every point the DG
operator ever samples:

- `A` carries no slicing boost at all (`flux_coeff_r`/`flux_coeff_z`/`flux_coeff_z2` do not
  take `H`), so its region has no meaning — the same callables serve every block. This is what
  makes the global bundle safe: `A` is the ONLY coefficient `element_coefficients` evaluates
  ON a face, i.e. at a slicing transition radius, where the pointwise classifier necessarily
  returns the closed-`t` value and a per-region boost would be ambiguous.
- `beta` and `gamma[1]` do depend on the boost and are evaluated only at the interior
  Legendre–Gauss volume nodes, which are strictly inside a block and hence strictly inside one
  slicing slab — so the radius dispatch reproduces the block's own region exactly.
- `gamma[2]` is boost-free as well (`gamma_z_m`/`gamma_z2_m`).

The whole equality is gated element-by-element on the production lattice in
test_domain_blocks.jl.
"""
function domain_coefficients(dom::BlockDomain, orbit::CircularOrbit, config::RunConfig)
    (config.orbit.impose_equatorial_symmetry ? :z2 : :z) === dom.chart || throw(ArgumentError(
        "domain_coefficients: config chart $(config.orbit.impose_equatorial_symmetry ?
         :z2 : :z) does not match the domain chart $(dom.chart)"))
    v = vtu_slicing(orbit, config)
    fv = fluxform_callables(orbit, config, :v)
    ft = fluxform_callables(orbit, config, :t)
    fu = fluxform_callables(orbit, config, :u)
    pick(r) = (d = vtu_domain(v, r); d === :v ? fv : d === :t ? ft : fu)
    return (A = ft.A,
            beta = (r, y) -> pick(r).beta(r, y),
            gamma = ((r, y) -> pick(r).gamma[1](r, y), ft.gamma[2]))
end

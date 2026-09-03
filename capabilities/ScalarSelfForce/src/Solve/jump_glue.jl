# Solve/jump_glue.jl — dg-jump-conditions (ESSF node, schedule W9-A, GLOBAL_DAG row 40): the
# binding of the CRP's two INTERIOR jump conditions to the EllipticCore ModifyBoundaryData
# seam, i.e. the `MeshOperator` `hooks` table of the block domain of arXiv:2606.04998
# (main.tex:596-608, §III.A.5 "Imposing jump conditions").
#
# The paper imposes both jumps by modifying the RECEIVED exterior data, never the assembly:
#
#   * WORLDTUBE (eq:jump_worldtube, main.tex:600-604)
#         Ψ^ext  →  Ψ^ext  ± Ψ^𝒫 ,     (n_i𝓕^i)_ext  →  (n_i𝓕^i)_ext ∓ n_i𝔄^{ij}∂_jΨ^𝒫 ,
#     TRIAL-FIELD INDEPENDENT ⇒ a [`FixedDataHook`](@ref) per mortar: it is dropped from the
#     linearized application and reaches the right-hand side through `operator_on_zero`
#     (spectre `ModifyBoundaryData::apply`, Equations.cpp:70-91, harvested in
#     DgOperator.hpp:1081-1130). The DATA are the sealed row-38 [`worldtube_jump_data`](@ref).
#   * VTU JUNCTIONS (eq:jump_vt/eq:jump_tu, main.tex:337-343, realized main.tex:606-608)
#         (n_i𝓕^i)_ext  →  (n_i𝓕^i)_ext − imΩ ⟨⟨Ψ_m⟩⟩ ,   ⟨⟨Ψ_m⟩⟩ = (Ψ^int + Ψ^ext)/2 ,
#     the SAME sign on both sides of the boundary and at both transition radii ("the change
#     in sign on either side is absorbed by the opposite face normal", main.tex:608).
#     LINEAR-HOMOGENEOUS in the trial field ⇒ a [`LinearizedHook`](@ref): the term is IN the
#     GMRES/Schwarz/multigrid matrix (spectre `ModifyBoundaryData::apply_linearized`,
#     Equations.cpp:92-116, applied in DgOperator.hpp:708-728). That asymmetry — vtu in the
#     matrix, worldtube in the fixed sources — is this node's oracle.
#
# ZERO new physics and ZERO re-derived coefficients. In particular the junction datum needs
# NO chart conversion: the sealed `A^rr = Δ/(r²+a²)` IS `dr/dr_*`, so `𝓕^r = A^rr ∂_rΨ` IS
# `∂_{r_*}Ψ` (main.tex:354-358; the [`VtuJunction`](@ref) docstring states it), and the whole
# datum is the sealed `vtu_junctions(v, orbit, m).coeff = −imΩ`. This file decides only WHICH
# faces carry WHICH sealed hook, and composes the two when a face carries both.
#
# WHY `Solve/` AND NOT `Worldtube/`. This is the `hooks` twin of `Solve/bc_glue.jl`'s `bcs`
# table: it is solver wiring that composes THREE sealed sources — the worldtube data
# (`Worldtube/worldtube.jl`), the junction data (`MModeOperator/vtu.jl`) and the block lattice
# (`Domain/blocks.jl`) — into one EllipticCore table, and `worldtube.jl` is scoped by its own
# header to "NO vtu-junction wiring". DESIGN.md §2.2 lists row 40 on the `Worldtube` line;
# the bc_glue precedent (row 41) makes `Solve/` the symmetric home and keeps the worldtube
# module a pure-data node.
#
# SCOPE BOUNDARIES (predecessors own the substance; later nodes own the rest):
#   - Worldtube/worldtube.jl (row 38) — `worldtube_jump_data` (one `FaceState` PER MORTAR in
#     the receiver's mortar order), `worldtube_jump_signs`, the fixed source;
#   - MModeOperator/vtu.jl (row 28) — `vtu_slicing`, `vtu_junctions`, `VtuJunction`;
#   - Domain/blocks.jl (row 39) — `block_domain`, `block_region` (the yaml `NullSlicingBlocks`
#     classification this file's junction rule tests), `domain_coefficients`;
#   - EllipticCore matrixfree.jl (rows 15 + 40 commit 1) — the per-MORTAR `hooks` seam,
#     `mortar_points`, `MeshInterface`; boundary_hooks.jl (row 14) — the hook types;
#   - Solve/bc_glue.jl (row 41) — `boundary_conditions`, the `bcs` table;
#   - solver-stack-ssf (row 42) owns the driver: the RHS b = M.S - operator_on_zero(op)
#     with S the SWITCHED volume evaluator (see the jump_hooks docstring), and the solve;
#     amr-driver (row 43) rebuilds this table per level (recipe in the same docstring).

export vtu_junction_faces, vtu_junction_hook, jump_hooks

# The junction modification as a NAMED callable (rather than a closure) so the hook is
# type-stable, printable, and identifiable in an assembled table.
struct VtuJunctionFlux{T<:Number}
    coeff::T
end
(f::VtuJunctionFlux)(int::FaceState, ext::FaceState) =
    FaceState(ext.psi, ext.nflux .+ f.coeff .* (int.psi .+ ext.psi) ./ 2)

"""
    vtu_junction_hook(j::VtuJunction) -> LinearizedHook
    vtu_junction_hook(coeff::Number)  -> LinearizedHook

The vtu slicing-junction condition as an EllipticCore [`LinearizedHook`](@ref): the received
exterior state of a face lying on a transition radius is modified by

    (n·𝓕)^ext  →  (n·𝓕)^ext + coeff · (Ψ^int + Ψ^ext)/2 ,      coeff = −imΩ ,

with `Ψ^ext` untouched (the mode function is continuous across the interface, main.tex:335).
This is eq:jump_vt/eq:jump_tu `[∂_{r_*}Ψ_m] = −imΩ Ψ_m` in the DG realization the paper
prints (main.tex:606-608) and spectre implements verbatim
(`ScalarSelfForce::ModifyBoundaryData::apply_linearized`, Equations.cpp:105-107:
`n_dot_flux_remote -= i·m·Ω·(field_local + field_remote)/2`). Three properties make it a
`LinearizedHook` rather than a `FixedDataHook`:

- it is LINEAR-HOMOGENEOUS in the trace data, so it acts in EVERY operator application and
  belongs to the GMRES matrix (`build_matrix` changes when it is added — the row-40 oracle);
- it contributes NOTHING to `operator_on_zero` (zero trial data ⇒ zero shift), so it never
  reaches the right-hand side;
- the SAME hook serves both sides of the interface and both transition radii — the sign
  change is absorbed by the opposite face normal (main.tex:608; spectre: "the signs are all
  the same"). No orientation bookkeeping exists here to get wrong.

NO chart conversion is applied: the sealed radial flux coefficient `A^rr = Δ/(r²+a²)` IS
`dr/dr_*`, so the DG's `φ = n·𝓕 = ±A^rr∂_rΨ` already IS `±∂_{r_*}Ψ` (main.tex:354-358, the
[`VtuJunction`](@ref) docstring), and `coeff` is the printed `−imΩ` unmodified. At `m = 0`
the coefficient is exactly zero and the hook is the identity (gated bitwise on
`build_matrix`).

Consistency with the sealed prescribed-jump table (EllipticCore boundary_hooks.jl): imposing
`[Ψ] = J`, `[𝓕] = K` needs `(Δψ, Δφ) = (∓J, +K)`; here `J = 0` and `K = coeff·⟨⟨Ψ⟩⟩`, giving
exactly the `Δφ` above on both receivers.
"""
vtu_junction_hook(coeff::Number) = LinearizedHook(VtuJunctionFlux(coeff))
vtu_junction_hook(j::VtuJunction) = vtu_junction_hook(j.coeff)

# The vtu region membership spectre keys the junction on: `null_slicing_blocks` is exactly the
# set of v and u blocks (yaml:39-45), i.e. the blocks whose boost H ≠ 0. A face is a junction
# face iff its owner and its neighbour(s) disagree on this predicate
# (Equations.cpp:97-101) — which is equivalent to, and cheaper and more robust than, testing
# the face radius against the transition radii, because block ids survive h-coarsening and AMR
# while element ids do not.
_jg_null_sliced(dom::BlockDomain, block::Integer) = block_region(dom, block) !== :t

# mortar count of one INTERIOR face of `mesh`, with the mesh-boundness diagnostics this node's
# consumers need (a `WorldtubeData` built on another level's mesh is the likely mistake).
function _jg_nmortars(mesh::Mesh2D, face::Tuple{ElementId,Integer,Integer}, who::AbstractString)
    el, dd, ss = face
    d, sidx = Int(dd), Int(ss)
    (d == 1 || d == 2) || throw(ArgumentError("$who: face $face has dim $d, expected 1 or 2"))
    (sidx == 1 || sidx == 2) ||
        throw(ArgumentError("$who: face $face has sidx $sidx, expected 1 or 2"))
    haskey(mesh.index, el) || throw(ArgumentError(
        "$who: $el is not an element of this mesh — jump data are MESH-BOUND (element ids " *
        "change under coarsening and AMR), so a WorldtubeData must be rebuilt for the mesh " *
        "the operator is built on"))
    iface = mesh.faces[mesh.index[el]][d][sidx]
    iface === nothing && throw(ArgumentError(
        "$who: face $face is a DOMAIN boundary, not an interior interface — external " *
        "conditions go in the `bcs` table (boundary_conditions), not in `hooks`"))
    return length(iface.nbs)
end

"""
    vtu_junction_faces(dom::BlockDomain, orbit::CircularOrbit, config::RunConfig,
                       m = config.orbit.m_mode_number; mesh = dom.mesh)
        -> Dict{Tuple{ElementId,Int,Int},VtuJunction}

Every interior face of `mesh` that lies on a vtu slicing transition, mapped to the sealed
[`VtuJunction`](@ref) datum it carries — the `vt` junction at `r_v` (`bounds_r[2]`) between
the v and t columns and the `tu` junction at `r_u` (`bounds_r[5]`) between t and u
(eq:jump_vt/eq:jump_tu). BOTH receivers of each interface appear: each is a receiver of
exterior data and needs its own hook.

The rule is spectre's, not a radius comparison: a face is a junction face exactly when its
owner block and its neighbour block(s) DISAGREE about `NullSlicingBlocks` membership, i.e.
about `block_region(dom, ·) !== :t` (`Equations.cpp:97-101`). Block ids survive h-coarsening
and AMR while element ids do not, so the same rule serves every level mesh. As a cross-check
the face's radial coordinate ([`mortar_points`](@ref)) must EQUAL one of the two transition
radii — exactly, since both come from the same `RunConfig` numbers and block-map endpoints
are exact by the `Mesh2D` contract; a face that crosses the slicing without sitting on a
transition radius is an `ArgumentError` rather than a silent mis-assignment. A face that
straddles the slicing only PARTIALLY (some mortars crossing, some not) is likewise refused —
it cannot happen on a block lattice, where a face's mortars all abut one neighbour block.

`mesh` defaults to `dom.mesh` but may be any mesh built from `dom`'s block lattice — a
coarser [`coarsen_mesh`](@ref) level or an AMR-refined one (the per-level `build_op` contract
of [`MultigridHierarchy`](@ref), the same freedom [`boundary_conditions`](@ref) has). `m`
must equal the config's `MModeNumber`: the junction coefficient is `−imΩ` and the operator's
own coefficients come from the same config.

At `m = 0` every coefficient is exactly zero (the slicing degenerates — `h` drops out of the
phase map), so the table is still returned and its hooks are the identity.
"""
function vtu_junction_faces(dom::BlockDomain, orbit::CircularOrbit, config::RunConfig,
                            m::Integer = config.orbit.m_mode_number;
                            mesh::Mesh2D = dom.mesh)
    Int(m) == config.orbit.m_mode_number || throw(ArgumentError(
        "vtu_junction_faces: m = $m disagrees with the config's MModeNumber " *
        "$(config.orbit.m_mode_number) — the junction coefficient −imΩ and the operator " *
        "coefficients must belong to the same mode (use crp_paper_config(m_mode = $m))"))
    v = vtu_slicing(orbit, config)                 # validates orbit ≡ config point
    jc = vtu_junctions(v, orbit, Int(m))
    out = Dict{Tuple{ElementId,Int,Int},typeof(jc.vt)}()
    for e in eachindex(mesh.elements), d in 1:2, sidx in 1:2
        iface = mesh.faces[e][d][sidx]
        iface === nothing && continue
        el = mesh.elements[e]
        here = _jg_null_sliced(dom, el.block)
        ncross = count(k -> _jg_null_sliced(dom, mesh.elements[iface.nbs[k]].block) != here,
                       1:length(iface.nbs))
        ncross == 0 && continue
        ncross == length(iface.nbs) || throw(ArgumentError(
            "vtu_junction_faces: face $((el, d, sidx)) crosses the vtu slicing on only " *
            "$ncross of its $(length(iface.nbs)) mortars — the slicing regions are whole " *
            "block columns (main.tex:323), so this mesh is not built from $(dom.nr)×$(dom.nz) " *
            "lattice blocks"))
        xf, _ = mortar_points(mesh, (el, d, sidx))
        out[(el, d, sidx)] = xf == jc.vt.r ? jc.vt : xf == jc.tu.r ? jc.tu :
            throw(ArgumentError(
                "vtu_junction_faces: face $((el, d, sidx)) crosses the vtu slicing at " *
                "coordinate $xf, which is neither transition radius (r_v = $(jc.vt.r), " *
                "r_u = $(jc.tu.r)) — the transitions must BE radial block bounds " *
                "(main.tex:323, enforced by block_domain)"))
    end
    return out
end

"""
    jump_hooks(mesh::Mesh2D, jump_data::AbstractDict, junctions::AbstractDict)
        -> Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}
    jump_hooks(dom::BlockDomain, wt::WorldtubeData, orbit::CircularOrbit, config::RunConfig,
               m = config.orbit.m_mode_number; mesh = dom.mesh) -> Dict{…}

The `MeshOperator(…; hooks = …)` table of the CRP problem: both jump conditions of
main.tex:596-608 bound onto the EllipticCore seam, PER MORTAR.

Every entry is a `Vector{BoundaryDataHook}` with one hook per mortar of that face, in the
receiver's mortar order ([`mortar_points`](@ref) / `MeshInterface`) — one uniform typed
representation, no per-face/per-mortar flag. That granularity is the point of this node: the
production lattice pre-refines the worldtube column one level above its left neighbour
(yaml:83-88), so `2^L` of its worldtube faces are 2:1 with genuinely different `Ψ^𝒫` on their
two halves, which the pre-row-40 per-FACE seam could not express at all.

Composition per face:

| face carries | hook |
|---|---|
| worldtube only | `FixedDataHook(Δ_k)` |
| junction only | `vtu_junction_hook(j)` (the same object on every mortar) |
| both | `ComposedHook(vtu_junction_hook(j), FixedDataHook(Δ_k))` |

The composed order is the sealed W5-B one: the LINEARIZED stage acts on the RAW received data
FIRST, then the fixed shift (spectre's vtu convention, Equations.cpp:105-107 — `⟨⟨Ψ⟩⟩` is the
average of the untransformed traces). On the CRP lattice the third row never occurs: the
worldtube block spans `[r₀−w, r₀+w] = [5.2, 14.8]` while the transitions sit at
`r₀−1.5w = 2.8` and `r₀+3w = 24.4` (main.tex:545-548 — the worldtube column lies strictly
inside the t slab), which the node tests assert on the production lattice; the composition is
built generically anyway so an AMR or a re-tuned lattice cannot silently drop a condition.

`jump_data` is the sealed [`worldtube_jump_data`](@ref) — face key ↦ one `FaceState` per
mortar — and `junctions` the [`vtu_junction_faces`](@ref) table; both are validated against
`mesh` (a face that is not an interior face of `mesh`, or whose datum has the wrong number of
mortars, is an `ArgumentError`, because jump data are MESH-BOUND and a stale `WorldtubeData`
is the likely mistake). The `BlockDomain` method builds `junctions` itself and reads the
worldtube data off `wt`; because `WorldtubeData` is mesh-bound, a caller passing a non-default
`mesh` (a multigrid level, an AMR level) must have REBUILT `wt` for that mesh — the generic
`worldtube_jump_data(mesh, isinside, values)` method is the seam for that, and the mismatch
is caught rather than mis-applied.

PER-LEVEL RECIPE (the multigrid / AMR contract, stated ONCE, here). Both tables are functions
of the level mesh; the `WorldtubeData` is MESH-BOUND and is rebuilt with the same `mesh` kwarg:

    wt_l    = worldtube_data(dom, orbit, config, m; mesh = mesh_l)
    hooks_l = jump_hooks(dom, wt_l, orbit, config, m; mesh = mesh_l)
    bcs_l   = boundary_conditions(dom, orbit, config; mesh = mesh_l)

The junction half needs no rebuild logic at all — it is recomputed from the block lattice, and
block ids survive coarsening and AMR while element ids do not. A stale `wt` is REFUSED, not
mis-applied: the mortar count AND the per-mortar point count of every datum are checked against
`mesh` (h- and p-refinement respectively), and the `BlockDomain` method additionally requires
`wt.dom === dom`, because block ids are lattice-local and a domain at another spin, radius or
chart carries the same ids.

RIGHT-HAND SIDE (node `solver-stack-ssf` owns it; recorded here so the two halves cannot
drift): `b = M·S - operator_on_zero(op)` with the LINEARIZED [`apply!`](@ref) as the matrix —
into which the junction hooks enter and the worldtube hooks do not. The VOLUME fixed source
`S` handed to that assembly is **row 42's switched evaluator**
([`effective_source_m_switched`](@ref)), NOT `wt.source`: the direct route's roundoff grows
like λ^-2 and is unusable below λ ≈ 1e-2 at Float64, while the AMR driver refines AT the
puncture eight times. Assemble it through the generic
`worldtube_source(mesh, isinside, source, T)` / `worldtube_massive_source(...)` seam —
switched inside the tube, exactly zero outside. The MORTAR data this file binds are a
different object and stay on the DIRECT [`effective_source_bundle`](@ref): Ψ^𝒫 and its two
derivatives are trace data on the tube boundary, far from the worldline, and are never
switched.
"""
function jump_hooks(mesh::Mesh2D, jump_data::AbstractDict, junctions::AbstractDict)
    out = Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}()
    for (face0, ms) in jump_data
        face = (face0[1], Int(face0[2]), Int(face0[3]))
        nm = _jg_nmortars(mesh, face, "jump_hooks")
        length(ms) == nm || throw(ArgumentError(
            "jump_hooks: the worldtube jump datum of face $face has $(length(ms)) mortar " *
            "states but the face has $nm mortar(s) on this mesh — jump data are MESH-BOUND " *
            "(worldtube_jump_data returns one FaceState per mortar in the receiver's mortar " *
            "order); rebuild the WorldtubeData for this mesh " *
            "(worldtube_data(dom, orbit, config, m; mesh = mesh))"))
        # …and each mortar's datum must have the mortar's POINT COUNT. h-topology alone does
        # not pin it: a p-refined mesh (node amr-driver raises the order EVERY step) has the
        # identical face/mortar structure with more points per mortar, so a stale datum would
        # pass the count above and die inside the flux with a raw DimensionMismatch.
        mp = mortar_points(mesh, face)[2]
        for k in 1:nm
            length(ms[k].psi) == length(mp[k]) || throw(ArgumentError(
                "jump_hooks: mortar $k of face $face carries $(length(ms[k].psi)) points but " *
                "that mortar has $(length(mp[k])) on this mesh — the jump data were built on " *
                "a mesh of different polynomial order; rebuild the WorldtubeData for this " *
                "mesh (worldtube_data(dom, orbit, config, m; mesh = mesh))"))
        end
        jn = get(junctions, face, nothing)
        out[face] = if jn === nothing
            BoundaryDataHook[FixedDataHook(s) for s in ms]
        else
            vh = vtu_junction_hook(jn)             # W5-B: linearized stage FIRST
            BoundaryDataHook[ComposedHook(vh, FixedDataHook(s)) for s in ms]
        end
    end
    for (face0, jn) in junctions
        face = (face0[1], Int(face0[2]), Int(face0[3]))
        haskey(out, face) && continue
        nm = _jg_nmortars(mesh, face, "jump_hooks")
        vh = vtu_junction_hook(jn)
        out[face] = BoundaryDataHook[vh for _ in 1:nm]
    end
    return out
end

function jump_hooks(dom::BlockDomain, wt::WorldtubeData, orbit::CircularOrbit,
                    config::RunConfig, m::Integer = config.orbit.m_mode_number;
                    mesh::Mesh2D = dom.mesh)
    Int(m) == wt.m || throw(ArgumentError(
        "jump_hooks: m = $m disagrees with the WorldtubeData's mode $(wt.m) — the jump data " *
        "would belong to a different equation"))
    wt.dom === dom || throw(ArgumentError(
        "jump_hooks: the WorldtubeData was built on a DIFFERENT BlockDomain object. Block " *
        "ids are lattice-local, so a domain at another spin, radius, chart or refinement " *
        "level has the same ids and would be accepted silently; identity is required"))
    junctions = vtu_junction_faces(dom, orbit, config, m; mesh = mesh)
    return jump_hooks(mesh, worldtube_jump_data(wt), junctions)
end

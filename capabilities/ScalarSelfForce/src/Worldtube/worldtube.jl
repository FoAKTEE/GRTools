# Worldtube/worldtube.jl — the WORLDTUBE FORMULATION (ESSF node `worldtube-formulation`,
# schedule W8-C, GLOBAL_DAG row 38) of arXiv:2606.04998 §II.D as PURE DATA for the DG layer.
#
# The paper solves for the REGULAR field inside a rectangular worldtube around the particle
# and for the RETARDED field outside (main.tex:513-523, eq:jump_worldtube):
#
#     △_m Ψ^ℛ_m = S^eff_m   (interior),      △_m Ψ_m = 0   (exterior),
#     Ψ^ℛ_m     = Ψ_m − Ψ^𝒫_m               (across the worldtube interface).
#
# Both statements are imposed WITHOUT touching the DG assembly (main.tex:596-608): the
# interface condition by modifying the RECEIVED exterior trace data on every element that
# receives across the worldtube boundary,
#
#     Ψ^ext_m        →  Ψ^ext_m        ±  Ψ^𝒫_m ,
#     (n_i𝓕^i_m)_ext →  (n_i𝓕^i_m)_ext ∓  n_i 𝔄^{ij}_m ∂_j Ψ^𝒫_m ,
#
# and the interior equation by a FIXED SOURCE +S^eff_m at every quadrature node inside the
# tube, zero outside (spectre `Tags::FixedSource<Tags::MMode>`).
#
# THIS FILE IS DATA ONLY. It owns four things and nothing else:
#   (1) the regularized-region partition of a [`BlockDomain`](@ref) — which elements are
#       inside, and which interior faces form the worldtube boundary;
#   (2) the per-face, PER-MORTAR jump data (the ±Ψ^𝒫 / ∓n·𝓕(Ψ^𝒫) `FaceState`s) at the
#       points the EllipticCore mortar exchange actually evaluates;
#   (3) the fixed-source vector (+`effective_source_m` inside, exactly zero outside) in the
#       nodal and the massive `M·S` form the DG residual `R(Ψ) = M·S` consumes;
#   (4) the sign convention itself, as ONE normative statement ([`worldtube_jump_signs`](@ref)).
# There is NO solve, NO operator construction, and NO vtu-junction wiring here: binding the
# data into `MeshOperator`'s `hooks`/`bcs` tables together with the vtu junction hooks is node
# `dg-jump-conditions`, and the RHS composition `b = M·S − operator_on_zero(op)` belongs to the
# solver driver (`solver-stack-ssf`).
#
# ZERO new physics: every value comes from the sealed [`effective_source_bundle`](@ref) (ONE
# call per evaluation point, supplying Ψ^𝒫_m, ∂_rΨ^𝒫_m, ∂_zΨ^𝒫_m and S^eff_m together) closed
# with the sealed flux coefficients [`flux_coeff_r`](@ref) / [`flux_coeff_z`](@ref) /
# [`flux_coeff_z2`](@ref). The z² equatorial chart uses the PRINTED chain rule of
# main.tex:389-394 (𝓕^{z²} = A^{z²}∂_wΨ with ∂_w = ∂_z/(2z)), exactly as `z2.jl` states it —
# no coefficient is re-derived here.
#
# COUPLING NOTE (deliberate, gated): the per-mortar geometry of a face
# ([`worldtube_face_points`](@ref)) is reconstructed from `Mesh2D`'s precomputed face topology
# (`mesh.faces[e][d][sidx]`) with the same arithmetic `MeshOperator` uses to place its mortars
# (EllipticCore `matrixfree.jl`), so the data land exactly on the points the hook sees. The
# manufactured gates of test_worldtube.jl run the assembled data THROUGH `MeshOperator`, which
# is what pins that agreement.

export WorldtubeData, worldtube_data, inside_worldtube, worldtube_elements,
       worldtube_boundary_faces, worldtube_face_points, worldtube_jump_signs,
       worldtube_jump_values, worldtube_jump_data, worldtube_hooks,
       worldtube_source, worldtube_massive_source

# ============================== the sign convention ========================================

"""
    worldtube_jump_signs(inside::Bool, sidx::Integer) -> (psi = ±1, flux = ±1)

**The normative statement of the worldtube jump signs** (eq:jump_worldtube realized as
main.tex:600-608; spectre `ScalarSelfForce::ModifyBoundaryData::apply`, Equations.cpp:80-91).
Everything else in this package that signs worldtube data refers here.

`inside` says whether the RECEIVING element solves for the regular field `Ψ^ℛ_m` (it lies in
the worldtube), and `sidx` is the receiving face's side (1 ↔ `ξ_d = −1`, 2 ↔ `ξ_d = +1`, the
EllipticCore face key). The returned integers are the factors applied to the SEALED puncture
values at the face points to build the additive [`FixedDataHook`](@ref) shift of the RECEIVED
exterior [`FaceState`](@ref) (SENDER convention, dg_operator.jl N3):

    Δψ = psi  · Ψ^𝒫_m ,        Δφ = flux · 𝓕^d(Ψ^𝒫_m) ,

with `𝓕^d = 𝔄^{dd}∂_dΨ` the flux component along the face normal dimension.

DERIVATION (each step independently gated, see the test file):
  * `psi = −1` when the receiver is inside, `+1` when it is outside. The received datum is
    the SENDER's variable; converting it to the receiver's variable subtracts Ψ^𝒫 going into
    the tube and adds it coming out. This is spectre's
    `sign = field_is_regularized ? −1 : +1` with `field += sign·Ψ^𝒫`.
  * `flux = s_in`, the side sign (`−1` for `sidx == 1`, `+1` for `sidx == 2`) of the INSIDE
    element's face — i.e. `Δφ = n^{in}_i𝓕^i(Ψ^𝒫)`, the same value on BOTH sides of the
    interface. When the receiver is inside `s_in` is its own side sign; when it is outside the
    inside element sits on the opposite side, so `s_in = −s`. This is the paper's remark that
    "the change in sign on either side of the boundary is absorbed by the opposite face
    normal" (main.tex:608), and spectre's `n_dot_flux -= sign·(n·𝓕(Ψ^𝒫))` with the singular
    flux built on the RECEIVER's own mortar normal (`InitializeEffectiveSource.hpp:205-215`).
  * Consistency with the sealed EllipticCore prescribed-jump table (boundary_hooks.jl,
    gated at 9.5e-15): imposing `[Ψ] = J` and `[𝓕] = K` across a face requires
    `(Δψ, Δφ) = (−J, +K)` on the lower element and `(+J, +K)` on the upper one — ψ shifts
    OPPOSITELY, φ identically. With `u = Ψ^ℛ` inside and `u = Ψ` outside, `J = ±Ψ^𝒫` and
    `K = ±𝓕(Ψ^𝒫)` with the signs fixed by which side is inside, reproducing this table.

The four values are `(true, 1) → (−1, −1)`, `(true, 2) → (−1, +1)`, `(false, 1) → (+1, +1)`,
`(false, 2) → (+1, −1)`.
"""
function worldtube_jump_signs(inside::Bool, sidx::Integer)
    (sidx == 1 || sidx == 2) || throw(ArgumentError(
        "worldtube_jump_signs: sidx must be 1 (side −1) or 2 (side +1), got $sidx"))
    s = sidx == 1 ? -1 : 1
    return (psi = inside ? -1 : 1, flux = inside ? s : -s)
end

# ============================== the pointwise physics kernel ===============================

"""
    worldtube_jump_values(p::Puncture, m::Integer, chart::Symbol, x1, x2; phi_p = 0)
        -> (psi, flux_r, flux_q, source)

The four worldtube quantities at ONE point of the discretization chart, from ONE sealed
[`effective_source_bundle`](@ref) call:

- `psi`    — `Ψ^𝒫_m`, the m-mode puncture (the jump `[Ψ_m] = ±Ψ^𝒫_m` of eq:jump_worldtube);
- `flux_r` — `𝓕^r(Ψ^𝒫_m) = A^{rr}∂_rΨ^𝒫_m` with `A^{rr}` the sealed [`flux_coeff_r`](@ref);
- `flux_q` — the ANGULAR flux component in the chart's own angular coordinate;
- `source` — `S^eff_m`, the fixed source of `△_mΨ^ℛ_m = S^eff_m` inside the tube.

`chart` selects the angular coordinate, matching `BlockDomain.chart`:

- `:z`  — `x2 = z = cosθ`, `flux_q = A^{zz}∂_zΨ^𝒫_m` ([`flux_coeff_z`](@ref));
- `:z2` — `x2 = w = z²` (the equatorial-symmetry half domain), `flux_q = A^{z²}∂_wΨ^𝒫_m`
  with the PRINTED chain rule of main.tex:389-394 as `z2.jl` states it: the coefficient is
  the sealed [`flux_coeff_z2`](@ref)`(orbit, r, √w) = 4z²A^{zz}` and the derivative slot is
  reinterpreted as `∂_w = ∂_z/(2z)`. `source` is chart-independent (`△_m` in the `w` chart
  equals `△_m` in the `z` chart on even fields — the sealed row-31 identity), so it is the
  same number in both.

`phi_p` is the particle's BOYER–LINDQUIST azimuth; the horizon-regular field-point rotation is
applied INTERNALLY by the sealed bundle and must never be applied twice.

DOMAIN (loud, never a silent NaN — node `amr-driver`'s puncture snap can otherwise place a
collocation point on the worldline): the sealed `DomainError`s at `ϱ² = 0` (the worldline
circle) and `|z| = 1` (the poles) propagate, and the `:z2` chart additionally refuses `w = 0`,
where its `∂_w = ∂_z/(2z)` Jacobian is singular. `w = 0` is the equatorial DOMAIN boundary of
the half-domain chart — never a worldtube face (the tube's angular boundary sits at
`w = dz² > 0`), and the flux coefficient `A^{z²}` vanishes there identically anyway
(main.tex:399-400), so no worldtube datum is lost. An unknown `chart` is an `ArgumentError`.
NOTE the `w = 0` refusal is on the WHOLE bundle, including the chart-independent `source`
slot, which is perfectly regular there: the guard is deliberately coarse because the four
slots come from ONE evaluation. It is unreachable on this discretization — Legendre-Gauss
nodes are strictly interior, so neither a volume node nor a mortar point ever lands on
`w = 0` — but a Lobatto/endpoint quadrature rule (or an equator observer) would have to take
the regular limit `∂_wΨ = lim_{z→0} ∂_zΨ/(2z)` instead, or call
[`effective_source_m`](@ref) directly for the source alone.

Pure and allocation-free (the sealed bundle is 0-byte on the Float64 path); type-generic —
BigFloat orbits and coordinates give `Complex{BigFloat}` slots at the ambient precision.
"""
function worldtube_jump_values(p::Puncture, m::Integer, chart::Symbol, x1::Real, x2::Real;
                               phi_p::Real = 0)
    o = p.orbit
    if chart === :z
        b = effective_source_bundle(p, m, x1, x2; phi_p = phi_p)
        return (psi = b.psi,
                flux_r = flux_coeff_r(o, x1) * b.dpsi_r,
                flux_q = flux_coeff_z(o, x1, x2) * b.dpsi_z,
                source = b.source)
    elseif chart === :z2
        z = sqrt(x2)
        b = effective_source_bundle(p, m, x1, z; phi_p = phi_p)
        z > 0 || throw(DomainError(x2,
            "worldtube_jump_values: the z² chart needs w > 0 — its angular chart Jacobian \
             ∂_w = ∂_z/(2z) is singular at the equator w = 0, which is the half domain's own \
             boundary (where A^{z²} vanishes identically, main.tex:399-400) and never a \
             worldtube face"))
        return (psi = b.psi,
                flux_r = flux_coeff_r(o, x1) * b.dpsi_r,
                flux_q = flux_coeff_z2(o, x1, z) * (b.dpsi_z / (2 * z)),
                source = b.source)
    end
    throw(ArgumentError("worldtube_jump_values: unknown chart $chart (expected :z or :z2)"))
end

# ============================== the region partition =======================================

"""
    inside_worldtube(dom::BlockDomain, el::ElementId) -> Bool
    inside_worldtube(wt::WorldtubeData, el::ElementId) -> Bool

Does element `el` lie in the regularized region, i.e. does it solve for `Ψ^ℛ_m = Ψ_m − Ψ^𝒫_m`?
The region is the WHOLE block containing the particle ([`worldtube_blocks`](@ref) — the red
rectangle of fig:domain, main.tex:529), matching spectre's `puncture_in_element`, which tests
block-logical containment of the puncture (`InitializeEffectiveSource.hpp:130-140`).
"""
inside_worldtube(dom::BlockDomain, el::ElementId) = el.block in worldtube_blocks(dom)

"""
    worldtube_elements(dom::BlockDomain) -> Vector{ElementId}
    worldtube_elements(wt::WorldtubeData) -> Vector{ElementId}

The elements inside the worldtube, in the mesh's global layout order (so their slices of a
global nodal vector are `mesh.offsets`-contiguous per element). At the production lattice this
is the `2^{L+1} × 2^{L+1}` element grid of the pre-refined source column (yaml:86-88).
"""
worldtube_elements(dom::BlockDomain) =
    filter(el -> inside_worldtube(dom, el), dom.mesh.elements)

"""
    worldtube_boundary_faces(mesh::Mesh2D, isinside) -> Vector{Tuple{ElementId,Int,Int}}
    worldtube_boundary_faces(dom::BlockDomain) -> Vector{Tuple{ElementId,Int,Int}}

The worldtube boundary: EXACTLY the interior element faces whose owner and whose neighbour(s)
fall on opposite sides of the partition, as `(element, dim, sidx)` keys (`sidx` 1 ↔ side −1,
2 ↔ side +1) in element-layout order. Both the inside-owned and the outside-owned view of each
interface appears — each is a receiver of exterior data and needs its own jump datum.

`isinside(::ElementId) -> Bool` is the region predicate (the `BlockDomain` method uses
[`inside_worldtube`](@ref)). A face all of whose neighbours are on the owner's side is not a
worldtube face; a face that would straddle the partition PARTIALLY (some mortars crossing and
some not) is an `ArgumentError` — the regularized region is a union of whole blocks, so it
cannot happen, and if it ever did the per-face hook seam could not express it.

Domain-boundary faces are never included: at the production lattice the worldtube block's
equatorial faces are `:equator` DOMAIN boundaries (node `dg-boundary-conditions`' table), not
worldtube boundaries, while in the full-z chart the equator is an interior face INSIDE the tube
and carries no jump.
"""
function worldtube_boundary_faces(mesh::Mesh2D, isinside)
    out = Tuple{ElementId,Int,Int}[]
    for e in eachindex(mesh.elements), d in 1:2, sidx in 1:2
        iface = mesh.faces[e][d][sidx]
        iface === nothing && continue
        here = isinside(mesh.elements[e])::Bool
        n = count(k -> isinside(mesh.elements[iface.nbs[k]])::Bool != here,
                  1:length(iface.nbs))
        n == 0 && continue
        n == length(iface.nbs) || throw(ArgumentError(
            "worldtube_boundary_faces: face $((mesh.elements[e], d, sidx)) straddles the " *
            "regularized region only partially ($n of $(length(iface.nbs)) mortars cross) — " *
            "the region must be a union of whole blocks (main.tex:529)"))
        push!(out, (mesh.elements[e], d, sidx))
    end
    return out
end
worldtube_boundary_faces(dom::BlockDomain) =
    worldtube_boundary_faces(dom.mesh, el -> inside_worldtube(dom, el))

"""
    worldtube_face_points(mesh::Mesh2D, face) -> (xf, mortar_points)
    worldtube_face_points(dom::BlockDomain, face) -> (xf, mortar_points)

The physical evaluation points of one interior `face = (element, dim, sidx)`: `xf` is the
face's coordinate along its NORMAL dimension, and `mortar_points[k]` the transverse physical
coordinates of the face's `k`-th mortar, in the receiver's mortar order — exactly the points
at which EllipticCore's `MeshOperator` evaluates the exchanged mortar data, so a jump datum
built here lands on the hook's own grid.

A conforming face has one mortar whose points ARE the element's face collocation points
(bitwise: the mortar carries `max` of the two abutting orders over the shared transverse
segment, which is the element's own segment). A 2:1 h-nonconforming face seen from the COARSE
side has two mortars, each covering half the face at the FINE neighbour's resolution — which is
why the jump data are per-mortar (CRP main.tex:678-686).

Throws `ArgumentError` for an element outside `mesh` or for a domain-boundary face.

IMPLEMENTATION DEBT (for node `dg-jump-conditions`): the mortar geometry is reconstructed here
from `Mesh2D`'s precomputed face topology (`mesh.faces[e][d][sidx]`, an UNDOCUMENTED field —
EllipticCore exposes no public neighbour or mortar query) with the same arithmetic
`MeshOperator` uses to place its mortars. The right fix, alongside the per-mortar hook seam, is
to promote a public `mortar_points(mesh, face)` accessor into EllipticCore and have this
function call it; likewise the z² puncture-flux chain rule of [`worldtube_jump_values`](@ref)
belongs in a sealed helper next to [`flux_coeff_z2`](@ref) rather than in this file. Until
then the manufactured gates of test_worldtube.jl pin the agreement by running the assembled
data THROUGH `MeshOperator`, where a divergence in that arithmetic shows up at O(1).
"""
function worldtube_face_points(mesh::Mesh2D{RT}, face::Tuple{ElementId,Integer,Integer}) where {RT}
    el, dd, ss = face
    d, sidx = Int(dd), Int(ss)
    (d == 1 || d == 2) ||
        throw(ArgumentError("worldtube_face_points: dim must be 1 or 2, got $d"))
    (sidx == 1 || sidx == 2) ||
        throw(ArgumentError("worldtube_face_points: sidx must be 1 or 2, got $sidx"))
    haskey(mesh.index, el) ||
        throw(ArgumentError("worldtube_face_points: $el is not an element of this mesh"))
    e = mesh.index[el]
    iface = mesh.faces[e][d][sidx]
    iface === nothing && throw(ArgumentError(
        "worldtube_face_points: face $face is a DOMAIN boundary, not an interior interface"))
    blk = mesh.blocks[mesh.blockindex[el.block]]
    lo_d, hi_d = logical_interval(RT, el, d)
    xf = RT(map_forward(blockmap(blk, d), sidx == 1 ? lo_d : hi_d))
    tmap = blockmap(blk, 3 - d)
    pts = Vector{Vector{RT}}(undef, length(iface.nbs))
    for k in eachindex(pts)
        Nm = max(mesh.orders[e][3-d], mesh.orders[iface.nbs[k]][3-d])
        lo, hi = logical_interval(RT, iface.span[k])
        pts[k] = RT[map_forward(tmap, (lo * (1 - x) + hi * (1 + x)) / 2)
                    for x in refelement(RT, Nm).nodes]
    end
    return (xf, pts)
end
worldtube_face_points(dom::BlockDomain, face::Tuple{ElementId,Integer,Integer}) =
    worldtube_face_points(dom.mesh, face)

# ============================== jump-data assembly =========================================

_wt_eval(values, d::Int, xf, t) = d == 1 ? values(xf, t) : values(t, xf)

"""
    worldtube_jump_data(mesh::Mesh2D, isinside, values)
        -> Dict{Tuple{ElementId,Int,Int}, Vector{FaceState}}
    worldtube_jump_data(wt::WorldtubeData) -> Dict{…}

The worldtube jump data: for every face of [`worldtube_boundary_faces`](@ref), one
[`FaceState`](@ref) PER MORTAR (receiver mortar order) holding the additive shift
`(Δψ, Δφ)` of the RECEIVED exterior trace data — the eq:jump_worldtube realization of
main.tex:600-608, signed by [`worldtube_jump_signs`](@ref) and evaluated at
[`worldtube_face_points`](@ref).

`values(x1, x2)` supplies the puncture data at a physical point and must return, in its first
three positions, `(Ψ^𝒫, 𝓕^1(Ψ^𝒫), 𝓕^2(Ψ^𝒫))` — the value and the two flux components of the
chart. The face's normal dimension selects which flux component enters `Δφ`. For the CRP
problem this callable is [`worldtube_jump_values`](@ref) closed over the puncture, the mode
number and the chart, which is what [`worldtube_data`](@ref) does; a manufactured `values` is
what the node's DG-residual gate uses.

The data are SENDER-convention shifts (dg_operator.jl N3) and are TRIAL-FIELD INDEPENDENT, so
they belong in a [`FixedDataHook`](@ref) and reach the right-hand side through
`operator_on_zero` — never into the linearized matrix ([`worldtube_hooks`](@ref)).
"""
function worldtube_jump_data(mesh::Mesh2D{RT}, isinside, values) where {RT}
    faces = worldtube_boundary_faces(mesh, isinside)
    isempty(faces) &&
        return Dict{Tuple{ElementId,Int,Int},Vector{FaceState{Complex{RT}}}}()
    xf0, pts0 = worldtube_face_points(mesh, faces[1])
    v0 = _wt_eval(values, faces[1][2], xf0, pts0[1][1])
    T = promote_type(typeof(v0[1]), typeof(v0[2]), typeof(v0[3]))
    out = Dict{Tuple{ElementId,Int,Int},Vector{FaceState{T}}}()
    for face in faces
        el, d, sidx = face
        sg = worldtube_jump_signs(isinside(el)::Bool, sidx)
        xf, pts = worldtube_face_points(mesh, face)
        ms = Vector{FaceState{T}}(undef, length(pts))
        for (k, xt) in enumerate(pts)
            np = length(xt)
            dpsi = Vector{T}(undef, np)
            dflx = Vector{T}(undef, np)
            for i in 1:np
                v = _wt_eval(values, d, xf, xt[i])
                dpsi[i] = sg.psi * v[1]
                dflx[i] = sg.flux * (d == 1 ? v[2] : v[3])
            end
            ms[k] = FaceState{T}(dpsi, dflx)
        end
        out[face] = ms
    end
    return out
end

"""
    worldtube_hooks(jump_data::AbstractDict; acknowledge_omitted = nothing)
        -> (hooks::Dict{Tuple{ElementId,Int,Int},BoundaryDataHook}, omitted::Vector)
    worldtube_hooks(wt::WorldtubeData; acknowledge_omitted = nothing) -> (hooks, omitted)

The worldtube jump data packed as EllipticCore [`FixedDataHook`](@ref)s for the
`MeshOperator(…; hooks = …)` table: trial-field-independent additive shifts of the received
exterior data, identity in the LINEARIZED application, reaching the right-hand side through
`operator_on_zero` (`b = M·S − operator_on_zero(op)`; spectre
`impose_inhomogeneous_boundary_conditions_on_source`). Returns the table AND the list of faces
it could NOT express, so an omission can never pass unnoticed.

**DISCHARGED — use [`jump_hooks`](@ref) for any production table.** This packer binds ONE hook
object per FACE, which the EllipticCore seam then applies to every mortar of that face, while
spectre's `ModifyBoundaryData` is keyed by `mortar_id`. A conforming face packs exactly, and so
does an h-nonconforming face whose mortars happen to carry identical data; a 2:1 face whose
mortars carry DIFFERENT Ψ^𝒫 cannot be expressed at all. Node `dg-jump-conditions` closed that
gap on both sides: `MeshOperator` now accepts one hook PER MORTAR (EllipticCore matrixfree.jl,
[`mortar_points`](@ref)), and `Solve/jump_glue.jl`'s [`jump_hooks`](@ref) consumes
[`worldtube_jump_data`](@ref) directly — already per mortar, in the receiver's mortar order —
and composes it with the vtu junction hooks. `worldtube_hooks` stays EXPORTED as the per-FACE
packer for conforming-face and DIAGNOSTIC use (the figure drivers and the tests below); it is
NOT the production path.

The production lattice HAS such faces, and their number GROWS — which is why this packer is
diagnostics-only there:

| initial level `L` | unpackable worldtube faces (both charts) |
|---|---|
| 0 | 1 |
| 1 | 2 |
| L | `2^L` |

because the worldtube column is pre-refined to `(L+1, L+1)` while its left neighbour column
stays at `(L, L)` (yaml:78-93), so every one of that column's `2^L` transverse element faces
meets two finer ones. Node `amr-driver` then refines the puncture column eight more times, so
the count keeps growing. **`worldtube_hooks` is therefore not production-wirable on that
lattice; nodes `solver-stack-ssf` and `amr-driver` build their tables with
[`jump_hooks`](@ref).**

`acknowledge_omitted` is the only way to get a partial table, and it is DIAGNOSTICS-ONLY:

- `nothing` (the default) — any unpackable face is an `ArgumentError` naming the faces.
- a collection of face keys — must list EXACTLY the unpackable faces (as a set), otherwise an
  `ArgumentError` names the difference. Those faces are then omitted from the table.

An omitted face carries NO jump: the operator there is bitwise the no-hook operator, so the
reconstructed field is discontinuous across it and the error GROWS with resolution (measured
on the production lattice: the split face's owner sits at 3.4e-2 at `N = 4` and 8.9e-2 at
`N = 10`, versus 2.2e-3 -> 3.9e-4 when the same face IS served per mortar). **The solver driver
must never pass `acknowledge_omitted` — it must call [`jump_hooks`](@ref) instead.**
"""
function worldtube_hooks(jump_data::AbstractDict; acknowledge_omitted = nothing)
    unpackable = Tuple{ElementId,Int,Int}[]
    for (face, ms) in jump_data
        length(ms) > 1 &&
            !all(k -> ms[k].psi == ms[1].psi && ms[k].nflux == ms[1].nflux, 2:length(ms)) &&
            push!(unpackable, face)
    end
    sort!(unpackable; by = f -> (f[1].block, f[1].segments[1].level, f[1].segments[1].index,
                                 f[1].segments[2].level, f[1].segments[2].index, f[2], f[3]))
    if !isempty(unpackable)
        if acknowledge_omitted === nothing
            throw(ArgumentError(
                "worldtube_hooks: $(length(unpackable)) worldtube face(s) are " *
                "h-nonconforming with DIFFERENT Ψ^P jump data on their mortars, which the " *
                "sealed per-FACE FixedDataHook seam cannot express (EllipticCore " *
                "matrixfree.jl applies one hook object to every mortar of a face, while " *
                "spectre's ModifyBoundaryData is keyed by mortar_id): $unpackable. This is " *
                "the node dg-jump-conditions BLOCKER — extend MeshOperator to bind a hook " *
                "or fixed delta per mortar index / neighbour id and consume " *
                "worldtube_jump_data directly (it is already per-mortar). To obtain a " *
                "DIAGNOSTIC partial table, pass acknowledge_omitted = that exact list; the " *
                "solver driver must never do so, because an omitted face carries no jump " *
                "at all and its error grows with resolution"))
        end
        ack = Set(Tuple{ElementId,Int,Int}[(f[1], Int(f[2]), Int(f[3]))
                                           for f in acknowledge_omitted])
        ack == Set(unpackable) || throw(ArgumentError(
            "worldtube_hooks: acknowledge_omitted must list EXACTLY the faces that cannot " *
            "be packed. Missing from it: $(sort(collect(setdiff(Set(unpackable), ack)); by = string)); " *
            "listed but packable: $(sort(collect(setdiff(ack, Set(unpackable))); by = string))"))
    elseif acknowledge_omitted !== nothing && !isempty(acknowledge_omitted)
        throw(ArgumentError(
            "worldtube_hooks: acknowledge_omitted lists $(collect(acknowledge_omitted)) but " *
            "every worldtube face packs exactly — nothing is omitted"))
    end
    omit = Set(unpackable)
    out = Dict{Tuple{ElementId,Int,Int},BoundaryDataHook}()
    for (face, ms) in jump_data
        face in omit && continue
        out[face] = FixedDataHook(ms[1])
    end
    return (out, unpackable)
end

# ============================== fixed-source assembly ======================================

"""
    worldtube_source(mesh::Mesh2D, isinside, source, T) -> Vector{T}
    worldtube_source(wt::WorldtubeData) -> Vector

The DG fixed source as a global nodal vector in the mesh layout: `source(x1, x2)` evaluated at
every collocation node of every INSIDE element and EXACTLY zero everywhere else — the
`△_mΨ^ℛ_m = S^eff_m` / `△_mΨ_m = 0` split of eq:jump_worldtube, and spectre's
`InitializeEffectiveSource.hpp:150-176` (`fixed_sources` = the effective source inside, a
zeroed `Variables` outside). For the CRP problem [`worldtube_data`](@ref) supplies
`source = +`[`effective_source_m`](@ref) — the SEALED sign of the row-36 contract, the same as
spectre's `Tags::FixedSource<Tags::MMode>`.

`T` is the solver scalar type (`ComplexF64` on the production path). The DG residual is the
MASSIVE form `R(Ψ) = M·S`, so the right-hand side the driver builds is
[`worldtube_massive_source`](@ref), not this vector.
"""
worldtube_source(mesh::Mesh2D, isinside, source, ::Type{T}) where {T<:Number} =
    _wt_source_massive(mesh, isinside, source, T)[1]

# ONE pass building both the nodal source and its mass-weighted twin: the per-element
# matrices are built once (they are the only reason this needs geometry at all), which
# matters on the refined meshes node `amr-driver` re-assembles at every level.
function _wt_source_massive(mesh::Mesh2D, isinside, source, ::Type{T}) where {T<:Number}
    nodal = zeros(T, mesh.ndofs)
    massive = zeros(T, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        isinside(el)::Bool || continue
        em = element_matrices(T, mesh.blocks[mesh.blockindex[el.block]], el, mesh.orders[e])
        off = mesh.offsets[e]
        nx, ny = mesh.orders[e]
        for j in 1:ny, i in 1:nx
            nodal[off + tensor_lin(i, j, nx)] = source(em.x[1][i], em.x[2][j])
        end
        rng = off+1 : off+nx*ny
        massive[rng] .= em.M * view(nodal, rng)      # M is diagonal (Gauss collocation)
    end
    return (nodal, massive)
end

"""
    worldtube_massive_source(mesh::Mesh2D, isinside, source, T) -> Vector{T}
    worldtube_massive_source(wt::WorldtubeData) -> Vector

The fixed source in the MASSIVE form `M·S` the DG residual consumes: [`worldtube_source`](@ref)
multiplied element-by-element by the diagonal mass matrix (spectre's `apply_mass_matrix` on the
fixed sources when `Massive: True`, yaml:170). The inhomogeneous problem is then solved as

    A_lin·u = worldtube_massive_source(wt) − operator_on_zero(op) ,

with `op` carrying the worldtube [`worldtube_hooks`](@ref) (whose lift `operator_on_zero`
harvests) and `A_lin` the linearized [`apply!`](@ref) matrix.
"""
worldtube_massive_source(mesh::Mesh2D, isinside, source, ::Type{T}) where {T<:Number} =
    _wt_source_massive(mesh, isinside, source, T)[2]

# ============================== the assembled bundle =======================================

"""
    WorldtubeData{T<:Number,RT<:AbstractFloat,PT<:AbstractFloat}

The complete worldtube formulation of one m-mode problem on one [`BlockDomain`](@ref), as pure
data (build with [`worldtube_data`](@ref)). Fields:

- `dom::BlockDomain{RT}`, `orbit::CircularOrbit{PT}`, `puncture::Puncture{PT}`, `m::Int`,
  `phi_p::RT`, `chart::Symbol` — the problem this data belongs to (`chart` is `dom.chart`;
  `PT` is the orbit's own real type, normally `== RT`).
- `elements::Vector{ElementId}` — the inside elements ([`worldtube_elements`](@ref)).
- `jumps` — the per-face, per-mortar jump `FaceState`s ([`worldtube_jump_data`](@ref)).
- `source`, `massive_source` — the nodal and `M·S` fixed-source vectors
  ([`worldtube_source`](@ref), [`worldtube_massive_source`](@ref)).

`T` is the data scalar type (`Complex{RT}`). Nothing here is a solver object: node
`dg-jump-conditions` binds `jumps` into a `MeshOperator` alongside the vtu junction hooks, and
`solver-stack-ssf` composes the right-hand side.
"""
struct WorldtubeData{T<:Number,RT<:AbstractFloat,PT<:AbstractFloat}
    dom::BlockDomain{RT}
    orbit::CircularOrbit{PT}
    puncture::Puncture{PT}
    m::Int
    phi_p::RT
    chart::Symbol
    elements::Vector{ElementId}
    jumps::Dict{Tuple{ElementId,Int,Int},Vector{FaceState{T}}}
    source::Vector{T}
    massive_source::Vector{T}
end

"""
    worldtube_data(dom::BlockDomain, orbit::CircularOrbit, config::RunConfig,
                   m = config.orbit.m_mode_number; phi_p = 0) -> WorldtubeData

Assemble the whole worldtube formulation of eq:jump_worldtube on `dom` in ONE pass: the
regularized-region partition, the per-face/per-mortar Ψ^𝒫 jump data, and the fixed source
(nodal and massive). Every puncture value comes from ONE sealed
[`effective_source_bundle`](@ref) call per evaluation point through
[`worldtube_jump_values`](@ref).

`m` must equal `config.orbit.m_mode_number` (an `ArgumentError` otherwise): the DG operator's
own coefficients come from `domain_coefficients(dom, orbit, config)`, which reads the mode
number from the same config, and jump data for a different mode would silently belong to a
different equation. The mesh itself is m-independent (main.tex:543), so a different mode is a
different `crp_paper_config(m_mode = …)` over the same lattice. The domain's angular chart must
likewise match `config.orbit.impose_equatorial_symmetry`.

SLICING (inherited, not re-decided): the worldtube column lies strictly inside the vtu `t` slab
(main.tex:545-548), where the boost `H = 0` — which is exactly the slab the sealed effective
source is defined on (`puncture_psi_chart`'s binding t-slab contract). The jump data therefore
need no boost and carry none.

`mesh` defaults to `dom.mesh` but may be ANY mesh built from `dom`'s block lattice — a coarser
[`coarsen_mesh`](@ref) level or an AMR-refined one. A `WorldtubeData` is MESH-BOUND (element
ids, mortar counts and the nodal source length all belong to one mesh), so a per-level solver
rebuilds it for each level and pairs it with the matching hook table:

    wt_l    = worldtube_data(dom, orbit, config, m; mesh = mesh_l)
    hooks_l = jump_hooks(dom, wt_l, orbit, config, m; mesh = mesh_l)

The region predicate stays block-based, so `inside_worldtube`, `elements` and the sign table are
correct on every level; only the geometry is re-evaluated.

Cost is one puncture evaluation per inside collocation node plus one per worldtube mortar point
(20 elements × `P²` + ~12·`P` at the production lattice); the result is plain data and can be
cached across solver iterations for a fixed mesh.
"""
function worldtube_data(dom::BlockDomain{RT}, orbit::CircularOrbit{PT}, config::RunConfig,
                        m::Integer = config.orbit.m_mode_number;
                        phi_p::Real = 0, mesh::Mesh2D = dom.mesh) where {RT,PT}
    Int(m) == config.orbit.m_mode_number || throw(ArgumentError(
        "worldtube_data: m = $m disagrees with the config's MModeNumber " *
        "$(config.orbit.m_mode_number) — the operator coefficients are built from the same " *
        "config, so the jump data would belong to a different equation (the lattice is " *
        "m-independent, main.tex:543: use crp_paper_config(m_mode = $m))"))
    cchart = config.orbit.impose_equatorial_symmetry ? :z2 : :z
    cchart === dom.chart || throw(ArgumentError(
        "worldtube_data: the config's angular chart ($cchart, from " *
        "ImposeEquatorialSymmetry) does not match the domain chart $(dom.chart)"))
    mm = Int(m)
    p = Puncture(orbit)
    chart = dom.chart
    pp = RT(phi_p)
    isinside(el::ElementId) = inside_worldtube(dom, el)
    vals(x1, x2) = worldtube_jump_values(p, mm, chart, x1, x2; phi_p = pp)
    jumps = worldtube_jump_data(mesh, isinside, vals)
    T = isempty(jumps) ? Complex{RT} : eltype(first(values(jumps))[1].psi)
    src, msrc = _wt_source_massive(mesh, isinside, (x1, x2) -> vals(x1, x2).source, T)
    els = filter(isinside, mesh.elements)
    return WorldtubeData{T,RT,PT}(dom, orbit, p, mm, pp, chart, els, jumps, src, msrc)
end

inside_worldtube(wt::WorldtubeData, el::ElementId) = inside_worldtube(wt.dom, el)
worldtube_elements(wt::WorldtubeData) = wt.elements
worldtube_jump_data(wt::WorldtubeData) = wt.jumps
worldtube_hooks(wt::WorldtubeData; kwargs...) = worldtube_hooks(wt.jumps; kwargs...)
worldtube_source(wt::WorldtubeData) = wt.source
worldtube_massive_source(wt::WorldtubeData) = wt.massive_source

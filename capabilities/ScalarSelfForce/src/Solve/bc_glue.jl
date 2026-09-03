# Solve/bc_glue.jl — dg-boundary-conditions (ESSF node, schedule W8-D, GLOBAL_DAG row 41):
# the binding of the CRP's PHYSICAL boundary conditions to the EllipticCore ghost seam, i.e.
# the `MeshOperator` `bcs` table of the block domain of arXiv:2606.04998 (main.tex:396-401
# + 612-618).
#
# The paper imposes boundary conditions "through the DG numerical flux by a choice of Ψ_ext
# and (n_i F^i)_ext" (main.tex:614), and the CRP domain has exactly two kinds of external
# boundary:
#
#   * HORIZON / EQUATOR / POLES — "on boundaries where the perpendicular flux vanishes,
#     n_i 𝔄^i_m = 0, neither of the two boundary terms contributes, as Ψ_ext is ignored and
#     (n_i F^i)_ext = 0. This occurs on the horizon boundary and on angular boundaries"
#     (main.tex:614-615; the regularity discussion main.tex:396-401). NO boundary data are
#     imposed there — the numerical flux must reduce to the INTERIOR flux. That is the
#     [`NoneGhost`](@ref) of this file: spectre's `ScalarSelfForce::BoundaryConditions::None`
#     ("applies no boundary condition at all … nothing but regularity at the horizon"),
#     whose empty `apply` leaves the exterior at the sealed convention N7 transmissive state.
#   * THE OUTER RADIUS R — the second-order Bayliss–Turkel condition for m ≠ 0 (the sealed
#     `bayliss-turkel-bc` node's [`BaylissTurkelGhost`](@ref) via [`outer_boundary_ghost`](@ref)),
#     and `None` again for m = 0, which is the yaml's own Jinja2 branch (yaml:115-125).
#
# ZERO new physics: no coefficient, no Robin datum and no ghost formula is derived here. This
# file decides only WHICH sealed ghost each boundary face gets, and at WHICH points the
# sealed Robin coefficient is evaluated.
#
# WHY THE `:horizon` CASE IS SOUND: `block_domain` (sealed, row 39) REFUSES a lattice whose
# inner radial bound is not the outer horizon, so `A^rr = Δ/(r²+a²)` is exactly zero on every
# `:horizon` face (gated). Away from r₊ the same coefficient measures 5.9e-2 and is NEGATIVE
# inside the horizon — there "no boundary condition" would be a silent error, which is why
# the precondition is enforced upstream rather than assumed here.
#
# SCOPE BOUNDARIES (predecessors own the substance; later nodes own the rest of the wiring):
#   - MModeOperator/bayliss_turkel.jl (row 32) — `BaylissTurkelGhost`, `bt_robin_coefficient`,
#     `outer_boundary_ghost` and the m = 0 `nothing` branch;
#   - Domain/blocks.jl (row 39) — `block_domain` (+ the horizon precondition),
#     `boundary_face_classes`, `block_lattice_index`, `domain_coefficients`;
#   - EllipticCore boundary_hooks.jl (row 14) — the `GhostBoundary` seam and conventions
#     N1-N7; matrixfree.jl (row 15) — `boundary_faces` and the `bcs` contract;
#   - dg-jump-conditions (row 40) owns the worldtube/vtu `hooks` table (INTERIOR faces);
#     solver-stack-ssf (row 42) owns the driver that assembles operator + source + solver.

export NoneGhost, boundary_conditions

"""
    NoneGhost() <: EllipticCore.GhostBoundary

The "no boundary condition at all" ghost — the faithful Julia port of spectre's
`ScalarSelfForce::BoundaryConditions::None` (`.../Scalar/BoundaryConditions/None.hpp`:
"Applies no boundary condition at all. Used to impose nothing but regularity at the horizon
in horizon-penetrating coordinates or at angular boundaries"), which is what the CRP's input
file selects with `BoundaryConditions[r].Lower: None`, `BoundaryConditions[z²]: None` and —
for m = 0 — `BoundaryConditions[r].Upper: None` (yaml:113, 127, 115-125).

spectre's `None::apply` and `None::apply_linearized` are both EMPTY: the boundary values are
left at the interior trace, and the elliptic DG operator then forms the exterior state as
`u_ext = 2u_bc − u_int`, `(n·F)_ext = −2(n·F)_bc + (n·F)_int`
(`Elliptic/DiscontinuousGalerkin/DgOperator.hpp:475-508`). With `u_bc = u_int` and
`(n·F)_bc = (n·F)_int` that is exactly the SENDER-convention (N3) state

    ghost_state(NoneGhost(), int) = FaceState(ψ_int, −φ_int) ,

the *transmissive* exterior of the sealed `dg_operator.jl` convention N7. Through ANY
numerical flux of the [`PenaltyFlux`](@ref)/[`CentralFlux`](@ref) family it gives

    Ψ*    = (ψ_int + ψ_int)/2      = ψ_int          (bitwise; the penalty jump is zero),
    (nF)* = (φ_int − (−φ_int))/2   = φ_int          (bitwise),

so BOTH lifted face terms of `strong_residual!` — `MLDᵀ(s𝔄(Ψ* − ψ_int))` and
`−ML((nF)* − φ_int)` — are the exact zero vector: the face contributes nothing to the
residual and imposes nothing on the solution. The ghost is linear-homogeneous, so
`ghost_state_linearized` IS `ghost_state` (nothing to move to the right-hand side by
operator-on-zero) and the face stays inert inside every GMRES/Schwarz/multigrid application.

PHYSICS. Where the CRP uses it, that inertness is not an approximation but the paper's own
statement (main.tex:614-615): the perpendicular flux coefficient vanishes identically there
— `A^rr = Δ/(r²+a²) = 0` at r = r₊ and `A^{z²} = 4w(1−w)/(r²+a²) = 0` at the equator w = 0
and the poles w = 1 (`A^zz = (1−z²)/(r²+a²) = 0` at z = ±1 in the full chart) — so the
elliptic operator is degenerate on those faces and admits no boundary data at all: the
polynomial trial space selects the regular solution by itself (main.tex:396-401). At the
m = 0 OUTER face the flux coefficient does NOT vanish, and `NoneGhost` there realizes a
different, deliberate statement: no condition is imposed on the outgoing field at R because
the m = 0 mode is static (k = mΩ = 0, nothing radiates) and the Bayliss–Turkel data
−αβ₀/γ^r_0 is ill-defined (γ^r_0 vanishes identically at a = 0). Measured consequence: the
m = 0 solution at the extraction radius changes by 1.1e-10 relative when R = 1e4 is replaced
by 1e5 (test gate (d)).

Distinct from the sealed [`DirichletGhost`](@ref)/[`NeumannGhost`](@ref), which cannot
express this state with trial-independent data: `DirichletGhost(g)` would need `g = ψ_int`
and `NeumannGhost(h)` would need `h = φ_int`. At a face where 𝔄 vanishes a
`DirichletGhost` of ANY data happens to produce the same RESIDUAL (its data ride the
vanishing 𝔄 — "Ψ_ext is ignored"), but `NeumannGhost(0)` does not (measured 7.7e-1 …
5.5e0 relative deviation), and at the m = 0 outer face `DirichletGhost(0)` is an O(1)
physical error (it drives the contamination gate from 8.3e-8 to 1.7e-6).
"""
struct NoneGhost <: GhostBoundary end

EllipticCore.ghost_state(::NoneGhost, int::FaceState) = FaceState(int.psi, -int.nflux)
EllipticCore.ghost_state_linearized(bc::NoneGhost, int::FaceState) = ghost_state(bc, int)

# `mesh` must be a mesh of THIS domain's lattice, not merely of one with the same shape.
# Block ids alone do not establish that — two `BlockDomain`s built from different configs
# carry the same 5×2 ids — so the block MAPS are checked against the domain's own bounds.
# Endpoints are exact by the mesh2d contract (`block_domain` builds the maps from these very
# numbers), so the comparison is exact equality, and it costs one pass over the ~10 blocks.
# Without it, the `mkop` closure of a solver driver could silently pick up the wrong domain.
function _bc_require_same_lattice(dom::BlockDomain, mesh::Mesh2D)
    for blk in mesh.blocks
        i, j = block_lattice_index(dom, blk.id)    # throws for a block outside the lattice
        for (d, bounds, k) in ((1, dom.bounds_r, i), (2, dom.bounds_z, j))
            mp = blockmap(blk, d)
            lo = map_forward(mp, -one(eltype(bounds)))
            hi = map_forward(mp, one(eltype(bounds)))
            (lo == bounds[k+1] && hi == bounds[k+2]) || throw(ArgumentError(
                "boundary_conditions: block $(blk.id) of the given mesh spans [$lo, $hi] in " *
                "dimension $d, but lattice cell ($i, $j) of this BlockDomain spans " *
                "[$(bounds[k+1]), $(bounds[k+2])] — the mesh belongs to a DIFFERENT domain"))
        end
    end
    return nothing
end

# The class of ONE domain-boundary face of a mesh built from `dom`'s lattice. The lattice is
# a full rectangle with every interior block face connected, so a face without a neighbour
# lies on the perimeter, and the class then follows from the face's dimension and side
# through the SHARED [`boundary_class_rule`](@ref) — the one-line rule `block_domain` itself
# uses to build `boundary_face_classes(dom)`, CALLED here rather than duplicated (that the
# two agree on `dom.mesh` and on a coarsened mesh is gated). Stated as a FUNCTION rather than a lookup
# because the multigrid hierarchy hands `build_op` COARSER meshes whose element ids do not
# exist in `dom.boundary_classes`. A face that is not on the lattice perimeter — i.e. a mesh
# that is not derived from this domain — is refused rather than classified by default.
function _bc_face_class(dom::BlockDomain, face::Tuple{ElementId,Int,Int})
    el, d, sidx = face
    i, j = block_lattice_index(dom, el.block)      # throws for a block outside the lattice
    seg = el.segments[d]
    atedge = sidx == 1 ? seg.index == 0 : seg.index == (Int64(1) << seg.level) - 1
    perimeter = d == 1 ? (sidx == 1 ? i == 0 : i == dom.nr - 1) :
                         (sidx == 1 ? j == 0 : j == dom.nz - 1)
    (atedge && perimeter) || throw(ArgumentError(
        "boundary_conditions: face $face is not on the perimeter of the " *
        "$(dom.nr)×$(dom.nz) lattice (lattice cell ($i, $j), segment $(seg)) — the mesh is " *
        "not derived from this BlockDomain"))
    return boundary_class_rule(dom.chart, d, sidx)
end

# The outer-boundary ghost of ONE face: the sealed Robin data evaluated at exactly the points
# `element_coefficients` samples 𝔄 on that face (convention N2 — the face's normal-dimension
# block-map endpoint × the element's transverse volume nodes), so the ghost's coefficients
# and the operator's `Aface` are the SAME numbers, bitwise. `element_matrices` is rebuilt
# here rather than the formula duplicated: it is the object `element_coefficients` itself
# consumes, so the sampling cannot drift, and only the handful of outer faces pay for it.
function _bc_outer_ghost(dom::BlockDomain{RT}, mesh::Mesh2D, config::RunConfig, fcs,
                         face::Tuple{ElementId,Int,Int}) where {RT}
    el, d, sidx = face
    (d == 1 && sidx == 2) || throw(ArgumentError(
        "boundary_conditions: the outer radial boundary must be a (dim 1, side +1) face, " *
        "got $face"))
    e = mesh.index[el]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    em = element_matrices(RT, blk, el, mesh.orders[e])
    _, hi = logical_interval(RT, el, d)
    xf = RT(map_forward(blockmap(blk, d), hi))
    xt = em.x[3-d]
    ghost = outer_boundary_ghost(config,
                                 [fcs.A[1](xf, t) for t in xt],
                                 [fcs.beta(xf, t) for t in xt],
                                 [fcs.gamma[1](xf, t) for t in xt])
    return ghost === nothing ? NoneGhost() : ghost
end

"""
    boundary_conditions(dom::BlockDomain, orbit::CircularOrbit, config::RunConfig;
                        mesh = dom.mesh)
        -> Dict{Tuple{ElementId,Int,Int},GhostBoundary}

The `bcs` table of the CRP block domain: the physical boundary condition of every
domain-boundary face of `mesh`, in the exact shape [`MeshOperator`](@ref) demands (its
constructor requires a [`GhostBoundary`](@ref) for EVERY face of
[`boundary_faces`](@ref)`(mesh)` and rejects any other key). Everything is driven by
`config` — the m-mode number, the equatorial-symmetry chart and the Sommerfeld order all
enter through the sealed [`RunConfig`](@ref) alone:

- `:horizon` (r = r₊), `:equator` (z² = 0) and `:pole` (z² = 1, or z = ±1) →
  [`NoneGhost`](@ref). The perpendicular flux coefficient vanishes identically on those
  faces (main.tex:614-615; `block_domain` enforces the horizon placement that makes it
  true), so no boundary data exist to impose and the numerical flux reduces to the interior
  flux, bitwise.
- `:outer` (r = R) → the sealed [`outer_boundary_ghost`](@ref)`(config, A^rr, β_m, γ^r_m)`
  evaluated at that face's own quadrature points: a [`BaylissTurkelGhost`](@ref) carrying the
  pointwise [`bt_robin_coefficient`](@ref) for m ≠ 0 (main.tex:616-618), and
  [`NoneGhost`](@ref) for m = 0, where the yaml selects no outer condition at all
  (yaml:115-125) — see the `NoneGhost` docstring for what that realizes physically.

`mesh` defaults to `dom.mesh` but may be ANY mesh built from `dom`'s block lattice — in
particular the coarser meshes a [`MultigridHierarchy`](@ref) hands its `build_op` factory,
which must rebuild the boundary conditions for that level's own `boundary_faces`
(boundary-condition tables are keyed by element id and cannot be reused across levels). The
face class is therefore recomputed from the lattice geometry rather than looked up in
[`boundary_face_classes`](@ref)`(dom)` (whose keys are `dom.mesh`'s element ids); on
`dom.mesh` the two agree exactly, and a mesh whose faces are not on this lattice's perimeter
raises `ArgumentError` instead of being classified by default. The idiomatic solver wiring is

    dom  = block_domain(orbit, config)
    fcs  = domain_coefficients(dom, orbit, config)
    mkop = msh -> MeshOperator(ComplexF64, msh, fcs.A, fcs.beta, fcs.gamma;
                               C = config.dg.penalty_parameter,
                               bcs = boundary_conditions(dom, orbit, config; mesh = msh))
    op   = mkop(dom.mesh)
    gmres(op, b; Minv! = MultigridV(MultigridHierarchy(op, mkop)))

The closure argument is named `msh`, not `m`: `m` is the mode number everywhere else in this
package. A mesh belonging to a DIFFERENT `BlockDomain` — even one with the same lattice
shape — is refused, because the block maps are checked against `dom.bounds_r`/`dom.bounds_z`;
that is exactly the mistake such a closure could otherwise make silently.

WHAT GATES WHAT (read this before citing a number from this node). The reduced-ODE
contamination gate does NOT discriminate the Bayliss-Turkel wiring: replacing the outer
[`BaylissTurkelGhost`](@ref) by [`NoneGhost`](@ref) passes the same 1e-6 gate with a 12.6x
margin, and over 18 measured cells (two independent sources x m in {1, 2, 5} x P in
{4, 6, 8}) the ratio BT/None lies in [0.86, 4.20] — on this u-sliced problem the outgoing
field is non-oscillatory at R, so an absorbing outer condition buys nothing measurable
there. The BT wiring is justified by the paper (main.tex:403-411 + 614-618) and gated by the
sealed `bayliss-turkel-bc` node (the spectre Sommerfeld hard regressions, the O(1/R^2)
reflection scale) together with this node's bitwise checks — the ghost's `robin` IS the
sealed [`bt_robin_coefficient`](@ref) broadcast and the alpha it consumes IS the operator's
own `Aface`. A gate that discriminated the two conditions would need an oscillatory
(t-sliced or ingoing) probe.

ONE MEASURED CAVEAT FOR THE DRIVER, at m = 0: with no outer condition anywhere the operator
is about 4x worse conditioned than at m = 2 (cond 1.506e6 vs 3.824e5 at P = 6; 2.705e5 vs
6.471e4 at P = 4), and the sealed [`gmres`](@ref) can then report `converged = false,
breakdown = :stagnation` at an absolute residual of 3.8e-13 while the iterate is in fact
correct to 4.5e-13 relative to `direct_solve`. An m = 0 solve must not be gated on
`res.converged` alone.

Type-generic: the geometry type of `dom` selects the coordinate precision and the sealed
coefficient callables carry it into the Robin data (a `BlockDomain{BigFloat}` yields
`BaylissTurkelGhost{Complex{BigFloat}}`). The returned dictionary is mutable and its value
type is the abstract `GhostBoundary`, so a caller may override individual faces (the
figure/experiment drivers do).
"""
function boundary_conditions(dom::BlockDomain{RT}, orbit::CircularOrbit,
                             config::RunConfig; mesh::Mesh2D = dom.mesh) where {RT}
    mesh isa Mesh2D{RT} || throw(ArgumentError(
        "boundary_conditions: the mesh geometry type $(typeof(mesh)) does not match the " *
        "domain's Mesh2D{$RT}"))
    _bc_require_same_lattice(dom, mesh)
    fcs = domain_coefficients(dom, orbit, config)      # validates the chart against `config`
    bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}()
    for face in boundary_faces(mesh)
        bcs[face] = _bc_face_class(dom, face) === :outer ?
                    _bc_outer_ghost(dom, mesh, config, fcs, face) : NoneGhost()
    end
    return bcs
end

# test_mortar_hooks.jl — gates for the PER-MORTAR widening of the ModifyBoundaryData seam
# and for the public mortar-geometry accessor (DAG node: dg-jump-conditions, elliptic-ssf
# GLOBAL_DAG row 40, commit 1 — the EllipticCore half). Written BEFORE the implementation
# (TDD red run captured in evidence/dg-jump-conditions/tests_red_ec.log).
#
# WHY. `MeshOperator` bound ONE `BoundaryDataHook` per FACE key `(ElementId, d, sidx)` and
# applied it to every mortar of that face, while spectre's `ModifyBoundaryData` is keyed by
# `mortar_id` (Elliptic/Systems/SelfForce/Scalar/Equations.cpp:70-91 + DgOperator.hpp:708-728).
# The CRP production lattice pre-refines the worldtube column one level above its left
# neighbour (yaml:83-88), so its 2:1 worldtube faces carry genuinely DIFFERENT Ψ^P on their two
# mortars — the row-38 `worldtube-formulation` node landed that as a [HOLE] ("EllipticCore binds
# one hook per FACE … row 40 must widen the seam"). This file gates the widening:
#
#   * a face key may now carry EITHER a single hook (unchanged semantics, BITWISE-identical
#     results — gate (a1)) OR a collection of hooks, one per mortar in the RECEIVER's mortar
#     order (`mesh.faces[e][d][sidx].nbs` order, `MORTAR_LOWER` before `MORTAR_UPPER`);
#   * a length mismatch is an `ArgumentError` naming the face (seam-validation gate);
#   * `mortar_points(mesh, face)` publishes the mortar quadrature geometry the operator itself
#     places, so a caller building per-mortar data lands on the operator's own points
#     (gate (a5)) instead of re-deriving it from the previously undocumented `Mesh2D.faces`.
#
# ORACLE GATES (constants calibrated PRE-RED in evidence/dg-jump-conditions/calib_ec1.jl on a
# throwaway prototype of the mortar placement — NOT the final API):
#   (a1) per-mortar ≡ per-face when every mortar carries the same hook object: `apply!`,
#        `apply_full!`, `operator_on_zero` and `build_matrix` agree BITWISE (`==`), on a
#        conforming face (1 mortar) and on a 2:1 face (2 mortars).
#   (a2) NODE GATE — a 2:1 face whose prescribed jump depends on the TRANSVERSE coordinate:
#        the manufactured piecewise solution `[Ψ] = J(y)`, `[∂Ψ] = K(y)` is reproduced to
#        ≤ 1e-10 (calibrated on the CONFORMING twin, which the sealed per-face seam already
#        expresses exactly: solve 6.0e-15, residual 3.1e-13). FALSIFIERS: the same table with
#        the two mortars' data SWAPPED, and with the hooks removed, miss by O(1e2) in the
#        residual (calibrated: uniform-mortar-1 1.07e+02, uniform-mortar-2 1.20e+02, no-hook
#        2.41e+02 relative) — the per-mortar granularity is load-bearing.
#   (a3) fixed data stay OUT of the matrix and are harvested exactly by `operator_on_zero`:
#        the linearized `apply!` and `build_matrix` are BITWISE the unhooked ones; the lift is
#        additive over the mortars (superposition, ≤ 1e-13 rel) and the affine decomposition
#        `A_full(u) = A_lin(u) + A_full(0)` holds (≤ 1e-13 rel).
#   (a4) per-mortar LINEARIZED hooks DO enter the matrix: `build_matrix` changes when one
#        mortar carries a `LinearizedHook`, and differs again when both do — the row-40 oracle
#        (the vtu junction term) in miniature.
#   (a5) `mortar_points(mesh, face)` equals the placement `MeshOperator` uses, BITWISE: the
#        operator evaluates 𝔄 ONCE per interface mortar at the mortar's physical points and
#        stores it as `Mortar.Am`, so an injective coefficient probe reads the placement back
#        out of the built operator. Both dimensions, both sides of every interface.
#
# Helpers are deliberately UNtyped where they touch the new API so this file LOADS before the
# implementation exists (each @testset then errors individually in the red run).

using Test
using LinearAlgebra
using EllipticCore

const PMH_N = 6                       # piecewise-cubic exact at N ≥ 4
const PMH_XM = 0.5                    # the jump interface abscissa
const PMH_BID = 77
const PMH_BETA = 0.3 + 0.2im

pmh_blk() = Block(PMH_BID, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
pmh_one(x, y) = 1.0
pmh_zero(x, y) = 0.0
pmh_beta(x, y) = PMH_BETA

# ---- manufactured piecewise solution with a TRANSVERSE-DEPENDENT jump ---------------------
# u_L is a complex cubic in x plus a quadratic in y; across x = x_m the field jumps by J(y)
# and its x-derivative by K(y), so the two halves of a 2:1 face carry DIFFERENT data — the
# production worldtube situation in miniature. (−∇² + β)u = S with S built from the exact
# Laplacian, so the discrete problem is exact for N ≥ 4 (Gauss collocation, 𝔄 ≡ 1).
pmh_uL(x, y) = (0.8 + 0.5im) * x^3 - (0.4 - 0.2im) * x + 0.3 + (0.25 + 0.1im) * y^2
pmh_J(y) = (0.35 - 0.15im) + (0.22 + 0.31im) * y + (-0.18 + 0.12im) * y^2
pmh_K(y) = (0.6 + 0.25im) + (0.41 - 0.13im) * y
pmh_uR(x, y) = pmh_uL(x, y) + pmh_J(y) + pmh_K(y) * (x - PMH_XM)
pmh_uex(x, y) = x < PMH_XM ? pmh_uL(x, y) : pmh_uR(x, y)
pmh_lapL(x, y) = 6 * (0.8 + 0.5im) * x + 2 * (0.25 + 0.1im)
pmh_lapR(x, y) = pmh_lapL(x, y) + 2 * (-0.18 + 0.12im)
pmh_S(x, y) = (x < PMH_XM ? -pmh_lapL(x, y) : -pmh_lapR(x, y)) + PMH_BETA * pmh_uex(x, y)

# ---- meshes -------------------------------------------------------------------------------
# 2:1 model: ONE coarse left element (full transverse extent) facing TWO fine right elements,
# so the left element's +x face carries two mortars covering different halves of y.
function pmh_mesh_split(N = PMH_N)
    els = [ElementId(PMH_BID, SegmentId(1, 0), SegmentId(0, 0)),
           ElementId(PMH_BID, SegmentId(1, 1), SegmentId(1, 0)),
           ElementId(PMH_BID, SegmentId(1, 1), SegmentId(1, 1))]
    return Mesh2D(pmh_blk(), els, (N, N))
end
# conforming twin: 2 × 2 elements — every interface has exactly ONE mortar, so the sealed
# per-FACE seam already expresses the transverse-dependent jump (the (a2) calibration).
function pmh_mesh_conf(N = PMH_N)
    els = [ElementId(PMH_BID, SegmentId(1, ix), SegmentId(1, iy)) for iy in 0:1 for ix in 0:1]
    return Mesh2D(pmh_blk(), els, (N, N))
end

# ---- element-level helpers ----------------------------------------------------------------
pmh_em(mesh, e) = element_matrices(Float64, mesh.blocks[mesh.blockindex[mesh.elements[e].block]],
                                   mesh.elements[e], mesh.orders[e])

function pmh_facedata(mesh, el, d, sidx, f)
    e = mesh.index[el]
    b = mesh.blocks[mesh.blockindex[el.block]]
    em = pmh_em(mesh, e)
    lo, hi = logical_interval(Float64, el, d)
    xf = map_forward(blockmap(b, d), sidx == 1 ? lo : hi)
    xt = em.x[3-d]
    return d == 1 ? ComplexF64[f(xf, t) for t in xt] : ComplexF64[f(t, xf) for t in xt]
end

pmh_bcs(mesh, f = pmh_uex) = Dict{Tuple{ElementId,Int,Int},GhostBoundary}(
    k => DirichletGhost(pmh_facedata(mesh, k[1], k[2], k[3], f)) for k in boundary_faces(mesh))
pmh_zerobcs(mesh) = Dict{Tuple{ElementId,Int,Int},GhostBoundary}(
    k => DirichletGhost(zeros(ComplexF64, mesh.orders[mesh.index[k[1]]][3-k[2]]))
    for k in boundary_faces(mesh))

function pmh_nodal(mesh, f)
    v = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        em = pmh_em(mesh, e)
        nx, ny = mesh.orders[e]
        v[mesh.offsets[e]+1:mesh.offsets[e]+nx*ny] .=
            vec([f(em.x[1][i], em.x[2][j]) for i in 1:nx, j in 1:ny])
    end
    return v
end

function pmh_rhs(mesh)
    b = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        em = pmh_em(mesh, e)
        nx, ny = mesh.orders[e]
        rng = mesh.offsets[e]+1 : mesh.offsets[e]+nx*ny
        b[rng] .= em.M * vec([pmh_S(em.x[1][i], em.x[2][j]) for i in 1:nx, j in 1:ny])
    end
    return b
end

# the interior faces that lie ON the jump interface x = x_m, as (face, receiver-is-left)
function pmh_jump_faces(mesh)
    out = Tuple{Tuple{ElementId,Int,Int},Bool}[]
    for e in eachindex(mesh.elements), sidx in 1:2
        mesh.faces[e][1][sidx] === nothing && continue
        el = mesh.elements[e]
        b = mesh.blocks[mesh.blockindex[el.block]]
        lo, hi = logical_interval(Float64, el, 1)
        xf = map_forward(blockmap(b, 1), sidx == 1 ? lo : hi)
        xf == PMH_XM || continue
        push!(out, ((el, 1, sidx), sidx == 2))
    end
    return out
end

# the exact per-mortar jump datum: (Δψ, Δφ) = (−J, +K) on the LEFT receiver, (+J, +K) on the
# RIGHT one (the sealed prescribed-jump sign table, boundary_hooks.jl)
pmh_datum(isleft, ys) = FixedDataHook(
    ComplexF64[(isleft ? -1 : 1) * pmh_J(y) for y in ys], ComplexF64[pmh_K(y) for y in ys])

# hook tables built THROUGH the new public accessor (undefined in the red run)
function pmh_hooks(mesh; reverse_mortars = false)
    h = Dict{Tuple{ElementId,Int,Int},Any}()
    for (face, isleft) in pmh_jump_faces(mesh)
        _, pts = mortar_points(mesh, face)
        ms = [pmh_datum(isleft, p) for p in pts]
        reverse_mortars && (ms = reverse(ms))
        h[face] = ms
    end
    return h
end

function pmh_op(mesh, hooks; bcs = pmh_bcs(mesh))
    return MeshOperator(ComplexF64, mesh, (pmh_one, pmh_one), pmh_beta, (pmh_zero, pmh_zero);
                        C = 1.5, hooks = hooks, bcs = bcs)
end

# relative residual of the EXACT field through the FULL (affine) operator
function pmh_residual(mesh, op)
    ue = pmh_nodal(mesh, pmh_uex)
    b = pmh_rhs(mesh)
    y = similar(ue)
    apply_full!(y, op, ue)
    return maximum(abs, y .- b) / maximum(abs, b)
end

function pmh_solve_err(mesh, op)
    ue = pmh_nodal(mesh, pmh_uex)
    u = direct_solve(op, pmh_rhs(mesh) .- operator_on_zero(op))
    return maximum(abs, u .- ue) / maximum(abs, ue)
end

@testset "per-mortar boundary-data hooks + mortar_points (dg-jump-conditions)" begin

@testset "(a5) mortar_points IS the placement MeshOperator uses (bitwise)" begin
    @test Base.isexported(EllipticCore, :mortar_points)
    mesh = pmh_mesh_split()
    # injective coefficient probes: MeshOperator evaluates 𝔄 ONCE per interface mortar at that
    # mortar's physical points and stores it as `Mortar.Am` (matrixfree.jl W5-A contract), so
    # 𝔄 = the transverse coordinate reads the placement back, and 𝔄 = the normal coordinate
    # reads the face abscissa.
    probe_t = MeshOperator(ComplexF64, mesh, ((x, y) -> y, (x, y) -> x), pmh_zero,
                           (pmh_zero, pmh_zero); C = 1.5, bcs = pmh_zerobcs(mesh))
    probe_n = MeshOperator(ComplexF64, mesh, ((x, y) -> x, (x, y) -> y), pmh_zero,
                           (pmh_zero, pmh_zero); C = 1.5, bcs = pmh_zerobcs(mesh))
    nifaces = 0
    for e in eachindex(mesh.elements), d in 1:2, sidx in 1:2
        iface = mesh.faces[e][d][sidx]
        iface === nothing && continue
        nifaces += 1
        face = (mesh.elements[e], d, sidx)
        xf, pts = mortar_points(mesh, face)
        @test length(pts) == length(iface.nbs)
        mt = probe_t.flux[e].faces[d][sidx].mortars
        mn = probe_n.flux[e].faces[d][sidx].mortars
        for k in eachindex(pts)
            @test pts[k] == mt[k].Am           # bitwise: same arithmetic, same points
            @test all(==(xf), mn[k].Am)
        end
    end
    # 5 interior faces: the coarse +x face (2 mortars), the two fine −x faces, and the two
    # sides of the y-interface between the fine pair
    @test nifaces == 5

    # the coarse side of the 2:1 face carries TWO mortars, MORTAR_LOWER before MORTAR_UPPER;
    # the two fine sides carry one each, and the halves TILE the coarse face
    coarse = (mesh.elements[1], 1, 2)
    xf, pts = mortar_points(mesh, coarse)
    @test xf == PMH_XM
    @test length(pts) == 2
    @test mesh.faces[1][1][2].msize == [MORTAR_LOWER, MORTAR_UPPER]
    @test all(<(0.5), pts[1]) && all(>(0.5), pts[2])
    # both sides of an interface see the SAME physical points, bitwise (mortar (σ, Am) parity)
    @test mortar_points(mesh, (mesh.elements[2], 1, 1))[2][1] == pts[1]
    @test mortar_points(mesh, (mesh.elements[3], 1, 1))[2][1] == pts[2]
    @test mortar_points(mesh, (mesh.elements[2], 1, 1))[1] == xf

    # a CONFORMING face's single mortar is the element's own face collocation points
    cmesh = pmh_mesh_conf()
    ce = cmesh.index[ElementId(PMH_BID, SegmentId(1, 0), SegmentId(1, 0))]
    cxf, cpts = mortar_points(cmesh, (cmesh.elements[ce], 1, 2))
    @test cpts[1] == pmh_em(cmesh, ce).x[2]
    @test cxf == PMH_XM

    # errors: domain-boundary face, unknown element, bad dim/side
    @test_throws ArgumentError mortar_points(mesh, (mesh.elements[1], 1, 1))
    @test_throws ArgumentError mortar_points(mesh, (ElementId(999, SegmentId(0, 0),
                                                             SegmentId(0, 0)), 1, 2))
    @test_throws ArgumentError mortar_points(mesh, (mesh.elements[1], 3, 2))
    @test_throws ArgumentError mortar_points(mesh, (mesh.elements[1], 1, 3))

    # Mesh2D.faces / MeshInterface are PUBLIC-READ (the row-38 node had to reach into them)
    @test isdefined(EllipticCore, :MeshInterface)
    @test mesh.faces[1][1][2] isa EllipticCore.MeshInterface
    @test fieldnames(EllipticCore.MeshInterface) == (:nbs, :msize, :nbmsize, :span, :nbslot)
end

@testset "(a1) per-mortar ≡ per-face when the mortars carry ONE hook object (bitwise)" begin
    rng_delta(n, s) = FixedDataHook(ComplexF64[cis(s * k) / (k + 1) for k in 1:n],
                                   ComplexF64[cis(-s * k) * 0.7 for k in 1:n])
    fω(i, e) = FaceState(e.psi, e.nflux .- im * 0.9 .* (i.psi .+ e.psi) ./ 2)

    for (mesh, face) in ((pmh_mesh_conf(),
                          (ElementId(PMH_BID, SegmentId(1, 0), SegmentId(1, 0)), 1, 2)),
                         (pmh_mesh_split(), (ElementId(PMH_BID, SegmentId(1, 0),
                                                       SegmentId(0, 0)), 1, 2)))
        _, pts = mortar_points(mesh, face)
        nm = length(pts)
        for h in (rng_delta(length(pts[1]), 0.3), LinearizedHook(fω),
                  ComposedHook(LinearizedHook(fω), rng_delta(length(pts[1]), 1.1)))
            opf = pmh_op(mesh, Dict(face => h))
            opm = pmh_op(mesh, Dict(face => [h for _ in 1:nm]))
            u = pmh_nodal(mesh, (x, y) -> cis(3x) * (1 + y) + im * x * y)
            yf = similar(u); yl = similar(u)
            apply!(yf, opf, u); apply!(yl, opm, u)
            @test yf == yl                                     # bitwise
            apply_full!(yf, opf, u); apply_full!(yl, opm, u)
            @test yf == yl
            @test operator_on_zero(opf) == operator_on_zero(opm)
            @test build_matrix(opf) == build_matrix(opm)
        end
        @test nm == (length(mesh.elements) == 3 ? 2 : 1)
    end
end

@testset "(a2) NODE GATE: 2:1 face, transverse-dependent prescribed jump [Ψ]=J(y), [∂Ψ]=K(y)" begin
    # calibration reference: the CONFORMING twin, which the sealed per-FACE seam expresses
    # exactly (measured solve 6.0e-15, residual 3.1e-13)
    cmesh = pmh_mesh_conf()
    copp = pmh_op(cmesh, pmh_hooks(cmesh))
    @test pmh_solve_err(cmesh, copp) ≤ 1e-10
    @test pmh_residual(cmesh, copp) ≤ 1e-10

    # THE GATE: the same manufactured problem on the 2:1 mesh, expressible ONLY per mortar
    mesh = pmh_mesh_split()
    op = pmh_op(mesh, pmh_hooks(mesh))
    err = pmh_solve_err(mesh, op)
    res = pmh_residual(mesh, op)
    @test err ≤ 1e-10
    @test res ≤ 1e-10

    # FALSIFIER 1 — swap the two mortars' data on the coarse face (calibrated O(1e2))
    opr = pmh_op(mesh, pmh_hooks(mesh; reverse_mortars = true))
    resr = pmh_residual(mesh, opr)
    @test resr > 1e0
    @test resr > 1e6 * res
    # FALSIFIER 2 — no hooks at all (calibrated 2.41e+02)
    opn = pmh_op(mesh, Dict{Tuple{ElementId,Int,Int},Any}())
    @test pmh_residual(mesh, opn) > 1e0

    println("mortar-hooks (a2) 2:1 transverse jump: solve ", err, " residual ", res,
            " | falsifiers reversed ", resr, " no-hook ", pmh_residual(mesh, opn),
            " | conforming twin ", pmh_solve_err(cmesh, copp))
end

@testset "(a3) per-mortar FIXED data: out of the matrix, harvested by operator_on_zero" begin
    mesh = pmh_mesh_split()
    face = (mesh.elements[1], 1, 2)
    _, pts = mortar_points(mesh, face)
    n1, n2 = length(pts[1]), length(pts[2])
    h1 = FixedDataHook(ComplexF64[0.3 + 0.1im * k for k in 1:n1],
                       ComplexF64[-0.2im * k for k in 1:n1])
    h2 = FixedDataHook(ComplexF64[1.7 - 0.4im * k for k in 1:n2],
                       ComplexF64[0.9 + 0.05im * k for k in 1:n2])
    op0 = pmh_op(mesh, Dict{Tuple{ElementId,Int,Int},Any}(); bcs = pmh_zerobcs(mesh))
    opb = pmh_op(mesh, Dict(face => [h1, h2]); bcs = pmh_zerobcs(mesh))
    op1 = pmh_op(mesh, Dict(face => [h1, NoHook()]); bcs = pmh_zerobcs(mesh))
    op2 = pmh_op(mesh, Dict(face => [NoHook(), h2]); bcs = pmh_zerobcs(mesh))

    u = pmh_nodal(mesh, (x, y) -> cis(2x - y) + im * x^2)
    y0 = similar(u); yb = similar(u)
    apply!(y0, op0, u); apply!(yb, opb, u)
    @test y0 == yb                                     # fixed data are NOT in the matrix
    @test build_matrix(op0) == build_matrix(opb)

    # the lift is additive over the mortars (each mortar's datum enters its own lift only)
    zb = operator_on_zero(opb)
    z1 = operator_on_zero(op1)
    z2 = operator_on_zero(op2)
    @test maximum(abs, zb .- z1 .- z2) ≤ 1e-13 * maximum(abs, zb)
    @test maximum(abs, zb) > 1e-3                       # non-vacuous
    @test maximum(abs, z1 .- z2) > 1e-3                 # the two mortars differ genuinely

    # affine decomposition through the per-mortar path
    yfull = similar(u)
    apply_full!(yfull, opb, u)
    @test maximum(abs, yfull .- yb .- zb) ≤ 1e-13 * maximum(abs, yfull)
end

@testset "(a4) per-mortar LINEARIZED hooks enter build_matrix (the row-40 oracle)" begin
    mesh = pmh_mesh_split()
    face = (mesh.elements[1], 1, 2)
    fω(i, e) = FaceState(e.psi, e.nflux .- im * 0.9 .* (i.psi .+ e.psi) ./ 2)
    lh = LinearizedHook(fω)
    A0 = build_matrix(pmh_op(mesh, Dict{Tuple{ElementId,Int,Int},Any}();
                             bcs = pmh_zerobcs(mesh)))
    A1 = build_matrix(pmh_op(mesh, Dict(face => [lh, NoHook()]); bcs = pmh_zerobcs(mesh)))
    A2 = build_matrix(pmh_op(mesh, Dict(face => [lh, lh]); bcs = pmh_zerobcs(mesh)))
    Af = build_matrix(pmh_op(mesh, Dict(face => lh); bcs = pmh_zerobcs(mesh)))
    @test A1 != A0
    @test A2 != A1
    @test A2 == Af                                       # per-face ≡ replicated per-mortar
    @test maximum(abs, A1 .- A0) > 1e-3
    @test maximum(abs, A2 .- A1) > 1e-3
    # a pure-linearized per-mortar hook set keeps the operator linear-HOMOGENEOUS
    @test all(iszero, operator_on_zero(pmh_op(mesh, Dict(face => [lh, NoHook()]);
                                              bcs = pmh_zerobcs(mesh))))
end

@testset "seam validation: per-mortar length, element type, face keys" begin
    mesh = pmh_mesh_split()
    coarse = (mesh.elements[1], 1, 2)       # 2 mortars
    fine = (mesh.elements[2], 1, 1)         # 1 mortar
    h = FixedDataHook(zeros(ComplexF64, PMH_N), zeros(ComplexF64, PMH_N))
    bcs = pmh_zerobcs(mesh)
    # too few / too many entries: an ArgumentError naming the face
    for bad in ([h], [h, h, h], BoundaryDataHook[])
        e = @test_throws ArgumentError pmh_op(mesh, Dict(coarse => bad); bcs = bcs)
        @test occursin("mortar", e.value.msg)
        @test occursin(string(mesh.elements[1]), e.value.msg)
    end
    @test_throws ArgumentError pmh_op(mesh, Dict(fine => [h, h]); bcs = bcs)
    # a non-hook element is refused, not silently converted
    @test_throws ArgumentError pmh_op(mesh, Dict(coarse => [h, 1.0]); bcs = bcs)
    @test_throws ArgumentError pmh_op(mesh, Dict(coarse => 1.0); bcs = bcs)
    # exact lengths are accepted
    @test pmh_op(mesh, Dict(coarse => [h, h]); bcs = bcs) isa MeshOperator
    @test pmh_op(mesh, Dict(fine => [h]); bcs = bcs) isa MeshOperator
    # a hook on a DOMAIN-boundary face is refused, and the message points at `bcs`
    eb = @test_throws ArgumentError pmh_op(mesh, Dict((mesh.elements[1], 1, 1) => [h]);
                                           bcs = bcs)
    @test occursin("bcs", eb.value.msg) && occursin("DOMAIN-BOUNDARY", eb.value.msg)
    # a malformed key is an ArgumentError, not a BoundsError from indexing mesh.faces
    for bad in ((mesh.elements[1], 3, 2), (mesh.elements[1], 1, 0), (mesh.elements[1], 0, 1))
        e3 = @test_throws ArgumentError pmh_op(mesh, Dict(bad => h); bcs = bcs)
        @test occursin("dim", e3.value.msg)
    end
    e4 = @test_throws ArgumentError pmh_op(mesh,
             Dict((ElementId(999, SegmentId(0, 0), SegmentId(0, 0)), 1, 2) => h); bcs = bcs)
    @test occursin("not an element of this mesh", e4.value.msg)
    # a TUPLE of hooks is NOT the per-mortar spelling (documented): it is refused, not
    # silently taken as a single hook
    @test_throws ArgumentError pmh_op(mesh, Dict(coarse => (h, h)); bcs = bcs)
end

@testset "BigFloat inheritance through the per-mortar seam" begin
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        blk = Block(PMH_BID, LinearMap1D(big(0.0), big(1.0)), LinearMap1D(big(0.0), big(1.0)))
        els = [ElementId(PMH_BID, SegmentId(1, 0), SegmentId(0, 0)),
               ElementId(PMH_BID, SegmentId(1, 1), SegmentId(1, 0)),
               ElementId(PMH_BID, SegmentId(1, 1), SegmentId(1, 1))]
        mesh = Mesh2D(blk, els, (4, 4))
        face = (els[1], 1, 2)
        xf, pts = mortar_points(mesh, face)
        @test xf isa BigFloat && precision(xf) == 256
        @test eltype(pts[1]) == BigFloat
        hks = [FixedDataHook(Complex{BigFloat}[big(0.1) * k for k in 1:length(p)],
                             Complex{BigFloat}[big(0.2) * k for k in 1:length(p)])
               for p in pts]
        bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}(
            k => DirichletGhost(zeros(Complex{BigFloat}, 4)) for k in boundary_faces(mesh))
        op = MeshOperator(Complex{BigFloat}, mesh, ((x, y) -> one(x), (x, y) -> one(x)),
                          (x, y) -> Complex{BigFloat}(0.3, 0.2),
                          ((x, y) -> zero(x), (x, y) -> zero(x));
                          C = 1.5, hooks = Dict(face => hks), bcs = bcs)
        z = operator_on_zero(op)
        @test eltype(z) == Complex{BigFloat}
        @test precision(real(z[1])) == 256
        @test maximum(abs, z) > 1e-6
    finally
        setprecision(BigFloat, old)
    end
end

end # testset per-mortar boundary-data hooks

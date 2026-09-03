# test_worldtube.jl — ESSF node `worldtube-formulation` (schedule W8-C, GLOBAL_DAG row 38),
# written BEFORE the implementation (TDD). The node delivers `Worldtube/worldtube.jl`: the
# worldtube FORMULATION of arXiv:2606.04998 §II.D as PURE DATA for the DG layer —
#
#   △_m Ψ^ℛ_m = S^eff_m   (interior),   △_m Ψ_m = 0   (exterior),
#   Ψ^ℛ_m = Ψ_m − Ψ^𝒫_m   (across the worldtube interface)     [eq:jump_worldtube,
#                                                                main.tex:513-523]
#
# realized, per main.tex:596-608, by MODIFYING THE RECEIVED EXTERIOR DATA on every element
# that receives across the worldtube boundary,
#
#   Ψ^ext_m       →  Ψ^ext_m       ±  Ψ^𝒫_m
#   (n_i𝓕^i_m)ext →  (n_i𝓕^i_m)ext ∓  n_i 𝔄^{ij}_m ∂_j Ψ^𝒫_m
#
# plus the FIXED SOURCE +S^eff_m at every quadrature node inside the worldtube and zero
# outside (spectre `Tags::FixedSource<Tags::MMode>`; `InitializeEffectiveSource.hpp:150-176`).
#
# ORACLES (read-only mirrors, transcribed by hand below — never edited):
#   * `main.tex:505-523` (eq:jump_worldtube, the interior/exterior split), `main.tex:529`
#     (fig:domain — the regularized region is the red RECTANGLE, i.e. a whole block),
#     `main.tex:545-548` (the worldtube column sits strictly inside the vtu `t` slab),
#     `main.tex:596-608` (the ±Ψ^𝒫 / ∓n·𝓕(Ψ^𝒫) exterior-data modification and the remark
#     that "the change in sign on either side of the boundary is absorbed by the opposite
#     face normal").
#   * `ref-code/spectre` @ d3cf5af0:
#     - `src/Elliptic/Systems/SelfForce/Scalar/Equations.cpp:70-91` — `ModifyBoundaryData::
#       apply`: `sign = field_is_regularized ? -1 : +1`, then `field += sign·Ψ^𝒫` and
#       `n_dot_flux -= sign·(n·𝓕(Ψ^𝒫))`, with `n_dot_flux` the RECEIVED remote datum (sender
#       convention, EllipticCore N3) and the singular-field flux built on the RECEIVER's own
#       mortar normal (`Actions/InitializeEffectiveSource.hpp:205-215`). Both signs are
#       transcribed into the sign table gated in (2c).
#     - `Actions/InitializeEffectiveSource.hpp:130-176` — the regularized region is the whole
#       BLOCK containing the puncture (`puncture_in_element` → `block_logical_coordinates`),
#       the fixed source is the effective source inside and IDENTICALLY ZERO outside, and it
#       is mass-multiplied (`apply_mass_matrix`) for the massive DG formulation.
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here):
#   * `effective-source-direct` — `effective_source_bundle` (source, psi, dpsi_r, dpsi_z) and
#     `puncture_psi_chart`; the `phi_p` = Boyer–Lindquist azimuth convention, the internal
#     field-point rotation, the `t`-slab-only slicing contract, the ϱ² = 0 / |z| = 1
#     `DomainError` guards, and the `+effective_source_m` fixed-source sign.
#   * `first-order-flux-form` / `equatorial-symmetry-z2` — `flux_coeff_r`, `flux_coeff_z`,
#     `flux_coeff_z2` (𝔄 is DIAGONAL, main.tex:575) and the printed z² chain rule
#     𝓕^{z²} = A^{z²}∂_wΨ with ∂_w = ∂_z/(2z) (main.tex:389-394).
#   * `domain-blocks-map` — `BlockDomain`, `worldtube_blocks`, `boundary_class`,
#     `domain_coefficients`; `vtu-slicing-jumps` — the `t` region of the worldtube column.
#   * `boundary-data-hooks` / `matrix-free-operator` (EllipticCore) — `FixedDataHook`
#     (SENDER-convention additive shift, identity in the LINEARIZED application, reaching the
#     RHS through `operator_on_zero`), `MeshOperator`, `apply_full!`, the mortar exchange.
#
# THIS node's gates (all constants calibrated PRE-RED in the evidence dir,
# `calibrate_wt1.jl` / `calibrate_wt2.jl` / `calibrate_wt3.jl`):
#  (1) region partition — against an INDEPENDENT hand enumeration of the production lattice in
#      both the 5×2 (r, z²) and the 5×3 (r, z) chart: 11 faces / 12 mortars and 15 / 16, with
#      exactly one 2:1 SPLIT face (the coarse left neighbour column).
#  (2) jump data — bitwise equality of the ψ slot with the sealed `puncture_psi_chart`, the
#      flux slot against the sealed coefficient callables, the two-sided sender-normal
#      consistency (BITWISE, calibrated 0.0 over all 12 mortar pairs) and the EXACT interface
#      cancellation of the reconstructed field.
#  (3) manufactured split (the GLOBAL_DAG row-38 gate) — a two-region 2D mesh with affine
#      blocks and complex coefficients, smooth Ψ^𝒫_man/Ψ_man, the SHIPPED assembly (face
#      enumeration, mortar geometry, sign table, `FixedDataHook` packing, massive source):
#      the DG residual of the piecewise field vanishes to ≤ 1e-10 of the source scale
#      (calibrated 1.3e-12 conforming, 1.9e-13 on a 2:1 h-nonconforming worldtube boundary)
#      while each of the three sign flips breaks it by ≥ 0.29 relative.
#  (4) fixed source — bitwise `effective_source_m` inside, EXACTLY zero outside, massive form
#      == M·S per element.
#  (5) charts + precision — both charts, the printed z² chain rule, BigFloat(256).
#  (6) guards — the sealed ϱ² = 0 / |z| = 1 throws, the z² equator chart-Jacobian throw, the
#      construction contracts, and the puncture-snap CORNER property on the production mesh.
#  (7) CRP-domain DG residual with the assembled jumps + sources: an ENVELOPE with a p-trend
#      and a no-hook falsifier (the puncture is log-singular, main.tex:522-523, so this is a
#      discretization measurement, not a machine-zero identity).
using Test
using ScalarSelfForce
using ScalarSelfForce.EllipticCore
using ScalarSelfForce.EllipticCore: refelement, element_matrices

# ---- independent transcriptions (parse-free literals) -------------------------------------
# The production lattice (anc/ScalarSelfForce.yaml:55-110 at the header example yaml:9-12):
# radial bounds {2, 2.8, 5.2, 14.8, 24.4, 1e4}, z² rows {0, 0.25, 1}, per-radial-column
# refinement (L, L) (L, L) (L+1, L+1) (L, L+1) (L+1, L). At L = 0 the worldtube column (i = 2)
# therefore holds 2×2 elements, its LEFT neighbour column (i = 1) a single element (a 2:1
# h-nonconforming interface), its RIGHT neighbour column (i = 3) 1×2 and the polar row above
# it 2×2 — all conforming. That fixes the worldtube-boundary faces by hand:
const WT_BLK52 = 2                                   # 5×2 chart: lattice cell (i, j) = (2, 0)
const WT_BLK53 = 7                                   # 5×3 chart: cell (2, 1) = 2 + 5·1
wt_el(b, lx, ix, ly, iy) = ElementId(b, SegmentId(lx, ix), SegmentId(ly, iy))

# 5×2, L = 0: 6 inside-owned + 5 outside-owned worldtube-boundary faces (element, dim, sidx)
const WT_FACES52 = Set([
    (wt_el(2, 1, 0, 1, 0), 1, 1), (wt_el(2, 1, 0, 1, 1), 1, 1),   # r = 5.2, inside
    (wt_el(2, 1, 1, 1, 0), 1, 2), (wt_el(2, 1, 1, 1, 1), 1, 2),   # r = 14.8, inside
    (wt_el(2, 1, 0, 1, 1), 2, 2), (wt_el(2, 1, 1, 1, 1), 2, 2),   # z² = 0.25, inside
    (wt_el(1, 0, 0, 0, 0), 1, 2),                                 # r = 5.2, outside (SPLIT)
    (wt_el(3, 0, 0, 1, 0), 1, 1), (wt_el(3, 0, 0, 1, 1), 1, 1),   # r = 14.8, outside
    (wt_el(7, 1, 0, 1, 0), 2, 1), (wt_el(7, 1, 1, 1, 0), 2, 1)])  # z² = 0.25, outside
# 5×3, L = 0: the worldtube block 7 spans z ∈ [−0.5, 0.5] (2×2 elements), so BOTH angular
# faces are worldtube boundaries and the equator z = 0 is an interior face INSIDE the tube.
const WT_FACES53 = Set([
    (wt_el(7, 1, 0, 1, 0), 1, 1), (wt_el(7, 1, 0, 1, 1), 1, 1),   # r = 5.2, inside
    (wt_el(7, 1, 1, 1, 0), 1, 2), (wt_el(7, 1, 1, 1, 1), 1, 2),   # r = 14.8, inside
    (wt_el(7, 1, 0, 1, 0), 2, 1), (wt_el(7, 1, 1, 1, 0), 2, 1),   # z = −0.5, inside
    (wt_el(7, 1, 0, 1, 1), 2, 2), (wt_el(7, 1, 1, 1, 1), 2, 2),   # z = +0.5, inside
    (wt_el(6, 0, 0, 0, 0), 1, 2),                                 # r = 5.2, outside (SPLIT)
    (wt_el(8, 0, 0, 1, 0), 1, 1), (wt_el(8, 0, 0, 1, 1), 1, 1),   # r = 14.8, outside
    (wt_el(2, 1, 0, 1, 1), 2, 2), (wt_el(2, 1, 1, 1, 1), 2, 2),   # z = −0.5, outside
    (wt_el(12, 1, 0, 1, 0), 2, 1), (wt_el(12, 1, 1, 1, 0), 2, 1)])  # z = +0.5, outside

# `block_domain` requires the inner radial bound to BE r₊ (main.tex:529), which the yaml's
# inner_radius = 2 satisfies only at a = 0 — every spin variation moves it with the horizon.
wt_cfg(; spin = 0.0, r0 = 10.0, m = 2, kwargs...) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, kwargs...)
wt_orbit(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)

# a RunConfig twin with ImposeEquatorialSymmetry = false (crp_paper_config hard-codes the
# yaml's True and exposes no keyword — the same field-wise reconstruction test_domain_blocks
# uses for the 5×3 chart)
function wt_full_z(cfg::RunConfig{T}) where {T}
    OC = typeof(cfg.orbit)
    oc = OC((f === :impose_equatorial_symmetry ? false : getfield(cfg.orbit, f)
             for f in fieldnames(OC))...)
    return RunConfig{T}((f === :orbit ? oc : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

# transmissive (no-op) exterior of dg_operator.jl convention N7: both face corrections vanish
# identically, so a CRP-domain residual measures ONLY the volume terms and the worldtube jumps
# — no boundary-condition contamination (node `dg-boundary-conditions` owns the real ghosts;
# this is a TEST instrument).
struct WTTransmissive <: EllipticCore.GhostBoundary end
EllipticCore.ghost_state(::WTTransmissive, int::FaceState) = FaceState(int.psi, -int.nflux)
EllipticCore.ghost_state_linearized(::WTTransmissive, int::FaceState) =
    FaceState(int.psi, -int.nflux)

# independent enumeration of the worldtube-boundary faces straight off the Mesh2D topology
function wt_expected_faces(dom)
    mesh = dom.mesh
    wtb = Set(worldtube_blocks(dom))
    isin(e) = mesh.elements[e].block in wtb
    out = Tuple{ElementId,Int,Int}[]
    for e in eachindex(mesh.elements), d in 1:2, sidx in 1:2
        iface = mesh.faces[e][d][sidx]
        iface === nothing && continue
        cross = [isin(iface.nbs[k]) != isin(e) for k in 1:length(iface.nbs)]
        any(cross) || continue
        @test all(cross)                     # a face never straddles the partition partially
        push!(out, (mesh.elements[e], d, sidx))
    end
    return out
end

# ---- (3) manufactured singular/regular split, driven through the SHIPPED assembly ---------
# `variant` ∈ (:conforming, :twoone); `flip` ∈ (:none, :psi, :flux, :source) injects a sign
# error into the DATA (never into the assembly) to calibrate the falsifiers.
function wt_manufactured(variant::Symbol; flip::Symbol = :none, Npt = 6, C = 1.5)
    blk = Block(90, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    if variant === :conforming
        ids = vec([ElementId(90, SegmentId(2, ix), SegmentId(2, iy)) for ix in 0:3, iy in 0:3])
        insideset = Set([ElementId(90, SegmentId(2, 1), SegmentId(2, 1))])
    else                                  # a 2:1 h-nonconforming worldtube boundary
        ids = [ElementId(90, SegmentId(1, 0), SegmentId(0, 0)),
               ElementId(90, SegmentId(1, 1), SegmentId(1, 0)),
               ElementId(90, SegmentId(1, 1), SegmentId(1, 1))]
        insideset = Set(ids[2:3])
    end
    mesh = Mesh2D(blk, ids, (Npt, Npt))
    isinside(el::ElementId) = el in insideset

    A1, A2 = 1.3, 0.7                     # 𝔄 diagonal and real, as in the physical operator
    bet, g1, g2 = 0.6 - 0.4im, 0.25 + 0.15im, -0.35 + 0.05im
    Af = ((x, y) -> A1, (x, y) -> A2)
    bf = (x, y) -> bet
    gf = ((x, y) -> g1, (x, y) -> g2)

    U(x, y)   = (0.8 + 0.5im) * x^3 + (0.3 - 0.2im) * x^2 * y + (0.4 + 0.1im) * y^3 -
                (0.6 - 0.3im) * x * y + (0.2 + 0.7im) * x + 0.35 - 0.15im
    Ux(x, y)  = 3(0.8 + 0.5im) * x^2 + 2(0.3 - 0.2im) * x * y - (0.6 - 0.3im) * y + (0.2 + 0.7im)
    Uy(x, y)  = (0.3 - 0.2im) * x^2 + 3(0.4 + 0.1im) * y^2 - (0.6 - 0.3im) * x
    Uxx(x, y) = 6(0.8 + 0.5im) * x + 2(0.3 - 0.2im) * y
    Uyy(x, y) = 6(0.4 + 0.1im) * y
    # the manufactured puncture. In the 2:1 variant it depends on the NORMAL coordinate only,
    # so the split face's two mortars carry IDENTICAL data and one FixedDataHook is exact —
    # the sealed hook seam binds per face, not per mortar.
    P, Px, Py, Pxx, Pyy = if variant === :twoone
        ((x, y) -> (0.45 - 0.25im) * x^3 - (0.7 + 0.2im) * x^2 + (0.5 + 0.6im) * x + 0.2 + 0.1im,
         (x, y) -> 3(0.45 - 0.25im) * x^2 - 2(0.7 + 0.2im) * x + (0.5 + 0.6im),
         (x, y) -> 0.0 + 0.0im,
         (x, y) -> 6(0.45 - 0.25im) * x - 2(0.7 + 0.2im),
         (x, y) -> 0.0 + 0.0im)
    else
        ((x, y) -> (0.45 - 0.25im) * x^3 - (0.7 + 0.2im) * x * y^2 + (0.5 + 0.6im) * y + 0.2 + 0.1im,
         (x, y) -> 3(0.45 - 0.25im) * x^2 - (0.7 + 0.2im) * y^2,
         (x, y) -> -2(0.7 + 0.2im) * x * y + (0.5 + 0.6im),
         (x, y) -> 6(0.45 - 0.25im) * x,
         (x, y) -> -2(0.7 + 0.2im) * x)
    end
    LU(x, y) = -A1 * Uxx(x, y) - A2 * Uyy(x, y) + bet * U(x, y) + g1 * Ux(x, y) + g2 * Uy(x, y)
    LP(x, y) = -A1 * Pxx(x, y) - A2 * Pyy(x, y) + bet * P(x, y) + g1 * Px(x, y) + g2 * Py(x, y)

    sp = flip === :psi ? -1 : 1
    sf = flip === :flux ? -1 : 1
    ss = flip === :source ? -1 : 1
    vals(x, y) = (psi = sp * P(x, y), flux_1 = sf * A1 * Px(x, y), flux_2 = sf * A2 * Py(x, y))
    jd = worldtube_jump_data(mesh, isinside, vals)
    hooks = worldtube_hooks(jd)
    # b = M·(background source) + the worldtube fixed source −△Ψ^𝒫_man inside
    b = worldtube_massive_source(mesh, el -> true, LU, ComplexF64) .+
        worldtube_massive_source(mesh, isinside, (x, y) -> -ss * LP(x, y), ComplexF64)

    ems = [element_matrices(ComplexF64, blk, mesh.elements[e], mesh.orders[e])
           for e in eachindex(mesh.elements)]
    rng(e) = mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])
    u = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        f = isinside(el) ? (x, y) -> U(x, y) - P(x, y) : U
        u[rng(e)] .= vec([f(ems[e].x[1][i], ems[e].x[2][j])
                          for i in 1:mesh.orders[e][1], j in 1:mesh.orders[e][2]])
    end
    bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}()
    for (el, d, sidx) in boundary_faces(mesh)
        e = mesh.index[el]
        lo, hi = logical_interval(Float64, el, d)
        xf = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
        f = isinside(el) ? (x, y) -> U(x, y) - P(x, y) : U
        bcs[(el, d, sidx)] = DirichletGhost(ComplexF64[d == 1 ? f(xf, t) : f(t, xf)
                                                       for t in ems[e].x[3-d]])
    end
    op = MeshOperator(ComplexF64, mesh, Af, bf, gf; C = C, hooks = hooks, bcs = bcs)
    y = similar(u)
    apply_full!(y, op, u)
    return maximum(abs, y .- b) / maximum(abs, b), length(jd), sum(length, values(jd))
end

# ---- (7) CRP-domain residual of the reconstructed-ZERO field ------------------------------
function wt_crp_residual(Npt; spin = 0.0, m = 2, hooks_on = true)
    cfg = wt_cfg(spin = spin, m = m)
    orb = wt_orbit(cfg)
    dom = block_domain(orb, cfg; orders = (Npt, Npt))
    mesh = dom.mesh
    p = Puncture(orb)
    wt = worldtube_data(dom, orb, cfg)
    coef = domain_coefficients(dom, orb, cfg)
    b = worldtube_massive_source(wt)
    u = zeros(eltype(b), mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        inside_worldtube(dom, el) || continue
        em = element_matrices(ComplexF64, mesh.blocks[mesh.blockindex[el.block]], el,
                              mesh.orders[e])
        rng = mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])
        u[rng] .= vec([-worldtube_jump_values(p, m, dom.chart, em.x[1][i], em.x[2][j]).psi
                       for i in 1:mesh.orders[e][1], j in 1:mesh.orders[e][2]])
    end
    hooks = hooks_on ? worldtube_hooks(wt; skip_split = true) :
                       Dict{Tuple{ElementId,Int,Int},BoundaryDataHook}()
    bcs = Dict(f => WTTransmissive() for f in boundary_faces(mesh))
    op = MeshOperator(ComplexF64, mesh, coef.A, coef.beta, coef.gamma;
                      C = cfg.dg.penalty_parameter, hooks = hooks, bcs = bcs)
    y = similar(u)
    apply_full!(y, op, u)
    r = y .- b
    # the one element that owns an unservable SPLIT worldtube face carries no jump data at
    # all through the sealed per-face seam, so it is excluded and reported separately
    split = Set{Int}()
    for (el, d, sidx) in worldtube_boundary_faces(dom)
        e = mesh.index[el]
        length(mesh.faces[e][d][sidx].nbs) > 1 && push!(split, e)
    end
    rin = 0.0; rout = 0.0; rsplit = 0.0
    for e in eachindex(mesh.elements)
        v = maximum(abs, view(r, mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])))
        if e in split
            rsplit = max(rsplit, v)
        elseif inside_worldtube(dom, mesh.elements[e])
            rin = max(rin, v)
        else
            rout = max(rout, v)
        end
    end
    return (rin = rin, rout = rout, rsplit = rsplit, nsplit = length(split))
end

@testset "Worldtube formulation: region, jump data, fixed source (worldtube-formulation)" begin

    # =====================================================================================
    @testset "(1) regularized-region partition of the CRP lattice" begin
        cfg = wt_cfg()
        orb = wt_orbit(cfg)
        dom = block_domain(orb, cfg)
        wt = worldtube_data(dom, orb, cfg)
        @test wt isa WorldtubeData

        # the region is the whole block containing the particle (spectre puncture_in_element)
        @test worldtube_blocks(dom) == [WT_BLK52]
        els = worldtube_elements(wt)
        @test length(els) == 4                                   # (L+1, L+1) ⇒ 2×2 elements
        @test Set(els) == Set(vec([wt_el(2, 1, ix, 1, iy) for ix in 0:1, iy in 0:1]))
        @test els == filter(e -> e.block == WT_BLK52, dom.mesh.elements)   # layout order
        @test worldtube_elements(dom) == els
        for el in dom.mesh.elements
            @test inside_worldtube(dom, el) == (el.block == WT_BLK52)
            @test inside_worldtube(wt, el) == (el.block == WT_BLK52)
        end

        # the worldtube boundary = exactly the interior faces between inside and outside
        faces = worldtube_boundary_faces(dom)
        @test Set(faces) == WT_FACES52
        @test length(faces) == 11
        @test Set(faces) == Set(wt_expected_faces(dom))
        @test faces == sort(faces, by = f -> (dom.mesh.index[f[1]], f[2], f[3]))
        # none of them is a domain boundary; the tube's equator faces ARE domain boundaries
        bf = Set(boundary_faces(dom.mesh))
        @test isempty(intersect(Set(faces), bf))
        for ix in 0:1
            f = (wt_el(2, 1, ix, 1, 0), 2, 1)
            @test f in bf && boundary_class(dom, f) === :equator
        end
        # mortar counts: one 2:1 SPLIT face (the coarse left neighbour), the rest conforming
        jd = worldtube_jump_data(wt)
        @test Set(keys(jd)) == WT_FACES52
        @test sum(length, values(jd)) == 12
        @test count(v -> length(v) == 2, values(jd)) == 1
        @test length(jd[(wt_el(1, 0, 0, 0, 0), 1, 2)]) == 2

        # the same in the full-z 5×3 chart
        cfz = wt_full_z(cfg)
        domz = block_domain(orb, cfz)
        wtz = worldtube_data(domz, orb, cfz)
        @test worldtube_blocks(domz) == [WT_BLK53]
        @test length(worldtube_elements(wtz)) == 4
        @test Set(worldtube_boundary_faces(domz)) == WT_FACES53
        @test length(worldtube_boundary_faces(domz)) == 15
        @test sum(length, values(worldtube_jump_data(wtz))) == 16
        # the equator is an INTERIOR face of the tube in the full chart — no jump there
        @test !((wt_el(7, 1, 0, 1, 0), 2, 2) in Set(worldtube_boundary_faces(domz)))
    end

    # =====================================================================================
    @testset "(2) Ψ^𝒫 jump data on the worldtube faces" begin
        cfg = wt_cfg()
        orb = wt_orbit(cfg)
        dom = block_domain(orb, cfg)
        p = Puncture(orb)
        m = cfg.orbit.m_mode_number
        wt = worldtube_data(dom, orb, cfg)
        jd = worldtube_jump_data(wt)

        # (2a) every slot reproduces the SEALED puncture bundle at the mortar points
        worst_psi = 0.0
        worst_flx = 0.0
        npts = 0
        for (face, mortars) in jd
            el, d, sidx = face
            inside = inside_worldtube(dom, el)
            sg = worldtube_jump_signs(inside, sidx)
            xf, pts = worldtube_face_points(dom, face)
            @test length(pts) == length(mortars)
            for (k, xt) in enumerate(pts)
                fs = mortars[k]
                @test fs isa FaceState
                @test length(fs.psi) == length(xt) == length(fs.nflux)
                for (i, t) in enumerate(xt)
                    r, w = d == 1 ? (xf, t) : (t, xf)
                    z = sqrt(w)
                    b = puncture_psi_chart(p, m, r, z)
                    Fn = d == 1 ? flux_coeff_r(orb, r) * b.dr :
                                  flux_coeff_z2(orb, r, z) * b.dz / (2z)
                    worst_psi = max(worst_psi, abs(fs.psi[i] - sg.psi * b.psi))
                    worst_flx = max(worst_flx, abs(fs.nflux[i] - sg.flux * Fn) / abs(Fn))
                    npts += 1
                end
            end
        end
        println("worldtube-formulation jump data vs the sealed puncture bundle over $npts ",
                "mortar points: ψ ", worst_psi, " (bitwise), flux rel ", worst_flx)
        @test worst_psi == 0.0
        @test worst_flx ≤ 1e-13
        @test npts == 12 * 4

        # (2b) two-sided consistency: ψ flips, n·𝓕 does NOT (the sender-normal convention of
        # main.tex:606-608); calibrated BITWISE over all 12 mortar pairs
        mesh = dom.mesh
        pairs = 0
        for face in worldtube_boundary_faces(dom)
            el, d, sidx = face
            iface = mesh.faces[mesh.index[el]][d][sidx]
            for k in 1:length(iface.nbs)
                nface = (mesh.elements[iface.nbs[k]], d, 3 - sidx)
                a = jd[face][k]
                b = jd[nface][iface.nbslot[k]]
                @test a.psi == -b.psi                          # BITWISE
                @test a.nflux == b.nflux                       # BITWISE
                pairs += 1
            end
        end
        @test pairs == 12

        # (2c) the SIGN TABLE (spectre Equations.cpp:80-91): the receiver solving for Ψ^ℛ
        # SUBTRACTS Ψ^𝒫, the one solving for Ψ adds it; the flux shift is n_in·𝓕(Ψ^𝒫) with
        # n_in the INSIDE element's outward normal, hence the SAME on both sides
        @test worldtube_jump_signs(true, 1) == (psi = -1, flux = -1)
        @test worldtube_jump_signs(true, 2) == (psi = -1, flux = +1)
        @test worldtube_jump_signs(false, 1) == (psi = +1, flux = +1)
        @test worldtube_jump_signs(false, 2) == (psi = +1, flux = -1)
        for (face, mortars) in jd
            el, d, sidx = face
            xf, pts = worldtube_face_points(dom, face)
            r, w = d == 1 ? (xf, pts[1][1]) : (pts[1][1], xf)
            ref = puncture_psi_chart(p, m, r, sqrt(w)).psi
            @test mortars[1].psi[1] == (inside_worldtube(dom, el) ? -ref : ref)
        end

        # (2d) EXACT interface cancellation of the reconstructed field: with u = −Ψ^𝒫 inside
        # and u = 0 outside the reconstructed Ψ = u + Ψ^𝒫 vanishes identically, so on every
        # worldtube mortar point the hook-shifted exterior datum equals the matched trace of
        # the receiver's own variable — every DG interface correction vanishes, BITWISE.
        for (face, mortars) in jd
            el, d, sidx = face
            inside = inside_worldtube(dom, el)
            s = sidx == 1 ? -1 : 1
            xf, pts = worldtube_face_points(dom, face)
            for (k, xt) in enumerate(pts), (i, t) in enumerate(xt)
                r, w = d == 1 ? (xf, t) : (t, xf)
                z = sqrt(w)
                b = puncture_psi_chart(p, m, r, z)
                Fn = d == 1 ? flux_coeff_r(orb, r) * b.dr :
                              flux_coeff_z2(orb, r, z) * b.dz / (2z)
                own_psi = inside ? -b.psi : zero(b.psi)          # ψ_int
                own_phi = inside ? s * (-Fn) : zero(Fn)          # φ_int (receiver normal)
                sent_psi = inside ? zero(b.psi) : -b.psi
                sent_phi = inside ? zero(Fn) : (-s) * (-Fn)      # sender convention (N3)
                @test sent_psi + mortars[k].psi[i] == own_psi
                @test sent_phi + mortars[k].nflux[i] == -own_phi
            end
        end

        # (2e) other (a, m) points and a nonzero particle azimuth
        for (spin, mm, phip) in ((0.5, 1, 0.0), (0.5, 2, 0.31), (0.0, 1, 0.0))
            c2 = wt_cfg(spin = spin, m = mm)
            o2 = wt_orbit(c2)
            d2 = block_domain(o2, c2)
            p2 = Puncture(o2)
            w2 = worldtube_data(d2, o2, c2; phi_p = phip)
            j2 = worldtube_jump_data(w2)
            @test length(j2) == 11
            face = (wt_el(2, 1, 0, 1, 0), 1, 1)
            xf, pts = worldtube_face_points(d2, face)
            ref = puncture_psi_chart(p2, mm, xf, sqrt(pts[1][1]); phi_p = phip).psi
            @test j2[face][1].psi[1] == -ref
            @test ref != puncture_psi_chart(p2, mm, xf, sqrt(pts[1][1]); phi_p = phip + 0.7).psi
        end
    end

    # =====================================================================================
    @testset "(3) manufactured singular/regular split — DG residual (row-38 gate)" begin
        # A two-region 2D mesh, affine blocks, complex β/γ; smooth manufactured Ψ_man and
        # Ψ^𝒫_man of degree ≤ 3 (inside the N = 6 element space, so the strong-form DG
        # residual of the exact piecewise field is M·S up to roundoff). The SHIPPED assembly
        # (face enumeration, mortar geometry, sign table, hook packing, massive source) must
        # annihilate it.
        rel, nf, nm = wt_manufactured(:conforming)
        fals = [wt_manufactured(:conforming; flip = f)[1] for f in (:psi, :flux, :source)]
        println("worldtube-formulation manufactured conforming split (", nf, " faces / ",
                nm, " mortars): |R−b|∞/|b|∞ = ", rel, " ; falsifiers ", fals)
        @test nf == 8 && nm == 8                            # 4 inside-owned + 4 outside-owned
        @test rel ≤ 1e-10                                   # calibrated 1.3e-12
        @test minimum(fals) > 1e-3                          # calibrated 0.29

        rel2, nf2, nm2 = wt_manufactured(:twoone)
        fals2 = [wt_manufactured(:twoone; flip = f)[1] for f in (:psi, :flux, :source)]
        println("worldtube-formulation manufactured 2:1 split (", nf2, " faces / ", nm2,
                " mortars): |R−b|∞/|b|∞ = ", rel2, " ; falsifiers ", fals2)
        @test nf2 == 3 && nm2 == 4                          # the coarse face carries 2 mortars
        @test rel2 ≤ 1e-10                                  # calibrated 1.9e-13
        @test minimum(fals2) > 1e-3                         # calibrated 0.61
    end

    # =====================================================================================
    @testset "(4) fixed-source assembly" begin
        cfg = wt_cfg()
        orb = wt_orbit(cfg)
        dom = block_domain(orb, cfg)
        p = Puncture(orb)
        m = cfg.orbit.m_mode_number
        wt = worldtube_data(dom, orb, cfg)
        S = worldtube_source(wt)
        Mb = worldtube_massive_source(wt)
        mesh = dom.mesh
        @test length(S) == mesh.ndofs == length(Mb)
        nin = 0
        for e in eachindex(mesh.elements)
            el = mesh.elements[e]
            rng = mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])
            if inside_worldtube(dom, el)
                em = element_matrices(ComplexF64, mesh.blocks[mesh.blockindex[el.block]],
                                      el, mesh.orders[e])
                exp = vec([effective_source_m(p, m, em.x[1][i], sqrt(em.x[2][j]))
                           for i in 1:mesh.orders[e][1], j in 1:mesh.orders[e][2]])
                @test S[rng] == exp                             # BITWISE, i-fastest layout
                @test Mb[rng] == em.M * exp                     # spectre apply_mass_matrix
                nin += length(rng)
            else
                @test all(iszero, S[rng])                       # EXACTLY zero outside
                @test all(iszero, Mb[rng])
            end
        end
        @test nin == 4 * 16                                     # 4 elements × 4×4 points
        @test !all(iszero, S)
    end

    # =====================================================================================
    @testset "(5) charts, precision, and the z² chain rule" begin
        cfg = wt_cfg()
        orb = wt_orbit(cfg)
        p = Puncture(orb)
        m = cfg.orbit.m_mode_number
        # the z² angular flux is the printed chain rule 𝓕^{z²} = A^{z²}∂_wΨ, ∂_w = ∂_z/(2z)
        for (r, z) in ((7.3, 0.31), (12.5, 0.5), (6.1, 0.12))
            vz = worldtube_jump_values(p, m, :z, r, z)
            vw = worldtube_jump_values(p, m, :z2, r, z^2)
            b = puncture_psi_chart(p, m, r, z)
            @test vz.psi == b.psi
            @test vw.psi == b.psi
            @test vz.flux_r == flux_coeff_r(orb, r) * b.dr
            @test vw.flux_r == vz.flux_r
            @test vz.flux_q == flux_coeff_z(orb, r, z) * b.dz
            @test isapprox(vw.flux_q, 2z * vz.flux_q; rtol = 1e-14)
            @test vz.source == effective_source_m(p, m, r, z)
            @test vw.source == vz.source
        end
        # full-z domain data
        cfz = wt_full_z(cfg)
        domz = block_domain(orb, cfz)
        wtz = worldtube_data(domz, orb, cfz)
        for (face, mortars) in worldtube_jump_data(wtz)
            el, d, sidx = face
            xf, pts = worldtube_face_points(domz, face)
            sg = worldtube_jump_signs(inside_worldtube(domz, el), sidx)
            for (k, xt) in enumerate(pts), (i, t) in enumerate(xt)
                r, z = d == 1 ? (xf, t) : (t, xf)
                b = puncture_psi_chart(p, m, r, z)
                Fn = d == 1 ? flux_coeff_r(orb, r) * b.dr : flux_coeff_z(orb, r, z) * b.dz
                @test mortars[k].psi[i] == sg.psi * b.psi
                @test mortars[k].nflux[i] == sg.flux * Fn
            end
        end
        # BigFloat(256) inheritance
        setprecision(BigFloat, 256) do
            rp = CircularOrbit(big"0.5", big"10.0").rp
            cb = crp_paper_config(spin = big"0.5", orbital_radius = big"10.0", m_mode = 2,
                                  inner_radius = rp, outer_radius = big"1.0e4",
                                  worldtube_radius = big"4.8", worldtube_dz2 = big"0.25")
            ob = wt_orbit(cb)
            db = block_domain(ob, cb)
            wb = worldtube_data(db, ob, cb)
            @test wb isa WorldtubeData{Complex{BigFloat},BigFloat}
            Sb = worldtube_source(wb)
            @test eltype(Sb) == Complex{BigFloat}
            @test all(v -> precision(real(v)) == 256, Sb)
            fs = first(values(worldtube_jump_data(wb)))[1]
            @test eltype(fs.psi) == Complex{BigFloat}
            @test precision(real(fs.psi[1])) == 256
            @test length(worldtube_boundary_faces(db)) == 11
            @test eltype(worldtube_massive_source(wb)) == Complex{BigFloat}
        end
    end

    # =====================================================================================
    @testset "(6) domain guards and the puncture-snap corner property" begin
        cfg = wt_cfg()
        orb = wt_orbit(cfg)
        dom = block_domain(orb, cfg)
        p = Puncture(orb)
        m = cfg.orbit.m_mode_number
        # the sealed ϱ² = 0 / |z| = 1 guards propagate (no silent NaN into the fixed source)
        @test_throws DomainError worldtube_jump_values(p, m, :z, orb.r0, 0.0)
        @test_throws DomainError worldtube_jump_values(p, m, :z2, orb.r0, 0.0)
        @test_throws DomainError worldtube_jump_values(p, m, :z, 9.0, 1.0)
        # the z² chart Jacobian ∂_w = ∂_z/(2z) is singular at the equator w = 0, which is a
        # DOMAIN boundary of the half-domain chart and never a worldtube face
        @test_throws DomainError worldtube_jump_values(p, m, :z2, 9.0, 0.0)
        @test worldtube_jump_values(p, m, :z, 9.0, 0.0).flux_q isa Complex   # legal in z
        @test_throws ArgumentError worldtube_jump_values(p, m, :bogus, 9.0, 0.1)

        # node `amr-driver`'s puncture snap must place the particle on an element CORNER:
        # Legendre-Gauss nodes are strictly interior, so no volume node and no mortar point
        # can coincide with it. Asserted on the production mesh, where r₀ = 10 IS an element
        # bound of the worldtube column and the equator w = 0 IS a block bound.
        mesh = dom.mesh
        dmin = Inf
        for e in eachindex(mesh.elements)
            el = mesh.elements[e]
            inside_worldtube(dom, el) || continue
            em = element_matrices(Float64, mesh.blocks[mesh.blockindex[el.block]], el,
                                  mesh.orders[e])
            for j in 1:mesh.orders[e][2], i in 1:mesh.orders[e][1]
                dmin = min(dmin, hypot(em.x[1][i] - orb.r0, sqrt(em.x[2][j])))
            end
        end
        println("worldtube-formulation closest production volume node to the puncture: ", dmin)
        @test dmin > 1e-2
        for face in worldtube_boundary_faces(dom)
            xf, pts = worldtube_face_points(dom, face)
            for xt in pts, t in xt
                r, w = face[2] == 1 ? (xf, t) : (t, xf)
                @test hypot(r - orb.r0, sqrt(w)) > 1e-2
            end
        end

        # construction contracts
        @test_throws ArgumentError worldtube_data(dom, orb, cfg, 3)      # m ≠ config's m
        @test_throws ArgumentError worldtube_data(dom, orb, wt_full_z(cfg))   # chart mismatch
        # the sealed FixedDataHook seam binds per FACE, not per mortar: the production
        # lattice's 2:1 worldtube face carries genuinely different data on its two mortars, so
        # `worldtube_hooks` REFUSES it and names the extension point (node dg-jump-conditions)
        wt = worldtube_data(dom, orb, cfg)
        @test_throws ArgumentError worldtube_hooks(wt)
        hk = worldtube_hooks(wt; skip_split = true)
        @test length(hk) == 10
        @test all(h -> h isa FixedDataHook, values(hk))
        @test !haskey(hk, (wt_el(1, 0, 0, 0, 0), 1, 2))
        # allocation-free per-point evaluation (the sealed bundle is 0-byte)
        worldtube_jump_values(p, m, :z2, 9.0, 0.1)
        @test (@allocated worldtube_jump_values(p, m, :z2, 9.0, 0.1)) == 0
    end

    # =====================================================================================
    @testset "(7) CRP-domain residual with the assembled jumps and sources" begin
        # u = −Ψ^𝒫 inside, 0 outside ⇒ the reconstructed Ψ ≡ 0 solves the exterior equation
        # exactly and the interior equation with source S^eff. The DG residual of that field
        # is then PURELY the discretization error of the log-singular puncture
        # (Ψ^𝒫_m ~ log λ, main.tex:522-523), which is why this gate is an ENVELOPE with a
        # p-trend and a falsifier, not a machine-zero identity.
        a4 = wt_crp_residual(4)
        a8 = wt_crp_residual(8)
        n4 = wt_crp_residual(4; hooks_on = false)
        println("worldtube-formulation CRP residual N=4: in ", a4.rin, " out ", a4.rout,
                " (split-owner ", a4.rsplit, ", ", a4.nsplit, " excluded) | N=8: in ",
                a8.rin, " out ", a8.rout, " | no-hook N=4: in ", n4.rin, " out ", n4.rout)
        @test a4.nsplit == 1
        @test a4.rin ≤ 2e-2 && a8.rin < a4.rin              # calibrated 8.21e-3 → 2.93e-3
        @test a4.rout ≤ 2e-3 && a8.rout < a4.rout           # calibrated 4.84e-4 → 1.45e-4
        @test n4.rin / a4.rin > 20                          # FALSIFIER: calibrated 71.6
        @test n4.rout / a4.rout > 100                       # calibrated 1331
    end
end

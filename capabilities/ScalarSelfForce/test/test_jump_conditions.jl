# test_jump_conditions.jl — gates for Solve/jump_glue.jl (ESSF node `dg-jump-conditions`,
# schedule W9-A, GLOBAL_DAG row 40, commit 2 — the ScalarSelfForce half). Written BEFORE the
# implementation (TDD red run in evidence/dg-jump-conditions/tests_red_ssf.log).
#
# WHAT THE NODE DELIVERS. The `MeshOperator(…; hooks = …)` table of the CRP problem: the two
# jump conditions of arXiv:2606.04998 §III.A.5 (main.tex:596-608) bound onto the EllipticCore
# ModifyBoundaryData seam, over the per-MORTAR widening this node's commit 1 added.
#
#   * WORLDTUBE (eq:jump_worldtube, main.tex:600-604) — Ψ^ext → Ψ^ext ± Ψ^𝒫,
#     (n·𝓕)^ext → (n·𝓕)^ext ∓ n·𝔄·∂Ψ^𝒫. Trial-field INDEPENDENT ⇒ `FixedDataHook`s, which
#     never enter the linearized matrix and reach the RHS through `operator_on_zero`. The
#     DATA are the sealed row-38 `worldtube_jump_data(wt)` (one `FaceState` per mortar); this
#     node only binds them.
#   * VTU JUNCTIONS (eq:jump_vt/eq:jump_tu, main.tex:337-343, realized main.tex:606-608) —
#     (n·𝓕)^ext → (n·𝓕)^ext − imΩ⟨⟨Ψ⟩⟩ with ⟨⟨Ψ⟩⟩ = (Ψ^int + Ψ^ext)/2, the SAME sign on both
#     sides of both transition radii (spectre `ModifyBoundaryData::apply_linearized`,
#     Equations.cpp:92-116). LINEAR-HOMOGENEOUS in the trial field ⇒ a `LinearizedHook`, i.e.
#     the term IS in the GMRES matrix. No coefficient is re-derived: `𝓕^r = A^rr ∂_rΨ` IS
#     `∂_{r_*}Ψ` (A^rr = dr/dr_*, main.tex:354-358 — the sealed `VtuJunction` docstring), so
#     the sealed `vtu_junctions(v, orbit, m).coeff = −imΩ` is the whole datum.
#   * COMPOSITION — a face carrying both is `ComposedHook(LinearizedHook(vtu),
#     FixedDataHook(worldtube))`: the W5-B ordering, the LINEARIZED stage acting on the RAW
#     received data FIRST. On the CRP lattice this never happens (gated below: the worldtube
#     block spans [5.2, 14.8] and the transitions sit at 2.8 / 24.4), but the composition is
#     built generically anyway.
#
# ORACLE GATES (calibrated PRE-RED in evidence/dg-jump-conditions/calib_ssf{1,2,3,4}.jl):
#  (b1) 2:1 two-region prescribed-jump reproduction through `jump_hooks` on a manufactured
#       mini-domain: |R_full(u) − b|∞/|b|∞ ≤ 1e-10 (the GLOBAL_DAG row-40 gate), with the
#       REVERSED mortar order as an O(1) falsifier.
#  (b2) a manufactured vtu solution — the sealed test_vtu.jl construction Ψ_dom = Ψ̃e^{−imΩh}
#       with the analytic branch derivatives — solves the modified-flux system SPECTRALLY on a
#       radial 3-block mesh spanning BOTH junctions (a = 0, z-uniform field, so the operator
#       reduces exactly to the radial ODE): calibrated residual 5.67e-2 → 2.21e-4 and solve
#       4.96e-4 → 2.75e-7 over N = 4..10 (fitted log10 slopes −0.403 / −0.552 per N).
#       FALSIFIERS at N = 8: no hooks (solve 1.51e-1), sign-flipped coefficient (2.91e-1),
#       zero coefficient (identical to no hooks, bitwise).
#  (b3) THE ROW-40 ORACLE on the production lattice: `build_matrix` CHANGES when the vtu hooks
#       are added (the junction term is in the GMRES matrix) and is BITWISE UNCHANGED when the
#       worldtube hooks are added, whose entire effect is in `operator_on_zero`.
#  (b4) the CRP table at L ∈ {0, 1}, both charts, a ∈ {0, 0.5}, m ∈ {1, 2}: it covers exactly
#       the worldtube ∪ junction faces, the per-mortar seam serves the 2:1 worldtube faces
#       (1 at L = 0, 2 at L = 1), worldtube ∩ junction = ∅, and the operator over
#       `domain_coefficients` + `boundary_conditions` applies FINITELY. The smooth-puncture
#       control on the production mesh — split-owner element INCLUDED — converges SPECTRALLY:
#       calibrated 4.45e-4 → 5.89e-9 over N = 4..10 (ratios 42.7/42.4/41.7), against 1.47e-1
#       with the mortar order reversed and 1.32e+01 with no hooks.
#  (b5) the REAL puncture: the row-38 envelope is reproduced BITWISE inside, and the
#       previously excluded split-owner element now falls with N like the rest — calibrated
#       2.160e-3 (N = 4) → 3.890e-4 (N = 8), versus 3.374e-2 flat when that face's jump is
#       dropped (all the sealed per-FACE seam could do).
#  (b6) BigFloat(256) construction; the per-LEVEL callable on a coarsened mesh (the row-42
#       contract); m = 0, where the junction coefficient vanishes identically.
#
# Helpers are jc_-prefixed (sibling test files share the runtests.jl include scope) and
# deliberately UNtyped where they touch the new API, so this file LOADS before it exists.

using Test
using LinearAlgebra
using ScalarSelfForce
using EllipticCore

# ---- configs ------------------------------------------------------------------------------
# `block_domain` REQUIRES the inner radial bound to BE r₊ (main.tex:529), which the yaml's
# inner_radius = 2 satisfies only at a = 0 — every spin moves it with the horizon.
jc_cfg(; spin = 0.0, r0 = 10.0, m = 2, L = 0) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, L = L)
jc_orb(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)
# the 5×3 full-z twin (crp_paper_config hard-codes the yaml's ImposeEquatorialSymmetry = true)
function jc_full_z(cfg::RunConfig{T}) where {T}
    OC = typeof(cfg.orbit)
    oc = OC((f === :impose_equatorial_symmetry ? false : getfield(cfg.orbit, f)
             for f in fieldnames(OC))...)
    return RunConfig{T}((f === :orbit ? oc : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

# transmissive (convention-N7) ghost: both face corrections vanish identically, so a residual
# measures ONLY the volume terms and the jump conditions — no boundary-condition contamination
struct JCTransmissive <: EllipticCore.GhostBoundary end
EllipticCore.ghost_state(::JCTransmissive, i::FaceState) = FaceState(i.psi, -i.nflux)
EllipticCore.ghost_state_linearized(::JCTransmissive, i::FaceState) =
    FaceState(i.psi, -i.nflux)

jc_nmortars(mesh, f) = length(mesh.faces[mesh.index[f[1]]][f[2]][f[3]].nbs)

# =========================== (b1) the manufactured mini-domain ==============================
# ONE coarse outside element facing TWO fine inside elements: the 2:1 worldtube boundary in
# miniature, with a manufactured puncture that depends on the TRANSVERSE coordinate so the
# split face's two mortars carry genuinely DIFFERENT data.
function jc_mini(Npt; reversed = false, hooks_on = true)
    blk = Block(91, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    ids = [ElementId(91, SegmentId(1, 0), SegmentId(0, 0)),
           ElementId(91, SegmentId(1, 1), SegmentId(1, 0)),
           ElementId(91, SegmentId(1, 1), SegmentId(1, 1))]
    mesh = Mesh2D(blk, ids, (Npt, Npt))
    inside = Set(ids[2:3])
    isin(el::ElementId) = el in inside
    A1, A2 = 1.3, 0.7
    bet, g1, g2 = 0.6 - 0.4im, 0.25 + 0.15im, -0.35 + 0.05im
    Af = ((x, y) -> A1, (x, y) -> A2)
    U(x, y)   = (0.8 + 0.5im)*x^3 + (0.3 - 0.2im)*x^2*y + (0.4 + 0.1im)*y^3 -
                (0.6 - 0.3im)*x*y + (0.2 + 0.7im)*x + 0.35 - 0.15im
    Ux(x, y)  = 3(0.8 + 0.5im)*x^2 + 2(0.3 - 0.2im)*x*y - (0.6 - 0.3im)*y + (0.2 + 0.7im)
    Uy(x, y)  = (0.3 - 0.2im)*x^2 + 3(0.4 + 0.1im)*y^2 - (0.6 - 0.3im)*x
    Uxx(x, y) = 6(0.8 + 0.5im)*x + 2(0.3 - 0.2im)*y
    Uyy(x, y) = 6(0.4 + 0.1im)*y
    P(x, y)   = (0.45 - 0.25im)*y^3 - (0.7 + 0.2im)*x*y^2 + (0.5 + 0.6im)*y +
                (0.3 - 0.1im)*x^2 + 0.2
    Px(x, y)  = -(0.7 + 0.2im)*y^2 + 2(0.3 - 0.1im)*x
    Py(x, y)  = 3(0.45 - 0.25im)*y^2 - 2(0.7 + 0.2im)*x*y + (0.5 + 0.6im)
    Pxx(x, y) = 2(0.3 - 0.1im) + 0.0im
    Pyy(x, y) = 6(0.45 - 0.25im)*y - 2(0.7 + 0.2im)*x
    LU(x, y) = -A1*Uxx(x, y) - A2*Uyy(x, y) + bet*U(x, y) + g1*Ux(x, y) + g2*Uy(x, y)
    LP(x, y) = -A1*Pxx(x, y) - A2*Pyy(x, y) + bet*P(x, y) + g1*Px(x, y) + g2*Py(x, y)
    vals(x, y) = (psi = P(x, y), f1 = A1 * Px(x, y), f2 = A2 * Py(x, y))
    jd0 = worldtube_jump_data(mesh, isin, vals)
    jd = reversed ? Dict(f => reverse(ms) for (f, ms) in jd0) : jd0
    hooks = hooks_on ?
        jump_hooks(mesh, jd, Dict{Tuple{ElementId,Int,Int},VtuJunction{Float64}}()) :
        Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}()
    b = worldtube_massive_source(mesh, el -> true, LU, ComplexF64) .+
        worldtube_massive_source(mesh, isin, (x, y) -> -LP(x, y), ComplexF64)
    u = zeros(ComplexF64, mesh.ndofs)
    bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}()
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        em = element_matrices(ComplexF64, blk, el, mesh.orders[e])
        f = isin(el) ? (x, y) -> U(x, y) - P(x, y) : U
        nx, ny = mesh.orders[e]
        u[mesh.offsets[e]+1:mesh.offsets[e]+nx*ny] .=
            vec([f(em.x[1][i], em.x[2][j]) for i in 1:nx, j in 1:ny])
    end
    for (el, d, sidx) in boundary_faces(mesh)
        e = mesh.index[el]
        em = element_matrices(ComplexF64, blk, el, mesh.orders[e])
        lo, hi = logical_interval(Float64, el, d)
        xf = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
        f = isin(el) ? (x, y) -> U(x, y) - P(x, y) : U
        bcs[(el, d, sidx)] = DirichletGhost(ComplexF64[d == 1 ? f(xf, t) : f(t, xf)
                                                       for t in em.x[3-d]])
    end
    op = MeshOperator(ComplexF64, mesh, Af, (x, y) -> bet, ((x, y) -> g1, (x, y) -> g2);
                      C = 1.5, hooks = hooks, bcs = bcs)
    y = similar(u)
    apply_full!(y, op, u; threaded = false)
    return (err = maximum(abs, y .- b) / maximum(abs, b),
            nfaces = length(jd), nmort = sum(length, values(jd)), mesh = mesh, jd = jd)
end

# =========================== (b2) the manufactured vtu problem ==============================
# Ψ_dom(r) = Ψ̃(r_*) e^{−imΩ h_dom(r_*)} with h_dom linear (main.tex:324-330), so Ψ is
# continuous and [∂_{r_*}Ψ] = −imΩΨ at BOTH transitions BY CONSTRUCTION (the sealed
# test_vtu.jl branch construction). a = 0 and a z-independent Ψ make the 2D operator reduce
# to the radial ODE exactly: F^z ≡ 0 and β_m is z-free at a = 0.
const JC_RV, JC_RU = 4.0, 12.0
const JC_RIN, JC_ROUT = 3.0, 20.0
const JC_ZLO, JC_ZHI = 0.2, 0.8
const JC_KAP = 1 / 17

function jc_vtu_case(N; m = 2, coeff = nothing, hooks_on = true)
    orb = CircularOrbit(0.0, 10.0)
    vsl = vtu_slicing(orb, JC_RV, JC_RU)
    jn = vtu_junctions(vsl, orb, m)
    Om = orb.Omega
    p0, p1, p2 = complex(0.83, -0.29), complex(0.31, 0.57)/40, complex(-0.12, 0.23)/1600
    Pt(x) = p0 + p1*x + p2*x^2
    Pt1(x) = p1 + 2p2*x
    til(x)  = Pt(x) * cis(JC_KAP*x)
    til1(x) = (Pt1(x) + im*JC_KAP*Pt(x)) * cis(JC_KAP*x)
    til2(x) = (2p2 + 2im*JC_KAP*Pt1(x) - JC_KAP^2*Pt(x)) * cis(JC_KAP*x)
    alp(d) = d === :v ? -1 : d === :u ? 1 : 0
    hc(d) = d === :v ? vsl.rstar_v : d === :u ? -vsl.rstar_u : 0.0
    E(d, x) = cis(-m*Om*(alp(d)*x + hc(d)))
    psi(d, x) = til(x) * E(d, x)
    dpsi(d, x) = (til1(x) - im*m*Om*alp(d)*til(x)) * E(d, x)
    d2psi(d, x) = (til2(x) - 2im*m*Om*alp(d)*til1(x) - (m*Om*alp(d))^2*til(x)) * E(d, x)
    reg(r) = vtu_domain(vsl, r)
    Psi(r, z) = psi(reg(r), rstar(orb, r))
    function Src(r, z)
        d = reg(r); x = rstar(orb, r); Arr = flux_coeff_r(orb, r)
        return -d2psi(d, x)/Arr + beta_m(orb, m, r, z, alp(d), 0)*psi(d, x) +
               gamma_r_m(orb, m, r, alp(d))*dpsi(d, x)/Arr
    end
    fv = fluxform_callables(orb, m, vsl, :v)
    ft = fluxform_callables(orb, m, vsl, :t)
    fu = fluxform_callables(orb, m, vsl, :u)
    pick(r) = (d = reg(r); d === :v ? fv : d === :t ? ft : fu)
    cA, cB = ft.A, (r, z) -> pick(r).beta(r, z)
    cG = ((r, z) -> pick(r).gamma[1](r, z), ft.gamma[2])

    bnds = [JC_RIN, JC_RV, JC_RU, JC_ROUT]
    blocks = [Block(20 + i, LinearMap1D(bnds[i], bnds[i+1]), LinearMap1D(JC_ZLO, JC_ZHI))
              for i in 1:3]
    els = ElementId[]
    for (i, nx) in enumerate((0, 1, 0)), ix in 0:(1<<nx)-1
        push!(els, ElementId(20 + i, SegmentId(nx, ix), SegmentId(0, 0)))
    end
    conns = [BlockConnection(20 + i, 1, 21 + i) for i in 1:2]
    mesh = Mesh2D(blocks, els, (N, N); connections = conns)
    # v and u are the null-sliced blocks; the two junction faces are the block seams
    jfaces = Dict{Tuple{ElementId,Int,Int},VtuJunction{Float64}}()
    ns(b) = b == 21 || b == 23
    for e in eachindex(mesh.elements), d in 1:2, sidx in 1:2
        iface = mesh.faces[e][d][sidx]
        iface === nothing && continue
        el = mesh.elements[e]
        any(k -> ns(mesh.elements[iface.nbs[k]].block) != ns(el.block),
            1:length(iface.nbs)) || continue
        xf, _ = mortar_points(mesh, (el, d, sidx))
        jfaces[(el, d, sidx)] = xf == jn.vt.r ? jn.vt : jn.tu
    end
    coeff !== nothing && (jfaces = Dict(f => VtuJunction(j.r, j.rstar, coeff)
                                        for (f, j) in jfaces))
    hooks = hooks_on ?
        jump_hooks(mesh, Dict{Tuple{ElementId,Int,Int},Vector{FaceState{ComplexF64}}}(),
                   jfaces) :
        Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}()

    u = zeros(ComplexF64, mesh.ndofs)
    b = zeros(ComplexF64, mesh.ndofs)
    bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}()
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        em = element_matrices(ComplexF64, blk, el, mesh.orders[e])
        nx, ny = mesh.orders[e]
        rng = mesh.offsets[e]+1 : mesh.offsets[e]+nx*ny
        u[rng] .= vec([Psi(em.x[1][i], em.x[2][j]) for i in 1:nx, j in 1:ny])
        b[rng] .= em.M * vec([Src(em.x[1][i], em.x[2][j]) for i in 1:nx, j in 1:ny])
    end
    for (el, d, sidx) in boundary_faces(mesh)
        e = mesh.index[el]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        em = element_matrices(ComplexF64, blk, el, mesh.orders[e])
        lo, hi = logical_interval(Float64, el, d)
        xf = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
        bcs[(el, d, sidx)] = d == 2 ? JCTransmissive() :
            DirichletGhost(ComplexF64[Psi(xf, t) for t in em.x[3-d]])
    end
    op = MeshOperator(ComplexF64, mesh, cA, cB, cG; C = 1.5, hooks = hooks, bcs = bcs)
    y = similar(u)
    apply_full!(y, op, u; threaded = false)
    us = direct_solve(op, b .- operator_on_zero(op))
    return (res = maximum(abs, y .- b) / maximum(abs, b),
            serr = maximum(abs, us .- u) / maximum(abs, u),
            op = op, mesh = mesh, njf = length(jfaces))
end

jc_slope(Ns, ys) = (length(Ns)*sum(Ns .* log10.(ys)) - sum(Ns)*sum(log10.(ys))) /
                   (length(Ns)*sum(Ns .^ 2) - sum(Ns)^2)

# ================= (b4)/(b5) the production lattice ========================================
jc_smooth_P(r, w)  = (0.021 - 0.013im)*(r-10)^3 + (0.35 + 0.24im)*w^3 +
                     (-0.11 + 0.07im)*(r-10)*w^2 + (0.4 - 0.2im)*(r-10)*w + (0.7 + 0.3im)
jc_smooth_Pr(r, w) = 3*(0.021 - 0.013im)*(r-10)^2 + (-0.11 + 0.07im)*w^2 + (0.4 - 0.2im)*w
jc_smooth_Pw(r, w) = 3*(0.35 + 0.24im)*w^2 + 2*(-0.11 + 0.07im)*(r-10)*w + (0.4 - 0.2im)*(r-10)

# The SMOOTH-PUNCTURE control on the production mesh: the identical construction as the real
# puncture but with an exactly representable degree-3 manufactured Ψ^𝒫, so what is measured is
# purely the INTERFACE consistency of the assembled table (the split-owner element INCLUDED).
function jc_smooth_control(Npt; spin = 0.0, m = 2, L = 0, mode = :on)
    cfg = jc_cfg(spin = spin, m = m, L = L)
    orb = jc_orb(cfg)
    dom = block_domain(orb, cfg; orders = (Npt, Npt))
    mesh = dom.mesh
    coef = domain_coefficients(dom, orb, cfg)
    isin(el) = inside_worldtube(dom, el)
    vals(r, w) = (psi = jc_smooth_P(r, w), f1 = coef.A[1](r, w)*jc_smooth_Pr(r, w),
                  f2 = coef.A[2](r, w)*jc_smooth_Pw(r, w))
    jd0 = worldtube_jump_data(mesh, isin, vals)
    jd = mode === :reversed ? Dict(f => reverse(ms) for (f, ms) in jd0) : jd0
    jf = vtu_junction_faces(dom, orb, cfg, m; mesh = mesh)
    hooks = mode === :none ? Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}() :
                             jump_hooks(mesh, jd, jf)
    u = zeros(ComplexF64, mesh.ndofs)
    b = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        isin(el) || continue
        blk = mesh.blocks[mesh.blockindex[el.block]]
        em = element_matrices(ComplexF64, blk, el, mesh.orders[e])
        c = element_coefficients(em, blk, el, coef.A, coef.beta, coef.gamma)
        nx, ny = mesh.orders[e]
        ue = vec([-jc_smooth_P(em.x[1][i], em.x[2][j]) for i in 1:nx, j in 1:ny])
        st = element_state(em, c, ue)
        rng = mesh.offsets[e]+1 : mesh.offsets[e]+nx*ny
        u[rng] .= ue
        b[rng] .= em.M * (c.beta .* st.psi .+ c.gamma[1] .* st.dpsi[1] .+
                          c.gamma[2] .* st.dpsi[2]) .- em.MD[1]*st.F[1] .- em.MD[2]*st.F[2]
    end
    bcs = Dict(f => JCTransmissive() for f in boundary_faces(mesh))
    op = MeshOperator(ComplexF64, mesh, coef.A, coef.beta, coef.gamma;
                      C = cfg.dg.penalty_parameter, hooks = hooks, bcs = bcs)
    y = similar(u)
    apply_full!(y, op, u; threaded = false)
    return maximum(abs, y .- b)
end

# The REAL puncture: the row-38 gate-(7) instrument, with the split-owner element now SERVED.
# `mode` ∈ (:served, :dropsplit, :none) — `:dropsplit` omits the split face's jump entirely,
# which is all the sealed per-FACE seam of row 38 could do.
function jc_crp_residual(Npt; spin = 0.0, m = 2, L = 0, mode = :served)
    cfg = jc_cfg(spin = spin, m = m, L = L)
    orb = jc_orb(cfg)
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
    jd = worldtube_jump_data(wt)
    mode === :dropsplit && (jd = Dict(f => ms for (f, ms) in jd if
                                      length(ms) == 1 || all(k -> ms[k].psi == ms[1].psi,
                                                             2:length(ms))))
    hooks = mode === :none ? Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}() :
            jump_hooks(mesh, jd, vtu_junction_faces(dom, orb, cfg, m; mesh = mesh))
    bcs = Dict(f => JCTransmissive() for f in boundary_faces(mesh))
    op = MeshOperator(ComplexF64, mesh, coef.A, coef.beta, coef.gamma;
                      C = cfg.dg.penalty_parameter, hooks = hooks, bcs = bcs)
    y = similar(u)
    apply_full!(y, op, u)
    r = y .- b
    split = Set(mesh.index[f[1]] for f in worldtube_boundary_faces(dom)
                if jc_nmortars(mesh, f) > 1)
    rin = 0.0; rout = 0.0; rsplit = 0.0
    for e in eachindex(mesh.elements)
        rng = mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])
        v = maximum(abs, view(r, rng))
        e in split && (rsplit = max(rsplit, v))
        inside_worldtube(dom, mesh.elements[e]) ? (rin = max(rin, v)) : (rout = max(rout, v))
    end
    return (rin = rin, rout = rout, rsplit = rsplit, nsplit = length(split))
end

@testset "DG jump conditions: worldtube + vtu hook table (dg-jump-conditions)" begin

@testset "exports and the seam contract" begin
    for s in (:vtu_junction_faces, :vtu_junction_hook, :jump_hooks)
        @test Base.isexported(ScalarSelfForce, s)
    end
    # the junction hook realizes main.tex:606-608 / spectre Equations.cpp:105-107 EXACTLY:
    # (n·F)_ext -> (n·F)_ext − imΩ·(Ψ_int + Ψ_ext)/2, and it is LINEARIZED (in the matrix)
    orb = CircularOrbit(0.5, 10.0)
    v = vtu_slicing(orb, 3.0, 26.0)
    jn = vtu_junctions(v, orb, 2)
    h = vtu_junction_hook(jn.vt)
    @test h isa LinearizedHook
    @test vtu_junction_hook(jn.vt.coeff) isa LinearizedHook
    int = FaceState(ComplexF64[1.0 + 0.5im, -0.25im], ComplexF64[0.3, -0.7 + 0.1im])
    ext = FaceState(ComplexF64[0.2 - 0.1im, 0.9], ComplexF64[-1.1, 0.05 + 0.2im])
    out = modify_boundary_data(h, int, ext)
    @test out.psi == ext.psi
    @test out.nflux == ext.nflux .- im * 2 * orb.Omega .* (int.psi .+ ext.psi) ./ 2
    @test modify_boundary_data_linearized(h, int, ext).nflux == out.nflux   # in the MATRIX
    # linear-homogeneous in the trace data (the LinearizedHook contract)
    z = FaceState(zeros(ComplexF64, 2), zeros(ComplexF64, 2))
    @test all(iszero, modify_boundary_data(h, z, z).nflux)
end

@testset "(b1) NODE GATE: 2:1 two-region prescribed jump through jump_hooks" begin
    a = jc_mini(6)
    @test a.nfaces == 3 && a.nmort == 4          # 2 fine-owned faces + 1 coarse (2 mortars)
    @test a.err ≤ 1e-10                          # GLOBAL_DAG row-40 gate
    r = jc_mini(6; reversed = true)
    n = jc_mini(6; hooks_on = false)
    @test r.err > 1e-3                           # the per-mortar ORDER is load-bearing
    @test n.err > 1e-3
    @test r.err > 1e6 * a.err
    println("jump-conditions (b1) mini-domain 2:1 jump: ", a.err,
            " | reversed ", r.err, " | no hooks ", n.err)
end

@testset "(b2) NODE GATE: manufactured vtu solution solves the modified-flux system" begin
    Ns = [4, 6, 8, 10]
    rs = [jc_vtu_case(N) for N in Ns]
    @test all(r -> r.njf == 4, rs)                # 2 junctions × 2 receivers, conforming
    res = [r.res for r in rs]
    ser = [r.serr for r in rs]
    @test issorted(res; rev = true) && issorted(ser; rev = true)
    @test res[3] ≤ 5e-3                           # N = 8 residual (calibrated 1.72e-3)
    @test ser[3] ≤ 1e-5                           # N = 8 solve  (calibrated 1.22e-6)
    @test jc_slope(Ns, res) ≤ -0.30                # SPECTRAL (calibrated −0.403 per N)
    @test jc_slope(Ns, ser) ≤ -0.40                # SPECTRAL (calibrated −0.552 per N)
    # FALSIFIERS at N = 8
    n8 = jc_vtu_case(8; hooks_on = false)
    s8 = jc_vtu_case(8; coeff = -vtu_junctions(vtu_slicing(CircularOrbit(0.0, 10.0),
                                                           JC_RV, JC_RU),
                                               CircularOrbit(0.0, 10.0), 2).vt.coeff)
    z8 = jc_vtu_case(8; coeff = 0.0 + 0.0im)
    @test n8.serr > 1e-2                          # calibrated 1.51e-1
    @test s8.serr > 1e-2                          # calibrated 2.91e-1 (sign flip)
    @test z8.serr == n8.serr                      # a zero coefficient IS no junction term
    @test n8.serr > 1e4 * ser[3]
    println("jump-conditions (b2) vtu N=4,6,8,10 residual ", res, " solve ", ser,
            " | slopes ", (jc_slope(Ns, res), jc_slope(Ns, ser)),
            " | N=8 falsifiers none ", n8.serr, " signflip ", s8.serr)
end

@testset "(b3) THE ROW-40 ORACLE: vtu in the matrix, worldtube in the fixed sources" begin
    cfg = jc_cfg(spin = 0.0, m = 2)
    orb = jc_orb(cfg)
    dom = block_domain(orb, cfg; orders = (4, 4))
    mesh = dom.mesh
    coef = domain_coefficients(dom, orb, cfg)
    wt = worldtube_data(dom, orb, cfg)
    jf = vtu_junction_faces(dom, orb, cfg, 2; mesh = mesh)
    empty_jd = Dict{Tuple{ElementId,Int,Int},Vector{FaceState{ComplexF64}}}()
    empty_jf = Dict{Tuple{ElementId,Int,Int},VtuJunction{Float64}}()
    bcs = Dict(f => JCTransmissive() for f in boundary_faces(mesh))
    mk(h) = MeshOperator(ComplexF64, mesh, coef.A, coef.beta, coef.gamma;
                         C = cfg.dg.penalty_parameter, hooks = h, bcs = bcs)
    op0 = mk(Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}())
    opv = mk(jump_hooks(mesh, empty_jd, jf))
    opw = mk(jump_hooks(mesh, worldtube_jump_data(wt), empty_jf))
    opb = mk(jump_hooks(dom, wt, orb, cfg, 2))
    A0, Av, Aw, Ab = build_matrix(op0), build_matrix(opv), build_matrix(opw), build_matrix(opb)
    # the VTU term IS the matrix change (it enters GMRES); the worldtube term is NOT
    @test Aw == A0                                       # bitwise
    @test Av != A0
    @test Ab == Av                                       # only the vtu part is in the matrix
    @test maximum(abs, Av .- A0) > 1e-6
    # …and the worldtube term is ALL of the fixed-source lift, the vtu term NONE of it
    @test all(iszero, operator_on_zero(op0))
    @test all(iszero, operator_on_zero(opv))
    zw = operator_on_zero(opw)
    @test maximum(abs, zw) > 1e-6
    @test operator_on_zero(opb) == zw
    println("jump-conditions (b3) row-40 oracle: |Av−A0|∞ ", maximum(abs, Av .- A0),
            " |Aw−A0|∞ ", maximum(abs, Aw .- A0), " |A_full(0)|∞ ", maximum(abs, zw))
end

@testset "(b4) the CRP production table: structure, coverage, finiteness" begin
    for L in (0, 1), chart in (:z2, :z), (a, m) in ((0.0, 2), (0.5, 1))
        cfg0 = jc_cfg(spin = a, m = m, L = L)
        cfg = chart === :z2 ? cfg0 : jc_full_z(cfg0)
        orb = jc_orb(cfg)
        dom = block_domain(orb, cfg; orders = (4, 4))
        mesh = dom.mesh
        wt = worldtube_data(dom, orb, cfg)
        jf = vtu_junction_faces(dom, orb, cfg, m)
        wf = Set(keys(worldtube_jump_data(wt)))
        tab = jump_hooks(dom, wt, orb, cfg, m)
        # the table covers EXACTLY the worldtube ∪ junction faces, one hook PER MORTAR
        @test Set(keys(tab)) == union(wf, Set(keys(jf)))
        @test all(f -> length(tab[f]) == jc_nmortars(mesh, f), keys(tab))
        # on the CRP lattice a worldtube face is NEVER a junction face: the tube spans
        # [5.2, 14.8] and the transitions sit at 2.8 / 24.4 (main.tex:545-548)
        @test isempty(intersect(wf, Set(keys(jf))))
        @test all(f -> dom.bounds_r[2] < mortar_points(mesh, f)[1] < dom.bounds_r[5],
                  [f for f in wf if f[2] == 1])
        # the junction faces sit exactly ON the transition radii
        @test Set(mortar_points(mesh, f)[1] for f in keys(jf)) ==
              Set([dom.bounds_r[2], dom.bounds_r[5]])
        @test all(j -> j.coeff == -im * m * orb.Omega, values(jf))
        # the per-mortar seam now serves the 2:1 worldtube faces (1 at L = 0, 2 at L = 1)
        nsplit = count(f -> jc_nmortars(mesh, f) > 1, wf)
        @test nsplit == 2^L
        @test all(f -> length(tab[f]) == jc_nmortars(mesh, f),
                  [f for f in wf if jc_nmortars(mesh, f) > 1])
        # …and the whole table builds an operator over the sealed coefficients + BCs that
        # applies FINITELY (the row-42 wiring)
        coef = domain_coefficients(dom, orb, cfg)
        op = MeshOperator(ComplexF64, mesh, coef.A, coef.beta, coef.gamma;
                          C = cfg.dg.penalty_parameter, hooks = tab,
                          bcs = boundary_conditions(dom, orb, cfg))
        u = ComplexF64[cis(0.7k) / (1 + k) for k in 1:mesh.ndofs]
        y = similar(u)
        apply!(y, op, u)
        @test all(isfinite, y)
        apply_full!(y, op, u)
        @test all(isfinite, y)
        @test all(isfinite, operator_on_zero(op))
    end
end

@testset "(b4) smooth-puncture control on the production mesh — SPLIT-OWNER INCLUDED" begin
    c = [jc_smooth_control(N) for N in (4, 6, 8, 10)]
    cr = [jc_smooth_control(N; mode = :reversed) for N in (4, 8)]
    cn = [jc_smooth_control(N; mode = :none) for N in (4, 8)]
    println("jump-conditions (b4) smooth control N=4,6,8,10: ", c,
            " | ratios ", [c[i] / c[i+1] for i in 1:3],
            " | reversed(4,8) ", cr, " | none(4,8) ", cn)
    @test c[1] ≤ 1e-3                              # calibrated 4.446e-4
    @test all(i -> c[i] / c[i+1] > 10, 1:3)        # SPECTRAL: calibrated 42.7 / 42.4 / 41.7
    @test c[4] ≤ 1e-8                              # calibrated 5.887e-9
    @test minimum(cr) / c[1] > 1e2                 # calibrated 330
    @test minimum(cn) / c[1] > 1e3                 # calibrated 2.97e4
    # the same control at L = 1 (TWO split worldtube faces) and at a = 0.5, m = 1
    c1 = jc_smooth_control(4; L = 1)
    @test c1 ≤ 1e-4                                # calibrated 2.613e-5
    @test jc_smooth_control(4; L = 1, mode = :none) / c1 > 1e3      # calibrated 6.0e5
    ca = jc_smooth_control(6; spin = 0.5, m = 1)
    @test ca ≤ 1e-4                                # calibrated 9.387e-6
    @test jc_smooth_control(6; spin = 0.5, m = 1, mode = :none) / ca > 1e3
    println("jump-conditions (b4) control L=1 ", c1, " | a=0.5,m=1 ", ca)
end

@testset "(b5) the REAL puncture: the split-owner element now falls with N" begin
    s4, s8 = jc_crp_residual(4), jc_crp_residual(8)
    d4, d8 = jc_crp_residual(4; mode = :dropsplit), jc_crp_residual(8; mode = :dropsplit)
    n4 = jc_crp_residual(4; mode = :none)
    println("jump-conditions (b5) N=4 in ", s4.rin, " out ", s4.rout, " split ", s4.rsplit,
            " | N=8 in ", s8.rin, " out ", s8.rout, " split ", s8.rsplit,
            " | dropped split ", (d4.rsplit, d8.rsplit), " | no hooks ", (n4.rin, n4.rout))
    @test s4.nsplit == 1
    # the row-38 envelope is reproduced BITWISE inside (the worldtube data are unchanged)
    @test s4.rin ≤ 1e-2 && s8.rin ≤ 5e-3           # calibrated 8.214e-3 / 2.930e-3
    @test s4.rin / s8.rin > 2
    # THE ROW-40 PRODUCT: the split-owner element is served and converges like the rest
    @test s4.rsplit ≤ 5e-3                         # calibrated 2.160e-3
    @test s8.rsplit ≤ 1e-3                         # calibrated 3.890e-4
    @test s4.rsplit / s8.rsplit > 3                # calibrated 5.55
    # …whereas dropping that one face's jump (all the sealed per-FACE seam could do) leaves a
    # residual an order larger that does NOT fall with N
    @test d4.rsplit > 1e-2                         # calibrated 3.374e-2
    @test d4.rsplit / s4.rsplit > 5
    @test d8.rsplit > d4.rsplit                     # it GROWS with N: an un-imposed jump, not
                                                   # truncation (calibrated 3.374e-2 → 7.012e-2)
    @test n4.rin > 1e-1 && n4.rout > 1e-1          # calibrated 0.588 / 0.645
end

@testset "(b6) BigFloat(256), the per-level callable, and m = 0" begin
    # --- m = 0: the junction coefficient vanishes identically (h drops out of the phase map)
    cfg0 = jc_cfg(m = 0)
    orb0 = jc_orb(cfg0)
    dom0 = block_domain(orb0, cfg0; orders = (4, 4))
    jf0 = vtu_junction_faces(dom0, orb0, cfg0, 0)
    @test !isempty(jf0)
    @test all(j -> iszero(j.coeff), values(jf0))
    coef0 = domain_coefficients(dom0, orb0, cfg0)
    bcs0 = Dict(f => JCTransmissive() for f in boundary_faces(dom0.mesh))
    mk0(h) = MeshOperator(ComplexF64, dom0.mesh, coef0.A, coef0.beta, coef0.gamma;
                          C = cfg0.dg.penalty_parameter, hooks = h, bcs = bcs0)
    empty_jd = Dict{Tuple{ElementId,Int,Int},Vector{FaceState{ComplexF64}}}()
    A0 = build_matrix(mk0(Dict{Tuple{ElementId,Int,Int},Vector{BoundaryDataHook}}()))
    Am0 = build_matrix(mk0(jump_hooks(dom0.mesh, empty_jd, jf0)))
    @test Am0 == A0                                # the m = 0 junction term is the identity

    # --- the per-LEVEL callable (the row-42 / multigrid build_op contract): a COARSER mesh of
    # the same lattice gets its own table; the WorldtubeData is mesh-bound and must be rebuilt
    cfg = jc_cfg(m = 2, L = 1)
    orb = jc_orb(cfg)
    dom = block_domain(orb, cfg; orders = (4, 4))
    co = coarsen_mesh(dom.mesh)
    @test co !== nothing
    cmesh = co[1]
    @test length(cmesh.elements) < length(dom.mesh.elements)
    jfc = vtu_junction_faces(dom, orb, cfg, 2; mesh = cmesh)
    @test !isempty(jfc)
    coefc = domain_coefficients(dom, orb, cfg)
    isin(el) = inside_worldtube(dom, el)
    valsc(r, w) = (psi = jc_smooth_P(r, w), f1 = coefc.A[1](r, w)*jc_smooth_Pr(r, w),
                   f2 = coefc.A[2](r, w)*jc_smooth_Pw(r, w))
    jdc = worldtube_jump_data(cmesh, isin, valsc)
    tabc = jump_hooks(cmesh, jdc, jfc)
    @test all(f -> length(tabc[f]) == jc_nmortars(cmesh, f), keys(tabc))
    opc = MeshOperator(ComplexF64, cmesh, coefc.A, coefc.beta, coefc.gamma;
                       C = cfg.dg.penalty_parameter, hooks = tabc,
                       bcs = boundary_conditions(dom, orb, cfg; mesh = cmesh))
    uc = ComplexF64[cis(0.3k) for k in 1:cmesh.ndofs]
    yc = similar(uc)
    apply!(yc, opc, uc)
    @test all(isfinite, yc)
    # a STALE table (built on the fine mesh) is REFUSED, not silently misapplied
    wt_fine = worldtube_data(dom, orb, cfg)
    @test_throws ArgumentError jump_hooks(cmesh, worldtube_jump_data(wt_fine), jfc)

    # THE PER-LEVEL RECIPE the docstring advertises, in one call each: worldtube_data takes the
    # same `mesh` kwarg, so a level rebuild needs no reconstruction of the region predicate
    wt_c = worldtube_data(dom, orb, cfg, 2; mesh = cmesh)
    @test length(wt_c.source) == cmesh.ndofs
    @test Set(wt_c.elements) == Set(filter(el -> inside_worldtube(dom, el), cmesh.elements))
    @test !isempty(wt_c.elements) && length(wt_c.elements) < length(wt_fine.elements)
    tab_l = jump_hooks(dom, wt_c, orb, cfg, 2; mesh = cmesh)
    @test all(f -> length(tab_l[f]) == jc_nmortars(cmesh, f), keys(tab_l))
    op_l = MeshOperator(ComplexF64, cmesh, coefc.A, coefc.beta, coefc.gamma;
                        C = cfg.dg.penalty_parameter, hooks = tab_l,
                        bcs = boundary_conditions(dom, orb, cfg; mesh = cmesh))
    apply_full!(yc, op_l, uc)
    @test all(isfinite, yc)
    @test maximum(abs, operator_on_zero(op_l)) > 1e-6      # the level's own worldtube lift
    # …and the fine-mesh datum is refused through the BlockDomain method too
    @test_throws ArgumentError jump_hooks(dom, wt_fine, orb, cfg, 2; mesh = cmesh)

    # M1: a P-REFINED mesh has the IDENTICAL face/mortar topology but more points per mortar,
    # so the mortar COUNT alone cannot catch a stale datum — the per-mortar POINT COUNT must.
    # (Node amr-driver p-refines every step.)
    pmesh = Mesh2D(dom.mesh.blocks, dom.mesh.elements, (6, 6);
                   connections = dom.mesh.connections)
    jf_p = vtu_junction_faces(dom, orb, cfg, 2; mesh = pmesh)
    @test Set(keys(jf_p)) == Set(keys(vtu_junction_faces(dom, orb, cfg, 2)))   # same topology
    jd_fine = worldtube_jump_data(wt_fine)
    @test all(f -> length(jd_fine[f]) == jc_nmortars(pmesh, f), keys(jd_fine))  # counts MATCH
    ep = @test_throws ArgumentError jump_hooks(pmesh, jd_fine, jf_p)
    @test occursin("points", ep.value.msg)
    # the same data on their OWN order do pass
    @test jump_hooks(dom.mesh, jd_fine, vtu_junction_faces(dom, orb, cfg, 2)) isa AbstractDict

    # M1: a WorldtubeData from a DIFFERENT BlockDomain (identical block ids) is refused
    cfg_b = jc_cfg(spin = 0.5, m = 2, L = 1)
    orb_b = jc_orb(cfg_b)
    dom_b = block_domain(orb_b, cfg_b; orders = (4, 4))
    @test dom_b.block_ids == dom.block_ids                 # same ids, different physics
    ed = @test_throws ArgumentError jump_hooks(dom, worldtube_data(dom_b, orb_b, cfg_b),
                                               orb, cfg, 2; mesh = dom_b.mesh)
    @test occursin("DIFFERENT BlockDomain", ed.value.msg)

    # --- BigFloat(256): the junction datum and the hook inherit the ambient precision
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        orbB = CircularOrbit(big(0.0), big(10.0))
        cfgB = crp_paper_config(; spin = big(0.0), orbital_radius = big(10.0), m_mode = 2,
                                inner_radius = orbB.rp, outer_radius = big(1.0e4),
                                worldtube_radius = big(4.8), worldtube_dz2 = big(0.25))
        @test cfgB isa RunConfig{BigFloat}
        domB = block_domain(orbB, cfgB; orders = (4, 4))
        jfB = vtu_junction_faces(domB, orbB, cfgB, 2)
        @test !isempty(jfB)
        jB = first(values(jfB))
        @test jB isa VtuJunction{BigFloat}
        @test precision(real(jB.coeff)) == 256
        hB = vtu_junction_hook(jB)
        intB = FaceState(Complex{BigFloat}[big(1.0), big(0.5)],
                         Complex{BigFloat}[big(0.25), big(-0.75)])
        extB = FaceState(Complex{BigFloat}[big(-0.3), big(0.9)],
                         Complex{BigFloat}[big(1.5), big(0.1)])
        outB = modify_boundary_data_linearized(hB, intB, extB)
        @test eltype(outB.nflux) == Complex{BigFloat}
        @test precision(real(outB.nflux[1])) == 256
        @test outB.nflux == extB.nflux .+ jB.coeff .* (intB.psi .+ extB.psi) ./ 2
        # and the whole table constructs in BigFloat
        valsB(r, w) = (psi = Complex{BigFloat}(jc_smooth_P(Float64(r), Float64(w))),
                       f1 = zero(Complex{BigFloat}), f2 = zero(Complex{BigFloat}))
        jdB = worldtube_jump_data(domB.mesh, el -> inside_worldtube(domB, el), valsB)
        tabB = jump_hooks(domB.mesh, jdB, jfB)
        @test all(f -> length(tabB[f]) == jc_nmortars(domB.mesh, f), keys(tabB))
        @test Set(keys(tabB)) == union(Set(keys(jdB)), Set(keys(jfB)))
    finally
        setprecision(BigFloat, old)
    end
end

end # testset dg-jump-conditions

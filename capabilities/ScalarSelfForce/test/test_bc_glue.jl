# test_bc_glue.jl — ESSF node `dg-boundary-conditions` (schedule W8-D, GLOBAL_DAG row 41),
# written BEFORE the implementation (TDD). The node delivers `Solve/bc_glue.jl`: the binding
# of the CRP's PHYSICAL boundary conditions to the EllipticCore ghost seam, i.e. the
# `MeshOperator` `bcs` table of the block domain — `boundary_conditions(dom, orbit, config)`
# covering exactly `boundary_faces(mesh)` for the domain's own mesh AND for every coarser
# mesh a multigrid hierarchy hands its `build_op`.
#
# ORACLE (the DAG row's gate column) — the READ-ONLY mirrors, transcribed by hand below:
#   * `main.tex:396-401` — "The vanishing of F^r_m at the horizon means that we do not have to
#     explicitly impose boundary conditions there. Similarly, factoring out the sin^{|m|}θ
#     behavior causes F^z_m to vanish at the poles such that the z boundaries do not require
#     explicit boundary treatment either." (regularity conditions, no explicit BCs).
#   * `main.tex:614-616` — "On external boundaries we impose boundary conditions through the
#     DG numerical flux by a choice of Ψ_ext and (n_i F^i)_ext. Note that on boundaries where
#     the perpendicular flux vanishes, n_i 𝔄^i_m = 0, neither of the two boundary terms
#     contributes, as Ψ_ext is ignored and (n_i F^i)_ext = 0. This occurs on the horizon
#     boundary and on angular boundaries. On the last remaining boundary at finite outer
#     radius, R = 1e4, we impose the second-order Bayliss-Turkel condition."
#   * `anc/ScalarSelfForce.yaml:113` `BoundaryConditions[r].Lower: None`, `:127`
#     `BoundaryConditions[z²]: None`, `:114-125` the Jinja2 `{% if m_mode == 0 %} None {% else
#     %} Sommerfeld` branch on `BoundaryConditions[r].Upper`.
#   * `ref-code/spectre .../Scalar/BoundaryConditions/None.{hpp,cpp}` — "Applies no boundary
#     condition at all. Used to impose nothing but regularity at the horizon in
#     horizon-penetrating coordinates or at angular boundaries"; `apply` and
#     `apply_linearized` are both EMPTY, and spectre's elliptic DG then forms the exterior as
#     `u_ext = 2 u_bc − u_int`, `(n.F)_ext = −2 (n.F)_bc + (n.F)_int`
#     (`Elliptic/DiscontinuousGalerkin/DgOperator.hpp:475-508`). With the empty `apply`
#     (u_bc = u_int, (n.F)_bc = (n.F)_int) that is EXACTLY `(ψ_int, −φ_int)` — the
#     transmissive exterior of the sealed dg_operator.jl convention N7.
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here):
#   * `bayliss-turkel-bc` — `BaylissTurkelGhost`, `bt_robin_coefficient`,
#     `outer_boundary_ghost` (incl. the m = 0 ⇒ `nothing` yaml branch). This node NEVER
#     recomputes a Robin coefficient; it only decides WHICH face gets one and at WHICH points.
#   * `domain-blocks-map` — `block_domain` (which ENFORCES bounds_r[1] == r₊, the precondition
#     the `:horizon` no-data contract rests on), `boundary_face_classes`/`boundary_class`,
#     `block_lattice_index`, `domain_coefficients`.
#   * `boundary-data-hooks` / `matrix-free-operator` / `multigrid-hierarchy` (EllipticCore) —
#     `GhostBoundary`, `ghost_state`, `DirichletGhost`, `NeumannGhost`, `MeshOperator`,
#     `boundary_faces`, `apply!`, `direct_solve`, `MultigridHierarchy`, `MultigridV`, `gmres`.
#
# THIS node's gates (constants calibrated PRE-RED in the evidence dir,
# calibrate_bc_glue{,2,3,4,5}.jl — never tuned post hoc):
#  (a) COVERAGE — `boundary_conditions` key set == `boundary_faces(mesh)` for the 5×2 and 5×3
#      layouts at L ∈ {0, 1} (18/20/36/40 faces), the class split matches the sealed
#      `boundary_face_classes`, only `:outer` carries a `BaylissTurkelGhost`, and a
#      `MeshOperator` over `domain_coefficients` constructs and applies.
#  (b) NO-OP EQUIVALENCE at horizon/equator/poles — `Aface` is EXACTLY 0 there (measured 0.0,
#      against 9.998e-1 on the outer face), so through the sealed `PenaltyFlux` the None ghost
#      gives (Ψ*, (nF)*) == (ψ_int, φ_int) BITWISE and both lifted face terms are exactly the
#      zero vector; at the operator level None == Dirichlet-with-ANY-data (bitwise, 0.0
#      measured — the face is deaf to Dirichlet data because the penalty and symmetrization
#      terms both carry the vanishing 𝔄), while the FALSIFIER `NeumannGhost(0)` — the only
#      sealed ghost that touches the flux — leaks 7.7e-1 (horizon) / 5.5e0 (equator) / 5.0e0
#      (pole) of the local residual scale per face, and 1.3e-2 of ‖A·u‖_∞ at the operator
#      level (3.5e-3 … 1.3e-2 over the 5×2/5×3, L ∈ {0, 1}, m ∈ {2, 5} calibration sweep —
#      the per-face and the global normalization differ by which element the maximum sits
#      on; calibrate_bc_glue5.jl). Float64 and BigFloat(256).
#  (c) REDUCED-ODE CONTAMINATION — a = 0 makes the operator z-independent and a z-uniform
#      field reduces it to the radial ODE (measured per-element z-uniformity 1.2e-13); with a
#      smooth compactly-supported source in the worldtube column, the solution at the
#      extraction column changes by 8.3e-8 (m = 2) / 1.0e-7 (m = 1) relative when R = 1e4 is
#      replaced by 1e5 — against the 1e-6 node gate. That production-resolution number is
#      DISCRETIZATION-limited: it falls to 2.0e-9 at P = 6 and 9.1e-11 at P = 8, i.e. below
#      the BT reflection scale (m+1)/(2mΩ²R³) = 7.5e-10 of the sealed row-32 gate. FALSIFIER:
#      a `DirichletGhost(0)` outer face gives 1.7e-6 at P = 4 and GROWS with resolution
#      (3.1e-6 at P = 6, 5.7e-6 at P = 8) — a genuine BC error, not a truncation error.
#      The solve runs through the sealed stack `gmres ∘ MultigridV(MultigridHierarchy(op,
#      build_op))` (11 iterations, agreeing with `direct_solve` to 2.5e-14), which is what
#      exercises `boundary_conditions` on the COARSE level meshes.
#  (d) m = 0 — `outer_radial_bc === nothing` ⇒ every face gets the None ghost; the operator
#      builds, applies and solves (cond 2.7e5, residual 7.9e-13) and its answer at the
#      extraction column is insensitive to R at 1.1e-10 (the static m = 0 mode does not
#      radiate, so "no outer condition" costs nothing).
#  (e) TYPE / CONFIG FLOW — BigFloat(256) construction inherits precision into the Robin
#      coefficients; the robin vector is BITWISE the sealed `bt_robin_coefficient` broadcast
#      over the operator's OWN face-point coefficients; the m-dependence and the config's
#      Sommerfeld order reach the ghost through the `RunConfig` alone.
#
# NO SCOPE CREEP: no worldtube/vtu jump hooks (`hooks` stays empty — nodes 38/40), no solver
# driver, no AMR, no effective source. The operator assembled here is the pure block-lattice
# operator of the sealed `domain_coefficients` bundle.
using Test
using LinearAlgebra
using ScalarSelfForce
using ScalarSelfForce.EllipticCore
using ScalarSelfForce.EllipticCore: element_matrices, element_coefficients

# ---------------------------------------------------------------------------------------
# helpers (test-side only)
# ---------------------------------------------------------------------------------------
bcg_orbit(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)

# a RunConfig twin with ImposeEquatorialSymmetry = false (the 5×3 full-z lattice);
# `crp_paper_config` hard-codes the yaml's True and exposes no keyword (domain-blocks-map
# caveat 6), so the layout is reached by a field-wise reconstruction of the sealed structs.
function bcg_full_z(cfg::RunConfig{T}) where {T}
    OC = typeof(cfg.orbit)
    oc = OC((f === :impose_equatorial_symmetry ? false : getfield(cfg.orbit, f)
             for f in fieldnames(OC))...)
    return RunConfig{T}((f === :orbit ? oc : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

# the face-point sampling `element_coefficients` uses (conventions N1/N2): the face's
# normal-dimension block-map endpoint × the transverse volume nodes
function bcg_face_points(mesh, el::ElementId, d::Int, sidx::Int)
    e = mesh.index[el]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    RT = typeof(mesh).parameters[1]
    em = element_matrices(RT, blk, el, mesh.orders[e])
    lo, hi = logical_interval(RT, el, d)
    return RT(map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)), em.x[3 - d]
end

# per-element matrices + coefficients of one element of a mesh
function bcg_element(::Type{T}, mesh, el::ElementId, fcs) where {T}
    e = mesh.index[el]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    em = element_matrices(T, blk, el, mesh.orders[e])
    return em, element_coefficients(em, blk, el, fcs.A, fcs.beta, fcs.gamma)
end

# the residual contribution of ONE face (the two lifted terms of eq:dgres_strong_matrix,
# sealed convention N5) — transcribed here so the no-op claim is checked against the
# ASSEMBLY, not against another ghost
function bcg_face_contribution(em, c, st, d::Int, sidx::Int, ext::FaceState, C)
    fm = em.faces[d][sidx]
    Af = c.Aface[d][sidx]
    int = st.faces[d][sidx]
    s = sidx == 1 ? -1 : 1
    ψstar, φstar = numerical_flux(PenaltyFlux(C), fm, Af, int, ext)
    return fm.MLDT * ((s .* Af) .* (ψstar .- int.psi)) .- fm.ML * (φstar .- int.nflux)
end

# massive DG right-hand side b = M·S of a nodal source callable S(x1, x2)
function bcg_massive_rhs(::Type{T}, mesh, S) where {T}
    b = zeros(T, mesh.ndofs)
    for (e, el) in enumerate(mesh.elements)
        blk = mesh.blocks[mesh.blockindex[el.block]]
        em = element_matrices(T, blk, el, mesh.orders[e])
        s = T[S(x1, x2) for x2 in em.x[2] for x1 in em.x[1]]
        b[mesh.offsets[e]+1 : mesh.offsets[e]+prod(mesh.orders[e])] = em.M * s
    end
    return b
end

# smooth compactly-supported radial bump, support (6, 13) strictly inside the worldtube
# column (5.2, 14.8) — identical in the R and 10R domains, zero in the wavezone column
bcg_bump(r, w) = (6.0 < r < 13.0) ? ComplexF64(((r - 6.0) * (13.0 - r) / 12.25)^4) :
                                    ComplexF64(0)

bcg_meshop(mesh, cfg, fcs, bcs) =
    MeshOperator(ComplexF64, mesh, fcs.A, fcs.beta, fcs.gamma;
                 C = cfg.dg.penalty_parameter, bcs = bcs)

# solve the manufactured radial problem on the CRP lattice of `cfg`, returning the solution
# and the dof range of the extraction (worldtube) column
function bcg_solve(cfg; outer = :production, stack = :direct)
    orb = bcg_orbit(cfg)
    dom = block_domain(orb, cfg)
    fcs = domain_coefficients(dom, orb, cfg)
    bcs = boundary_conditions(dom, orb, cfg)
    if outer !== :production
        for (f, _) in collect(bcs)
            boundary_class(dom, f) === :outer || continue
            nf = length(bcg_face_points(dom.mesh, f[1], f[2], f[3])[2])
            bcs[f] = outer === :dirichlet0 ? DirichletGhost(zeros(ComplexF64, nf)) :
                     outer === :none ? NoneGhost() :
                     throw(ArgumentError("unknown outer variant $outer"))
        end
    end
    op = bcg_meshop(dom.mesh, cfg, fcs, bcs)
    b = bcg_massive_rhs(ComplexF64, dom.mesh, bcg_bump)
    u = if stack === :direct
        direct_solve(op, b)
    else
        hier = MultigridHierarchy(op, m -> bcg_meshop(m, cfg, fcs,
                                          boundary_conditions(dom, orb, cfg; mesh = m)))
        res = gmres(op, b; Minv! = MultigridV(hier), abstol = 0, reltol = 1e-12,
                    maxiter = 200)
        @test res.converged
        res.x
    end
    idx = Int[]
    for (e, el) in enumerate(dom.mesh.elements)
        i, _ = block_lattice_index(dom, el.block)
        i == 2 && append!(idx,
                          dom.mesh.offsets[e]+1 : dom.mesh.offsets[e]+prod(dom.mesh.orders[e]))
    end
    return (dom = dom, op = op, b = b, u = u, idx = idx)
end

# relative change of the extraction-column solution when R → 10R (the contamination oracle)
function bcg_contamination(; m, P = 4, L = 0, outer = :production, stack = :direct)
    base(R) = crp_paper_config(spin = 0.0, orbital_radius = 10.0, m_mode = m,
                               outer_radius = R, P = P, L = L)
    a = bcg_solve(base(1.0e4); outer = outer, stack = stack)
    b = bcg_solve(base(1.0e5); outer = outer, stack = stack)
    # columns 0..3 are the SAME elements in both domains (only the wavezone map changes),
    # so the two solution vectors are directly comparable on the extraction dofs
    @test a.dom.mesh.elements == b.dom.mesh.elements
    return maximum(abs, a.u[a.idx] .- b.u[b.idx]) / maximum(abs, a.u[a.idx]), a
end

const BCG_CFG = crp_paper_config()                       # production 5×2, m = 2, a = 0

@testset "Physical BCs → EllipticCore ghost seam (dg-boundary-conditions)" begin

    # -----------------------------------------------------------------------------------
    @testset "(0) API surface" begin
        for name in (:NoneGhost, :boundary_conditions)
            @test name in names(ScalarSelfForce)
        end
        @test NoneGhost() isa GhostBoundary
        @test NoneGhost() == NoneGhost()                 # stateless singleton semantics
    end

    # -----------------------------------------------------------------------------------
    # (a) coverage: exactly the boundary faces, the sealed classes, and a working operator
    # -----------------------------------------------------------------------------------
    @testset "(a) coverage of boundary_faces + MeshOperator wiring" begin
        for L in (0, 1), sym in (true, false)
            cfg0 = crp_paper_config(L = L)
            cfg = sym ? cfg0 : bcg_full_z(cfg0)
            orb = bcg_orbit(cfg)
            dom = block_domain(orb, cfg)
            bcs = boundary_conditions(dom, orb, cfg)
            bf = boundary_faces(dom.mesh)
            @test bcs isa AbstractDict
            @test Set(keys(bcs)) == Set(bf)              # EXACTLY, no more and no fewer
            @test length(bcs) == length(bf)
            @test valtype(bcs) <: GhostBoundary
            # the class split is the sealed table's, and only :outer carries physical data
            classes = boundary_face_classes(dom)
            @test Set(keys(classes)) == Set(bf)
            for f in bf
                if classes[f] === :outer
                    @test bcs[f] isa BaylissTurkelGhost  # m = 2 ≠ 0 ⇒ the BT ghost
                else
                    @test classes[f] in (:horizon, :equator, :pole)
                    @test bcs[f] isa NoneGhost
                end
            end
            # the wavezone column carries (L+1, L) levels, so 2^L outer faces per angular row
            nouter = count(f -> classes[f] === :outer, bf)
            @test nouter == (sym ? 2 : 3) * 2^L
            # the operator the CRP solves with: ONE coefficient bundle, no hooks
            fcs = domain_coefficients(dom, orb, cfg)
            op = bcg_meshop(dom.mesh, cfg, fcs, bcs)
            u = ComplexF64[cis(0.3k) / (1 + k) for k in 1:dom.mesh.ndofs]
            y = similar(u)
            @test apply!(y, op, u) === y
            @test all(isfinite, y)
            @test !iszero(y)
        end
    end

    @testset "(a2) the per-level meshes a multigrid hierarchy builds" begin
        cfg = BCG_CFG
        orb = bcg_orbit(cfg)
        dom = block_domain(orb, cfg)
        fcs = domain_coefficients(dom, orb, cfg)
        op = bcg_meshop(dom.mesh, cfg, fcs, boundary_conditions(dom, orb, cfg))
        nlevels = 0
        hier = MultigridHierarchy(op, m -> begin
            nlevels += 1
            b = boundary_conditions(dom, orb, cfg; mesh = m)
            # the contract the hierarchy demands: keyed by THIS mesh's own boundary faces
            @test Set(keys(b)) == Set(boundary_faces(m))
            @test any(g -> g isa BaylissTurkelGhost, values(b))
            @test any(g -> g isa NoneGhost, values(b))
            bcg_meshop(m, cfg, fcs, b)
        end)
        @test length(hier.levels) >= 2
        @test nlevels == length(hier.levels) - 1
        # a coarse mesh whose blocks are foreign to `dom` must be REFUSED, not guessed at
        foreign = Mesh2D(Block(99, LinearMap1D(1.0, 2.0), LinearMap1D(0.0, 1.0)),
                         [ElementId(99, SegmentId(0, 0), SegmentId(0, 0))], (4, 4))
        @test_throws ArgumentError boundary_conditions(dom, orb, cfg; mesh = foreign)
    end

    # -----------------------------------------------------------------------------------
    # (b) the no-data faces: the flux coefficient vanishes and the seam is a no-op
    # -----------------------------------------------------------------------------------
    @testset "(b1) 𝔄 vanishes EXACTLY on horizon/equator/pole faces" begin
        for sym in (true, false)
            cfg = sym ? BCG_CFG : bcg_full_z(BCG_CFG)
            orb = bcg_orbit(cfg)
            dom = block_domain(orb, cfg)
            fcs = domain_coefficients(dom, orb, cfg)
            outer_scale = 0.0
            for f in boundary_faces(dom.mesh)
                el, d, sidx = f
                _, c = bcg_element(ComplexF64, dom.mesh, el, fcs)
                Af = c.Aface[d][sidx]
                if boundary_class(dom, f) === :outer
                    outer_scale = max(outer_scale, maximum(abs, Af))
                else
                    @test all(iszero, Af)                # EXACTLY zero (main.tex:614)
                    # the sealed callables say the same thing at the same points
                    xf, xt = bcg_face_points(dom.mesh, el, d, sidx)
                    vals = d == 1 ? [fcs.A[1](xf, t) for t in xt] :
                                    [fcs.A[2](t, xf) for t in xt]
                    @test all(iszero, vals)
                    @test vals == Af
                end
            end
            @test outer_scale > 0.9                      # calibrated 9.998e-1
        end
    end

    @testset "(b2) the None ghost makes the numerical flux the INTERIOR flux, bitwise" begin
        for TT in (ComplexF64, Complex{BigFloat})
            setprecision(BigFloat, 256) do
                cfg = BCG_CFG
                orb = bcg_orbit(cfg)
                dom = block_domain(orb, cfg)
                fcs = domain_coefficients(dom, orb, cfg)
                C = cfg.dg.penalty_parameter
                nchecked = 0
                for f in boundary_faces(dom.mesh)
                    boundary_class(dom, f) === :outer && continue
                    el, d, sidx = f
                    em, c = bcg_element(TT, dom.mesh, el, fcs)
                    n = prod(dom.mesh.orders[dom.mesh.index[el]])
                    u = TT[cis(TT === ComplexF64 ? 0.7k : big(0.7) * k) / (1 + k) for k in 1:n]
                    st = element_state(em, c, u)
                    int = st.faces[d][sidx]
                    ext = ghost_state(NoneGhost(), int)
                    # spectre None.cpp through DgOperator.hpp:475-508 == convention N7
                    @test ext.psi == int.psi
                    @test ext.nflux == -int.nflux
                    # linear-homogeneous ⇒ the linearized ghost IS the full one, bitwise
                    lin = ghost_state_linearized(NoneGhost(), int)
                    @test lin.psi == ext.psi && lin.nflux == ext.nflux
                    # (Ψ*, (nF)*) == (ψ_int, φ_int) bitwise ⇒ both lifted terms vanish
                    ψs, φs = numerical_flux(PenaltyFlux(C), em.faces[d][sidx],
                                            c.Aface[d][sidx], int, ext)
                    @test ψs == int.psi
                    @test φs == int.nflux
                    @test all(iszero, bcg_face_contribution(em, c, st, d, sidx, ext, C))
                    nchecked += 1
                end
                @test nchecked == 16                     # 2 horizon + 7 equator + 7 pole
            end
        end
    end

    @testset "(b3) operator-level equivalence and the flux-touching falsifier" begin
        cfg = BCG_CFG
        orb = bcg_orbit(cfg)
        dom = block_domain(orb, cfg)
        fcs = domain_coefficients(dom, orb, cfg)
        base = boundary_conditions(dom, orb, cfg)
        nodata = [f for f in boundary_faces(dom.mesh) if boundary_class(dom, f) !== :outer]
        swap(g) = begin
            b = copy(base)
            for f in nodata
                b[f] = g(length(bcg_face_points(dom.mesh, f[1], f[2], f[3])[2]))
            end
            b
        end
        op0 = bcg_meshop(dom.mesh, cfg, fcs, base)
        u = ComplexF64[cis(0.9k) / (1 + sqrt(k)) for k in 1:dom.mesh.ndofs]
        y0 = similar(u); apply!(y0, op0, u)
        scale = maximum(abs, y0)
        # Dirichlet data of ANY value is INVISIBLE at these faces (main.tex:615 "Ψ_ext is
        # ignored") — both the symmetrization and the penalty carry the vanishing 𝔄
        for gd in (nf -> DirichletGhost(zeros(ComplexF64, nf)),
                   nf -> DirichletGhost(ComplexF64[cis(0.4k) for k in 1:nf]))
            opd = bcg_meshop(dom.mesh, cfg, fcs, swap(gd))
            yd = similar(u); apply!(yd, opd, u)
            @test maximum(abs, yd .- y0) <= 1e-14 * scale   # NODE GATE
            @test yd == y0                                  # in fact BITWISE (calibrated 0.0)
        end
        # FALSIFIER: the only sealed ghost that touches the FLUX changes the operator far
        # above the gate — the no-op is a property of the None ghost, not of the faces alone
        opn = bcg_meshop(dom.mesh, cfg, fcs, swap(nf -> NeumannGhost(zeros(ComplexF64, nf))))
        yn = similar(u); apply!(yn, opn, u)
        leak = maximum(abs, yn .- y0) / scale
        @info "bc-glue gate (b3) NeumannGhost(0) leakage $(leak) (Dirichlet: bitwise 0)"
        @test leak > 1e-3               # calibrated 1.28e-2 here, 3.5e-3 … 1.3e-2 over the
                                        # 5×2/5×3, L ∈ {0,1}, m ∈ {2,5} sweep — 11 orders
                                        # above the 1e-14 no-op gate
    end

    # -----------------------------------------------------------------------------------
    # (c) the reduced-ODE contamination oracle at the outer boundary
    # -----------------------------------------------------------------------------------
    @testset "(c) radial reduced-ODE contamination at R = 1e4 vs 1e5" begin
        cont2, run2 = bcg_contamination(m = 2)
        cont1, _ = bcg_contamination(m = 1)
        @info "bc-glue gate (c) contamination (P=4, L=0): m=1 $(cont1)  m=2 $(cont2)"
        @test cont2 < 1e-6                               # NODE GATE (calibrated 8.3e-8)
        @test cont1 < 1e-6                               # calibrated 1.0e-7
        # the a = 0 operator is z-independent, so the solution IS the radial ODE's
        dev = 0.0
        for (e, el) in enumerate(run2.dom.mesh.elements)
            nx, ny = run2.dom.mesh.orders[e]
            U = reshape(run2.u[run2.dom.mesh.offsets[e]+1 : run2.dom.mesh.offsets[e]+nx*ny],
                        nx, ny)
            dev = max(dev, maximum(abs, U .- U[:, 1]) / maximum(abs, U))
        end
        @info "bc-glue gate (c) per-element z-uniformity $(dev)"
        @test dev < 1e-11                                # calibrated 1.2e-13
        # the production number is TRUNCATION-limited: raising p drives it below the
        # BT reflection scale (m+1)/(2mΩ²R³) of the sealed row-32 gate
        cont2p6, _ = bcg_contamination(m = 2, P = 6)
        refl = 3 / (2 * 2 * bcg_orbit(BCG_CFG).Omega^2 * 1.0e4^3)
        @info "bc-glue gate (c) P=6 contamination $(cont2p6); BT reflection scale $(refl)"
        @test cont2p6 < cont2 / 10                       # calibrated 2.0e-9 (41× better)
        # FALSIFIER: a wrong outer condition breaks the gate and gets WORSE with resolution
        bad4, _ = bcg_contamination(m = 2, outer = :dirichlet0)
        bad6, _ = bcg_contamination(m = 2, P = 6, outer = :dirichlet0)
        @info "bc-glue gate (c) falsifier Dirichlet(0) outer: P=4 $(bad4)  P=6 $(bad6)"
        @test bad4 > 1e-6                                # calibrated 1.7e-6
        @test bad6 > bad4                                # calibrated 3.1e-6 — a BC error
    end

    @testset "(c2) the sealed gmres ∘ MultigridV stack drives the same answer" begin
        cont, run = bcg_contamination(m = 2, stack = :mg)
        @info "bc-glue gate (c2) contamination through gmres ∘ MultigridV $(cont)"
        @test cont < 1e-6
        # and the preconditioned Krylov solution IS the direct one
        ud = direct_solve(run.op, run.b)
        @test norm(run.u .- ud) <= 1e-12 * norm(ud)      # calibrated 2.5e-14
    end

    # -----------------------------------------------------------------------------------
    # (d) the m = 0 branch: the yaml's `None` outer boundary
    # -----------------------------------------------------------------------------------
    @testset "(d) m = 0 — no outer boundary condition" begin
        cfg = crp_paper_config(spin = 0.0, m_mode = 0)
        @test cfg.domain.outer_radial_bc === nothing     # yaml:115-125 Jinja2 branch
        orb = bcg_orbit(cfg)
        dom = block_domain(orb, cfg)
        bcs = boundary_conditions(dom, orb, cfg)
        @test Set(keys(bcs)) == Set(boundary_faces(dom.mesh))
        @test all(g -> g isa NoneGhost, values(bcs))     # incl. the OUTER faces
        @test !any(g -> g isa BaylissTurkelGhost, values(bcs))
        # the operator builds, applies and solves — and its answer at the extraction column
        # is insensitive to the outer radius: the m = 0 mode is static and does not radiate
        cont, run = bcg_contamination(m = 0)
        @info "bc-glue gate (d) m = 0 contamination $(cont)"
        @test cont < 1e-6                                # calibrated 1.1e-10
        y = similar(run.u); apply!(y, run.op, run.u)
        @test norm(y .- run.b) <= 1e-10 * norm(run.b)    # calibrated 7.9e-13
        # the m = 0 lattice is identical to the m ≠ 0 one — only the outer ghost differs
        b2 = boundary_conditions(dom, orb, crp_paper_config(spin = 0.0, m_mode = 2))
        @test Set(keys(b2)) == Set(keys(bcs))
        @test count(g -> g isa BaylissTurkelGhost, values(b2)) == 2
    end

    # -----------------------------------------------------------------------------------
    # (e) type genericity and RunConfig flow
    # -----------------------------------------------------------------------------------
    @testset "(e) BigFloat(256), the Robin data, and the config knobs" begin
        cfg = BCG_CFG
        orb = bcg_orbit(cfg)
        dom = block_domain(orb, cfg)
        fcs = domain_coefficients(dom, orb, cfg)
        bcs = boundary_conditions(dom, orb, cfg)
        nouter = 0
        for f in boundary_faces(dom.mesh)
            boundary_class(dom, f) === :outer || continue
            el, d, sidx = f
            @test d == 1 && sidx == 2
            xf, xt = bcg_face_points(dom.mesh, el, d, sidx)
            want = bt_robin_coefficient.(Ref(cfg.domain.outer_radial_bc),
                                         [fcs.A[1](xf, t) for t in xt],
                                         [fcs.beta(xf, t) for t in xt],
                                         [fcs.gamma[1](xf, t) for t in xt])
            @test bcs[f].robin == want                   # BITWISE the sealed broadcast
            # and the α it used is the operator's OWN face coefficient
            _, c = bcg_element(ComplexF64, dom.mesh, el, fcs)
            @test [fcs.A[1](xf, t) for t in xt] == c.Aface[1][2]
            nouter += 1
        end
        @test nouter == 2
        # the m-dependence flows through the RunConfig alone
        r(m) = begin
            c = crp_paper_config(spin = 0.0, m_mode = m)
            d = block_domain(bcg_orbit(c), c)
            b = boundary_conditions(d, bcg_orbit(c), c)
            first(g.robin for g in values(b) if g isa BaylissTurkelGhost)
        end
        @test r(1) != r(2) && r(2) != r(3)
        # BigFloat(256) construction: precision is inherited end to end
        setprecision(BigFloat, 256) do
            cb = crp_paper_config(spin = big(0.0), orbital_radius = big(10.0), m_mode = 2,
                                  inner_radius = big(2.0), outer_radius = big(1.0e4),
                                  worldtube_radius = big(4.8), worldtube_dz2 = big(0.25))
            @test cb isa RunConfig{BigFloat}
            ob = bcg_orbit(cb)
            db = block_domain(ob, cb)
            bb = boundary_conditions(db, ob, cb)
            @test Set(keys(bb)) == Set(boundary_faces(db.mesh))
            bt = first(g for g in values(bb) if g isa BaylissTurkelGhost)
            @test eltype(bt.robin) === Complex{BigFloat}
            @test precision(real(bt.robin[1])) == 256
            @test all(g -> g isa NoneGhost,
                      [bb[f] for f in boundary_faces(db.mesh)
                       if boundary_class(db, f) !== :outer])
            # the None ghost is precision-preserving as well
            int = FaceState(Complex{BigFloat}[big(1.0) + 0im], Complex{BigFloat}[big(2.0) + 0im])
            g = ghost_state(NoneGhost(), int)
            @test eltype(g.psi) === Complex{BigFloat}
            @test g.nflux == -int.nflux
        end
        # the config's Sommerfeld ORDER reaches the ghost: order 1 under hyperboloidal
        # slicing is spectre's plain-outflow robin = 0 (Sommerfeld.cpp:43)
        SB = typeof(cfg.domain.outer_radial_bc)
        bc1 = SB((f === :order ? 1 : getfield(cfg.domain.outer_radial_bc, f)
                  for f in fieldnames(SB))...)
        DC = typeof(cfg.domain)
        dom1 = DC((f === :outer_radial_bc ? bc1 : getfield(cfg.domain, f)
                   for f in fieldnames(DC))...)
        cfg1 = RunConfig{Float64}((f === :domain ? dom1 : getfield(cfg, f)
                                   for f in fieldnames(RunConfig{Float64}))...)
        b1 = boundary_conditions(dom, orb, cfg1)
        for f in boundary_faces(dom.mesh)
            boundary_class(dom, f) === :outer || continue
            @test b1[f] isa BaylissTurkelGhost
            @test all(iszero, b1[f].robin)
        end
    end
end

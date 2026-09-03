# test_assemble.jl — operator assembly (ESSF node `fig-matrix`, GLOBAL_DAG row 49, commit 1;
# the assembly half of the row-42 `solver-stack-ssf` contract).
#
# WHAT THE NODE DELIVERS. `Solve/assemble.jl`: the ONE composition of the sealed CRP pieces
# into the production `MeshOperator` — `domain_coefficients` (row 39) + `jump_hooks` (row 40)
# + `boundary_conditions` (row 41) + `config.dg.penalty_parameter` over EllipticCore's
# `MeshOperator` (row 15) — plus the three things a multi-level consumer needs around it:
# `operator_builder` (the `MultigridHierarchy` `build_op` factory, taking the caller's
# per-level `worldtube` bundle factory — rebuilding that bundle is row 38's
# `worldtube_data(...; mesh)`, deliberately NOT re-implemented here) and `coarsest_operator`
# (the multigrid COARSEST level — what main.tex:622-626 Fig. 2 plots),
# with `linearized_rhs(op, wt)` = `worldtube_massive_source(wt) − operator_on_zero(op)`.
#
# WHY THESE GATES.
#  (a) BITWISE identity with the manual spelling. The whole value of the helper is that it is
#      not a *different* operator: `build_matrix`, `operator_on_zero` and `apply!` must all be
#      bit-for-bit equal to the hand-composed `MeshOperator`. Anything weaker would let the
#      helper drift from what the sealed node tests measure.
#  (b) The row-40 ORACLE reproduced THROUGH the helper: the vtu junction term IS in the
#      matrix (|A_full − A_novtu|_∞ = 9.0137e-3 at the production point) while the worldtube
#      term is BITWISE absent from it and lives entirely in `operator_on_zero`. That asymmetry
#      is what the figure's red rectangle (worldtube) vs red dotted lines (vtu) mean, so the
#      figure driver must be able to rely on it at the level it plots.
#  (c) The COARSEST LEVEL is the yaml's, not a guess: `InitialCoarseLevels: L` (yaml:188)
#      counts the coarser levels only, so the hierarchy has `L + 1` levels and at the
#      production `L = 0` the coarsest level IS the initial 20-element / 320-dof mesh. Both
#      that default and the fully-coarsened block structure (10/160 on the 5×2 lattice,
#      15/240 on the 5×3) are pinned, and the returned operator is gated BITWISE against the
#      sealed `MultigridHierarchy`'s own `levels[end]`.
#  (d) Mesh-boundness is enforced, not assumed: the per-level bundle works on a level mesh,
#      a stale one is refused, a foreign mesh is refused.
#  (e) The parameter envelope the CRP needs: a ∈ {0, 0.5}, m ∈ {0, 1, 2}, L ∈ {0, 1}, both
#      angular charts, and BigFloat(256) precision inheritance.

using Test
using LinearAlgebra
using ScalarSelfForce
using EllipticCore

# ---- configs (as_-prefixed: sibling test files share the runtests.jl include scope) --------
# `block_domain` REQUIRES the inner radial bound to BE r₊, which the yaml's inner_radius = 2
# satisfies only at a = 0 — every spin moves it with the horizon.
as_cfg(; spin = 0.0, r0 = 10.0, m = 2, L = 0) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, L = L)
as_orb(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)
# the 5×3 full-z twin (crp_paper_config hard-codes the yaml's ImposeEquatorialSymmetry = true)
function as_full_z(cfg::RunConfig{T}) where {T}
    OC = typeof(cfg.orbit)
    oc = OC((f === :impose_equatorial_symmetry ? false : getfield(cfg.orbit, f)
             for f in fieldnames(OC))...)
    return RunConfig{T}((f === :orbit ? oc : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

# the hand-composed operator this node must reproduce bitwise
function as_manual(T, dom, wt, orb, cfg, m; mesh = dom.mesh)
    fcs = domain_coefficients(dom, orb, cfg)
    return MeshOperator(T, mesh, fcs.A, fcs.beta, fcs.gamma;
                        C = cfg.dg.penalty_parameter,
                        hooks = jump_hooks(dom, wt, orb, cfg, m; mesh = mesh),
                        bcs = boundary_conditions(dom, orb, cfg; mesh = mesh))
end

as_probe(n) = ComplexF64[cis(0.37k) / (1 + 0.01k) for k in 1:n]

# the per-level WorldtubeData factory the multi-level helpers take (row 38's constructor —
# this node deliberately owns no WorldtubeData rebuild of its own)
as_wtf(dom, orb, cfg, m) = mesh -> worldtube_data(dom, orb, cfg, m; mesh = mesh)

@testset "operator assembly: production MeshOperator + multigrid levels (fig-matrix)" begin

    # ============================ (a) the production composition ===========================
    cfg = as_cfg()
    orb = as_orb(cfg)
    dom = block_domain(orb, cfg; orders = (4, 4))
    wt = worldtube_data(dom, orb, cfg, 2)

    @test length(dom.mesh.elements) == 20        # yaml:78-110 RefinedLevels at L = 0
    @test dom.mesh.ndofs == 320                  # 20 elements × 4×4 grid points

    op = assemble_operator(ComplexF64, dom, wt, orb, cfg, 2)
    opm = as_manual(ComplexF64, dom, wt, orb, cfg, 2)
    A = build_matrix(op)
    Am = build_matrix(opm)
    @test size(A) == (320, 320)
    @test A == Am                                                     # BITWISE
    @test operator_on_zero(op) == operator_on_zero(opm)               # BITWISE
    u = as_probe(dom.mesh.ndofs)
    y1 = similar(u); y2 = similar(u)
    apply!(y1, op, u); apply!(y2, opm, u)
    @test y1 == y2                                                    # BITWISE
    @test all(isfinite, y1)
    # the mode number defaults to the config's
    @test build_matrix(assemble_operator(ComplexF64, dom, wt, orb, cfg)) == A
    # …and a disagreeing one is refused by the sealed guards
    @test_throws ArgumentError assemble_operator(ComplexF64, dom, wt, orb, cfg, 1)

    # ============================ (b) the row-40 oracle through the helper ==================
    fcs = domain_coefficients(dom, orb, cfg)
    bcs = boundary_conditions(dom, orb, cfg)
    mk(h) = MeshOperator(ComplexF64, dom.mesh, fcs.A, fcs.beta, fcs.gamma;
                         C = cfg.dg.penalty_parameter, hooks = h, bcs = bcs)
    nojd = Dict{Tuple{ElementId,Int,Int},Vector{FaceState{ComplexF64}}}()
    jf = vtu_junction_faces(dom, orb, cfg, 2)
    A_novtu = build_matrix(mk(jump_hooks(dom.mesh, worldtube_jump_data(wt),
                                         Dict{Tuple{ElementId,Int,Int},VtuJunction{Float64}}())))
    A_nowt  = build_matrix(mk(jump_hooks(dom.mesh, nojd, jf)))
    dvtu = maximum(abs, A .- A_novtu)
    dwt  = maximum(abs, A .- A_nowt)
    @info "assemble (b) row-40 oracle at the coarsest level: |A−A_novtu|∞ = $dvtu  " *
          "|A−A_nowt|∞ = $dwt  |operator_on_zero|∞ = $(maximum(abs, operator_on_zero(op)))"
    @test dvtu > 1e-6                        # the junction term IS in the GMRES matrix
    @test isapprox(dvtu, 9.013699828708553e-3; rtol = 1e-9)
    @test A == A_nowt                        # the worldtube term is BITWISE absent from it
    @test dwt == 0
    @test maximum(abs, operator_on_zero(op)) > 1e-3   # …and lives entirely in the RHS lift

    # ============================ (c) linearized_rhs ========================================
    b = linearized_rhs(op, wt)
    @test b == worldtube_massive_source(wt) .- operator_on_zero(op)
    @test length(b) == dom.mesh.ndofs
    @test all(isfinite, b)
    # the general EllipticCore fold is the same function (one arithmetic, two methods)
    @test linearized_rhs(op, worldtube_massive_source(wt)) == b
    # skipping the fold is an O(1) error, not a small one
    @test maximum(abs, b .- worldtube_massive_source(wt)) > 1e-3

    # ============================ (d) the coarsest multigrid level =========================
    # yaml:188 InitialCoarseLevels = L counts the COARSER levels only ⇒ L + 1 levels total,
    # so at the production L = 0 the coarsest level IS the initial mesh: the 20-element,
    # 320-dof lattice the paper's Fig. 2 numbers 1…20 (main.tex:622-626).
    @test cfg.multigrid.initial_coarse_levels == 0
    wtf = as_wtf(dom, orb, cfg, 2)
    c0 = coarsest_operator(ComplexF64, dom, orb, cfg, 2; worldtube = wtf)
    @test c0.nlevels == 1
    @test c0.mesh === dom.mesh
    @test length(c0.mesh.elements) == 20
    @test c0.mesh.ndofs == 320
    @test build_matrix(c0.op) == A
    @test worldtube_massive_source(c0.wt) == worldtube_massive_source(wt)

    # coarsening all the way down instead: the irreducible block structure, one element/block
    cfull = coarsest_operator(ComplexF64, dom, orb, cfg, 2; worldtube = wtf, maxlevels = nothing)
    @test cfull.nlevels == 2
    @test length(cfull.mesh.elements) == 10          # 5×2 lattice
    @test cfull.mesh.ndofs == 160
    yb = similar(as_probe(160)); apply!(yb, cfull.op, as_probe(160))
    @test all(isfinite, yb)
    @test size(build_matrix(cfull.op)) == (160, 160)
    @test_throws ArgumentError coarsest_operator(ComplexF64, dom, orb, cfg, 2;
                                                 worldtube = wtf, maxlevels = 0)

    # BITWISE identity with the sealed MultigridHierarchy's own coarsest level
    bop = operator_builder(ComplexF64, dom, orb, cfg, 2; worldtube = wtf)
    hier = MultigridHierarchy(op, bop; maxlevels = nothing)
    @test length(hier.levels) == 2
    @test hier.levels[end].mesh.elements == cfull.mesh.elements
    @test build_matrix(hier.levels[end].op) == build_matrix(cfull.op)
    @test operator_on_zero(hier.levels[end].op) == operator_on_zero(cfull.op)
    h1 = MultigridHierarchy(op, bop; maxlevels = cfg.multigrid.initial_coarse_levels + 1)
    @test length(h1.levels) == 1
    @test build_matrix(h1.levels[end].op) == A

    # ============================ (e) mesh-boundness ========================================
    cmesh = cfull.mesh
    wtc = cfull.wt
    @test length(worldtube_massive_source(wtc)) == cmesh.ndofs
    @test worldtube_elements(wtc) == filter(el -> inside_worldtube(dom, el), cmesh.elements)
    @test !isempty(worldtube_jump_data(wtc))
    @test build_matrix(assemble_operator(ComplexF64, dom, wtc, orb, cfg, 2; mesh = cmesh)) ==
          build_matrix(cfull.op)
    # a STALE WorldtubeData on a level mesh is REFUSED, never silently misapplied
    @test_throws ArgumentError assemble_operator(ComplexF64, dom, wt, orb, cfg, 2; mesh = cmesh)
    # a FOREIGN mesh (blocks that are not this lattice's) is refused by the sealed guards
    fmesh = Mesh2D(Block(91, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0)),
                   [ElementId(91, SegmentId(0, 0), SegmentId(0, 0))], (4, 4))
    # the sealed jump_hooks mesh-boundness guard fires first: the fine mesh's element ids are
    # not elements of `fmesh`, so the datum cannot be bound to it
    @test_throws ArgumentError assemble_operator(ComplexF64, dom, wt, orb, cfg, 2; mesh = fmesh)

    # ============================ (f) the parameter envelope ===============================
    for a in (0.0, 0.5), mm in (0, 1, 2)
        c = as_cfg(spin = a, m = mm)
        o = as_orb(c)
        d = block_domain(o, c; orders = (4, 4))
        w = worldtube_data(d, o, c, mm)
        oo = assemble_operator(ComplexF64, d, w, o, c, mm)
        @test build_matrix(oo) == build_matrix(as_manual(ComplexF64, d, w, o, c, mm))
        uu = as_probe(d.mesh.ndofs); yy = similar(uu)
        apply!(yy, oo, uu)
        @test all(isfinite, yy)
        @test all(isfinite, linearized_rhs(oo, w))
        # m = 0 degenerates the junction coefficient to exactly zero, but the table is built
        mm == 0 && @test all(j -> iszero(j.coeff), values(vtu_junction_faces(d, o, c, 0)))
    end

    # L = 1: one initial coarse level ⇒ the coarsest level is the L = 0 lattice's size
    cfg1 = as_cfg(L = 1)
    orb1 = as_orb(cfg1)
    dom1 = block_domain(orb1, cfg1; orders = (4, 4))
    wt1 = worldtube_data(dom1, orb1, cfg1, 2)
    @test length(dom1.mesh.elements) == 80
    @test cfg1.multigrid.initial_coarse_levels == 1
    c1 = coarsest_operator(ComplexF64, dom1, orb1, cfg1, 2; worldtube = as_wtf(dom1, orb1, cfg1, 2))
    @test c1.nlevels == 2
    @test length(c1.mesh.elements) == 20
    @test c1.mesh.ndofs == 320
    y1c = similar(as_probe(320)); apply!(y1c, c1.op, as_probe(320))
    @test all(isfinite, y1c)

    # the 5×3 full-z chart
    cfgz = as_full_z(as_cfg())
    orbz = as_orb(cfgz)
    domz = block_domain(orbz, cfgz; orders = (4, 4))
    wtz = worldtube_data(domz, orbz, cfgz, 2)
    @test domz.chart === :z
    @test length(domz.mesh.elements) == 30
    opz = assemble_operator(ComplexF64, domz, wtz, orbz, cfgz, 2)
    @test build_matrix(opz) == build_matrix(as_manual(ComplexF64, domz, wtz, orbz, cfgz, 2))
    cz = coarsest_operator(ComplexF64, domz, orbz, cfgz, 2;
                           worldtube = as_wtf(domz, orbz, cfgz, 2), maxlevels = nothing)
    @test length(cz.mesh.elements) == 15
    @test cz.mesh.ndofs == 240

    # ============================ (g) BigFloat(256) inheritance ============================
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        orbB = CircularOrbit(big(0.0), big(10.0))
        cfgB = crp_paper_config(; spin = big(0.0), orbital_radius = big(10.0), m_mode = 2,
                                inner_radius = orbB.rp, outer_radius = big(1.0e4),
                                worldtube_radius = big(4.8), worldtube_dz2 = big(0.25))
        @test cfgB isa RunConfig{BigFloat}
        domB = block_domain(orbB, cfgB; orders = (4, 4))
        # the fully coarsened level keeps the BigFloat cost to the 10-element mesh
        cB = coarsest_operator(Complex{BigFloat}, domB, orbB, cfgB, 2;
                               worldtube = as_wtf(domB, orbB, cfgB, 2), maxlevels = nothing)
        @test length(cB.mesh.elements) == 10
        @test cB.mesh.ndofs == 160
        uB = Complex{BigFloat}[big(1) / (k + 1) for k in 1:cB.mesh.ndofs]
        yB = similar(uB)
        apply!(yB, cB.op, uB)
        @test all(isfinite, yB)
        @test eltype(yB) === Complex{BigFloat}
        @test precision(real(yB[1])) == 256
        bB = linearized_rhs(cB.op, cB.wt)
        @test eltype(bB) === Complex{BigFloat}
        @test all(isfinite, bB)
    finally
        setprecision(BigFloat, old)
    end

    # ============================ (h) API surface + documentation ==========================
    for f in (:assemble_operator, :operator_builder, :coarsest_operator, :linearized_rhs)
        @test isdefined(ScalarSelfForce, f)
        @test f in names(ScalarSelfForce)
    end
    d_as = string(@doc assemble_operator)
    @test occursin("domain_coefficients", d_as) && occursin("jump_hooks", d_as) &&
          occursin("boundary_conditions", d_as)
    d_co = string(@doc coarsest_operator)
    @test occursin("InitialCoarseLevels", d_co)
    d_ob = string(@doc operator_builder)
    @test occursin("MESH-BOUND", d_ob) && occursin("worldtube_data", d_ob)
end

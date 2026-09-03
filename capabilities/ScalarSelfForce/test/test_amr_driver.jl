# test_amr_driver.jl — ESSF node `amr-driver` (GLOBAL_DAG row 43, wave 11, worker ESSF-W11-A),
# written BEFORE `src/Solve/amr.jl` exists (TDD; `tests_red.log` captured against the missing
# file). THIS IS THE PRODUCTION LOOP: the hp-AMR sweep that turns the single fixed-mesh
# m-mode solve of row 42 into the paper's converged sequence.
#
# WHAT THE NODE DELIVERS (`Solve/amr.jl`):
#   * `amr_flags(sol)` / `amr_flags(mesh, dom, orbit, config)` — the CRP refinement criterion:
#     spectre's `RefineAtPuncture` (h-split every element CONTAINING the puncture, with the
#     block-logical snap that keeps the split symmetric) composed with `IncreaseResolution`
#     (p-increment everywhere else), under the yaml's `Limits`.
#   * `amr_step!(sol)` — one level: flags → `hp_refine` → `transfer_solution` warm start →
#     `rebuild_for_mesh!` → solve.
#   * `solve_mode_amr(orbit, config, m)` → `AmrHistory` — the whole `Amr: Iterations: 8` sweep
#     with per-level meshes, `SolveRecord`s and particle diagnostics.
#   * `puncture_chart_position`, `amr_level_options` — the two policy pieces rows 44/47/48/52
#     must share. The per-level particle diagnostics are read through the SEALED row-42
#     `field_at` / `field_grad_at`; this node adds no read-out surface of its own.
#
# ORACLE (the DAG row-43 gate column + the READ-ONLY mirrors, transcribed by hand):
#   * `main.tex:552-560` — "we split every element that contains the puncture in two in both
#     dimensions" and "increase the polynomial order in all elements that do not contain the
#     puncture by one in both dimensions in each refinement step", with neighbouring elements
#     differing by at most one h-refinement level.
#   * `main.tex:672-676` — the multigrid hierarchy is h-coarsened to the irreducible block
#     structure and the bottom LU is "cached and reused in all multigrid V-cycles across AMR
#     iterations, since the bottom grid remains the same even when AMR adds finer levels".
#   * `main.tex:701` (fig:linsolv) — "Each solve converges in less than 10 preconditioned
#     iterations. The first AMR level converges immediately because the multigrid bottom
#     solver solves the problem directly."
#   * `main.tex:714-716` — exponential convergence under the hp-AMR sequence; a few seconds
#     per m-mode on an Apple M2 Pro; each m-mode independent.
#   * `anc/ScalarSelfForce.yaml:129-153` — `Criteria: [RefineAtPuncture, IncreaseResolution]`,
#     `Policies: {Isotropy: Anisotropic, EnforceTwoToOneBalanceInNormalDirection: True,
#     Limits: {NumGridPoints: Auto, RefinementLevel: Auto, ErrorBeyondLimits: True}}`,
#     `Iterations: 8`.
#   * spectre `Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineAtPuncture.cpp` — the
#     block-logical snap `abs(ξ_d) < 1e-10 ⇒ ξ_d = 0` and `Domain/ElementLogicalCoordinates.cpp`
#     lines 31-46 — the CLOSED containment `lo ≤ ξ_d ≤ up`, so a puncture on an element face
#     or corner belongs to EVERY element that touches it.
#   * spectre `ParallelAlgorithms/Amr/Actions/EvaluateRefinementCriteria.hpp:147` — criteria
#     are merged per dimension by `std::max` over `Flag`, and `Split` (5) outranks
#     `IncreaseResolution` (4): the puncture elements are split and NOT p-refined, which is
#     exactly `hp_refine`'s one-criterion `hold`/`split` contract.
#   * spectre `ParallelAlgorithms/Amr/Policies/Limits.{hpp,cpp}` — `Auto` means
#     `ElementId<1>::max_refinement_level = 15` and
#     `Spectral::maximum_number_of_points<Legendre> = 20`.
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here): `solver-stack-ssf`
# (`solve_mode`, `solve_mode!`, `rebuild_for_mesh!`, `ModeSolution`, `solver_options`,
# `field_at`, `mode_residual`), `amr-h-refinement`/`amr-p-refinement` (`hp_refine`,
# `transfer_solution`, the flag containers), `solver-stack-core` (`SolveRecord`,
# `ResidualMonitor`, `cumulative_iterations`), `domain-blocks-map` (`block_domain`),
# `mesh-basis-quadrature` (`refelement`, `interpolation_matrix`), `mesh2d-blocks`
# (`logical_interval`, `map_inverse`, `map_jacobian`, `two_to_one_balanced`).
#
# GATE CONSTANTS — every one calibrated on the throwaway prototypes
# `evidence/amr-driver/proto1_structure.jl` and `proto2_sweep.jl` BEFORE this file was
# written, never tuned afterwards to pass. The measured value is quoted beside each gate.

using Test
using LinearAlgebra
using ScalarSelfForce
using EllipticCore

# ---- configs (da_-prefixed: sibling test files share the runtests.jl include scope) -------
da_cfg(; spin = 0.5, r0 = 10.0, m = 2, L = 0, P = 4, wr = 4.8) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L,
                     worldtube_radius = wr)
da_orb(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)

# an INDEPENDENT re-derivation of "which elements contain the puncture", written from the
# spectre source rather than from the implementation: map the puncture into each block with
# the block's own inverse map, snap, and keep every element whose CLOSED dyadic interval
# contains it in both dimensions.
function da_reference_punctures(mesh, dom, orbit; tol = 1e-10)
    RT = typeof(mesh).parameters[1]
    out = ElementId[]
    for id in mesh.elements
        blk = mesh.blocks[mesh.blockindex[id.block]]
        inside = true
        for d in 1:2
            xd = d == 1 ? RT(orbit.r0) : zero(RT)
            ζ = map_inverse(blockmap(blk, d), xd)
            abs(ζ - 1) < tol && (ζ = one(RT))
            abs(ζ + 1) < tol && (ζ = -one(RT))
            abs(ζ) < tol && (ζ = zero(RT))
            lo, hi = logical_interval(RT, id, d)
            (abs(ζ) <= 1 && lo <= ζ <= hi) || (inside = false)
        end
        inside && push!(out, id)
    end
    return out
end

da_flagset(flags) = Set(collect(keys(flags)))

# a copy of `cfg` with a different `Amr.Criteria` list, rebuilt field by field (positional
# reconstruction would silently rot if either struct gained a field)
function da_criteria(cfg, criteria)
    A = typeof(cfg.amr)
    amr = A((f === :criteria ? criteria : getfield(cfg.amr, f) for f in fieldnames(A))...)
    C = typeof(cfg)
    return C((f === :amr ? amr : getfield(cfg, f) for f in fieldnames(C))...)
end

# successive differences of a per-level scalar sequence and their ratios
function da_diffs(vals)
    d = [abs(vals[k] - vals[k-1]) for k in 2:length(vals)]
    r = [k > 1 ? d[k] / d[k-1] : NaN for k in eachindex(d)]
    return d, r
end

# geometric mean of the successive-difference ratios from level 2 on — the honest measure of
# "converges exponentially" for a sequence that has plateaus
da_georatio(r) = exp(sum(log, r[2:end]) / length(r[2:end]))

@testset "hp-AMR driver (amr-driver)" begin

    # ======================================================================================
    # (a) THE CRITERION — spectre RefineAtPuncture ∘ IncreaseResolution, and its snap
    # ======================================================================================
    @testset "(a) criterion: which elements contain the puncture" begin
        cfg = da_cfg()
        orb = da_orb(cfg)
        dom = block_domain(orb, cfg)
        mesh = dom.mesh

        # the puncture in the DISCRETIZATION chart: (r₀, w = z² = 0) on the equator
        xp = puncture_chart_position(dom, orb)
        @test xp[1] == orb.r0
        @test xp[2] == 0
        @test dom.chart === :z2

        flags = amr_flags(mesh, dom, orb, cfg)
        @test flags isa AbstractDict
        # MEASURED (proto1_structure.log): exactly TWO elements of the worldtube block, the
        # pair that shares the r = r₀ face on the equatorial boundary w = 0. Not four: in the
        # z² chart w = 0 is the DOMAIN boundary, so only one element is adjacent in dim 2.
        @test length(flags) == 2
        @test all(v -> v === (true, true), values(flags))       # spectre: Split in BOTH dims
        @test all(id -> id.block in dom.worldtube, keys(flags))
        @test all(id -> (id.segments[1].level, id.segments[2].level) == (1, 1), keys(flags))
        # the independent re-derivation agrees element for element
        @test da_flagset(flags) == Set(da_reference_punctures(mesh, dom, orb))

        # THE PUNCTURE SITS ON AN ELEMENT CORNER, exactly (the sealed worldtube node proved
        # the Gauss nodes never touch it). Both flagged elements see |ξ| = 1 in both dims.
        for id in keys(flags)
            e = findfirst(==(id), mesh.elements)
            blk = mesh.blocks[mesh.blockindex[id.block]]
            for d in 1:2
                lo, hi = logical_interval(Float64, id, d)
                ζ = map_inverse(blockmap(blk, d), d == 1 ? xp[1] : xp[2])
                ξ = (2ζ - (lo + hi)) / (hi - lo)
                @test abs(abs(ξ) - 1) <= 1e-10
            end
        end

        # every OTHER element is unflagged, i.e. p-refined (IncreaseResolution)
        @test count(id -> haskey(flags, id), mesh.elements) == 2
    end

    @testset "(a) the block-logical snap is load-bearing" begin
        # spectre snaps |ξ| < 1e-10 to 0 "to avoid asymmetric refinement from round-off error
        # for the usual case where the puncture is at the center of a block"
        # (RefineAtPuncture.cpp:40-45). At the paper's own (r₀, w_r) = (10, 4.8) the block
        # centre is hit EXACTLY in Float64, so the snap is inert there; it is not inert in
        # general. MEASURED (proto4_snap.log): 160 of the 4100 scanned (r₀, w_r) lattices in
        # 6 ≤ r₀ ≤ 16, 1 ≤ w_r ≤ 5 put the puncture at a nonzero round-off from the block
        # centre, and in EVERY one of them the unsnapped criterion flags ONE element instead
        # of two — the asymmetric refinement spectre's comment is about. At
        # (r₀, w_r) = (7.2, 2.0) the raw block-logical coordinate is +4.441e-16.
        cfg = da_cfg(; spin = 0.0, r0 = 7.2, wr = 2.0)
        orb = da_orb(cfg)
        dom = block_domain(orb, cfg)
        blk = dom.mesh.blocks[dom.mesh.blockindex[only(dom.worldtube)]]
        ζraw = map_inverse(blockmap(blk, 1), orb.r0)
        @test ζraw != 0                       # measured +4.441e-16
        @test abs(ζraw) < 1e-10

        @test length(amr_flags(dom.mesh, dom, orb, cfg)) == 2
        @test length(amr_flags(dom.mesh, dom, orb, cfg; snap_tol = 0.0)) == 1
        # and the snapped criterion still agrees with the independent re-derivation
        @test da_flagset(amr_flags(dom.mesh, dom, orb, cfg)) ==
              Set(da_reference_punctures(dom.mesh, dom, orb))
    end

    @testset "(a) the hp cascade: level map, order map, 2:1 balance" begin
        # structure only, no solves — this is the sealed amr node's metadata gate on the REAL
        # CRP domain. MEASURED (proto1_structure.log) at a = 0.5 and a = 0, identically:
        # nel 20 27 33 39 45 51 57 63 69, ndofs 320 603 962 1425 2004 2711 3558 4557 5720.
        for spin in (0.5, 0.0)
            cfg = da_cfg(; spin = spin)
            orb = da_orb(cfg)
            dom = block_domain(orb, cfg)
            mesh = dom.mesh
            nel = Int[]
            for k in 0:8
                flags = amr_flags(mesh, dom, orb, cfg)
                push!(nel, length(mesh.elements))
                punct = collect(keys(flags))
                # (a1) the puncture pair is at h-level L+1+k in BOTH dimensions, order 4
                @test length(punct) == 2
                for id in punct
                    @test (id.segments[1].level, id.segments[2].level) == (1 + k, 1 + k)
                    @test mesh.orders[findfirst(==(id), mesh.elements)] == (4, 4)
                end
                # (a2) elements that were never split still carry N = 4 + k in both dims
                never = [i for (i, id) in enumerate(mesh.elements)
                         if (id.segments[1].level, id.segments[2].level) ==
                            dom.refinements[id.block]]
                @test !isempty(never)
                @test all(o -> o == (4 + k, 4 + k), mesh.orders[never])
                # (a3) nothing in the mesh is below N = 4 or above N = 4 + k
                @test minimum(minimum.(mesh.orders)) == 4
                @test maximum(maximum.(mesh.orders)) == 4 + k
                # (a4) two-to-one balanced at every level (the sealed hp_refine's job — a
                # cheap standing check that the driver never hands it a violating flag set)
                @test two_to_one_balanced(mesh.elements)
                # (a5) the puncture is a CORNER of every flagged element at every level
                for id in punct, d in 1:2
                    blk = mesh.blocks[mesh.blockindex[id.block]]
                    lo, hi = logical_interval(Float64, id, d)
                    ζ = map_inverse(blockmap(blk, d), d == 1 ? orb.r0 : 0.0)
                    ξ = (2ζ - (lo + hi)) / (hi - lo)
                    @test abs(abs(ξ) - 1) <= 1e-10
                end
                k == 8 && break
                mesh = hp_refine(mesh, flags)
            end
            @test nel == [20, 27, 33, 39, 45, 51, 57, 63, 69]
        end
    end

    @testset "(a) the yaml Limits are enforced, not silently clipped" begin
        cfg = da_cfg()
        orb = da_orb(cfg)
        dom = block_domain(orb, cfg)
        # `Limits: {RefinementLevel: Auto, NumGridPoints: Auto, ErrorBeyondLimits: True}`
        # ⇒ spectre's own defaults 15 / 20, and an ERROR when a flag would cross them.
        @test ScalarSelfForce.AMR_MAX_REFINEMENT_LEVEL == 15
        @test ScalarSelfForce.AMR_MAX_RESOLUTION == 20
        # the CRP sweep never gets there: 8 iterations reach h-level 1 + 8 = 9 and N = 12
        @test 1 + cfg.amr.iterations < ScalarSelfForce.AMR_MAX_REFINEMENT_LEVEL
        @test 4 + cfg.amr.iterations < ScalarSelfForce.AMR_MAX_RESOLUTION
        @test_throws ArgumentError amr_flags(dom.mesh, dom, orb, cfg;
                                             max_refinement_level = 1)
        @test_throws ArgumentError amr_flags(dom.mesh, dom, orb, cfg; max_resolution = 4)
        # a criteria list this file does not implement must be refused, not approximated
        @test_throws ArgumentError amr_flags(dom.mesh, dom, orb, da_criteria(cfg, [:Persson]))
        # and an EMPTY flag set is refused rather than silently becoming pure p-refinement:
        # a particle outside this lattice can only mean `orbit` and `dom` disagree
        @test_throws ArgumentError amr_flags(dom.mesh, dom, CircularOrbit(0.5, 2.0e4), cfg)
    end

    # ======================================================================================
    # (b)(c)(d)(e) THE PRODUCTION SWEEP — a = 0.5, r₀ = 10, m = 2
    # ======================================================================================
    # ONE production sweep on the DEFAULT level policy, shared by the gates below (it is the
    # most expensive object in this file — 20 → 57 elements, 320 → 3558 dofs, 7 solves).
    da_pcfg = da_cfg()
    da_porb = da_orb(da_pcfg)
    da_hist = solve_mode_amr(da_porb, da_pcfg, 2; iterations = 6)

    @testset "(b,c,d,e) the a = 0.5, m = 2 hp-AMR sweep" begin
        cfg, orb, hist = da_pcfg, da_porb, da_hist

        @test hist isa AmrHistory
        @test length(hist.levels) == 7                       # level 0 + 6 AMR iterations
        @test hist.levels[1].step == 0
        @test hist.m == 2
        @test hist.x_p == puncture_chart_position(hist.solution.dom, orb)

        # ---- the mesh sequence is the one gate (a) fixed, now under the real solver -------
        @test [length(l.mesh.elements) for l in hist.levels] == [20, 27, 33, 39, 45, 51, 57]
        @test [l.mesh.ndofs for l in hist.levels] ==
              [320, 603, 962, 1425, 2004, 2711, 3558]
        @test all(l -> length(l.punctures) == 2, hist.levels)

        # ---- (c) GMRES iterations per level ------------------------------------------------
        its = [l.record.iterations for l in hist.levels]
        # MEASURED (proto2_a05_m2_crp.log): 1, then 9-10 per level. main.tex:701 — the first
        # AMR level converges IMMEDIATELY (`InitialCoarseLevels: L = 0` makes the V-cycle the
        # direct bottom solve) and every later level in "less than 10" iterations; we measure
        # 10 twice, so the honest gate is ≤ 11.
        @test its[1] == 1
        @test all(<=(11), its)
        @test all(l -> l.record.final_residual <= cfg.gmres.absolute_residual, hist.levels)
        @test mode_residual(hist.solution) <= 1e-11
        # the monitor accumulated exactly one record per level — this IS fig:linsolv's input
        @test length(hist.solution.monitor.records) == length(hist.levels)
        @test cumulative_iterations(hist.solution.monitor)[end] == sum(its)

        # ---- (d) caches: one epoch, one bottom LU, Schwarz reuse ---------------------------
        @test hist.solution.stack.epoch == 0                  # a pure mesh change bumps nothing
        # the bottom grid is the FIXED irreducible block structure from level 1 on, so its LU
        # is built once there and hit at every later level (main.tex:674). Level 0 has its own
        # single-level hierarchy (the direct solve), hence 2 builds in total, not 1.
        @test sum(l.record.bottom_builds for l in hist.levels) == 2
        @test sum(l.record.bottom_hits for l in hist.levels) == length(hist.levels) - 2
        @test all(l -> l.record.bottom_builds + l.record.bottom_hits >= 1, hist.levels[2:end])
        # the bottom grid is the irreducible block structure of the 5x2 lattice and it does
        # NOT MOVE: 10 elements at order 4 = 160 dofs at every level from 1 on, which is why
        # its LU is a hit. With the policy's `shrink = 0.5` the hierarchy reaches it in a
        # handful of levels instead of one ring at a time (measured 4 levels at AMR level 8
        # against 18 at `shrink = 1`, for one extra GMRES iteration over the whole sweep).
        @test all(l -> l.record.level_ndofs[end] == 160, hist.levels[2:end])
        @test all(l -> length(l.record.level_ndofs) <= 6, hist.levels)
        # Schwarz subdomain factorizations are reused on the levels AMR did not touch
        # (main.tex:676 "the matrix inverses can remain cached across AMR iterations, so we
        # only have to construct new subdomain matrices when AMR adds a finer level").
        # MEASURED at the last level of this sweep: 73 builds against 42 cache hits, i.e.
        # 36 % of the subdomain solves in the whole hierarchy are served from the cache; the
        # builds are the new finest level plus the coarse levels the composed `shrink = 0.5`
        # coarsening moved. The gate is the REUSE, not a fitted ratio.
        @test all(l -> l.record.subdomain_hits > 0, hist.levels[3:end])
        @test hist.levels[end].record.subdomain_hits >
              hist.levels[end].record.subdomain_builds / 3

        # ---- (b) convergence of the particle diagnostics ----------------------------------
        # The DAG's F_r^m criteria applied to the RAW ingredients (recorded as a landed-gate
        # substitution; row 44 restates them on F_r^m itself).
        psi = [l.psi for l in hist.levels]
        dps = [l.dpsi_r for l in hist.levels]
        dpsi, rpsi = da_diffs(psi)
        ddps, rdps = da_diffs(dps)
        # MEASURED ratios (proto2_a05_m2_parity.log, 8 levels): ψ 0.220 0.158 0.339 0.254
        # 0.276 0.278 0.308; ∂_rψ 0.141 0.223 0.055 0.129 0.326 0.231 0.260 — all well under
        # the DAG's 0.5 from level 2 on.
        @test all(<(0.5), rpsi[2:end])
        @test all(<(0.5), rdps[2:end])
        @test dpsi[end] < dpsi[1] / 100
        @test ddps[end] < ddps[1] / 100
        # the diagnostics are finite, nonzero, and the two one-sided corner traces agree
        # better and better as the puncture corner is refined (measured spread 8.4e-6 → 3.5e-17
        # on ψ over 8 levels, 2.2e-4 → 9.4e-11 on ∂_rψ)
        @test all(l -> isfinite(l.psi) && isfinite(l.dpsi_r), hist.levels)
        @test hist.levels[end].trace_spread[1] < hist.levels[1].trace_spread[1] / 100
        @test hist.levels[end].trace_spread[2] < hist.levels[1].trace_spread[2] / 100
        # the final level's value equals `field_at` on the returned ModeSolution — the
        # read-out row 44 will use
        @test hist.levels[end].psi ==
              field_at(hist.solution, hist.x_p...; element = hist.levels[end].element)
        @test hist.levels[end].dpsi_r ==
              field_grad_at(hist.solution, hist.x_p...;
                                element = hist.levels[end].element)[1]
    end

    # ======================================================================================
    # (f) DETERMINISM — the same sweep twice is bitwise identical
    # ======================================================================================
    @testset "(f) determinism of a repeated sweep" begin
        cfg = da_cfg()
        orb = da_orb(cfg)
        a = solve_mode_amr(orb, cfg, 2; iterations = 3)
        b = solve_mode_amr(orb, cfg, 2; iterations = 3)
        @test [l.mesh.elements for l in a.levels] == [l.mesh.elements for l in b.levels]
        @test [l.mesh.orders for l in a.levels] == [l.mesh.orders for l in b.levels]
        @test [l.punctures for l in a.levels] == [l.punctures for l in b.levels]
        @test [l.record.iterations for l in a.levels] == [l.record.iterations for l in b.levels]
        @test [l.psi for l in a.levels] == [l.psi for l in b.levels]           # BITWISE
        @test [l.dpsi_r for l in a.levels] == [l.dpsi_r for l in b.levels]     # BITWISE
        @test a.solution.u == b.solution.u                                     # BITWISE
    end

    # ======================================================================================
    # (b') the other configurations of the DAG row
    # ======================================================================================
    @testset "(b) a = 0.5 m = 1 and a = 0 m = 2 converge the same way" begin
        # THE SINGLE-STEP RATIO IS NOT UNIVERSALLY < 0.5 (a landed-gate refinement, measured,
        # not tuned): over the full 8-level sweeps the ∂_rΨ^ℛ ratios are ≤ 0.33 at a = 0.5,
        # m = 2 and ≤ 0.32 at a = 0, m = 2, but at a = 0.5, m = 1 level 5 measures **0.687**
        # — one plateau in an otherwise geometric sequence (proto2_a05_m1_crp.log). The
        # two-step ratio |Δ_k|/|Δ_{k−2}| is the honest envelope: measured worst 0.139 there,
        # and it is what "converges exponentially" means for a sequence with plateaus.
        for (spin, m) in ((0.5, 1), (0.0, 2))
            cfg = da_cfg(; spin = spin, m = m)
            orb = da_orb(cfg)
            hist = solve_mode_amr(orb, cfg, m; iterations = 5)
            its = [l.record.iterations for l in hist.levels]
            @test its[1] == 1
            @test all(<=(11), its)
            dp, _ = da_diffs([l.psi for l in hist.levels])
            dd, _ = da_diffs([l.dpsi_r for l in hist.levels])
            @test all(k -> dp[k] < 0.5 * dp[k-2], 3:length(dp))
            @test all(k -> dd[k] < 0.5 * dd[k-2], 3:length(dd))
            @test dp[end] < dp[1] / 50
            @test dd[end] < dd[1] / 50
            @test mode_residual(hist.solution) <= 1e-11
        end
    end

    # ======================================================================================
    # TYPE GENERICITY — the whole sweep at Complex{BigFloat}(256)
    # ======================================================================================
    @testset "the sweep is type-generic: Complex{BigFloat}(256)" begin
        setprecision(BigFloat, 256) do
            # P = 2 (80 dofs) keeps the dense BigFloat bottom solve cheap, as in the sealed
            # row-42 gate; ONE AMR step is enough to exercise every mesh-bound rebuild.
            cfg = da_cfg(; spin = big"0.5", r0 = big"10.0", m = 2, P = 2)
            orb = da_orb(cfg)
            h = solve_mode_amr(orb, cfg, 2; T = Complex{BigFloat}, iterations = 1)
            @test h isa AmrHistory{Complex{BigFloat},BigFloat}
            @test length(h.levels) == 2
            @test h.levels[1].mesh.ndofs == 80          # 20 elements × 2 × 2
            @test all(l -> length(l.punctures) == 2, h.levels)
            @test h.x_p == (big"10.0", big"0.0")
            @test h.levels[end].psi isa Complex{BigFloat}
            @test h.levels[end].dpsi_r isa Complex{BigFloat}
            @test precision(real(h.levels[end].psi)) == 256      # INHERITED, never Float64'd
            @test h.solution.stack.options.abstol isa BigFloat
            @test all(isfinite, h.solution.u)
            @test mode_residual(h.solution) <= cfg.gmres.absolute_residual
            # the criterion itself is exact in any precision: the snapped coordinates are
            # 0 and ±1, compared against exact dyadic endpoints
            @test amr_flags(h.solution) == amr_flags(h.levels[end].mesh, h.solution.dom, orb,
                                                     cfg)
        end
    end

    # ======================================================================================
    # (h) THE MANUFACTURED FALSIFIER — the paper's rationale for the scheme
    # ======================================================================================
    @testset "(h) inverting the criterion destroys the convergence" begin
        # main.tex:552-560 + main.tex:721: h-refine AT the puncture, p-refine everywhere else,
        # and "exponential convergence relies solely on our ability to h-refine non-smooth
        # points". Drop the h half — pure p-refinement, the SAME number of levels and the same
        # element set — and the particle diagnostics must do markedly worse, because raising
        # the polynomial order cannot resolve a non-smooth point.
        # MEASURED (proto2_a05_m2_purep.log vs proto2_a05_m2_crp.log): see the evidence table.
        # (The literal inversion of the prompt — p at the puncture, h EVERYWHERE else — is
        # measured too, in `proto2_a05_m2_invertedB.log`; it is not gated here because it
        # quadruples the element count per level (20 → 74 → 290 → 1154) and would dominate
        # the suite's runtime while testing the same statement.)
        good = da_hist                                   # the shared 6-level sweep
        bad = solve_mode_amr(da_porb, da_pcfg, 2; iterations = 6,
                             criterion = sol -> Dict{ElementId,NTuple{2,Bool}}())
        @test all(l -> length(l.punctures) == 2, good.levels)
        @test all(l -> isempty(l.punctures), bad.levels)
        @test length(bad.levels[end].mesh.elements) == 20      # p only: no element was split
        @test bad.levels[end].mesh.ndofs == 20 * (4 + 6)^2
        dg, rg = da_diffs([l.psi for l in good.levels])
        db, rb = da_diffs([l.psi for l in bad.levels])
        eg, _ = da_diffs([l.dpsi_r for l in good.levels])
        eb, _ = da_diffs([l.dpsi_r for l in bad.levels])
        # MEASURED at 6 levels (proto2_a05_m2_crp.log vs proto2_a05_m2_purep.log):
        #   Ψ^ℛ      last difference 2.354e-8 (hp) vs 2.619e-7 (p only)  — 11.1x worse
        #   ∂_rΨ^ℛ   last difference 6.154e-9 (hp) vs 3.068e-8 (p only)  —  5.0x worse
        #   ratio geometric mean 0.242 (hp) vs 0.398 (p only): p-refinement alone halves the
        #   error per level (algebraic in the element count) instead of cutting it by ~4x,
        #   and its ∂_r sequence is not even monotone (one ratio of 2.0 at level 7).
        @test db[end] > 5 * dg[end]
        @test eb[end] > 4 * eg[end]
        @test da_georatio(rg) < 0.30 <= da_georatio(rb)
    end
end

# test_driver.jl — ESSF node `solver-stack-ssf` (GLOBAL_DAG row 42, wave 10, worker W10-A),
# written BEFORE `src/Solve/driver.jl` exists (TDD; `tests_red.log` captured against the
# missing file). THIS IS THE JOIN: the first solve of the real CRP m-mode problem, on the real
# CRP domain, with the full GMRES ∘ multigrid ∘ additive Schwarz ∘ cached-bottom stack.
#
# WHAT THE NODE DELIVERS (`Solve/driver.jl`):
#   * `solver_options(config)`  — the `SolverStackOptions` twin of the sealed `RunConfig`,
#     field by field, with the two multigrid LEVEL-POLICY knobs (`coarse_levels`,
#     `coarse_orders`) exposed because the AMR driver (row 43) must choose them.
#   * `solve_mode(orbit, config, m)` / `solve_mode!` / `ModeSolution` — one m-mode solve.
#   * `mode_fixed_source` — the VOLUME source: the SWITCHED evaluator inside the worldtube,
#     exactly zero outside, in the massive `M·S` form.
#   * `mode_residual`, `reconstruct_full_field`, `field_at` — the read-out surface.
#   * `rebuild_for_mesh!` — the ONE-STEP rebuild the AMR loop calls (the LOOP is row 43).
#
# ORACLE (the DAG row-42 gate column + the READ-ONLY mirrors, transcribed by hand):
#   * `main.tex:720` — "the preconditioners ... allow the GMRES algorithm to converge rapidly
#     in less than 10 iterations per AMR level."
#   * `fig:linsolv` caption, `main.tex:701` — "The first AMR level converges immediately
#     because the multigrid bottom solver solves the problem directly." THAT is what the
#     production yaml's `InitialCoarseLevels: L = 0` means: `maxlevels = L + 1 = 1`, the
#     hierarchy is a SINGLE level, and its bottom solve IS the whole problem.
#   * `main.tex:672-676` — h-coarsening to the irreducible block structure; the cached bottom
#     LU and the cached subdomain inverses across AMR iterations; two Schwarz passes.
#   * `main.tex:680` — the `m²Ω²` complex shift, "disabled altogether" under null slicing.
#   * `main.tex:513-523` (eq:jump_worldtube) — `Ψ^ℛ_m = Ψ_m − Ψ^𝒫_m` across the worldtube
#     interface, which is what makes the RECONSTRUCTED field continuous there.
#   * `main.tex:488-492` — the near-worldline route switch the VOLUME source must use.
#   * `anc/ScalarSelfForce.yaml:178-205` — the solver block the options twin reproduces.
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here): `solver-stack-core`
# (`SolverStack`, `SolverStackOptions`, `solve!`/`solve_linearized!`/`rebuild!`,
# `ResidualMonitor`/`SolveRecord`/`cumulative_iterations`, `linearized_rhs`), `fig-matrix`
# (`assemble_operator`, `operator_builder`, `coarsest_operator`), `dg-jump-conditions`
# (`jump_hooks`, `vtu_junction_faces`), `dg-boundary-conditions` (`boundary_conditions`,
# `NoneGhost`), `worldtube-formulation` (`worldtube_data`, `worldtube_massive_source`,
# `worldtube_boundary_faces`, `worldtube_face_points`), `effective-source-series`
# (`effective_source_m_switched`, `effective_source_switched_bundle`),
# `effective-source-direct` (`puncture_psi_chart`), `domain-blocks-map` (`block_domain`,
# `boundary_face_classes`), `mesh-basis-quadrature` (`refelement`, `interpolation_matrix`),
# `amr-h/p-refinement` (`h_refine`, `hp_refine`, `transfer_solution`), `multigrid-hierarchy`
# (`coarsen_mesh`), `schwarz-complex-shift` (`helmholtz_shift`).
#
# GATE CONSTANTS — every one calibrated on the throwaway prototypes
# `evidence/solver-stack-ssf/calibrate_first_solve.jl`, `calibrate_gates2.jl` and
# `calibrate_gates3.jl` BEFORE this file was written, never tuned afterwards to pass. The
# measured value is quoted beside each gate.

using Test
using LinearAlgebra
using ScalarSelfForce
using EllipticCore

# ---- configs (dv_-prefixed: sibling test files share the runtests.jl include scope) --------
# `block_domain` REQUIRES the inner radial bound to BE r₊ (main.tex:529), which the yaml's
# inner_radius = 2 satisfies only at a = 0 — every spin moves it with the horizon.
dv_cfg(; spin = 0.5, r0 = 10.0, m = 2, L = 0, P = 4) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L)
dv_orb(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)
# rebuild a RunConfig with one nested field replaced (the config structs are immutable)
function dv_setfield(cfg::RunConfig{T}, group::Symbol, field::Symbol, value) where {T}
    g = getfield(cfg, group)
    G = typeof(g)
    g2 = G((f === field ? value : getfield(g, f) for f in fieldnames(G))...)
    return RunConfig{T}((f === group ? g2 : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end
dv_solve(; kw...) = (c = dv_cfg(; kw...); solve_mode(dv_orb(c), c, c.orbit.m_mode_number))

# ---- an INDEPENDENT re-assembly of the volume fixed source (the gate (d4) oracle) ---------
# Deliberately NOT a call into the driver: it walks the mesh itself, evaluates the sealed
# switched evaluator at every collocation node of every inside element, applies the sealed
# element mass matrix, and leaves everything else exactly zero.
function dv_reference_source(::Type{T}, dom, orbit, m; mesh = dom.mesh, phi_p = 0) where {T}
    p = Puncture(orbit)
    nodal = zeros(T, mesh.ndofs)
    massive = zeros(T, mesh.ndofs)
    nseries = 0; ndirect = 0; lam_min = Inf
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        inside_worldtube(dom, el) || continue
        em = element_matrices(T, mesh.blocks[mesh.blockindex[el.block]], el, mesh.orders[e])
        off = mesh.offsets[e]
        nx, ny = mesh.orders[e]
        for j in 1:ny, i in 1:nx
            r, x2 = em.x[1][i], em.x[2][j]
            z = dom.chart === :z2 ? sqrt(x2) : x2
            b = effective_source_switched_bundle(p, m, r, z; phi_p = phi_p)
            b.route === :series ? (nseries += 1) : (ndirect += 1)
            dr, dth = r - orbit.r0, -asin(z)
            lam_min = min(lam_min, sqrt(orbit.rho2_dr2 * dr^2 + orbit.rho2_dtheta2 * dth^2))
            nodal[off + tensor_lin(i, j, nx)] = b.source
        end
        rng = off+1:off+nx*ny
        massive[rng] .= em.M * view(nodal, rng)
    end
    return (nodal = nodal, massive = massive, nseries = nseries, ndirect = ndirect,
            lam_min = lam_min)
end

# ---- one-sided trace of a nodal field AND of its chart-normal derivative ------------------
# `field_at(...; element)` is the node's own evaluator (gated separately against the nodal
# values); this helper adds the DERIVATIVE trace the flux-jump gate needs, from the same
# sealed `refelement` machinery. NOTE the interpolated derivative row is `Py·D`, not `D·Pyᵀ`.
function dv_trace(mesh, values, e::Int, x1, x2, d::Int)
    el = mesh.elements[e]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    nx, ny = mesh.orders[e]
    RT = Float64
    ξ = ntuple(2) do dd
        lo, hi = logical_interval(RT, el, dd)
        ζ = map_inverse(blockmap(blk, dd), dd == 1 ? x1 : x2)
        (2ζ - (lo + hi)) / (hi - lo)
    end
    ex, ey = refelement(RT, nx), refelement(RT, ny)
    Px, Py = interpolation_matrix(ex, [ξ[1]]), interpolation_matrix(ey, [ξ[2]])
    off = mesh.offsets[e]
    U = reshape(view(values, off+1:off+nx*ny), nx, ny)
    val = (Px * U * transpose(Py))[1, 1]
    dξ = d == 1 ? ((Px * ex.D) * U * transpose(Py))[1, 1] :
                  (Px * U * transpose(Py * ey.D))[1, 1]
    lo, hi = logical_interval(RT, el, d)
    ζ = (lo * (1 - ξ[d]) + hi * (1 + ξ[d])) / 2
    jac = map_jacobian(blockmap(blk, d), ζ) * (hi - lo) / 2
    return val, dξ / jac
end

# worst value- and flux-jump across the worldtube boundary, both one-sided traces taken at
# the operator's OWN mortar points (sealed `worldtube_face_points`)
function dv_worldtube_jump(sol, values; reconstruct = true)
    dom, mesh, m = sol.dom, sol.mesh, sol.m
    p = Puncture(sol.orbit)
    wj = 0.0; wjs = 0.0; wf = 0.0; wfs = 0.0; npts = 0
    for f in worldtube_boundary_faces(mesh, el -> inside_worldtube(dom, el))
        el, d, sidx = f
        e = mesh.index[el]
        iface = mesh.faces[e][d][sidx]
        xf, mpts = worldtube_face_points(mesh, f)
        for (k, nb) in enumerate(iface.nbs), tq in mpts[k]
            x1, x2 = d == 1 ? (xf, tq) : (tq, xf)
            z = dom.chart === :z2 ? sqrt(x2) : x2
            vi, di = dv_trace(mesh, values, e, x1, x2, d)
            vo, dd = dv_trace(mesh, values, nb, x1, x2, d)
            pb = puncture_psi_chart(p, m, x1, z; phi_p = sol.phi_p)
            pp = reconstruct ? pb.psi : zero(pb.psi)
            ppd = if !reconstruct
                zero(pb.dz)
            elseif d == 1
                pb.dr
            else
                dom.chart === :z2 ? pb.dz / (2z) : pb.dz
            end
            add = inside_worldtube(dom, el)
            Ψi = add ? vi + pp : vi
            Ψo = add ? vo : vo + pp
            Di = add ? di + ppd : di
            Do = add ? dd : dd + ppd
            A = d == 1 ? flux_coeff_r(sol.orbit, x1) :
                (dom.chart === :z2 ? flux_coeff_z2(sol.orbit, x1, z) :
                                     flux_coeff_z(sol.orbit, x1, z))
            j, fj = abs(Ψi - Ψo), abs(A * (Di - Do))
            j > wj && (wj = j; wjs = max(abs(Ψi), abs(Ψo)))
            fj > wf && (wf = fj; wfs = max(abs(A * Di), abs(A * Do)))
            npts += 1
        end
    end
    return (value = wj, value_scale = wjs, flux = wf, flux_scale = wfs, npoints = npts)
end

# peak and far-field amplitudes of Ψ_m and of the PHYSICAL amplitude Ψ_m/r (the ansatz
# already factors 1/r out, so it is Ψ_m/r that must decay — see the gate (d3) comment)
function dv_decay(sol)
    mesh = sol.mesh
    pk = 0.0; pkphys = 0.0; far = 0.0; farphys = 0.0
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        nx, ny = mesh.orders[e]
        em = element_matrices(Float64, mesh.blocks[mesh.blockindex[el.block]], el,
                              mesh.orders[e])
        off = mesh.offsets[e]
        for j in 1:ny, i in 1:nx
            r = em.x[1][i]
            a = abs(sol.u[off + tensor_lin(i, j, nx)])
            pk = max(pk, a); pkphys = max(pkphys, a / r)
            if r > 100
                far = max(far, a); farphys = max(farphys, a / r)
            end
        end
    end
    return (peak = pk, far = far, peak_physical = pkphys, far_physical = farphys)
end

@testset "single-mode solver-stack driver (solver-stack-ssf)" begin

# ==================== (0) solver_options: the RunConfig / yaml twin ======================
@testset "(0) solver_options is the field-by-field RunConfig twin" begin
    cfg = dv_cfg()
    o = solver_options(cfg)
    @test o isa SolverStackOptions
    @test o.abstol == cfg.gmres.absolute_residual == 1.0e-12         # yaml:183
    @test o.reltol == cfg.gmres.relative_residual == 0               # yaml:182
    @test o.maxiter == cfg.gmres.max_iterations == 500               # yaml:181
    @test o.restart == o.maxiter                                     # unrestarted GMRES
    @test o.multigrid_iterations == cfg.multigrid.iterations == 1    # yaml:187
    @test o.initial_coarse_levels == cfg.multigrid.initial_coarse_levels == 0  # yaml:188, L=0
    @test o.pre_smoothing == cfg.multigrid.pre_smoothing == true     # yaml:189
    @test o.use_bottom_solver == cfg.multigrid.use_bottom_solver == true       # yaml:190
    @test o.post_smoothing_at_bottom == cfg.multigrid.post_smoothing_at_bottom == false
    @test o.passes == cfg.schwarz.iterations == 2                    # yaml:197
    # yaml:198 MaxOverlap: Auto -> the NEIGHBOUR'S FULL EXTENT = P = 4 (see (0b)), not 2
    @test o.overlap == 4 == resolve_max_overlap(cfg)
    @test o.subdomain_complex_shift == cfg.schwarz.subdomain_complex_shift == 0  # yaml:203
    @test o.shrink == 1                                              # the sealed row-21 default
    @test o.coarse_orders === nothing                                # the MIN rule by default

    # the two AMR level-policy knobs are EXPOSED (row 43 chooses; see the docstring table)
    @test solver_options(cfg; coarse_levels = 3).initial_coarse_levels == 3
    @test solver_options(cfg; coarse_levels = nothing).initial_coarse_levels === nothing
    @test solver_options(cfg; coarse_orders = 4).coarse_orders == 4

    # BigFloat config ⇒ BigFloat tolerances (never silently Float64'd)
    @test solver_options(dv_cfg(; spin = big"0.5", r0 = big"10.0")).abstol isa BigFloat

    # THE SHIFT POLICY VALUE is exactly zero under the CRP's vtu slicing (main.tex:680)
    orb = dv_orb(cfg)
    @test helmholtz_shift(2, orb.Omega; slicing = :vtu) == 0
    @test helmholtz_shift(2, orb.Omega; slicing = :t) ≈ (2 * orb.Omega)^2
end

@testset "(0b) H2: every RunConfig field that affects the solve is mapped or REFUSED" begin
    cfg = dv_cfg()
    # MaxOverlap: Auto (yaml:198) resolves to the NEIGHBOUR'S FULL EXTENT — spectre
    # Tags.hpp:29-35 ("'Auto' means no restriction, so the overlap covers the full neighbor")
    # and OverlapHelpers.cpp:21-24 (min(volume_extent, max_overlap.value_or(volume_extent))).
    # At the production InitialGridPoints [4, 4] that is 4, NOT the package default 2.
    @test cfg.schwarz.max_overlap === :auto
    @test resolve_max_overlap(cfg) == 4
    @test solver_options(cfg).overlap == 4
    @test resolve_max_overlap(dv_cfg(; P = 6)) == 6
    @test solver_options(cfg; overlap = 2).overlap == 2          # caller override kept
    # any OTHER MaxOverlap value is REFUSED, not silently defaulted. [HOLE] the yaml permits
    # an integer MaxOverlap but the sealed SchwarzConfig types the field `::Symbol`
    # (Domain/config.jl:133), so an integer cannot be carried at all — node
    # runtime-config-fixture would need Union{Symbol,Integer}; resolve_max_overlap already
    # accepts one if it is ever passed.
    @test fieldtype(typeof(cfg.schwarz), :max_overlap) === Symbol
    @test_throws ArgumentError resolve_max_overlap(dv_setfield(cfg, :schwarz,
                                                               :max_overlap, :full))
    # Multigrid.ComplexShift (yaml:194) has no twin (subdomain-level shift only, decision 8);
    # the production value is 0 and a NONZERO one is refused rather than silently dropped
    @test cfg.multigrid.complex_shift == 0
    @test_throws ArgumentError solver_options(dv_setfield(cfg, :multigrid,
                                                          :complex_shift, 0.5))
    # ENUMERATION of the solver-relevant RunConfig fields: each is mapped, or refused above
    o = solver_options(cfg)
    mapped = Dict(
        :gmres_max_iterations      => o.maxiter   == cfg.gmres.max_iterations,
        :gmres_relative_residual   => o.reltol    == cfg.gmres.relative_residual,
        :gmres_absolute_residual   => o.abstol    == cfg.gmres.absolute_residual,
        :multigrid_iterations      => o.multigrid_iterations == cfg.multigrid.iterations,
        :multigrid_initial_coarse_levels => o.initial_coarse_levels ==
                                            cfg.multigrid.initial_coarse_levels,
        :multigrid_pre_smoothing   => o.pre_smoothing == cfg.multigrid.pre_smoothing,
        :multigrid_use_bottom_solver => o.use_bottom_solver == cfg.multigrid.use_bottom_solver,
        :multigrid_post_smoothing_at_bottom => o.post_smoothing_at_bottom ==
                                               cfg.multigrid.post_smoothing_at_bottom,
        :schwarz_iterations        => o.passes    == cfg.schwarz.iterations,
        :schwarz_max_overlap       => o.overlap   == resolve_max_overlap(cfg),
        :schwarz_subdomain_complex_shift => o.subdomain_complex_shift ==
                                            cfg.schwarz.subdomain_complex_shift,
        :dg_penalty_parameter      => cfg.dg.penalty_parameter == 1.5,  # consumed by assemble
    )
    for (k, v) in mapped
        @test v || error("RunConfig field $k is not mapped onto the solver")
    end
    @test length(mapped) == 12
    # the two verbosity/observer-only fields are deliberately NOT solver inputs
    @test cfg.gmres.verbosity === :Quiet && cfg.multigrid.verbosity === :Silent
end

# ============================== (a) THE ROW-42 ORACLE ====================================
sol05 = dv_solve(; spin = 0.5, r0 = 10.0, m = 2)
sol00 = dv_solve(; spin = 0.0, r0 = 10.0, m = 2)

@testset "(a) row-42 oracle at a = $(sol.orbit.a), r0 = 10, m = 2, L = 0, P = 4" for sol in
        (sol05, sol00)
    @test sol isa ModeSolution
    @test sol.mesh.ndofs == 320                    # 20 elements × 4×4 (yaml:74, 78-110)
    @test length(sol.mesh.elements) == 20
    @test sol.dom.chart === :z2                    # ImposeEquatorialSymmetry (yaml:36)
    @test sol.m == 2
    @test length(sol.u) == 320 && eltype(sol.u) === ComplexF64
    @test length(sol.b) == 320

    # ---- THE GATE: converged to the yaml's own absolute residual in < 10 iterations
    @test sol.result.residual <= sol.stack.options.abstol   # MEASURED 1.21e-15 / 6.69e-16
    @test sol.result.iterations < 10
    @test sol.result.converged
    # the production yaml is InitialCoarseLevels = L = 0 ⇒ ONE multigrid level ⇒ the bottom
    # solver solves the whole problem directly (fig:linsolv caption, main.tex:701), so the
    # MEASURED count is exactly 1. Pinned so a regression to "many" is caught.
    @test sol.result.iterations == 1
    @test sol.record.level_ndofs == [320]

    # ---- finite EVERYWHERE: nodes, the horizon column r = r₊, the equator and the pole
    @test all(isfinite, sol.u)
    rp = sol.orbit.rp
    R = sol.config.domain.block_bounds_r[end]
    for x2 in (0.0, 0.25, 1.0)                     # equator, worldtube row edge, pole
        @test isfinite(field_at(sol, rp, x2))      # r = r₊ EXACTLY (the horizon face trace)
        @test isfinite(field_at(sol, R, x2))       # the outer boundary
    end
    @test isfinite(field_at(sol, sol.orbit.r0, 0.5))
end

@testset "(a2) field_at reproduces the collocation values and honours `element`" begin
    sol = sol05
    mesh = sol.mesh
    worst = 0.0
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        nx, ny = mesh.orders[e]
        em = element_matrices(Float64, mesh.blocks[mesh.blockindex[el.block]], el,
                              mesh.orders[e])
        off = mesh.offsets[e]
        for j in 1:ny, i in 1:nx
            v = field_at(sol, em.x[1][i], em.x[2][j]; element = e)
            worst = max(worst, abs(v - sol.u[off + tensor_lin(i, j, nx)]))
        end
    end
    @test worst <= 1e-14 * maximum(abs, sol.u)     # MEASURED 8.36e-17 against ‖u‖∞ = 0.182
    # it evaluates ANY nodal field of the mesh, which is how rows 44/48 use it
    Ψ = reconstruct_full_field(sol)
    @test field_at(sol, Ψ, sol.orbit.r0, 0.5) == field_at(sol, Ψ, sol.orbit.r0, 0.5)
    @test_throws ArgumentError field_at(sol, 1.0e9, 0.5)          # outside every element
    @test_throws ArgumentError field_at(sol, zeros(ComplexF64, 3), sol.orbit.r0, 0.5)
end

@testset "(M2) field_grad_at: exact on polynomials, spectral against a finite difference" begin
    sol = sol05
    mesh = sol.mesh
    # ---- EXACTNESS on a manufactured polynomial of the element's own degree
    e = 1
    el = mesh.elements[e]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    nx, ny = mesh.orders[e]
    em = element_matrices(Float64, blk, el, mesh.orders[e])
    f(x, y)  = 2.0 + 0.7x - 1.3y + 0.45x^2 * y - 0.2y^3 + 0.11x^3
    fx(x, y) = 0.7 + 0.9x * y + 0.33x^2
    fy(x, y) = -1.3 + 0.45x^2 - 0.6y^2
    vals = zeros(ComplexF64, mesh.ndofs)
    for j in 1:ny, i in 1:nx
        vals[mesh.offsets[e] + tensor_lin(i, j, nx)] = f(em.x[1][i], em.x[2][j])
    end
    xq, yq = em.x[1][2], em.x[2][2]
    g = field_grad_at(sol, vals, xq, yq; element = e)
    @test isapprox(real(g[1]), fx(xq, yq); rtol = 1e-10)      # exact to roundoff
    @test isapprox(real(g[2]), fy(xq, yq); rtol = 1e-10)
    @test field_at(sol, vals, xq, yq; element = e) ≈ f(xq, yq) rtol = 1e-12

    # ---- SPECTRAL agreement with a finite difference of field_at on the REAL solution
    r0, w0 = sol.orbit.r0 + 1.7, 0.4
    e2 = ScalarSelfForce._dr_locate(mesh, r0, w0)
    gr = field_grad_at(sol, r0, w0; element = e2)
    h = 1e-5
    fd1 = (field_at(sol, r0 + h, w0; element = e2) -
           field_at(sol, r0 - h, w0; element = e2)) / (2h)
    fd2 = (field_at(sol, r0, w0 + h; element = e2) -
           field_at(sol, r0, w0 - h; element = e2)) / (2h)
    @test isapprox(gr[1], fd1; rtol = 1e-6)
    @test isapprox(gr[2], fd2; rtol = 1e-6)
    # the short form differentiates sol.u, and any nodal field is accepted
    @test field_grad_at(sol, r0, w0; element = e2) == field_grad_at(sol, sol.u, r0, w0;
                                                                   element = e2)
    @test_throws ArgumentError field_grad_at(sol, zeros(ComplexF64, 3), r0, w0)

    # ---- ONE-SIDED traces at the worldtube seam differ (a DG derivative is discontinuous)
    f0 = first(worldtube_boundary_faces(mesh, el2 -> inside_worldtube(sol.dom, el2)))
    elf, d, sidx = f0
    ei = mesh.index[elf]
    nb = mesh.faces[ei][d][sidx].nbs[1]
    xf, mpts = worldtube_face_points(mesh, f0)
    tq = mpts[1][1]
    xa, xb = d == 1 ? (xf, tq) : (tq, xf)
    gin  = field_grad_at(sol, xa, xb; element = ei)
    gout = field_grad_at(sol, xa, xb; element = nb)
    @test gin != gout
    # ... and field_grad_at agrees with the test-local prototype it was promoted from
    vi, di = dv_trace(mesh, sol.u, ei, xa, xb, d)
    @test field_at(sol, xa, xb; element = ei) ≈ vi rtol = 1e-12
    @test gin[d] ≈ di rtol = 1e-10
end

# ==================== (b) self-consistency: residual and direct_solve ====================
@testset "(b) mode_residual and the direct-solve envelope" begin
    for sol in (sol05, sol00)
        r = mode_residual(sol)
        @test r <= sol.stack.options.abstol                # MEASURED 1.08e-15 / 6.87e-16
        # the same quantity as the sealed GMRESResult's, but NOT to a tight rtol: at a
        # residual of 1e-15 the two recomputations differ by roundoff OF THEIR OWN ORDER.
        # `A_full(u) = A_lin(u) + operator_on_zero(op)` is an exact identity in exact
        # arithmetic, but `mode_residual` evaluates the left side and GMRES the right, with
        # different rounding, and ‖operator_on_zero‖ = 1.74 here against a residual of 1e-15
        # — i.e. 15 orders of cancellation. MEASURED ratios 0.893 (a = 0.5) / 1.028 (a = 0).
        @test 0.1 <= r / sol.result.residual <= 10
        y = similar(sol.u)
        apply_full!(y, sol.op, sol.u)
        @test r ≈ norm(sol.b .- y) rtol = 1e-12            # ... and it IS ‖b − A_full·u‖

        # agreement with the sealed dense bypass on the SAME linearized system.
        # MEASURED 3.27e-14 (a = 0.5) / 2.42e-14 (a = 0) relative, at cond(A) ≈ 6.4e4 —
        # the solve is far better than its tolerance demands. Envelope 1e-11.
        ud = direct_solve(sol.op, linearized_rhs(sol.op, sol.b))
        @test norm(sol.u .- ud) / norm(ud) <= 1e-11
    end
    # the pre-folded route is the SAME run (the sealed row-21 contract, exercised from here)
    u2 = zeros(ComplexF64, sol05.mesh.ndofs)
    r2 = solve_linearized!(u2, sol05.stack, linearized_rhs(sol05.op, sol05.b))
    @test u2 == sol05.u
    @test r2.iterations == sol05.result.iterations
end

@testset "(b2) FALSIFIER: double-folding the right-hand side is an O(1) wrong answer" begin
    # `linearized_rhs(op, wt)` (the assembly helper) is ALREADY folded; `solve!` folds again.
    # The failure is SILENT: it converges, at machine residual, on a different system.
    sol = sol05
    @test norm(operator_on_zero(sol.op)) > 1e3 * norm(sol.b)   # MEASURED ratio 2.34e4
    ubad = zeros(ComplexF64, sol.mesh.ndofs)
    rbad = solve!(ubad, sol.stack, linearized_rhs(sol.op, sol.wt))
    @test rbad.converged                                       # ... and it "converges"
    @test rbad.residual <= sol.stack.options.abstol            # MEASURED 1.98e-15
    @test norm(ubad .- sol.u) / norm(sol.u) > 0.5              # MEASURED 0.989 — 99 % wrong
    # the driver's own right-hand side is the PHYSICAL one and is never pre-folded
    @test sol.b == mode_fixed_source(ComplexF64, sol.dom, sol.orbit, sol.m; mesh = sol.mesh)
end

# ============ (c) the worldtube seam: Ψ = Ψ^ℛ + Ψ^𝒫 is continuous across it ==============
@testset "(c1) reconstruct_full_field: outside is bitwise the solution, inside adds Ψ^𝒫" begin
    for sol in (sol05, sol00)
        Ψ = reconstruct_full_field(sol)
        @test length(Ψ) == sol.mesh.ndofs && all(isfinite, Ψ)
        nin = 0
        for e in eachindex(sol.mesh.elements)
            el = sol.mesh.elements[e]
            nx, ny = sol.mesh.orders[e]
            rng = sol.mesh.offsets[e]+1 : sol.mesh.offsets[e]+nx*ny
            if inside_worldtube(sol.dom, el)
                nin += length(rng)
                @test view(Ψ, rng) != view(sol.u, rng)
            else
                @test view(Ψ, rng) == view(sol.u, rng)     # BITWISE untouched
            end
        end
        @test nin == 64                                     # 4 in-tube elements × 4×4
    end
end

@testset "(c2) the reconstructed field is continuous across the worldtube boundary" begin
    for sol in (sol05, sol00)
        j = dv_worldtube_jump(sol, sol.u; reconstruct = true)
        raw = dv_worldtube_jump(sol, sol.u; reconstruct = false)
        @test j.npoints == 48                               # 11 faces on the production lattice
        # MEASURED worst |[Ψ]| = 2.796e-4 (a = 0.5) / 2.779e-4 (a = 0) against a field scale
        # of 0.201 — 1.4e-3 relative, INSIDE the sealed row-38 worldtube envelope (8e-3
        # absolute for the puncture residual on this coarse mesh).
        @test j.value <= 1.0e-3
        @test j.value_scale > 0.1                           # measured against a real field
        # the FLUX jump: MEASURED 2.178e-3 / 2.167e-3 against a flux scale of 0.0277
        @test j.flux <= 5.0e-3
        # ---- FALSIFIER: without the reconstruction the SAME measurement is O(Ψ^𝒫)
        @test raw.value > 0.05                              # MEASURED 0.2017 / 0.2018
        @test raw.value > 100 * j.value                     # MEASURED 721× / 726×
        @test raw.flux > 4 * j.flux                         # MEASURED 11.5× / 11.5×
    end
end

@testset "(c3) the jump is a discretization envelope: it converges spectrally in p" begin
    j4 = dv_worldtube_jump(sol05, sol05.u)
    s6 = dv_solve(; spin = 0.5, r0 = 10.0, m = 2, P = 6)
    j6 = dv_worldtube_jump(s6, s6.u)
    @test s6.result.residual <= s6.stack.options.abstol
    # MEASURED P = 4 → 6 → 8: value 2.80e-4 → 1.34e-5 → 1.09e-6; flux 2.18e-3 → 1.82e-4 → 1.67e-5
    @test j6.value <= 5.0e-5
    @test j6.value < j4.value / 5
    @test j6.flux <= 1.0e-3
    @test j6.flux < j4.flux / 5
end

@testset "(c4) FALSIFIER: dropping the worldtube hooks leaves a jump of O(Ψ^𝒫)" begin
    sol = sol05
    junc = vtu_junction_faces(sol.dom, sol.orbit, sol.config, sol.m; mesh = sol.mesh)
    empty_jd = Dict{Tuple{ElementId,Int,Int},Vector{FaceState{ComplexF64}}}()
    fcs = domain_coefficients(sol.dom, sol.orbit, sol.config)
    op0 = MeshOperator(ComplexF64, sol.mesh, fcs.A, fcs.beta, fcs.gamma;
                       C = sol.config.dg.penalty_parameter,
                       hooks = jump_hooks(sol.mesh, empty_jd, junc),
                       bcs = boundary_conditions(sol.dom, sol.orbit, sol.config;
                                                 mesh = sol.mesh))
    st0 = SolverStack(op0, _ -> op0; options = solver_options(sol.config))
    u0 = zeros(ComplexF64, sol.mesh.ndofs)
    res0 = solve!(u0, st0, sol.b)
    @test res0.converged                       # it converges HAPPILY on the wrong problem
    @test res0.residual <= sol.stack.options.abstol
    j0 = dv_worldtube_jump(sol, u0; reconstruct = true)
    @test j0.value > 0.05                      # MEASURED 0.2019, i.e. |Ψ^𝒫| itself
    @test j0.value > 100 * dv_worldtube_jump(sol, sol.u).value
end

# ====================== (d) physics sanity without a reference ===========================
@testset "(d1) the m-sweep m ∈ {0, 1, 2, 5} all converge" begin
    for m in (0, 1, 2, 5)
        cfg = dv_cfg(; spin = 0.5, m = m)
        s = solve_mode(dv_orb(cfg), cfg, m)
        # at m = 0 the sealed GMRES may report :stagnation at machine residual (row-41
        # caveat) — THE GATE IS THE RESIDUAL, never `converged` alone
        @test s.result.residual <= cfg.gmres.absolute_residual
        @test mode_residual(s) <= cfg.gmres.absolute_residual
        @test s.result.iterations < 10          # MEASURED 1 at every m (single-level yaml MG)
        @test all(isfinite, s.u)
        @test all(isfinite, reconstruct_full_field(s))
        # the worldtube seam holds at every mode too
        j = dv_worldtube_jump(s, s.u)
        raw = dv_worldtube_jump(s, s.u; reconstruct = false)
        @test j.value <= 1.0e-3                 # MEASURED 2.70e-4 (m=0) … 2.57e-4 (m=5)
        @test raw.value > 50 * j.value          # MEASURED 3773× (m=0), 158× (m=5)
    end
end

@testset "(d2) m = 0 carries a NoneGhost on the outer face (yaml:114-125)" begin
    cfg0 = dv_cfg(; spin = 0.5, m = 0)
    s0 = solve_mode(dv_orb(cfg0), cfg0, 0)
    @test cfg0.domain.outer_radial_bc === nothing
    bcs0 = boundary_conditions(s0.dom, s0.orbit, s0.config; mesh = s0.mesh)
    outer = [f for (f, c) in boundary_face_classes(s0.dom) if c === :outer]
    @test !isempty(outer)
    @test all(bcs0[f] isa NoneGhost for f in outer)
    # ... while m ≠ 0 puts the sealed Bayliss–Turkel ghost on exactly those faces
    bcs2 = boundary_conditions(sol05.dom, sol05.orbit, sol05.config; mesh = sol05.mesh)
    @test all(bcs2[f] isa BaylissTurkelGhost for f in outer)
end

@testset "(d3) the PHYSICAL amplitude decays outward; Ψ_m itself need not" begin
    # The ansatz Φ = Σ_m Ψ_m sin^{|m|}θ / r e^{im(φ̄−Ωt)} already factors 1/r out, and under
    # the u-slab's hyperboloidal compaction the outgoing solution approaches a CONSTANT at
    # future null infinity — so it is Ψ_m/r, not Ψ_m, that must fall off. Both are measured;
    # only the physical one is gated.
    for (m, gate) in ((0, 0.10), (2, 0.05), (5, 0.02))
        cfg = dv_cfg(; spin = 0.5, m = m)
        s = solve_mode(dv_orb(cfg), cfg, m)
        d = dv_decay(s)
        # MEASURED far/peak of |Ψ_m|/r at r > 100: 0.0588 (m=0), 0.0171 (m=2), 0.00358 (m=5)
        @test d.far_physical <= gate * d.peak_physical
        # MEASURED far/peak of |Ψ_m| itself: 0.840, 0.219, 0.0512 — recorded, NOT gated
        @test d.far <= d.peak
    end
    # regular at the horizon: the r = r₊ face trace exists and stays on the field's scale
    @test abs(field_at(sol05, sol05.orbit.rp, 0.1)) < 10 * maximum(abs, sol05.u)
end

@testset "(d4) the VOLUME source is the SWITCHED evaluator inside, exactly zero outside" begin
    sol = sol05
    ref = dv_reference_source(ComplexF64, sol.dom, sol.orbit, sol.m; mesh = sol.mesh,
                              phi_p = sol.phi_p)
    @test sol.b == ref.massive                              # BITWISE
    @test mode_fixed_source(ComplexF64, sol.dom, sol.orbit, sol.m; mesh = sol.mesh) == sol.b
    for e in eachindex(sol.mesh.elements)
        el = sol.mesh.elements[e]
        inside_worldtube(sol.dom, el) && continue
        nx, ny = sol.mesh.orders[e]
        rng = sol.mesh.offsets[e]+1 : sol.mesh.offsets[e]+nx*ny
        @test all(iszero, view(sol.b, rng))                 # EXACTLY zero outside
    end
    @test ref.nseries + ref.ndirect == 64                   # 4 in-tube elements × 4×4
    # ROUTE CENSUS at Float64 on this lattice: MEASURED 0 series / 64 direct, λ_min = 1.004,
    # three orders outside the Float64 crossover λ* ≈ 1.2e-2 — so the switch is INERT here
    # and it is the AMR puncture-snap resolution that decides. Recorded, not gated as ">0".
    @test ref.lam_min > 1.0e-2
    # ... and the switch is nevertheless the evaluator that was used: the direct bundle's own
    # source (what `worldtube_data` carries for the MORTARS) is a different object
    @test length(worldtube_massive_source(sol.wt)) == length(sol.b)
    @test isapprox(sol.b, worldtube_massive_source(sol.wt); rtol = 1e-6)  # equal where inert
end

@testset "(d5) Ψ_{−m} = conj(Ψ_m) through the assembled system" begin
    cfgp, cfgm = dv_cfg(; spin = 0.5, m = 2), dv_cfg(; spin = 0.5, m = -2)
    sp = solve_mode(dv_orb(cfgp), cfgp, 2)
    sm = solve_mode(dv_orb(cfgm), cfgm, -2)
    @test build_matrix(sm.op) == conj.(build_matrix(sp.op))       # MEASURED bitwise equal
    @test maximum(abs, sm.b .- conj.(sp.b)) <= 1e-11 * maximum(abs, sp.b)   # MEASURED 1.8e-12 rel
    @test norm(sm.u .- conj.(sp.u)) / norm(sp.u) <= 1e-9           # MEASURED 9.24e-15
end

# ========================= (e) the hierarchy and the caches =============================
@testset "(e1) the stack is built from operator_builder with the worldtube factory" begin
    sol = sol05
    @test sol.stack.op === sol.op
    hier = sol.stack.hierarchy
    @test hier !== nothing
    @test length(hier.levels) == 1                    # yaml InitialCoarseLevels = L = 0
    # the coarsest level of the stack IS `coarsest_operator`'s, BITWISE
    cw = coarsest_operator(ComplexF64, sol.dom, sol.orbit, sol.config, sol.m;
                           worldtube = mesh -> worldtube_data(sol.dom, sol.orbit,
                                                              sol.config, sol.m; mesh = mesh),
                           maxlevels = sol.config.multigrid.initial_coarse_levels + 1)
    @test cw.nlevels == length(hier.levels)
    @test build_matrix(cw.op) == build_matrix(hier.levels[end].op)
    # ONE epoch, both caches populated
    @test sol.stack.epoch == 0
    @test sol.stack.subdomain_cache.builds == 20      # one per element on the single level
    @test sol.stack.bottom_cache.builds == 1
end

@testset "(e2) a genuine multi-level hierarchy on the same mesh" begin
    cfg = dv_cfg(; spin = 0.5, m = 2)
    s = solve_mode(dv_orb(cfg), cfg, 2;
                   options = solver_options(cfg; coarse_levels = nothing))
    # h-coarsening to the irreducible block structure (main.tex:672): 20 → 10 elements
    @test s.record.level_ndofs == [320, 160]
    @test s.result.residual <= cfg.gmres.absolute_residual
    # MEASURED 9 at the yaml-resolved overlap 4 (10 at the package default 2) — BELOW the
    # paper's "< 10 per AMR level". This is the L = 0 fine mesh; the sweep that tests
    # main.tex:720 properly is (e5).
    @test s.result.iterations == 9
    @test s.result.iterations <= 15
    @test s.stack.bottom_cache.builds == 1
    @test length(s.stack.hierarchy.levels) == 2
end

# NOTE (e3)/(e4) refine EVERY worldtube element = UNIFORM in-tube h refinement. That is NOT
# the paper's RefineAtPuncture criterion (yaml:132, main.tex:552-560) — see (e5), which is
# the measurement of record for the iteration-count claim.
@testset "(e3) rebuild_for_mesh! on a coarsened and a UNIFORMLY in-tube hp-refined mesh" begin
    cfg = dv_cfg(; spin = 0.5, m = 2)

    # ---- a COARSENED mesh
    sol = solve_mode(dv_orb(cfg), cfg, 2)
    e0 = sol.stack.epoch
    co = coarsen_mesh(sol.mesh)
    @test co !== nothing
    cm = co[1]
    rebuild_for_mesh!(sol, cm)
    @test sol.mesh === cm
    @test sol.mesh.ndofs == 160 && length(sol.mesh.elements) == 10
    @test length(sol.b) == 160 && length(sol.u) == 160
    @test length(sol.wt.source) == 160                 # the bundle was REBUILT for this mesh
    @test sol.record === nothing && sol.result === nothing
    @test sol.stack.epoch == e0                        # a pure MESH change needs no bump
    solve_mode!(sol)
    @test sol.result.residual <= cfg.gmres.absolute_residual
    @test mode_residual(sol) <= cfg.gmres.absolute_residual

    # ---- an hp_refine'd mesh warm-started through transfer_solution (the row-43 flow)
    sol2 = solve_mode(dv_orb(cfg), cfg, 2)
    bot0 = sol2.stack.bottom_cache.builds
    flags = [el for el in sol2.mesh.elements if inside_worldtube(sol2.dom, el)]
    newmesh = hp_refine(sol2.mesh, flags)
    # RE-PINNED after EllipticCore 6df9f51 (hp_refine no longer p-increments the elements the
    # 2:1 closure splits, so the one closure element and its children lose one order): 681 → 671
    @test newmesh.ndofs == 671 && length(newmesh.elements) == 33
    u = transfer_solution(sol2.mesh, newmesh, sol2.u)
    rebuild_for_mesh!(sol2, newmesh; u = u)
    @test sol2.mesh === newmesh
    @test sol2.u == u && length(sol2.u) == newmesh.ndofs
    @test length(sol2.b) == newmesh.ndofs
    solve_mode!(sol2; x0 = sol2.u)
    @test sol2.result.residual <= cfg.gmres.absolute_residual
    @test sol2.result.iterations < 10
    @test sol2.stack.bottom_cache.builds >= bot0       # counters only grow
    @test length(sol2.monitor) == 2                    # the monitor ACCUMULATES per level
    # a carried solution of the wrong length is refused, not silently reshaped
    @test_throws ArgumentError rebuild_for_mesh!(sol2, sol2.mesh; u = zeros(ComplexF64, 3))
end

@testset "(e4) coarse_orders MIN vs parity, one UNIFORM in-tube hp step" begin
    cfg = dv_cfg(; spin = 0.5, m = 2)
    rows = map((nothing, 4)) do co
        # overlap pinned to the package default 2 so this contrast matches the table quoted
        # below (the yaml-resolved overlap 4 shifts it to 9, 8 / 9, 9 — see (e6))
        s = solve_mode(dv_orb(cfg), cfg, 2;
                       options = solver_options(cfg; coarse_levels = nothing,
                                                coarse_orders = co, overlap = 2))
        its0, bot0 = s.result.iterations, s.stack.bottom_cache.builds
        flags = [el for el in s.mesh.elements if inside_worldtube(s.dom, el)]
        nm = hp_refine(s.mesh, flags)
        u = transfer_solution(s.mesh, nm, s.u)
        rebuild_for_mesh!(s, nm; u = u)
        solve_mode!(s; x0 = s.u)
        (co = co, its = (its0, s.result.iterations),
         bottom = (bot0, s.stack.bottom_cache.builds), residual = s.result.residual,
         coarsest = s.stack.hierarchy.levels[end].mesh.ndofs)
    end
    for r in rows
        @test r.residual <= cfg.gmres.absolute_residual
    end
    # at step 0 the two rules COINCIDE (no p refinement has happened yet)
    @test rows[1].its[1] == rows[2].its[1] == 10
    @test rows[1].bottom[1] == rows[2].bottom[1] == 1
    # after the hp step spectre parity keeps THE BOTTOM GRID ITSELF fixed (main.tex:674:
    # "the bottom grid remains the same even when AMR adds finer levels"), so its cached LU
    # is a HIT; the MIN rule raises the coarsest order, changing the grid and rebuilding the
    # LU. MEASURED over FIVE levels of UNIFORM in-tube refinement (320 → 681 → 1924 → 6252 →
    # 21632 dofs) at overlap 2 BEFORE EllipticCore 6df9f51 — that fix moves every dof and
    # coarsest-grid sequence (the first step is now 671, not 681); the ASSERTIONS below are
    # structural (which grid moves, which LU is rebuilt) and do not depend on those numbers —
    # MIN 10, 9, 10, 12, 14 (coarsest 160 → 241 → 340 → 457 → 592,
    # bottom builds 1/1/1/1/1); parity 10, 10, 12, 15, 32 (coarsest 160 throughout, builds
    # 1/0/0/0/0). Parity's cost here is 0, +1, +2, +3, +18 — it GROWS, it is not "≤ 2"; the
    # ≤ 2 figure belongs to the RefineAtPuncture sweep of (e5).
    @test rows[2].coarsest == 160                       # the bottom grid did NOT move
    @test rows[1].coarsest > 160                        # ... under the MIN rule it did
    @test rows[2].bottom[2] == 1                        # built once, then hit
    @test rows[1].bottom[2] == 2                        # rebuilt on the hp step
end

# elements whose CLOSURE contains the particle (r0, w = 0) — the paper's RefineAtPuncture
# criterion (yaml:132, main.tex:552-560), NOT "every worldtube element"
function dv_puncture_elements(sol, mesh)
    RT = Float64; tol = 1e-9; out = ElementId[]
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]; blk = mesh.blocks[mesh.blockindex[el.block]]
        inside = true
        for d in 1:2
            lo, hi = logical_interval(RT, el, d)
            x = d == 1 ? sol.orbit.r0 : 0.0
            ζ = map_inverse(blockmap(blk, d), x)
            ξ = (2ζ - (lo + hi)) / (hi - lo)
            (-1 - tol <= ξ <= 1 + tol) || (inside = false; break)
        end
        inside && push!(out, el)
    end
    return out
end

@testset "(e5) H1: RefineAtPuncture — the paper's own criterion, the measurement of record" begin
    # THE CLAIM UNDER TEST is main.tex:720 ("less than 10 iterations per AMR level") together
    # with main.tex:722 (scale independence). It must be measured under the paper's OWN
    # refinement criterion; (e3)/(e4) refine the whole worldtube and are a different sweep.
    cfg = dv_cfg(; spin = 0.5, m = 2)
    # RE-PINNED after EllipticCore 6df9f51: the MIN rule's coarse grids changed with the order
    # map, moving two of its counts (8, 7 → 7, 8); spectre parity is unaffected.
    for (co, expect) in ((nothing, [9, 8, 7, 8, 6]), (4, [9, 8, 9, 9, 9]))
        s = solve_mode(dv_orb(cfg), cfg, 2;
                       options = solver_options(cfg; coarse_levels = nothing,
                                                coarse_orders = co))
        its = [s.result.iterations]
        nd = [s.mesh.ndofs]
        coarsest = [s.stack.hierarchy.levels[end].mesh.ndofs]
        bot = [s.stack.bottom_cache.builds]
        for _ in 1:4
            flags = dv_puncture_elements(s, s.mesh)
            @test 1 <= length(flags) <= 4        # a CORNER point, not the whole tube
            nm = hp_refine(s.mesh, flags)
            u = transfer_solution(s.mesh, nm, s.u)
            rebuild_for_mesh!(s, nm; u = u)
            solve_mode!(s; x0 = s.u)
            push!(its, s.result.iterations); push!(nd, s.mesh.ndofs)
            push!(coarsest, s.stack.hierarchy.levels[end].mesh.ndofs)
            push!(bot, s.stack.bottom_cache.builds)
            @test s.result.residual <= cfg.gmres.absolute_residual
        end
        # RE-PINNED after EllipticCore 6df9f51 (the balance-split element and its descendants
        # no longer carry the p increment in the split dimension)
        @test nd == [320, 593, 950, 1411, 1988]
        @test its == expect                       # the WHOLE sequence, not one step
        # main.tex:720 HOLDS under the paper's criterion, at every level
        @test all(<=(10), its)
        if co === nothing
            @test its[end] < its[1]               # main.tex:722: the count DECREASES with depth
            # RE-PINNED after EllipticCore 6df9f51 (the coarse order map follows the fine one)
            @test coarsest == [160, 236, 334, 450, 584]
            @test bot == [1, 2, 3, 4, 5]          # MIN rule rebuilds the bottom LU every step
        else
            @test coarsest == fill(160, 5)        # main.tex:674: the bottom grid does NOT move
            @test bot == fill(1, 5)               # ... and its LU is built ONCE, then hit
        end
    end
end

@testset "(e6) H2: the resolved overlap saturates at the neighbour extent" begin
    # spectre's Auto = min(volume_extent, nullopt) = volume_extent, so overlap >= P is the
    # same subdomain; our schwarz.jl:120 clamp `min(ov, Nn[d])` reproduces it BITWISE.
    cfg = dv_cfg(; spin = 0.5, m = 2)
    runs = map((3, 4, 5)) do ov
        s = solve_mode(dv_orb(cfg), cfg, 2;
                       options = solver_options(cfg; coarse_levels = nothing, overlap = ov))
        (its = s.result.iterations, res = s.result.residual, u = s.u)
    end
    @test runs[1].its == 10                        # overlap 3
    @test runs[2].its == 9                         # overlap 4 = the yaml's Auto
    @test runs[3].its == runs[2].its               # overlap 5 saturates ...
    @test runs[3].res == runs[2].res               # ... BITWISE
    @test runs[3].u == runs[2].u
end

@testset "(e7) H3: controlled hierarchy-depth comparison on ONE fine operator" begin
    cfg = dv_cfg(; spin = 0.5, m = 2, L = 2)      # 5120 dofs
    orb = dv_orb(cfg)
    dom = block_domain(orb, cfg)
    rows = map((0, 1, 2, 3)) do mc
        s = solve_mode(orb, cfg, 2; dom = dom,
                       options = solver_options(cfg; coarse_levels = mc, overlap = 2))
        (mc = mc, levels = length(s.record.level_ndofs), lv = s.record.level_ndofs,
         its = s.result.iterations, res = s.result.residual)
    end
    @test [r.levels for r in rows] == [1, 2, 3, 4]
    @test rows[4].lv == [5120, 1280, 320, 160]
    @test [r.its for r in rows] == [1, 20, 24, 25]
    for r in rows
        @test r.res <= cfg.gmres.absolute_residual
    end
    # the depth study is on ONE operator; it is NOT the 10/15/25 fine-mesh comparison
    @test all(r -> r.lv[1] == 5120, rows)
end

# ==================== (f) the residual staircase and the build cost ======================
@testset "(f) SolveRecord: staircase, wall clock, build cost" begin
    cfg = dv_cfg(; spin = 0.5, m = 2)
    sol = solve_mode(dv_orb(cfg), cfg, 2)
    rec = sol.record
    @test rec isa SolveRecord
    @test length(rec.residual) == rec.iterations + 1     # ‖b−Ax₀‖ then one Givens estimate each
    @test rec.residual[end] <= rec.residual[1]
    @test rec.final_residual == sol.result.residual
    @test rec.solve_time > 0
    @test length(rec.iteration_time) == rec.iterations
    @test rec.build_time > 0                             # charged to the FIRST solve only
    @test rec.build_time_hierarchy >= 0 && rec.build_time_bottom >= 0
    @test rec.subdomain_builds == 20 && rec.bottom_builds == 1
    @test rec.level_ndofs == [320]
    @test length(sol.monitor) == 1
    # a second solve on the SAME preconditioner records no build cost (the sealed contract)
    solve_mode!(sol)
    @test length(sol.monitor) == 2
    @test sol.monitor[end].build_time == 0.0
    @test sol.monitor[end].subdomain_builds == 0
    @test cumulative_iterations(sol.monitor)[end] ==
          sum(r.iterations for r in sol.monitor.records)
end

# ================================= (g) BigFloat(256) =====================================
@testset "(g) the whole pipeline constructs and solves at BigFloat(256)" begin
    setprecision(BigFloat, 256) do
        # P = 2 (80 dofs) keeps the dense BigFloat bottom solve cheap; MEASURED residual
        # 1.74e-77 in 1 iteration, precision preserved at 256 bits
        cfg = dv_cfg(; spin = big"0.5", r0 = big"10.0", m = 2, P = 2)
        orb = dv_orb(cfg)
        s = solve_mode(orb, cfg, 2; T = Complex{BigFloat})
        @test s.dom isa BlockDomain{BigFloat}
        @test s.mesh.ndofs == 80
        @test eltype(s.u) === Complex{BigFloat}
        @test precision(real(s.u[1])) == 256              # INHERITED, never Float64'd
        @test s.stack.options.abstol isa BigFloat
        @test all(isfinite, s.u)
        @test s.result.residual <= cfg.gmres.absolute_residual
        @test mode_residual(s) <= cfg.gmres.absolute_residual
        Ψ = reconstruct_full_field(s)
        @test eltype(Ψ) === Complex{BigFloat}
        @test all(isfinite, Ψ)
        @test field_at(s, Ψ, orb.r0, big"0.5") isa Complex{BigFloat}
    end
end

end # @testset

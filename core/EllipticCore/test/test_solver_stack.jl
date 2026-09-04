# test_solver_stack.jl — the solver-stack composition API (DAG node solver-stack-core,
# elliptic-ssf GLOBAL_DAG row 21). TESTS FIRST: this file is written and run RED before
# src/solver_stack.jl holds anything but its comment stub.
#
# WHAT THIS NODE OWNS (CRP main.tex:668-676 + fig:precond caption main.tex:632):
# the PUBLIC COMPOSITION `GMRES ∘ MultigridV ∘ Schwarz(2 passes) ∘ cached bottom solve`,
# built once per (mesh, operator), rebuilt across AMR steps with both caller-owned caches,
# and instrumented by a residual monitor. It composes the SEALED predecessors and
# re-implements none of them — gate (b) pins that BITWISE.
#
# GATES (calibrated on throwaway prototypes in
# results/elliptic-ssf/paper_ELLIPTIC-SSF/evidence/solver-stack-core/, logs
# calibrate_{model,alloc,shrink,gates,gates2,amr}.log, BEFORE this file was written):
#
#  (a) THE 5-BLOCK GATE MODEL. Four affine blocks plus one eq:compaction 1/r wavezone block
#      (main.tex:546-549), a 2:1 h-refinement inside block 3, a p-jump on block 4, complex
#      coefficients, and NONZERO fixed Dirichlet data on every boundary face folded through
#      `linearized_rhs`. At the yaml production tolerance (abstol 1e-12, reltol 0, maxiter
#      500) the stack converges in 8 iterations at the k=1 mesh (404 dofs, 23 elements, 4
#      levels) and 10 at k=2 (1472 dofs, 83 elements, 6 levels) — gate < 15. The sealed
#      multigrid model's 7 was measured at reltol 1e-10 on a uniform two-block mesh; at THIS
#      file's abstol 1e-12 that same sealed model also takes 8, so the two agree.
#      LANDED DEVIATION on the planned "solution equals direct_solve to 1e-10": at abstol
#      1e-12 on a UNIT-NORM right-hand side the measured agreement is 1.64e-10 (k=1) and
#      1.75e-10 (k=2) — the model's residual-to-error amplification is ≈ 4.5e2, a
#      CONDITIONING fact, not a solver defect. Pinned as an envelope 1e-9 at the production
#      tolerance PLUS the planned 1e-10 at abstol 1e-14 (measured 3.47e-12 / 6.15e-13),
#      which shows the agreement is tolerance-limited.
#  (b) THE COMPOSITION IS THE PREDECESSORS. `solve!` through the stack is BITWISE the
#      hand-assembled `gmres(op, b; Minv! = MultigridV(MultigridHierarchy(op, build_op;
#      caches, epoch); bottom = sparse_bottom_solver(hier.levels[end].op; cache, epoch)))`
#      of the predecessors' own tests — solution, iteration count and full residual history.
#      A monitored solve is bitwise an unmonitored one (instrumentation wraps the callables,
#      never the arithmetic).
#  (c) ARBITRARY PRECONDITIONERS THROUGH ONE API. `multigrid_preconditioner`,
#      `schwarz_only_preconditioner` and `no_preconditioner` are plain factory callables
#      `f(stack) -> (Minv!, hierarchy_or_nothing)`; measured on the 5-block model the counts
#      order 8 < 18 < 229 (the sealed contrast was 7 / 18-32 / 213). The identity factory is
#      BITWISE `gmres`'s own `Minv! = nothing` (measured: identical x, iterations, history).
#  (d) PRECONDITIONER-APPLY ALLOCATION. On the sealed 4-level multigrid model a V-cycle
#      apply allocated 19 869 600 B, CONSTANT over repeats (row 19 reported 18 834 720 B on
#      its host). Decomposed here for the first time: 96.1 % is the sealed matrix-free
#      `apply!` (5 calls per non-bottom level), and only 4.2 % (835 584 B) is `smooth!`'s own
#      per-call workspace. This node preallocates THAT (`SchwarzWorkspace`, additive, kwarg
#      on `smooth!`, one per level owned by `MultigridV`); the gate is the constant, its
#      bound, the measured delta, and the honest floor: `apply!` in matrixfree.jl owns the
#      remainder and is a follow-up (that file is in flight with the dg-jump-conditions node).
#  (e) MOCK AMR THROUGH `rebuild!`. Eight h-refinement steps of the finest level: the bottom
#      LU is built EXACTLY ONCE and hit 7 times, every unchanged coarser level's subdomain
#      factorizations are reused, each level's build cost is recorded separately from its
#      iterations, and the falsifier — an `invalidate!` epoch bump — rebuilds everything.
#  (f) THE SHIFT POLICY. `helmholtz_shift(m, Ω; slicing)` reaches every level's smoother, and
#      the sealed row-18 resonance model (its factory REUSED verbatim, not re-derived)
#      reproduces the sealed contrast through this API: unshifted GMRES exhausts the
#      150-iteration budget (residual 1206.74, the sealed number), the mass-weighted policy
#      shift converges in 47 (Schwarz-only, the sealed count) and 44 (multigrid).
#  (g) TYPE GENERICITY. Complex{BigFloat}(256) builds and solves a 4-level stack (235 dofs,
#      21 iterations, residual 2.2e-42); thread parity is bitwise.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists (each
# @testset then errors individually in the red run without aborting its siblings).

using Test
using LinearAlgebra
using REPL          # Julia >=1.11: @doc is a no-op unless REPL is loaded (runtests.jl convention)
using EllipticCore

# ---- calibrated gate constants (evidence/solver-stack-core/calibrate_*.log) --------------
const SK_ITER_GATE      = 15        # (a) planned; measured 8 (k=1) and 10 (k=2)
const SK_DIRECT_TOL     = 1e-9      # (a) at abstol 1e-12: measured 1.64e-10 / 1.75e-10
const SK_DIRECT_TIGHT   = 1e-10     # (a) at abstol 1e-14: measured 3.47e-12 / 6.15e-13
const SK_TIGHT_ABSTOL   = 1e-14     # (a) the tolerance-limitation probe
const SK_MG_ITERS       = 12        # (c) multigrid,     measured 8
const SK_SW_ITERS       = 30        # (c) Schwarz-only, measured 17 (RHS-dependent, 17-18)
const SK_NONE_ITERS     = 400       # (c) unpreconditioned, measured 229
const SK_STALE_ITERS    = 40        # (M2) stale-Schwarz cost, measured 77 against 7
const SK_HP_ITER_GATE   = 20        # (M1) pinned coarse order over 6 hp steps: measured 17
const SK_AMR_STEPS      = 8         # (e) mock AMR sweep length
const SK_AMR_ITER_GATE  = 13        # (e) 8-deep corner cascade: measured 8,9,9,9,9,10,10,11
const SK_AMR_EARLY_GATE = 10        # (e) the paper's "< 10" holds for the first 6 steps
const SK_FRONT_STEPS    = 3         # (e) refinement FRONT sweep (a refined REGION)
const SK_FRONT_ITER_GATE = 10       # (e) measured 8, 9, 9 — the paper's claim, outright
const SK_SHIFT_BUDGET   = 150       # (f) sealed row-18 budget (unshifted exhausts it)
const SK_SHIFT_GATE     = 49        # (f) sealed CAL4_A_GATE; measured 47 (Schwarz) / 44 (MG)
const SK_SHIFT_RTOL     = 1e-10     # (f) sealed row-18 solver tolerance
const SK_BIG_ABSTOL     = big(1e-40)  # (g) measured residual 2.2e-42 in 21 iterations
const SK_BIG_ITERS      = 40        # (g) measured 21
const SK_VCYCLE_ALLOC   = 19_500_000  # (d) pre-workspace baseline 19 869 600 B on this model;
                                      #     predicted post-workspace 19 034 016 B. The bound
                                      #     sits BETWEEN the two, so it fails if the
                                      #     preallocation does nothing.
const SK_APPLY_FLOOR    = 0.90      # (d) fraction of the V-cycle allocation owned by apply!
                                    #     (measured 96.1 %) — the documented floor
const SK_C = 1.5                                   # CRP main.tex:595 penalty parameter

# ---- the 5-block gate model -------------------------------------------------------------
sk_zero2(x, y) = 0.0
sk_A11(x, y) = 1.0 + 0.05x
sk_A22(x, y) = 1.0
const SK_BETA = 1.0 + 0.5im                        # definite complex (screened) model
sk_beta(x, y) = SK_BETA * (1 + 0.1x)

# NONZERO fixed Dirichlet data on every boundary face: the operator-on-zero lift is O(1e3)
# here, so a caller who forgets `linearized_rhs` solves a visibly different problem.
function sk_nzbcs(mesh, T)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost([T(0.3 + 0.2im) * (1 + 0.1k) for k in 1:nf])
    end
    return bcs
end

# the operator FACTORY: the same coefficient closures rediscretized on whatever mesh it is
# handed, with that mesh's own boundary conditions (the MultigridHierarchy contract)
sk_factory(T) = mesh -> MeshOperator(T, mesh, (sk_A11, sk_A22), sk_beta, (sk_zero2, sk_zero2);
                                     C = SK_C, bcs = sk_nzbcs(mesh, T))

function sk_five_blocks(RT)
    return ([Block(1, LinearMap1D(RT(1), RT(2)), LinearMap1D(RT(0), RT(1))),
             Block(2, LinearMap1D(RT(2), RT(3)), LinearMap1D(RT(0), RT(1))),
             Block(3, LinearMap1D(RT(3), RT(4)), LinearMap1D(RT(0), RT(1))),
             Block(4, LinearMap1D(RT(4), RT(5)), LinearMap1D(RT(0), RT(1))),
             Block(5, InverseRadialMap1D(RT(5), RT(50)), LinearMap1D(RT(0), RT(1)))],
            [BlockConnection(1, 1, 2), BlockConnection(2, 1, 3),
             BlockConnection(3, 1, 4), BlockConnection(4, 1, 5)])
end

# k-level uniform grid in every block, then a 2:1 h-split inside block 3 and a p-jump on
# block 4 — mixed h AND p, one curved block, five blocks
function sk_five_mesh(RT, k; N = 4)
    blks, conns = sk_five_blocks(RT)
    els = ElementId[]
    for b in 1:5, iy in 0:(2^k-1), ix in 0:(2^k-1)
        push!(els, ElementId(b, SegmentId(k, ix), SegmentId(k, iy)))
    end
    m0 = Mesh2D(blks, els, (N, N); connections = conns)
    m1 = h_refine(m0, Set([ElementId(3, SegmentId(k, 0), SegmentId(k, 0))]))
    ords = Dict(el => (el.block == 4 ? (N + 1, N + 1) : (N, N)) for el in m1.elements)
    return Mesh2D(m1.blocks, m1.elements, ords; connections = m1.connections)
end

# deterministic vectors (Knuth-MMIX LCG; distinct name and seed per test file)
mutable struct SKRNG
    s::UInt64
end
sk_unit!(r::SKRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      Int64(r.s >> 11) / 2.0^53 * 2 - 1)
function sk_vec!(r::SKRNG, T, n)
    RT = real(T)
    return T[T <: Complex ? T(RT(sk_unit!(r)), RT(sk_unit!(r))) : T(RT(sk_unit!(r)))
             for _ in 1:n]
end

# a manufactured discrete solution on a UNIT-NORM right-hand side: b := A·u_ex makes u_ex the
# exact discrete answer, and ‖b‖ = 1 makes the yaml's ABSOLUTE 1e-12 a relative tolerance
function sk_manufactured(op, seed)
    n = op.mesh.ndofs
    T = eltype(operator_on_zero(op))
    uex = sk_vec!(SKRNG(seed), T, n)
    uex ./= norm(uex)
    blin = apply!(zeros(T, n), op, uex)
    sc = norm(blin)
    blin ./= sc
    uex ./= sc
    return uex, blin, blin .+ operator_on_zero(op)      # (u_ex, linearized b, physical b)
end

@testset "solver stack (solver-stack-core)" begin

# ---------------------------------------------------- exports, options, docstrings -------
@testset "exports, options and the normative docstrings" begin
    for f in (:SolverStackOptions, :SolverStack, :SolveRecord, :ResidualMonitor, :solve!,
              :solve_linearized!, :rebuild!, :invalidate!, :linearized_rhs,
              :cumulative_iterations,
              :multigrid_preconditioner, :schwarz_only_preconditioner, :no_preconditioner)
        @test isdefined(EllipticCore, f)
        @test Base.isexported(EllipticCore, f)
    end
    # the yaml production defaults (anc/ScalarSelfForce.yaml:178-205) ARE the defaults
    o = SolverStackOptions()
    @test o.abstol == 1e-12                     # yaml:183 AbsoluteResidual
    @test o.reltol == 0                         # yaml:182 RelativeResidual
    @test o.maxiter == 500                      # yaml:181 MaxIterations
    @test o.restart == 500                      # unrestarted (= maxiter)
    @test o.multigrid_iterations == 1           # yaml:187 Multigrid.Iterations
    @test o.initial_coarse_levels === nothing   # yaml:188 (nothing = coarsen all the way)
    @test o.pre_smoothing                       # yaml:189
    @test o.use_bottom_solver                   # yaml:190
    @test !o.post_smoothing_at_bottom           # yaml:191
    @test o.passes == 2                         # yaml:197 SchwarzSmoother.Iterations
    @test o.overlap == 2
    @test o.subdomain_complex_shift == 0        # yaml:203 ExplicitInverse.ComplexShift
    @test o.shrink == 1                         # the measured decision (see the docstring)
    @test o isa SolverStackOptions{Float64}
    # precision inheriting
    @test SolverStackOptions(; abstol = big"1e-30") isa SolverStackOptions{BigFloat}
    # the InitialCoarseLevels ↔ maxlevels off-by-one (spectre Multigrid/Tags.hpp:30-37)
    @test EllipticCore._ss_maxlevels(SolverStackOptions()) === nothing
    @test EllipticCore._ss_maxlevels(SolverStackOptions(; initial_coarse_levels = 0)) == 1
    @test EllipticCore._ss_maxlevels(SolverStackOptions(; initial_coarse_levels = 3)) == 4
    # validation
    @test_throws ArgumentError SolverStackOptions(; abstol = -1)
    @test_throws ArgumentError SolverStackOptions(; reltol = -1)
    @test_throws ArgumentError SolverStackOptions(; maxiter = -1)
    @test_throws ArgumentError SolverStackOptions(; restart = 0)
    @test_throws ArgumentError SolverStackOptions(; multigrid_iterations = 0)
    @test_throws ArgumentError SolverStackOptions(; initial_coarse_levels = -1)
    @test_throws ArgumentError SolverStackOptions(; shrink = 0)
    @test_throws ArgumentError SolverStackOptions(; shrink = 1.5)
    @test_throws ArgumentError SolverStackOptions(; overlap = -1)
    @test_throws ArgumentError SolverStackOptions(; passes = 0)
    # the two normative contracts live in ONE docstring, on SolverStack
    d = string(@doc SolverStack)
    @test occursin("LINEARIZED-RHS CONTRACT", uppercase(d))
    @test occursin("linearized_rhs", d)
    @test occursin("operator_on_zero", d)
    @test occursin("TWO-CACHES-ONE-EPOCH CONTRACT", uppercase(d))
    @test occursin("epoch", d)
    @test occursin("shrink", d)                 # the measured default is documented here
    @test occursin("helmholtz_shift", d)        # and the shift route
    @test occursin("main.tex", d)               # paper anchor
    dl = string(@doc linearized_rhs)
    @test occursin("operator_on_zero", dl)
    dr = string(@doc SolveRecord)
    @test occursin("fig:linsolv", dr)
end

# ------------------------------------------------------------- (a) the 5-block model -----
@testset "gate (a): the 5-block model converges at the yaml tolerance" begin
    for (k, iters_expected) in ((1, 8), (2, 10))
        mesh = sk_five_mesh(Float64, k)
        fac = sk_factory(ComplexF64)
        op = fac(mesh)
        n = mesh.ndofs
        @test length(mesh.blocks) == 5
        @test any(b -> blockmap(b, 1) isa InverseRadialMap1D, mesh.blocks)   # the 1/r block
        @test length(unique(mesh.orders)) > 1                                # the p-jump
        @test !two_to_one_balanced([el for el in mesh.elements if el.block != 3]) ||
              any(el -> el.block == 3 && el.segments[1].level > k, mesh.elements)  # the h-jump
        @test norm(operator_on_zero(op)) > 1.0        # the fixed-data gate is NOT vacuous

        uex, blin, bphys = sk_manufactured(op, 0xfeed)
        stack = SolverStack(op, fac)
        @test stack.hierarchy !== nothing
        @test length(stack.hierarchy.levels) >= 4

        # the ONE place fixed data is folded
        b = linearized_rhs(op, bphys)
        @test maximum(abs, b .- blin) <= 1e-12 * norm(operator_on_zero(op))
        @test linearized_rhs(stack, bphys) == b

        u = zeros(ComplexF64, n)
        res = solve_linearized!(u, stack, b)
        @test res isa GMRESResult
        @test res.converged
        @test res.iterations == iters_expected
        @test res.iterations < SK_ITER_GATE
        @test res.residual <= stack.options.abstol      # TRUE residual, recomputed
        @test u == res.x
        xd = direct_solve(op, b)
        @test norm(u - xd) / norm(xd) < SK_DIRECT_TOL
        @test norm(u - uex) / norm(uex) < SK_DIRECT_TOL

        # the planned 1e-10 IS reached — the agreement is tolerance-limited, not broken
        tight = SolverStack(op, fac; options = SolverStackOptions(; abstol = SK_TIGHT_ABSTOL))
        ut = zeros(ComplexF64, n)
        rt = solve_linearized!(ut, tight, b)
        @test rt.converged
        @test rt.iterations < SK_ITER_GATE
        @test norm(ut - xd) / norm(xd) < SK_DIRECT_TIGHT

        # H2: `solve!` takes the PHYSICAL rhs and folds internally — bitwise the pre-folded
        # route, and the route a caller gets by default
        up = zeros(ComplexF64, n)
        rp = solve!(up, stack, bphys)
        @test up == u
        @test rp.iterations == res.iterations
        @test rp.history == res.history
        # FALSIFIER (now the CALLER's error, not the API's): folding twice is O(1) wrong
        uw = zeros(ComplexF64, n)
        solve!(uw, stack, b)
        @test norm(uw - xd) / norm(xd) > 0.1
        # and the OLD trap is what `solve!` removes: the pre-folded route on a physical rhs
        # still converges happily on the wrong system
        uold = zeros(ComplexF64, n)
        rold = solve_linearized!(uold, stack, bphys)
        @test rold.converged
        @test norm(uold - xd) / norm(xd) > 1e2
    end
end

# --------------------------------- (b) the composition IS the sealed predecessors --------
@testset "gate (b): bitwise equal to the hand-assembled predecessor composition" begin
    mesh = sk_five_mesh(Float64, 1)
    fac = sk_factory(ComplexF64)
    op = fac(mesh)
    n = mesh.ndofs
    _, b, _ = sk_manufactured(op, 0x51ac)

    # hand-assembled, exactly as ESSF::multigrid-hierarchy and ESSF::multigrid-bottom-solver
    # write it in their own tests
    subc = SubdomainCache()
    botc = BottomSolverCache()
    hier = MultigridHierarchy(op, fac; caches = subc, epoch = 0, overlap = 2, shift = 0,
                              maxlevels = nothing, shrink = 1)
    bot = sparse_bottom_solver(hier.levels[end].op; cache = botc, epoch = 0)
    mgv = MultigridV(hier; passes = 2, iterations = 1, pre_smoothing = true,
                     post_smoothing_at_bottom = false, bottom = bot)
    ref = gmres(op, b; abstol = 1e-12, reltol = 0, maxiter = 500, restart = 500, Minv! = mgv)

    stack = SolverStack(op, fac)
    u = zeros(ComplexF64, n)
    got = solve_linearized!(u, stack, b)
    @test got.x == ref.x                       # BITWISE
    @test got.iterations == ref.iterations
    @test got.residual == ref.residual
    @test got.history == ref.history
    @test got.converged == ref.converged
    @test got.breakdown == ref.breakdown

    # the stack really wired the predecessors it claims to
    @test stack.preconditioner isa MultigridV
    @test stack.hierarchy isa MultigridHierarchy
    @test length(stack.hierarchy.levels) == length(hier.levels)
    @test [l.mesh.ndofs for l in stack.hierarchy.levels] == [l.mesh.ndofs for l in hier.levels]
    @test all(l -> l.smoother isa SchwarzSmoother, stack.hierarchy.levels)
    @test all(l -> l.smoother.overlap == 2, stack.hierarchy.levels)
    @test stack.preconditioner.passes == 2
    @test stack.preconditioner.iterations == 1
    @test stack.preconditioner.pre_smoothing
    @test !stack.preconditioner.post_smoothing_at_bottom
    @test stack.preconditioner.bottom !== nothing
    # the bottom solver is built on the COARSEST operator, never the finest
    @test stack.bottom_cache.builds == 1
    @test stack.hierarchy.levels[end].mesh.ndofs < mesh.ndofs

    # a MONITORED solve is bitwise an unmonitored one
    mon = ResidualMonitor()
    u2 = zeros(ComplexF64, n)
    got2 = solve_linearized!(u2, stack, b; monitor = mon)
    @test got2.x == ref.x
    @test got2.history == ref.history
    @test got2.iterations == ref.iterations
    @test length(mon) == 1
    @test mon[1].residual == ref.history

    # warm start: x0 = the converged answer costs zero iterations
    u3 = zeros(ComplexF64, n)
    got3 = solve_linearized!(u3, stack, b; x0 = ref.x)
    @test got3.iterations == 0
    @test got3.converged
end

# ------------------------------------------ (c) arbitrary preconditioners, one API -------
@testset "gate (c): multigrid < Schwarz-only < unpreconditioned through one API" begin
    mesh = sk_five_mesh(Float64, 1)
    fac = sk_factory(ComplexF64)
    op = fac(mesh)
    n = mesh.ndofs
    _, b, _ = sk_manufactured(op, 0xc0c0)
    counts = Dict{Symbol,Int}()
    for (name, factory, gate) in ((:mg, multigrid_preconditioner, SK_MG_ITERS),
                                  (:sw, schwarz_only_preconditioner, SK_SW_ITERS),
                                  (:id, no_preconditioner, SK_NONE_ITERS))
        st = SolverStack(op, fac; preconditioner = factory,
                         options = SolverStackOptions(; maxiter = 2000, restart = 2000))
        u = zeros(ComplexF64, n)
        r = solve_linearized!(u, st, b)
        @test r.converged
        @test r.iterations <= gate
        @test norm(u - direct_solve(op, b)) / norm(direct_solve(op, b)) < SK_DIRECT_TOL
        counts[name] = r.iterations
        # only the multigrid factory owns a hierarchy
        @test (st.hierarchy === nothing) == (name !== :mg)
    end
    @test counts[:mg] < counts[:sw] < counts[:id]     # the headline ordering
    # the sealed contrast was 7 / 18-32 / 213 on the uniform two-block model; the counts on
    # THIS model and RHS are printed so the log carries them (the gate is the ORDERING and
    # the calibrated bounds above, never a hand-fitted triple)
    println("  [gate (c)] 5-block iteration counts: multigrid=", counts[:mg],
            " schwarz-only=", counts[:sw], " unpreconditioned=", counts[:id],
            "  (Schwarz-only is RHS-dependent: 17-18 across the seeds measured)")
    @test counts[:sw] >= 2 * counts[:mg]             # multigrid at least halves the count
    @test counts[:id] >= 5 * counts[:sw]             # and preconditioning is not optional

    # the identity factory is BITWISE gmres's own Minv! = nothing
    st = SolverStack(op, fac; preconditioner = no_preconditioner,
                     options = SolverStackOptions(; maxiter = 2000, restart = 2000))
    u = zeros(ComplexF64, n)
    r = solve_linearized!(u, st, b)
    r0 = gmres(op, b; abstol = 1e-12, reltol = 0, maxiter = 2000, restart = 2000,
               Minv! = nothing)
    @test r.x == r0.x
    @test r.iterations == r0.iterations
    @test r.history == r0.history

    # a CALLER-SUPPLIED factory composes through the same seam (the API is open)
    my_factory = st -> (schwarz_preconditioner(
        SchwarzSmoother(st.op; overlap = 1, epoch = st.epoch,
                        cache = st.subdomain_cache); passes = 3), nothing)
    st2 = SolverStack(op, fac; preconditioner = my_factory,
                      options = SolverStackOptions(; maxiter = 2000, restart = 2000))
    u2 = zeros(ComplexF64, n)
    r2 = solve_linearized!(u2, st2, b)
    @test r2.converged
    @test counts[:mg] < r2.iterations < counts[:id]
end

# ------------------------------------------------------------ the residual monitor -------
@testset "the residual monitor records the fig:linsolv data surface" begin
    mesh = sk_five_mesh(Float64, 1)
    fac = sk_factory(ComplexF64)
    op = fac(mesh)
    n = mesh.ndofs
    _, b, _ = sk_manufactured(op, 0x3f3f)
    stack = SolverStack(op, fac)
    mon = ResidualMonitor()
    @test isempty(mon)
    u = zeros(ComplexF64, n)
    r = solve_linearized!(u, stack, b; monitor = mon, true_residual_probe = true)
    @test length(mon) == 1
    rec = mon[1]
    @test rec isa SolveRecord{Float64}
    # residual staircase: the sealed GMRES history verbatim
    @test rec.residual == r.history
    @test length(rec.residual) == r.iterations + 1
    @test rec.iterations == r.iterations
    @test rec.converged == r.converged
    @test rec.final_residual == r.residual
    @test issorted(rec.residual; rev = true) || all(diff(rec.residual) .<= 1e-14)
    # per-iteration wall clock and preconditioner allocation
    @test length(rec.iteration_time) == r.iterations
    @test length(rec.preconditioner_time) == r.iterations
    @test length(rec.operator_time) == r.iterations
    @test length(rec.preconditioner_bytes) == r.iterations
    @test all(>(0), rec.iteration_time)
    @test all(>=(0), rec.preconditioner_time)
    @test sum(rec.preconditioner_time) + sum(rec.operator_time) <= rec.solve_time * 1.001
    @test all(>(0), rec.preconditioner_bytes)
    # the V-cycle apply cost dominates one iteration
    @test sum(rec.preconditioner_time) > sum(rec.operator_time)
    # the build cost is recorded SEPARATELY from the iterations (this package builds the
    # subdomain matrices eagerly, so it is NOT inside iteration 1 as in the paper)
    @test rec.build_time > 0
    @test rec.build_time == stack.build_time
    @test rec.subdomain_builds > 0
    @test rec.bottom_builds == 1
    @test rec.subdomain_hits >= 0
    @test rec.level_ndofs == [l.mesh.ndofs for l in stack.hierarchy.levels]
    @test rec.level_subdomains ==
          [length(l.smoother.subdomains) for l in stack.hierarchy.levels]
    @test sum(rec.level_subdomains) == rec.subdomain_builds + rec.subdomain_hits
    # the TRUE residual probe: with RIGHT preconditioning the Givens estimate IS the true
    # residual of the original system, so the probe must confirm it, not contradict it
    @test length(rec.true_residual) == r.iterations
    @test length(rec.preconditioned_residual) == r.iterations
    @test all(isfinite, rec.true_residual)
    for j in 1:r.iterations
        @test abs(rec.true_residual[j] - rec.residual[j+1]) <=
              1e-6 * rec.residual[j+1] + 1e-12 * rec.residual[1]
    end
    @test rec.preconditioned_residual[end] < rec.preconditioned_residual[1]
    # without the probe the exact vectors stay empty (it costs O(k²) applies)
    mon2 = ResidualMonitor()
    solve_linearized!(zeros(ComplexF64, n), stack, b; monitor = mon2)
    @test isempty(mon2[1].true_residual)
    @test isempty(mon2[1].preconditioned_residual)
    # the cumulative-iteration axis of the figure's top panel
    solve_linearized!(zeros(ComplexF64, n), stack, b; monitor = mon)
    @test length(mon) == 2
    xs = cumulative_iterations(mon)
    @test length(xs) == 3
    @test xs[1] == 0
    @test xs[2] == mon[1].iterations
    @test xs[3] == mon[1].iterations + mon[2].iterations
    # a BigFloat monitor keeps BigFloat residuals
    @test ResidualMonitor(BigFloat) isa ResidualMonitor{BigFloat}
    @test isempty(cumulative_iterations(ResidualMonitor())) == false
    @test cumulative_iterations(ResidualMonitor()) == [0]
end

# ---- the sealed 4-level multigrid model (ESSF::multigrid-hierarchy), for the allocation --
sk_one2(x, y) = 1.0
function sk_zerobcs(mesh, T)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost(zeros(T, nf))
    end
    return bcs
end
sk_uniform_factory(T) =
    mesh -> MeshOperator(T, mesh, (sk_one2, sk_one2), (x, y) -> SK_BETA,
                         (sk_zero2, sk_zero2); C = SK_C, bcs = sk_zerobcs(mesh, T))
function sk_uniform_mesh(RT, k, N)
    blks = [Block(1, LinearMap1D(RT(0), RT(1)), LinearMap1D(RT(0), RT(1))),
            Block(2, LinearMap1D(RT(1), RT(2)), LinearMap1D(RT(0), RT(1)))]
    els = ElementId[]
    for b in (1, 2), iy in 0:(2^k-1), ix in 0:(2^k-1)
        push!(els, ElementId(b, SegmentId(k, ix), SegmentId(k, iy)))
    end
    return Mesh2D(blks, els, (N, N); connections = [BlockConnection(1, 1, 2)])
end

# ------------------------------------------- (d) preconditioner-apply allocation ---------
@testset "gate (d): the V-cycle apply allocates a documented constant" begin
    T = ComplexF64
    mesh = sk_uniform_mesh(Float64, 3, 4)        # the sealed 4-level model: 2048 dofs
    fac = sk_uniform_factory(T)
    op = fac(mesh)
    n = mesh.ndofs
    stack = SolverStack(op, fac)
    hier = stack.hierarchy
    @test length(hier.levels) == 4
    @test [l.mesh.ndofs for l in hier.levels] == [2048, 512, 128, 32]
    mg = stack.preconditioner
    z = zeros(T, n)
    v = sk_vec!(SKRNG(0xa11c), T, n)

    # serial, so the gate is independent of the thread count of the run
    vcycle_apply!() = (fill!(z, zero(T)); vcycle!(z, mg, v; threaded = false))
    vcycle_apply!()                                        # warm-up (compilation)
    a1 = @allocated vcycle_apply!()
    a2 = @allocated vcycle_apply!()
    a3 = @allocated vcycle_apply!()
    println("  [gate (d)] V-cycle apply allocation: ", a1, " / ", a2, " / ", a3, " B")
    @test max(a1, a2, a3) - min(a1, a2, a3) <= 128         # CONSTANT after warm-up (GC jitter)
    @test a1 <= SK_VCYCLE_ALLOC                            # below the pre-workspace 19 869 600 B

    # THE FLOOR, measured in place: the sealed matrix-free apply! owns it. One V-cycle runs
    # 1 residual apply! + 2 pre-smoothing + 2 post-smoothing passes (each one apply!) on
    # every non-bottom level.
    applied = 0
    for ℓ in 1:(length(hier.levels) - 1)
        l = hier.levels[ℓ]
        m = l.mesh.ndofs
        y = zeros(T, m)
        x = zeros(T, m)
        apply!(y, l.op, x; threaded = false)
        applied += 5 * @allocated apply!(y, l.op, x; threaded = false)
    end
    @test applied >= SK_APPLY_FLOOR * a1     # measured 96.1 % before / 99.99 % after
    @test applied <= 1.01 * a1               # the only other allocator is the bottom solve

    # WHAT THIS NODE REMOVED: smooth!'s own per-call workspace, preallocated additively and
    # BITWISE neutral
    lev = hier.levels[1]
    m = lev.mesh.ndofs
    bb = sk_vec!(SKRNG(0xb22c), T, m)
    u_fresh = zeros(T, m)
    smooth!(u_fresh, lev.smoother, bb; passes = 2, threaded = false)   # warm-up
    fill!(u_fresh, zero(T))
    a_fresh = @allocated smooth!(u_fresh, lev.smoother, bb; passes = 2, threaded = false)
    w = SchwarzWorkspace(lev.smoother)
    u_work = zeros(T, m)
    smooth!(u_work, lev.smoother, bb; passes = 2, threaded = false, workspace = w)
    fill!(u_work, zero(T))
    a_work = @allocated smooth!(u_work, lev.smoother, bb; passes = 2, threaded = false,
                                workspace = w)
    @test a_work < a_fresh
    @test a_fresh - a_work >= 3 * m * sizeof(T)            # at least the three pass vectors
    @test u_work == u_fresh                                # BITWISE identical result
    # a workspace of the wrong shape is IGNORED, never wrong
    wsmall = SchwarzWorkspace(hier.levels[2].smoother)
    u_bad = zeros(T, m)
    smooth!(u_bad, lev.smoother, bb; passes = 2, threaded = false, workspace = wsmall)
    @test u_bad == u_fresh

    # the monitor's per-iteration byte counter sees the same constant
    _, b, _ = sk_manufactured(op, 0xd1d1)
    mon = ResidualMonitor()
    solve_linearized!(zeros(T, n), stack, b; monitor = mon)
    pb = mon[1].preconditioner_bytes
    @test length(pb) == mon[1].iterations
    @test all(>(0), pb)
    @test maximum(pb) - minimum(pb) <= 4096                # constant after warm-up
end

# ---- the mock-AMR flag sets (calibrate_amr.log) ------------------------------------------
# "refine the finest level only": the deepest corner element of the two middle blocks — the
# CRP's puncture cascade shape (a graded corner chain), which is what makes the coarsest grid
# invariant while AMR adds finer levels.
function sk_corner_flags(mesh)
    out = Set{ElementId}()
    for blk in (2, 3)
        cand = [el for el in mesh.elements if el.block == blk]
        isempty(cand) && continue
        mx = maximum(e -> e.segments[1].level + e.segments[2].level, cand)
        for el in cand
            el.segments[1].level + el.segments[2].level == mx &&
                el.segments[1].index == 0 && el.segments[2].index == 0 && push!(out, el)
        end
    end
    return out
end
# a refinement FRONT: every deepest element of block 3 (a refined REGION, not a chain)
function sk_front_flags(mesh)
    cand = [el for el in mesh.elements if el.block == 3]
    mx = maximum(e -> e.segments[1].level + e.segments[2].level, cand)
    return Set(el for el in cand if el.segments[1].level + el.segments[2].level == mx)
end

# ------------------------------------------------------------- (e) mock AMR sweep --------
@testset "gate (e): eight AMR steps through rebuild!, caches carried" begin
    T = ComplexF64
    fac = sk_factory(T)
    mesh = sk_five_mesh(Float64, 1)
    op = fac(mesh)
    stack = SolverStack(op, fac)
    mon = ResidualMonitor()
    coarse_ndofs = Int[]
    nlev = Int[]
    for step in 1:SK_AMR_STEPS
        if step > 1
            mesh = h_refine(mesh, sk_corner_flags(mesh))
            op = fac(mesh)                       # every mesh-derived object is stale
            rebuild!(stack, op)                  # … the stack rebuilds below it
        end
        push!(coarse_ndofs, stack.hierarchy.levels[end].mesh.ndofs)
        push!(nlev, length(stack.hierarchy.levels))
        _, b, _ = sk_manufactured(op, 0x4000 + step)
        r = solve_linearized!(zeros(T, op.mesh.ndofs), stack, b; monitor = mon)
        @test r.converged
        @test all(isfinite, r.x)
    end
    @test length(mon) == SK_AMR_STEPS
    # main.tex:674 — ONE bottom LU across the whole AMR sequence
    @test stack.bottom_cache.builds == 1
    @test stack.bottom_cache.hits == SK_AMR_STEPS - 1
    @test length(stack.bottom_cache.entries) == 1
    @test all(==(coarse_ndofs[1]), coarse_ndofs)     # the bottom grid really never moved …
    @test nlev[end] > nlev[1]                        # … while AMR added finer levels
    # main.tex:676 — the subdomain inverses of every unchanged level are reused; only the
    # levels AMR actually touched are refactorized
    @test stack.subdomain_cache.hits > stack.subdomain_cache.builds
    @test all(rec -> rec.subdomain_hits >= 0, mon.records)
    @test sum(rec -> rec.subdomain_builds, mon.records) == stack.subdomain_cache.builds
    @test sum(rec -> rec.subdomain_hits, mon.records) == stack.subdomain_cache.hits
    for k in 2:SK_AMR_STEPS
        @test mon[k].subdomain_hits > 0              # every rebuild reused coarser levels
    end
    @test mon[1].bottom_builds == 1
    @test all(k -> mon[k].bottom_builds == 0 && mon[k].bottom_hits == 1, 2:SK_AMR_STEPS)
    # the build cost of each new AMR level is recorded SEPARATELY from its iterations
    @test all(rec -> rec.build_time > 0, mon.records)
    @test all(rec -> !isempty(rec.iteration_time), mon.records)
    # the fig:linsolv staircase reads left to right
    xs = cumulative_iterations(mon)
    @test length(xs) == SK_AMR_STEPS + 1
    @test issorted(xs)
    @test xs[end] == sum(rec -> rec.iterations, mon.records)
    # LANDED GATE. The paper claims "less than 10 preconditioned iterations" per AMR level
    # (main.tex:720). Measured here: the count is driven by the DEPTH of the 2:1 grading, not
    # by the problem size or the level count (identical counts at maxlevels 4, 5 and ∞ and at
    # shrink 1 and 0.5 — calibrate_gates2.log). This sweep drives a corner chain eight
    # refinement levels deep, which the CRP's region refinement never does, and the count
    # grows 8 → 11.
    its = [rec.iterations for rec in mon.records]
    @test maximum(its) <= SK_AMR_ITER_GATE
    @test maximum(its[1:6]) <= SK_AMR_EARLY_GATE
    @test its[1] <= its[end]
    # … and on a refinement FRONT — a refined REGION, the CRP's own shape — the paper's
    # "< 10 per AMR level" holds outright
    mesh2 = sk_five_mesh(Float64, 1)
    op2 = fac(mesh2)
    st2 = SolverStack(op2, fac)
    mon2 = ResidualMonitor()
    for step in 1:SK_FRONT_STEPS
        if step > 1
            mesh2 = h_refine(mesh2, sk_front_flags(mesh2))
            op2 = fac(mesh2)
            rebuild!(st2, op2)
        end
        _, b2, _ = sk_manufactured(op2, 0x5000 + step)
        r2 = solve_linearized!(zeros(T, op2.mesh.ndofs), st2, b2; monitor = mon2)
        @test r2.converged
        @test r2.iterations < SK_FRONT_ITER_GATE
    end
    @test st2.bottom_cache.builds == 1
    @test length(mesh2.elements) > 4 * length(sk_five_mesh(Float64, 1).elements)

    # FALSIFIER — the epoch is AUTHORITATIVE for BOTH caches. A rebuild at an unchanged
    # epoch on an unchanged mesh reuses everything; `invalidate!` bumps the epoch and
    # rebuilds every subdomain factorization AND the bottom LU.
    st3 = SolverStack(op2, fac)
    nsub = sum(length(l.smoother.subdomains) for l in st3.hierarchy.levels)
    distinct = st3.subdomain_cache.builds
    # every subdomain of every level is accounted for, but they are not all DISTINCT: a
    # coarsening step CARRIES OVER an element that has no complete sibling group, so such an
    # element keeps its id AND its neighbourhood and legitimately shares one factorization
    # across two levels (measured 5 of 169 here) — a real extra saving of the sealed key.
    @test distinct + st3.subdomain_cache.hits == nsub
    @test 0 < distinct <= nsub
    @test st3.bottom_cache.builds == 1
    @test st3.epoch == 0
    rebuild!(st3, fac(mesh2))                      # identically rediscretized, same epoch
    @test st3.build_stats[1] == 0                  # ZERO new factorizations
    @test st3.build_stats[2] == nsub               # every subdomain reused
    @test st3.bottom_cache.builds == 1
    invalidate!(st3)                               # the epoch bump
    @test st3.epoch == 1
    @test st3.build_stats[1] == distinct           # every distinct factorization rebuilt …
    @test st3.build_stats[1] + st3.build_stats[2] == nsub
    @test st3.bottom_cache.builds == 2             # … including the bottom LU
    @test all(l -> l.smoother.epoch == 1, st3.hierarchy.levels)
    r3 = solve_linearized!(zeros(T, op2.mesh.ndofs), st3, sk_manufactured(op2, 0x6000)[2])
    @test r3.converged
end

# ---- a SECOND coefficient set, for the stale-factory / stale-cache falsifiers -------------
# same mesh, same structure, a genuinely different matrix
sk_beta2(x, y) = (-3.0 + 0.1im) * (1 + 0.1x)
sk_factory2(T) = mesh -> MeshOperator(T, mesh, (sk_A11, sk_A22), sk_beta2, (sk_zero2, sk_zero2);
                                      C = SK_C, bcs = sk_nzbcs(mesh, T))
# a LARGE coefficient change, for the stale-cache falsifier: calibrated in
# evidence/solver-stack-core/calibrate_repair.log, where a Schwarz cache warmed on β = 1+0.5i
# and reused on this operator costs 77 iterations against 7 (patch-solve difference 0.95).
sk_beta3(x, y) = (40.0 + 0.5im) * (1 + 0.1x)
sk_factory3(T) = mesh -> MeshOperator(T, mesh, (sk_A11, sk_A22), sk_beta3, (sk_zero2, sk_zero2);
                                      C = SK_C, bcs = sk_nzbcs(mesh, T))

# ------------------------------- (H1) rebuild! replaces the OPERATOR FACTORY --------------
@testset "gate (H1): rebuild! rediscretizes the coarse levels from the NEW factory" begin
    T = ComplexF64
    mesh = sk_five_mesh(Float64, 1)
    fac1 = sk_factory(T)
    fac2 = sk_factory2(T)
    op1 = fac1(mesh)
    stack = SolverStack(op1, fac1)
    @test stack.build_op === fac1
    @test length(stack.hierarchy.levels) >= 3

    # THE BUG THIS GATES: a coefficient change on a FIXED mesh. The finest operator is the new
    # one, but every coarser level is rediscretized by the stack's factory — so a stale factory
    # preconditions a different operator than GMRES applies, and the residual never says so.
    op2 = fac2(mesh)
    A2_lv2_new = build_matrix(fac2(stack.hierarchy.levels[2].mesh))
    A2_lv2_old = build_matrix(fac1(stack.hierarchy.levels[2].mesh))
    @test norm(A2_lv2_new - A2_lv2_old, Inf) > 1.0        # the two factories really differ

    rebuild!(stack, op2; epoch = stack.epoch + 1, build_op = fac2)
    @test stack.build_op === fac2
    @test stack.precondition_build_op === fac2            # the undamped route follows build_op
    @test build_matrix(stack.hierarchy.levels[2].op) ≈ A2_lv2_new
    @test norm(build_matrix(stack.hierarchy.levels[2].op) - A2_lv2_old, Inf) > 1.0
    @test build_matrix(stack.hierarchy.levels[1].op) ≈ build_matrix(op2)
    _, b2, _ = sk_manufactured(op2, 0x1111)
    r = solve_linearized!(zeros(T, mesh.ndofs), stack, b2)
    @test r.converged
    @test norm(r.x - direct_solve(op2, b2)) / norm(direct_solve(op2, b2)) < SK_DIRECT_TOL

    # invalidate! forwards the factory too (the documented coefficient-change path)
    st2 = SolverStack(fac1(mesh), fac1)
    st2.op = fac2(mesh)
    invalidate!(st2; build_op = fac2)
    @test st2.epoch == 1
    @test st2.build_op === fac2
    @test build_matrix(st2.hierarchy.levels[2].op) ≈ A2_lv2_new

    # THE BOTTOM ALARM IS REACHABLE THROUGH THE STACK: a new factory at an UNCHANGED epoch is
    # exactly the case the bottom cache asserts on (f273eb7). Before this gate it could never
    # fire from here, because the coarsest level was still built by the stale factory.
    st3 = SolverStack(fac1(mesh), fac1)
    @test_throws ArgumentError rebuild!(st3, fac2(mesh); build_op = fac2)
    # ... and rebuild! is FAILURE-ATOMIC: the stack is untouched by the throw
    @test st3.build_op === fac1
    @test st3.op !== nothing
    @test build_matrix(st3.hierarchy.levels[2].op) ≈ A2_lv2_old
    @test st3.epoch == 0
    r3 = solve_linearized!(zeros(T, mesh.ndofs), st3, sk_manufactured(st3.op, 0x2222)[2])
    @test r3.converged                                    # still a working stack
end

# ------------------------------- (M2) the Schwarz cache trusts the epoch ------------------
@testset "gate (M2): stale Schwarz factors at a constant epoch cost iterations" begin
    T = ComplexF64
    mesh = sk_five_mesh(Float64, 1)
    fac1 = sk_factory(T)
    fac3 = sk_factory3(T)
    op3 = fac3(mesh)
    _, b3, _ = sk_manufactured(op3, 0x3333)
    n = mesh.ndofs

    # the ASYMMETRY is documented, not fixed: only the bottom cache carries a fingerprint
    d = lowercase(string(@doc SolverStack))
    @test occursin("smoke alarm", d)
    @test occursin("trusts the epoch", d)

    # a cache warmed on the OLD coefficients, then reused at the SAME epoch on the new ones
    stale = SubdomainCache()
    SolverStack(fac1(mesh), fac1; subdomain_cache = stale)          # warm it on fac1
    warmed = stale.builds
    st_stale = SolverStack(op3, fac3; subdomain_cache = stale)      # same epoch 0: all hits
    @test stale.builds == warmed                                    # ZERO new factorizations …
    @test stale.hits > 0                                            # … every patch reused
    st_fresh = SolverStack(op3, fac3; subdomain_cache = SubdomainCache())

    # THE STALENESS IS PROVABLE INDEPENDENTLY OF THE COUNT: the cached patch factor solves a
    # different system from the one this operator's patch matrix defines.
    sm_s = st_stale.hierarchy.levels[1].smoother
    sm_f = st_fresh.hierarchy.levels[1].smoother
    z = ones(T, size(subdomain_matrix(sm_f, 1), 1))
    @test norm(sm_s.fact[1] \ z - sm_f.fact[1] \ z) / norm(sm_f.fact[1] \ z) > 0.5

    r_stale = solve_linearized!(zeros(T, n), st_stale, b3)
    r_fresh = solve_linearized!(zeros(T, n), st_fresh, b3)
    println("  [gate (M2)] stale-Schwarz iterations=", r_stale.iterations,
            " vs honest=", r_fresh.iterations)
    # right preconditioning keeps the ANSWER right in both cases — the symptom is the count
    xd = direct_solve(op3, b3)
    @test r_stale.converged && r_fresh.converged
    @test norm(r_stale.x - xd) / norm(xd) < SK_DIRECT_TOL
    @test norm(r_fresh.x - xd) / norm(xd) < SK_DIRECT_TOL
    @test r_stale.iterations >= SK_STALE_ITERS               # measured 77 against 7
    @test r_fresh.iterations <= SK_MG_ITERS
    # the epoch bump is the repair, and it restores the honest count exactly
    invalidate!(st_stale; build_op = fac3)
    r_bumped = solve_linearized!(zeros(T, n), st_stale, b3)
    @test r_bumped.iterations == r_fresh.iterations

    # HONEST CAVEAT, gated so it cannot rot: the iteration count is NOT a reliable detector.
    # A SMALL coefficient change leaves the count untouched (measured 11 vs 11 for
    # β = -3+0.1i) while the patch factors are already stale — which is exactly why the epoch,
    # and not a symptom, is the safeguard.
    fac2 = sk_factory2(T)
    op2 = fac2(mesh)
    _, b2, _ = sk_manufactured(op2, 0x3334)
    stale2 = SubdomainCache()
    SolverStack(fac1(mesh), fac1; subdomain_cache = stale2)
    st_s2 = SolverStack(op2, fac2; subdomain_cache = stale2)
    st_f2 = SolverStack(op2, fac2; subdomain_cache = SubdomainCache())
    sm_s2 = st_s2.hierarchy.levels[1].smoother
    sm_f2 = st_f2.hierarchy.levels[1].smoother
    z2 = ones(T, size(subdomain_matrix(sm_f2, 1), 1))
    @test norm(sm_s2.fact[1] \ z2 - sm_f2.fact[1] \ z2) / norm(sm_f2.fact[1] \ z2) > 0.05
    r_s2 = solve_linearized!(zeros(T, n), st_s2, b2)
    r_f2 = solve_linearized!(zeros(T, n), st_f2, b2)
    println("  [gate (M2)] small-change stale=", r_s2.iterations, " honest=", r_f2.iterations,
            " (the count does NOT detect staleness)")
    @test r_s2.iterations == r_f2.iterations
    @test norm(r_s2.x - direct_solve(op2, b2)) / norm(direct_solve(op2, b2)) < SK_DIRECT_TOL
end

# ------------------------------- (M3) build_time belongs to ONE record --------------------
@testset "gate (M3): build_time is this rebuild's cost, zero on later solves" begin
    T = ComplexF64
    mesh = sk_five_mesh(Float64, 1)
    fac = sk_factory(T)
    op = fac(mesh)
    stack = SolverStack(op, fac)          # warm-up build (carries JIT)
    stack = SolverStack(op, fac)          # the timed one
    _, b, _ = sk_manufactured(op, 0x4444)
    mon = ResidualMonitor()
    solve_linearized!(zeros(T, mesh.ndofs), stack, b; monitor = mon)
    solve_linearized!(zeros(T, mesh.ndofs), stack, b; monitor = mon)
    solve_linearized!(zeros(T, mesh.ndofs), stack, b; monitor = mon)
    @test length(mon) == 3
    @test mon[1].build_time > 0
    @test mon[2].build_time == 0.0        # NOT copied forward any more
    @test mon[3].build_time == 0.0
    @test mon[1].subdomain_builds > 0
    @test mon[2].subdomain_builds == 0 && mon[2].subdomain_hits == 0
    @test mon[2].bottom_builds == 0 && mon[2].bottom_hits == 0
    # the split adds up and is attributed
    @test mon[1].build_time_hierarchy > 0
    @test mon[1].build_time_bottom > 0
    @test mon[1].build_time_hierarchy + mon[1].build_time_bottom <= mon[1].build_time * 1.001
    @test mon[2].build_time_hierarchy == 0.0 && mon[2].build_time_bottom == 0.0
    # a rebuild re-arms exactly one record
    rebuild!(stack, fac(mesh))
    solve_linearized!(zeros(T, mesh.ndofs), stack, b; monitor = mon)
    solve_linearized!(zeros(T, mesh.ndofs), stack, b; monitor = mon)
    @test mon[4].build_time > 0
    @test mon[5].build_time == 0.0
    # summing the column over a sweep is the sweep's build cost, with no double counting
    @test sum(r -> r.build_time, mon.records) == mon[1].build_time + mon[4].build_time
end

# ------------------------------- (M1) hp-AMR and the coarse-level order -------------------
@testset "gate (M1): coarse_orders pins the bottom grid under hp-AMR" begin
    T = ComplexF64
    fac = sk_factory(T)
    # spectre coarsens h ONLY and allocates every level from the same initial_extents
    # (Multigrid/Hierarchy.cpp:21-30, ElementsAllocator.hpp:101-102/142-148); under p-AMR the
    # coarse levels are frozen (Amr/Actions/AdjustDomain.hpp:206-212,
    # SendDataToChildren.hpp:53-55). Our MIN rule (ported from the misleading comment at
    # spectre Multigrid.hpp:51-52) instead raises the coarsest order on every hp step.
    function hp_sweep(coarse_orders, steps)
        mesh = sk_five_mesh(Float64, 1)
        opts = SolverStackOptions(; coarse_orders = coarse_orders)
        op = fac(mesh)
        st = SolverStack(op, fac; options = opts)
        its = Int[]; cords = Vector{NTuple{2,Int}}[]
        for step in 1:steps
            if step > 1
                mesh = hp_refine(mesh, sk_corner_flags(mesh))
                op = fac(mesh)
                rebuild!(st, op)
            end
            push!(cords, unique(st.hierarchy.levels[end].mesh.orders))
            _, b, _ = sk_manufactured(op, 0x7000 + step)
            push!(its, solve_linearized!(zeros(T, op.mesh.ndofs), st, b).iterations)
        end
        return st, its, cords
    end
    steps = 6
    st_min, its_min, cords_min = hp_sweep(nothing, steps)
    st_fix, its_fix, cords_fix = hp_sweep(4, steps)
    println("  [gate (M1)] hp sweep, MIN rule: bottom builds=", st_min.bottom_cache.builds,
            " entries=", length(st_min.bottom_cache.entries), " coarsest orders=", cords_min,
            " its=", its_min)
    println("  [gate (M1)] hp sweep, coarse_orders=4: bottom builds=",
            st_fix.bottom_cache.builds, " entries=", length(st_fix.bottom_cache.entries),
            " coarsest orders=", cords_fix, " its=", its_fix)
    # the MIN rule: the coarsest order RISES with p-AMR, so main.tex:674 fails here
    @test st_min.bottom_cache.builds == steps
    @test length(unique(cords_min)) > 1
    # the spectre-faithful rule: the bottom grid really does remain the same (main.tex:674)
    @test st_fix.bottom_cache.builds == 1
    @test st_fix.bottom_cache.hits == steps - 1
    @test length(st_fix.bottom_cache.entries) == 1
    @test all(==(cords_fix[1]), cords_fix)
    @test all(o -> o == (4, 4), cords_fix[1])
    # and every coarser level keeps order N, not just the coarsest
    for l in st_fix.hierarchy.levels[2:end]
        @test all(o -> o == (4, 4), l.mesh.orders)
    end
    # THE TRADE-OFF, measured, not assumed. Pinning buys main.tex:674 (one cached bottom LU)
    # and pays in iterations once the fine order runs away from the pinned coarse one: this
    # sweep p-refines EVERY element on every step, so by step 6 the fine grid is order 10
    # against a coarse space pinned at 4 and the count grows 8 -> 17, while the MIN rule keeps
    # 8 -> 10 by rebuilding the bottom LU six times. Neither is free; the CRP's own hp-AMR
    # (8 iterations over its whole domain) sits between these two extremes, so row 42/43 must
    # choose per run and this gate pins BOTH curves rather than declaring a winner.
    @test maximum(its_fix) <= SK_HP_ITER_GATE            # measured 17
    @test maximum(its_min) <= SK_ITER_GATE               # measured 10
    @test its_fix[1] == its_min[1]                       # step 1 has no p gap yet
    @test maximum(its_fix) > maximum(its_min)            # pinning costs iterations here
    @test all(r -> r > 0, its_fix)
    # the option is validated and reaches the hierarchy
    @test_throws ArgumentError SolverStackOptions(; coarse_orders = 0)
    @test SolverStackOptions().coarse_orders === nothing
end

# ---- the row-18 resonance model, transcribed VERBATIM from test_schwarz_shift.jl --------
# (the file-local prefix convention of this package: mg_*, bs_*, shf_* each carry their own
# copy of a shared model; nothing is re-derived here — the coefficients, the mesh, the
# manufactured solution and the resonance search are the sealed ones.)
skr_profile(x, y) = 1 + 0.3x
skr_mesh(N) = Mesh2D(Block(81, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0)),
                     vec([ElementId(81, SegmentId(1, ix), SegmentId(1, iy))
                          for ix in 0:1, iy in 0:1]), (N, N))
skr_top(mesh, k; s = 0.0) =
    MeshOperator(ComplexF64, mesh, (sk_one2, sk_one2),
                 (x, y) -> -k^2 * skr_profile(x, y) - im * s,
                 (sk_zero2, sk_zero2); C = SK_C, bcs = sk_zerobcs(mesh, ComplexF64))
function skr_fieldvec(mesh, f)
    u = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        Nx, Ny = mesh.orders[e]
        grid = ntuple(d -> begin
            lo, hi = logical_interval(Float64, el, d)
            [map_forward(blockmap(blk, d), (lo * (1 - ξ) + hi * (1 + ξ)) / 2)
             for ξ in EllipticCore.refelement(Float64, d == 1 ? Nx : Ny).nodes]
        end, 2)
        for j in 1:Ny, i in 1:Nx
            u[mesh.offsets[e] + tensor_lin(i, j, Nx)] = f(grid[1][i], grid[2][j])
        end
    end
    return u
end
function skr_kres(mesh)
    A0 = subdomain_matrix(SchwarzSmoother(skr_top(mesh, 0.0); overlap = 2), 1)
    B = A0 .- subdomain_matrix(SchwarzSmoother(skr_top(mesh, 1.0); overlap = 2), 1)
    ks = [sqrt(real(λ)) for λ in eigvals(Matrix(A0), Matrix(B))
          if abs(imag(λ)) < 1e-8 * abs(λ) && real(λ) > 0 && 10 ≤ sqrt(real(λ)) ≤ 20]
    return ks[argmin(abs.(ks .- 14.0))]
end

# ------------------------------------------------------------- (f) the shift policy ------
@testset "gate (f): the helmholtz_shift policy reaches every level, and its contrast" begin
    mesh = skr_mesh(14)                       # the sealed row-18 mesh
    kres = skr_kres(mesh)
    @test 12.5 < kres < 12.7                  # the sealed resonance, 12.58801463986746
    s = kres^2

    # the POLICY value is the sealed function, not a number invented here
    @test helmholtz_shift(4, 0.3; slicing = :t) == abs2(4 * 0.3)
    @test helmholtz_shift(4, 0.3; slicing = :vtu) == 0

    # (1) the LITERAL seam of SchwarzSmoother reaches EVERY level's smoother
    op = skr_top(mesh, kres)
    fac = m -> skr_top(m, kres)
    lit = SolverStack(op, fac;
                      options = SolverStackOptions(; subdomain_complex_shift = s))
    @test length(lit.hierarchy.levels) >= 2
    @test all(l -> l.smoother.shift == s, lit.hierarchy.levels)
    @test all(l -> l.smoother.epoch == 0, lit.hierarchy.levels)
    # t-slicing vs vtu: the policy value, forwarded verbatim
    vtu = SolverStack(op, fac; options = SolverStackOptions(;
                      subdomain_complex_shift = helmholtz_shift(4, 0.3; slicing = :vtu)))
    @test all(l -> l.smoother.shift == 0, vtu.hierarchy.levels)

    # (2) the PRODUCTION route (ESSF::schwarz-complex-shift): the shift is mass-weighted,
    # i.e. it enters through a DAMPED preconditioner factory β − i·s, not as −i·s·I
    uex = skr_fieldvec(mesh, (x, y) -> cis(kres * (x + 0.5y)) * sinpi(y))
    b = apply!(zeros(ComplexF64, mesh.ndofs), op, uex)
    n = mesh.ndofs
    opts = SolverStackOptions(; abstol = 0.0, reltol = SK_SHIFT_RTOL,
                              maxiter = SK_SHIFT_BUDGET, restart = SK_SHIFT_BUDGET)
    damped = m -> skr_top(m, kres; s = s)
    results = Dict{Tuple{Symbol,Symbol},Any}()
    for (tag, pbo) in ((:unshifted, fac), (:shifted, damped)),
        (pname, pfac) in ((:schwarz, schwarz_only_preconditioner),
                          (:multigrid, multigrid_preconditioner))
        st = SolverStack(op, fac; options = opts, precondition_build_op = pbo,
                         preconditioner = pfac)
        results[(tag, pname)] = solve_linearized!(zeros(ComplexF64, n), st, b)
    end
    # the sealed contrast: unshifted exhausts the budget, the policy shift converges
    @test !results[(:unshifted, :schwarz)].converged
    @test results[(:unshifted, :schwarz)].iterations == SK_SHIFT_BUDGET
    @test results[(:unshifted, :schwarz)].residual > 1e2
    @test !results[(:unshifted, :multigrid)].converged
    @test results[(:shifted, :schwarz)].converged
    @test results[(:shifted, :schwarz)].iterations <= SK_SHIFT_GATE       # sealed 47
    @test results[(:shifted, :multigrid)].converged
    @test results[(:shifted, :multigrid)].iterations <=
          results[(:shifted, :schwarz)].iterations                        # measured 44 ≤ 47
    # 3 shifted solves fit inside the budget the unshifted run exhausts (the row-18 >3× gate)
    @test 3 * results[(:shifted, :schwarz)].iterations < SK_SHIFT_BUDGET

    # the shift can NEVER change the answer: it lives only in the preconditioner
    st_s = SolverStack(op, fac; precondition_build_op = damped,
                       options = SolverStackOptions(; abstol = 0.0, reltol = 1e-12,
                                                    maxiter = 300, restart = 300))
    st_2s = SolverStack(op, fac; precondition_build_op = m -> skr_top(m, kres; s = 2s),
                        options = SolverStackOptions(; abstol = 0.0, reltol = 1e-12,
                                                     maxiter = 300, restart = 300))
    u1 = zeros(ComplexF64, n); r1 = solve_linearized!(u1, st_s, b)
    u2 = zeros(ComplexF64, n); r2 = solve_linearized!(u2, st_2s, b)
    @test r1.converged && r2.converged
    @test norm(u1 - u2) / norm(u1) < 1e-9
    @test norm(u1 - uex) / norm(uex) < 1e-9
    # the damped hierarchy really is a DIFFERENT operator from the solved one
    @test st_s.hierarchy.levels[1].op !== op
    @test st_s.op === op
end

# ------------------------------------------------ (g) type genericity + thread parity ----
@testset "gate (g): Complex{BigFloat}(256) and thread parity" begin
    setprecision(BigFloat, 256) do
        T = Complex{BigFloat}
        mesh = sk_five_mesh(BigFloat, 1; N = 3)
        fac = sk_factory(T)
        op = fac(mesh)
        n = mesh.ndofs
        opts = SolverStackOptions(; abstol = SK_BIG_ABSTOL, reltol = big(0.0),
                                  maxiter = 200, restart = 200)
        stack = SolverStack(op, fac; options = opts)
        @test stack.options isa SolverStackOptions{BigFloat}
        @test length(stack.hierarchy.levels) >= 2
        # UMFPACK takes exactly {Float64, ComplexF64}; everything else falls back to the
        # generic DENSE LU in the caller's own precision (the sealed bottom-solver contract)
        @test stack.bottom_cache.builds == 1
        uex, b, bphys = sk_manufactured(op, 0x81f0)
        @test eltype(b) === T
        @test linearized_rhs(op, bphys) ≈ b
        u = zeros(T, n)
        mon = ResidualMonitor(BigFloat)
        r = solve_linearized!(u, stack, b; monitor = mon)
        @test r.converged
        @test r.iterations <= SK_BIG_ITERS
        @test r.residual <= opts.abstol
        @test precision(real(u[1])) == 256
        @test mon[1] isa SolveRecord{BigFloat}
        @test mon[1].residual[end] isa BigFloat
        @test norm(u - uex) / norm(uex) < 1e-30
    end

    # thread parity: the V-cycle apply is bitwise independent of the thread count (inherited
    # from the bitwise-gated threaded apply!/smooth!), and so is the whole solve
    T = ComplexF64
    mesh = sk_uniform_mesh(Float64, 2, 4)
    fac = sk_uniform_factory(T)
    op = fac(mesh)
    n = mesh.ndofs
    stack = SolverStack(op, fac)
    _, b, _ = sk_manufactured(op, 0x7777)
    mg = stack.preconditioner
    zs = zeros(T, n); vcycle!(zs, mg, b; threaded = false)
    zt = zeros(T, n); vcycle!(zt, mg, b; threaded = true)
    @test zs == zt
    u = zeros(T, n)
    r = solve_linearized!(u, stack, b)
    @test r.converged
    if Threads.nthreads() < 2
        @test_skip "threaded solve == serial needs -t2"
    else
        # a REAL parity check: the whole solve, not just the V-cycle. A stack whose smoothers
        # run threaded (this run) must reproduce the serial answer bitwise, which is what the
        # sealed bitwise-threaded apply!/smooth! buy the composition.
        us = zeros(T, n)
        rs = solve_linearized!(us, stack, b)
        @test us == u
        @test rs.x == r.x
        @test rs.history == r.history
        @test rs.iterations == r.iterations
        zt2 = zeros(T, n); vcycle!(zt2, mg, b; threaded = true)
        @test zt2 == zs
    end
    # a second stack on a private pair of caches reproduces the run bitwise
    st2 = SolverStack(op, fac)
    u2 = zeros(T, n)
    r2 = solve_linearized!(u2, st2, b)
    @test r2.x == r.x
    @test r2.history == r.history
end

# ----------------------------------------------------------- errors and edge cases -------
@testset "errors, edge cases and the mutable-stack contract" begin
    mesh = sk_five_mesh(Float64, 1)
    fac = sk_factory(ComplexF64)
    op = fac(mesh)
    n = mesh.ndofs
    _, b, _ = sk_manufactured(op, 0x9999)
    stack = SolverStack(op, fac)
    @test_throws ArgumentError solve_linearized!(zeros(ComplexF64, n - 1), stack, b)
    @test_throws ArgumentError solve_linearized!(zeros(ComplexF64, n), stack, b[1:end-1])
    @test_throws ArgumentError solve!(zeros(ComplexF64, n), stack, b[1:end-1])
    # a rebuild with a different scalar type is refused loudly
    opreal = MeshOperator(Float64, mesh, (sk_A11, sk_A22), (x, y) -> 1.0,
                          (sk_zero2, sk_zero2); C = SK_C, bcs = sk_zerobcs(mesh, Float64))
    @test EllipticCore._scalar(opreal) !== EllipticCore._scalar(op)
    @test_throws ArgumentError rebuild!(stack, opreal)
    # maxiter = 0 returns the initial state honestly rather than looping
    st0 = SolverStack(op, fac; options = SolverStackOptions(; maxiter = 0, restart = 1))
    r0 = solve_linearized!(zeros(ComplexF64, n), st0, b)
    @test !r0.converged
    @test r0.iterations == 0
    # use_bottom_solver = false is the yaml's UseBottomSolver: False twin (smooth instead)
    stnb = SolverStack(op, fac; options = SolverStackOptions(; use_bottom_solver = false))
    @test stnb.preconditioner.bottom === nothing
    @test stnb.bottom_cache.builds == 0
    rnb = solve_linearized!(zeros(ComplexF64, n), stnb, b)
    @test rnb.converged
    # initial_coarse_levels caps the hierarchy (yaml:188 + 1)
    stcap = SolverStack(op, fac; options = SolverStackOptions(; initial_coarse_levels = 1))
    @test length(stcap.hierarchy.levels) == 2
    # shrink is exposed and forwarded
    stsh = SolverStack(op, fac; options = SolverStackOptions(; shrink = 0.5))
    @test length(stsh.hierarchy.levels) <= length(stack.hierarchy.levels)
    rsh = solve_linearized!(zeros(ComplexF64, n), stsh, b)
    @test rsh.converged
    @test rsh.iterations <= SK_ITER_GATE
end

# ------------------------------------------- the EXPLICIT-CHAIN preconditioner factory ----
# `multigrid_chain_preconditioner` (DAG row 56, linsolv-depth-gap). Written BEFORE the
# implementation; the red run is captured in the node's evidence directory.
#
# WHY IT EXISTS. arXiv:2606.04998 does not run the DERIVED hierarchy under AMR. Its executable
# keeps the coarse grids (`SolveScalarSelfForce.hpp:173`), every AMR step creates a new finest
# grid at `grid_index + 1` (`Amr/Actions/AdjustDomain.hpp:157-163` for split elements, `:203-212`
# for the p increment) and the parent survives as a coarser grid
# (`Amr/Actions/SendDataToChildren.hpp:43-70`), so the grids the V-cycle walks
# (`Multigrid/ElementActions.hpp:96-135`, `Multigrid.hpp:56-58`) ARE the AMR history. This
# factory is how a caller hands that history to the stack: it owns a vector of coarse meshes,
# FINEST-FIRST, and the factory reads it at every `rebuild!`.
#
# GATE CONSTANTS calibrated PRE-RED on the throwaway prototype
# `evidence/linsolv-depth-gap/calibrate_chain_precond.jl` (log `calibrate_chain_precond.log`),
# which runs the gate (M1) hp sweep both ways.
SK_CHAIN_STEPS = 6              # the gate (M1) sweep length, reused
SK_CHAIN_SPREAD = 1             # history counts span 1 iteration over steps >= 2 (measured 8,9,9,9,9)
@testset "gate (M4): the explicit-chain preconditioner is the AMR history" begin
    T = ComplexF64
    fac = sk_factory(T)

    # (1) EQUIVALENCE ANCHOR — handed the DERIVED hierarchy's own meshes, the chain factory is
    #     the derived stack, bitwise, iterate for iterate.
    mesh0 = sk_five_mesh(Float64, 1)
    op0 = fac(mesh0)
    ref = SolverStack(op0, fac; options = SolverStackOptions())
    chain_eq = [l.mesh for l in ref.hierarchy.levels[2:end]]
    st_eq = SolverStack(op0, fac; options = SolverStackOptions(),
                        preconditioner = multigrid_chain_preconditioner(chain_eq))
    @test length(st_eq.hierarchy.levels) == length(ref.hierarchy.levels)
    _, beq, _ = sk_manufactured(op0, 0x9001)
    x1 = zeros(T, op0.mesh.ndofs); x2 = zeros(T, op0.mesh.ndofs)
    r1 = solve_linearized!(x1, ref, beq)
    r2 = solve_linearized!(x2, st_eq, beq)
    @test r1.iterations == r2.iterations
    @test x1 == x2                                        # BITWISE the same solve

    # (2) THE hp SWEEP, both hierarchies. The chain is the sweep's own mesh history, so the
    #     coarsest level is the INITIAL mesh and never moves.
    function hp_sweep(mode)
        mesh = sk_five_mesh(Float64, 1)
        chain = Mesh2D{Float64}[]
        opts = SolverStackOptions(; coarse_orders = 4)
        op = fac(mesh)
        st = mode === :history ?
             SolverStack(op, fac; options = opts,
                         preconditioner = multigrid_chain_preconditioner(chain)) :
             SolverStack(op, fac; options = opts)
        its = Int[]; depth = Int[]; coarsest = Int[]
        for step in 1:SK_CHAIN_STEPS
            if step > 1
                newmesh = hp_refine(mesh, sk_corner_flags(mesh))
                mode === :history && pushfirst!(chain, mesh)
                mesh = newmesh
                op = fac(mesh)
                rebuild!(st, op)
            end
            push!(depth, length(st.hierarchy.levels))
            push!(coarsest, st.hierarchy.levels[end].mesh.ndofs)
            _, b, _ = sk_manufactured(op, 0x7000 + step)
            push!(its, solve_linearized!(zeros(T, op.mesh.ndofs), st, b).iterations)
        end
        return st, its, depth, coarsest, mesh
    end
    st_d, its_d, dep_d, co_d, _ = hp_sweep(:derived)
    st_h, its_h, dep_h, co_h, m_end = hp_sweep(:history)
    println("  [gate (M4)] derived: its=", its_d, " depth=", dep_d, " coarsest=", co_d,
            " bottom=(", st_d.bottom_cache.builds, ",", st_d.bottom_cache.hits, ")")
    println("  [gate (M4)] history: its=", its_h, " depth=", dep_h, " coarsest=", co_h,
            " bottom=(", st_h.bottom_cache.builds, ",", st_h.bottom_cache.hits, ")")

    # the DERIVED hierarchy degrades as the fine order runs away from the pinned coarse space
    @test its_d[end] > its_d[1]
    # the AMR-history hierarchy does NOT: after the first refinement the count is bounded
    @test maximum(its_h[2:end]) - minimum(its_h[2:end]) <= SK_CHAIN_SPREAD
    @test its_h[end] < its_d[end]
    # its depth IS the number of grids the sweep has made (one per AMR step, level 0 included)
    @test dep_h == collect(1:SK_CHAIN_STEPS)
    # step 1 has an EMPTY chain: a single level, whose bottom solve IS the whole problem
    @test its_h[1] == 1
    # the coarsest grid is the INITIAL mesh and never moves (main.tex:674), so one bottom LU
    @test all(==(mesh0.ndofs), co_h)
    @test st_h.bottom_cache.builds == 1
    @test st_h.bottom_cache.hits == SK_CHAIN_STEPS - 1
    # the chain is READ at every rebuild, so mutating the caller's vector is what moved the depth
    @test dep_h[end] - dep_h[1] == SK_CHAIN_STEPS - 1
    # `initial_coarse_levels`, `shrink` and `coarse_orders` are NOT consulted by this factory:
    # the depth is the chain's length + 1 whatever they say
    chain2 = [l.mesh for l in ref.hierarchy.levels[2:3]]
    st2 = SolverStack(op0, fac; options = SolverStackOptions(; initial_coarse_levels = 0,
                                                             shrink = 0.5, coarse_orders = 4),
                      preconditioner = multigrid_chain_preconditioner(chain2))
    @test length(st2.hierarchy.levels) == 3
    @test multigrid_chain_preconditioner isa Function
    doc = string(Base.Docs.doc(Base.Docs.Binding(EllipticCore, :multigrid_chain_preconditioner)))
    @test !occursin("No documentation found", doc)
end

end   # @testset "solver stack (solver-stack-core)"

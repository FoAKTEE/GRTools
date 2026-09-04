# solver_stack.jl — the CRP solver-stack COMPOSITION API (main.tex:668-676, fig:precond
# caption main.tex:632): right-preconditioned complex GMRES ∘ one multigrid V-cycle ∘ two
# additive-Schwarz passes per level ∘ a cached direct bottom solve, built ONCE per
# (mesh, operator) and rebuilt across AMR steps with every cache the sealed contracts allow,
# plus the residual monitor that records the per-AMR-level GMRES staircase and its wall clock
# (the fig:linsolv data surface, main.tex:701).
# DAG node: solver-stack-core (elliptic-ssf GLOBAL_DAG row 21). Built TDD.
#
# THIS FILE COMPOSES; IT DOES NOT RE-IMPLEMENT. Every numerical step is a sealed predecessor:
#   gmres.jl        `gmres(applyA!, b; Minv!, …)`   — ESSF::complex-gmres
#   multigrid.jl    `MultigridHierarchy`, `MultigridV` — ESSF::multigrid-hierarchy
#   bottom_solve.jl `sparse_bottom_solver`          — ESSF::multigrid-bottom-solver
#   schwarz.jl      `SchwarzSmoother`, `SubdomainCache` — ESSF::schwarz-smoother
#   schwarz_shift.jl `helmholtz_shift`, `schwarz_preconditioner` — ESSF::schwarz-complex-shift
# The adopted preconditioner seam is the `Minv!` CALLABLE of `gmres` (DESIGN.md §4.3 landed
# deviation): there is no `Preconditioner`/`ldiv!` type hierarchy and none is introduced here.
# The composition is gated BITWISE against the hand-assembled call in test_solver_stack.jl.
#
# The `preconditioner` seam of `SolverStack` is a plain FACTORY CALLABLE `f(stack) -> (Minv!,
# hierarchy_or_nothing)`; three are shipped (`multigrid_preconditioner`,
# `schwarz_only_preconditioner`, `no_preconditioner`) and a caller may pass its own.

export SolverStackOptions, SolverStack, SolveRecord, ResidualMonitor,
       solve!, solve_linearized!, rebuild!, invalidate!, linearized_rhs,
       cumulative_iterations,
       multigrid_preconditioner, multigrid_chain_preconditioner,
       schwarz_only_preconditioner, no_preconditioner

# ---- options ----------------------------------------------------------------------------

"""
    SolverStackOptions(; abstol = 1e-12, reltol = 0, maxiter = 500, restart = maxiter,
                         multigrid_iterations = 1, initial_coarse_levels = nothing,
                         pre_smoothing = true, use_bottom_solver = true,
                         post_smoothing_at_bottom = false, shrink = 1,
                         overlap = 2, passes = 2, subdomain_complex_shift = 0)

Plain configuration record for [`SolverStack`](@ref). Every field mirrors a key of the CRP's
own solver input `anc/ScalarSelfForce.yaml:178-205`, under the yaml's NAME, and every default
is the yaml's production value — so the `ScalarSelfForce` driver fills this struct field by
field from the sealed `RunConfig` (`GmresConfig`, `MultigridConfig`, `SchwarzConfig`) without
`EllipticCore` ever depending on that package (DESIGN.md §1: this core stratum carries zero
physics and no capability dependency).

`LinearSolver.Gmres` (yaml:178-184)
  - `abstol` — `ConvergenceCriteria.AbsoluteResidual` (yaml:183), default `1e-12`.
  - `reltol` — `ConvergenceCriteria.RelativeResidual` (yaml:182), default `0`.
  - `maxiter` — `ConvergenceCriteria.MaxIterations` (yaml:181), default `500`.
  - `restart` — the Arnoldi-basis cap per restart cycle. No yaml twin (spectre's GMRES is
    full-memory); the default `maxiter` reproduces unrestarted GMRES exactly.

`LinearSolver.Multigrid` (yaml:186-194)
  - `multigrid_iterations` — `Iterations` (yaml:187), default `1` V-cycle per GMRES iteration.
  - `initial_coarse_levels` — `InitialCoarseLevels` (yaml:188), the number of COARSER levels;
    `nothing` (the default) coarsens all the way to the irreducible block structure. The
    hierarchy's `maxlevels` counts TOTAL levels, so the twin is `initial_coarse_levels + 1`
    ([`MultigridHierarchy`](@ref)).
  - `pre_smoothing` — `PreSmoothing` (yaml:189), default `true`.
  - `use_bottom_solver` — `UseBottomSolver` (yaml:190), default `true`; `false` smooths at the
    bottom instead ([`MultigridV`](@ref)'s `bottom = nothing`).
  - `post_smoothing_at_bottom` — `PostSmoothingAtBottom` (yaml:191), default `false`.
  - `shrink` — the level policy of [`MultigridHierarchy`](@ref); no yaml twin. Default `1`
    (never skip a coarsening step) — see the `shrink` paragraph of [`SolverStack`](@ref) for
    the measured tradeoff table that fixed that default.
  - `coarse_orders` — the coarse-level polynomial policy of [`MultigridHierarchy`](@ref); no
    yaml twin. Default `nothing` (the MIN rule). Set it to the run's initial `P` for an
    `hp_refine` sweep: spectre's multigrid coarsens h ONLY and allocates every level from the
    same `initial_extents`, so its coarse levels keep a fixed order while p-AMR raises the fine
    grid — which is what makes main.tex:674 ("the bottom grid remains the same even when AMR
    adds finer levels") true. Under the MIN rule an hp step raises the coarsest order too and
    the cached bottom LU is rebuilt every step (measured 6/6). It is a TRADE, not a free win:
    on a sweep that p-refines every element six times the pinned coarse space falls six
    orders behind the fine grid and the count grows 8 -> 17, against 8 -> 10 for the MIN
    rule that rebuilds the bottom LU each step. Choose per run; both curves are gated.
    See the `coarse_orders`
    paragraph of [`MultigridHierarchy`](@ref) for the spectre file:line citations.
  - the yaml's multigrid-level `ComplexShift` (yaml:194) has NO field: the complex shift is
    subdomain-level only in this package (adopted decision 8) and the CRP's production value
    is `0.0` anyway.

`LinearSolver.SchwarzSmoother` (yaml:196-205)
  - `passes` — `Iterations` (yaml:197), default `2` smoothing passes per level.
  - `overlap` — the subdomain overlap in grid points. The yaml says `MaxOverlap: Auto`
    (yaml:198); `2` is this package's [`SchwarzSmoother`](@ref) default.
  - `subdomain_complex_shift` — `SubdomainSolver.ExplicitInverse.ComplexShift` (yaml:203),
    default `0`, forwarded to every level's smoother as the LITERAL `−i·s·I` seam. READ the
    shift paragraph of [`SolverStack`](@ref) before setting it: on this package's massive Gauss
    formulation the literal seam is NOT the policy route.

Numeric fields promote to one real type `RT`; `BigFloat` tolerances stay `BigFloat`.
"""
struct SolverStackOptions{RT<:AbstractFloat}
    abstol::RT
    reltol::RT
    maxiter::Int
    restart::Int
    multigrid_iterations::Int
    initial_coarse_levels::Union{Nothing,Int}
    pre_smoothing::Bool
    use_bottom_solver::Bool
    post_smoothing_at_bottom::Bool
    shrink::RT
    coarse_orders::Union{Nothing,Int}
    overlap::Int
    passes::Int
    subdomain_complex_shift::RT
end

function SolverStackOptions(; abstol::Real = 1e-12, reltol::Real = 0,
                            maxiter::Integer = 500, restart::Integer = maxiter,
                            multigrid_iterations::Integer = 1,
                            initial_coarse_levels::Union{Nothing,Integer} = nothing,
                            pre_smoothing::Bool = true, use_bottom_solver::Bool = true,
                            post_smoothing_at_bottom::Bool = false, shrink::Real = 1,
                            coarse_orders::Union{Nothing,Integer} = nothing,
                            overlap::Integer = 2, passes::Integer = 2,
                            subdomain_complex_shift::Real = 0)
    RT = float(promote_type(typeof(abstol), typeof(reltol), typeof(shrink),
                            typeof(subdomain_complex_shift)))
    abstol >= 0 || throw(ArgumentError("SolverStackOptions: abstol must be ≥ 0, got $abstol"))
    reltol >= 0 || throw(ArgumentError("SolverStackOptions: reltol must be ≥ 0, got $reltol"))
    maxiter >= 0 ||
        throw(ArgumentError("SolverStackOptions: maxiter must be ≥ 0, got $maxiter"))
    restart >= 1 ||
        throw(ArgumentError("SolverStackOptions: restart must be ≥ 1, got $restart"))
    multigrid_iterations >= 1 ||
        throw(ArgumentError("SolverStackOptions: multigrid_iterations must be ≥ 1, " *
                            "got $multigrid_iterations"))
    initial_coarse_levels === nothing || initial_coarse_levels >= 0 ||
        throw(ArgumentError("SolverStackOptions: initial_coarse_levels must be ≥ 0 " *
                            "(it counts the COARSER levels), got $initial_coarse_levels"))
    0 < shrink <= 1 ||
        throw(ArgumentError("SolverStackOptions: shrink must be in (0, 1], got $shrink"))
    overlap >= 0 ||
        throw(ArgumentError("SolverStackOptions: overlap must be ≥ 0, got $overlap"))
    passes >= 1 ||
        throw(ArgumentError("SolverStackOptions: passes must be ≥ 1, got $passes"))
    coarse_orders === nothing || coarse_orders >= 1 ||
        throw(ArgumentError("SolverStackOptions: coarse_orders must be ≥ 1, " *
                            "got $coarse_orders"))
    icl = initial_coarse_levels === nothing ? nothing : Int(initial_coarse_levels)
    cor = coarse_orders === nothing ? nothing : Int(coarse_orders)
    return SolverStackOptions{RT}(RT(abstol), RT(reltol), Int(maxiter), Int(restart),
                                  Int(multigrid_iterations), icl, pre_smoothing,
                                  use_bottom_solver, post_smoothing_at_bottom, RT(shrink),
                                  cor, Int(overlap), Int(passes),
                                  RT(subdomain_complex_shift))
end

# yaml InitialCoarseLevels (COARSER levels) → MultigridHierarchy maxlevels (TOTAL levels)
_ss_maxlevels(o::SolverStackOptions) =
    o.initial_coarse_levels === nothing ? nothing : o.initial_coarse_levels + 1

# ---- the monitor ------------------------------------------------------------------------

"""
    SolveRecord{RT}

What one [`solve!`](@ref) through a [`SolverStack`](@ref) recorded — one entry of a
[`ResidualMonitor`](@ref), and exactly the data surface `fig:linsolv` plots (main.tex:701:
"the GMRES iterations needed to solve consecutive AMR levels … the bottom panel shows the
walltime for each iteration. The first iteration on each new AMR level is slower than the
following because it builds and caches subdomain matrices").

Residuals (all in `RT = real(T)`):
  - `residual` — the sealed [`GMRESResult`](@ref)`.history`: `residual[1] = ‖b − A x₀‖` then
    one Givens estimate per inner iteration (`length == iterations + 1`). With RIGHT
    preconditioning this estimate is the residual of the ORIGINAL system, so this vector IS
    the figure's staircase.
  - `true_residual`, `preconditioned_residual` — per iteration `‖b − A xⱼ‖` and `‖M⁻¹(b −
    A xⱼ)‖`, RECOMPUTED, not estimated. EMPTY unless `solve!` was called with
    `true_residual_probe = true`: the sealed `gmres` does not expose its intermediate
    iterates, so the probe re-solves with `maxiter = 1, 2, …` and costs O(k²) applies. It is a
    diagnostic (and the gate that the Givens estimate does not drift), never the production path.

Wall clock (seconds, `Float64` at every scalar type):
  - `iteration_time[j]` — from the start of iteration `j`'s preconditioner apply to the start
    of iteration `j+1`'s; the LAST entry runs to the end of the solve and therefore also
    carries the back-substitution and the final residual recomputation.
  - `preconditioner_time[j]`, `operator_time[j]` — the V-cycle apply and the `A·z` apply
    inside iteration `j`; `iteration_time - (preconditioner_time + operator_time)` is the
    Arnoldi/Givens work.
  - `preconditioner_bytes[j]` — bytes allocated by iteration `j`'s preconditioner apply,
    from `Base.gc_bytes()` differences. That counter is PROCESS-WIDE, so on a threaded run or
    beside any concurrent allocation it is an upper bound, not an isolated measurement; take
    the per-apply figure from a quiet serial run (as the allocation gate does) and read this
    field as a monitor, not as a benchmark.

Build side (the "first iteration on each new AMR level is slower" story — see the
`build_time` paragraph of [`ResidualMonitor`](@ref)):
  - `build_time` — wall clock of THIS solve's own preceding [`rebuild!`](@ref): the factory
    call that built the hierarchy, its Schwarz factorizations and the bottom solve. It is
    charged to exactly ONE record — the first solve after that rebuild — and is `0.0` on every
    later solve over the same preconditioner, so summing the column over an AMR sweep gives
    the sweep's true build cost with no double counting. `build_time_hierarchy` and
    `build_time_bottom` split it (hierarchy construction, which includes every Schwarz patch
    factorization, versus the coarsest-level assembly and LU); they are `0.0` for a
    preconditioner factory that builds neither.
    `subdomain_builds`/`subdomain_hits`/`bottom_builds`/`bottom_hits` are that same rebuild's
    cache-counter deltas — how many subdomain factorizations and bottom LUs it actually
    constructed versus reused — and are likewise zero on the later records.
    WARM UP BEFORE TIMING: the first rebuild in a process carries JIT (measured 9.9 s against
    0.25-0.67 s warm on the same mesh).
  - `level_ndofs`, `level_subdomains` — per multigrid level (finest first); empty for a
    preconditioner without a hierarchy.

Summary: `iterations`, `converged`, `final_residual` (the recomputed true residual) and
`solve_time`.
"""
struct SolveRecord{RT<:Real}
    residual::Vector{RT}
    true_residual::Vector{RT}
    preconditioned_residual::Vector{RT}
    iteration_time::Vector{Float64}
    preconditioner_time::Vector{Float64}
    operator_time::Vector{Float64}
    preconditioner_bytes::Vector{Int}
    iterations::Int
    converged::Bool
    final_residual::RT
    solve_time::Float64
    build_time::Float64
    build_time_hierarchy::Float64
    build_time_bottom::Float64
    level_ndofs::Vector{Int}
    level_subdomains::Vector{Int}
    subdomain_builds::Int
    subdomain_hits::Int
    bottom_builds::Int
    bottom_hits::Int
end

"""
    ResidualMonitor()            # Float64 residuals
    ResidualMonitor(RT)          # e.g. ResidualMonitor(BigFloat)

Recorder for a sequence of [`solve!`](@ref)s through one [`SolverStack`](@ref) — one
[`SolveRecord`](@ref) appended per solve, in call order. Across an AMR sweep that is one
record per AMR level, which is precisely `fig:linsolv`: the top panel plots each record's
`residual` against a cumulative-iteration axis ([`cumulative_iterations`](@ref) gives the
level boundaries), the bottom panel plots each record's `iteration_time`.

THE FIRST-ITERATION BUILD COST IS RECORDED SEPARATELY, not inside iteration 1. The CRP builds
its subdomain matrices lazily, so its first iteration on a new AMR level absorbs that cost
(main.tex:701); this package builds them EAGERLY in [`SchwarzSmoother`](@ref) at hierarchy
construction, so the cost lands in [`rebuild!`](@ref) and is reported as `build_time` together
with the `subdomain_builds`/`bottom_builds` counters that explain it. A caller reproducing the
figure adds `build_time` to the level's first iteration; a caller measuring the SOLVER reads
`iteration_time` alone.

Ported idea (not code) from the SpecCore `elliptic/solvers/residual_monitor.jl` recorder,
itself modelled on spectre's `LinearSolver/Gmres/ResidualMonitor.hpp` reduction: a plain
struct filled synchronously by the driving solver instead of a parallel component.

The recorder is passive: `solve!(u, stack, b; monitor = mon)` appends; nothing else mutates it.
"""
struct ResidualMonitor{RT<:Real}
    records::Vector{SolveRecord{RT}}
end
ResidualMonitor() = ResidualMonitor{Float64}(SolveRecord{Float64}[])
ResidualMonitor(::Type{RT}) where {RT<:Real} = ResidualMonitor{RT}(SolveRecord{RT}[])

Base.length(mon::ResidualMonitor) = length(mon.records)
Base.getindex(mon::ResidualMonitor, i) = mon.records[i]
Base.lastindex(mon::ResidualMonitor) = lastindex(mon.records)
Base.isempty(mon::ResidualMonitor) = isempty(mon.records)

"""
    cumulative_iterations(mon::ResidualMonitor) -> Vector{Int}

The cumulative-inner-iteration x-positions at which the recorded solves start:
`[0, m₁, m₁+m₂, …]` with `mᵢ = mon.records[i].iterations`. The `fig:linsolv` top panel draws
record `i`'s `residual` starting at `cumulative_iterations(mon)[i]`, which is what makes
consecutive AMR levels read as one left-to-right staircase (main.tex:701, "the lines in the
top panel, from left to right, show the GMRES iterations needed to solve consecutive AMR
levels"). `length == length(mon) + 1`; the last entry is the total iteration count.
"""
function cumulative_iterations(mon::ResidualMonitor)
    xs = Int[0]
    acc = 0
    for r in mon.records
        acc += r.iterations
        push!(xs, acc)
    end
    return xs
end

# ---- the stack --------------------------------------------------------------------------

"""
    SolverStack(op::MeshOperator, build_op; options = SolverStackOptions(),
                preconditioner = multigrid_preconditioner, precondition_build_op = build_op,
                subdomain_cache = SubdomainCache(), bottom_cache = BottomSolverCache(),
                epoch = 0)

The CRP linear-solver stack (main.tex:668-676, `fig:precond`): right-preconditioned complex
GMRES whose preconditioner is ONE multigrid V-cycle, with two additive-Schwarz passes on every
finer level and a direct, cached solve at the bottom. Built ONCE per (mesh, operator) and
carried across AMR steps with [`rebuild!`](@ref); solved with [`solve!`](@ref).

    stack = SolverStack(op, build_op)              # production defaults = the yaml
    res   = solve!(u, stack, b_phys)               # a GMRESResult

`op` is the finest [`MeshOperator`](@ref) — the matrix GMRES applies. `build_op(mesh) ->
MeshOperator` is the caller's operator FACTORY, rediscretizing the SAME problem on any mesh it
is handed with that mesh's own boundary conditions (the [`MultigridHierarchy`](@ref) contract).
`options` is the yaml twin ([`SolverStackOptions`](@ref)).

# THE LINEARIZED-RHS CONTRACT (normative — read before using this API)

Everything below GMRES is LINEARIZED-only. [`apply!`](@ref) is the linear-homogeneous matrix
(fixed hook data and inhomogeneous boundary data are NOT in it), the Schwarz smoother and the
V-cycle inherit that, and the multigrid coarse levels solve for a CORRECTION from a zero guess,
so no fixed data may enter them. The right-hand side therefore has to carry the
[`operator_on_zero`](@ref) lift — and [`solve!`](@ref) APPLIES THAT LIFT FOR YOU:

    solve!(u, stack, b_phys)                        # folds internally; the safe default
    solve_linearized!(u, stack, linearized_rhs(op, b_phys))   # the pre-folded route

The two are bitwise the same run (gated). The fold costs one full-operator apply, about 1 % of
a single V-cycle. [`linearized_rhs`](@ref) is the ONE place fixed data is folded; the failure it
guards against is not loud — a physical right-hand side handed to the linearized solver
CONVERGES, `converged = true`, on a different problem (measured 3.4e3 relative error on the
5-block gate model), which is why the safe route is the default and the pre-folded one carries
its own name instead of a flag. Folding TWICE (handing `linearized_rhs(op, b)` to `solve!`) is
an O(1) wrong answer and is the caller's error; both directions are gated.

# THE TWO-CACHES-ONE-EPOCH CONTRACT (normative)

The stack holds exactly TWO caller-owned caches, and they are what makes main.tex:674-676 true
("the LU decomposition of A^bottom is cached and reused in all multigrid V-cycles across AMR
iterations … also here, the matrix inverses can remain cached across AMR iterations, so we only
have to construct new subdomain matrices when AMR adds a finer level"):

  - `subdomain_cache::SubdomainCache` — ONE shared cache for every level's Schwarz
    factorizations. Its key `(element id, subdomain shape, epoch, shift)` cannot collide across
    levels or AMR steps, so an unchanged coarser level is all hits.
  - `bottom_cache::BottomSolverCache` — the coarsest-level LU.

ONE `epoch` GOVERNS BOTH, AND IT IS AUTHORITATIVE. No cache can see a caller's coefficient,
hook, boundary-condition or penalty closures, so the caller — not the cache — decides when a
factorization is stale. BUMP THE EPOCH ON ANY MATRIX-AFFECTING CHANGE (coefficients, hooks,
boundary conditions, penalty, shift) through [`invalidate!`](@ref) or
[`rebuild!`](@ref)`(stack, op; epoch = e)`. A pure MESH change (an AMR step) needs no bump: the
mesh is in both keys. Never rely on a cache "noticing" a coefficient change — the bottom
cache's operator fingerprint is a CONSISTENCY ASSERTION, not a detector: at an unchanged epoch
with a matching structural key but a different operator it raises, naming this rule.

ONLY THE BOTTOM CACHE CARRIES A SMOKE ALARM; THE SCHWARZ CACHE TRUSTS THE EPOCH. The bottom
cache stores an operator fingerprint beside each factorization and raises when a structural-key
hit at an unchanged epoch meets a different operator; the Schwarz cache has no such check by
design — one fingerprint per patch would cost an operator apply per subdomain, and there are
hundreds — so at a constant epoch a coefficient change silently reuses stale patch factors.
Right preconditioning keeps the ANSWER correct in that case; what degrades is the preconditioner,
and the symptom is an iteration count, not a wrong result (gated as a stack-level falsifier: the
stale-Schwarz solve costs iterations and the epoch bump restores them). The epoch is the
safeguard — bump it.

Both caches GROW without bound and are never evicted: the caller owns their lifetime. Drop the
stack (or hand `rebuild!` a fresh pair) to release the factorizations of an AMR sequence.
NEITHER CACHE IS THREAD-SAFE, and neither is the stack: one stack, one thread at a time —
exactly the [`MultigridV`](@ref) rule. Pass your own caches to share them across stacks of the
same operator family; the defaults are fresh, private, and still reused across this stack's own
rebuilds. Measured on the 5-block AMR sweep: eight refinement steps of the finest level build
the bottom LU ONCE (seven hits) and reuse every subdomain of every unchanged level; an
`invalidate!` rebuilds all of them.

# THE PRECONDITIONER SEAM

The stack is MUTABLE and its callable fields (`build_op`, `precondition_build_op`,
`make_preconditioner`, `preconditioner`) are untyped, because [`rebuild!`](@ref) replaces them:
the per-solve cost is one dynamic call per GMRES iteration against a V-cycle, which is noise.

`preconditioner` is a plain factory callable `f(stack) -> (Minv!, hierarchy_or_nothing)` — the
`Minv!` right-preconditioner hook of [`gmres`](@ref) (DESIGN.md §4.3: no `Preconditioner`/
`ldiv!` type hierarchy exists in this package and none is wanted). Three are shipped and any
caller may pass another: [`multigrid_preconditioner`](@ref) (the default and the CRP stack),
[`schwarz_only_preconditioner`](@ref) (the smoother alone) and [`no_preconditioner`](@ref)
(the identity). Measured on the 5-block gate model, the three order as the paper's stack
predicts: multigrid < Schwarz-only < none.

# THE COMPLEX SHIFT

`options.subdomain_complex_shift` is forwarded to every level's smoother as the LITERAL
`−i·s·I` seam of [`SchwarzSmoother`](@ref). On this package's MASSIVE Gauss formulation that is
NOT the policy route: `ESSF::schwarz-complex-shift` measured that an operator-level damping
`β → β − i·s` discretizes to `−i·s·M` with the diagonal Gauss mass `M ~ 1e-3`, so the literal
unweighted seam at the policy value `s = helmholtz_shift(m, Ω)` degenerates every patch solve
and FAILS. The production route is `precondition_build_op`: hand the stack a second factory
that rediscretizes the DAMPED operator `β − i·s`, and the whole preconditioner (levels,
smoothers and bottom solve) is built on it while GMRES keeps applying the true `op`. Right
preconditioning leaves the residual, the convergence test and the converged answer those of the
original system, so the shift can never change the answer — only the iteration count.

    s     = helmholtz_shift(m, Ω; slicing = :t)      # (mΩ)²; 0 under vtu/hyperboloidal
    stack = SolverStack(op, build_op;
                        precondition_build_op = mesh -> damped_op(mesh, s))

# THE `shrink` DEFAULT

`options.shrink` is [`MultigridHierarchy`](@ref)'s level policy and defaults to `1` (every
coarsening step becomes a level). Measured here on the three model families, at GMRES abstol
`1e-12` on a unit-norm right-hand side — the stack's figure of merit is the ITERATION COUNT,
not the stationary contraction:

| model                        | shrink | levels | Σndofs/ndofs | ρ      | GMRES its | t_solve |
|------------------------------|--------|--------|--------------|--------|-----------|---------|
| uniform 8×8 two-block        | 1 … ¼  | 4      | 1.328        | 0.1716 | 8         | 0.183 s |
| graded puncture base 2, 4 rd | 1      | 12     | 3.393        | 0.1731 | 9         | 0.090 s |
| graded puncture base 2, 4 rd | ½      | 5      | 1.786        | 0.1731 | 9         | 0.057 s |
| graded puncture base 3, 5 rd | 1      | 16     | 2.671        | 0.1702 | 9         | 0.416 s |
| graded puncture base 3, 5 rd | ½      | 6      | 1.620        | 0.1724 | 9         | 0.172 s |
| anisotropic 44 el            | 1      | 8      | 2.182        | 0.1679 | 9         | 0.097 s |
| anisotropic 44 el            | ½      | 5      | 1.727        | 0.2218 | 9         | 0.086 s |
| anisotropic 157 el two-block | 1      | 10     | 2.446        | 0.1970 | 9         | 0.523 s |
| anisotropic 157 el two-block | ½      | 6      | 1.554        | 0.2973 | 11        | 0.402 s |

Reading: on a UNIFORM mesh `shrink` does nothing (each step already halves the problem). On the
GRADED puncture cascade — the CRP's own refinement pattern — `shrink = 0.5` cuts the hierarchy
cost ratio to ≈1.6-1.8 and the V-cycle to less than half the wall clock at an IDENTICAL
iteration count, and it is the right choice there. On ANISOTROPIC meshes it costs contraction
(0.168 → 0.222, 0.197 → 0.297) and, on the larger of them, iterations (9 → 11) — past the
paper's own "less than 10 preconditioned iterations per AMR level" (main.tex:720). `shrink = 0.25`
loses on every family (ρ up to 0.50, 11-13 iterations). The CRP's AMR policy is
`Isotropy: Anisotropic` (`anc/ScalarSelfForce.yaml:145`), so the DEFAULT stays the conservative
`1`, which never loses an iteration on any family measured; `0.5` is exposed and documented for
callers whose meshes are graded.

Fields: `op`, `build_op`, `precondition_build_op`, `options`, `subdomain_cache`,
`bottom_cache`, `epoch`, `preconditioner` (the live `Minv!`), `hierarchy` (the live
[`MultigridHierarchy`](@ref) or `nothing`), `build_time` with its `build_time_hierarchy` /
`build_time_bottom` split and `build_stats` (the last [`rebuild!`](@ref)'s wall clock and its
`(subdomain builds, hits, bottom builds, hits)` deltas), and `build_pending` (whether the next
[`solve!`](@ref) still owes that rebuild's cost to its [`SolveRecord`](@ref)). The stack is
MUTABLE and single-application at a time, exactly like the [`MultigridV`](@ref) it holds:
build one per thread.
"""
mutable struct SolverStack{RT<:AbstractFloat}
    op::MeshOperator
    build_op::Any
    precondition_build_op::Any
    options::SolverStackOptions{RT}
    make_preconditioner::Any
    subdomain_cache::SubdomainCache
    bottom_cache::BottomSolverCache
    epoch::Int
    preconditioner::Any
    hierarchy::Any
    build_time::Float64
    build_time_hierarchy::Float64
    build_time_bottom::Float64
    build_stats::NTuple{4,Int}
    build_pending::Bool
end

function SolverStack(op::MeshOperator, build_op;
                     options::SolverStackOptions = SolverStackOptions(),
                     preconditioner = multigrid_preconditioner,
                     precondition_build_op = build_op,
                     subdomain_cache::SubdomainCache = SubdomainCache(),
                     bottom_cache::BottomSolverCache = BottomSolverCache(),
                     epoch::Integer = 0)
    RT = typeof(options).parameters[1]
    stack = SolverStack{RT}(op, build_op, precondition_build_op, options, preconditioner,
                            subdomain_cache, bottom_cache, Int(epoch), nothing, nothing,
                            0.0, 0.0, 0.0, (0, 0, 0, 0), false)
    return rebuild!(stack, op; epoch = epoch)
end

"""
    rebuild!(stack::SolverStack, op = stack.op; epoch = stack.epoch,
             build_op = stack.build_op, precondition_build_op = …) -> stack

Rebuild the stack's preconditioner. This is BOTH the AMR step (a new mesh, same factory) and
the coefficient step (the same mesh, a NEW factory) — and the two are not interchangeable:

  - AMR: `rebuild!(stack, build_op(new_mesh))`. Both caches are reused, so an unchanged coarser
    level costs zero factorizations and the bottom LU survives as long as the coarsest grid and
    operator do. No `epoch` bump (the mesh is in both keys).
  - COEFFICIENTS / hooks / boundary conditions / penalty / shift: hand `rebuild!` (or
    [`invalidate!`](@ref)) the NEW `build_op` **and** bump the `epoch`. THE FACTORY IS WHAT THE
    COARSE LEVELS ARE REDISCRETIZED FROM: leaving it stale silently preconditions with the old
    operator on every level below the finest, which the residual never reveals because right
    preconditioning keeps the answer correct — it only costs iterations. The stack therefore
    replaces `build_op` here rather than capturing it once at construction.

`precondition_build_op` defaults to `build_op` when the stack was UNDAMPED
(`precondition_build_op === build_op` at the previous rebuild), so the undamped route follows
the new factory automatically; a stack carrying a separate damped factory (the complex-shift
route) keeps it, and a caller whose damping changed must pass the new one explicitly.

FAILURE-ATOMIC: `op`, `epoch`, the factories, the preconditioner and the hierarchy are all
committed only after the factory returns. If it throws, the stack is left exactly as it was and
the exception propagates. The one thing that is NOT rolled back is the caller-owned caches: a
factorization built before the throw stays in them. That is safe — every entry is keyed by
`(element, shape, epoch, shift)` or by the bottom key, so it is a correct entry, merely built
early — but it means `stack.build_stats` and the cache counters may move on a failed rebuild.

Records this rebuild's wall clock in `stack.build_time` (split into `build_time_hierarchy` and
`build_time_bottom`) and its cache-counter deltas in `stack.build_stats`; the NEXT
[`solve!`](@ref) copies them into its [`SolveRecord`](@ref) and later solves record `0.0`
(see [`SolveRecord`](@ref)).

Every mesh-derived object of the previous step is stale after an AMR step (the sealed AMR
contract), so hand `rebuild!` the operator built on the NEW mesh.
"""
function rebuild!(stack::SolverStack, op::MeshOperator = stack.op;
                  epoch::Integer = stack.epoch,
                  build_op = stack.build_op,
                  precondition_build_op =
                      stack.precondition_build_op === stack.build_op ? build_op :
                      stack.precondition_build_op)
    _scalar(op) === _scalar(stack.op) ||
        throw(ArgumentError("rebuild!: the new operator has scalar type $(_scalar(op)), " *
                            "the stack was built for $(_scalar(stack.op))"))
    sc, bc = stack.subdomain_cache, stack.bottom_cache
    b0 = (sc.builds, sc.hits, bc.builds, bc.hits)
    saved = (stack.op, stack.build_op, stack.precondition_build_op, stack.epoch)
    stack.op = op
    stack.build_op = build_op
    stack.precondition_build_op = precondition_build_op
    stack.epoch = Int(epoch)
    stack.build_time_hierarchy = 0.0
    stack.build_time_bottom = 0.0
    local pc, hier, dt
    try
        t0 = time_ns()
        pc, hier = stack.make_preconditioner(stack)
        dt = (time_ns() - t0) / 1e9
    catch
        stack.op, stack.build_op, stack.precondition_build_op, stack.epoch = saved
        rethrow()
    end
    stack.preconditioner = pc
    stack.hierarchy = hier
    stack.build_time = dt
    stack.build_stats = (sc.builds - b0[1], sc.hits - b0[2], bc.builds - b0[3], bc.hits - b0[4])
    stack.build_pending = true
    return stack
end

"""
    invalidate!(stack::SolverStack; epoch = stack.epoch + 1, build_op = stack.build_op,
                precondition_build_op = …) -> stack

Bump the caches' coefficient `epoch` and rebuild — the explicit invalidator of the
two-caches-one-epoch contract ([`SolverStack`](@ref)). Call it whenever the operator's
coefficients, hooks, boundary conditions or penalty changed without the mesh changing, and pass
the NEW `build_op` along with the new finest `op`: the epoch alone invalidates the CACHES, the
factory is what the coarse levels are rediscretized from ([`rebuild!`](@ref)).

    op = new_build_op(mesh)
    invalidate!(stack; build_op = new_build_op)      # after stack.op = op, or:
    rebuild!(stack, op; epoch = stack.epoch + 1, build_op = new_build_op)
"""
invalidate!(stack::SolverStack; epoch::Integer = stack.epoch + 1,
            build_op = stack.build_op,
            precondition_build_op =
                stack.precondition_build_op === stack.build_op ? build_op :
                stack.precondition_build_op) =
    rebuild!(stack, stack.op; epoch = epoch, build_op = build_op,
             precondition_build_op = precondition_build_op)

# ---- the shipped preconditioner factories -----------------------------------------------

# the operator the PRECONDITIONER is built on: `op` itself unless the caller supplied a
# separate factory (the mass-weighted complex-shift route — see the SolverStack docstring)
_ss_precond_op(stack::SolverStack) =
    stack.precondition_build_op === stack.build_op ? stack.op :
    stack.precondition_build_op(stack.op.mesh)::MeshOperator

"""
    multigrid_preconditioner(stack::SolverStack) -> (MultigridV, MultigridHierarchy)

The CRP preconditioner and the default `preconditioner` factory of [`SolverStack`](@ref): one
[`MultigridV`](@ref) V-cycle over a [`MultigridHierarchy`](@ref) of the stack's operator, with
`options.passes` Schwarz passes per level and the coarsest level solved directly by
[`sparse_bottom_solver`](@ref) (`options.use_bottom_solver = false` smooths there instead). The
hierarchy takes the stack's shared [`SubdomainCache`](@ref) and the bottom solver its
[`BottomSolverCache`](@ref), both at the stack's `epoch`; the bottom solver is built from
`hier.levels[end].op`, never from the finest operator.
"""
function multigrid_preconditioner(stack::SolverStack)
    o = stack.options
    th = time_ns()
    hier = MultigridHierarchy(_ss_precond_op(stack), stack.precondition_build_op;
                              caches = stack.subdomain_cache, epoch = stack.epoch,
                              overlap = o.overlap, shift = o.subdomain_complex_shift,
                              maxlevels = _ss_maxlevels(o), shrink = o.shrink,
                              coarse_orders = o.coarse_orders)
    return _ss_from_hierarchy(stack, hier, th)
end

"""
    multigrid_chain_preconditioner(meshes::AbstractVector) -> factory

`preconditioner` factory for [`SolverStack`](@ref) that builds the EXPLICIT-CHAIN hierarchy
[`MultigridHierarchy`](@ref)`(op, build_op, meshes)` instead of deriving the coarse levels by
h-coarsening the current mesh. `meshes` is the CALLER'S OWN vector of coarse meshes,
FINEST-FIRST and not repeating `op.mesh`; the factory READS IT AT EVERY [`rebuild!`](@ref), so
mutating it (`pushfirst!` a mesh the sweep has just left behind) is how the hierarchy follows an
AMR sweep. An empty vector is the single-level hierarchy whose bottom solve is the whole problem.

`options.initial_coarse_levels`, `options.shrink` and `options.coarse_orders` are NOT consulted
— the caller supplies the levels — while `subdomain_complex_shift`, `passes`,
`multigrid_iterations`, `pre_smoothing`, `post_smoothing_at_bottom` and `use_bottom_solver` mean
exactly what they mean for [`multigrid_preconditioner`](@ref).

`overlap` overrides `options.overlap` with any of [`MultigridHierarchy`](@ref)'s overlap
policies — an integer, one integer per level, or a callable of the level's own mesh. USE THE
PER-LEVEL FORM WHEN THE CHAIN IS AN AMR HISTORY: the subdomain cache is keyed on
`(element, subdomain shape, epoch, shift)` and the shape carries the overlap, so a scalar taken
from the FINEST mesh changes the key of every coarser level at every rebuild and rebuilds
factorizations that main.tex:676 says are reused. A callable that returns each level's own
maximum order leaves the frozen coarse levels bit-identical — the same subdomains either way
(the smoother clamps per neighbour), only the key differs.

    chain = Mesh2D{Float64}[]                                   # coarse meshes, finest first
    stack = SolverStack(op, build_op; preconditioner = multigrid_chain_preconditioner(chain))
    …
    pushfirst!(chain, sol.mesh)                                  # leaving this mesh behind
    rebuild!(stack, build_op(finer_mesh))                        # the hierarchy grew by one

WHY THE FORM EXISTS, and what supplying an AMR history measures: see
`docs/src/core/ellipticcore.md`. The reference keeps its coarse grids, so its
multigrid levels ARE its AMR history; the measured before/after is recorded there
and in the node's evidence.
"""
multigrid_chain_preconditioner(meshes::AbstractVector; overlap = nothing) =
    function (stack::SolverStack)
        o = stack.options
        th = time_ns()
        hier = MultigridHierarchy(_ss_precond_op(stack), stack.precondition_build_op, meshes;
                                  caches = stack.subdomain_cache, epoch = stack.epoch,
                                  overlap = overlap === nothing ? o.overlap : overlap,
                                  shift = o.subdomain_complex_shift)
        return _ss_from_hierarchy(stack, hier, th)
    end

# the common tail of every hierarchy-based factory: charge the hierarchy build, take (or hit)
# the cached bottom LU through the hierarchy-aware constructor of
# ESSF::multigrid-bottom-solver — it selects `hier.levels[end].op` itself, so the coarsest
# operator can never be confused with the finest one at this call site — and wrap the V-cycle.
function _ss_from_hierarchy(stack::SolverStack, hier, t_hierarchy::UInt64)
    o = stack.options
    stack.build_time_hierarchy = (time_ns() - t_hierarchy) / 1e9
    tb = time_ns()
    bottom = o.use_bottom_solver ?
             sparse_bottom_solver(hier; cache = stack.bottom_cache, epoch = stack.epoch) :
             nothing
    stack.build_time_bottom = (time_ns() - tb) / 1e9
    mg = MultigridV(hier; passes = o.passes, iterations = o.multigrid_iterations,
                    pre_smoothing = o.pre_smoothing,
                    post_smoothing_at_bottom = o.post_smoothing_at_bottom, bottom = bottom)
    return mg, hier
end

"""
    schwarz_only_preconditioner(stack::SolverStack) -> (Minv!, nothing)

`preconditioner` factory building the additive-Schwarz smoother of the FINEST level alone
(`options.passes` passes through the sealed [`schwarz_preconditioner`](@ref) adapter) — no
hierarchy, no bottom solve. The row-18 configuration, and the middle term of the stack's
iteration-count contrast (multigrid < Schwarz-only < none).
"""
function schwarz_only_preconditioner(stack::SolverStack)
    o = stack.options
    sm = SchwarzSmoother(_ss_precond_op(stack); overlap = o.overlap,
                         shift = o.subdomain_complex_shift, epoch = stack.epoch,
                         cache = stack.subdomain_cache)
    return schwarz_preconditioner(sm; passes = o.passes), nothing
end

"""
    no_preconditioner(stack::SolverStack) -> (Minv!, nothing)

`preconditioner` factory for unpreconditioned GMRES: `Minv! = M⁻¹ = I`, written as an explicit
`copyto!` callable rather than `gmres`'s `Minv! = nothing` so that the identity stack is
instrumented by the [`ResidualMonitor`](@ref) exactly like the others. The two are BITWISE the
same run (gated): the preconditioned directions are exact copies of the Arnoldi vectors.
"""
function no_preconditioner(::SolverStack)
    return (z::AbstractVector, v::AbstractVector) -> copyto!(z, v), nothing
end

# ---- the right-hand side ----------------------------------------------------------------

"""
    linearized_rhs(op::MeshOperator, b_phys) -> b
    linearized_rhs(stack::SolverStack, b_phys) -> b

`b_phys − operator_on_zero(op)`: the ONE place that folds an operator's fixed boundary/hook
data into its right-hand side. The whole stack below GMRES is LINEARIZED-only, so this is the
form [`solve!`](@ref) requires — see the linearized-RHS contract of [`SolverStack`](@ref).
Costs one full-operator apply; call it once per operator and never fold twice.
"""
linearized_rhs(op::MeshOperator, b_phys::AbstractVector) = b_phys .- operator_on_zero(op)
linearized_rhs(stack::SolverStack, b_phys::AbstractVector) = linearized_rhs(stack.op, b_phys)

# ---- solve ------------------------------------------------------------------------------

"""
    solve!(u, stack::SolverStack, b_phys; x0 = nothing, monitor = nothing,
           true_residual_probe = false) -> GMRESResult

Solve `stack.op · u = b_phys` for the PHYSICAL right-hand side, writing the solution into `u`
and returning the sealed [`GMRESResult`](@ref) unchanged (`x`, `iterations`, `residual`,
`converged`, `breakdown`, `history`).

THE FOLD IS DONE FOR YOU. Everything below GMRES is LINEARIZED-only, so the operator's fixed
boundary/hook data has to be lifted out of the right-hand side ([`linearized_rhs`](@ref));
`solve!` does exactly that, once, and then calls [`solve_linearized!`](@ref):

    solve!(u, stack, b_phys) ≡ solve_linearized!(u, stack, linearized_rhs(stack.op, b_phys))

(bitwise, gated). It costs ONE extra full-operator apply — about 1 % of a single V-cycle — and
it removes the failure this API previously invited: a physical right-hand side handed to the
linearized solver converges happily on the WRONG system (measured 3.4e3 relative error on the
5-block gate model, `converged = true`). Callers that already hold a folded right-hand side use
[`solve_linearized!`](@ref); DOUBLE-FOLDING (passing `linearized_rhs(op, b)` to `solve!`) is an
O(1) wrong answer and is the caller's error, gated as a falsifier.

`x0` warm-starts the iteration (default zero).

`monitor::ResidualMonitor` appends one [`SolveRecord`](@ref) — see [`solve_linearized!`](@ref).
"""
function solve!(u::AbstractVector, stack::SolverStack, b_phys::AbstractVector;
                x0::Union{Nothing,AbstractVector} = nothing,
                monitor::Union{Nothing,ResidualMonitor} = nothing,
                true_residual_probe::Bool = false)
    length(b_phys) == stack.op.mesh.ndofs ||
        throw(ArgumentError("solve!: b has length $(length(b_phys)), the mesh has " *
                            "$(stack.op.mesh.ndofs) dofs"))
    return solve_linearized!(u, stack, linearized_rhs(stack.op, b_phys); x0 = x0,
                             monitor = monitor, true_residual_probe = true_residual_probe)
end

"""
    solve_linearized!(u, stack::SolverStack, b; x0 = nothing, monitor = nothing,
                      true_residual_probe = false) -> GMRESResult

The pre-folded route: solve `apply!(·, stack.op, u) = b` with right-preconditioned GMRES over
the stack's preconditioner, for a `b` that ALREADY carries the operator-on-zero correction
(`b = linearized_rhs(op, b_phys)`). Prefer [`solve!`](@ref), which folds for you; use this one
when the folded right-hand side is already in hand (an AMR sweep that reuses it, or a caller
assembling `b` in linearized form to begin with) and the extra apply is worth avoiding.

`x0` warm-starts the iteration (default zero).

`monitor::ResidualMonitor` appends one [`SolveRecord`](@ref) with the per-iteration residual
staircase, wall clock and preconditioner allocation, plus — on the FIRST solve after a
[`rebuild!`](@ref) — that rebuild's build cost and cache-counter deltas. With
`monitor = nothing` (the default) the call is literally
`gmres(stack.op, b; Minv! = stack.preconditioner, …)` and nothing is measured; instrumentation
wraps the two callables and does not touch the arithmetic, so a monitored solve is BITWISE an
unmonitored one (gated). `true_residual_probe = true` additionally fills the record's
`true_residual`/`preconditioned_residual` by re-solving with `maxiter = 1, 2, …` — exact,
O(k²) applies, diagnostic only.
"""
function solve_linearized!(u::AbstractVector, stack::SolverStack, b::AbstractVector;
                           x0::Union{Nothing,AbstractVector} = nothing,
                           monitor::Union{Nothing,ResidualMonitor} = nothing,
                           true_residual_probe::Bool = false)
    o = stack.options
    n = stack.op.mesh.ndofs
    length(u) == n ||
        throw(ArgumentError("solve_linearized!: u has length $(length(u)), the mesh has " *
                            "$n dofs"))
    length(b) == n ||
        throw(ArgumentError("solve_linearized!: b has length $(length(b)), the mesh has " *
                            "$n dofs"))
    M! = stack.preconditioner
    if monitor === nothing
        res = gmres(stack.op, b; abstol = o.abstol, reltol = o.reltol, maxiter = o.maxiter,
                    restart = o.restart, Minv! = M!, x0 = x0)
        copyto!(u, res.x)
        return res
    end
    return _ss_solve_monitored!(u, stack, b, x0, monitor, true_residual_probe)
end

function _ss_solve_monitored!(u, stack::SolverStack, b, x0, monitor::ResidualMonitor{RT},
                              probe::Bool) where {RT}
    o = stack.options
    M! = stack.preconditioner
    marks = Float64[]                 # timestamp at the start of each iteration
    ptime = Float64[]
    otime = Float64[]
    pbytes = Int[]
    inflight = Ref(false)
    wrapM = M! === nothing ? nothing : function (z::AbstractVector, v::AbstractVector)
        push!(marks, time_ns() / 1e9)
        inflight[] = true
        g0 = Base.gc_bytes()
        t0 = time_ns()
        M!(z, v)
        push!(ptime, (time_ns() - t0) / 1e9)
        push!(pbytes, Int(Base.gc_bytes() - g0))
        return z
    end
    wrapA = function (y::AbstractVector, v::AbstractVector)
        t0 = time_ns()
        apply!(y, stack.op, v)
        dt = (time_ns() - t0) / 1e9
        if inflight[]
            push!(otime, dt)
            inflight[] = false
        end
        return y
    end
    t_start = time_ns()
    res = gmres(wrapA, b; abstol = o.abstol, reltol = o.reltol, maxiter = o.maxiter,
                restart = o.restart, Minv! = wrapM, x0 = x0)
    t_solve = (time_ns() - t_start) / 1e9
    push!(marks, time_ns() / 1e9)
    copyto!(u, res.x)

    itime = [marks[j+1] - marks[j] for j in 1:(length(marks) - 1)]
    truer = RT[]
    pcr = RT[]
    probe && _ss_probe_residuals!(truer, pcr, stack, b, x0, res.iterations)
    hier = stack.hierarchy
    lnd = hier === nothing ? Int[] : [l.mesh.ndofs for l in hier.levels]
    lns = hier === nothing ? Int[] : [length(l.smoother.subdomains) for l in hier.levels]
    # THE BUILD COST BELONGS TO ONE RECORD. `rebuild!` arms `build_pending`; the first solve
    # after it takes the cost, every later solve on the same preconditioner records 0.
    bt, bth, btb, bs = stack.build_pending ?
        (stack.build_time, stack.build_time_hierarchy, stack.build_time_bottom,
         stack.build_stats) : (0.0, 0.0, 0.0, (0, 0, 0, 0))
    stack.build_pending = false
    push!(monitor.records,
          SolveRecord{RT}(RT.(res.history), truer, pcr, itime, ptime, otime, pbytes,
                          res.iterations, res.converged, RT(res.residual), t_solve,
                          bt, bth, btb, lnd, lns, bs[1], bs[2], bs[3], bs[4]))
    return res
end

# EXACT per-iteration residuals: the sealed `gmres` does not expose its intermediate iterates,
# so iterate j is reproduced by re-solving with maxiter = j (abstol = reltol = 0 forces exactly
# j Arnoldi steps). O(k²) applies — diagnostic only, never the production path.
function _ss_probe_residuals!(truer::Vector{RT}, pcr::Vector{RT}, stack::SolverStack,
                              b, x0, k::Int) where {RT}
    o = stack.options
    M! = stack.preconditioner
    r = similar(b)
    z = similar(b)
    for j in 1:k
        rj = gmres(stack.op, b; abstol = 0, reltol = 0, maxiter = j, restart = o.restart,
                   Minv! = M!, x0 = x0)
        apply!(r, stack.op, rj.x)
        @. r = b - r
        push!(truer, RT(norm(r)))
        if M! === nothing
            push!(pcr, RT(norm(r)))
        else
            M!(z, r)
            push!(pcr, RT(norm(z)))
        end
    end
    return nothing
end

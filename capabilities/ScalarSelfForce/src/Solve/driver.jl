# Solve/driver.jl — THE SINGLE-MODE SOLVER-STACK DRIVER (ESSF node `solver-stack-ssf`,
# GLOBAL_DAG row 42, wave 10). This is the JOIN of the reproduction: the first solve of the
# real CRP m-mode problem, on the real CRP domain, with the full GMRES ∘ multigrid ∘ additive
# Schwarz ∘ cached-bottom stack of arXiv:2606.04998 §III.
#
# WHAT THIS FILE OWNS, AND ONLY THIS:
#   1. `solver_options`    — the `SolverStackOptions` twin of the sealed `RunConfig`.
#   2. `mode_fixed_source` — the VOLUME fixed source of the DG system (the SWITCHED evaluator).
#   3. `ModeSolution` / `solve_mode` / `solve_mode!` — one m-mode solve on ONE fixed mesh.
#   4. `rebuild_for_mesh!` — the ONE-STEP rebuild an hp-AMR sweep performs between levels.
#      THE LOOP ITSELF (criterion, iteration count, stopping rule) IS NODE `amr-driver`.
#   5. `mode_residual` / `reconstruct_full_field` / `field_at` — the read-out surface that
#      rows 43, 44, 47 and the figure drivers consume.
#
# EVERYTHING NUMERICAL IS A SEALED PREDECESSOR, called and never re-derived:
#   block_domain (row 39) · worldtube_data / worldtube_massive_source (row 38) ·
#   assemble_operator / operator_builder (row 49) · jump_hooks (row 40) ·
#   boundary_conditions (row 41) · effective_source_m_switched (row 37) ·
#   puncture_psi_chart (row 36) · SolverStack / solve! / rebuild! / ResidualMonitor (row 21) ·
#   refelement / interpolation_matrix (row 8) · helmholtz_shift (row 18).
# NO new physics, NO AMR criterion, NO self-force extraction, NO figure.

export solver_options, ModeSolution, solve_mode, solve_mode!, rebuild_for_mesh!,
       mode_fixed_source, mode_residual, reconstruct_full_field, field_at

# ============================== the options twin ==========================================

"""
    solver_options(config::RunConfig;
                   coarse_levels = config.multigrid.initial_coarse_levels,
                   coarse_orders = nothing, shrink = 1) -> SolverStackOptions

The [`SolverStackOptions`](@ref) twin of the CRP run configuration, FIELD BY FIELD. Every
value is read from the sealed `RunConfig`, which is itself the literal twin of
`anc/ScalarSelfForce.yaml`, so the solver stack is configured by the paper's own input file
and `EllipticCore` keeps its zero dependency on this package:

| `SolverStackOptions` | `RunConfig` field | yaml |
|---|---|---|
| `abstol` | `gmres.absolute_residual` | 183 (`1e-12`) |
| `reltol` | `gmres.relative_residual` | 182 (`0`) |
| `maxiter` | `gmres.max_iterations` | 181 (`500`) |
| `restart` | — (`= maxiter`; spectre's GMRES is full-memory) | — |
| `multigrid_iterations` | `multigrid.iterations` | 187 (`1` V-cycle per GMRES iteration) |
| `initial_coarse_levels` | `multigrid.initial_coarse_levels` | 188 (`L`) |
| `pre_smoothing` | `multigrid.pre_smoothing` | 189 |
| `use_bottom_solver` | `multigrid.use_bottom_solver` | 190 |
| `post_smoothing_at_bottom` | `multigrid.post_smoothing_at_bottom` | 191 |
| `passes` | `schwarz.iterations` | 197 (`2` Schwarz passes per level) |
| `overlap` | — (`2`; the yaml says `MaxOverlap: Auto`) | 198 |
| `subdomain_complex_shift` | `schwarz.subdomain_complex_shift` | 203 (`0`) |

The yaml's multigrid-level `ComplexShift` (194) has no twin: the complex shift is
subdomain-level only in this stack (adopted decision 8) and the production value is `0`.

THE SHIFT POLICY VALUE for this problem is `helmholtz_shift(m, orbit.Omega; slicing = :vtu)`,
which is EXACTLY ZERO under the CRP's piecewise-null slicing (main.tex:680 — the shift "can
be reduced or disabled altogether when hyperboloidal or null slicing is used"). Under
t-slicing it would be `(mΩ)²`, and it would then be applied through the stack's
`precondition_build_op` factory (a rediscretized `β − i·s`), NEVER through the literal
`subdomain_complex_shift` seam, which the sealed row-18 calibration measured to FAIL at the
policy value on this package's massive Gauss formulation ([`helmholtz_shift`](@ref)).
[`solve_mode`](@ref) passes the vtu value through and therefore builds an UNDAMPED
preconditioner, exactly as the paper's production run does.

# THE TWO LEVEL-POLICY KNOBS (exposed on purpose — the AMR driver must choose)

`coarse_levels` is the yaml's `InitialCoarseLevels`: the number of COARSER multigrid levels.
The stack's `maxlevels` is `coarse_levels + 1`. At the production point `L = 0` that is ONE
level, so the multigrid V-cycle IS the direct bottom solve of the whole problem — which is
precisely what the paper reports for its first AMR level (`fig:linsolv`, main.tex:701: "The
first AMR level converges immediately because the multigrid bottom solver solves the problem
directly"). MEASURED at `a = 0.5`, `r₀ = 10`, `m = 2`, `L = 0`, `P = 4`: **1 GMRES iteration**
at a true residual of `1.2e-15`.

AN AMR SWEEP MUST RAISE IT. main.tex:674 keeps the bottom grid fixed while "AMR adds finer
levels"; a driver that left `coarse_levels` at `L` would instead let the bottom grid follow
the fine grid upward. Pass `coarse_levels = nothing` to h-coarsen all the way to the
irreducible block structure (main.tex:672) — the fixed bottom grid — or raise it by one per
h-step. MEASURED on the production lattice at `L = 0`: `nothing` gives 2 levels
(320 → 160 dofs) and **10** GMRES iterations; at `L = 1`, 3 levels and 15; at `L = 2`,
4 levels and 25. Under the CRP's own hp-AMR sweep (h at the puncture, p elsewhere) the counts
stay flat instead: 10, 9, 10, 12 over four steps. See this node's evidence table.

`coarse_orders` is the coarse-level POLYNOMIAL policy of [`MultigridHierarchy`](@ref):
`nothing` (the default) p-coarsens by the MIN rule; an integer pins every coarser level at
that order. spectre coarsens `h` ONLY and allocates every level from the same initial
extents, so `coarse_orders = P` is the spectre-faithful setting and is what keeps the cached
bottom LU alive across an `hp_refine` sweep (main.tex:674). It is a TRADE, not a free win —
the pinned coarse space falls behind a p-refined fine grid and costs iterations. The default
here stays `nothing` because at `L = 0`, before any hp step, the two rules COINCIDE; row 43
chooses per sweep with the measured table in front of it.

`shrink` stays at the sealed row-21 default `1`, the only value that never costs an iteration
on the anisotropic meshes the CRP's `Isotropy: Anisotropic` policy (yaml:145) produces.

Tolerances keep the config's own real type, so a `RunConfig{BigFloat}` yields `BigFloat`
tolerances and the solve is never silently Float64'd.
"""
function solver_options(config::RunConfig;
                        coarse_levels::Union{Nothing,Integer} =
                            config.multigrid.initial_coarse_levels,
                        coarse_orders::Union{Nothing,Integer} = nothing,
                        shrink::Real = 1)
    return SolverStackOptions(; abstol = config.gmres.absolute_residual,
                              reltol = config.gmres.relative_residual,
                              maxiter = config.gmres.max_iterations,
                              multigrid_iterations = config.multigrid.iterations,
                              initial_coarse_levels = coarse_levels,
                              pre_smoothing = config.multigrid.pre_smoothing,
                              use_bottom_solver = config.multigrid.use_bottom_solver,
                              post_smoothing_at_bottom =
                                  config.multigrid.post_smoothing_at_bottom,
                              shrink = shrink,
                              coarse_orders = coarse_orders,
                              overlap = 2,
                              passes = config.schwarz.iterations,
                              subdomain_complex_shift =
                                  config.schwarz.subdomain_complex_shift)
end

# ============================ the volume fixed source =====================================

"""
    mode_fixed_source(T, dom::BlockDomain, orbit::CircularOrbit, m;
                      mesh = dom.mesh, phi_p = 0) -> Vector{T}

The VOLUME fixed source of the CRP linear system, in the massive `M·S` form the DG residual
consumes: `+S^eff_m` at every collocation node of every element INSIDE the worldtube and
EXACTLY zero everywhere else, assembled through the sealed generic seam
[`worldtube_massive_source`](@ref)`(mesh, isinside, source, T)`.

# THE SWITCHED-SOURCE / DIRECT-MORTAR RULE (normative — the one statement for this package)

THE VOLUME SOURCE IS SWITCHED. The evaluator here is [`effective_source_m_switched`](@ref),
never the bare direct route: the direct assembly of `S^eff_m` cancels four orders of
singularity, so its roundoff grows like `λ⁻²` and it is unusable below `λ ≈ 1e-2` at
`Float64` (main.tex:488-492) — and the AMR driver refines AT the puncture eight times,
driving the innermost quadrature nodes exactly there. The switch picks whichever of the
near-worldline series and the direct route has the smaller predicted error, so the volume
source stays accurate all the way in.

THE WORLDTUBE MORTAR DATA ARE A DIFFERENT OBJECT AND STAY DIRECT. `Ψ^𝒫` and its two
derivatives on the tube boundary are TRACE data, evaluated far from the worldline where the
direct route is accurate and the truncated series is not; [`worldtube_data`](@ref) builds
them from the direct [`effective_source_bundle`](@ref), and this function never touches them.
Consequently `mode_fixed_source(...)` and `worldtube_massive_source(wt)` are DIFFERENT
vectors wherever the switch fires, and only the former may be used as the right-hand side.

MEASURED ROUTE CENSUS at the production point (`a = 0.5`, `r₀ = 10`, `m = 2`, `L = 0`,
`P = 4`, `Float64`): 0 of the 64 in-tube quadrature nodes select the series — the closest
node sits at `λ = 1.004`, three orders outside the Float64 crossover `λ* ≈ 1.2e-2`, and it is
still 0 of 4096 after three `hp_refine` steps (`λ_min = 0.333`). On THIS lattice the switch is
therefore inert; it is the AMR puncture-snap resolution, not the initial mesh, that decides,
and the rule is written here so the source cannot silently become the roundoff-limited one
when it does.

The physical right-hand side of the solve is this vector: [`solve!`](@ref) folds the
operator's own fixed-data lift (`− operator_on_zero(op)`, which carries the worldtube jump
data) internally, so it must NOT be pre-folded here.

`phi_p` is the particle's Boyer–Lindquist azimuth; the field-point horizon-regular rotation
is internal to the sealed evaluator and must never be applied twice.
"""
function mode_fixed_source(::Type{T}, dom::BlockDomain, orbit::CircularOrbit, m::Integer;
                           mesh::Mesh2D = dom.mesh, phi_p::Real = 0) where {T<:Number}
    p = Puncture(orbit)
    mm = Int(m)
    z2 = dom.chart === :z2
    src(x1, x2) = effective_source_m_switched(p, mm, x1, z2 ? sqrt(x2) : x2; phi_p = phi_p)
    return worldtube_massive_source(mesh, el -> inside_worldtube(dom, el), src, T)
end

# ================================ the solved mode =========================================

"""
    ModeSolution{T,RT}

One solved m-mode problem of arXiv:2606.04998: the block domain, the mesh-bound worldtube
bundle, the assembled operator, the solver stack, the physical right-hand side and the
solution, together with the [`ResidualMonitor`](@ref) whose records are the `fig:linsolv`
data surface. Build with [`solve_mode`](@ref); re-solve with [`solve_mode!`](@ref); move to a
new mesh with [`rebuild_for_mesh!`](@ref).

Fields:
  - `config`, `orbit`, `m`, `phi_p` — the problem.
  - `dom::BlockDomain{RT}`, `mesh::Mesh2D{RT}` — the CRP block lattice and the CURRENT mesh
    (an AMR level need not be `dom.mesh`).
  - `wt::WorldtubeData` — the bundle for `mesh`. MESH-BOUND: rebuilt on every mesh change.
  - `op::MeshOperator` — the operator on `mesh`. `stack::SolverStack` — the
    GMRES ∘ multigrid ∘ Schwarz ∘ bottom stack over `op` and the per-level factory.
  - `b::Vector{T}` — the PHYSICAL right-hand side ([`mode_fixed_source`](@ref)); never
    pre-folded, because [`solve!`](@ref) folds `operator_on_zero` itself.
  - `u::Vector{T}` — THE SOLUTION, and note what it is: `Ψ^ℛ_m` INSIDE the worldtube and
    `Ψ_m` outside it (eq:jump_worldtube, main.tex:513-523). Use
    [`reconstruct_full_field`](@ref) for `Ψ_m` everywhere.
  - `monitor::ResidualMonitor{RT}` — one [`SolveRecord`](@ref) per solve, in call order;
    across an AMR sweep that is one record per level, i.e. `fig:linsolv` itself.
  - `record`, `result` — the last solve's [`SolveRecord`](@ref) and sealed
    [`GMRESResult`](@ref); both `nothing` before the first solve and after a rebuild.

# ONE STACK, ONE MODE, ONE THREAD (normative)

Each m-mode is an INDEPENDENT problem (main.tex:716) and gets its OWN `ModeSolution` with its
OWN stack and its own two caches. Nothing is shared between modes: the Schwarz and bottom
caches are keyed by `(element, subdomain shape, epoch, shift)` and would be correct but
useless across modes, and neither cache nor the stack is thread-safe. The process-parallel
fan-out (row 47) therefore builds one `ModeSolution` PER PROCESS, never one per thread.

# ONE EPOCH, BUMPED BY THE CALLER (normative)

The stack's two caches share ONE `epoch` and it is AUTHORITATIVE (the sealed row-21
two-caches-one-epoch contract). A pure MESH change needs no bump — the mesh is in both keys —
which is why [`rebuild_for_mesh!`](@ref) does not bump it, and why an unchanged coarser level
costs zero factorizations while the bottom LU survives as long as the coarsest grid and
operator do. Any matrix-affecting change the mesh does NOT carry (a different `m`, orbit,
penalty, hook or boundary-condition policy) is a different problem and wants a NEW
`ModeSolution`, not an epoch bump on this one: a different mode is a different `RunConfig`.
"""
mutable struct ModeSolution{T<:Number,RT<:AbstractFloat}
    config::RunConfig{RT}
    orbit::CircularOrbit
    m::Int
    phi_p::RT
    dom::BlockDomain{RT}
    mesh::Mesh2D{RT}
    wt::WorldtubeData
    op::MeshOperator
    stack::SolverStack
    b::Vector{T}
    u::Vector{T}
    monitor::ResidualMonitor{RT}
    record::Union{Nothing,SolveRecord{RT}}
    result::Any
end

"""
    solve_mode(orbit::CircularOrbit, config::RunConfig, m = config.orbit.m_mode_number;
               T = ComplexF64, options = solver_options(config),
               dom = block_domain(orbit, config), mesh = dom.mesh, phi_p = 0,
               x0 = nothing, monitor = ResidualMonitor(real(T))) -> ModeSolution

Solve ONE m-mode of the CRP elliptic problem on a FIXED mesh with the full solver stack — the
join of arXiv:2606.04998 §II-III. This is the single home of the sequence

    dom   = block_domain(orbit, config)                                     # row 39
    wtf   = mesh -> worldtube_data(dom, orbit, config, m; mesh = mesh)      # row 38
    op    = assemble_operator(T, dom, wtf(mesh), orbit, config, m; mesh)    # row 49
    stack = SolverStack(op, operator_builder(T, dom, orbit, config, m; worldtube = wtf);
                        options = solver_options(config))                   # row 21
    b     = mode_fixed_source(T, dom, orbit, m; mesh)                       # the SWITCHED source
    solve!(u, stack, b)                                                     # folds internally

and nothing more. The argument order is this package's own sealed convention
`(…, orbit, config, m)` — the one `worldtube_data`, `jump_hooks`, `assemble_operator` and
`vtu_junction_faces` already use.

`mesh` may be any mesh of `dom`'s block lattice; `dom` is accepted so a caller that already
built the lattice (the AMR driver, the parallel fan-out, a figure driver) does not rebuild
it. `x0` warm-starts GMRES. `m` must equal the config's `MModeNumber` — the sealed guards say
so, and a mismatch would silently solve a different equation.

GATE ON THE RESIDUAL, NOT ON `converged`: `sol.result.residual` is the sealed
[`GMRESResult`](@ref)'s RECOMPUTED true residual, while at `m = 0` the sealed GMRES can
report `:stagnation` at machine residual (the row-41 caveat). [`mode_residual`](@ref) is the
independent recomputation against the full operator table.

THE PRECONDITIONER IS UNDAMPED, deliberately: `helmholtz_shift(m, orbit.Omega;
slicing = :vtu)` is exactly zero for the CRP's null slicing (main.tex:680), so there is no
damped factory to hand over and `precondition_build_op === build_op`. See
[`solver_options`](@ref) for what a t-sliced caller would do instead.

MEASURED at the row-42 oracle point (`a = 0.5`, `r₀ = 10`, `m = 2`, `L = 0`, `P = 4`,
`z²` chart, 320 dofs): converged in **1** preconditioned GMRES iteration at a true residual
of `1.2e-15` — the paper's own first-AMR-level behaviour (main.tex:701), because the yaml's
`InitialCoarseLevels: L = 0` makes the V-cycle a direct bottom solve. At `a = 0`:
`6.7e-16`. With a genuine two-level hierarchy on the same mesh (`coarse_levels = nothing`):
10 iterations at `2.5e-13`.
"""
function solve_mode(orbit::CircularOrbit, config::RunConfig{RT},
                    m::Integer = config.orbit.m_mode_number;
                    T::Type = ComplexF64,
                    options::SolverStackOptions = solver_options(config),
                    dom::BlockDomain = block_domain(orbit, config),
                    mesh::Mesh2D = dom.mesh,
                    phi_p::Real = 0,
                    x0::Union{Nothing,AbstractVector} = nothing,
                    monitor::ResidualMonitor = ResidualMonitor(real(T))) where {RT}
    mm = Int(m)
    wtf = mh -> worldtube_data(dom, orbit, config, mm; mesh = mh, phi_p = phi_p)
    wt = wtf(mesh)
    op = assemble_operator(T, dom, wt, orbit, config, mm; mesh = mesh)
    build_op = operator_builder(T, dom, orbit, config, mm; worldtube = wtf)
    stack = SolverStack(op, build_op; options = options)
    b = mode_fixed_source(T, dom, orbit, mm; mesh = mesh, phi_p = phi_p)
    sol = ModeSolution{T,RT}(config, orbit, mm, RT(phi_p), dom, mesh, wt, op, stack, b,
                             zeros(T, mesh.ndofs), monitor, nothing, nothing)
    return solve_mode!(sol; x0 = x0)
end

"""
    solve_mode!(sol::ModeSolution; x0 = nothing, true_residual_probe = false) -> sol

Run (or re-run) the solve of `sol` in place: `solve!(sol.u, sol.stack, sol.b)` with the
PHYSICAL right-hand side — the sealed [`solve!`](@ref) applies the `operator_on_zero` fold
itself, once, so `sol.b` is never pre-folded and never folded twice (folding twice is an
`O(1)` wrong answer and the caller's error).

Appends one [`SolveRecord`](@ref) to `sol.monitor` and stores it in `sol.record` alongside
the [`GMRESResult`](@ref) in `sol.result`. The FIRST solve after a build or a
[`rebuild_for_mesh!`](@ref) carries that rebuild's `build_time` and cache-counter deltas;
later solves on the same preconditioner record `0.0`, so summing the column over an AMR
sweep gives its true build cost with no double counting.

`x0` warm-starts the iteration — the AMR driver passes the transferred previous level's
solution, and `x0 = sol.u` is safe because the sealed `gmres` copies its initial guess.
`true_residual_probe` additionally fills the record's EXACT per-iteration residuals by
re-solving with `maxiter = 1, 2, …`; that is `O(k²)` operator applies and is a diagnostic,
never the production path.
"""
function solve_mode!(sol::ModeSolution; x0::Union{Nothing,AbstractVector} = nothing,
                     true_residual_probe::Bool = false)
    res = solve!(sol.u, sol.stack, sol.b; x0 = x0, monitor = sol.monitor,
                 true_residual_probe = true_residual_probe)
    sol.result = res
    sol.record = sol.monitor.records[end]
    return sol
end

"""
    rebuild_for_mesh!(sol::ModeSolution, mesh::Mesh2D; u = nothing, options = nothing,
                      epoch = sol.stack.epoch) -> sol

Move `sol` to a NEW mesh of the same block lattice in one step — the per-level rebuild an
hp-AMR sweep performs between refinements. THE LOOP ITSELF IS NODE `amr-driver` (row 43):
this function is the step it calls, and it owns no criterion, no iteration count and no
stopping rule.

Everything mesh-bound is rebuilt, in the order the sealed contracts require:

1. `wt = worldtube_data(dom, orbit, config, m; mesh)` — a [`WorldtubeData`](@ref) is
   MESH-BOUND (one `FaceState` per mortar, a nodal source in that mesh's layout), and a stale
   one is REFUSED by the sealed `jump_hooks` guard rather than mis-applied.
2. `op = assemble_operator(…; mesh)` — hooks and boundary conditions recomputed from the new
   mesh's own faces (both tables are keyed by element id, which AMR changes).
3. `rebuild!(stack, op; build_op = …)` — THE FACTORY IS HANDED OVER, because it is what the
   coarse levels are rediscretized from; leaving it stale silently preconditions with the old
   operator on every level below the finest, which the residual never reveals.
4. `b = mode_fixed_source(…; mesh)` — the switched volume source on the new nodes.

NO EPOCH BUMP: a pure mesh change is already in both cache keys (the sealed
two-caches-one-epoch contract), so an unchanged coarser level costs zero factorizations and
the bottom LU survives as long as the coarsest grid and operator do. Pass `epoch` explicitly
only if something the mesh does not carry also changed.

`u` is the solution to carry over; the row-43 flow is

    flags = criterion(sol.mesh, sol.u)
    mesh  = hp_refine(sol.mesh, flags)
    u     = transfer_solution(sol.mesh, mesh, sol.u)
    rebuild_for_mesh!(sol, mesh; u = u)
    solve_mode!(sol; x0 = sol.u)

Omitted, the solution restarts from zero. `options` replaces the stack's options — the seam
for raising `coarse_levels` as AMR adds finer levels, and for choosing `coarse_orders` (see
[`solver_options`](@ref)); omitted, the current options are kept.

`sol.record` and `sol.result` are CLEARED: they describe a solve on the old mesh. The monitor
is KEPT, so its records accumulate one per AMR level — which is exactly `fig:linsolv`.
"""
function rebuild_for_mesh!(sol::ModeSolution{T}, mesh::Mesh2D;
                           u::Union{Nothing,AbstractVector} = nothing,
                           options::Union{Nothing,SolverStackOptions} = nothing,
                           epoch::Integer = sol.stack.epoch) where {T}
    u === nothing || length(u) == mesh.ndofs || throw(ArgumentError(
        "rebuild_for_mesh!: the carried solution has length $(length(u)) but the new mesh " *
        "has $(mesh.ndofs) dofs — transfer it first, transfer_solution(old, new, u)"))
    wtf = mh -> worldtube_data(sol.dom, sol.orbit, sol.config, sol.m; mesh = mh,
                               phi_p = sol.phi_p)
    wt = wtf(mesh)
    op = assemble_operator(T, sol.dom, wt, sol.orbit, sol.config, sol.m; mesh = mesh)
    build_op = operator_builder(T, sol.dom, sol.orbit, sol.config, sol.m; worldtube = wtf)
    options === nothing || (sol.stack.options = options)
    rebuild!(sol.stack, op; epoch = epoch, build_op = build_op)
    sol.mesh = mesh
    sol.wt = wt
    sol.op = op
    sol.b = mode_fixed_source(T, sol.dom, sol.orbit, sol.m; mesh = mesh, phi_p = sol.phi_p)
    sol.u = u === nothing ? zeros(T, mesh.ndofs) : convert(Vector{T}, u)
    sol.record = nothing
    sol.result = nothing
    return sol
end

# ================================== read-out ==============================================

"""
    mode_residual(sol::ModeSolution) -> Real

The DG residual of the stored solution against the FULL operator table, `‖b − A_full(u)‖₂`
with `b` the PHYSICAL right-hand side — the independent check that the solve solved THIS
problem, whatever GMRES reported.

It equals the sealed [`GMRESResult`](@ref)`.residual` to roundoff, because the operator is
affine: `A_full(u) = A_lin(u) + operator_on_zero(op)` exactly, and that is precisely the fold
[`solve!`](@ref) applies. A discrepancy therefore means the right-hand side and the operator
belong to DIFFERENT discretizations — a stale `WorldtubeData`, a mesh mismatch, a
double-folded right-hand side — which is exactly the failure this function exists to catch.
MEASURED at the oracle point: `1.08e-15` against `result.residual = 1.21e-15`.

Gate on this and on `sol.result.residual`; never on `converged` alone.
"""
function mode_residual(sol::ModeSolution)
    y = similar(sol.u)
    apply_full!(y, sol.op, sol.u)
    return norm(sol.b .- y)
end

"""
    reconstruct_full_field(sol::ModeSolution) -> Vector

`Ψ_m` on the WHOLE domain as one nodal vector: the stored solution outside the worldtube
(where it already IS the retarded field) and `Ψ^ℛ_m + Ψ^𝒫_m` inside it, with `Ψ^𝒫_m` taken
from the sealed [`puncture_psi_chart`](@ref) at each in-tube collocation node
(eq:jump_worldtube, `Ψ^ℛ_m = Ψ_m − Ψ^𝒫_m`, main.tex:513-523).

Outside the tube the result is `sol.u` BITWISE — nothing is recomputed there.

THE RECONSTRUCTION IS WHAT MAKES THE FIELD CONTINUOUS across the worldtube boundary, and that
is the sharpest available check that the whole jump machinery is wired correctly. MEASURED at
the oracle point, both one-sided traces taken at the operator's OWN mortar points: the worst
jump of the reconstructed field is `2.80e-4` against a field scale of `0.201`, while the same
measurement on the RAW solution is `0.2017` — the full size of `Ψ^𝒫`, a factor **722** larger
— and with the worldtube hooks removed from the operator the reconstructed jump is `0.2019`
even though GMRES still converges happily. The jump is a DISCRETIZATION envelope and
converges spectrally in `p`: `2.80e-4 → 1.34e-5 → 1.09e-6` at `P = 4, 6, 8`; the one-sided
flux jump behaves the same way (`2.18e-3 → 1.82e-4 → 1.67e-5`, against `2.51e-2` unreconstructed).

This is the field `fig:modes` plots outside the regularized region (it plots `Ψ^ℛ_m` inside
it, i.e. `sol.u`); the self-force extraction of row 44 also uses the REGULAR field `sol.u` at
the particle, not this one.
"""
function reconstruct_full_field(sol::ModeSolution{T}) where {T}
    mesh, dom = sol.mesh, sol.dom
    p = Puncture(sol.orbit)
    z2 = dom.chart === :z2
    out = copy(sol.u)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        inside_worldtube(dom, el) || continue
        em = element_matrices(T, mesh.blocks[mesh.blockindex[el.block]], el, mesh.orders[e])
        off = mesh.offsets[e]
        nx, ny = mesh.orders[e]
        for j in 1:ny, i in 1:nx
            x1, x2 = em.x[1][i], em.x[2][j]
            out[off + tensor_lin(i, j, nx)] +=
                puncture_psi_chart(p, sol.m, x1, z2 ? sqrt(x2) : x2; phi_p = sol.phi_p).psi
        end
    end
    return out
end

"""
    field_at(sol::ModeSolution, x1, x2; element = nothing) -> T
    field_at(sol::ModeSolution, values::AbstractVector, x1, x2; element = nothing) -> T

Evaluate a nodal field of `sol`'s mesh at the DISCRETIZATION-CHART point `(x1, x2)` by the
containing element's own tensor-product Lagrange interpolation — the sealed
[`interpolation_matrix`](@ref) of [`refelement`](@ref) in each dimension, never a re-derived
basis. The short form evaluates `sol.u`; pass [`reconstruct_full_field`](@ref)`(sol)` to
evaluate `Ψ_m` instead, or any other nodal vector of the same mesh.

COORDINATES ARE THE DISCRETIZATION CHART, not `(r, θ)`: `x1` is the horizon-penetrating
Boyer–Lindquist radius `r`, and `x2` is the domain's own angular coordinate — `w = z²` under
the production `ImposeEquatorialSymmetry` (`sol.dom.chart === :z2`, main.tex:389-394) and
`z = cos θ` otherwise. Read `sol.dom.chart` rather than assuming.

`element` selects WHICH element's polynomial is used, as an index into `sol.mesh.elements`.
A DG field is discontinuous across element faces, so a point ON a face has two legitimate
values and this keyword is how a caller asks for a specific one-sided trace (the worldtube
seam, an extraction radius). Omitted, the first element containing the point is used, with a
`16 eps` logical tolerance so a face or corner point resolves rather than falling through.

Interpolation is exact for the element's own polynomial degree and reproduces a collocation
node's value exactly (MEASURED worst deviation `8.4e-17` against `‖u‖_∞ = 0.18` over all 320
nodes); evaluation at a block-boundary coordinate such as `r = r₊` or `r = R` is the face
trace and is finite. Extrapolation outside the chosen element is possible and is the caller's
responsibility. A point in no element is an `ArgumentError`.
"""
field_at(sol::ModeSolution, x1::Real, x2::Real; element = nothing) =
    field_at(sol, sol.u, x1, x2; element = element)

function field_at(sol::ModeSolution, values::AbstractVector, x1::Real, x2::Real;
                  element = nothing)
    mesh = sol.mesh
    length(values) == mesh.ndofs || throw(ArgumentError(
        "field_at: the field has length $(length(values)) but the mesh has " *
        "$(mesh.ndofs) dofs"))
    e = element === nothing ? _dr_locate(mesh, x1, x2) : Int(element)
    (1 <= e <= length(mesh.elements)) || throw(ArgumentError(
        "field_at: element index $e is outside 1:$(length(mesh.elements))"))
    return _dr_interp(mesh, values, e, _dr_logical(mesh, e, x1, x2))
end

# ---- element location + tensor-product interpolation (sealed refelement machinery) --------

# logical coordinates of a physical point within element `e` (outside the element ⇒ |ξ| > 1)
function _dr_logical(mesh::Mesh2D{RT}, e::Int, x1::Real, x2::Real) where {RT}
    el = mesh.elements[e]
    blk = mesh.blocks[mesh.blockindex[el.block]]
    return ntuple(2) do d
        lo, hi = logical_interval(RT, el, d)
        ζ = map_inverse(blockmap(blk, d), RT(d == 1 ? x1 : x2))
        (2ζ - (lo + hi)) / (hi - lo)
    end
end

# the first element containing the point; a face point belongs to both its neighbours, and
# `field_at`'s `element` keyword is how a caller picks one (see that docstring)
function _dr_locate(mesh::Mesh2D{RT}, x1::Real, x2::Real) where {RT}
    tol = 16 * eps(RT)
    for e in eachindex(mesh.elements)
        ξ = _dr_logical(mesh, e, x1, x2)
        (-1 - tol <= ξ[1] <= 1 + tol && -1 - tol <= ξ[2] <= 1 + tol) && return e
    end
    throw(ArgumentError(
        "field_at: ($x1, $x2) lies in no element of this mesh. The coordinates are the " *
        "DISCRETIZATION chart: x1 = r (horizon-penetrating Boyer–Lindquist) and x2 = z² " *
        "under equatorial symmetry, z = cosθ otherwise — check `sol.dom.chart`"))
end

function _dr_interp(mesh::Mesh2D{RT}, values::AbstractVector, e::Int, ξ) where {RT}
    nx, ny = mesh.orders[e]
    Px = interpolation_matrix(refelement(RT, nx), [ξ[1]])
    Py = interpolation_matrix(refelement(RT, ny), [ξ[2]])
    off = mesh.offsets[e]
    acc = zero(eltype(values))
    @inbounds for j in 1:ny, i in 1:nx
        acc += Px[1, i] * Py[1, j] * values[off + tensor_lin(i, j, nx)]
    end
    return acc
end

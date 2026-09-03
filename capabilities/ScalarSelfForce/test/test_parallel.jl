# test_parallel.jl — ESSF node `mmode-parallel-driver` (GLOBAL_DAG row 47, wave 12, worker
# ESSF-W12-B), written BEFORE `src/Solve/parallel.jl` exists (TDD; `tests_red.log` captured
# against the missing file). THIS IS THE PRODUCTION FAN-OUT: the paper's m = 0…20 sweep run as
# one INDEPENDENT PROCESS per mode, under the mission's binding compute policy.
#
# WHAT THE NODE DELIVERS (`Solve/parallel.jl`):
#   * `solve_modes_parallel(orbit, config, ms)` → `ModeSet` — the process fan-out: a pool of
#     `min(budget ÷ threads_per_worker, length(ms))` worker processes, each with pinned BLAS
#     and its own stack/caches, one `solve_mode_amr` per mode, results collected in `ms` order.
#   * `solve_modes_serial(orbit, config, ms)` → `ModeSet` — the SAME loop in-process. It is the
#     BITWISE ORACLE of the parallel path and the only thing gate (b) can compare against.
#   * `ModeSet` / `ModeResult` / `ModeFailure` / `ModeRunRecord` / `LoadSample` — the F-ready
#     per-mode payload (the `AmrLevel` record of every level + the final `sol.u`) plus the run
#     record rows 50/51/53 need to size their own pools.
#   * `ModePool` / `mode_pool` — a guarded, REUSABLE pool of worker processes, because a fresh
#     process pays ~200 s of JIT before it can solve anything (MEASURED below).
#   * `system_load`, `mode_pool_size`, `mode_config` — the three policy pieces: the live
#     whole-system reading, the pool arithmetic, and the per-mode `RunConfig`.
#
# ORACLE (the DAG row-47 gate column + the READ-ONLY mirrors + the binding compute policy):
#   * `main.tex:716-720` — "Each m-mode is independent … each mode solves in a few seconds",
#     and the m = 0…20 sweep; `main.tex:684-695` (fig:convergence, fig:modes) — modes up to
#     m = 20. `CRP_PAPER_MMODES = 0:20` is the sealed twin of that sweep.
#   * the `compute-usage` skill (BINDING owner policy): at most 24 CPUs for our jobs; whole
#     system below 90 % before AND after; at most 4 threads per Julia process; for an
#     embarrassingly parallel sweep "size the worker pool as min(24, modes) and leave it
#     there; never auto-scale up"; GPU banned.
#   * DESIGN.md §5.1 — one worker PROCESS per m, single-threaded BLAS and Julia per worker,
#     `nworkers × threads ≤ 24`, poll the total system load to stay < 90 %.
#   * the sealed row-43 API contract (`evidence/amr-driver/api_contract_rows44_47_48_52.md`)
#     — one `AmrHistory` per PROCESS, never per thread; `OPENBLAS_NUM_THREADS=1` is worth a
#     factor 8.9 and is what makes a run bitwise reproducible (the diagnostics move by
#     4.2e-12 relative between BLAS thread counts); the per-mode cost is m-independent to
#     ~10 %, so `min(24, nmodes)` is the right pool size.
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here): `amr-driver` (`solve_mode_amr`,
# `AmrHistory`, `AmrLevel`, `amr_flags`, `amr_level_options`), `solver-stack-ssf`
# (`ModeSolution`), `runtime-config-fixture` (`crp_paper_config`, `CRP_PAPER_MMODES`,
# `RunConfig`), `orbits-adapters` (`CircularOrbit`).
#
# GATE CONSTANTS — every one calibrated on the throwaway prototypes
# `evidence/mmode-parallel-driver/proto1_fanout.jl` and `proto2_threads.jl` BEFORE this file
# was written, never tuned afterwards to pass. The measured value is quoted beside each gate.
#
# WALL-TIME POLICY OF THIS FILE: the test itself must obey the compute policy it gates. Every
# pool it starts is sized from the LIVE load through the same `mode_pool_size` + guard the
# production path uses, the AMR sweeps are 0-2 levels (not the paper's 8), and the full
# 21 × 8 run is the EVIDENCE ARTEFACT, not a unit test.

using Test
using LinearAlgebra
using ScalarSelfForce
using EllipticCore
using Distributed        # the leak check only: `nprocs()` must return to its start value

# ---- configs (pa_-prefixed: sibling test files share the runtests.jl include scope) --------
pa_cfg(; spin = 0.5, r0 = 10.0, m = 2, L = 0, P = 4, wr = 4.8) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L,
                     worldtube_radius = wr)
pa_orb(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)

# recursive field-by-field comparison — Julia's `==` on composite types falls back to `===`,
# so the `mode_config` gate needs its own structural equality.
function pa_same(a, b)
    a === b && return true
    typeof(a) === typeof(b) || return false
    if a isa Number || a isa AbstractString || a isa Symbol || a isa Nothing
        return a == b
    elseif a isa Tuple || a isa AbstractArray
        length(a) == length(b) || return false
        return all(pa_same(x, y) for (x, y) in zip(a, b))
    else
        fs = fieldnames(typeof(a))
        isempty(fs) && return a == b
        return all(pa_same(getfield(a, f), getfield(b, f)) for f in fs)
    end
end

# every bitwise-comparable slot of one mode's payload
pa_bitwise(x::ModeResult, y::ModeResult) =
    x.m == y.m && length(x.u) == length(y.u) && all(x.u .=== y.u) &&
    length(x.levels) == length(y.levels) &&
    all(lx.psi === ly.psi && lx.dpsi_r === ly.dpsi_r && lx.dpsi_w === ly.dpsi_w &&
        lx.element == ly.element && lx.mesh.ndofs == ly.mesh.ndofs &&
        lx.record.iterations == ly.record.iterations &&
        lx.record.final_residual === ly.record.final_residual
        for (lx, ly) in zip(x.levels, y.levels))

pa_rel(x, y) = abs(x - y) / max(abs(y), floatmin(Float64))

# how many processes this session had BEFORE the node touched anything: every pool this file
# starts must be given back (a leaked worker is 2.5 GB of RSS and a core the owner did not lend).
pa_nprocs0 = Distributed.nprocs()

@testset "mmode-parallel-driver (row 47)" begin

# ============================ (0) the policy arithmetic ===================================
# Pure, launch-free gates on the pool sizing: the owner policy is `min(24, modes)` and
# `nworkers × threads_per_worker ≤ 24`, ALWAYS, with every clamp reported (never silent).
@testset "(0) pool sizing obeys the compute budget" begin
    @test PARALLEL_CPU_BUDGET == 24
    @test PARALLEL_MAX_THREADS_PER_WORKER == 4
    @test PARALLEL_LOAD_GUARD == 0.90

    # the paper's own sweep: 21 modes, one single-threaded process each
    p = mode_pool_size(21)
    @test p.nworkers == 21
    @test p.threads_per_worker == 1
    @test !p.clamped
    @test p.nworkers * p.threads_per_worker <= PARALLEL_CPU_BUDGET

    # more modes than the budget: the pool stops at the budget, it never auto-scales up
    @test mode_pool_size(100).nworkers == 24
    @test mode_pool_size(100).clamped == false   # the DEFAULT is the budget, not a clamp
    # fewer modes than the budget: one worker per mode, no idle processes
    @test mode_pool_size(3).nworkers == 3

    # threads per worker divide the budget
    @test mode_pool_size(100; threads_per_worker = 4).nworkers == 6
    @test mode_pool_size(100; threads_per_worker = 3).nworkers == 8
    for nth in 1:4, nm in (1, 5, 21, 100)
        q = mode_pool_size(nm; threads_per_worker = nth)
        @test q.nworkers * q.threads_per_worker <= PARALLEL_CPU_BUDGET
        @test q.nworkers >= 1
    end

    # an EXPLICIT over-budget request is clamped and REPORTED (gate c)
    q = mode_pool_size(100; nworkers = 64)
    @test q.nworkers == 24 && q.requested == 64 && q.clamped
    q = mode_pool_size(21; nworkers = 40, threads_per_worker = 2)
    @test q.nworkers == 12 && q.clamped && q.nworkers * q.threads_per_worker <= 24
    # a thread count beyond the owner's 4-per-process limit is clamped and reported
    q = mode_pool_size(21; threads_per_worker = 16)
    @test q.threads_per_worker == 4 && q.clamped
    # a budget beyond the owner's 24 is clamped and reported
    q = mode_pool_size(100; budget = 128)
    @test q.budget == 24 && q.nworkers == 24 && q.clamped
    # asking for more workers than modes wastes processes: clamped down
    q = mode_pool_size(3; nworkers = 20)
    @test q.nworkers == 3 && q.clamped
end

# ============================== (1) the live load reading =================================
@testset "(1) system_load reads the LIVE whole-system usage" begin
    s = system_load()
    @test s isa LoadSample
    @test s.ncpu == Sys.CPU_THREADS
    @test s.cpu_fraction >= 0.0
    @test isfinite(s.cpu_fraction)
    @test s.loadavg1 >= 0.0
    @test s.source in (:ps, :loadavg)
    @test s.time >= 0.0
    # a second reading, relative to a caller-supplied origin, is later
    t0 = time()
    s2 = system_load(t0)
    @test s2.time >= 0.0
    # this host is 256-core and shared; the reading must be a FRACTION, not a percentage
    @test s.cpu_fraction < 100.0
end

# ============================ (2) the per-mode config =====================================
# The fan-out sends ONE template config; every worker rebuilds the RunConfig for its own m.
# Only two things may move: `orbit.m_mode_number` and the m = 0 no-outer-BC yaml branch.
@testset "(2) mode_config is the yaml's own per-m config" begin
    cfg = pa_cfg(; m = 2)
    for m in (0, 1, 2, 5, 20)
        c = mode_config(cfg, m)
        @test c.orbit.m_mode_number == m
        @test pa_same(c, pa_cfg(; m = m))          # identical to building it from scratch
    end
    # the m = 0 branch (yaml:114-125): no outer radial BC
    @test mode_config(cfg, 0).domain.outer_radial_bc === nothing
    @test mode_config(cfg, 3).domain.outer_radial_bc isa SommerfeldBC
    @test mode_config(cfg, 3).domain.outer_radial_bc.m_mode_number == 3
    # a config already at m keeps every field
    @test pa_same(mode_config(cfg, 2), cfg)
    # NON-paper configs are preserved too (only the two m-dependent slots move)
    other = pa_cfg(; spin = -0.9, r0 = 6.0, m = 4, P = 6, L = 1, wr = 1.2)
    c = mode_config(other, 7)
    @test c.orbit.m_mode_number == 7
    @test pa_same(c.domain.block_bounds_r, other.domain.block_bounds_r)
    @test c.domain.initial_grid_points == (6, 6)
    @test c.multigrid.initial_coarse_levels == other.multigrid.initial_coarse_levels
    @test_throws ArgumentError mode_config(cfg, -1)
end

# ========================== (3) the pre-launch load guard ==================================
# The guard reads the LIVE system usage and REFUSES to launch when the planned share would
# push the whole system past `load_guard`. It must name the numbers, and it must never launch
# a process on the way to the refusal.
@testset "(3) the guard refuses to launch on a loaded host" begin
    cfg = pa_cfg(; m = 1)
    orb = pa_orb(cfg)
    err = try
        solve_modes_parallel(orb, cfg, [1]; iterations = 0, load_guard = 1.0e-6)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("load", lowercase(msg))
    @test occursin("1.0e-06", msg) || occursin("1.0e-6", msg) || occursin("1e-06", msg)
    @test occursin(string(Sys.CPU_THREADS), msg)          # the core count
    # the SERIAL path is guarded by the same rule (one process' share)
    @test_throws ArgumentError solve_modes_serial(orb, cfg, [1]; iterations = 0,
                                                  load_guard = 1.0e-6)
    # a guard above 1 is meaningless — the owner policy caps it at 0.90
    @test_throws ArgumentError solve_modes_parallel(orb, cfg, [1]; iterations = 0,
                                                    load_guard = 1.5)
    @test_throws ArgumentError solve_modes_parallel(orb, cfg, Int[]; iterations = 0)
    @test_throws ArgumentError solve_modes_parallel(orb, cfg, [1, 1]; iterations = 0)
    # the POOL is guarded by the same rule, before a single process is started
    @test_throws ArgumentError mode_pool(4; load_guard = 1.0e-6)
    @test_throws ArgumentError mode_pool(4; load_guard = 1.5)
    @test_throws ArgumentError mode_pool(0)
    # nothing above may leak a worker process
    @test Distributed.nprocs() == pa_nprocs0
end

# ============ (4) the sweep: completion, order, isolation, budget, load ==================
# ONE POOL for every parallel gate below. MEASURED (proto1/proto3): a FRESH worker process
# pays 166 s of JIT for the level-0 solve and 45 s more for the first AMR step, against 1.3 s
# for a warm 2-level sweep — so a file that started a pool per gate would spend a quarter of
# an hour recompiling the same stack. `mode_pool` is the seam that fixes it, and it is the
# same seam the figure drivers of rows 50/51/53 need.
#
# The pool is sized from the OWNER POLICY, not from a wish: 4 worker processes, each with the
# SAME Julia thread count as this process (so gate (5) compares like with like), and a budget
# that makes the over-request visible. 4 × 4 = 16 ≤ 24 CPUs in the worst case.
pa_nth = min(Threads.nthreads(), PARALLEL_MAX_THREADS_PER_WORKER)
pa_pool = @test_logs (:warn,) match_mode = :any mode_pool(64; threads_per_worker = pa_nth,
                                                          budget = 4 * pa_nth)
try

@testset "(4) the fan-out completes, in order, isolated, inside the budget" begin
    cfg = pa_cfg(; m = 0)
    orb = pa_orb(cfg)
    ms = [0, 1, 2, 5, 20]          # m = 0 exercises the no-outer-BC branch, m = 20 the top mode
    logfile = joinpath(mktempdir(), "load_fanout.log")
    set = solve_modes_parallel(orb, cfg, ms; iterations = 1, pool = pa_pool,
                               log = logfile, sample_interval = 1.0)
    # (a) every requested mode, collected in the REQUESTED order, no failures
    @test set.ms == ms
    @test isempty(set.failures)
    @test length(set.results) == length(ms)
    @test [r.m for r in set.results] == ms
    @test all(length(r.levels) == 2 for r in set.results)      # level 0 + 1 AMR step
    @test all(r.levels[end].record.final_residual < 1.0e-10 for r in set.results)
    @test all(isfinite(abs(r.levels[end].dpsi_r)) for r in set.results)
    @test all(r.x_p[1] == cfg.orbit.orbital_radius for r in set.results)
    @test set[0] === set.results[1]
    @test set[20].m == 20
    @test haskey(set, 5) && !haskey(set, 4)
    @test_throws KeyError set[4]
    # the per-mode config really was rebuilt per mode: m = 0 has no outer BC, m = 20 does
    @test mode_config(cfg, 0).domain.outer_radial_bc === nothing
    @test mode_config(cfg, 20).domain.outer_radial_bc.m_mode_number == 20

    # (c) the budget is ENFORCED, and the over-request was clamped and logged
    @test set.run.nworkers == 4
    @test set.run.requested_workers == 64
    @test set.run.clamped
    @test set.run.threads_per_worker == pa_nth
    @test set.run.nworkers * set.run.threads_per_worker <= PARALLEL_CPU_BUDGET
    @test set.run.budget <= PARALLEL_CPU_BUDGET
    @test !set.run.owned_pool
    @test set.run.blas_threads == 1

    # (f) NO SHARED STATE. 5 modes on 4 processes: every mode built its own Schwarz subdomain
    # matrices and its own bottom LU from scratch, and no mode's level-0 solve reused a cache
    # entry — including the two modes that ran in the SAME process one after the other.
    @test all(r.ospid != getpid() for r in set.results)
    @test all(r.worker != 1 for r in set.results)
    @test length(unique(r.ospid for r in set.results)) == 4     # one process ran two modes
    # a mode that had INHERITED a cache would show `builds == 0` and pure hits; every mode,
    # including the second one on the process that ran two, builds both caches from scratch.
    @test all(r.levels[1].record.subdomain_builds > 0 for r in set.results)
    @test all(r.levels[1].record.bottom_builds > 0 for r in set.results)
    @test all(r.blas_threads == 1 for r in set.results)
    @test all(r.julia_threads == pa_nth for r in set.results)

    # (d) the load monitor sampled the LIVE system during the run and it stayed under the guard
    @test length(set.run.samples) >= 2
    @test all(s.cpu_fraction < PARALLEL_LOAD_GUARD for s in set.run.samples)
    @test set.run.max_load == maximum(s.cpu_fraction for s in set.run.samples)
    @test set.run.max_load < PARALLEL_LOAD_GUARD
    @test set.run.guard_before isa LoadSample && set.run.guard_after isa LoadSample
    @test set.run.guard_before.ncpu == Sys.CPU_THREADS
    @test set.run.planned_share ≈ set.run.nworkers * set.run.threads_per_worker / Sys.CPU_THREADS
    @test isfile(logfile)
    txt = read(logfile, String)
    @test occursin("cpu_fraction", txt)
    @test occursin(string(PARALLEL_CPU_BUDGET), txt)
    @test count(==('\n'), txt) >= length(set.run.samples) + 3   # 3 header lines

    # (e) the cost record exists and is coherent
    @test set.run.wall_time > 0
    @test set.run.startup_time == 0.0            # a REUSED pool pays no startup
    @test all(r.wall_time > 0 for r in set.results)
    @test sum(r.wall_time for r in set.results) > 0
    @test set.run.parallel
    txt = sprint(show, MIME"text/plain"(), set)
    @test occursin("ModeSet", txt) && occursin("parallel", txt)
end

# =================== (5) bitwise identity against the serial oracle ========================
# GATE (b): the fan-out changes WHERE a mode is solved, never WHAT it computes. At the same
# BLAS thread count and the same Julia thread count the process result is BIT-IDENTICAL to the
# in-process one — the final `sol.u` and every AmrHistory diagnostic.
# MEASURED on the prototype: worker-vs-worker and worker-vs-master both bitwise `true`.
@testset "(5) parallel == serial, bitwise" begin
    cfg = pa_cfg(; m = 0)
    orb = pa_orb(cfg)
    ms = [0, 1, 2]
    par = solve_modes_parallel(orb, cfg, ms; iterations = 2, pool = pa_pool)
    ser = solve_modes_serial(orb, cfg, ms; iterations = 2)
    @test isempty(par.failures) && isempty(ser.failures)
    @test [r.m for r in ser.results] == ms
    @test !ser.run.parallel && ser.run.nworkers == 1 && !ser.run.owned_pool
    @test all(r.worker == 1 for r in ser.results)
    @test all(r.ospid == getpid() for r in ser.results)
    @test all(r.julia_threads == Threads.nthreads() for r in ser.results)
    for (p, s) in zip(par.results, ser.results)
        @test pa_bitwise(p, s)
    end
    # the serial oracle also records the load it caused, through the same guard and monitor
    @test !isempty(ser.run.samples)
    @test ser.run.max_load < PARALLEL_LOAD_GUARD
    @test ser.run.guard_before isa LoadSample && ser.run.guard_after isa LoadSample
end

# ====================== (6) the BLAS thread count is the determinism knob ==================
# GATE (b), second half: a DIFFERENT BLAS thread count does NOT reproduce bitwise — it moves
# the answer at ROUNDOFF only. MEASURED on the prototype at 2 AMR levels: `u` not bitwise,
# max relative Δu 3.2e-12, ΔΨ^R 7.5e-13, Δ∂_rΨ^R 2.1e-13 (the row-43 node measured 4.2e-12 on
# the 8-level sweep). This is why the pin is normative and why every worker sets it EXPLICITLY
# instead of trusting the environment.
@testset "(6) a different BLAS thread count breaks bitwise identity at roundoff only" begin
    cfg = pa_cfg(; m = 1)
    orb = pa_orb(cfg)
    one_thread = solve_modes_parallel(orb, cfg, [1]; iterations = 2, pool = pa_pool,
                                      blas_threads = 1)
    two_thread = solve_modes_parallel(orb, cfg, [1]; iterations = 2, pool = pa_pool,
                                      blas_threads = 2)
    @test isempty(one_thread.failures) && isempty(two_thread.failures)
    a, b = one_thread.results[1], two_thread.results[1]
    @test a.blas_threads == 1
    @test b.blas_threads == 2
    @test one_thread.run.blas_threads == 1 && two_thread.run.blas_threads == 2
    # not bitwise …
    @test !all(a.u .=== b.u)
    # … but the difference is roundoff, three orders below the level-to-level AMR change
    @test maximum(pa_rel.(a.u, b.u)) < 1.0e-9
    @test pa_rel(a.levels[end].dpsi_r, b.levels[end].dpsi_r) < 1.0e-9
    @test pa_rel(a.levels[end].psi, b.levels[end].psi) < 1.0e-9
    @test a.levels[end].record.iterations == b.levels[end].record.iterations
    # … and the same pinned count reproduces itself bitwise, which is the point of the pin
    again = solve_modes_parallel(orb, cfg, [1]; iterations = 2, pool = pa_pool,
                                 blas_threads = 1)
    @test pa_bitwise(a, again.results[1])
end

# ============================== (7) the failure path =======================================
# GATE (g): a mode that throws is reported BY ID and never swallowed; every other mode still
# returns. The failure is injected through the row-43 `criterion` seam — the only
# caller-supplied callable on the solve path — so it fires INSIDE the remote solve.
@testset "(7) one failing mode is reported by id, the others still return" begin
    cfg = pa_cfg(; m = 0)
    orb = pa_orb(cfg)
    boom = sol -> sol.m == 1 ? error("injected failure in mode $(sol.m)") : amr_flags(sol)
    set = @test_logs (:warn,) match_mode = :any solve_modes_parallel(
        orb, cfg, [0, 1, 2]; iterations = 0, pool = pa_pool, criterion = boom)
    @test length(set.failures) == 1
    @test set.failures[1].m == 1
    @test occursin("injected failure", set.failures[1].message)
    @test set.failures[1].worker != 1
    @test [r.m for r in set.results] == [0, 2]
    @test set[0].m == 0 && set[2].m == 2
    @test_throws KeyError set[1]
    @test occursin("FAILED", sprint(show, MIME"text/plain"(), set))
    # the serial oracle reports the same failure the same way, by id, and keeps going
    sset = @test_logs (:warn,) match_mode = :any solve_modes_serial(
        orb, cfg, [0, 1, 2]; iterations = 0, criterion = boom)
    @test [f.m for f in sset.failures] == [1]
    @test occursin("injected failure", sset.failures[1].message)
    @test [r.m for r in sset.results] == [0, 2]
end

finally
    close(pa_pool)
end
@test !isopen(pa_pool)
@test length(pa_pool) == 0
@test Distributed.nprocs() == pa_nprocs0          # no worker process leaked

end # testset mmode-parallel-driver

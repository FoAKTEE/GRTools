# parallel.jl — elliptic-ssf node `mmode-parallel-driver` (GLOBAL_DAG row 47, DESIGN §5.1).
#
# THE PROCESS FAN-OUT over the paper's m = 0…20 sweep (main.tex:716 "Each m-mode is
# independent … each mode solves in a few seconds"). Five products, in dependency order:
#
#   1. `system_load` / `LoadSample`  — the LIVE whole-system reading the owner's compute
#      policy is written in terms of (Σ `ps -eo pcpu` / ncpu, plus the load average).
#   2. `mode_pool_size`              — the pool arithmetic: `min(budget ÷ threads, modes)`,
#      every clamp reported, `nworkers × threads_per_worker ≤ 24` ALWAYS.
#   3. `mode_config`                 — the per-mode `RunConfig` (only `MModeNumber` and the
#      m = 0 no-outer-BC yaml branch may move).
#   4. `ModePool`                    — a guarded, reusable pool of worker processes.
#   5. `solve_modes_parallel` / `solve_modes_serial` → `ModeSet` — the fan-out and its
#      in-process BITWISE ORACLE.
#
# WHY PROCESSES AND NOT THREADS (sealed, row 42/43): one `AmrHistory` per PROCESS. The solver
# stack and its two caches (Schwarz subdomain matrices, bottom LU) are neither thread-safe nor
# reusable across modes, BigFloat GC contention makes `-t16` net slower than `-t4`, and each
# mode fits comfortably in one core.
#
# DEPENDENCY-BUDGET DECISION (orchestrator-approved, wave 12; mirrors adopted decision 3 on
# SparseArrays): the STDLIB `Distributed` is added to the ScalarSelfForce budget, and
# `test/hygiene.jl` is widened in the same commit. It is a stdlib, not a registry package, and
# it is the KISS choice: the alternative — spawning `julia --project` children with `Base.run`
# and writing results to files — buys zero new dependencies at the price of a hand-rolled
# binary serialization of `AmrLevel`/`Mesh2D`/`SolveRecord` plus child-stderr parsing for the
# failure path, i.e. more code and a second thing to get wrong. `Distributed` already
# transfers exactly the payload this node needs and surfaces a remote throw as an exception
# object carrying its worker id.

# ==========================================================================================
#  THE COMPUTE POLICY (the ONE normative statement; quoted from the binding `compute-usage`
#  skill, which is owner policy for this host — 256 cores, shared with other users)
# ==========================================================================================
#   * CPU: "use no more than 50 % of all CPUs in principle; in practice use middle,
#     non-occupied cores and AT MOST 24 CPUs for our jobs."
#   * Whole-system guard: "before AND after launching jobs, total CPU usage of ALL users must
#     stay BELOW 90 %."
#   * Threads: "prefer N ≤ 4 threads for BigFloat-heavy work … scale BigFloat via processes,
#     not threads."
#   * Sweeps: "for embarrassingly parallel sweeps (e.g. independent m-modes), size the worker
#     pool as min(24, modes) and leave it there; NEVER auto-scale up."
#   * GPU: "banned by admin. Never."
# Everything below enforces those five sentences mechanically.

"The owner's hard CPU cap for this mission's jobs (`compute-usage` skill). Never exceeded."
const PARALLEL_CPU_BUDGET = 24

"The owner's cap on Julia threads per process (`compute-usage` skill: `N ≤ 4`)."
const PARALLEL_MAX_THREADS_PER_WORKER = 4

"The owner's whole-system CPU guard (`compute-usage` skill: total usage of ALL users < 90 %)."
const PARALLEL_LOAD_GUARD = 0.90

# ================================= the live load reading ==================================

"""
    LoadSample

One reading of the whole-system CPU usage, the quantity the owner's 90 % guard is written in
terms of. Fields:

  - `time::Float64` — seconds since the run's origin (`0` is the pre-launch guard sample).
  - `cpu_fraction::Float64` — total CPU usage of ALL users as a fraction of the machine,
    measured as the BUSY fraction of a short `/proc/stat` window (`1 − Δidle/Δtotal`). THE
    guarded quantity, and bounded by construction to `[0, 1]`.
  - `ps_fraction::Float64` — the `compute-usage` skill's own shell formula,
    `Σ ps -eo pcpu / (100 · ncpu)`, recorded beside it. MEASURED to agree with `cpu_fraction`
    on a quiet host (0.0456 against 0.0400) but NOT bounded: `ps` reports each process'
    LIFETIME average, so a burst of short-lived processes that each averaged ~100 % sums far
    above the machine — this node measured **2.868** (287 % of a 256-core host) during a
    compile storm, which is why it is no longer the guarded quantity.
  - `loadavg1::Float64` — the 1-minute load average divided by `ncpu` (`uptime`'s number). It
    lags by up to a minute, so it guards nothing; it is what a human reads.
  - `ncpu::Int` — `Sys.CPU_THREADS`.
  - `source::Symbol` — `:procstat` normally; `:ps` if `/proc/stat` was unreadable; `:loadavg`
    if neither was. The two fallbacks are WEAKER readings, and the source travels with the
    sample so a run can be audited.

See [`system_load`](@ref).
"""
struct LoadSample
    time::Float64
    cpu_fraction::Float64
    ps_fraction::Float64
    loadavg1::Float64
    ncpu::Int
    source::Symbol
end

function _pf_ps_total()
    try
        out = read(`ps -eo pcpu --no-headers`, String)
        s = 0.0
        n = 0
        for line in eachsplit(out, '\n')
            v = tryparse(Float64, strip(line))
            if v !== nothing
                s += v
                n += 1
            end
        end
        return n == 0 ? nothing : s
    catch
        return nothing
    end
end

# the cumulative jiffy counters of /proc/stat's aggregate `cpu` line: (busy, total).
# NOTE: procfs reports size 0, and `read(path, String)` trusts `filesize` and returns "" —
# MEASURED here, and the reason this streams the first line off an open handle instead.
function _pf_procstat()
    try
        line = open(readline, "/proc/stat")
        startswith(line, "cpu ") || return nothing
        f = [something(tryparse(Int64, w), Int64(0)) for w in eachsplit(strip(line)) if
             tryparse(Int64, w) !== nothing]
        length(f) >= 4 || return nothing
        total = sum(f)
        idle = f[4] + (length(f) >= 5 ? f[5] : 0)     # idle + iowait
        return (total - idle, total)
    catch
        return nothing
    end
end

"The width of the `/proc/stat` window `system_load` averages the busy fraction over."
const PARALLEL_LOAD_WINDOW = 0.15

"""
    system_load(t0::Real = time(); window = PARALLEL_LOAD_WINDOW) -> LoadSample

The LIVE whole-system CPU usage, read fresh at every call — the pre-launch guard, the
during-run monitor and the post-run check all go through this one function. `t0` is the origin
the sample's `time` is measured from.

`cpu_fraction` is the BUSY fraction of a `window`-second `/proc/stat` interval,
`1 − Δ(idle + iowait)/Δtotal`, over the aggregate `cpu` line: every user's time, not only ours
— this host is shared and the guard is a whole-system guard — and bounded to `[0, 1]` by
construction, which is what makes it safe to compare against `PARALLEL_LOAD_GUARD`.

`ps_fraction` records the `compute-usage` skill's own shell formula
`Σ ps -eo pcpu / (100 · ncpu)` beside it. WHY BOTH, AND WHY THE SWITCH: `ps` prints each
process' LIFETIME average, so the sum is not an instantaneous utilisation and is not bounded by
the machine. On a quiet host the two agree (MEASURED 0.0456 against 0.0400); during a compile
storm on this 256-core box the `ps` sum reached **2.868**, i.e. 287 % of the machine, and would
have tripped a 0.90 gate that reality never crossed. The skill's number is kept because it is
the number a human checks at the shell, but the guard is taken on the bounded one.

`source` is `:procstat` normally, `:ps` if `/proc/stat` could not be read, `:loadavg` if
neither was available — the last two are weaker readings and travel labelled.
"""
function system_load(t0::Real = time(); window::Real = PARALLEL_LOAD_WINDOW)
    ncpu = Sys.CPU_THREADS
    la = Sys.loadavg()[1] / ncpu
    total = _pf_ps_total()
    psf = total === nothing ? NaN : total / (100 * ncpu)
    a = _pf_procstat()
    if a !== nothing
        sleep(window)
        b = _pf_procstat()
        if b !== nothing && b[2] > a[2]
            frac = clamp((b[1] - a[1]) / (b[2] - a[2]), 0.0, 1.0)
            return LoadSample(time() - Float64(t0), frac, psf, la, ncpu, :procstat)
        end
    end
    t = time() - Float64(t0)
    return total === nothing ? LoadSample(t, la, psf, la, ncpu, :loadavg) :
                               LoadSample(t, psf, psf, la, ncpu, :ps)
end

# The pre-launch refusal. `planned_share` is the fraction of the machine THIS job is about to
# take (`nworkers × threads_per_worker / ncpu`). Refuses with every number in the message.
function _pf_guard(planned_share::Real, load_guard::Real, t0::Real, what::AbstractString)
    (load_guard isa Real && 0 < load_guard <= PARALLEL_LOAD_GUARD) || throw(ArgumentError(
        "$what: load_guard must lie in (0, $PARALLEL_LOAD_GUARD] — the owner's whole-system " *
        "guard is $(PARALLEL_LOAD_GUARD) (compute-usage skill), got $load_guard"))
    s = system_load(t0)
    total = s.cpu_fraction + planned_share
    total <= load_guard || throw(ArgumentError(
        "$what: REFUSING TO LAUNCH — the whole-system CPU guard would be exceeded. " *
        "live usage $(round(s.cpu_fraction; digits = 4)) (source :$(s.source), " *
        "ps sum $(round(s.ps_fraction; digits = 4)), " *
        "loadavg1/ncpu $(round(s.loadavg1; digits = 4)), ncpu $(s.ncpu)) " *
        "+ planned share $(round(planned_share; digits = 4)) = $(round(total; digits = 4)) " *
        "> load_guard $(load_guard). Wait for the host to quieten or shrink the pool."))
    return s
end

# ================================= the pool arithmetic =====================================

"""
    mode_pool_size(nmodes; nworkers = nothing, threads_per_worker = 1,
                   budget = PARALLEL_CPU_BUDGET)
        -> (; nworkers, requested, threads_per_worker, budget, clamped)

The owner's pool rule, as arithmetic: **size the pool as `min(budget ÷ threads_per_worker,
nmodes)` and leave it there; never auto-scale up.**

`nworkers === nothing` (the default) asks for that rule. An explicit `nworkers` is a REQUEST:
it is clamped to `min(budget ÷ threads_per_worker, nmodes)` — a pool bigger than the budget
would break owner policy, and one bigger than the mode count would only start idle processes.
`threads_per_worker` is clamped to `PARALLEL_MAX_THREADS_PER_WORKER` and `budget` to
`PARALLEL_CPU_BUDGET`; the returned `clamped` flag is `true` whenever ANY of the three was
reduced, and the callers log it — a clamp is never silent.

The post-condition `nworkers × threads_per_worker ≤ PARALLEL_CPU_BUDGET` holds for every
input, including hostile ones.

    julia> mode_pool_size(21)                       # the paper's sweep
    (nworkers = 21, requested = 21, threads_per_worker = 1, budget = 24, clamped = false)

    julia> mode_pool_size(100; nworkers = 64)       # an over-request, clamped and reported
    (nworkers = 24, requested = 64, threads_per_worker = 1, budget = 24, clamped = true)
"""
function mode_pool_size(nmodes::Integer; nworkers = nothing, threads_per_worker::Integer = 1,
                        budget::Integer = PARALLEL_CPU_BUDGET)
    nmodes >= 1 || throw(ArgumentError("mode_pool_size: nmodes must be ≥ 1, got $nmodes"))
    threads_per_worker >= 1 || throw(ArgumentError(
        "mode_pool_size: threads_per_worker must be ≥ 1, got $threads_per_worker"))
    budget >= 1 || throw(ArgumentError("mode_pool_size: budget must be ≥ 1, got $budget"))
    nth = min(Int(threads_per_worker), PARALLEL_MAX_THREADS_PER_WORKER)
    bud = min(Int(budget), PARALLEL_CPU_BUDGET)
    cap = max(1, bud ÷ nth)
    req = nworkers === nothing ? min(cap, Int(nmodes)) : Int(nworkers)
    req >= 1 || throw(ArgumentError("mode_pool_size: nworkers must be ≥ 1, got $req"))
    n = clamp(req, 1, min(cap, Int(nmodes)))
    clamped = n != req || nth != Int(threads_per_worker) || bud != Int(budget)
    return (nworkers = n, requested = req, threads_per_worker = nth, budget = bud,
            clamped = clamped)
end

# ================================== the per-mode config ====================================

# one field of an immutable struct replaced, everything else carried over verbatim
_pf_replace(x::T, name::Symbol, v) where {T} =
    T((f === name ? v : getfield(x, f) for f in fieldnames(T))...)

"""
    mode_config(config::RunConfig, m::Integer) -> RunConfig

The `RunConfig` of mode `m`, derived from a template config. The fan-out ships ONE config to
the workers and each worker rebuilds its own mode's config with this function, so the sweep is
described by a single object and no yaml parameter can drift between modes.

EXACTLY TWO things depend on `m` in the CRP input (`anc/ScalarSelfForce.yaml`), and this
function moves exactly those two:

  - `orbit.m_mode_number` — `Background.CircularOrbit.MModeNumber` (yaml:29);
  - `domain.outer_radial_bc` — the Jinja2 branch `{% if m_mode == 0 %} None {% else %}
    Sommerfeld …` (yaml:114-125): `nothing` at `m = 0`, otherwise the order-2 hyperboloidal
    Sommerfeld condition carrying the same `m`.

Every other field — block bounds, refinement, DG, GMRES, multigrid, Schwarz, AMR, observers —
is carried over verbatim, so a non-paper template keeps its customizations. If the template is
itself an `m = 0` config (no outer BC to copy) the condition is synthesized from the template's
own `(M, a, r₀)` with the yaml's `HyperboloidalSlicing: True` and `Order: 2`.

`m` is the FOLDED mode index of the package-wide `m ≥ 0` ×2 real-part convention; a negative
`m` is an error, not a retrograde orbit (that is `a < 0`).
"""
function mode_config(config::RunConfig{RT}, m::Integer) where {RT}
    m >= 0 || throw(ArgumentError(
        "mode_config: m must be ≥ 0 (the package's folded m ≥ 0 ×2 real-part convention — " *
        "a retrograde orbit is a < 0, not m < 0), got $m"))
    mm = Int(m)
    mm == config.orbit.m_mode_number && return config
    orbit = _pf_replace(config.orbit, :m_mode_number, mm)
    bc = config.domain.outer_radial_bc
    newbc = if mm == 0
        nothing
    elseif bc === nothing
        SommerfeldBC{RT}(config.orbit.black_hole_mass, config.orbit.black_hole_spin,
                         config.orbit.orbital_radius, mm, true, 2)
    else
        _pf_replace(bc, :m_mode_number, mm)
    end
    domain = _pf_replace(config.domain, :outer_radial_bc, newbc)
    return _pf_replace(_pf_replace(config, :orbit, orbit), :domain, domain)
end

# ==================================== the run record =======================================

"""
    ModeResult{T,RT}

One mode's TRANSFERABLE result: everything of an [`AmrHistory`](@ref) that crosses a process
boundary, and everything rows 44/50/51/53 read.

  - `m::Int`, `x_p::NTuple{2,RT}` — the mode and the puncture in the discretization chart.
  - `levels::Vector{AmrLevel{T,RT}}` — the sweep's per-level record VERBATIM (mesh, puncture
    elements, `SolveRecord`, and the particle diagnostics `psi`, `dpsi_r`, `dpsi_w`), level 0
    first. `levels[end]` carries the converged values row 44 turns into `F_r^m`.
  - `u::Vector{T}` — the FINAL level's solution vector, `hist.solution.u`.
  - `wall_time::Float64` — seconds this mode took inside its process, INCLUDING the first
    call's JIT when it was the first mode that process solved (see [`ModePool`](@ref)).
  - `worker::Int`, `ospid::Int` — the Distributed id (`1` = in-process) and the OS pid.
  - `blas_threads::Int`, `julia_threads::Int` — what the producing process actually had. These
    are not decoration: a run is bitwise reproducible only against another run with the SAME
    two numbers (MEASURED: the diagnostics move by ~4e-12 relative between BLAS thread counts).

WHAT IS DELIBERATELY NOT HERE: the live `ModeSolution` — its operator closures and its two
solver caches (Schwarz subdomain matrices, bottom LU) are hundreds of megabytes, are
process-local by row 42's own rule, and are useless to a different process. A caller that
needs the live object calls the sealed `solve_mode_amr` in its own process. MEASURED payload
of a full 8-level mode: **597 KB**, against ~2.5 GB of process memory.
"""
struct ModeResult{T<:Number,RT<:AbstractFloat}
    m::Int
    x_p::NTuple{2,RT}
    levels::Vector{AmrLevel{T,RT}}
    u::Vector{T}
    wall_time::Float64
    worker::Int
    ospid::Int
    blas_threads::Int
    julia_threads::Int
end

ModeResult(h::AmrHistory{T,RT}, wall_time::Real) where {T,RT} =
    ModeResult{T,RT}(h.m, h.x_p, h.levels, h.solution.u, Float64(wall_time),
                     Distributed.myid(), getpid(), BLAS.get_num_threads(), Threads.nthreads())

"""
    ModeFailure

One mode that threw: its `m`, the rendered exception `message` (the remote backtrace included
when it came from a worker), and the `worker` it happened on (`0` when the id is unknown).
A failure is NEVER swallowed — it is `@warn`ed as it happens, kept here, and printed by
[`ModeSet`](@ref)'s `show`; a caller that requires every mode must assert
`isempty(set.failures)`.
"""
struct ModeFailure
    m::Int
    message::String
    worker::Int
end

"""
    ModeRunRecord

What a sweep DID, as opposed to what it computed: the pool it used, the budget arithmetic that
sized it, the guard decisions, and the load-monitor samples that are this node's DAG evidence.

  - `parallel::Bool`, `nworkers::Int`, `threads_per_worker::Int`, `blas_threads::Int`;
  - `budget::Int`, `requested_workers::Int`, `clamped::Bool` — the [`mode_pool_size`](@ref)
    arithmetic as applied, with `clamped` recording that a request was reduced (and logged);
  - `owned_pool::Bool` — whether the sweep started and stopped its own [`ModePool`](@ref);
  - `load_guard::Float64`, `planned_share::Float64` — the guard and the fraction of the
    machine this run planned to take;
  - `guard_before::LoadSample`, `guard_after::LoadSample` — the `compute-usage` skill's
    BEFORE and AFTER readings;
  - `max_load::Float64` — the largest `cpu_fraction` over every sample of the run. THE gate:
    it must stay under `PARALLEL_LOAD_GUARD`;
  - `startup_time::Float64` — seconds spent starting the pool (spawn + `using
    ScalarSelfForce`); `0` when an existing pool was reused;
  - `wall_time::Float64` — seconds for the whole call;
  - `samples::Vector{LoadSample}` — the monitor's series (also written to `log`).
"""
struct ModeRunRecord
    parallel::Bool
    nworkers::Int
    threads_per_worker::Int
    blas_threads::Int
    budget::Int
    requested_workers::Int
    clamped::Bool
    owned_pool::Bool
    load_guard::Float64
    planned_share::Float64
    guard_before::LoadSample
    guard_after::LoadSample
    max_load::Float64
    startup_time::Float64
    wall_time::Float64
    samples::Vector{LoadSample}
end

"""
    ModeSet{T,RT}

A whole m-mode sweep: the `orbit`, the template `config`, the requested modes `ms`, the
successful `results` (in `ms` order, one [`ModeResult`](@ref) each), the `failures`, and the
[`ModeRunRecord`](@ref).

`set[m]` is the result of mode `m` (`KeyError` if it failed or was never requested);
`haskey(set, m)` asks first. Build it with [`solve_modes_parallel`](@ref) or its bitwise
oracle [`solve_modes_serial`](@ref).
"""
struct ModeSet{T<:Number,RT<:AbstractFloat}
    orbit::CircularOrbit
    config::RunConfig{RT}
    ms::Vector{Int}
    results::Vector{ModeResult{T,RT}}
    failures::Vector{ModeFailure}
    run::ModeRunRecord
end

function Base.getindex(s::ModeSet, m::Integer)
    i = findfirst(r -> r.m == Int(m), s.results)
    i === nothing && throw(KeyError(Int(m)))
    return s.results[i]
end
Base.haskey(s::ModeSet, m::Integer) = any(r -> r.m == Int(m), s.results)

function Base.show(io::IO, ::MIME"text/plain", s::ModeSet)
    r = s.run
    println(io, "ModeSet: ", length(s.results), " of ", length(s.ms), " modes, a = ",
            s.orbit.a, ", r₀ = ", s.orbit.r0,
            r.parallel ? " [parallel]" : " [serial]")
    @printf(io, "  %d worker(s) × %d thread(s), BLAS %d, budget %d (requested %d%s) — %.1f s wall, %.1f s startup\n",
            r.nworkers, r.threads_per_worker, r.blas_threads, r.budget, r.requested_workers,
            r.clamped ? ", CLAMPED" : "", r.wall_time, r.startup_time)
    @printf(io, "  system load: before %.4f, max %.4f, after %.4f (guard %.2f, %d samples)\n",
            r.guard_before.cpu_fraction, r.max_load, r.guard_after.cpu_fraction,
            r.load_guard, length(r.samples))
    println(io, "    m  levels   ndofs  iters   wall[s]  worker      Re Ψ^R(x_p)     ∂_r Ψ^R(x_p)")
    for x in s.results
        l = x.levels[end]
        @printf(io, "  %3d  %6d  %6d  %5d  %8.2f  %6d  %+13.6e  %+13.6e\n",
                x.m, length(x.levels), l.mesh.ndofs, l.record.iterations, x.wall_time,
                x.ospid, real(l.psi), real(l.dpsi_r))
    end
    for f in s.failures
        println(io, "  FAILED m = ", f.m, " on worker ", f.worker, ": ",
                first(f.message, 200))
    end
end

# ================================== the load monitor =======================================

mutable struct _PfMonitor
    t0::Float64
    interval::Float64
    io::Union{Nothing,IO}
    own_io::Bool
    samples::Vector{LoadSample}
    running::Bool
    task::Union{Nothing,Task}
end

function _pf_monitor(t0::Real, interval::Real, log)
    interval > 0 || throw(ArgumentError("sample_interval must be > 0, got $interval"))
    io = log === nothing ? nothing : (log isa IO ? log : open(String(log), "w"))
    mon = _PfMonitor(Float64(t0), Float64(interval), io, !(log === nothing || log isa IO),
                     LoadSample[], false, nothing)
    if io !== nothing
        println(io, "# ScalarSelfForce mmode-parallel-driver load monitor")
        println(io, "# owner compute policy (compute-usage skill): ≤ $PARALLEL_CPU_BUDGET CPUs for our jobs, ",
                "whole system < $(PARALLEL_LOAD_GUARD), ncpu = ", Sys.CPU_THREADS)
        println(io, "t_seconds\tcpu_fraction\tps_fraction\tloadavg1_per_cpu\tsource")
        flush(io)
    end
    return mon
end

function _pf_sample!(mon::_PfMonitor, s::LoadSample = system_load(mon.t0))
    push!(mon.samples, s)
    if mon.io !== nothing
        @printf(mon.io, "%.3f\t%.6f\t%.6f\t%.6f\t%s\n",
                s.time, s.cpu_fraction, s.ps_fraction, s.loadavg1, s.source)
        flush(mon.io)
    end
    return s
end

# The monitor runs as a TASK on the master. It samples while the master is parked on the
# workers' channels (the whole fan-out), which is exactly when the machine is loaded. It does
# NOT sample while the master itself computes, so `solve_modes_serial` drives it by hand
# between modes instead of starting it.
function _pf_start!(mon::_PfMonitor)
    mon.running = true
    mon.task = @async begin
        while mon.running
            try
                _pf_sample!(mon)
            catch
            end
            slept = 0.0
            while mon.running && slept < mon.interval
                sleep(min(0.2, mon.interval))
                slept += min(0.2, mon.interval)
            end
        end
    end
    return mon
end

function _pf_stop!(mon::_PfMonitor)
    mon.running = false
    if mon.task !== nothing
        try
            wait(mon.task)
        catch
        end
        mon.task = nothing
    end
    return mon
end

_pf_close!(mon::_PfMonitor) = (mon.own_io && mon.io !== nothing && close(mon.io); nothing)

# ====================================== the pool ===========================================

"""
    ModePool

A live pool of worker PROCESSES for the m-mode fan-out: `pids` (Distributed ids),
`threads_per_worker`, the budget arithmetic that sized it (`requested`, `budget`, `clamped`),
the `project` the workers were started in, and `startup_time`.

WHY IT IS A FIRST-CLASS OBJECT. A fresh Julia process pays the WHOLE solver stack's JIT before
it can solve anything: MEASURED on this host, **166 s** for the level-0 solve and a further
**45 s** for the first AMR step, against **15.9 s** for a warm 8-level sweep. A sweep that
starts and stops its own pool therefore pays ~211 s of compilation no matter how small it is,
while a REUSED pool pays it once and then runs each further sweep at full speed. Callers that
run several sweeps — the figure drivers of rows 50/51/53 — must hold one pool across all of
them. Build it with [`mode_pool`](@ref) and always `close` it (the `do`-block form does that
for you).
"""
mutable struct ModePool
    pids::Vector{Int}
    threads_per_worker::Int
    budget::Int
    requested::Int
    clamped::Bool
    project::String
    startup_time::Float64
    open::Bool
end

Base.isopen(p::ModePool) = p.open
Base.length(p::ModePool) = length(p.pids)

function Base.show(io::IO, p::ModePool)
    print(io, "ModePool(", length(p.pids), " worker", length(p.pids) == 1 ? "" : "s",
          " × ", p.threads_per_worker, " thread", p.threads_per_worker == 1 ? "" : "s",
          p.open ? "" : ", CLOSED", ")")
end

function _pf_project()
    p = Base.active_project()
    p === nothing && throw(ArgumentError(
        "mode_pool: this session has no active project, so a worker could not load " *
        "ScalarSelfForce — pass `project = \"/path/to/Project.toml\"` explicitly"))
    return p
end

"""
    mode_pool(nworkers; threads_per_worker = 1, blas_threads = 1,
              budget = PARALLEL_CPU_BUDGET, load_guard = PARALLEL_LOAD_GUARD,
              project = Base.active_project(), exeflags = String[]) -> ModePool
    mode_pool(f, nworkers; kwargs...)

Start a guarded pool of `nworkers` worker processes, each launched with
`--project=<project> -t<threads_per_worker>` and `OPENBLAS_NUM_THREADS=<blas_threads>`, and
each loading `ScalarSelfForce`.

The request goes through [`mode_pool_size`](@ref), so `nworkers × threads_per_worker` never
exceeds `PARALLEL_CPU_BUDGET`; a reduced request is `@warn`ed, never silently applied. Before
a single process is started the LIVE whole-system usage is read and the pool is REFUSED with
an `ArgumentError` naming every number if `current + nworkers × threads / ncpu` would cross
`load_guard`.

`close(pool)` removes the workers; the `do`-block form does it even on a throw:

    mode_pool(21) do pool
        conv = solve_modes_parallel(orbit, config, 0:20; pool = pool)
        modes = solve_modes_parallel(orbit, config, 0:20; pool = pool, iterations = 4)
    end

BLAS pinning is applied INSIDE each worker at every solve (see [`solve_modes_parallel`](@ref)),
so one pool serves sweeps with different `blas_threads`; the environment variable set here only
sizes OpenBLAS's pool at process start.
"""
function mode_pool(nworkers::Integer; threads_per_worker::Integer = 1,
                   blas_threads::Integer = 1, budget::Integer = PARALLEL_CPU_BUDGET,
                   load_guard::Real = PARALLEL_LOAD_GUARD,
                   project::AbstractString = _pf_project(),
                   exeflags::AbstractVector{<:AbstractString} = String[])
    blas_threads >= 1 || throw(ArgumentError("mode_pool: blas_threads must be ≥ 1, got $blas_threads"))
    p = mode_pool_size(typemax(Int); nworkers = nworkers,
                       threads_per_worker = threads_per_worker, budget = budget)
    p.clamped && @warn "mode_pool: the compute budget clamped the request (compute-usage skill)" requested_workers = p.requested granted_workers = p.nworkers requested_threads = threads_per_worker threads_per_worker = p.threads_per_worker requested_budget = budget budget = p.budget cpu_budget = PARALLEL_CPU_BUDGET
    share = p.nworkers * p.threads_per_worker / Sys.CPU_THREADS
    _pf_guard(share, load_guard, time(),
              "mode_pool($(p.nworkers) × $(p.threads_per_worker) threads)")
    t = time()
    pids = Distributed.addprocs(p.nworkers;
        exeflags = String[String("--project=" * project),
                          String("-t" * string(p.threads_per_worker)), exeflags...],
        env = ["OPENBLAS_NUM_THREADS" => string(blas_threads),
               "JULIA_NUM_THREADS" => string(p.threads_per_worker),
               "JULIA_PKG_OFFLINE" => "true"],
        topology = :master_worker)
    try
        Distributed.remotecall_eval(Main, pids, :(using ScalarSelfForce))
    catch e
        Distributed.rmprocs(pids; waitfor = 60)
        rethrow(e)
    end
    return ModePool(collect(Int, pids), p.threads_per_worker, p.budget, p.requested,
                    p.clamped, String(project), time() - t, true)
end

function Base.close(p::ModePool)
    if p.open
        p.open = false
        isempty(p.pids) || Distributed.rmprocs(p.pids; waitfor = 60)
        empty!(p.pids)
    end
    return p
end

function mode_pool(f::Function, nworkers::Integer; kwargs...)
    p = mode_pool(nworkers; kwargs...)
    try
        return f(p)
    finally
        close(p)
    end
end

# =================================== the mode solve ========================================

# THE function every worker runs. Top-level (not a closure) so it serializes by name, and the
# ONE place the BLAS pin is applied — the environment is best effort, this is authoritative.
function _pf_solve(orbit::CircularOrbit, config::RunConfig, m::Integer;
                   iterations::Integer, T::Type, criterion, level_options,
                   blas_threads::Integer)
    old = BLAS.get_num_threads()
    try
        BLAS.set_num_threads(blas_threads)
        cfg = mode_config(config, m)
        t = time()
        h = solve_mode_amr(orbit, cfg, m; iterations = iterations, T = T,
                           criterion = criterion, level_options = level_options)
        return ModeResult(h, time() - t)
    finally
        BLAS.set_num_threads(old)
    end
end

function _pf_modes(ms)
    v = collect(Int, ms)
    isempty(v) && throw(ArgumentError(
        "solve_modes_*: the mode list is empty (the paper's sweep is CRP_PAPER_MMODES = 0:20)"))
    length(unique(v)) == length(v) || throw(ArgumentError(
        "solve_modes_*: the mode list has duplicates, which would solve the same problem " *
        "twice and break the `set[m]` lookup: $v"))
    return v
end

_pf_msg(e) = try
        first(sprint(showerror, e), 8000)
    catch
        first(string(e), 8000)
    end

_pf_failed_pid(e) = e isa Distributed.RemoteException ? e.pid : 0

# ==================================== the two drivers ======================================

"""
    solve_modes_serial(orbit, config, ms = CRP_PAPER_MMODES; iterations = config.amr.iterations,
                       T = ComplexF64, criterion = amr_flags, level_options = amr_level_options,
                       blas_threads = 1, load_guard = PARALLEL_LOAD_GUARD,
                       log = nothing, sample_interval = 5.0) -> ModeSet

The m-mode sweep run IN PROCESS, one mode after another — the BITWISE ORACLE of
[`solve_modes_parallel`](@ref) and the only thing its determinism gate can be measured
against. Same arguments, same [`ModeSet`](@ref), same failure handling; the only differences
are that `run.parallel` is `false`, every `ModeResult` carries `worker = 1` and this process'
`ospid`, and the load monitor is driven by hand between modes (a busy master never yields to a
sampling task).

It is also the right driver for a SINGLE mode, or for a caller that is itself already one
worker of somebody else's fan-out.
"""
function solve_modes_serial(orbit::CircularOrbit, config::RunConfig{RT},
                            ms = CRP_PAPER_MMODES;
                            iterations::Integer = config.amr.iterations,
                            T::Type = ComplexF64,
                            criterion = amr_flags,
                            level_options = amr_level_options,
                            blas_threads::Integer = 1,
                            load_guard::Real = PARALLEL_LOAD_GUARD,
                            log = nothing,
                            sample_interval::Real = 5.0) where {RT}
    blas_threads >= 1 || throw(ArgumentError(
        "solve_modes_serial: blas_threads must be ≥ 1, got $blas_threads"))
    mlist = _pf_modes(ms)
    t0 = time()
    share = 1 / Sys.CPU_THREADS
    before = _pf_guard(share, load_guard, t0, "solve_modes_serial(1 process)")
    mon = _pf_monitor(t0, sample_interval, log)
    _pf_sample!(mon, before)
    results = ModeResult{T,RT}[]
    failures = ModeFailure[]
    old = BLAS.get_num_threads()
    after = before
    try
        BLAS.set_num_threads(blas_threads)
        for m in mlist
            _pf_sample!(mon)
            try
                push!(results, _pf_solve(orbit, config, m; iterations = iterations, T = T,
                                         criterion = criterion, level_options = level_options,
                                         blas_threads = blas_threads))
            catch e
                @warn "solve_modes_serial: mode m = $m FAILED; the other modes still run" exception = (e, catch_backtrace())
                push!(failures, ModeFailure(m, _pf_msg(e), Distributed.myid()))
            end
        end
    finally
        BLAS.set_num_threads(old)
        after = _pf_sample!(mon)          # the skill's mandatory AFTER reading, throw or not
        _pf_close!(mon)
    end
    rec = ModeRunRecord(false, 1, Threads.nthreads(), Int(blas_threads),
                        PARALLEL_CPU_BUDGET, 1, false, false, Float64(load_guard), share,
                        before, after, maximum(s -> s.cpu_fraction, mon.samples),
                        0.0, time() - t0, mon.samples)
    return ModeSet{T,RT}(orbit, config, mlist, results, failures, rec)
end

"""
    solve_modes_parallel(orbit, config, ms = CRP_PAPER_MMODES;
                         iterations = config.amr.iterations, T = ComplexF64,
                         criterion = amr_flags, level_options = amr_level_options,
                         pool = nothing, nworkers = nothing, threads_per_worker = 1,
                         blas_threads = 1, budget = PARALLEL_CPU_BUDGET,
                         load_guard = PARALLEL_LOAD_GUARD, log = nothing,
                         sample_interval = 5.0, project = Base.active_project(),
                         exeflags = String[]) -> ModeSet

THE m-MODE FAN-OUT of arXiv:2606.04998: the paper's `m = 0…20` sweep (main.tex:716 "Each
m-mode is independent", `CRP_PAPER_MMODES`) solved as ONE INDEPENDENT PROCESS PER MODE, under
the owner's compute policy, with the results collected in `ms` order.

Each mode runs the sealed [`solve_mode_amr`](@ref) on its own [`mode_config`](@ref), in its own
process, with its own solver stack and its own two caches — nothing is shared, which is the
whole reason the fan-out is over PROCESSES and not threads (row 42's rule: the stack and its
caches are neither thread-safe nor reusable across modes). What comes back is a
[`ModeResult`](@ref) per mode: every `AmrLevel` and the final `sol.u`.

# THE COMPUTE POLICY

The binding `compute-usage` skill is owner policy for this 256-core SHARED host, and this
section is the package's SINGLE normative statement of it.

  - **at most `PARALLEL_CPU_BUDGET` = 24 CPUs for our jobs.** `nworkers × threads_per_worker ≤
    24` is enforced for every input by [`mode_pool_size`](@ref); a larger request is CLAMPED
    and `@warn`ed, never silently honoured and never exceeded.
  - **the whole system, all users, must stay below `PARALLEL_LOAD_GUARD` = 90 %.** Before a
    single process is started, [`system_load`](@ref) reads the LIVE usage and the call is
    REFUSED with an `ArgumentError` naming every number if `live + planned > load_guard`.
    During the run a monitor samples the system every `sample_interval` seconds into
    `run.samples` (and into `log`, if given, as a flushed table); `run.max_load` is the
    largest reading and `run.guard_after` is the skill's mandatory after-check.
  - **at most `PARALLEL_MAX_THREADS_PER_WORKER` = 4 Julia threads per process**, clamped the
    same way. The default is 1: an m-mode fits in one core.
  - **the pool is sized `min(24 ÷ threads, nmodes)` and LEFT THERE.** It is computed once,
    before launch, and this function never scales it up.
  - **GPU is banned.** Nothing here touches one.

# BLAS

Every worker calls `BLAS.set_num_threads(blas_threads)` (default 1) around its own solve —
authoritative, not merely inherited from the environment. This is not a micro-optimization:
the subdomain and bottom solves are dense and TINY, and the host's default OpenBLAS pool
spends more time in barriers than in arithmetic, costing a MEASURED factor **8.9** on the
8-level sweep (15.4 s pinned against 136 s unpinned, identical iteration counts). It is also
what makes the fan-out reproducible: results are bitwise identical to
[`solve_modes_serial`](@ref) only at the SAME `blas_threads` and the same Julia thread count
(between BLAS thread counts the diagnostics move by ~4e-12 relative — roundoff, but not
bitwise).

# COST, AND WHY A POOL IS WORTH REUSING

A fresh process must JIT the whole solver stack before it can solve anything: MEASURED **166 s**
for the level-0 solve and **45 s** more for the first AMR step, against **15.9 s** for a warm
8-level sweep. So the paper's 21 modes × 8 levels take ~211 s + 15.9 s on 21 workers against
~211 s + 21 × 15.9 s in one process — a ~2.3× speed-up, JIT-bound, NOT 21×. The speed-up
approaches the worker count only when there are many more modes than workers (a figure driver
sweeping several `(a, r₀)` points) or when the pool is REUSED: pass a [`ModePool`](@ref) via
`pool` and the 211 s is paid once for all sweeps. With `pool` given, `nworkers`,
`threads_per_worker`, `budget`, `project` and `exeflags` are the pool's, not this call's.

# FAILURES

A mode that throws is reported BY ID and never swallowed: it is `@warn`ed with its remote
backtrace, recorded in `set.failures`, printed by `show`, and every OTHER mode still returns.
A caller that requires the complete sweep asserts `isempty(set.failures)`.

    set = solve_modes_parallel(CircularOrbit(0.5, 10.0), crp_paper_config(; spin = 0.5))
    isempty(set.failures) || error("modes failed: ", [f.m for f in set.failures])
    Fr_ingredients = [(set[m].levels[end].psi, set[m].levels[end].dpsi_r) for m in 0:20]
"""
function solve_modes_parallel(orbit::CircularOrbit, config::RunConfig{RT},
                              ms = CRP_PAPER_MMODES;
                              iterations::Integer = config.amr.iterations,
                              T::Type = ComplexF64,
                              criterion = amr_flags,
                              level_options = amr_level_options,
                              pool::Union{Nothing,ModePool} = nothing,
                              nworkers = nothing,
                              threads_per_worker::Integer = 1,
                              blas_threads::Integer = 1,
                              budget::Integer = PARALLEL_CPU_BUDGET,
                              load_guard::Real = PARALLEL_LOAD_GUARD,
                              log = nothing,
                              sample_interval::Real = 5.0,
                              project::AbstractString = _pf_project(),
                              exeflags::AbstractVector{<:AbstractString} = String[]) where {RT}
    blas_threads >= 1 || throw(ArgumentError(
        "solve_modes_parallel: blas_threads must be ≥ 1, got $blas_threads"))
    mlist = _pf_modes(ms)
    if pool === nothing
        p = mode_pool_size(length(mlist); nworkers = nworkers,
                           threads_per_worker = threads_per_worker, budget = budget)
    else
        isopen(pool) || throw(ArgumentError(
            "solve_modes_parallel: the supplied ModePool is closed"))
        p = (nworkers = length(pool), requested = pool.requested,
             threads_per_worker = pool.threads_per_worker, budget = pool.budget,
             clamped = pool.clamped)
    end
    pool === nothing && p.clamped &&
        @warn "solve_modes_parallel: the compute budget clamped the request (compute-usage skill)" requested_workers = p.requested granted_workers = p.nworkers requested_threads = threads_per_worker threads_per_worker = p.threads_per_worker requested_budget = budget budget = p.budget cpu_budget = PARALLEL_CPU_BUDGET
    t0 = time()
    share = p.nworkers * p.threads_per_worker / Sys.CPU_THREADS
    before = _pf_guard(share, load_guard, t0,
                       "solve_modes_parallel($(p.nworkers) × $(p.threads_per_worker) threads)")
    mon = _pf_monitor(t0, sample_interval, log)
    _pf_sample!(mon, before)
    results = ModeResult{T,RT}[]
    failures = ModeFailure[]
    own = pool === nothing
    thepool = pool
    startup = 0.0
    after = before
    try
        _pf_start!(mon)
        if own
            ts = time()
            thepool = mode_pool(p.nworkers; threads_per_worker = p.threads_per_worker,
                                blas_threads = blas_threads, budget = p.budget,
                                load_guard = load_guard, project = project,
                                exeflags = exeflags)
            startup = time() - ts
        end
        raw = Distributed.pmap(Distributed.WorkerPool(copy(thepool.pids)), mlist;
                               on_error = identity) do m
            _pf_solve(orbit, config, m; iterations = iterations, T = T,
                      criterion = criterion, level_options = level_options,
                      blas_threads = blas_threads)
        end
        for (m, r) in zip(mlist, raw)
            if r isa ModeResult
                push!(results, r)
            else
                @warn "solve_modes_parallel: mode m = $m FAILED; the other modes still run" worker = _pf_failed_pid(r) exception = r
                push!(failures, ModeFailure(m, _pf_msg(r), _pf_failed_pid(r)))
            end
        end
    finally
        # the pool is ALWAYS given back and the AFTER reading is ALWAYS taken, throw or not:
        # a leaked worker is 2.5 GB of RSS and a core the owner did not lend.
        _pf_stop!(mon)
        own && thepool !== nothing && close(thepool)
        after = _pf_sample!(mon)
        _pf_close!(mon)
    end
    rec = ModeRunRecord(true, p.nworkers, p.threads_per_worker, Int(blas_threads),
                        p.budget, p.requested, p.clamped, own, Float64(load_guard), share,
                        before, after, maximum(s -> s.cpu_fraction, mon.samples),
                        startup, time() - t0, mon.samples)
    return ModeSet{T,RT}(orbit, config, mlist, results, failures, rec)
end

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
# WAVE-12 FOLLOW-UP (second-pass review). Five things changed
# and every one of them is a CORRECTNESS fix, not a polish pass:
#   * `pmap(…; on_error = identity)` does NOT catch a worker DEATH — MEASURED in
#     `evidence/mmode-parallel-driver/proto_death.jl`, a `ProcessExitedException` escapes pmap
#     and takes the WHOLE sweep down, so one dead worker lost all 21 modes. The fan-out now
#     schedules the modes itself (`_pf_dispatch!`), which is also what makes `timeout_s`
#     possible: `remotecall_fetch` blocks INDEFINITELY and the only way to stop a remote call
#     is to remove the worker.
#   * the CPU budget now counts `nworkers × max(threads_per_worker, blas_threads)`: a CPU spent
#     inside OpenBLAS is a CPU the owner lent us.
#   * the whole-system guard is STRICT (`< load_guard`, the owner's "must stay BELOW 90 %"),
#     takes the MAX of `PARALLEL_LOAD_PROBES` probes, and governs LAUNCH ONLY — `run.max_load`
#     is RECORDED, and a crossing is `@warn`ed after the fact, never silently gated.
#   * `close`/`isopen` are lifecycle-correct: processes are given back FIRST and the pids are
#     cleared only then (a failed teardown leaves them for a retry), and `isopen` asks the
#     CLUSTER whether the workers are alive instead of trusting a flag.
#   * the BLAS pin goes through the binding `NumericsCore.with_blas_threads` lease
#     (docs/src/topics/gpu-parallel.md §"The BLAS-guard rule"), not a hand-rolled
#     get/set/restore that races with any other task in the process.
#
# DEPENDENCY-BUDGET DECISION (mission-approved, wave 12; mirrors adopted decision 3 on
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

"The owner's whole-system CPU guard (`compute-usage` skill: total usage of ALL users < 90 %).
It governs the LAUNCH decision only — see [`solve_modes_parallel`](@ref)."
const PARALLEL_LOAD_GUARD = 0.90

"""How many `/proc/stat` windows the pre-launch guard reads before deciding. The guard takes
the MAX of them: one 0.15 s window can land in a lull between two bursts, and a launch
decision that is wrong is expensive (24 CPUs on a shared host), so the heuristic is
deliberately pessimistic. It is still a HEURISTIC — nothing can promise what the other users
of this box will do one second after we launch."""
const PARALLEL_LOAD_PROBES = 3

# ================================= the live load reading ==================================

"""
    LoadSample

One reading of the whole-system CPU usage, the quantity the owner's 90 % guard is written in
terms of. Fields:

  - `time::Float64` — seconds since the run's origin (`0` is the pre-launch guard sample).
  - `cpu_fraction::Float64` — total CPU usage of ALL users as a fraction of the machine,
    measured as the BUSY fraction of a short `/proc/stat` window (`1 − Δidle/Δtotal`). THE
    guarded quantity, and CLAMPED to `[0, 1]` on every source path — `/proc/stat` exists on
    Linux only, and the two fallbacks below are not bounded by the machine, so the clamp (not
    the estimator) is what makes the `≤ 1` claim true everywhere. A clamped reading reports
    `1.0`, i.e. a full machine, which is the conservative direction for a guard.
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

# The ONE place a `LoadSample` is built, and therefore the one place the `[0, 1]` bound is
# imposed. Pure: takes raw readings, so every source path can be gated in the test suite
# without a loaded machine (`busy === nothing` selects the ps fallback, both `nothing` the
# loadavg one). `ps_fraction` is recorded UNCLAMPED — it is the skill's shell number and its
# unboundedness is exactly the fact this node measured.
function _pf_load_sample(t::Real, ncpu::Integer, la::Real, ps_total, busy)
    nc = Int(ncpu)
    psf = ps_total === nothing ? NaN : Float64(ps_total) / (100 * nc)
    if busy !== nothing
        return LoadSample(Float64(t), clamp(Float64(busy), 0.0, 1.0), psf, Float64(la), nc,
                          :procstat)
    elseif ps_total !== nothing
        return LoadSample(Float64(t), clamp(psf, 0.0, 1.0), psf, Float64(la), nc, :ps)
    else
        return LoadSample(Float64(t), clamp(Float64(la), 0.0, 1.0), NaN, Float64(la), nc,
                          :loadavg)
    end
end

"""
    system_load(t0::Real = time(); window = PARALLEL_LOAD_WINDOW) -> LoadSample

The LIVE whole-system CPU usage, read fresh at every call — the pre-launch guard, the
during-run monitor and the post-run check all go through this one function. `t0` is the origin
the sample's `time` is measured from.

`cpu_fraction` is the BUSY fraction of a `window`-second `/proc/stat` interval,
`1 − Δ(idle + iowait)/Δtotal`, over the aggregate `cpu` line: every user's time, not only ours
— this host is shared and the guard is a whole-system guard — and clamped to `[0, 1]`, which
is what makes it safe to compare against `PARALLEL_LOAD_GUARD`.

WHAT IT IS NOT: an instantaneous truth about a 256-core shared box. It is a `window`-wide
average taken at ONE moment; the pre-launch guard therefore takes the max of
`PARALLEL_LOAD_PROBES` of them, and even that is a LAUNCH HEURISTIC, not a promise about the
next second.

`ps_fraction` records the `compute-usage` skill's own shell formula
`Σ ps -eo pcpu / (100 · ncpu)` beside it. WHY BOTH, AND WHY THE SWITCH: `ps` prints each
process' LIFETIME average, so the sum is not an instantaneous utilisation and is not bounded by
the machine. On a quiet host the two agree (MEASURED 0.0456 against 0.0400); during a compile
storm on this 256-core box the `ps` sum reached **2.868**, i.e. 287 % of the machine, and would
have tripped a 0.90 gate that reality never crossed. The skill's number is kept because it is
the number a human checks at the shell, but the guard is taken on the bounded one.

`source` is `:procstat` normally (Linux), `:ps` if `/proc/stat` could not be read, `:loadavg`
if neither was available — the last two are WEAKER readings, are clamped the same way, and
travel labelled so a run can be audited.
"""
function system_load(t0::Real = time(); window::Real = PARALLEL_LOAD_WINDOW)
    ncpu = Sys.CPU_THREADS
    la = Sys.loadavg()[1] / ncpu
    total = _pf_ps_total()
    a = _pf_procstat()
    if a !== nothing
        sleep(window)
        b = _pf_procstat()
        if b !== nothing && b[2] > a[2]
            return _pf_load_sample(time() - Float64(t0), ncpu, la, total,
                                   (b[1] - a[1]) / (b[2] - a[2]))
        end
    end
    return _pf_load_sample(time() - Float64(t0), ncpu, la, total, nothing)
end

# The owner's sentence is "total CPU usage of ALL users must stay BELOW 90 %", so the
# comparison is STRICT: a plan that lands exactly on the guard is refused, not admitted.
_pf_admits(live::Real, planned::Real, guard::Real) = (Float64(live) + Float64(planned)) < Float64(guard)

# The pre-launch refusal. `planned_share` is the fraction of the machine THIS launch is about
# to ADD (`nworkers × cpus_per_worker / ncpu`); it is ZERO for a sweep that reuses an
# already-running pool, whose processes are already inside the live reading and would
# otherwise be counted twice. Takes the MAX of `PARALLEL_LOAD_PROBES` windows and refuses with
# every number in the message.
function _pf_guard(planned_share::Real, load_guard::Real, t0::Real, what::AbstractString;
                   probes::Integer = PARALLEL_LOAD_PROBES)
    (load_guard isa Real && 0 < load_guard <= PARALLEL_LOAD_GUARD) || throw(ArgumentError(
        "$what: load_guard must lie in (0, $PARALLEL_LOAD_GUARD] — the owner's whole-system " *
        "guard is $(PARALLEL_LOAD_GUARD) (compute-usage skill), got $load_guard"))
    s = system_load(t0)
    for _ in 2:max(1, Int(probes))
        t = system_load(t0)
        t.cpu_fraction > s.cpu_fraction && (s = t)
    end
    total = s.cpu_fraction + planned_share
    _pf_admits(s.cpu_fraction, planned_share, load_guard) || throw(ArgumentError(
        "$what: REFUSING TO LAUNCH — the whole-system CPU guard would be exceeded. " *
        "live usage $(round(s.cpu_fraction; digits = 4)) (worst of $(max(1, Int(probes))) " *
        "probes, source :$(s.source), ps sum $(round(s.ps_fraction; digits = 4)), " *
        "loadavg1/ncpu $(round(s.loadavg1; digits = 4)), ncpu $(s.ncpu)) " *
        "+ planned share $(round(planned_share; digits = 4)) = $(round(total; digits = 4)) " *
        "≥ load_guard $(load_guard). Wait for the host to quieten or shrink the pool."))
    return s
end

# ================================= the pool arithmetic =====================================

"""
    mode_pool_size(nmodes; nworkers = nothing, threads_per_worker = 1, blas_threads = 1,
                   budget = PARALLEL_CPU_BUDGET)
        -> (; nworkers, requested, threads_per_worker, blas_threads, cpus_per_worker,
              budget, clamped)

The owner's pool rule, as arithmetic: **size the pool as `min(budget ÷ cpus_per_worker,
nmodes)` and leave it there; never auto-scale up.**

`cpus_per_worker = max(threads_per_worker, blas_threads)` is what one worker actually occupies.
Counting the BLAS threads is not pedantry: OpenBLAS runs inside the worker's own process, so a
24-worker pool at `blas_threads = 2` is 48 CPUs on the owner's 24-CPU budget, and the landed
first version of this function could not see that. The two counts do not ADD — BLAS work runs
on the worker's threads — so the cost of a worker is their MAX.

`nworkers === nothing` (the default) asks for the rule. An explicit `nworkers` is a REQUEST:
it is clamped to `min(budget ÷ cpus_per_worker, nmodes)` — a pool bigger than the budget would
break owner policy, and one bigger than the mode count would only start idle processes.
`threads_per_worker` and `blas_threads` are clamped to `PARALLEL_MAX_THREADS_PER_WORKER` and
`budget` to `PARALLEL_CPU_BUDGET`; the returned `clamped` flag is `true` whenever ANY of them
was reduced, and the callers log it — a clamp is never silent.

The post-condition `nworkers × cpus_per_worker ≤ PARALLEL_CPU_BUDGET` holds for every input,
including hostile ones.

    julia> mode_pool_size(21)                       # the paper's sweep
    (nworkers = 21, requested = 21, threads_per_worker = 1, blas_threads = 1,
     cpus_per_worker = 1, budget = 24, clamped = false)

    julia> mode_pool_size(100; nworkers = 64)       # an over-request, clamped and reported
    (nworkers = 24, requested = 64, …, clamped = true)

    julia> mode_pool_size(100; blas_threads = 2).nworkers      # 2 BLAS threads each: 12 × 2 = 24
    12
"""
function mode_pool_size(nmodes::Integer; nworkers = nothing, threads_per_worker::Integer = 1,
                        blas_threads::Integer = 1, budget::Integer = PARALLEL_CPU_BUDGET)
    nmodes >= 1 || throw(ArgumentError("mode_pool_size: nmodes must be ≥ 1, got $nmodes"))
    threads_per_worker >= 1 || throw(ArgumentError(
        "mode_pool_size: threads_per_worker must be ≥ 1, got $threads_per_worker"))
    blas_threads >= 1 || throw(ArgumentError(
        "mode_pool_size: blas_threads must be ≥ 1, got $blas_threads"))
    budget >= 1 || throw(ArgumentError("mode_pool_size: budget must be ≥ 1, got $budget"))
    nth = min(Int(threads_per_worker), PARALLEL_MAX_THREADS_PER_WORKER)
    nbl = min(Int(blas_threads), PARALLEL_MAX_THREADS_PER_WORKER)
    cpw = max(nth, nbl)
    bud = min(Int(budget), PARALLEL_CPU_BUDGET)
    cap = max(1, bud ÷ cpw)
    req = nworkers === nothing ? min(cap, Int(nmodes)) : Int(nworkers)
    req >= 1 || throw(ArgumentError("mode_pool_size: nworkers must be ≥ 1, got $req"))
    n = clamp(req, 1, min(cap, Int(nmodes)))
    clamped = n != req || nth != Int(threads_per_worker) || nbl != Int(blas_threads) ||
              bud != Int(budget)
    return (nworkers = n, requested = req, threads_per_worker = nth, blas_threads = nbl,
            cpus_per_worker = cpw, budget = bud, clamped = clamped)
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

# ================================== the per-mode stamp =====================================

"""
    ModeStamp{RT}

The m-DEPENDENT solver knobs one worker actually used, reported back with its result so a
caller can PROVE that the process fan-out and the in-process oracle solved the same problem —
a bitwise comparison of two answers to two different questions proves nothing.

  - `m::Int`, `spin::RT`, `orbital_radius::RT` — the problem.
  - `subdomain_complex_shift::RT`, `multigrid_complex_shift::RT` — the shifts the config
    carried into the Schwarz subdomain solves and the multigrid levels (yaml:203, 194).
  - `helmholtz_shift_policy::RT` — `helmholtz_shift(m, Ω; slicing = :vtu)`, the CRP policy
    value under the paper's piecewise-null slicing. It is EXACTLY ZERO (main.tex:680: under
    hyperboloidal slicing the shift "can be reduced or disabled altogether"), and the shipped
    vtu config sets `ComplexShift: 0.0`.
  - `outer_bc_m::Int`, `outer_bc_order::Int` — the outer radial condition's own mode number and
    order; `outer_bc_m = -1` records the `m = 0` yaml branch, which HAS no outer condition.
  - `solve_kwargs::String` — the extra [`solve_mode_amr`](@ref) keywords this mode was solved
    with, canonically rendered (keys sorted) — `""` when the sweep used the driver's own
    arguments and nothing else. A sweep run with `horizon_hold = true` (the row-58 cure) and
    one run without it are DIFFERENT problems, and this is where an archived artefact says
    which one it holds instead of leaving a reader to guess.

WHY THE SHIFT IS RECORDED PER MODE. The preconditioner this sweep runs is m-AGNOSTIC — shift
zero — while the operator's Helmholtz term scales as `(mΩ)²`. That is the paper's own
configuration, and it is the reason the measured GMRES count GROWS with `m` (see
[`solve_modes_parallel`](@ref)). Turning the shift on makes it far worse on this formulation
(MEASURED in the row-47 review at `s = (mΩ)²`: m = 10 goes 13 → 92 iterations, m = 20 goes
16 → 215), which is the round-1 calibration falsifier of [`helmholtz_shift`](@ref) showing up
again. Recording the value per mode keeps that claim auditable instead of remembered.
"""
struct ModeStamp{RT<:AbstractFloat}
    m::Int
    spin::RT
    orbital_radius::RT
    subdomain_complex_shift::RT
    multigrid_complex_shift::RT
    helmholtz_shift_policy::RT
    outer_bc_m::Int
    outer_bc_order::Int
    solve_kwargs::String
end

# Canonical rendering of a pass-through: keys SORTED, so two callers that wrote the same
# keywords in a different order produce the same stamp and compare equal.
function _pf_render_kwargs(nt::NamedTuple)
    isempty(nt) && return ""
    ks = sort(collect(keys(nt)))
    return string(NamedTuple{Tuple(ks)}(Tuple(getfield(nt, k) for k in ks)))
end

Base.:(==)(a::ModeStamp, b::ModeStamp) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(ModeStamp))

"""
    mode_stamp(orbit::CircularOrbit, config::RunConfig, m::Integer; solve_kwargs = (;))
        -> ModeStamp

The [`ModeStamp`](@ref) of mode `m` under `config`, computed from [`mode_config`](@ref) — the
same function the worker runs, so the master can check a returned stamp against what it would
itself have used.
"""
function mode_stamp(orbit::CircularOrbit, config::RunConfig{RT}, m::Integer;
                    solve_kwargs::NamedTuple = (;)) where {RT}
    cfg = mode_config(config, m)
    bc = cfg.domain.outer_radial_bc
    return ModeStamp{RT}(Int(m), RT(orbit.a), RT(orbit.r0),
                         RT(cfg.schwarz.subdomain_complex_shift),
                         RT(cfg.multigrid.complex_shift),
                         RT(helmholtz_shift(Int(m), orbit.Omega; slicing = :vtu)),
                         bc === nothing ? -1 : Int(bc.m_mode_number),
                         bc === nothing ? -1 : Int(bc.order),
                         _pf_render_kwargs(solve_kwargs))
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
  - `converged::Bool` — whether EVERY AMR level's linear solve met its tolerance. It is the
    roll-up of `levels[i].record.converged`, which stays available level by level. A stagnated
    GMRES still returns a `ModeResult` — it is a number, just not a converged one — and a
    caller that plots it without looking has published a wrong value, so the flag travels
    beside the answer rather than inside it.
  - `stamp::ModeStamp{RT}` — the m-dependent knobs the producing process used, so `parallel ==
    serial` can be asserted on the QUESTION as well as on the answer (see [`ModeStamp`](@ref)).

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
    converged::Bool
    stamp::ModeStamp{RT}
end

ModeResult(h::AmrHistory{T,RT}, wall_time::Real, stamp::ModeStamp{RT}) where {T,RT} =
    ModeResult{T,RT}(h.m, h.x_p, h.levels, h.solution.u, Float64(wall_time),
                     Distributed.myid(), getpid(), BLAS.get_num_threads(), Threads.nthreads(),
                     all(l -> l.record.converged, h.levels), stamp)

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

  - `parallel::Bool`, `nworkers::Int`, `threads_per_worker::Int`, `blas_threads::Int`,
    `cpus_per_worker::Int` — the shape of the run; `cpus_per_worker` is
    `max(threads_per_worker, blas_threads)`, the CPUs one worker actually occupies, and
    `nworkers × cpus_per_worker` is what the owner's 24-CPU budget is spent on;
  - `budget::Int`, `requested_workers::Int`, `clamped::Bool` — the [`mode_pool_size`](@ref)
    arithmetic as applied, with `clamped` recording that a request was reduced (and logged);
  - `owned_pool::Bool` — whether the sweep started and stopped its own [`ModePool`](@ref);
  - `load_guard::Float64` — the guard the LAUNCH decision was taken against;
  - `planned_share::Float64` — the fraction of the machine this run OCCUPIES
    (`nworkers × cpus_per_worker / ncpu`);
  - `guard_share::Float64` — the fraction the pre-launch guard ADDED to the live reading. It
    equals `planned_share` when this call creates the capacity, and is `0` when an existing
    pool was reused: those processes are already inside the live reading and counting them
    twice would refuse launches that are perfectly legal;
  - `guard_before::LoadSample`, `guard_after::LoadSample` — the `compute-usage` skill's
    BEFORE and AFTER readings;
  - `max_load::Float64` — the largest `cpu_fraction` over every sample of the run. It is
    RECORDED, NOT ENFORCED: the guard governs the launch decision, and nothing here can stop
    another user of a shared box afterwards. A crossing is `@warn`ed when the run ends;
  - `solve_kwargs::NamedTuple` — the extra [`solve_mode_amr`](@ref) keywords forwarded to
    every mode of this sweep (`(;)` by default);
  - `timeout_s::Union{Nothing,Float64}`, `timed_out::Bool` — the deadline, if any, and whether
    it fired (in which case the workers were REMOVED to cancel the outstanding modes);
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
    cpus_per_worker::Int
    budget::Int
    requested_workers::Int
    clamped::Bool
    owned_pool::Bool
    load_guard::Float64
    planned_share::Float64
    guard_share::Float64
    guard_before::LoadSample
    guard_after::LoadSample
    max_load::Float64
    solve_kwargs::NamedTuple
    timeout_s::Union{Nothing,Float64}
    timed_out::Bool
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
    @printf(io, "  %d worker(s) × %d thread(s), BLAS %d ⇒ %d CPU(s) each, budget %d (requested %d%s) — %.1f s wall, %.1f s startup%s\n",
            r.nworkers, r.threads_per_worker, r.blas_threads, r.cpus_per_worker, r.budget,
            r.requested_workers, r.clamped ? ", CLAMPED" : "", r.wall_time, r.startup_time,
            r.timed_out ? ", TIMED OUT" : "")
    @printf(io, "  system load: before %.4f, max %.4f, after %.4f (guard %.2f, %d samples)\n",
            r.guard_before.cpu_fraction, r.max_load, r.guard_after.cpu_fraction,
            r.load_guard, length(r.samples))
    println(io, "    m  levels   ndofs  iters   wall[s]  worker      Re Ψ^R(x_p)     ∂_r Ψ^R(x_p)")
    for x in s.results
        l = x.levels[end]
        @printf(io, "  %3d  %6d  %6d  %5d  %8.2f  %6d  %+13.6e  %+13.6e%s\n",
                x.m, length(x.levels), l.mesh.ndofs, l.record.iterations, x.wall_time,
                x.ospid, real(l.psi), real(l.dpsi_r),
                x.converged ? "" : "   NOT CONVERGED")
    end
    isempty(r.solve_kwargs) ||
        println(io, "  solve_kwargs: ", _pf_render_kwargs(r.solve_kwargs))
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
        println(io, "# the guard governs the LAUNCH decision only (max of $PARALLEL_LOAD_PROBES ",
                "probes, taken before the first process starts); the rows below are a RECORD of ",
                "what the shared machine then did, and are not enforced — a crossing is warned, ",
                "never silently gated. cpu_fraction = /proc/stat busy fraction, clamped to [0,1]; ",
                "ps_fraction = the skill's Σ ps -eo pcpu / (100·ncpu), recorded UNCLAMPED.")
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
`threads_per_worker`, the `blas_threads` it was BUDGETED for and the resulting
`cpus_per_worker`, the budget arithmetic that sized it (`requested`, `budget`, `clamped`), the
`project` the workers were started in, and `startup_time`.

WHY IT IS A FIRST-CLASS OBJECT. A fresh Julia process pays the WHOLE solver stack's JIT before
it can solve anything: MEASURED on this host, **166 s** for the level-0 solve and a further
**45 s** for the first AMR step, against **15.9 s** for a warm 8-level sweep. A sweep that
starts and stops its own pool therefore pays ~211 s of compilation no matter how small it is,
while a REUSED pool pays it once and then runs each further sweep at full speed. Callers that
run several sweeps — the figure drivers of rows 50/51/53 — must hold one pool across all of
them. Build it with [`mode_pool`](@ref) and always `close` it (the `do`-block form does that
for you).

LIVENESS, NOT A FLAG. `isopen(pool)` asks the CLUSTER whether the workers are still there
(`Distributed.procs()`), because a worker can die on its own — a segfault, an OOM kill, an
`exit` inside a solve — and a pool that believed a stale flag would hand a corpse out as if it
were a worker. `length(pool)` is the number of pids it currently holds; a sweep prunes the
dead ones (with a warning) before and after it runs.
"""
mutable struct ModePool
    pids::Vector{Int}
    threads_per_worker::Int
    blas_threads::Int
    cpus_per_worker::Int
    budget::Int
    requested::Int
    clamped::Bool
    project::String
    startup_time::Float64
    open::Bool
end

# the cluster's own answer, which is the only authority on whether a worker still exists
_pf_live(pids) = filter(w -> w in Distributed.procs(), pids)

Base.isopen(p::ModePool) = p.open && !isempty(_pf_live(p.pids))
Base.length(p::ModePool) = length(p.pids)

# drop the corpses and say so. Returns the pids that were removed.
function _pf_prune!(p::ModePool, what::AbstractString)
    dead = filter(w -> !(w in Distributed.procs()), p.pids)
    if !isempty(dead)
        filter!(w -> w in Distributed.procs(), p.pids)
        @warn "$what: worker process(es) are gone and were pruned from the pool" dead_workers = dead remaining = copy(p.pids)
    end
    return dead
end

function Base.show(io::IO, p::ModePool)
    print(io, "ModePool(", length(p.pids), " worker", length(p.pids) == 1 ? "" : "s",
          " × ", p.threads_per_worker, " thread", p.threads_per_worker == 1 ? "" : "s",
          ", BLAS ≤ ", p.blas_threads, isopen(p) ? "" : ", CLOSED", ")")
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

The request goes through [`mode_pool_size`](@ref), so `nworkers × max(threads_per_worker,
blas_threads)` never exceeds `PARALLEL_CPU_BUDGET`; a reduced request is `@warn`ed, never
silently applied. Before a single process is started the LIVE whole-system usage is read (the
worst of `PARALLEL_LOAD_PROBES` probes) and the pool is REFUSED with an `ArgumentError` naming
every number if `live + planned ≥ load_guard`.

`blas_threads` is the pool's BUDGETED BLAS ceiling: a later sweep may pin fewer threads, but
one that asks for more than the pool was sized for is refused rather than quietly overspending
the owner's CPU budget.

`close(pool)` removes the workers; the `do`-block form does it even on a throw:

    mode_pool(21) do pool
        conv = solve_modes_parallel(orbit, config, 0:20; pool = pool)
        modes = solve_modes_parallel(orbit, config, 0:20; pool = pool, iterations = 4)
    end

BLAS pinning is applied INSIDE each worker at every solve (see [`solve_modes_parallel`](@ref))
through the binding `NumericsCore.with_blas_threads` lease, so one pool serves sweeps with
different `blas_threads`; the environment variable set here only sizes OpenBLAS's pool at
process start.
"""
function mode_pool(nworkers::Integer; threads_per_worker::Integer = 1,
                   blas_threads::Integer = 1, budget::Integer = PARALLEL_CPU_BUDGET,
                   load_guard::Real = PARALLEL_LOAD_GUARD,
                   project::AbstractString = _pf_project(),
                   exeflags::AbstractVector{<:AbstractString} = String[])
    blas_threads >= 1 || throw(ArgumentError("mode_pool: blas_threads must be ≥ 1, got $blas_threads"))
    p = mode_pool_size(typemax(Int); nworkers = nworkers,
                       threads_per_worker = threads_per_worker,
                       blas_threads = blas_threads, budget = budget)
    p.clamped && @warn "mode_pool: the compute budget clamped the request (compute-usage skill)" requested_workers = p.requested granted_workers = p.nworkers requested_threads = threads_per_worker threads_per_worker = p.threads_per_worker requested_blas = blas_threads blas_threads = p.blas_threads cpus_per_worker = p.cpus_per_worker requested_budget = budget budget = p.budget cpu_budget = PARALLEL_CPU_BUDGET
    share = p.nworkers * p.cpus_per_worker / Sys.CPU_THREADS
    _pf_guard(share, load_guard, time(),
              "mode_pool($(p.nworkers) × $(p.cpus_per_worker) CPU(s))")
    t = time()
    pids = Distributed.addprocs(p.nworkers;
        exeflags = String[String("--project=" * project),
                          String("-t" * string(p.threads_per_worker)), exeflags...],
        env = ["OPENBLAS_NUM_THREADS" => string(p.blas_threads),
               "JULIA_NUM_THREADS" => string(p.threads_per_worker),
               "JULIA_PKG_OFFLINE" => "true"],
        topology = :master_worker)
    try
        Distributed.remotecall_eval(Main, pids, :(using ScalarSelfForce))
    catch e
        Distributed.rmprocs(pids; waitfor = 60)
        rethrow(e)
    end
    return ModePool(collect(Int, pids), p.threads_per_worker, p.blas_threads,
                    p.cpus_per_worker, p.budget, p.requested, p.clamped, String(project),
                    time() - t, true)
end

"""
    close(pool::ModePool) -> pool

Give the worker processes back. The processes are removed FIRST and the bookkeeping is cleared
only once they are gone, so a teardown that throws leaves the pids in place for a retry
instead of forgetting 2.5 GB of RSS per worker; `close` is therefore idempotent AND retryable,
and closing a pool whose workers already died is clean, not an error.
"""
function Base.close(p::ModePool)
    (p.open || !isempty(p.pids)) || return p
    try
        if !isempty(p.pids)
            alive = _pf_live(p.pids)
            isempty(alive) || Distributed.rmprocs(alive; waitfor = 60)
            empty!(p.pids)          # reached only when the removal returned
        end
    finally
        p.open = false
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

# ============================ the solve_mode_amr pass-through ==============================

"""The keywords the driver itself owns and passes explicitly. A caller that repeated one of
these in `solve_kwargs` would produce a duplicate-keyword call inside the worker, so they are
refused on the master instead."""
const PARALLEL_DRIVER_KWARGS = (:iterations, :T, :criterion, :level_options)

# every keyword `solve_mode_amr` actually accepts, read off its own method table so this can
# never drift from the sealed signature
function _pf_amr_kwargs()
    ks = Set{Symbol}()
    for mth in methods(solve_mode_amr), k in Base.kwarg_decl(mth)
        endswith(String(k), "...") || push!(ks, k)
    end
    return ks
end

# The caller-side gate. It runs BEFORE the load guard and before a single worker is started:
# a bad keyword that only surfaced inside the pool would cost a pool start and come back as
# `nworkers` identical remote MethodErrors with no useful attribution.
function _pf_check_solve_kwargs(nt::NamedTuple)
    isempty(nt) && return nothing
    owned = [k for k in keys(nt) if k in PARALLEL_DRIVER_KWARGS]
    isempty(owned) || throw(ArgumentError(
        "solve_modes_*: solve_kwargs may not carry $(join(owned, ", ")) — the driver passes " *
        "$(join(PARALLEL_DRIVER_KWARGS, ", ")) itself, and repeating one would be a duplicate " *
        "keyword inside the worker. Pass it as the driver's own argument instead."))
    allowed = _pf_amr_kwargs()
    bad = [k for k in keys(nt) if !(k in allowed)]
    isempty(bad) || throw(ArgumentError(
        "solve_modes_*: solve_kwargs carries $(join(bad, ", ")), which solve_mode_amr does " *
        "not accept. It accepts: $(join(sort(collect(allowed)), ", ")). " *
        "(Checked on the master, before any worker process is started.)"))
    return nothing
end

# THE function every worker runs. Top-level (not a closure) so it serializes by name, and the
# ONE place the BLAS pin is applied — the environment is best effort, this is authoritative.
# The pin goes through `NumericsCore.with_blas_threads`, the package-wide lease of
# docs/src/topics/gpu-parallel.md §"The BLAS-guard rule": BLAS' thread count is process-global,
# so a hand-rolled get/set/restore races with every other task in the same process (the load
# monitor's sampling task, a caller that threads over sweeps) and can restore the wrong value.
function _pf_solve(orbit::CircularOrbit, config::RunConfig, m::Integer;
                   iterations::Integer, T::Type, criterion, level_options,
                   blas_threads::Integer, solve_kwargs::NamedTuple = (;))
    return with_blas_threads(blas_threads) do
        cfg = mode_config(config, m)
        t = time()
        h = solve_mode_amr(orbit, cfg, m; iterations = iterations, T = T,
                           criterion = criterion, level_options = level_options,
                           solve_kwargs...)
        return ModeResult(h, time() - t,
                          mode_stamp(orbit, config, m; solve_kwargs = solve_kwargs))
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

# Who to blame for a failure. A remote THROW carries its worker id in the RemoteException; a
# remote DEATH carries it in `ProcessExitedException.worker_id` — unwrapping that is what turns
# "worker = 0" into an actionable pid (review finding 4).
function _pf_failed_pid(e)
    e isa Distributed.RemoteException && return e.pid
    e isa Distributed.ProcessExitedException && return e.worker_id
    if e isa CompositeException && !isempty(e.exceptions)
        return _pf_failed_pid(first(e.exceptions))
    end
    e isa TaskFailedException && return _pf_failed_pid(e.task.result)
    return 0
end

# ================================== the mode scheduler =====================================

# WHY THIS EXISTS INSTEAD OF `pmap`. `pmap(f, pool, items; on_error = identity)` turns a THROW
# inside `f` into a value, but it does NOT catch the death of the worker running it: MEASURED
# in `evidence/mmode-parallel-driver/proto_death.jl`, a worker that exits mid-element raises
# `ProcessExitedException` OUT of pmap and the whole sweep dies with it — one dead process
# would have lost all 21 modes of the paper's sweep. `pmap` also has no deadline, while
# `remotecall_fetch` blocks forever by construction.
#
# So the fan-out schedules its own modes: one task per mode, workers handed out through a
# channel, EVERY exception (throw, death, cancellation) caught at the mode that owns it, a
# corpse never handed out twice, and the queue unblocked when the last worker dies. On a
# deadline the outstanding modes are cancelled the only way a remote call can be cancelled —
# by removing the workers — and each is reported by its own id.
function _pf_dispatch!(call, out::Vector{Any}, pids::Vector{Int}, mlist::Vector{Int};
                       timeout_s = nothing, what::AbstractString = "solve_modes_parallel")
    avail = Channel{Int}(max(1, length(pids)))
    foreach(w -> put!(avail, w), pids)
    nlive = Threads.Atomic{Int}(length(pids))
    cancelled = Ref(false)
    body = @async begin
        @sync for (i, m) in enumerate(mlist)
            @async begin
                w = 0
                try
                    w = take!(avail)
                    out[i] = call(w, m)
                catch e
                    out[i] = e
                end
                if w != 0
                    if !cancelled[] && (w in Distributed.procs())
                        try
                            put!(avail, w)
                        catch
                        end
                    elseif Threads.atomic_sub!(nlive, 1) <= 1
                        try                      # the last worker is gone: unblock the queue
                            close(avail)
                        catch
                        end
                    end
                end
            end
        end
    end
    timed_out = false
    if timeout_s === nothing
        try
            wait(body)
        catch e
            @warn "$what: the mode scheduler itself threw" exception = e
        end
    else
        deadline = time() + Float64(timeout_s)
        while !istaskdone(body) && time() < deadline
            sleep(0.05)
        end
        if !istaskdone(body)
            timed_out = true
            cancelled[] = true
            live = _pf_live(pids)
            @warn "$what: TIMED OUT — cancelling the outstanding modes by removing their workers (a remote call cannot be interrupted any other way)" timeout_s = Float64(timeout_s) workers = live
            try
                close(avail)
            catch
            end
            isempty(live) || Distributed.rmprocs(live; waitfor = 60)
            try
                wait(body)
            catch e
                @warn "$what: the mode scheduler threw while cancelling" exception = e
            end
        end
    end
    for i in eachindex(out)
        isassigned(out, i) || (out[i] = ErrorException(
            "$what: mode m = $(mlist[i]) never ran — the sweep timed out after " *
            "$(timeout_s) s and its worker was removed"))
    end
    return timed_out
end

# ================================ the shared run scaffold ==================================
# The two drivers differ in ONE thing: where a mode is solved. Everything around it — the
# guard, the monitor, the collection, the teardown discipline and the record — is the same
# code, and lives here so the oracle cannot drift from the thing it is the oracle of.

# guard + monitor, in that order, before anything is launched
function _pf_open_run(what, guard_share, load_guard, t0, log, sample_interval)
    before = _pf_guard(guard_share, load_guard, t0, what)
    mon = _pf_monitor(t0, sample_interval, log)
    _pf_sample!(mon, before)
    return before, mon
end

# one raw entry per mode, in `mlist` order, into results/failures — a failure is never swallowed
function _pf_collect!(results, failures, mlist, raw, what)
    for (m, r) in zip(mlist, raw)
        if r isa ModeResult
            push!(results, r)
        else
            @warn "$what: mode m = $m FAILED; the other modes still run" worker = _pf_failed_pid(r) exception = r
            push!(failures, ModeFailure(m, _pf_msg(r), _pf_failed_pid(r)))
        end
    end
    return results, failures
end

# The teardown. EVERY step runs in its own try/catch so that a secondary failure cannot mask
# the primary exception and cannot skip the step after it: the `compute-usage` skill's AFTER
# reading is mandatory, and the landed version silently skipped it whenever the pool teardown
# threw.
function _pf_close_run!(mon, fallback::LoadSample, teardown, what)
    try
        _pf_stop!(mon)
    catch e
        @warn "$what: the load monitor did not stop cleanly" exception = e
    end
    try
        teardown()
    catch e
        @warn "$what: giving the worker processes back FAILED — they may be leaked; remove them by pid" exception = e
    end
    after = fallback
    try
        after = _pf_sample!(mon)          # the skill's mandatory AFTER reading, throw or not
    catch e
        @warn "$what: the mandatory AFTER load reading failed" exception = e
    end
    try
        _pf_close!(mon)
    catch e
        @warn "$what: the load-monitor log did not close cleanly" exception = e
    end
    return after
end

# `max_load` is RECORDED, not enforced — but it is never silently ignored either.
function _pf_max_load(mon, load_guard, what)
    ml = isempty(mon.samples) ? 0.0 : maximum(s -> s.cpu_fraction, mon.samples)
    ml < load_guard || @warn "$what: the whole-system load CROSSED the guard while the run was in flight. The guard governs the LAUNCH decision; this is a record of what the shared machine then did, and it is reported, not enforced." max_load = ml load_guard = load_guard samples = length(mon.samples)
    return ml
end

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
                            solve_kwargs::NamedTuple = (;),
                            log = nothing,
                            sample_interval::Real = 5.0) where {RT}
    blas_threads >= 1 || throw(ArgumentError(
        "solve_modes_serial: blas_threads must be ≥ 1, got $blas_threads"))
    _pf_check_solve_kwargs(solve_kwargs)
    mlist = _pf_modes(ms)
    what = "solve_modes_serial(1 process)"
    t0 = time()
    # one process, solving one mode at a time: it occupies its own core, plus however many
    # OpenBLAS may use inside it — the same `max` the pool arithmetic applies to a worker.
    cpw = max(1, Int(blas_threads))
    share = cpw / Sys.CPU_THREADS
    before, mon = _pf_open_run(what, share, load_guard, t0, log, sample_interval)
    results = ModeResult{T,RT}[]
    failures = ModeFailure[]
    after = before
    try
        with_blas_threads(blas_threads) do
            for m in mlist
                _pf_sample!(mon)
                try
                    push!(results, _pf_solve(orbit, config, m; iterations = iterations, T = T,
                                             criterion = criterion,
                                             level_options = level_options,
                                             blas_threads = blas_threads,
                                             solve_kwargs = solve_kwargs))
                catch e
                    @warn "$what: mode m = $m FAILED; the other modes still run" exception = (e, catch_backtrace())
                    push!(failures, ModeFailure(m, _pf_msg(e), Distributed.myid()))
                end
            end
        end
    finally
        after = _pf_close_run!(mon, before, () -> nothing, what)
    end
    rec = ModeRunRecord(false, 1, Threads.nthreads(), Int(blas_threads), cpw,
                        PARALLEL_CPU_BUDGET, 1, false, false, Float64(load_guard), share,
                        share, before, after, _pf_max_load(mon, load_guard, what),
                        solve_kwargs, nothing, false, 0.0, time() - t0, mon.samples)
    return ModeSet{T,RT}(orbit, config, mlist, results, failures, rec)
end

"""
    solve_modes_parallel(orbit, config, ms = CRP_PAPER_MMODES;
                         iterations = config.amr.iterations, T = ComplexF64,
                         criterion = amr_flags, level_options = amr_level_options,
                         pool = nothing, nworkers = nothing, threads_per_worker = 1,
                         blas_threads = 1, budget = PARALLEL_CPU_BUDGET,
                         load_guard = PARALLEL_LOAD_GUARD, timeout_s = nothing,
                         log = nothing, sample_interval = 5.0,
                         project = Base.active_project(), exeflags = String[]) -> ModeSet

THE m-MODE FAN-OUT of arXiv:2606.04998: the paper's `m = 0…20` sweep (main.tex:716 "Each
m-mode is independent", `CRP_PAPER_MMODES`) solved as ONE INDEPENDENT PROCESS PER MODE, under
the owner's compute policy, with the results collected in `ms` order.

Each mode runs the sealed [`solve_mode_amr`](@ref) on its own [`mode_config`](@ref), in its own
process, with its own solver stack and its own two caches — nothing is shared, which is the
whole reason the fan-out is over PROCESSES and not threads (row 42's rule: the stack and its
caches are neither thread-safe nor reusable across modes). What comes back is a
[`ModeResult`](@ref) per mode: every `AmrLevel`, the final `sol.u`, and the
[`ModeStamp`](@ref) of the knobs that mode actually used.

# THE COMPUTE POLICY

The binding `compute-usage` skill is owner policy for this 256-core SHARED host, and this
section is the package's SINGLE normative statement of it.

  - **at most `PARALLEL_CPU_BUDGET` = 24 CPUs for our jobs.** `nworkers × max(threads_per_worker,
    blas_threads) ≤ 24` is enforced for every input by [`mode_pool_size`](@ref); a larger
    request is CLAMPED and `@warn`ed, never silently honoured and never exceeded. A sweep on a
    SUPPLIED pool that asks for more BLAS threads than the pool was budgeted for is REFUSED.
  - **the whole system, all users, must stay below `PARALLEL_LOAD_GUARD` = 90 %.** Before a
    single process is started, [`system_load`](@ref) reads the LIVE usage — the worst of
    `PARALLEL_LOAD_PROBES` probes — and the call is REFUSED with an `ArgumentError` naming
    every number if `live + planned ≥ load_guard`. **The guard governs the LAUNCH decision
    only.** During the run a monitor samples the system every `sample_interval` seconds into
    `run.samples` (and into `log`, if given, as a flushed table); `run.max_load` is the largest
    reading and is a RECORD, not a gate — nothing here can stop another user of a shared box
    once we have launched, so a crossing is `@warn`ed at the end of the run instead of
    pretending to have prevented it. `run.guard_after` is the skill's mandatory after-check,
    and it is taken even when the sweep threw.
  - **at most `PARALLEL_MAX_THREADS_PER_WORKER` = 4 Julia threads per process**, clamped the
    same way. The default is 1: an m-mode fits in one core.
  - **the pool is sized `min(24 ÷ cpus_per_worker, nmodes)` and LEFT THERE.** It is computed
    once, before launch, and this function never scales it up.
  - **GPU is banned.** Nothing here touches one.

# BLAS

Every worker pins OpenBLAS to `blas_threads` (default 1) around its own solve through the
binding `NumericsCore.with_blas_threads` lease — authoritative, not merely inherited from the
environment. This is not a micro-optimization: the subdomain and bottom solves are dense and
TINY, and the host's default OpenBLAS pool spends more time in barriers than in arithmetic,
costing a MEASURED factor **8.9** on the 8-level sweep (15.4 s pinned against 136 s unpinned,
identical iteration counts). It is also what makes the fan-out reproducible: results are
bitwise identical to [`solve_modes_serial`](@ref) only at the SAME `blas_threads` and the same
Julia thread count (between BLAS thread counts the diagnostics move by ~4e-12 relative —
roundoff, but not bitwise).

# COST — THREE REGIMES, AND WHY A POOL IS WORTH REUSING

MEASURED on the current tree for the paper's whole `m = 0…20` sweep at `a = 0.5`, `r₀ = 10M`
through 9 grids (`evidence/mmode-parallel-driver/artefact_21modes_v3.json`; ONE run on a
256-core SHARED host, so treat the wall times as PRELIMINARY):

| regime | wall |
|---|---|
| COLD — a fresh 21-worker pool, half domain | 818.2 s (16.6 s pool start + 801.5 s) |
| SERIAL — one process, half domain | 864.7 s |
| **REUSED POOL** — the FULL-z sweep on those same workers | **282.6 s**, for meshes ~1.7× larger |

A fresh worker pays the whole solver stack's JIT before its first solve, so 21 of them pay it 21
times over and the COLD fan-out only just matches a single warm process (1.06× here; an earlier
tree measured 0.75×, an outright loss). The win is REUSE: hold ONE [`ModePool`](@ref) across
every `(a, r₀)` point of a figure and pass it as `pool`. The full-z sweep above is the
demonstration — a LARGER problem than the cold one, on an already-warm pool, in a third of the
time.

# THE GMRES COUNT DOES NOT GROW WITH m — SETTLED AT THE PAPER'S OWN DEPTH

An earlier pass of this node reported that the count at the finest AMR grid GREW with `m`
(9…20 over `m = 0…20`). That was measured before the all-dimension 2:1 closure, which moved the
iteration counts and not only the meshes, and it was WITHDRAWN. It is now settled by a re-run of
the whole 21-mode sweep at the paper's own 9 grids
(`evidence/mmode-parallel-driver/artefact_21modes_v3.json`), and the answer is that **it does
not grow**:

| | `m = 0…20`, iterations at the finest grid | max over EVERY level and mode |
|---|---|---|
| production HALF domain | 3 3 3 3 3 3 2 2 2 2 2 2 3 3 4 3 3 3 3 4 4 | **7** |
| FULL-z lattice | 2 2 2 2 2 2 2 2 2 3 2 2 2 2 3 3 2 3 3 3 3 | **7** |

**main.tex:720's "less than 10 iterations per AMR level" therefore HOLDS for EVERY mode of the
sweep, on both lattices** — not merely for the paper's one exemplary mode. `m = 2` runs
`[1, 7, 6, 5, 5, 4, 4, 4, 3]` per level and `m = 20` runs `[1, 5, 4, 5, 6, 4, 5, 5, 4]`: the
count falls with refinement depth at every `m`, which is the behaviour the paper describes.
The dof cascade is 320, 632, 1014, 1506, 2120, 2868, 3762, 4814, 6036 (72 elements) on the half
domain and 480 … 10292 (131 elements) on the full-z lattice. All 21 modes converged on both.

The CONFIGURATION fact that the earlier growth was attributed to still holds and is still
recorded per mode: the shipped vtu config sets `ComplexShift: 0` (yaml:194,203; main.tex:680),
so the preconditioner is m-agnostic while the operator's Helmholtz term scales as `(mΩ)²`. It
simply does not cost iterations at this configuration and depth. Each result's
[`ModeStamp`](@ref) carries the shift it used.

# FAILURES, DEATHS AND DEADLINES

A mode that throws is reported BY ID and never swallowed: it is `@warn`ed with its remote
backtrace, recorded in `set.failures`, printed by `show`, and every OTHER mode still returns.
The same is true of a worker that DIES — the failure carries the dead worker's id and the
corpse is pruned from the pool rather than handed out again.

`solve_kwargs` is forwarded VERBATIM to [`solve_mode_amr`](@ref) — it is how a caller reaches
a keyword the driver does not name, `horizon_hold = true` (the row-58 high-spin cure) above
all. It is validated on the MASTER against `solve_mode_amr`'s own method table before a single
worker starts, refused if it repeats a keyword the driver owns, and recorded in every
[`ModeStamp`](@ref) so an archived sweep says which keywords produced it. MEASURED at
`a = 0.5`, `r₀ = 10M`, `m = 2`: the first refinement step alone gives 632 dofs by default
against 614 with `horizon_hold = true`, so the keyword is visible in the returned mesh
(measured on the all-dimension 2:1 closure, GRTools 13ea24c; it was 593 against 575 before that
commit, and `test_parallel.jl` was re-pinned at 93fcc41).

`timeout_s` bounds the whole sweep. Without it a hung mode hangs the sweep forever, because
`remotecall_fetch` has no deadline; with it, the outstanding modes are cancelled the only way a
remote call can be cancelled — by REMOVING their workers — each is reported by id, and
`run.timed_out` is `true`. A supplied pool loses those workers, so check `isopen(pool)`
afterwards.

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
                              solve_kwargs::NamedTuple = (;),
                              timeout_s = nothing,
                              log = nothing,
                              sample_interval::Real = 5.0,
                              project::AbstractString = _pf_project(),
                              exeflags::AbstractVector{<:AbstractString} = String[]) where {RT}
    blas_threads >= 1 || throw(ArgumentError(
        "solve_modes_parallel: blas_threads must be ≥ 1, got $blas_threads"))
    timeout_s === nothing || timeout_s > 0 || throw(ArgumentError(
        "solve_modes_parallel: timeout_s must be positive or `nothing`, got $timeout_s"))
    _pf_check_solve_kwargs(solve_kwargs)      # on the MASTER, before any worker is started
    mlist = _pf_modes(ms)
    own = pool === nothing
    if own
        p = mode_pool_size(length(mlist); nworkers = nworkers,
                           threads_per_worker = threads_per_worker,
                           blas_threads = blas_threads, budget = budget)
        p.clamped &&
            @warn "solve_modes_parallel: the compute budget clamped the request (compute-usage skill)" requested_workers = p.requested granted_workers = p.nworkers requested_threads = threads_per_worker threads_per_worker = p.threads_per_worker requested_blas = blas_threads blas_threads = p.blas_threads cpus_per_worker = p.cpus_per_worker requested_budget = budget budget = p.budget cpu_budget = PARALLEL_CPU_BUDGET
    else
        _pf_prune!(pool, "solve_modes_parallel")
        isopen(pool) || throw(ArgumentError(
            "solve_modes_parallel: the supplied ModePool has no live worker process " *
            "(closed, or every worker died / was removed) — start a new one with mode_pool"))
        need = max(pool.threads_per_worker, Int(blas_threads))
        need <= pool.cpus_per_worker || throw(ArgumentError(
            "solve_modes_parallel: this sweep asks for blas_threads = $blas_threads, so each " *
            "of the pool's $(length(pool)) workers would occupy $need CPUs, but the pool was " *
            "budgeted for $(pool.cpus_per_worker) CPU(s) each — " *
            "$(length(pool) * need) CPUs against a budget of $(pool.budget) and an owner cap " *
            "of $PARALLEL_CPU_BUDGET. Start the pool with blas_threads = $blas_threads, or " *
            "pin fewer BLAS threads."))
        p = (nworkers = length(pool), requested = pool.requested,
             threads_per_worker = pool.threads_per_worker, blas_threads = pool.blas_threads,
             cpus_per_worker = pool.cpus_per_worker, budget = pool.budget,
             clamped = pool.clamped)
    end
    what = "solve_modes_parallel($(p.nworkers) × $(p.cpus_per_worker) CPU(s))"
    t0 = time()
    share = p.nworkers * p.cpus_per_worker / Sys.CPU_THREADS
    # A REUSED pool's processes are already running and are therefore already inside the live
    # reading: adding their share again would double count them and refuse legal launches.
    # Only a launch that CREATES capacity adds it.
    guard_share = own ? share : 0.0
    before, mon = _pf_open_run(what, guard_share, load_guard, t0, log, sample_interval)
    results = ModeResult{T,RT}[]
    failures = ModeFailure[]
    thepool = pool
    startup = 0.0
    timed_out = false
    after = before
    try
        _pf_start!(mon)
        if own
            ts = time()
            thepool = mode_pool(p.nworkers; threads_per_worker = p.threads_per_worker,
                                blas_threads = p.blas_threads, budget = p.budget,
                                load_guard = load_guard, project = project,
                                exeflags = exeflags)
            startup = time() - ts
        end
        raw = Vector{Any}(undef, length(mlist))
        timed_out = _pf_dispatch!(raw, copy(thepool.pids), mlist, timeout_s = timeout_s,
                                  what = what) do w, m
            Distributed.remotecall_fetch(_pf_solve, w, orbit, config, m;
                                         iterations = iterations, T = T,
                                         criterion = criterion,
                                         level_options = level_options,
                                         blas_threads = blas_threads,
                                         solve_kwargs = solve_kwargs)
        end
        _pf_prune!(thepool, what)          # a death shrinks the pool; never hand a corpse out
        _pf_collect!(results, failures, mlist, raw, what)
    finally
        # the pool is ALWAYS given back and the AFTER reading is ALWAYS taken, throw or not:
        # a leaked worker is 2.5 GB of RSS and a core the owner did not lend.
        after = _pf_close_run!(mon, before,
                               () -> (own && thepool !== nothing && close(thepool)), what)
    end
    rec = ModeRunRecord(true, p.nworkers, p.threads_per_worker, Int(blas_threads),
                        p.cpus_per_worker, p.budget, p.requested, p.clamped, own,
                        Float64(load_guard), share, guard_share, before, after,
                        _pf_max_load(mon, load_guard, what), solve_kwargs,
                        timeout_s === nothing ? nothing : Float64(timeout_s), timed_out,
                        startup, time() - t0, mon.samples)
    return ModeSet{T,RT}(orbit, config, mlist, results, failures, rec)
end

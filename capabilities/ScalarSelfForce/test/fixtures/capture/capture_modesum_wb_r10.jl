#!/usr/bin/env julia
# capture_modesum_wb_r10.jl — capture script for test/fixtures/modesum_wb_r10.jl
# (ESSF node `mode-sum-total-force`, GLOBAL_DAG row 46).
#
# WHAT IT RECORDS, AND WHY IT IS A RECORDING AND NOT AN ORACLE
# ------------------------------------------------------------
# The per-mode radial self-force `F_r^m` (and the SIGNED difference between the two highest AMR
# levels, the CRP's own resolution measure of main.tex:722) for the five Warburton–Barack
# comparison spins at `r0 = 10M`, `m = 0…20`, at the sealed production configuration. The
# oracle is THIS PACKAGE's own row-42/43/44/47 pipeline — the numbers are a RECORDING of it at
# one commit, so that the node's Warburton–Barack gate (an assertion about the mode SUM and its
# large-m tail, given per-mode forces that rows 44/45/47 already gate against spectre, against
# a BigFloat finite-difference oracle and against the printed Heffernan parameters) runs in
# seconds instead of ~30 CPU-minutes. `test_modesum.jl` gate (e) re-solves ONE mode per spin
# live at every suite run and requires this table to match, so the recording cannot drift away
# from the code unnoticed.
#
# The EXTERNAL oracle of the node is the printed Warburton–Barack Table III
# (arXiv:1003.1860 `Warburton.tex:620-645`), which is a different fixture
# (`test/fixtures/wb_reference.jl`) and is string-equal to the tex mirror.
#
# RUN (≈ 12 min wall on 6 worker processes, ≈ 30 CPU-minutes; the compute policy of
# `solve_modes_parallel` applies and is enforced by the package itself):
#
#     cd <GRTools>
#     JULIA_PKG_OFFLINE=true julia --project=capabilities/ScalarSelfForce \
#         capabilities/ScalarSelfForce/test/fixtures/capture/capture_modesum_wb_r10.jl \
#         > capabilities/ScalarSelfForce/test/fixtures/modesum_wb_r10.jl
#
# Options: `--workers N` (default 6), `--json PATH` (also dump the raw artefact),
# `--log PATH` (the load monitor log of the fan-out).

using ScalarSelfForce

const R0         = 10.0
const SPINS      = [-0.9, -0.5, 0.0, 0.5, 0.9]     # the printed WB columns the CRP compares at
const MMODES     = 0:20                            # CRP_PAPER_MMODES, main.tex:703
const LIVE_MODE  = 2                               # the mode test_modesum.jl (e) re-solves live

function parse_args(argv)
    o = Dict{String,Any}("workers" => 6, "json" => nothing, "log" => nothing)
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--workers"; o["workers"] = parse(Int, argv[i+1]); i += 2
        elseif a == "--json"; o["json"] = argv[i+1]; i += 2
        elseif a == "--log";  o["log"]  = argv[i+1]; i += 2
        else error("capture_modesum_wb_r10: unknown argument $a")
        end
    end
    return o
end

git_head() = try
    strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String))
catch
    "UNKNOWN"
end

function main(argv)
    o = parse_args(argv)
    base = crp_paper_config(; spin = 0.0, orbital_radius = R0, m_mode = 0,
                            inner_radius = CircularOrbit(0.0, R0).rp, P = 4, L = 0)
    iters = base.amr.iterations                    # the sealed default: iters refinements
    tables = Dict{Float64,Any}()
    # ONE pool across all five spins: a fresh worker pays its whole JIT before its first solve.
    mode_pool(o["workers"]; threads_per_worker = 1, blas_threads = 1) do pool
        for a in SPINS
            orb = CircularOrbit(a, R0)
            cfg = crp_paper_config(; spin = a, orbital_radius = R0, m_mode = 0,
                                   inner_radius = orb.rp, P = 4, L = 0)
            t0 = time()
            set = solve_modes_parallel(orb, cfg, MMODES; pool = pool, iterations = iters,
                                       log = o["log"])
            isempty(set.failures) ||
                error("capture_modesum_wb_r10: modes failed at a = $a: " *
                      string([f.m for f in set.failures]))
            fr, dl, ft, ftdl, red = per_mode(set)
            tables[a] = (Fr = fr, delta = dl, Ft = ft, ft_delta = ftdl, reduction = red,
                         wall = time() - t0,
                         nworkers = set.run.nworkers, max_load = set.run.max_load,
                         ndofs = set[last(MMODES)].levels[end].mesh.ndofs)
            @info "capture_modesum_wb_r10: spin done" a wall = time() - t0
        end
    end
    emit(stdout, tables, iters, git_head())
    if o["json"] !== nothing
        open(o["json"], "w") do io
            print(io, json_artefact(tables, iters, git_head()))
        end
    end
    return nothing
end

# The per-mode final value and the SIGNED difference between the two highest AMR levels,
# through the SEALED row-44 formula alone (`mode_self_force` on the sealed row-43 per-level
# diagnostics). This script deliberately does NOT call the mode-sum node it feeds: a recording
# produced by the code it is used to test would be worth nothing.
function per_mode(set)
    fr, dl, ft, ftdl = Float64[], Float64[], Float64[], Float64[]
    reds = Symbol[]
    for m in sort(collect(Int, set.ms))
        r = set[m]
        # the SEALED row-44 level route: one call, both components, no chart handling here
        h = mode_self_force_history(set.orbit, mode_config(set.config, m), m, r.levels;
                                    x_p = r.x_p)
        vr = [f.Fr for f in h]
        vt = [f.Ft for f in h]
        push!(fr, vr[end]);  push!(dl, vr[end] - vr[end - 1])
        push!(ft, vt[end]);  push!(ftdl, vt[end] - vt[end - 1])
        push!(reds, h[end].reduction)
    end
    length(unique(reds)) == 1 || error("capture_modesum_wb_r10: mixed reductions $(unique(reds))")
    return fr, dl, ft, ftdl, only(unique(reds))
end

vecrepr(v) = "[" * join((repr(x) for x in v), ", ") * "]"

function emit(io, tables, iters, sha)
    println(io, """
# modesum_wb_r10.jl — BAKED per-mode radial self-force table for the Warburton–Barack
# comparison points of ESSF node `mode-sum-total-force` (GLOBAL_DAG row 46). DO NOT EDIT BY
# HAND: regenerate with test/fixtures/capture/capture_modesum_wb_r10.jl.
#
# WHAT THIS IS: a RECORDING of this package's own m-mode pipeline (rows 42/43/44/47) at the
# commit named below — NOT an external oracle. The external oracle of the node is the printed
# Warburton–Barack Table III (arXiv:1003.1860 Warburton.tex:620-645), baked separately and
# string-equal to the tex in test/fixtures/wb_reference.jl. test_modesum.jl gate (e) re-solves
# MSWB_LIVE_MODE live for every spin and requires this table to reproduce, so the recording
# cannot silently drift from the code.
#
# CONFIGURATION (identical for every spin except the spin itself):
#   crp_paper_config(spin = a, orbital_radius = 10.0, m_mode = 0,
#                    inner_radius = CircularOrbit(a, 10.0).rp, P = 4, L = 0)
#   per-mode config from the sealed `mode_config` (only MModeNumber and the m = 0 outer-BC
#   branch move); solve_modes_parallel over m = 0:20, one process per mode,
#   OPENBLAS_NUM_THREADS = 1; solve_mode_amr with `iterations` refinement steps
#   ⇒ `iterations + 1` grids. Reduction: the one-sided `:trace` of the sealed per-level
#   diagnostics (a ModeResult ships no live ModeSolution, so the spectre-parity `:mean` is not
#   reachable from the fan-out; the two differ by the DG trace spread, far below the
#   resolution floor).
#
# MSWB_FR[a][k]    = F_r^m at m = k-1, final AMR level, normalized (M²/q²)F_r (q = M = 1).
# MSWB_DELTA[a][k] = the SIGNED difference between the final and the penultimate level —
#                    the CRP's resolution measure (main.tex:722). Its magnitude is the
#                    per-mode resolution error; its sign is what makes the total's own
#                    resolution error a coherent sum.
#
# STANDALONE: pure Base, no package or stdlib deps (it is loaded into a fresh module).
""")
    println(io, "const MSWB_COMMIT     = ", repr(sha))
    println(io, "const MSWB_CAPTURED   = ", repr(Libc.strftime("%Y-%m-%d", time())))
    println(io, "const MSWB_R0         = ", repr(R0))
    println(io, "const MSWB_MMODES     = 0:", last(MMODES))
    println(io, "const MSWB_ITERATIONS = ", iters, "        # AMR refinement steps")
    println(io, "const MSWB_LEVELS     = ", iters + 1, "        # grids = iterations + 1")
    println(io, "const MSWB_REDUCTION  = :", tables[SPINS[1]].reduction,
            "        # what the SEALED row-44 level route actually used")
    println(io, "const MSWB_LIVE_MODE  = ", LIVE_MODE)
    println(io, "const MSWB_SPINS      = ", vecrepr(SPINS))
    println(io)
    println(io, "const MSWB_FR = Dict{Float64,Vector{Float64}}(")
    for a in SPINS
        println(io, "    ", repr(a), " => ", vecrepr(tables[a].Fr), ",")
    end
    println(io, ")\n")
    println(io, "const MSWB_DELTA = Dict{Float64,Vector{Float64}}(")
    for a in SPINS
        println(io, "    ", repr(a), " => ", vecrepr(tables[a].delta), ",")
    end
    println(io, ")\n")
    println(io, "const MSWB_FT = Dict{Float64,Vector{Float64}}(")
    for a in SPINS
        println(io, "    ", repr(a), " => ", vecrepr(tables[a].Ft), ",")
    end
    println(io, ")\n")
    println(io, "const MSWB_FT_DELTA = Dict{Float64,Vector{Float64}}(")
    for a in SPINS
        println(io, "    ", repr(a), " => ", vecrepr(tables[a].ft_delta), ",")
    end
    println(io, ")")
    return nothing
end

function json_artefact(tables, iters, sha)
    io = IOBuffer()
    println(io, "{")
    println(io, "  \"commit\": ", repr(sha), ",")
    println(io, "  \"iterations\": ", iters, ",")
    println(io, "  \"r0\": ", R0, ",")
    println(io, "  \"spins\": [")
    for (i, a) in enumerate(SPINS)
        t = tables[a]
        println(io, "    {\"a\": ", a,
                ", \"wall_s\": ", repr(t.wall), ", \"nworkers\": ", t.nworkers,
                ", \"max_load\": ", repr(t.max_load), ", \"ndofs_final\": ", t.ndofs,
                ", \"Fr_modes\": ", vecrepr(t.Fr),
                ", \"deltas\": ", vecrepr(t.delta),
                ", \"Ft_modes\": ", vecrepr(t.Ft),
                ", \"Ft_deltas\": ", vecrepr(t.ft_delta), "}", i == length(SPINS) ? "" : ",")
    end
    println(io, "  ]")
    println(io, "}")
    return String(take!(io))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)

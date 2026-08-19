# Run every reproduction driver listed in reproduce/MANIFEST.toml, each in its own
# subprocess, and summarise pass / fail / skip.
#
#     julia --project=reproduce reproduce/run_all.jl [--check|--regen] [--tier fast|heavy|all]
#                                                    [--paper <id>] [--only <id>] [--threads <n>]
#
# Tier ceiling: `--tier fast` (default) runs fast rows only; `--tier heavy` / `all` also
# run the heavy ones.  Threads default to 4 and are capped at 8 (compute-usage policy).
# Exit status is non-zero when any driver fails (exit 1/2 or crash); a driver that exits 3
# (HOLE/skip, reason printed) counts as skipped, not failed.

include(joinpath(@__DIR__, "_lib", "repro_harness.jl"))
using .ReproHarness
using TOML
using Printf

const REPRO_ROOT = ReproHarness.REPRO_ROOT
const MANIFEST = joinpath(REPRO_ROOT, "MANIFEST.toml")

opts = ReproHarness.parse_args()
tier = opts.tier === :all ? :heavy : opts.tier          # run_all treats `all` as the heavy ceiling
if opts.tier === :all && !any(a -> startswith(a, "--tier"), ARGS)
    tier = :fast                                          # default ceiling for run_all is fast
end
isempty(opts.extra) || (println("run_all: unknown arguments ", opts.extra); exit(2))

manifest = TOML.parsefile(MANIFEST)
rows = get(manifest, "result", Any[])

# ── compute-usage guard (report; the cap is enforced by the thread clamp) ──────────────
ncpu = Sys.CPU_THREADS
load = Sys.loadavg()
@printf("[run_all] nproc=%d  load=%.2f/%.2f/%.2f (%.0f%% of cores, 1-min)  threads/driver=%d  tier<=%s  mode=%s\n",
        ncpu, load[1], load[2], load[3], 100 * load[1] / ncpu, opts.threads, tier, opts.mode)
if load[1] / ncpu > 0.9
    println("[run_all] WARNING: system load above 90% of cores — compute-usage policy says do not launch more work")
end

norm_paper(p) = ReproHarness.paper_dir(String(p))
selected = filter(rows) do r
    (opts.paper === nothing || norm_paper(r["paper"]) == norm_paper(opts.paper)) &&
    (opts.only === nothing || r["id"] == opts.only) &&
    ReproHarness.tier_rank(Symbol(get(r, "tier", "fast"))) <= ReproHarness.tier_rank(tier)
end
println("[run_all] $(length(selected))/$(length(rows)) manifest rows selected")
isempty(selected) && (println("[run_all] nothing to run"); exit(0))

julia = joinpath(Sys.BINDIR, Base.julia_exename())
modeflag = opts.mode === :regen ? "--regen" : "--check"

struct Outcome
    paper::String; id::String; kind::String; tier::String; status::String; code::Int; seconds::Float64
end
outcomes = Outcome[]

for r in selected
    script = joinpath(REPRO_ROOT, r["script"])
    id = r["id"]; paper = r["paper"]
    println("\n" * "="^100)
    println("[run_all] $(paper) :: $(id)  ($(get(r, "kind", "?")), tier=$(get(r, "tier", "fast")))  →  $(r["script"]) $(modeflag)")
    println("="^100)
    if !isfile(script)
        println("[run_all] FAIL: script not found: $(script)")
        push!(outcomes, Outcome(paper, id, get(r, "kind", "?"), get(r, "tier", "fast"), "FAIL(no-script)", 2, 0.0))
        continue
    end
    if opts.mode === :check && haskey(r, "expected") && !isfile(joinpath(REPRO_ROOT, r["expected"]))
        println("[run_all] FAIL: golden missing: $(r["expected"]) (run with --regen)")
        push!(outcomes, Outcome(paper, id, get(r, "kind", "?"), get(r, "tier", "fast"), "FAIL(no-golden)", 2, 0.0))
        continue
    end
    cmd = Cmd([julia, "--project=$(REPRO_ROOT)", "--threads=$(opts.threads)", "--startup-file=no",
               script, modeflag, "--tier", string(tier)])
    t0 = time()
    proc = run(ignorestatus(cmd); wait = true)
    dt = time() - t0
    code = proc.exitcode
    status = code == 0 ? (opts.mode === :regen ? "REGEN" : "PASS") :
             code == 3 ? "SKIP" :
             code == 2 ? "FAIL(harness)" :
             code == 1 ? "FAIL" : "FAIL(crash $(code))"
    push!(outcomes, Outcome(paper, id, get(r, "kind", "?"), get(r, "tier", "fast"), status, code, dt))
end

println("\n" * "="^100)
println("[run_all] summary")
println("="^100)
w1 = maximum(length(o.paper) for o in outcomes); w2 = maximum(length(o.id) for o in outcomes)
println(rpad("paper", w1), " | ", rpad("id", w2), " | kind   | tier  | status        | seconds")
println(repeat("-", w1), "-+-", repeat("-", w2), "-+--------+-------+---------------+--------")
for o in outcomes
    println(rpad(o.paper, w1), " | ", rpad(o.id, w2), " | ", rpad(o.kind, 6), " | ", rpad(o.tier, 5), " | ",
            rpad(o.status, 13), " | ", @sprintf("%7.1f", o.seconds))
end
npass = count(o -> o.code == 0, outcomes)
nskip = count(o -> o.code == 3, outcomes)
nfail = length(outcomes) - npass - nskip
println("[run_all] pass=$(npass) skip=$(nskip) fail=$(nfail) of $(length(outcomes))")
exit(nfail == 0 ? 0 : 1)

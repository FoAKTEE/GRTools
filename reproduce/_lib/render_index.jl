#!/usr/bin/env julia
# Render the reproduce/ index table into reproduce/README.md from MANIFEST.toml + the goldens.
#
#     julia --project=reproduce reproduce/_lib/render_index.jl [--check]
#
# The index is a GENERATED view, never hand-authored (research-admission contract: prose
# artifacts are rendered from the canonical record, corrected by editing the record and
# re-rendering).  The canonical record here is MANIFEST.toml plus the golden headers; this
# script rewrites everything between
#
#     <!-- BEGIN GENERATED INDEX ... -->  and  <!-- END GENERATED INDEX -->
#
# leaving the rest of README.md untouched.  `--check` rewrites nothing and exits 1 when the
# committed README does not match what would be rendered — the form used in CI and by the
# verification command of the GRR::reproduce-readme-render ledger row.
#
# Every number in the rendered block is derived, never typed: result counts and the
# status/tier/kind breakdowns come from the manifest, and the golden column is measured on
# the files themselves (present? how many registered values? how large?).

module RenderIndex

using Printf
using TOML

include(joinpath(@__DIR__, "repro_harness.jl"))
using .ReproHarness: json_read

const REPRO_ROOT = normpath(joinpath(@__DIR__, ".."))
const MANIFEST_PATH = joinpath(REPRO_ROOT, "MANIFEST.toml")
const README_PATH = joinpath(REPRO_ROOT, "README.md")
const BEGIN_MARK = "<!-- BEGIN GENERATED INDEX -->"
const END_MARK = "<!-- END GENERATED INDEX -->"

# Human labels for the papers, so the index is readable without opening each driver.
# Keyed by the manifest `paper` field; an unknown paper renders with an empty label rather
# than failing, and `--check` in CI will show the gap.
const PAPER_TITLES = Dict(
    "_selftest"      => "harness self-test (not a paper)",
    "gr-qc-0703028"  => "Berti, Cardoso & Starinets — QNM tables",
    "1010.0758"      => "Fujita & Hikida — Kerr geodesics / fundamental frequencies",
    "1010.0759"      => "Fujita, Hikida & Tagoshi — Teukolsky MST fluxes",
    "1301.1591"      => "Chrzanowski-Cohen-Kegeles metric reconstruction",
    "1811.09628"     => "Pound et al. — second-order self-force",
    "2303.16798"     => "Piovano et al. — generic linear-in-spin fluxes",
    "2506.20726"     => "scalar self-force / effective source",
    "2510.02433"     => "modified-gravity ringdown",
    "2601.19123"     => "hybrid EMRI waveform generation",
    "2603.12189"     => "post-Newtonian self-force comparison",
    "2603.13482"     => "elliptic scalar self-force",
    "2608.17003"     => "EMRI in a gas disc — spiral-wave torque (mission paper)",
)

const STATUS_ORDER = ["solid", "preliminary", "hole", "future"]
const STATUS_MARK = Dict("solid" => "solid", "preliminary" => "prelim",
                         "hole" => "HOLE", "future" => "future")

struct Row
    paper::String
    id::String
    kind::String
    tier::String
    status::String
    script::String
    expected::String
    depends::Vector{String}
    notes::String
end

function load_rows()
    manifest = TOML.parsefile(MANIFEST_PATH)
    rows = Row[]
    for r in get(manifest, "result", Any[])
        push!(rows, Row(String(r["paper"]), String(r["id"]), String(get(r, "kind", "?")),
                        String(get(r, "tier", "fast")), String(get(r, "status", "future")),
                        String(get(r, "script", "")), String(get(r, "expected", "")),
                        String[String(d) for d in get(r, "depends", String[])],
                        String(get(r, "notes", ""))))
    end
    return rows
end

"Measured facts about a row's golden: presence, registered-value count, size in KB."
struct Golden
    present::Bool
    nvalues::Int
    kbytes::Float64
    commit::String
end

function inspect_golden(row::Row)
    isempty(row.expected) && return Golden(false, 0, 0.0, "")
    path = joinpath(REPRO_ROOT, row.expected)
    isfile(path) || return Golden(false, 0, 0.0, "")
    txt = read(path, String)
    doc = json_read(txt)
    vals = (doc isa AbstractDict && haskey(doc, "values")) ? doc["values"] : Dict{String,Any}()
    commit = (doc isa AbstractDict) ? String(get(doc, "generated_at_commit", "")) : ""
    return Golden(true, length(vals), ncodeunits(txt) / 1024, commit)
end

count_by(rows, f) = begin
    d = Dict{String,Int}()
    for r in rows
        d[f(r)] = get(d, f(r), 0) + 1
    end
    d
end

"`3 solid · 1 prelim` — status counts in a fixed order, omitting the zeros."
function status_cell(rows)
    c = count_by(rows, r -> r.status)
    parts = String[]
    for s in STATUS_ORDER
        n = get(c, s, 0)
        n > 0 && push!(parts, "$(n) $(STATUS_MARK[s])")
    end
    isempty(parts) ? "—" : join(parts, " · ")
end

function kind_cell(rows)
    c = count_by(rows, r -> r.kind)
    parts = String[]
    for k in ("figure", "table", "check")
        n = get(c, k, 0)
        n > 0 && push!(parts, "$(n) $(k == "figure" ? "fig" : k == "table" ? "tab" : "chk")")
    end
    isempty(parts) ? "—" : join(parts, " · ")
end

function tier_cell(rows)
    nf = count(r -> r.tier == "fast", rows)
    nh = count(r -> r.tier == "heavy", rows)
    nh == 0 ? "$(nf) fast" : "$(nf) fast · $(nh) heavy"
end

paper_label(p) = startswith(p, "_") ? p : "arXiv:$(p)"
paper_dirname(p) = startswith(p, "_") ? p : "arXiv-$(p)"

function render(rows::Vector{Row})
    io = IOBuffer()
    papers = sort(unique(r.paper for r in rows); by = p -> (startswith(p, "_"), p))
    total_goldens = 0
    total_kb = 0.0
    missing_goldens = String[]        # golden absent on a row that claims to be reproduced — a defect
    pending_goldens = String[]        # golden absent on a future/hole row — the declared state

    println(io, BEGIN_MARK)
    println(io, "<!-- Generated by reproduce/_lib/render_index.jl from MANIFEST.toml + the goldens.")
    println(io, "     Do not edit by hand: edit the manifest (and the ledger the status mirrors),")
    println(io, "     then re-run `julia --project=reproduce reproduce/_lib/render_index.jl`. -->")
    println(io)

    # ── per-paper summary ────────────────────────────────────────────────────────────────
    println(io, "### By paper")
    println(io)
    println(io, "| paper | what | results | kinds | tier | status |")
    println(io, "|---|---|---:|---|---|---|")
    for p in papers
        rs = filter(r -> r.paper == p, rows)
        label = paper_label(p)
        dir = "`$(paper_dirname(p))/`"
        title = get(PAPER_TITLES, p, "")
        @printf(io, "| **%s**<br>%s | %s | %d | %s | %s | %s |\n",
                label, dir, title, length(rs), kind_cell(rs), tier_cell(rs), status_cell(rs))
    end
    rs_all = rows
    @printf(io, "| **total** | %d papers | **%d** | %s | %s | %s |\n",
            count(p -> !startswith(p, "_"), papers), length(rs_all),
            kind_cell(rs_all), tier_cell(rs_all), status_cell(rs_all))
    println(io)

    # ── per-result detail ────────────────────────────────────────────────────────────────
    println(io, "### By result")
    println(io)
    println(io, "`--only <id>` runs one row; `--paper <id>` runs one block.")
    println(io, "The **golden** column is measured on the file, not declared: `n values / size`,")
    println(io, "`(none)` for rows that intentionally store no golden.")
    println(io)
    println(io, "| paper | id | kind | tier | status | driver | golden | depends |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for p in papers
        for r in sort(filter(x -> x.paper == p, rows); by = x -> x.id)
            g = inspect_golden(r)
            if isempty(r.expected)
                gold = "(none)"
            elseif g.present
                total_goldens += 1
                total_kb += g.kbytes
                gold = @sprintf("%d / %.0f KB", g.nvalues, g.kbytes)
            elseif r.status in ("future", "hole")
                gold = "not yet generated"
                push!(pending_goldens, "$(r.paper)::$(r.id)")
            else
                gold = "**MISSING**"
                push!(missing_goldens, "$(r.paper)::$(r.id)")
            end
            deps = isempty(r.depends) ? "—" : join(("`" * d * "`" for d in r.depends), " ")
            @printf(io, "| %s | `%s` | %s | %s | %s | `%s` | %s | %s |\n",
                    paper_label(r.paper), r.id, r.kind, r.tier, STATUS_MARK[r.status],
                    r.script, gold, deps)
        end
    end
    println(io)

    # ── totals ───────────────────────────────────────────────────────────────────────────
    nfast = count(r -> r.tier == "fast", rows)
    c = count_by(rows, r -> r.status)
    @printf(io, "**%d registered results** across %d papers (+ the harness self-test): %s.\n",
            length(rows), count(p -> !startswith(p, "_"), papers), status_cell(rows))
    @printf(io, "%d run at the default `--tier fast` ceiling; %d goldens on disk, %.0f KB total ",
            nfast, total_goldens, total_kb)
    @printf(io, "(cap %d KB each).\n", 200)
    if !isempty(missing_goldens)
        println(io)
        @printf(io, "> **DEFECT — %d golden(s) missing on a row that claims to be reproduced**; ",
                length(missing_goldens))
        @printf(io, "`run_all --check` reports `FAIL(no-golden)` for: %s\n", join(missing_goldens, ", "))
    end
    holes = [r for r in rows if r.status in ("hole", "future")]
    if !isempty(holes)
        println(io)
        @printf(io, "> **Not yet reproduced (%d)** — declared open, not a defect: %s. ",
                length(holes),
                join(("`" * r.paper * "::" * r.id * "` (" * r.status *
                      ", " * r.tier * ((r.paper * "::" * r.id) in pending_goldens ? ", no golden yet" : "") * ")"
                      for r in holes), ", "))
        @printf(io, "%d of them are `heavy`, so the default `--tier fast` sweep does not select them.\n",
                count(r -> r.tier == "heavy", holes))
    end
    println(io, END_MARK)
    return String(take!(io))
end

"Splice `block` into README.md between the markers; returns the new file text."
function splice(readme::AbstractString, block::AbstractString)
    i = findfirst(BEGIN_MARK, readme)
    j = findfirst(END_MARK, readme)
    if i === nothing || j === nothing
        error("README.md has no generated-index markers; add\n\n$(BEGIN_MARK)\n$(END_MARK)\n\nwhere the index belongs")
    end
    first(i) < first(j) || error("README.md markers are out of order")
    # `block` ends with a newline after END_MARK, and the text after the old END_MARK starts
    # with one; keep exactly one so re-rendering an already-current README is a no-op.
    return readme[1:prevind(readme, first(i))] * rstrip(block, '\n') * readme[nextind(readme, last(j)):end]
end

function main(args = ARGS)
    check = "--check" in args
    extra = filter(a -> a != "--check", args)
    isempty(extra) || (println("render_index: unknown arguments ", extra); return 2)

    rows = load_rows()
    block = render(rows)
    readme = read(README_PATH, String)
    updated = splice(readme, block)

    if check
        if updated == readme
            @printf("[render_index] OK: README.md index matches MANIFEST.toml (%d rows)\n", length(rows))
            return 0
        end
        println("[render_index] STALE: reproduce/README.md does not match MANIFEST.toml + the goldens.")
        println("[render_index] re-run: julia --project=reproduce reproduce/_lib/render_index.jl")
        return 1
    end

    if updated == readme
        @printf("[render_index] unchanged: README.md already current (%d rows)\n", length(rows))
    else
        write(README_PATH, updated)
        @printf("[render_index] wrote reproduce/README.md index (%d rows, %d lines in block)\n",
                length(rows), count(==('\n'), block))
    end
    return 0
end

end # module

if abspath(PROGRAM_FILE) == (@__FILE__)
    exit(RenderIndex.main())
end

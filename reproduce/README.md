# reproduce/ — reproduction of published results

Every numerical result that GrTools claims to reproduce (a figure, a table, or a printed
number from a paper) lives here as a **thin driver** plus a **golden** (`expected/*.json`).
`run_all.jl` re-runs them and fails loudly when a number drifts. The registry is
`MANIFEST.toml`; the status column mirrors the mission knowledge ledger.

## Layout

```
reproduce/
  Project.toml, Manifest.toml   one pinned environment: GrTools (path ..) + Plots (+ stdlibs)
  MANIFEST.toml                 registry: [[result]] paper, id, kind, script, tier, expected, depends, status, notes
  run_all.jl                    iterate the manifest; one subprocess per driver; summary table; exit≠0 on failure
  _lib/repro_harness.jl         module ReproHarness: setup!, @repro_result, write_expected/compare_expected, CLI, JSON
  _style/                       presentation only: emri_style.jl, mma_style.jl, manuals/<paper>.md (VLM style manuals)
  _selftest/                    harness self-test (toy result, paper "_selftest")
  arXiv-<id>/                   one directory per paper
      figN.jl, tabN.jl, check_*.jl     thin drivers (zero physics inside)
      expected/<id>.json               goldens (≤ 200 KB each)
      derivations/                     optional: .tex / .wls symbolic derivations that back a check
      README.md                        which figure/table each driver reproduces, tolerances, gaps
```

Paper directories are named `reproduce/arXiv-<id>/` (`<id>` = the arXiv number, e.g.
`arXiv-2608.17003`). Names starting with `_` are harness-internal.

## Running

```bash
# once, offline is fine — the Manifest is pinned
JULIA_PKG_OFFLINE=true julia --project=reproduce -e 'using Pkg; Pkg.instantiate()'

julia --project=reproduce reproduce/run_all.jl --check                      # all fast rows (default tier)
julia --project=reproduce reproduce/run_all.jl --check --tier heavy         # fast + heavy rows
julia --project=reproduce reproduce/run_all.jl --check --paper 2608.17003   # one paper
julia --project=reproduce reproduce/run_all.jl --check --only fig03         # one result id
julia --project=reproduce reproduce/run_all.jl --regen --only fig03         # rewrite that golden
julia --project=reproduce reproduce/run_all.jl --check --threads 8          # threads per driver (max 8)

# a single driver directly
julia --project=reproduce -t4 reproduce/arXiv-2608.17003/fig03.jl --check
```

`run_all.jl` prints the compute-usage guard (`nproc`, load average) before launching and
runs each driver as `julia --project=reproduce --threads=<n> <script> --check|--regen
--tier <ceiling>`. Driver exit codes: `0` pass · `1` mismatch · `2` missing golden /
harness error · `3` HOLE or skip (reason printed; counted as skipped, never as pass).
Any `1`/`2`/crash makes `run_all.jl` exit non-zero.

## Determinism policy

* `ReproHarness.setup!()` pins `BLAS.set_num_threads(1)`; results must be invariant under
  the Julia thread count (`-t1` and `-t4` produce the same golden within `rtol`).
* No RNG in drivers; if an algorithm needs one, seed it explicitly and say so in `notes`.
* The environment is pinned: commit `reproduce/Manifest.toml` with every change; Julia
  **1.12.6** is the reference version (recorded in each golden header). After changing the
  dependency list of any GrTools package, re-run `julia --project=reproduce -e 'using Pkg;
  Pkg.resolve()'` (offline is fine) and commit the updated Manifest.
* Goldens carry provenance: `{paper, id, julia, generated_at_commit, script_sha256, values}`;
  `--check` prints a note when the driver script changed since the golden was generated.
* No absolute paths, no `/tmp` reads, no untracked inputs. External corpora that cannot be
  committed are located through the `Data` package data root (`GRTOOLS_DATADIR`) and the
  driver calls `ReproHarness.hole("...")` (exit 3) when they are absent — never fabricate.

## Adding a result

1. **Physics first, in a package.** Anything that computes a number belongs in a capability
   package (`capabilities/<P>` or `core/<P>`) with tests; a driver only calls package
   functions, arranges the values, and (for figures) plots. Zero physics in scripts. Never
   write the same function twice: promote figure-local helpers into the owning package.
2. Write the driver `reproduce/arXiv-<id>/<name>.jl` from the skeleton:
   ```julia
   include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl")); using .ReproHarness
   opts = ReproHarness.parse_args()        # --check | --regen | --tier | --paper | --only
   ReproHarness.setup!()                   # BLAS single thread + run-log header
   ReproHarness.require_tier(opts, :fast)  # or :heavy — skips (exit 3) above the requested ceiling
   using GrTools.Spacetime
   @repro_result "isco_a0.9" kerr_isco(0.9)
   # … figures: plot with _style/*.jl and save a ≤ 200 KB PNG beside the golden if useful
   exit(ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig03", script = @__FILE__, rtol = 1e-10))
   ```
   Per-name tolerances: `tolerances = Dict("name" => (rtol = 1e-6, atol = 1e-12))`.
3. Generate the golden with `--regen`, inspect the printed values against the paper, run
   `--check`, and add the `[[result]]` row to `MANIFEST.toml` (status from the ledger).
4. **Golden-size rule:** each `expected/*.json` ≤ 200 KB (`write_expected` refuses more).
   Reduce heavy figure data to the plotted arrays before registering; no PDFs; PNG
   reference renders only if ≤ 200 KB.
5. Tier: `fast` = deterministic and < 5 min at `-t4` (runs on every commit); `heavy` =
   anything longer (nightly / manual `--tier heavy`, with the compute-usage self-check).
6. Commit driver + golden + manifest row together; the commit body states what is
   reproduced and the verification command/output.

## Figures

Figure drivers use `_style/emri_style.jl` (matplotlib-rendered papers) or
`_style/mma_style.jl` (Mathematica-rendered papers) and the per-paper VLM style manual in
`_style/manuals/`; see `_style/README.md` for the recurring VLM procedure.

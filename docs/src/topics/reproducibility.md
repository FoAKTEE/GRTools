# [Reproducibility](@id topic-reproducibility)

*This page is the single home for the GrTools reproducibility policy: how the environment is
pinned, what determinism means here, where a number's provenance is recorded, and how to
regenerate every artifact in the tree. The repository-root `REPRODUCIBILITY.md` is a short
pointer to this page; [Fixtures](fixtures.md) and [Data artifacts](data-artifacts.md) hold
the per-artifact detail. Correct a claim by fixing the mechanism it describes, not the prose.*

Three concentric contracts are in force, and they fail independently:

| contract | question it answers | mechanism | breaks CI? |
|---|---|---|---|
| **environment** | does the tree still resolve and load the same code? | committed `Manifest.toml` per environment, pinned Julia | yes (`test`) |
| **fixtures** | do the numbers a *test* asserts still come from the oracle they claim? | `FIXTURES.toml` + tracked capture scripts | yes (`test`) |
| **published results** | do the numbers a *paper* claims still come out of this tree? | `reproduce/` drivers + goldens | yes (`reproduce`) |

## Environment pinning

* **Julia.** `1.12.6` is the reference version: `reproduce/Manifest.toml` is resolved against it and
  every golden header records the Julia version it was generated under. Every package declares
  `julia = "1.10"` as its `[compat]` floor, and `.github/workflows/CI.yml` runs the whole matrix on
  **both** `1.10` (the floor) and `1.12` (the pinned reproduction version). Neither entry is the
  floating `"1"` — a Julia release must be adopted deliberately, not inherited from a runner image.
* **Manifests.** Every environment in the tree ships a committed `Manifest.toml`: the umbrella, all
  18 packages, three test-only environments (`core/Data/test`, `core/Wigner/test`,
  `capabilities/PostNewtonian/test`), `docs/`, and `reproduce/`. Intra-suite dependencies resolve
  through relative `[sources]` paths (`../../core/NumericsCore`, `path = ".."`), so nothing inside
  the suite is fetched from a registry and a checkout is self-contained.
* **Offline.** Every environment instantiates with the network off:

  ```bash
  JULIA_PKG_OFFLINE=true julia --project=reproduce -e 'using Pkg; Pkg.instantiate()'
  ```

  If that fails, the Manifest is stale — re-resolve and commit it, never work around it by going
  online in a run that produces numbers.
* **After a dependency change.** Adding or removing a dep of any package invalidates the resolved
  set: re-run `Pkg.resolve()` in the affected environments (offline is fine) and commit the updated
  Manifests in the same commit as the `Project.toml` edit.
* **Weak dependencies.** `HDF5` and `JLD2` are `Data` weakdeps; an extension activates only for a
  package that declares the trigger in its *own* `[deps]`, and `TestKit.assert_dep_budget` enforces
  who may. See [Data artifacts](data-artifacts.md) — this is why a plain `using Perturbation` never
  drags an I/O stack into the closure.

## Determinism rules

1. **BLAS is single-threaded in reproduction runs.** `ReproHarness.setup!()` calls
   `BLAS.set_num_threads(1)` before any driver arithmetic, removing the largest source of
   run-to-run reassociation.
2. **Thread invariance.** A golden must be identical at `-t1` and `-t4` (within its declared
   `rtol`); the harness self-test is checked both ways. Threaded reductions in the tree are written
   gather-then-ordered-fold for exactly this reason — the `EllipticCore` Schwarz sweep is asserted
   *bitwise* equal to its serial counterpart, and `Kernels`/`matrixfree` carry the same parity gates.
3. **No unseeded randomness.** Drivers use no RNG. Where a test needs varying data it uses a fixed
   deterministic cycle, not `rand` (`test_qnm_fast.jl` band sizes are a `7*trial%21` cycle).
4. **No hidden inputs.** No absolute paths, no `/tmp` reads, no untracked files. External corpora
   are resolved through the `Data` root (`GRTOOLS_DATADIR`, `set_datadir!`) or a named
   `GRTOOLS_*` oracle root, and a driver that cannot find one calls `ReproHarness.hole("...")`
   (exit 3) — it never fabricates a substitute.
5. **Caches are keyed, not stateful.** Global memo tables (`_REFELEMENT_CACHE`, `_ANGMAT_CACHE`,
   `_LAMBDA_CACHE`, `_GSN_COEFF_CACHE`) run get/compute/store under a `ReentrantLock` and are keyed
   by their full input; mutation-tested concurrency gates fail if a lock is removed. A cached run
   and a cold run return the same numbers.
6. **Precision is inherited, never assumed.** Numerical code is type-generic `{T}`; a `BigFloat`
   input under the active `setprecision` keeps its full mantissa. The one Float64 ceiling on the
   angular side is *named* (`angular_precision_floor()`), not implicit — see
   [Precision](precision.md).
7. **Test discovery is explicit.** Four `runtests.jl` files (`EllipticCore`, `Perturbation`,
   `Wigner`, `Harmonics`) list their test files explicitly and assert the list matches the
   directory exactly, so a stray file in a dirty working tree can neither join nor silently leave a
   suite.

## Provenance: where a number's origin is recorded

* **Baked fixtures** — one `<pkg>/test/fixtures/FIXTURES.toml` per package, one row per artifact
  with `sha256`, the tracked `capture` script, the `oracle` and its commit, and a `status` from a
  fixed vocabulary. The roll-up, the status meanings and the honest hole list live on
  [Fixtures](fixtures.md).
* **Out-of-tree oracles** — reached through named roots (`GRTOOLS_FEWDATA_DIR`,
  `GRTOOLS_H1LORENZ_ROOT`, `GRTOOLS_REFPAPER_ROOT`, `GRTOOLS_EFFSOURCE_ROOT`,
  `GRTOOLS_QQNM_COEFF_ROOT`, …) via `TestKit.oracle_root`. When a root is absent the gate records a
  visible skip through `TestKit.oracle_or_skip` instead of quietly passing, and

  ```bash
  GRTOOLS_STRICT_ORACLES=1 julia --project=<pkg> -e 'using Pkg; Pkg.test()'
  ```

  promotes every such skip to a hard failure — the mode to run on a machine that *does* have the
  mirrors, and the honest way to ask "is anything hiding behind a skip?".
* **Re-running a capture safely.** Point the capture at a scratch directory and diff; never
  overwrite in place:

  ```bash
  export GRTOOLS_FIXTURE_OUT_DIR=/tmp/refresh
  julia --project=<oracle env> -t4 <pkg>/test/fixtures/capture/<capture>.jl
  diff /tmp/refresh/<fixture>.jl <pkg>/test/fixtures/<fixture>.jl
  ```

  A nonzero data difference is a **finding**, not a licence to refresh: record it and investigate
  which side moved.
* **Published-result goldens** — each `reproduce/arXiv-<id>/expected/<id>.json` carries
  `{paper, id, julia, generated_at_commit, script_sha256, values}`. `--check` prints a note when the
  driver script's hash no longer matches the one the golden was generated from.

## The `reproduce/` tree contract

```
reproduce/
  Project.toml, Manifest.toml   one pinned environment: GrTools (path ..) + Plots (+ stdlibs)
  MANIFEST.toml                 registry: [[result]] paper, id, kind, script, tier, expected, depends, status
  run_all.jl                    iterate the registry, one subprocess per driver, summary table
  _lib/repro_harness.jl         setup!, @repro_result, write/compare_expected, CLI, JSON
  _style/                       presentation only: plot styles + per-paper VLM manuals
  arXiv-<id>/                   thin drivers + expected/ goldens + README
```

Binding rules (the full text is `reproduce/README.md`):

* **Thin drivers, zero physics.** Anything that computes a number lives in a package with tests; a
  driver only calls package functions, arranges values, and plots. Never write the same function
  twice — a figure-local helper that earns its keep is promoted into the owning package first.
* **Tiers.** `fast` = deterministic and under ~5 min at `-t4`; it runs on every push and pull
  request. `heavy` = anything longer, run manually or nightly with `--tier heavy` after the
  compute-usage self-check.
* **Exit codes.** `0` pass · `1` mismatch · `2` missing golden / harness error · `3` documented
  HOLE or tier skip. A `3` counts as *skipped*, never as a pass; any `1`/`2`/crash fails the job.
* **Size.** Each golden ≤ 200 KB (`write_expected` refuses more); reduce heavy figure data to the
  plotted arrays before registering.
* **Status mirrors the ledger.** The `status` field of a `[[result]]` row copies the mission
  knowledge-ledger node status. Never hand-promote a row: promote the node, then the row.

CI runs the fast tier as its own job:

```yaml
reproduce:   # .github/workflows/CI.yml
  run: julia --project=reproduce reproduce/run_all.jl --check --tier fast
```

with `continue-on-error: false` — a drifted golden breaks the build.

## Regenerating everything

```bash
# 1. environments (offline)
for env in . docs reproduce core/* capabilities/*; do
  JULIA_PKG_OFFLINE=true julia --project=$env -e 'using Pkg; Pkg.instantiate()'
done

# 2. package test suites (what CI's `test` job runs)
julia --project=core/EllipticCore -t4 -e 'using Pkg; Pkg.test()'      # one package
GRTOOLS_STRICT_ORACLES=1 julia --project=capabilities/Perturbation -t4 -e 'using Pkg; Pkg.test()'

# 3. the seven self-checking examples (CI's `examples` job)
for f in examples/0*.jl; do julia --project=. "$f"; done

# 4. published results (CI's `reproduce` job)
julia --project=reproduce reproduce/run_all.jl --check --tier fast
julia --project=reproduce reproduce/run_all.jl --check --tier heavy      # + the long rows
julia --project=reproduce reproduce/run_all.jl --regen --only figNN      # rewrite ONE golden

# 5. documentation
julia --project=docs docs/make.jl
```

A golden is rewritten only with `--regen --only <id>`, after inspecting the printed values against
the paper, and the commit body states what moved and why.

## Known holes

Stated from the ledgers, not from optimism. None of these is hidden behind a passing test — each
is a visible skip, a `[[result]]` row that is not `solid`, or a `FIXTURES.toml` status.

* **Fixture capture scripts that no longer exist (35 rows, `hole-script-lost`).** The largest
  groups are the eight remaining `Perturbation/leg_*.jl` legacy-Teukolsky tables, four `Harmonics`
  spheroidal fixtures, ten `SelfForce` `h1/` corpus slices plus three `src_*.h5` sources, and four
  `Spacetime` reference tables. In every one the oracle is still resolvable, so each is a bounded
  regenerator task on the pattern of `regen_swsh_ref.jl` / `regen_carlson_ref.jl` — not a dead end.
  Four rows are `tracked-notrerunnable` by construction (they froze a pre-refactor code state) and
  three are `tracked-notrerun` (measured cost: the two `l = 3` `h1lorenz` modes need > 29 min solves;
  `rt_coeffs_reference.jl` did not finish a 20-minute sympy `collect`).
* **External corpora are not in the tree.** The FEW `.fewdata` amplitude/flux grids
  (`GRTOOLS_FEWDATA_DIR`), the H1Lorenz `h1` retarded-field corpus (`GRTOOLS_H1LORENZ_ROOT`), the
  `EffectiveSource` C oracle (`GRTOOLS_EFFSOURCE_ROOT`) and the QQNM coefficient tex mirror
  (`GRTOOLS_QQNM_COEFF_ROOT`) live outside the checkout. Suites that need them **skip visibly**
  (and fail under `GRTOOLS_STRICT_ORACLES=1`); the committed fixtures remain authoritative. A
  lazy-artifact plan for the >1 MB files is recorded but not implemented — nothing is deleted
  without the owner.
* **The arXiv:2608.17003 supplementary material is no longer in this repository.** The paper ships
  `paper-expressions.mx` beside a *Mathematica* notebook, and route R1 (a transcription of that
  notebook's output boxes) was once an in-tree oracle. On 2026-08-22 the notebook-derived oracle and
  its off-shell fixture were **removed from this repository** by owner decision
  (`EGD::jd-notebook-purge`); the shipped master system is derived in Julia (route R3) and gated by
  our own independent *Mathematica* derivation (route R2), which never read the notebook. The R1
  lineage, the notebook mirror and the full R1-vs-R2-vs-R3 equivalence evidence are archived in the
  mission repository (emri-disk `ref-code/arXiv-2608.17003-anc/`, `results/emri-gasdisk`
  `paper_EMRI-GASDISK`, `evidence/{R1-vs-R2-gate,jd-equivalence-gate,jd-notebook-purge}`), which is
  the canonical archive. Nothing in this checkout loads `.mx` or a notebook, so no Mathematica
  version dependence is left in the reproduction path.
* **MR19 (arXiv:1811.09628) figures.** Registered as `future`: the published figures need reduced
  data that has to be regenerated before drivers can land.
* **The heavy tier does not run in CI.** Roughly two fifths of the registry is `heavy` (2026-08-20:
  20 of 52 rows); CI checks the fast rows only. Heavy rows are re-checked manually, and a heavy-only
  regression can therefore survive a green build until the next nightly/manual sweep. The live split
  is `grep -c 'tier = "heavy"' reproduce/MANIFEST.toml` — the registry grows weekly, so trust the
  command, not this sentence.
* **Not every registered row is `solid`.** The registry carries `preliminary` and `future` rows
  alongside the solid ones; the `status` column is the honest state, and a `preliminary` row means
  the driver runs and its golden is stable, *not* that the physics claim is settled. It mirrors the
  mission knowledge ledger — promote the node first, then the row.
* **arXiv:2606.04998 (elliptic scalar self-force) has no drivers yet.** The `EllipticCore` solver
  stack is still being built; its figure drivers are a `future` node.

## Related pages

* [Fixtures](fixtures.md) — the per-artifact provenance ledger and its status vocabulary.
* [Data artifacts](data-artifacts.md) — corpora, the schema contract, the injected data root.
* [Precision](precision.md) — where BigFloat is inherited and where the Float64 ceilings are.
* [Conventions](conventions.md) — the sign/normalization conventions a number is stated in.

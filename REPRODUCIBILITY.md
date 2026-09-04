# Reproducibility

Every number this repository publishes — in a test, in an example, in a figure — is meant to come
back the same on another machine, from a checkout alone. This file is the short contract; the full
policy, with the mechanisms and the honest hole list, is
[`docs/src/topics/reproducibility.md`](docs/src/topics/reproducibility.md).

## Three contracts

| contract | question | mechanism | CI job |
|---|---|---|---|
| environment | does the tree resolve and load the same code? | committed `Manifest.toml` per environment, pinned Julia | `test` |
| fixtures | do the numbers a *test* asserts still come from the oracle they claim? | `<pkg>/test/fixtures/FIXTURES.toml` + tracked capture scripts | `test` |
| published results | do the numbers a *paper* claims still come out of this tree? | `reproduce/` drivers + goldens | `reproduce` |

## Environment

* Julia **1.12.6** is the reference version; `[compat]` floor is 1.10. CI runs both `1.10` and
  `1.12` — never the floating `"1"`.
* Every environment (umbrella, 18 packages, 3 test-only envs, `docs/`, `reproduce/`) ships a
  committed `Manifest.toml`. Intra-suite deps resolve through relative `[sources]` paths, so a
  checkout is self-contained.
* Everything instantiates offline:

  ```bash
  JULIA_PKG_OFFLINE=true julia --project=reproduce -e 'using Pkg; Pkg.instantiate()'
  ```

  A failure there means the Manifest is stale — re-resolve and commit it; do not go online for a
  run that produces numbers.

## Determinism

* `ReproHarness.setup!()` pins `BLAS.set_num_threads(1)`.
* Goldens must be invariant under the thread count (`-t1` == `-t4` within the row's `rtol`);
  threaded reductions are gather-then-ordered-fold and gated bitwise against serial.
* No unseeded randomness, no absolute paths, no `/tmp` reads, no untracked inputs. A driver that
  cannot find an external corpus calls `ReproHarness.hole(...)` (exit 3) — it never fabricates.
* Global caches are locked and fully keyed: a warm run and a cold run give the same numbers.
* Numerical code is type-generic `{T}`; BigFloat precision is inherited, and the one Float64
  ceiling on the angular side is named, not implicit.

## Provenance

* **Fixtures** — `<pkg>/test/fixtures/FIXTURES.toml` records `sha256`, capture script, oracle and
  its commit, and a status (`regenerated-identical`, `reverified-identical`, `hole-script-lost`, …).
  Roll-up: `docs/src/topics/fixtures.md`.
* **Out-of-tree oracles** — named `GRTOOLS_*` roots via `TestKit.oracle_root`; a missing root is a
  *visible* skip, and `GRTOOLS_STRICT_ORACLES=1` promotes every such skip to a hard failure.
* **Goldens** — `reproduce/arXiv-<id>/expected/<id>.json` carries
  `{paper, id, julia, generated_at_commit, script_sha256, values}`; `--check` warns when the driver
  changed since the golden was written.
* Re-run a capture into a scratch directory and **diff** — never overwrite in place:

  ```bash
  GRTOOLS_FIXTURE_OUT_DIR=/tmp/refresh julia --project=<oracle env> -t4 <pkg>/test/fixtures/capture/<x>.jl
  diff /tmp/refresh/<fixture>.jl <pkg>/test/fixtures/<fixture>.jl
  ```

  A data difference is a finding, not a licence to refresh.

## Regenerating

```bash
julia --project=core/EllipticCore -t4 -e 'using Pkg; Pkg.test()'   # one package suite
for f in examples/0*.jl; do julia --project=. "$f"; done           # the 7 self-checking examples
julia --project=reproduce reproduce/run_all.jl --check --tier fast # published results (CI job)
julia --project=reproduce reproduce/run_all.jl --check --tier heavy
julia --project=reproduce reproduce/run_all.jl --regen --only figNN  # rewrite ONE golden
julia --project=docs docs/make.jl                                  # documentation
```

## Known holes

Stated from the ledgers. None is hidden behind a passing test — each is a visible skip, a
non-`solid` registry row, or a `FIXTURES.toml` status.

* 35 fixture rows have **no surviving capture script** (`hole-script-lost`); their oracles are still
  resolvable, so each is a bounded regenerator task. 4 more are un-re-runnable by construction
  (they froze a pre-refactor code state) and 3 are deferred on measured cost.
* **External corpora are not in the tree**: FEW `.fewdata` grids, the H1Lorenz `h1` corpus, the
  `EffectiveSource` C oracle, the QQNM tex mirror. Suites needing them skip visibly; the committed
  fixtures remain authoritative.
* The arXiv:2608.17003 supplementary material (**`.mx`** binary and *Mathematica* notebook) is not
  in this repository at all since 2026-08-22 (`EGD::jd-notebook-purge`, owner decision): the shipped
  master system is derived in Julia (route R3) and gated by our own independent re-derivation
  (route R2). The notebook mirror, the route-R1 lineage and the R1-vs-R2-vs-R3 equivalence evidence
  are archived in the mission repository (emri-disk `paper_EMRI-GASDISK`).
* **MR19 (arXiv:1811.09628)** figures are `future` — they need reduced data regenerated first.
* The **heavy tier is not run by CI** (2026-08-20: 20 of 52 registered results; the live split is
  `grep -c 'tier = "heavy"' reproduce/MANIFEST.toml`); a heavy-only regression can survive a green
  build until the next manual sweep.
* **arXiv:2606.04998** has THREE drivers so far — `fig_matrix` (the DG operator sparsity at the
  coarsest grid level), `fig_domain` (Fig. 1, the DG computational domain over the hp-AMR mesh)
  and `fig_linsolv` (Fig. 6, the GMRES staircase). The remaining four figures need the
  `ScalarSelfForce` mode-sum / error-budget nodes, which are still being built.

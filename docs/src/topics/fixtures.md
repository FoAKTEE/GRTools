# Baked fixtures and their capture scripts

*This page is the single home for the GrTools baked-fixture provenance story: what a fixture is,
which capture script produced it, whether that script is tracked, and what happened when it was
re-run. The machine-readable ledger is one `FIXTURES.toml` per package under
`<pkg>/test/fixtures/`; this page is its roll-up. Correct a claim by editing the package's
`FIXTURES.toml`, not this page.*

## Why fixtures exist

Most GrTools test suites are **characterization** suites: a capability was ported from a heavy
legacy package (`emri-disk/packages/*`) or from an upstream mirror (`emri-disk/ref-code/*`), and
the port is gated against the original's *actual output*. Depending on those oracles at test time
would drag Mathematica, GSL/Boost, multi-gigabyte corpora and unmaintained Julia environments into
every CI run, so the oracle output is evaluated **once** on a fixed grid and **baked** into a
tracked fixture. The suites then depend only on the tree.

The cost of that trade is reproducibility: a baked number is only as good as the record of how it
was produced. Two failure modes are specific to it —

1. the capture script lives outside the repository, so the tree cannot regenerate its own fixtures;
2. the oracle changes (or is refactored away) and nobody notices the fixture has drifted.

## The provenance ledger

Every package that ships fixtures carries `<pkg>/test/fixtures/FIXTURES.toml`, one `[[fixture]]`
row per tracked artifact:

| field | meaning |
|---|---|
| `file` | package-relative path of the fixture |
| `sha256`, `bytes` | content pin at ledger-generation time |
| `capture` | the tracked capture script or regenerator (empty = none exists) |
| `oracle` | what was evaluated: legacy package / ref-code mirror / closed form, and which entry points |
| `oracle_commit` | commit of that oracle where it is a repository |
| `status` | see the vocabulary below |
| `note` | the honest detail — including what was *not* verified |

### Status vocabulary

| status | meaning |
|---|---|
| `regenerated-identical` | the capture was re-run on 2026-08-19 and the fixture data came back bit-identical (only capture-date/host comment lines may differ) |
| `reverified-identical` | the original capture is lost, but a **new key-driven regenerator** re-evaluated the oracle at every baked key and every value matched exactly |
| `tracked-generator` | an in-tree generator or regenerator already existed before this pass |
| `tracked-notrerunnable` | the capture script is now tracked, but the oracle *state* it captured no longer exists (it was a pre-refactor snapshot); re-running is impossible **by construction** |
| `tracked-notrerun` | the capture script is now tracked; the re-run was deferred (cost or missing corpus) — a `[FUTURE]` obligation |
| `hole-script-lost` | no capture script exists anywhere; values are frozen — a `[HOLE]` |
| `upstream-copy` | vendored upstream data; no generator applies (provenance is a README/PROVENANCE record) |
| `shared-input` | a generator *input* shared verbatim between the capture and the test, so inputs cannot drift by transcription |

## Where the capture scripts live

Tracked capture scripts sit in `<pkg>/test/fixtures/capture/`; tracked oracle generators that need
a non-Julia engine (Wolfram, python) sit in `<pkg>/test/oracles/`. **Nothing in either directory is
executed by any test suite** — `runtests.jl` includes only `test_*.jl` files in `test/` itself — so
tracking them adds no dependency.

Every tracked capture had its hard-coded absolute paths rewritten into ENV-overridable constants,
with defaults that reproduce the original capture host:

| variable | default | what it points at |
|---|---|---|
| `GRTOOLS_LEGACY_PACKAGES_ROOT` | `/data/haiyangw/claude/emri-disk/packages` | legacy oracle packages |
| `GRTOOLS_REFCODE_ROOT` | `/data/haiyangw/claude/emri-disk/ref-code` | upstream source mirrors |
| `GRTOOLS_TEUKEVOLUTION_ROOT` | `$GRTOOLS_REFCODE_ROOT/TeukEvolution.jl` | the TeukEvolution.jl working tree |
| `GRTOOLS_DGTEUK_ROOT` | `/data/haiyangw/claude/DG-Teulolsky` | the DG time-domain oracle repo |
| `GRTOOLS_EMRIDISK_ROOT` | `/data/haiyangw/claude/emri-disk` | the mission repo root |
| `GRTOOLS_ROOT` | this checkout | the GRTools tree |
| `GRTOOLS_FIXTURE_OUT_DIR` | the fixture's own directory | **set this to a scratch dir to re-run without touching the committed fixture** |

The safe re-run/diff idiom is therefore:

```bash
export GRTOOLS_FIXTURE_OUT_DIR=/tmp/refresh
julia --project=<the oracle's env> -t4 <pkg>/test/fixtures/capture/<capture>.jl
diff /tmp/refresh/<fixture>.jl <pkg>/test/fixtures/<fixture>.jl
```

A nonzero data difference is a **finding**, not a licence to overwrite: record it, do not refresh
the fixture, and investigate which side moved.

## Roll-up (2026-08-19)

| package | rows | regen-identical | reverified-identical | tracked-generator | tracked-notrerun | tracked-notrerunnable | hole-script-lost | upstream-copy | shared-input |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `capabilities/AccretionDisc` | 5 | 5 | – | – | – | – | – | – | – |
| `capabilities/EMRIWaveforms` | 2 | – | – | – | – | 1 | 1 | – | – |
| `capabilities/Harmonics` | 6 | 1 | 1 | – | – | – | 4 | – | – |
| `capabilities/Hybrid` | 1 | – | – | – | – | – | – | 1 | – |
| `capabilities/MetricReconstruction` | 2 | 1 | – | – | – | – | 1 | – | – |
| `capabilities/ModifiedGravity` | 3 | – | – | – | – | – | – | 3 | – |
| `capabilities/Perturbation` | 28 | 3 | – | 3 | 8 | – | 14 | – | – |
| `capabilities/PostNewtonian` | 1 | – | – | 1 | – | – | – | – | – |
| `capabilities/ScalarSelfForce` | 11 | – | – | 10 | – | – | – | 1 | – |
| `capabilities/SelfForce` | 35 | 7 | – | – | 2 | – | 13 | 13 | – |
| `capabilities/Spacetime` | 7 | 2 | – | – | – | – | 4 | 1 | – |
| `core/BigSpecial` | 1 | – | – | – | – | – | 1 | – | – |
| `core/Data` | 1 | – | – | – | – | 1 | – | – | – |
| `core/Kernels` | 5 | 1 | – | 1 | – | 2 | – | – | 1 |
| `core/NumericsCore` | 4 | 2 | 2 | – | – | – | – | – | – |
| `core/Wigner` | 1 | – | – | – | – | – | 1 | – | – |
| **total** | **113** | **22** | **3** | **15** | **10** | **4** | **39** | **19** | **1** |

25 fixture groups were re-derived from their oracle in this pass and **every one came back
bit-identical** — nothing verified here has drifted from the oracle it claims. Concretely:

* `dg_lgl_ref.jl`, `numericscore_fd_uniform_ref.jl`, `sphere_collocation_ref.jl`,
  `npghp_np0_fixture.jl`, `td_hphc_recon_upstream.jl` — byte-identical files;
* `spin_cg_reference.jl`, `d2hdr2_oracle.jl`, `spacetime_spinforce_ref.jl`, `fdspa_ref.jl` —
  identical data, differing only in a capture-date / host / repo-state comment line;
* the five `AccretionDisc` `h1-l{2..6}m2.h5` windows — every HDF5 dataset and attribute equal
  (max abs difference `0.0`); only HDF5 container metadata differs;
* `td_qnm_a0_l2_m2_n0.jl` — the re-run reproduced the intermediate text source byte-for-byte and
  its sha256 matches the value the fixture header records;
* `numericscore_swsh_ref.jl`, `spherical_harmonics_ref.jl` (128 values) and
  `numericscore_carlson_ref.jl` (49 literals) — re-evaluated against the live legacy oracle by the
  two new regenerators `capture/regen_swsh_ref.jl` and `capture/regen_carlson_ref.jl`;
* the `MetricReconstruction` KLC literals — all 22 scalars printed by the capture found verbatim in
  the test bodies;
* six of the eight `SelfForce` `h1lorenz/h1-l*.h5` modes (`l0m0, l1m0, l1m1, l2m0, l2m1, l2m2`) —
  every dataset and attribute equal; the two `l = 3` modes exceeded the re-run budget.

One mapping error was caught by the re-run and is worth recording, because a careless refresh would
have corrupted a fixture: `capture_h1.jl` originally wrote into `test/fixtures/h1/`, but that
directory was later **renamed** to `h1lorenz/` when the six-column SecondOrderRicci corpus took the
name `h1` (`test/lorenz_solver.jl:4-8`, `test/second_order_readcorpus.jl:6-9`). Diffing the re-run
against `h1/` shows a spurious structural mismatch (16 vs 24 rows, two extra datasets); against
`h1lorenz/` it is bit-identical. The tracked copy of the script writes to `h1lorenz/`.

## The holes, stated plainly

* **`hole-script-lost` (39 rows).** The largest groups are the twelve `Perturbation/leg_*.jl`
  legacy-Teukolsky fixtures (headers say "regenerate via `scratchpad gen.jl`"; no such file exists
  anywhere), the four remaining `Harmonics` spheroidal fixtures (`capture_swsh.jl`), the ten
  `SelfForce` `h1/` six-column corpus slices and the three `src_*.h5` source fixtures
  (`capture_ms.jl`), and the four `Spacetime` KerrGeodesics/SpinningOrbit reference tables. In every
  one of these the **oracle is still resolvable**, so each is a bounded regenerator task on the
  pattern already demonstrated by `regen_swsh_ref.jl` / `regen_carlson_ref.jl` — not a dead end.
* **`tracked-notrerunnable` (4 rows).** `s4a`/`s4b`/`s4c`-style captures froze a *pre-refactor*
  code state (the FEW-namespace mode-sum before the hoist into `Kernels`, the FEW HDF5 loaders
  before the move into `Data`, the inline `coupling(...)` spec in `r2.jl` before it became
  `R2_COUPLING_SPEC`). Their scripts are tracked for the record, but they cannot be re-run: the
  code they measured is the code the change removed. `fdspa_ref.jl` is the happy exception — its
  front door still exists and delegates to the hoisted kernel, and the re-run matched exactly.
* **`tracked-notrerun` (10 rows).** Scripts are now tracked; re-running them costs more than this
  pass budgeted: the six `td_*` DG time-domain fixtures and `td_hphc_upstream.jl` need full
  `tf ≈ 1000` / 600-step evolutions, the two `l = 3` `h1lorenz/` modes need long `H1Lorenz` solves
  (measured: > 29 min for the first one, versus ~10 s each for `l ≤ 2`), and
  `rt_coeffs_reference.jl` did not finish a 20-minute sympy `collect` over the full jet.

## Related pages

* [Data artifacts](data-artifacts.md) — the corpora the fixtures are sliced from, and the schema
  contract those slices obey.
* [Conventions](conventions.md) — the sign/normalization conventions a fixture's numbers are stated in.

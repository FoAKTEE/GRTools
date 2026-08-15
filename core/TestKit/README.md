# TestKit — shared test harness (test-only)

Tier-1 core, test-infrastructure package: a tolerance-recording
oracle layered on stdlib `Test`, riding the `Data` spine for fixtures.

## Scope

During construction the 281 **predecessor `runtests` remained the AUTHORITATIVE regression oracle**, run
UNCHANGED against the predecessor codes. TestKit is scoped to **NEW tests only** — characterization-gate
snapshots for fork relocations + any new physics. The ~1900 predecessor tolerance literals were NOT
transcribed into `@approx_oracle` (that transcription would itself have been the regression risk). This
new-test net + the predecessor runtests are what made additive-then-sweep de-duplication safe.

## What it owns

This package absorbs, across the predecessor corpus:
- every member's inline `rtol`/`atol` asserts (mechanically, behind `@approx_oracle`);
- the scattered test-fixture loaders (behind `load_fixture`, over the `Data` reader spine).

## Public API

- `@approx_oracle(value, ref; rtol, atol, label)` — `isapprox` that RECORDS its frozen tolerance in
  `oracle_records()` and returns the boolean. Wrap in `@test` to assert. The recording lets CI prove a
  new oracle's tolerance was frozen and never silently loosened.
- `TOLERANCE_CLASSES` / `tolerance_class(:machine|:spectral|:physics|:gpu)` — default `(rtol, atol)`
  bundles for NEW tests, ordered tightest→loosest.
- `assert_can_only_tighten(class_tol, original_tol)` — CI guard that a new tolerance may only TIGHTEN
  the legacy one (`class_tol <= original_tol`), throwing otherwise.
- `load_fixture(parts...; kwargs...)` — load a captured fixture as a numeric matrix, resolving the path
  under the active GrTools data root (`Data.datapath`) and reading via `Data.read_dat`. Configure the
  root with `Data.set_datadir!` / `ENV["GRTOOLS_DATADIR"]`.
- `characterization_snapshot(f, inputs)` / `characterization_gate(candidate, snap; ulp, rtol, atol)` —
  snapshot a reference function's full-precision output (BigFloat stays BigFloat), then gate a
  candidate to the recorded ULP. A perturbed candidate FAILS. This is the additive-then-sweep
  relocation gate.

## Deps

`NumericsCore`, `Data`, `Test`, `Printf`. Light by construction: there is NO legacy oracle (the harness
is greenfield, grounded by self-tests — the macro passes/fails correctly, the tolerance-class ordering,
the can-only-tighten guard, a `Data` fixture round-trip, and the gate catching a perturbation), so the
cold-compile stays fast.

## Grounding

Self-tests in `test/test_testkit.jl` (infra, not physics): `@approx_oracle` pass/fail + tol recording,
tolerance-class ordering (`:machine` strictly tighter than `:physics`), the `class_tol <= original_tol`
guard catching a loosened tol, a fixture round-trip through `Data`, and `characterization_gate`
detecting a 1-ULP and a gross perturbation.

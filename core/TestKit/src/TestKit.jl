"""
    TestKit (GrToolsTestKit)

Tier-1 core, TEST-ONLY shared harness for the GrTools capability packages (design §5.6,
`GRTOOLS_FILE_FUNCTION_DAG.md §GrToolsTestKit`). It layers a thin, tolerance-recording oracle on top of
stdlib `Test` and rides the `Data` spine for fixtures.

Scope (Rev 2 decision, design §5.6): the 281 existing legacy `runtests` remain the AUTHORITATIVE
regression oracle, run UNCHANGED against the legacy facades. This package is scoped to NEW tests only —
characterization-gate snapshots for fork relocations + any new physics — and is what makes the
additive-then-sweep dedup safe.

Public API:
  * `@approx_oracle(value, ref; rtol, atol, label)` — `isapprox` that RECORDS its (frozen) tolerance in
    [`oracle_records`](@ref) and returns the boolean (wrap in `@test`).
  * `TOLERANCE_CLASSES` / `tolerance_class` — `:machine`/`:spectral`/`:physics`/`:gpu` default bundles
    for NEW tests, ordered tightest→loosest.
  * `assert_can_only_tighten(class_tol, original_tol)` — CI guard that a new tol may only TIGHTEN.
  * `load_fixture(parts...)` — captured-fixture loader over `Data.datapath` + `Data.read_dat`.
  * `characterization_snapshot` / `characterization_gate` — snapshot a reference function's
    full-precision output, then gate a candidate to the recorded ULP (a perturbed output FAILS).
  * `assert_dep_budget(project_toml_path; allow, forbid, must_declare, weakdeps_allow)` — §7.7
    gates 1+5: TOML-static dep-budget / extension-trigger checker over a package's own Project.toml.
  * `assert_not_loaded(mods)` — §7.7 gate 2: runtime closure-load poison check over
    `Base.loaded_modules` (binding caveat: run AFTER `using <Pkg>`).

Absorbs (DAG crosswalk): every legacy member's inline `rtol`/`atol` asserts (mechanically, via
`@approx_oracle`) and the scattered test-fixture loaders. Deps: `NumericsCore`, `Data`, `Test`,
`Printf`, `TOML` (stdlib, for the §7.7 hygiene gates) — light by construction (no legacy oracle,
grounded by self-tests), so the cold-compile stays fast.
"""
module TestKit

include("oracle.jl")    # leaf: Test/Printf only — tolerance recording, classes, tighten-guard
include("fixtures.jl")  # rides Data: fixture loader + characterization gate
include("hygiene.jl")   # leaf: TOML/Test only — §7.7 gates 1/2/5 dep-budget + closure-load checkers

export @approx_oracle, oracle_records, reset_oracle_log!, OracleRecord,
    TOLERANCE_CLASSES, tolerance_class, assert_can_only_tighten,
    load_fixture,
    CharacterizationSnapshot, characterization_snapshot, characterization_gate,
    assert_dep_budget, assert_not_loaded

end # module TestKit

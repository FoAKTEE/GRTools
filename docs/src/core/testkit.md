# TestKit

Tier-1 core, **test-only** shared harness for the GrTools capability packages. It layers a
thin, tolerance-recording oracle on top of stdlib `Test` and rides the [Data](data.md) spine for fixtures. Nothing here computes physics: TestKit
exists so that every *new* numerical assertion in the monorepo records its frozen tolerance, every
fork relocation is gated against a full-precision snapshot of the legacy primitive, and every
package's dependency budget is CI-enforced.

## Overview — what the harness guarantees

Three invariants, one per source file:

1. **Frozen tolerances** (`oracle.jl`). `@approx_oracle` is an `isapprox` that *records* the
   `(rtol, atol)` actually used in an append-only log (`oracle_records()`), so a CI pass can prove a
   new oracle's tolerance was registered and never silently loosened. Four named tolerance classes
   order the default budgets tightest → loosest:

   | class | rtol | intended for |
   |---|---|---|
   | `:machine` | ``10^{-14}`` | bit-level / closed-form identities |
   | `:spectral` | ``10^{-10}`` | spectral / series-convergent numerics |
   | `:physics` | ``10^{-6}`` | physics-level agreement, accumulated error |
   | `:gpu` | ``10^{-4}`` | CPU-vs-GPU `sincos`/`fma` divergence (the only class wider than physics) |

   `assert_can_only_tighten(class_tol, original_tol)` enforces the migration rule
   ``\mathrm{tol}_{\text{new}} \le \mathrm{tol}_{\text{legacy}}`` — a relocated test may only
   *tighten* the legacy tolerance, never loosen it.

2. **Characterization gates** (`fixtures.jl`). `characterization_snapshot(f, inputs)` freezes a
   reference function's outputs at full precision (`BigFloat` stays `BigFloat` — no rounding), and
   `characterization_gate(candidate, snap; ulp, rtol, atol)` asserts a candidate reproduces every
   output to the recorded ULP (default `ulp = 0`, i.e. bit-exact for floats). A perturbed candidate
   FAILS. This is the *additive-then-sweep* relocation gate: snapshot the legacy primitive, relocate
   it, gate the new home, only then delete the fork. `load_fixture(parts...)` loads captured fixtures
   through `Data.datapath` + `Data.read_dat` (configure the root with `Data.set_datadir!` or
   `ENV["GRTOOLS_DATADIR"]`; there is deliberately no implicit fallback).

3. **Dependency-hygiene gates** (`hygiene.jl`). Two checkers shared by all 16 packages'
   `test/hygiene.jl`:
   - `assert_dep_budget(project_toml_path; allow, forbid, must_declare, weakdeps_allow)` —
     a TOML-*static* check of the package's own `Project.toml` (`[deps]`/`[weakdeps]`/`[sources]`;
     nothing is loaded). It enforces the closed `allow` budget, the `forbid` poison axes (trailing-`*`
     prefix globs, e.g. `"OrdinaryDiffEq*"`), the extension triggers (`must_declare` names must
     be *direct* hard `[deps]`, since a Julia extension activates only when the loading package
     declares its trigger), and — always — that no `[sources]` path points into a predecessor-code
     tree (`packages/` or `ref-code/`).
   - `assert_not_loaded(mods)` — a *runtime* closure-load check that no poison module is in
     `Base.loaded_modules`. Binding usage contract: call it **after** `using <Pkg>`, and list only
     poison names nothing in the test env legitimately loads.

## Scope decision

The 281 `runtests` of the predecessor codes served as the **authoritative regression oracle**
during construction, run unchanged against their own facades. TestKit is scoped to **new tests
only** — characterization-gate snapshots for fork relocations plus any new physics. The ~1900
predecessor tolerance literals are *not* transcribed into `@approx_oracle` (that transcription would
itself be the regression risk). The regression net is therefore: predecessor runtests (verbatim,
during construction) + new characterization gates, which is what makes additive-then-sweep dedup
safe.

## Module map

TestKit is a single module (no submodules); the internal structure is by file:

| file | depends on | owns |
|---|---|---|
| `src/oracle.jl` | `Test`, `Printf` (leaf) | `@approx_oracle`, `OracleRecord`, `oracle_records`, `reset_oracle_log!`, `TOLERANCE_CLASSES`, `tolerance_class`, `assert_can_only_tighten` |
| `src/fixtures.jl` | `oracle.jl`, `Data` | `load_fixture`, `CharacterizationSnapshot`, `characterization_snapshot`, `characterization_gate` |
| `src/hygiene.jl` | `TOML`, `Test` (leaf) | `assert_dep_budget`, `assert_not_loaded` |

## Dependencies and poison-axis declaration

Hard `[deps]`: `NumericsCore`, `Data`, `Test`, `Printf`, `TOML`. `Test` is a normal dependency here
precisely because this package *is* the shared harness; `TOML` (stdlib, zero cost) powers the static
dep-budget gate. There is no legacy oracle dependency — the harness is greenfield and grounded by
self-tests, so the cold-compile stays light.

TestKit eats its own dog food: its `test/hygiene.jl` runs the hygiene gates on itself with the closed
budget `allow = [NumericsCore, Data, Test, Printf, TOML]` and forbids (in `[deps]`/`[weakdeps]`)
`HDF5`, `JLD2`, `Symbolics`, `OrdinaryDiffEq*`, `CUDA`, `AMDGPU`, `KernelAbstractions`, `Kernels`,
`WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase`, `SpecialFunctions`, `HypergeometricFunctions`,
`QuadGK`; the closure-load check additionally asserts none of the heavy axes (`HDF5`, `JLD2`,
`Symbolics`, `OrdinaryDiffEq*`, `CUDA`, `AMDGPU`, `KernelAbstractions`, `WignerSymbols`) enter
`Base.loaded_modules` after `using TestKit`. In particular, Data's HDF5/JLD2 weakdep extensions are
*not* triggered by TestKit's fixture loader (`.dat` path only).

## GPU / data needs

None. No GPU code, no runtime I/O of its own — fixture reads go through the injected `Data` root
(`Data.set_datadir!` / `ENV["GRTOOLS_DATADIR"]`). The `:gpu` tolerance class exists for *other*
packages' CPU-vs-GPU gates (e.g. [Kernels](kernels.md)); TestKit itself never touches a device.

## Provenance

TestKit is the one core unit with **no function-level legacy provenance**: `@approx_oracle` had zero
usages in the legacy corpus, so the harness is 100% greenfield. What it absorbs from the
predecessor codes is *patterns*, not code:

| legacy source (all ~20 packages) | absorbed as | note |
|---|---|---|
| inline `rtol`/`atol` asserts (~1924 tolerance literals, ~1754 `isapprox`) | `@approx_oracle` + tolerance classes | pattern only — the literals are NOT transcribed; legacy runtests stay authoritative |
| scattered per-package test-fixture loaders | `load_fixture` over the `Data` spine | kills the verified `ref-code` run-output path escapes |
| ad-hoc "does the port match?" spot checks | `characterization_snapshot` / `characterization_gate` | full-precision, ULP-bounded, perturbation-detecting |
| (new, no legacy analogue) | `assert_dep_budget` / `assert_not_loaded` | CI hygiene gates, wired into all 16 packages |

Grounding: self-tests in `test/` (infra, not physics) — the macro passes/fails correctly and records
its tolerance, the class ordering (`:machine` strictly tighter than `:physics`), the
can-only-tighten guard catching a loosened tolerance, a fixture round-trip through `Data`, the gate
detecting 1-ULP and gross perturbations, and the hygiene gates on TestKit's own `Project.toml`.

## What changed from the predecessor codes

- **New package.** Nothing in the legacy corpus corresponds to TestKit; migrating users lose nothing
  and gain the frozen-tolerance log, the relocation gate, and the hygiene gates.
- **Tolerance discipline is now enforced, not conventional**: new tests must declare a class or a
  frozen tolerance, and CI can prove (via `oracle_records()` + `assert_can_only_tighten`) that no
  tolerance was loosened over a relocation.
- **Dependency budgets are now executable** (`assert_dep_budget`/`assert_not_loaded`), replacing
  prose poison-axis declarations with per-package CI gates — including the invariant that no
  `[sources]` entry may point into a predecessor-code tree (`packages/` or `ref-code/`).

## Public API

```@autodocs
Modules = [TestKit]
Order = [:module, :type, :constant, :macro, :function]
```

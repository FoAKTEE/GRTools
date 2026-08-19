# Contributing to GrTools

## Repository layout

GrTools is a monorepo of 18 independent Julia packages plus an umbrella:

- `core/` — physics-free numerics (NumericsCore, BigSpecial, Data, TestKit,
  Wigner, Kernels, EllipticCore). Dependency budgets are deliberate and tight:
  NumericsCore depends on LinearAlgebra only, and every other core package
  declares its allowed dependency set in its `Project.toml` header comment.
- `capabilities/` — physics packages built over core (Spacetime, Harmonics,
  Perturbation, MetricReconstruction, SelfForce, ScalarSelfForce,
  PostNewtonian, ModifiedGravity, EMRIWaveforms, Hybrid, AccretionDisc).
- Root `Project.toml` / `src/GrTools.jl` — the umbrella package.
- `docs/` — one Documenter site over the whole suite.
- `examples/` — runnable, self-contained example scripts.

Packages reference each other with relative-path `[sources]` entries, so a
plain clone resolves with no registry set-up.

## Working on a package

```bash
julia --project=capabilities/Perturbation
julia> using Pkg; Pkg.instantiate(); Pkg.test()
```

Each package's tests must pass standalone. Cross-package test dependencies go
through `[sources]` path entries in that package's `Project.toml`.

## Design rules

1. **Type-generic numerics.** Every numerical path is generic in the element
   type `T`; BigFloat inputs must inherit the active
   `setprecision(BigFloat, bits)` end-to-end. Never silently downcast — if a
   fast Float64 tier exists, route BigFloat callers to a reference path.
2. **Dependency budgets are binding.** Several packages enforce their allowed
   dependency set with a `test/hygiene.jl` gate; keep it green. Heavy
   dependencies (Symbolics, HDF5, GPU stacks) are confined to the packages
   that already carry them or to weak extensions (`ext/`).
3. **Tests are grounded.** New numerics land with characterization tests
   against a stated oracle (a paper value, a baked upstream fixture under
   `test/fixtures/`, or an exact identity), with tolerances via TestKit's
   tolerance classes. Fixtures record their capture provenance in a header
   comment; never depend on the upstream source live at test time.
4. **Conventions.** Geometric units `G = c = 1`, `M = 1` unless stated;
   Boyer–Lindquist coordinates; Kinnersley tetrad; fields `∝ e^{i(mφ − ωt)}`;
   Teukolsky eigenvalue `λ = E − 2mγ + γ²` with `γ = aω`. State any deviation
   at the call site and in the docs page.

## Documentation

Docs build locally with:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Every exported symbol needs a docstring (`checkdocs = :exports`); each package
has one page under `docs/src/core/` or `docs/src/capabilities/`.

## Commit messages

Commit titles follow `<type>(<scope>)[!]: <imperative summary>` with
`type ∈ {feat, fix, perf, refactor, docs, test, build, ci, style, revert, diag,
exp, chore, infra, notes}`. The body, when present, is a markdown list of typed
objects — `- why:`, `- finding:`, `- change:`, `- run:`, `- result:`,
`- verify:`, `- caveat:`, `- next:`, `- files:` — one kind per line, with every
`finding`/`result` carrying a claim tag `[SOLID|PRELIMINARY|HOLE|FUTURE]`.
Footers are git trailers (`Refs:`, `Verify:`, `Claim:`, `BREAKING CHANGE:`).

The grammar is enforced by the tracked hook [`.githooks/commit-msg`](.githooks/commit-msg).
Activate it once per clone:

```bash
git config core.hooksPath .githooks
git config commit.template .gitmessage
```

The hook **rejects** a commit whose title fails the grammar, and **rejects** any
message carrying AI/assistant attribution — `Co-Authored-By` trailers of any
kind, or a mention of Claude / Anthropic / "generated with" / the robot emoji.
Body and footer issues are warnings; `COMMIT_GATE_STRICT=1` escalates them to
errors. Run it standalone with `.githooks/commit-msg <message-file>`
(exit 0 = admit, 1 = reject).

## Pull requests

- One logical change per PR; include tests.
- Run the affected package's full suite (`Pkg.test`) before pushing, and the
  umbrella smoke test (`julia --project=. -e 'using Pkg; Pkg.test()'`) if you
  touched package layout or `[sources]`.

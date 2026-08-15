# SelfForce

Tier-5 capability of the GrTools monorepo (design §4.5): the **1st + 2nd-order gravitational
self-force (GSF) corpus**. It is HDF5 HOME #1 of the capability layer — the place where the h1
corpus (the first-order Lorenz-gauge mode data) is generated, persisted, and re-ingested.

## Physics this capability owns

- **Lorenz-gauge first-order field** — the Barack–Lousto–Sago (BLS) mode-by-mode Lorenz-gauge
  solution of the linearised Einstein equations for a point mass on a Kerr/Schwarzschild orbit:
  R1..R10 sources, in/out boundary-condition series, integration of the radial field equations,
  gauge-field assembly, and writing the h1 corpus.
- **Second-order Ricci source** — the 10 quadratic combinations `r2_1..r2_10` of the first-order
  field, their 3-harmonic angular coupling, the analytic mass/spin pieces, and the assembled
  second-order Ricci / effective source.
- **Punctures / effective source** — the static `h0` field, the singular-field puncture, the
  α-integration (terminating `2F1`), the β=π/2 mode mixing, the regularisation sums, and the
  retarded/singular subtraction giving the regular self-field on the worldline.

## Legacy packages absorbed

| Legacy package      | Lands in submodule      |
|---------------------|-------------------------|
| `H1Lorenz`          | `SelfForce.LorenzField` |
| `SecondOrderRicci`  | `SelfForce.SecondOrder` |
| `Punctures`         | `SelfForce.Puncture`    |

The radial RHS (`even_ode!`/`odd_ode!`/`monopole_ode!`/`odd_dipole_ode!`) is NOT owned here; its
canonical home is `MetricReconstruction.Lorenz.LorenzRadial`, which `LorenzField` consumes.

## Public API

Re-exported through the submodule namespaces:

- `SelfForce.LorenzField` — Lorenz-gauge solver + `Data`-backed h1-corpus writer.
- `SelfForce.SecondOrder` — `r2_1..r2_10`, `second_order_ricci`, `ricci_source`, `write_src`.
- `SelfForce.Puncture` — `h0` elliptic field, puncture / effective-source, worldline self-field.

(Concrete symbols land incrementally with each TDD module.)

## Dependencies

Core: `NumericsCore`, `BigSpecial`, `Wigner`, `Kernels`, `Data`. Capability: `Harmonics`,
`MetricReconstruction`. External: `HDF5` (the h1 corpus). Standard: `LinearAlgebra`, `Printf`.

All primitives are type-generic `{T}` so BigFloat precision (256..2048-bit) is inherited.

## Grounding (test characterization)

The legacy `H1Lorenz` / `Punctures` / `SecondOrderRicci` packages are HEAVY and are **never** live
test-deps (TEST-DEP WEIGHT RULE). Expected values are CAPTURED one-shot and baked as literal
fixtures: the compiled g++/GSL `r2_oracle.cpp` + `h1.cc` oracles (compiled once to capture, as the
legacy did), wolframscript / Mathematica for the puncture, analytic identities, and a round-trip
h1-corpus write/read (which fires the `Data` HDF5 weakdep extension). The test cold-compile stays
fast against the GrTools deps only.

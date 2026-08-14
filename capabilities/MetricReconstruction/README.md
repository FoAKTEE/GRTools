# MetricReconstruction

Tier-4 GrTool capability (design §4.4): reconstruction of the linear metric perturbation
`h_{μν}` of the Kerr spacetime from perturbed-field scalars. It fuses two legacy packages into
disjoint submodule namespaces so the gauge-dependent concept-name overlaps (`weyl_scalar`,
`to_BLS`) coexist without rename.

## Physics this capability owns

- **Radiation-gauge CCK reconstruction** (`MetricReconstruction.Radiation`): the
  Chrzanowski–Cohen–Kegeles / Hertz-potential procedure. From the Teukolsky ψ4 amplitudes it
  builds the Hertz potential, regenerates the `weyl_scalar`/`psi4` Weyl scalars, the 10
  reconstructed metric components, and the `reconstruct_kerr` / `to_BLS_kerr` / `closure_kerr`
  assembly.
- **Lorenz-gauge reconstruction** (`MetricReconstruction.Lorenz`): the Dolan-2024 a≠0
  Lorenz-gauge metric reconstruction for a circular Kerr orbit (angular/radial closure, component
  assembly, App-H jump conditions), plus `Lorenz.LorenzRadial` — the canonical home of the a=0
  first-order Lorenz-gauge radial field-equation RHS (`even_ode!`/`odd_ode!`/`monopole_ode!`/
  `odd_dipole_ode!`) that `SelfForce.LorenzField` and `AccretionDisc.Force` consume but do not own.

## Legacy packages absorbed

| Legacy package          | Absorbed into          |
|-------------------------|------------------------|
| `MetricReconstruction`  | `Radiation` submodule  |
| `KerrLorenzCirc`        | `Lorenz` submodule     |

KLC's hand-rolled Clebsch–Gordan routine is promoted to `Wigner.clebsch_gordan`.

## Public API

Exposed through the submodule namespaces to keep the disambiguated concept names:

- `MetricReconstruction.Radiation`: Hertz potential, `weyl_scalar`/`psi4`, the 10 metric
  components, `reconstruct_kerr`, `to_BLS_kerr`, `closure_kerr`.
- `MetricReconstruction.Lorenz`: Dolan-2024 reconstruction entrypoints, `weyl_scalar`/`to_BLS`
  (Lorenz variants), and `Lorenz.LorenzRadial.{even_ode!, odd_ode!, monopole_ode!, odd_dipole_ode!}`.

Concrete symbols are exported as each TDD module lands.

## Dependencies

- GrTool: `NumericsCore`, `Wigner`, `Spacetime`, `Harmonics`, `Perturbation`.
- External: `QuadGK`, plus the `LinearAlgebra`/`Printf` standard libraries.

## Precision

Every primitive is type-generic `{T}`, so BigFloat precision (256..2048-bit) is inherited from the
input element type.

## Testing / grounding

Tests are written first and grounded by characterization against the LIGHT live legacy oracles
`MetricReconstruction` and `KerrLorenzCirc` (verified runtime-HDF5-free; KLC's HDF5 is test-only)
at their frozen tolerances, supplemented by analytic identities (CCK reconstruction vs the
Teukolsky ψ4, metric-compatibility, Dolan App-H jump conditions). No heavy package is a live
test dependency.

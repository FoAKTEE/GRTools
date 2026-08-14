# Perturbation

Tier-3 vacuum/test perturbed-field solver for the Kerr spacetime — the **Teukolsky-radial
engine** and the **light-dep firewall** of GrTool (design §4.3).

## Physics this capability owns

The homogeneous solutions of the radial Teukolsky equation and everything built directly on
them for a *vacuum / test* (point-particle) perturbation:

- **Homogeneous radial solutions** In/Up via the MST (Mano–Suzuki–Takasugi) hypergeometric
  representation, plus a static ω=0 branch, confluent-Heun, a vendored numerical (DP5) backend,
  and the Sasaki–Nakamura form — all behind one dispatcher.
- **Renormalized angular momentum** ν (the MST continued-fraction root) and the MST series
  coefficients, generated through the BigFloat 256..2048-bit tier.
- **Point-particle amplitudes** Z/C and the **energy / angular-momentum / Carter fluxes**
  (circular / eccentric / spherical / generic orbits × spin weights).
- **MST scattering phase** + Starobinsky constant + B_refl/B_inc reflection/transmission.
- (Later round) spinning-secondary **spin-fluxes**.

## Legacy packages absorbed

Functions relocate to this physics home from (design §6 crosswalk rows dest=Perturbation):

| Legacy package | What moves in |
| --- | --- |
| `Teukolsky` | homogeneous radial In/Up, RenormalizedAngularMomentum ν, MST coefficients, amplitudes Z/C, point-particle fluxes |
| `ScatteringPhaseMST` | `starobinsky_constant`, `_brefl_over_binc`, `scattering_phase_numerical`, `mst_coefficient`, `scattering_phase_table` |
| `TeukolskySpinFluxes` (later round) | the 13 spinning-secondary `TeukolskySpinMode*` flux entrypoints |

## The firewall

This capability is **HDF5-free, JLD2-free, Symbolics-free, and external-ODE-free**. There are
zero `HDF5 / JLD2 / Symbolics / OrdinaryDiffEq / CUDA / AMDGPU / KernelAbstractions` edges in its
runtime closure (a closure-load gate test enforces this). "External-ODE-free" rests on the
**vendored DP5 + Sasaki–Nakamura** radial integrators — these must never be refactored onto
OrdinaryDiffEq. The only heavy dependency reaching a user is the
`HypergeometricFunctions`/`SpecialFunctions` closure carried transitively through `BigSpecial`.

## Precision

Every primitive is type-generic `{T}`, so BigFloat precision (256/512/1024/2048 bits) is
inherited. The BigFloat MST/continued-fraction generation tier stays on CPU (BigFloat-CPU /
Float64-GPU boundary); heavy outer flux sweeps scale via `Distributed` pmap + CPU Threads ≤20
under the BLAS-guard.

## Dependencies

- **GrTool:** `NumericsCore` (continued-fraction, interpolation, conventions, precision),
  `BigSpecial` (Γ / 2F1 / pochhammer / Kummer), `Spacetime` (geodesic constants, frequencies,
  roots), `Harmonics` (SWSH eigenvalue + S).
- **stdlib:** `LinearAlgebra`, `Printf`.
- **external:** none.

## Test-characterization oracle

Grounded test-only against the LIGHT live oracles `Teukolsky` (closure
KerrGeodesics/SWSH/SpecialFunctions/HypGeo — all light) and `ScatteringPhaseMST` (deps
Teukolsky/SWSH) at their frozen tolerances, plus analytic identities (ν continued-fraction root,
MST 2F1 connection, flux Newtonian-quadrupole 32/5 limit, Starobinsky–Chandrasekhar scattering
identity). No heavy package is ever a live test-dep.

## Public API

Built TDD across the BUILD rounds. The `Radial` submodule (the core everything downstream needs)
lands first; `AmplitudesFluxes`, `Scattering`, and `SpinFluxes` follow.

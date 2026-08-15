# AccretionDisc

Tier-6 EMRI-accretion-disc capability — the **top of the capability DAG**;
nothing depends on it. The Lindblad-torque / spiral-density-wave disc front-end.

## Physics owned

The fluid response of a thin, circular-equatorial accretion disc around a Kerr black hole to an
embedded EMRI secondary, reduced via the κ²-dispersion (Lin-Shu / Goldreich-Tremaine) density-wave
approximation to a one-dimensional radial master equation, and the resulting Lindblad / corotation
resonance torque.

### Submodules

- **`Modes`** — the relativistic disc-mode solver (from legacy `DiscModes`):
  - `kappar2` — the radial epicyclic frequency κ_r² from the **Spacetime geodesic shear** (a
    closed-form geometric quantity, independent of the metric perturbation 𝒲);
  - the height-averaged κ²-dispersion reduced operator `D = κ² − ω̃²`;
  - `B_m` mode functions and the `𝔉_m` forcing source;
  - the Lindblad / corotation resonance torque.

- **`Force`** — the circular-equatorial Kerr thin-disc torque front-end (from legacy
  `DiscEMRIForce`): the spinning-BH fluid solve / torque balance. **Consumes**
  `MetricReconstruction.Lorenz.closure_kerr` / `reconstruct_kerr` (NOT `SelfForce`, NOT
  `Kernels`) plus the a=0 Lorenz corpus. `HDF5` is required by **this submodule only** (captured
  disc-mode data, injected `datadir`).

## Provenance

| Legacy | Absorbed as |
|--------|-------------|
| `DiscModes` | `AccretionDisc.Modes` |
| `DiscEMRIForce` | `AccretionDisc.Force` |

## Dependencies

GrTools: `NumericsCore`, `Data`, `Spacetime`, `Harmonics`, `Perturbation`, `MetricReconstruction`.
External: `HDF5` (Force submodule only).

## Grounding

Built tests-first. The disc-mode leaves are grounded **without heavy live deps**: the geodesic
κ_r² closed form, the corotation regular-singular exponent (vs RT-A6 −0.0818, non-circular
validation), Newtonian disc-mode limits, and one-shot fixtures captured from the published
arXiv:2601.19123 disc-mode data and the `DiscModes` oracle. `DiscModes`/`DiscEMRIForce` are NOT
live `[targets].test` deps (they pull an HDF5/Teukolsky closure).

Precision: every primitive is type-generic `{T}` so BigFloat precision (256..2048-bit) is inherited.

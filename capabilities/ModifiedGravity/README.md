# ModifiedGravity

Beyond-GR (modified Teukolsky) source machinery — a thin, independently versioned GrTools research
capability (Tier-4).

## Physics owned

The perturbed-source machinery for an EMRI inspiraling into a Kerr black hole surrounded by matter
(Polcar & Witzany 2025, arXiv:2507.15720). A thin equatorial matter ring back-reacts on the Kerr
background; the reconstructed (0,1) metric assembles an **effective source** so the perturbed Weyl
scalar ψ4 obeys a Teukolsky equation with that source. This capability owns:

- the 12-column `.dat` ring-forcing **source tables** (POD + reader/writer),
- the per-mode **radial Hertz potential** front-end (`(l,m,n,k)`, with the load-bearing
  `kloc = (x==1 ? 0 : k)` convention),
- the closed-form **NP/GHP source kernels** (rational functions at equatorial Schwarzschild
  θ=π/2),
- the assembled **source amplitudes** + the `.dat` generator,
- the perturbed **(1,1) ring dynamics** (geodesic + density/continuity response, NP→BL tetrad map),
- the closed-form **ψ4←Hertz (CCK) reconstruction**.

## Provenance

- `ModifiedTeukolsky` (`source_table_io`, `hertz_potential`, `source_kernels`, `source_amplitudes`,
  `rhs_amplitudes`, `ring_geodesic`, `ring_functions`, `weyl_from_hertz`) → the six submodules below.

## Submodules / public API

| Submodule | Public surface |
|-----------|----------------|
| `SourceTables` | `SourceMode`, `read_source_table`, `write_source_table`, `write_source_mode`, `source_mode_line` |
| `HertzPotential` | `hertz_potential_mode` |
| `SourceKernels` | closed-form `Sourceace*`/`Sourcenorm*` NP/GHP kernels |
| `SourceAmplitudes` | `source_amplitudes`, `rhs_amplitudes` |
| `RingDynamics` | `ring_urho_amplitudes`, `ring_functions` |
| `WeylFromHertz` | `psi4_from_hertz` |

`[HOLE]`: `h44_from_hertz` (full reconstructed (4,4) metric component) is a documented hole — needs
the xAct PDNP machinery not yet ported.

## Dependencies

- **Runtime:** `Perturbation` (Teukolsky radial engine, reused for the Hertz front-end),
  `Harmonics` (SWSH), `NumericsCore` (interpolants / `dataschema` / conventions), `Data` (`.dat`
  table readers), `Printf`.
- **Self-contained:** does NOT depend on `MetricReconstruction` — the CCK ψ4←Hertz identity is
  hard-baked as a closed-form kernel rather than called through Radiation machinery.
- **Precision:** type-generic `{T}`; BigFloat inheritance from Perturbation/Harmonics preserved.
  HDF5-free + Symbolics-free.

## Test grounding

Test-only characterization against the predecessor `ModifiedTeukolsky` code at its frozen
tolerances (its outputs are baked as literal fixtures in the test files), the three shipped `.dat`
source-table oracles (round-trip via `Data`), and the analytic CCK reconstruction identity
(reldiff ~4e-9).

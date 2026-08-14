# ModifiedGravity

Beyond-GR (modified Teukolsky) source machinery — a thin, independently versioned Tier-4 research
capability (design §4.7).

## Overview

The perturbed-source machinery for an EMRI inspiraling into a Kerr black hole surrounded by matter
(Polcar & Witzany 2025, arXiv:2507.15720). A thin equatorial matter ring back-reacts on the Kerr
background; the reconstructed ``(0,1)`` metric assembles an **effective source** so the perturbed
Weyl scalar ``\psi_4`` obeys a Teukolsky equation with that source. This capability owns:

- the 12-column `.dat` ring-forcing **source tables** (`SourceMode` POD + reader/writer),
- the per-mode **radial Hertz potential** front-end over `(l,m,n,k)`, with the load-bearing
  `kloc = (x==1 ? 0 : k)` convention driving both the amplitude sign ``(-1)^{l+k_{\rm loc}+m}``
  and the minus-mode conjugation sign ``\sigma = (-1)^{l+k_{\rm loc}}``,
- the closed-form **NP/GHP source kernels** (rational functions at equatorial Schwarzschild
  ``\theta=\pi/2``, ported verbatim from the saved notebook output cells),
- the assembled per-mode **source amplitudes** and the summed-over-``l`` `.dat` generator
  (`rhs_amplitudes` — the exact pipeline that produced the shipped tables),
- the perturbed **(1,1) ring dynamics** (geodesic + density/continuity response, NP→BL tetrad map),
- native **GHP operators** ``\text{þ}, \text{þ}'`` on equatorial Schwarzschild (jet-calculus based),
- the closed-form ``\psi_4\leftarrow`` Hertz (CCK) reconstruction **and** the finite-``r``
  reconstructed ``(4,4)=\bar m\bar m`` metric component `h44_from_hertz` (the legacy `[HOLE]`,
  now filled natively — no xAct).

The Hertz front-end covers circular-equatorial (`e=0, x=1`) and spherical (`e=0, |x|<1`) orbits;
eccentric orbits throw (no migrated grounding yet).

## Submodule map

`ModifiedGravity` is a **flat namespace**: the design-§4.7 "submodules" are source files included
into the single top-level module (verified — `src/` declares no nested `module` blocks). The map
below is file → public surface:

| Source file (design-§4.7 area) | Public surface |
|---|---|
| `source_tables.jl` (SourceTables) | `SourceMode`, `read_source_table` (both re-exported from `Data`), `source_mode_line`, `write_source_mode`, `write_source_table` |
| `hertz_potential.jl` (HertzPotential) | `hertz_potential_mode` |
| `source_kernels.jl` (SourceKernels) | *internal*: 10 closed-form `_source_kernel_{l,n,m,mdag,norm}_{plus,minus}` NP/GHP kernels; `yv`-callable, so they carry no Harmonics dependency |
| `source_amplitudes.jl` (SourceAmplitudes) | `source_amplitudes` |
| `rhs_amplitudes.jl` (SourceAmplitudes) | `rhs_amplitudes` |
| `ring_dynamics.jl` (RingDynamics) | `ring_urho_amplitudes`, `ring_functions` |
| `ghp_operators.jl` (GHPOperators) | `schw_spincoeff`, `thorn`, `thorn_prime` |
| `weyl_from_hertz.jl` (WeylFromHertz) | `psi4_from_hertz`, `h44_from_hertz` |

The read side of the table I/O (`SourceMode` + `read_source_table`) was promoted into core
`Data` (shared with the FEW flux/amplitude `.dat` readers) and is re-exported here — their
docstrings render on the [Data](../core/data.md) page. This capability adds the write side,
reusing `Data._mma_sci` for byte-exact ``\omega`` reproduction.

## Dependencies and poison axis

Runtime `[deps]` (Project.toml): `Perturbation` (Teukolsky radial engine, reused by the Hertz
front-end), `Harmonics` (SWSH values at ``\pi/2``), `NumericsCore`, `Data` (`.dat` reader +
`_mma_sci`), `MetricReconstruction` (see below), `Printf`.

**Sanctioned MetricReconstruction edge (Rev-3).** The design-§3.2 dep line for ModifiedGravity was
`Perturbation, Harmonics, NumericsCore (+Printf, Data)`; `MetricReconstruction` is a **sanctioned
Rev-3 amendment** for the h44 GHP fill: `ghp_operators.jl` reuses
`MetricReconstruction.Lorenz`'s jet calculus (`jmul`/`jadd`/`jderiv`/…, the single source of truth
for radial-derivative Taylor towers) and `h44_from_hertz` reuses
`MetricReconstruction.Radiation.radial_jet` (the ``s=+2`` Teukolsky ODE tower). The CCK
``\psi_4\leftarrow`` Hertz identity itself remains hard-baked as a closed-form kernel
(`_psi4_hertz_kernel`) — it does not call Radiation machinery.

**Poison axis** — the binding §7.7 CI budget in
`capabilities/ModifiedGravity/test/hygiene.jl` (`TestKit.assert_dep_budget` +
closure-load gate `TestKit.assert_not_loaded`):

- *may carry*: `Printf`, `NumericsCore`, `Harmonics`, `Perturbation`, `MetricReconstruction`, `Data`;
- *may NOT carry (declared or loaded)*: `HDF5`, `JLD2`, `Symbolics`, `OrdinaryDiffEq*`, `CUDA`,
  `AMDGPU`, `KernelAbstractions`, `Kernels`, `WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase`.

The package is HDF5-free and Symbolics-free by CI gate, not by convention.

## GPU / data needs

- **GPU:** none. Parallelism is CPU threads/processes over the `(l,m,n,k)` mode axis (design
  §"per-capability parallelism").
- **Data:** the only data surface is the 12-column tab-separated `.dat` ring-forcing source tables
  (three shipped oracles in `test/data/`, e.g. `i1.0_p12.0_e0_r20.0.dat`), read through `Data`.
  No HDF5, no grids.
- **Precision:** the numeric surface is type-generic; BigFloat inheritance from
  Perturbation/Harmonics is preserved (see "What changed vs legacy").

## Legacy provenance (drawsFrom)

Design §6.7 crosswalk — 1:1 with legacy `ModifiedTeukolsky`, zero cross-legacy pulls:

| Legacy pkg / file | New home |
|---|---|
| ModifiedTeukolsky / `source_table_io.jl` | `source_tables.jl` (write side); read side + `SourceMode` → `Data` (re-exported) |
| ModifiedTeukolsky / `source_table_io.jl` / `_mma_sci` | → `Data` (reused here for ``\omega`` formatting) |
| ModifiedTeukolsky / `hertz_potential.jl` | `hertz_potential.jl` |
| ModifiedTeukolsky / `source_kernels.jl` | `source_kernels.jl` |
| ModifiedTeukolsky / `source_amplitudes.jl`, `rhs_amplitudes.jl` | `source_amplitudes.jl`, `rhs_amplitudes.jl` |
| ModifiedTeukolsky / `ring_geodesic.jl`, `ring_functions.jl` | `ring_dynamics.jl` |
| ModifiedTeukolsky / `weyl_from_hertz.jl` | `weyl_from_hertz.jl` |
| — (legacy `[HOLE]`, xAct-only upstream) | `ghp_operators.jl` + `h44_from_hertz` (new, native) |

Test grounding: frozen `ModifiedTeukolsky` characterization literals (no live legacy dep — every
expected value is baked), the three shipped `.dat` oracles (round-trip via `Data`), the analytic
CCK ``\psi_4`` reconstruction identity (reldiff ``\sim 4\times10^{-9}``), standalone GHP algebra
identities for ``\text{þ}/\text{þ}'``, and — for `h44` — the necessary asymptotic
``h_{44} \to 2\psi_4/\omega^2`` plus a finite-``r`` match to the independent wolframscript operator
oracle (`oracle/h44_from_hertz_oracle.wl`).

## What changed vs legacy

- **BigFloat leak fixed (design hole H-MG-BIGFLOAT).** Legacy `hertz_potential_mode` coerced
  `float(a)`/`float(ω)` and pinned `ComplexF64`, silently downcasting extended-precision inputs.
  The port threads the element type through, so BigFloat flows into
  `TeukolskyRadial`/`TeukolskyPointParticleMode*` unmolested. The load-bearing
  `kloc = (x==1 ? 0 : k)` convention and the ``(-1)^{l+k_{\rm loc}+m}`` sign are verbatim.
- **Table read side de-duplicated.** `SourceMode` + `read_source_table` moved to core `Data`
  (shared with FEW `.dat` readers); this package re-exports them and owns only the write side.
- **Kernels decoupled from Harmonics.** The 10 source kernels take a caller-supplied
  `yv(s,m) = {}_sY_{l,m}(\pi/2)` callable; `source_amplitudes` wires in `Harmonics`.
- **h44 hole closed, xAct-free.** `h44_from_hertz` assembles the reconstructed
  ``\bar m\bar m`` component natively from the recovered GHP operator (SC l.368):
  ``-(-1)^m\big[(3\varepsilon'+3\varepsilon'^{\dagger}+2\rho')\,\text{þ}'\Psi + \Psi B +
  \text{þ}'\text{þ}'\Psi\big]`` with ``\Psi = \mathrm{RadHZ}^{\dagger}/\xi^4`` — this is what
  motivates the sanctioned Rev-3 `MetricReconstruction` dependency (jet calculus + `radial_jet`).
- **Dead legacy deps dropped** (grep-verified in the design round): `KerrGeodesics`,
  `LinearAlgebra`, `DelimitedFiles` are gone from `[deps]`.

## Public API

```@autodocs
Modules = [ModifiedGravity]
Order = [:module, :type, :constant, :function]
```

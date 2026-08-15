# AccretionDisc

Tier-6 EMRI-accretion-disc capability (design §4.10) — the **top of the GrTools capability
DAG**; nothing depends on it. The Lindblad-torque / spiral-density-wave disc front-end.

## Physics owned

The fluid response of a thin, circular-equatorial accretion disc around a Kerr black hole to an
embedded EMRI secondary, reduced via the ``\kappa^2``-dispersion (Lin–Shu /
Goldreich–Tremaine) density-wave approximation to a one-dimensional radial master equation, and
the resulting Lindblad / corotation resonance torque.

The physicist's map of the pipeline:

1. **Disc background** — the radial epicyclic frequency ``\kappa_r^2`` comes from the
   `Spacetime` geodesic shear as a closed-form geometric quantity, *independent* of the metric
   perturbation ``\mathcal{W}``; circular-orbit ``\Omega_K``, sound speed, surface density and
   ISCO complete the relativistic (Dyson–D'Orazio) and Newtonian (Miranda–Rafikov) backgrounds.
2. **Mode reduction** — Gaussian height-averaging of the off-equatorial (``z=\cos\theta``)
   background to ``O(z^2)`` collapses the linearised fluid equations onto the height-averaged
   reduced operator ``D = \kappa^2 - \tilde\omega^2``, giving the RT-A1 radial master-PDE /
   reduced-ODE coefficients ``A_m``, ``B_m`` and the forcing source ``\mathfrak{F}_m``.
3. **Resonances and torque** — Wronskian plus Lindblad / corotation residuals locate the
   resonances (including the corotation regular-singular exponent); the Lindblad torque follows.
4. **EMRI forcing** — the embedded secondary's Lorenz-gauge metric perturbation
   ``h_{\mu\nu}`` (a=0 corpus, or the a≠0 Dolan-2024 Kerr reconstruction) is assembled into the
   fluid-frame effective force, projected onto the ``\mathcal{Q}``-source, and drives the forced
   master-ODE solve; torque-balance (F8) and flux cross-checks close the loop.

Every primitive is type-generic `{T}` so BigFloat precision (256–2048 bits) is inherited.

## Submodule map

Verified module paths: `AccretionDisc`, `AccretionDisc.Modes`, `AccretionDisc.Force`.

### `Modes` — relativistic disc-mode solver (from legacy `DiscModes`)

Runtime-HDF5-free.

| File | Contents |
|---|---|
| `Background.jl` | `RelativisticDisc`, `disc_OmegaK`, `kappar2` (geodesic-shear closed form), `delta_r`, `csound`, `eden0`, `risco_disc` |
| `Velocity.jl` | RT-A2 velocity-reconstruction operators `qinv_operator`, `q_operator` |
| `OffEquator.jl` | off-equatorial background to ``O(z^2)`` (`disc_OmegaK_z`, `alpha_z`, `Uvec_z`, `pressure_z`, `eden_z`), off-equator operators (`qinv_operator_z`, `q_operator_z`), the PMLA Gaussian `heightavg` |
| `MasterPDE.jl` | `master_radial_coeffs` — the RT-A1 radial master-PDE coefficients |
| `Resonances.jl` | `lindblad_residual`, `corotation_residual`, `wronskian`, `lindblad_locations` |
| `ReducedODE.jl` | `A_m`, `B_m_corotation`, `B_m_wkb`, `B_m_full`, `corotation_exponent`, `dispersion`, `reduced_operator` (``D = \kappa^2 - \tilde\omega^2``) |
| `Forcing.jl` | `lindblad_torque`, `forcing_Fm` (``\mathfrak{F}_m``), `laplace_b` |

### `Force` — circular-equatorial Kerr thin-disc torque front-end (from legacy `DiscEMRIForce`)

The full design-§4.10 pipeline, including the S5-5 feature-gap port
(Source / Fluid / Pipeline / KerrFrontend):

| File | Contents |
|---|---|
| `Loader.jl` | `H1Mode` (in-memory, distinct from the on-disk `Data.H1Schema`), `load_h1_mode`, `load_modes`, `available_modes` — the a=0 Lorenz BLS corpus loader over `Data`'s h1schema (injected `datadir`) |
| `Assemble.jl` | parity-aware tetrad→BLS assembly: `kerr_parity`, `tetmodes_from_bls`, `hmunu`, `hmunu_eq` |
| `Source.jl` | fluid-frame effective force → ``\mathcal{Q}``-projected source → radial forcing ``\mathfrak{F}(r)``: `disc_a0`, `disc`, `force_mode`, `Vvec`, `source_F`, plus the wing-precision `source_F_grid` array pass |
| `Fluid.jl` | fluid reconstruction ``(\delta U, \delta p, \delta\varepsilon)``: `FluidVars`, `reconstruct`, and the advected angular-momentum flux `Jadv` |
| `Pipeline.jl` | forced master-ODE solve: `ForcingInterp`, `build_forcing`, `forced_solve`; the F8 torque-balance / flux cross-checks `f8_torque_balance`, `f8_torque_balance_ilr`, `f8_flux_crosscheck`; the F5 source-structure driver `f5_source` |
| `KerrFrontend.jl` | the a≠0 metric source: `kerr_mode`, `kerr_modes` consume `MetricReconstruction.Lorenz.{closure_kerr, smatrices_kerr, reconstruct_kerr, to_BLS_kerr}` and pack the 10 Barack–Lousto–Sago ``\bar h^i(r)`` into the same `H1Mode{T}` container the a=0 corpus loader returns — so the downstream (spin-agnostic) Source/Fluid/Pipeline stages run unchanged |

The Kerr front-end **replaces the metric source only**: per ``(l,m)`` the ``(l+m)``-parity
selects the physical BLS fields (even ⇒ the 7 even fields, odd ⇒ the 3 odd fields, in loader
order), matching the corpus convention.

## Dependencies and poison-axis declaration

Declared `[deps]` (Project.toml): `LinearAlgebra`, `Printf`, and the GrTools packages
`NumericsCore`, `Data`, `Spacetime`, `Harmonics`, `Perturbation`, `MetricReconstruction`, plus
`HDF5` (**direct** — required by the `Force` submodule only, as the trigger that activates
`Data`'s HDF5 weakdep extension; design §5.4 Rev-2 extension-activation rule).

The binding CI dep budget lives in `test/hygiene.jl` (§7.7 gates 1 + 2 + 5):

- **may carry** (allow): `LinearAlgebra`, `Printf`, `NumericsCore`, `Data`, `Spacetime`,
  `Harmonics`, `Perturbation`, `MetricReconstruction`, `HDF5`;
- **must declare**: `HDF5` (the Data HDF5-ext trigger, gate 5);
- **may NOT carry** (forbid): `SelfForce`, `Kernels`, `Symbolics`, `JLD2`, `OrdinaryDiffEq*`,
  `CUDA`, `AMDGPU`, `KernelAbstractions`, `WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase`,
  `BigSpecial`, `Wigner`;
- a closure-load gate asserts after `using AccretionDisc` that none of `SelfForce`, `Kernels`,
  `Symbolics`, `JLD2`, `OrdinaryDiffEq(Verner)`, `CUDA`, `AMDGPU`, `KernelAbstractions`,
  `WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase` were loaded.

**§7.3 rule (L3 leak fix, binding).** The a=0 / a≠0 Lorenz-gauge radial physics is reached
*only* through `MetricReconstruction.Lorenz` (the canonical, HDF5-free, Kernels-free home
`MetricReconstruction.Lorenz.LorenzRadial` for the a=0 radial ODE systems) — **never** through
`SelfForce` and **never** through `Kernels`. A `SelfForce` or `Kernels` edge is the exact leak
§7.3 forbids: it would drag `KernelAbstractions` + the HDF5 writer + `Wigner` onto the disc
front-end. Because AccretionDisc is top-of-DAG, the one accepted cost — HDF5 reaching a
pure-disc-mode user through the fused capability — propagates to no other capability, and
`Modes` remains cleanly separable (zero `Force` deps, zero HDF5).

## GPU / data needs

- **GPU: none.** CPU-only; no `Kernels`/`KernelAbstractions`/`CUDA` anywhere in the closure
  (asserted by the hygiene gate).
- **`Modes`: no runtime data.** All backgrounds are closed-form profiles.
- **`Force`: HDF5-backed.** Reads the a=0 Lorenz-gauge ``h^{(1)}`` BLS corpus via `Data`'s
  h1schema — the sole HDF5 entry point — with an injected `datadir` (the legacy `H1DIR`
  environment path is gone). The test suite carries its own small captured fixtures
  (`test/fixtures/h1-l{2..6}m2.h5`).
- **Threading.** Per-``l`` / per-grid-point sweeps (the `kerr_mode(s)` reconstruction and
  BLS evaluation) are embarrassingly parallel: race-free preallocated slots under the
  `NumericsCore.blas_single_thread` guard, threaded == serial bit-identical (design §8). Pass
  `threaded=false` to `kerr_mode` when the caller already threads the outer ``l`` axis.

## Legacy provenance (drawsFrom, design §6.10)

| Legacy package / files / functions | New home |
|---|---|
| `DiscModes` — `newtonian_disc.jl`, `master_equation_cool_newt.jl`, `solver_newt.jl`, `heightaverage.jl`, `relativistic_disc.jl`, `off_equator.jl`, relativistic master-PDE, `lindblad_torque.jl`, `metric_recon.jl` (`bl_angular_eq`) | `AccretionDisc.Modes` |
| `DiscEMRIForce` — `loader.jl`, `assemble.jl`, `kerr_source.jl`, `source.jl`, `fluid.jl`, `pipeline.jl` | `AccretionDisc.Force` |
| `DiscModes` — `{kerr_geometry, metric_recon, relativistic_disc}.jl` geometry + `_risco` | → `Spacetime.Geometry` |
| `DiscModes` — `metric_recon.jl` / `Ylm_real`, `Ylm_eq`, `assoc_legendre` | → `Harmonics` |
| `DiscModes` — `lorenz_metric.jl` / `even_ode!`, `odd_ode!`, `monopole_ode!`, `odd_dipole_ode!` | → `MetricReconstruction.Lorenz.LorenzRadial` (canonical) |
| `DiscModes` — `lorenz_metric.jl` / `r_of_rs` | → `Spacetime.Geometry` |
| `DiscModes` — `integrator.jl` DP5; `DiscEMRIForce` — `numerics.jl` fornberg/cheb/bary | → `NumericsCore` (ODE / interpolation) |
| `DiscEMRIForce` — `loader.jl` `read_h1_mode` path handling | → `Data` h1schema (injected `datadir`) |

## What changed vs legacy

- **Fusion, with hoists.** `DiscModes` + `DiscEMRIForce` became the two namespaced submodules;
  their geometry (`kerr_metric`/`christoffel`/`_risco` forks), spherical harmonics, a=0 Lorenz
  radial ODE systems, vendored DP5 integrator and interpolation primitives moved to their
  canonical GrTools homes (table above) instead of being copied in.
- **NEW a≠0 front-end (S5-5).** `kerr_mode` / `kerr_modes` (`KerrFrontend.jl`) are the ported
  spinning-BH path: `MetricReconstruction.Lorenz` Dolan-2024 reconstruction → BLS
  ``\bar h^i(r)`` → `H1Mode{T}`, sharing one banded `closure_kerr` solve across ``l``. The
  legacy code reached this physics through a private Lorenz fork; the port consumes the single
  canonical home (§7.3).
- **Off-equator / height-average (S5-5).** `OffEquator.jl` ports the PMLA ``O(z^2)``
  off-equatorial background + Gaussian height-average; the connection now comes from the
  analytic `Spacetime.christoffel` (same ``\Gamma^\lambda{}_{\mu\nu}`` convention as the legacy
  finite-difference `kerr_christoffel`, collapsed). The legacy `F6_NOT5` env experiment knob is
  *not* ported — the validated default (t5 pressure term kept) is hardwired.
- **Data path.** The legacy `H1DIR=../../../ref-code/...` environment escape is replaced by
  `Data.datadir` injection / Artifacts resolution.
- **Test grounding without heavy live deps.** `DiscModes`/`DiscEMRIForce` are *not* live
  `[targets].test` deps (TEST-DEP WEIGHT RULE): ``\kappa_r^2`` is grounded by the geodesic
  closed form, the corotation exponent against RT-A6 (−0.0818), Newtonian limits, and one-shot
  fixtures captured from the arXiv:2601.19123 disc-mode data and the legacy oracles.
- **Dropped:** the vestigial `FastGaussQuadrature` dep (declared by legacy `DiscModes`, never
  used in src).
- **Precision:** every primitive is type-generic `{T}`, so BigFloat (256–2048-bit) is inherited
  end-to-end.

## Public API

```@autodocs
Modules = [AccretionDisc, AccretionDisc.Modes, AccretionDisc.Force]
Order = [:module, :type, :constant, :function]
```

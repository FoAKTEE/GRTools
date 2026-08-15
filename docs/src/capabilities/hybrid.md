# Hybrid

Tier-5 GrTools capability: the **WaSABI self-force / post-Newtonian hybrid**
EMRI waveform engine.

## Overview

Hybrid owns the SF-PN-hybrid inspiral-and-waveform pipeline. A self-force-accurate forcing
model is fused with a post-Newtonian completion to drive a spinning EMRI inspiral, and the
resulting trajectory is mapped to a binary gravitational waveform. This is the WaSABI-C
scheme: a 9-variable SF-PN-hybrid forcing (``F_\omega/F_\nu/F_m/F_\chi``) with
``\delta x``-coupled trajectory evolution and closed-form shift channels.

The physicist mental map, bottom-up:

1. **Interpolants over on-disk SF tables.** Wolfram-faithful interpolation (moving-window
   order-8 1D, tensor-product Hermite/Lagrange 2D/3D) turns the legacy `.m`/`.wl`
   self-force flux / local-invariant / amplitude grids into callable forcing channels.
2. **Forcing models.** Per-order adiabatic/post-adiabatic inspiral models —
   `0PASchwarz`, `1PAT1`, `1PAT1ea`/`1PAT1R` (spinning secondary, 7-variable state
   ``(\Omega, \varphi, \nu, M, \tilde\chi_1, \tilde\chi_2, \delta m)``), `0PAKerr` — plus
   their amplitude counterparts, and the WaSABI-C 9-variable SF-PN-hybrid forcing with its
   machine-generated ``\delta x`` closed forms.
3. **Evolution.** A `Vern9` ODE driver integrates ``d\Omega/dt = F/M``,
   ``d\varphi/dt = \Omega/M`` (N-var states) to the ``\Omega_{\mathrm{crit}}`` /
   ``r_0``-stop termination, single-trajectory or threaded batch.
4. **Waveform.** The trajectory is assembled into strain modes
   ``h_{\ell m}(t) = A_{\ell m}(\Omega(t))\, e^{-i m \varphi(t)}`` via the model registry.

Hybrid is a **top-of-DAG leaf**: nothing depends up on it, so the heavy
`OrdinaryDiffEqVerner`/`SciMLBase` precompile floor is confined here and never reaches the
core numeric/GPU waveform paths. It is also deliberately **outside the Symbolics
quarantine**: the legacy WaSABI source has zero Symbolics tokens — its
digamma was hand-rolled — so Hybrid is firewalled from `PostNewtonian` and never pays the
Symbolics load floor.

## Submodule map

| Submodule | Contents |
|---|---|
| `Hybrid.Interpolation` | WL-faithful interpolation foundation: `Interp1`/`interp1_order8` (order-8 moving-window), `Interp2Hermite`, `Interp3Lagrange`, analytic `derivative`/`derivative2`/`partial`, `interp_logexp`, and the Kerr grid-abscissa remap `kerr_grid_coords`. Bottom of the intra-package DAG. |
| `Hybrid.BindingEnergy` | First-law binding-energy / flux closed forms over the redshift interpolant: `E_FLx`, `E_dm`, `dEdΩ` and derivatives. Pure-algebra leaf consumed by the 1PA forcing. |
| `Hybrid.Models` | Inspiral forcing models (`Forcing0PASchwarz`, `Forcing1PAT1`, `Forcing1PAT1ea`/`Model1PAT1R`, `Forcing0PAKerr`) with their ODE right-hand sides (`inspiral_rhs!`, `kerr_inspiral_rhs!`, `_spin_eom!`), stop conditions (`stop_r0`, `omega_crit`), and amplitude models (`Amp0PASchwarz`, `Amp1PAT1`, `Amp1PAT1ea`, `Amp0PAKerr`). |
| `Hybrid.HybridC` | The WaSABI-C engine: `WaSABICModel` interpolant bundle, the machine-generated ``\delta x``/``d\delta x/d\chi``/forcing closed forms (`wc_δx`, `wc_forcing`, …), and the ``\delta x``-coupled trajectory evolution (`WaSABICSolution`). |
| `Hybrid.Drivers` | N-var inspiral evolution over `OrdinaryDiffEqVerner`/`SciMLBase`: `int_inspiral`, threaded `int_inspiral_batch`, `InspiralSolution`/`SpinInspiralSolution`, and the `*_at` trajectory evaluators. |
| `Hybrid.Registry` | The upstream 7-plugin name→constructor catalogue: `list_models`, `get_inspiral`, `get_amplitudes`. |
| `Hybrid.BinaryWaveform` | `BinaryInspiralModel`, `binary_inspiral`, `waveform`/`waveform_mode`/`mode_list`, `trajectory`. |
| `Hybrid.Introspection` | `cite_model` (citation metadata), `list_modes`, `forcing_terms`, `get_inspiral_equations`/`inspiral_evaluate`. |

Layering is acyclic: `Interpolation` feeds `BindingEnergy`, both feed the `Models`
forcing/amplitude leaves; `Drivers`, `HybridC`, `Registry`, `BinaryWaveform` and
`Introspection` sit on top.

## Dependencies and poison-axis budget

Declared `[deps]` (Project.toml):

- **GrTools:** `NumericsCore` (the canonical interpolation suite + `blas_single_thread`),
  `BigSpecial` (`PolyGamma` digamma), `Spacetime` (`kerr_isco`, `kepler_Ω`), `Data`
  (`.wl`/`.m` table readers).
- **External:** `OrdinaryDiffEqVerner`, `SciMLBase`, `JLD2`, plus stdlib
  `Printf`/`LinearAlgebra`.

The CI hygiene gate (`test/hygiene.jl` — hygiene gates 1, 2 and 5, enforced by
`TestKit.assert_dep_budget` / `assert_not_loaded`) pins this budget exactly:

- **May carry:** only the nine packages above (the `allow` list); the SciML floor is
  legitimate here because Hybrid is a top-of-DAG leaf.
- **Must declare:** `JLD2` — Hybrid is the *extension trigger* for `Data`'s JLD2
  weakdep-extension (`DataJLD2Ext`). A Julia package extension activates only when the
  loading package declares the trigger in its own `[deps]`, so if Hybrid dropped the
  direct `JLD2` entry, `Data`'s JLD2 reader would silently never load.
- **May not carry (poison axes):** `HDF5` (Hybrid touches `Data` only for `.m`/`.wl`
  artifacts; the HDF5 ext must stay dormant), `Symbolics` (the PN quarantine must not
  leak in), `WignerSymbols`, the GPU axis (`Kernels`, `KernelAbstractions`, `CUDA`,
  `AMDGPU`), the `OrdinaryDiffEq` monolith (Verner-only), `Roots`, `NLsolve`.
  Gate 2 additionally asserts none of these appear in the loaded-module set after
  `using Hybrid`.

The test target adds only `Test`, `Data` and the test-only `TestKit`; the legacy `WaSABI`
package is **never** a live test dependency — tests ground against wolframscript runs of
the WL sources, baked fixtures from the legacy `.m`/`.wl` data
(`test/fixtures/wasabi/`), and analytic identities.

## GPU / data needs

- **No GPU.** Hybrid is Float64 CPU-only; the GPU/Kernels axis is a forbidden dep (gate
  above). Parallelism is CPU threads over the trajectory axis (`int_inspiral_batch`,
  with BLAS pinned single-threaded via `NumericsCore.blas_single_thread`).
- **Data.** The forcing/amplitude models are built from on-disk SF tables
  (`InspiralModels/sf_data`, `AmplitudeModels/sf_amp_data`: 1SF/2SF fluxes, local
  invariants, leading-``S_2`` flux, stop conditions) read through `Data`'s WL
  InputForm readers (`read_wl_table`, `read_wl_grid2d`, `read_wl_grid3d`, …). JLD2 IO
  goes through `Data`'s JLD2 extension, triggered by Hybrid's direct `JLD2` dep.
  Citation metadata ships as package assets (`assets/MetadataModels/*.m`).

## Provenance

Hybrid absorbs exactly one legacy package — `WaSABI` — in full:

| Legacy source (WaSABI) | Lands in |
|---|---|
| `inspiral_{0PASchwarz,1PAT1,1PAT1ea,0PAKerr}.jl`, `amplitude_*.jl`, `interpolation.jl` (`interp_logexp`, `kerr_grid_coords`) | `Hybrid.Models` (interp foundation re-exported from `Hybrid.Interpolation`) |
| `inspiral_WaSABI-C{,_deps,_codegen}.jl` | `Hybrid.HybridC` |
| `binding_energy.jl` | `Hybrid.BindingEnergy` |
| `inspiral_driver.jl`, `inspiral.jl` | `Hybrid.Drivers` |
| `registry.jl` | `Hybrid.Registry` |
| `binary_inspiral.jl`, `waveform.jl` | `Hybrid.BinaryWaveform` |
| `utilities.jl` | `Hybrid.Introspection` |

## What changed from the predecessor codes

- **ODE stack swap.** The legacy `OrdinaryDiffEq` monolith dep was replaced by
  `OrdinaryDiffEqVerner` (`Vern9`-only) + `SciMLBase`; the monolith is a forbidden dep.
- **Digamma collapse.** WaSABI's hand-rolled digamma/polygamma fork was collapsed onto
  the canonical `BigSpecial.PolyGamma` (re-exported by `HybridC`).
- **WL-interp collapse onto NumericsCore.** The moving-window order-8
  interpolant, `Interp2Hermite`/`Interp3Lagrange` and their analytic
  derivatives/`partial`s now live in `NumericsCore` (`Interp1WL`,
  `interp1_movingwindow`, …); `Hybrid.Interpolation` keeps only the historical
  Hybrid-facing names plus the two composites NumericsCore does not own
  (`interp_logexp`, `kerr_grid_coords`). Type-generic `{T}`: BigFloat tables and queries
  are inherited, never hard-Float64'd.
- **JLD2 as a direct ext trigger.** `import JLD2` at the top of the module
  is load-bearing: it fires `Data`'s `DataJLD2Ext`. It is asserted by hygiene gate 5 and
  by a dedicated liveness test (`test/test_data_jld2_ext.jl`).
- **Upstream stubs filled.** The upstream declaration-only
  `GetInspiralEquations`/`InspiralEvaluate` now have functional bodies
  (`get_inspiral_equations`, `inspiral_evaluate`), assembled from the already-grounded
  `Models`/`HybridC` forcing primitives.
- **Data-path injection.** Model constructors take an explicit `sf_datadir` instead of
  the legacy hard-wired run-output paths.

## Public API

```@autodocs
Modules = [Hybrid, Hybrid.Interpolation, Hybrid.BindingEnergy, Hybrid.Models,
           Hybrid.Drivers, Hybrid.HybridC, Hybrid.Registry, Hybrid.BinaryWaveform,
           Hybrid.Introspection]
Order = [:module, :type, :constant, :function]
```

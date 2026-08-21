# AccretionDisc

Tier-6 EMRI-accretion-disc capability — the **top of the GrTools capability
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

Verified module paths: `AccretionDisc`, `AccretionDisc.Modes`, `AccretionDisc.Force`,
`AccretionDisc.SpiralWaves`.

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
| `PressurelessLindblad.jl` | the CLOSED-FORM pressureless Kerr Lindblad torque of arXiv:2510.03564 (the comparison baseline of arXiv:2608.17003 Fig. 9): `hkr_shear_argument` (``X``), `hkr_b0`/`hkr_b1`/`hkr_b2`, `hkr_lindblad_power_quadratic`/`_linear` (Eqs. 18a/18b), `hkr_lindblad_torque_rate` (Eq. 20), `hkr_lindblad_torque_mode`/`hkr_lindblad_torque` (2608 Eqs. 77/78), `hkr_mode_cutoff`/`hkr_mode_cutoff_powerlaw`/`HKR_M45_COEFFICIENT` (2608 Eqs. 79/80, 2510 Eq. 23), `hkr_torque_reversal_criterion`/`hkr_torque_reversal_radius` (``\mathcal{G}_1``, Eq. 25), and the `PowerLawSigma` surface-density profile. Distinct from the `Forcing.jl` Hirata GW-overlap kernel and the `ThinDisc.jl` mode-sum: different model, different evidence chain, `hkr_`-prefixed names. ``K_{0,1}`` comes from `SpecialFunctions` and is injectable (no BigFloat `besselk`). |

### `SpiralWaves` — relativistic spiral density waves and disc-EMRI torques (arXiv:2608.17003)

The paper's linear Lagrangian perturbation theory on the Page-Thorne equatorial line element,
end to end: disc equilibrium, adaptive corotation/WKB grid, the first-order master system (47),
the EMRI singular field (52)-(60), the implicit mode solver, the angular-momentum diagnostics
(72a)/(72b), the surveys, the image synthesis and the fully Newtonian limit pipeline.

| File | Contents |
|---|---|
| `background.jl` | `DiscBackground`, `rotation_profile` (exact hydrostatic balance (11)), `Omega`, `omega_d`, power laws (12), `gamma1` (22), `r_isco`, `F_J0`/`F_J1` |
| `grid.jl` | `build_grid` — the App.-C adaptive radial grid: `corotation_radius` (61), sinh corotation block, WKB wavenumber, smooth-min cell sizing |
| `master.jl` | `MasterInputs`, `master_inputs`, `master_matrices` (47), `master_matrices_newtonian` (50b-g), `delta_u_phi`, `landau_shift` (62)-(63), `dOmega_dr`/`dv_dr` |
| `master_derived.jl` | **GENERATED** — `_master_entries`, the entries of C and S. Emitted by route **R3** (see below) |
| `r2_derived.jl` | **GENERATED** — `_r2_delta_u_phi`, `_r2_amf` (72a), `_r2_torque_density` (72b). Same route |
| `source.jl` | the EMRI singular field (52)-(55), softening (59) and its COVARIANT cosine Fourier modes (60) (`singular_field_mode`, `source_modes`; `raise_modes` for the contravariant ones) |
| `solver.jl` | `solve_mode` — the implicit-midpoint global solve (64) with outgoing-wave / WKB-impedance boundary conditions |
| `diagnostics.jl` | `amf` (72a), `torque_density` (72b), `cumulative_torque`, `total_torque` (74), `lindblad_torque` (75), `balance_residual` (71), `m_cut` |
| `survey.jl` | threaded m-sweeps and configuration surveys |
| `image_synthesis.jl` | mode-sum reconstruction of the wake (Figs. 1-2) |
| `newtonian.jl` | the same scheme on a purely Newtonian background — the ``T_{\rm Newt}`` of Fig. 5 |

#### Route R3 — the master system is DERIVED, not transcribed

`master_derived.jl` and `r2_derived.jl` are generated by
`capabilities/AccretionDisc/derivation/` — a **standalone `Symbolics.jl` project** (its own `Project.toml`, **not** a dependency of `AccretionDisc`; the package dependency
budget is unchanged and `Symbolics` is on the package's forbidden-load list) that derives the master
system from the paper's own equations: line element (4), fluid equations (3), linearisation
(15)-(22), Fourier ansatz (45)-(46) and the reduction of M:604-608. It certifies every step either
symbolically or by exact `Rational{BigInt}` evaluation (residual exactly `0//1`), emits the two files
with a deterministic printer, and re-reads the emitted text to check it interprets exactly like the
symbolic expression.

**No external CAS is in the build path.** The authors' supplementary *Mathematica* notebook (route
R1) and an independent *Mathematica* derivation (route R2) are kept as ORACLES in
`test/SpiralWaves/fixtures/oracle/` and gate the derived sources forever
(`test/SpiralWaves/test_derived_equivalence.jl`).

**On-shell only.** The three routes agree wherever the background satisfies both algebraic
constraints — ``W = 1/\sqrt{1-v^2}`` (8) and the hydrostatic balance (11) — and carry different
off-shell continuations, so every route comparison is made on-shell. Measured: `3.5e-74`
(BigFloat 256) at the 40 projected fixture points, `2.6e-12` (Float64) against the 50-point physical
Mathematica fixture — where the notebook form itself only reaches `7.0e-10`, i.e. the derived form is
also better conditioned.

### `Force` — circular-equatorial Kerr thin-disc torque front-end (from legacy `DiscEMRIForce`)

The full pipeline, including the feature-gap port
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
`Data`'s HDF5 weakdep extension — the extension-activation rule).

The binding CI dep budget lives in `test/hygiene.jl` (hygiene gates 1, 2 and 5):

- **may carry** (allow): `LinearAlgebra`, `Printf`, `NumericsCore`, `Data`, `Spacetime`,
  `Harmonics`, `Perturbation`, `MetricReconstruction`, `HDF5`;
- **must declare**: `HDF5` (the Data HDF5-ext trigger, gate 5);
- **may NOT carry** (forbid): `SelfForce`, `Kernels`, `Symbolics`, `JLD2`, `OrdinaryDiffEq*`,
  `CUDA`, `AMDGPU`, `KernelAbstractions`, `WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase`,
  `BigSpecial`, `Wigner`;
- a closure-load gate asserts after `using AccretionDisc` that none of `SelfForce`, `Kernels`,
  `Symbolics`, `JLD2`, `OrdinaryDiffEq(Verner)`, `CUDA`, `AMDGPU`, `KernelAbstractions`,
  `WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase` were loaded.

**Canonical-home rule (binding).** The a=0 / a≠0 Lorenz-gauge radial physics is reached
*only* through `MetricReconstruction.Lorenz` (the canonical, HDF5-free, Kernels-free home
`MetricReconstruction.Lorenz.LorenzRadial` for the a=0 radial ODE systems) — **never** through
`SelfForce` and **never** through `Kernels`. A `SelfForce` or `Kernels` edge is exactly the leak
this rule forbids: it would drag `KernelAbstractions` + the HDF5 writer + `Wigner` onto the disc
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
  `NumericsCore.blas_single_thread` guard, threaded == serial bit-identical. Pass
  `threaded=false` to `kerr_mode` when the caller already threads the outer ``l`` axis.

## Provenance

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

## What changed from the predecessor codes

- **Fusion, with hoists.** `DiscModes` + `DiscEMRIForce` became the two namespaced submodules;
  their geometry (`kerr_metric`/`christoffel`/`_risco` forks), spherical harmonics, a=0 Lorenz
  radial ODE systems, vendored DP5 integrator and interpolation primitives moved to their
  canonical GrTools homes (table above) instead of being copied in.
- **NEW a≠0 front-end.** `kerr_mode` / `kerr_modes` (`KerrFrontend.jl`) are the ported
  spinning-BH path: `MetricReconstruction.Lorenz` Dolan-2024 reconstruction → BLS
  ``\bar h^i(r)`` → `H1Mode{T}`, sharing one banded `closure_kerr` solve across ``l``. The
  legacy code reached this physics through a private Lorenz fork; the port consumes the single
  canonical home.
- **Off-equator / height-average.** `OffEquator.jl` ports the PMLA ``O(z^2)``
  off-equatorial background + Gaussian height-average; the connection now comes from the
  analytic `Spacetime.christoffel` (same ``\Gamma^\lambda{}_{\mu\nu}`` convention as the legacy
  finite-difference `kerr_christoffel`, collapsed). The legacy `F6_NOT5` env experiment knob is
  *not* ported — the validated default (t5 pressure term kept) is hardwired.
- **Data path.** The legacy `H1DIR` environment-path escape is replaced by
  `Data.datadir` injection / Artifacts resolution.
- **Test grounding without heavy live deps.** `DiscModes`/`DiscEMRIForce` are *not* live
  `[targets].test` deps: ``\kappa_r^2`` is grounded by the geodesic
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

# SelfForce

Tier-5 capability: the **1st + 2nd-order gravitational self-force (GSF) corpus**,
fusing the earlier `H1Lorenz`, `Punctures`, and `SecondOrderRicci` research codes over the shared **h1
corpus** (the first-order Lorenz-gauge mode data, intrinsic HDF5 mediated by `Data.H1Schema`).
This is **HDF5 HOME #1** of the capability layer: the h1 corpus is generated, persisted, and
re-ingested here.

## Physics this capability owns

- **Lorenz-gauge first-order field** — the Barack–Lousto–Sago (BLS) mode-by-mode Lorenz-gauge
  solution of the linearised Einstein equations for a point mass on a circular Schwarzschild
  orbit: the ``R_1..R_{10}`` sources, the in/out boundary-condition series, integration of the
  radial field equations, gauge-field assembly, and the writer of the h1 corpus.
- **Second-order Ricci source** — the 10 quadratic combinations ``r^{(2)}_1..r^{(2)}_{10}`` of the
  first-order field, their 3-harmonic angular coupling, the analytic mass/spin pieces
  (`h1_mass`/`h1_spin`, `d2hdr2`), and the assembled second-order Ricci / effective source
  (`second_order_ricci`, `ricci_source`, `write_src`).
- **Punctures / effective source** — the static ``h^{(0)}`` field (AGM elliptic integrals), the
  singular-field puncture, the α-integration (terminating ``{}_2F_1``), the ``\beta=\pi/2``
  Wigner-D mode mixing, the digamma-tail regularisation sums, and the retarded/singular
  subtraction ``h^R = h^{\mathrm{ret}} - h^S`` giving the regular self-field on the worldline.

The radial field-equation RHS (`even_ode!`/`odd_ode!`/`monopole_ode!`/`odd_dipole_ode!`) is **not
owned here**: its single canonical, physics-justified home is
`MetricReconstruction.Lorenz.LorenzRadial`, which `LorenzField` consumes.

## Submodule map

| Submodule | Legacy origin | Contents |
|---|---|---|
| `SelfForce.LorenzField` | `H1Lorenz` | Orbit energetics (`OrbitalParams`, `orbit_energy`, `orbit_ang_mom`, `orbit_phi_freq`), coupled-set bookkeeping (`CSETS`, `select_coupled_set`, `gamma_to_lm`, `lm_to_gamma`, `mode_list`), the `R1..R10` sources, in/out boundary-condition series (`even_out_bcs` … `even_static_in_bcs`), `integrate_field_equations` + `calculate_scaling_coefficients` on the vendored NumericsCore DP5+GBS integrators, the static monopole, mode drivers (`NModeData`, `compute_mode`, `compute_modes`, `compute_one_mode`), tortoise maps `rs_of_r`/`r_of_rs`, and the h1-corpus writer `write_h1_mode`/`h1_schema`. |
| `SelfForce.SecondOrder` | `SecondOrderRicci` | h1-corpus reader `read_h1` + the `ZeroPadField` zero-pad shim, mode selection `selection_modes`, 3-harmonic coupling (front door to `Wigner.coupling`), analytic `h1_mass`/`h1_spin`/`d2hdr2`, the coupling-evaluating `r2_term` host wrapper over the `Kernels.r2_term` GPU host, and assembly `second_order_ricci` / `ricci_source`(`_mass`/`_spin`) / `write_src`. |
| `SelfForce.Puncture` | `Punctures` | Static `h0_static` (AGM `ellipK`/`ellipE`), worldline Christoffels (`christoffel`, `christoffel_dr`, `christoffel_dtheta`), α-integration (`exactint`, `modeint`, `alphaint`, `alphaintv1..4`, `alphaintt1..2`), regularised l-sums (`reg_sum`/`noreg_sum` → `RegResult`), field-equation second derivative `ddhret`, Martel–Poisson tensor basis (`tensor_harmonic`, `coord_component`, `coord_deriv`), puncture modes `hP`/`dhP` (+ the exact-rational helper `rat`), Wigner-D mode mixing `h_mode`, singular field `hS_mode`, subtraction `compute_hR`/`compute_ddhR`, and the worldline regular self-field `assemble_{hR1,dhR1,d2hR1}Data` + `covariant_first`/`covariant_second`. |

Re-exported-not-re-owned front doors: `a_il` (owner: `NumericsCore`), `coupling` (owner:
`Wigner`), `h1_parity_fields` (owner: `Data`) — single method tables, no duplication.

## Dependencies and poison-axis budget

Declared `[deps]` (Project.toml): core `NumericsCore`, `BigSpecial`, `Wigner`, `Kernels`, `Data`;
capability `Spacetime`, `Harmonics`, `MetricReconstruction`; external `HDF5` (hard),
`ForwardDiff` + `QuadGK` (Puncture); stdlib `LinearAlgebra`, `Printf` (Printf is flagged
vestigial — allowed for a green suite but has no `src/` usage).

The binding budget is enforced by the CI gates in `SelfForce/test/hygiene.jl`
(`TestKit.assert_dep_budget` + `TestKit.assert_not_loaded`):

- **Must declare** — `HDF5` **directly**. This is the extension trigger: `Data`'s
  HDF5 reader/writer bodies live in a weakdep extension that activates only when the *loading*
  package declares the trigger in its own `[deps]`. Omitting it would silently disable corpus IO.
- **May carry** — exactly the allow-list above; any other `[deps]` entry fails the gate.
- **Must NOT carry (poison axes)** — `Symbolics` (quarantined in PostNewtonian),
  `OrdinaryDiffEq*` (LorenzField uses the vendored NumericsCore DP5+GBS integrators — a
  CI-gated invariant), `WignerSymbols`, `JLD2`, `CUDA`/`AMDGPU`, direct `KernelAbstractions` (must arrive
  only through `Kernels`), `Roots`, `NLsolve`, `SciMLBase`. A closure-load gate re-asserts after
  `using SelfForce` that none of these modules appears in the loaded-module set.

## GPU and data needs

- **GPU** is confined to `SecondOrder`: the ``r^{(2)}`` quadratic forms run on the
  `Kernels.r2_term` KernelAbstractions host (one work-item per output axis, serial reduction —
  bit-identical to serial CPU). **Float64 on device only**; every CPU primitive in the package is
  type-generic `{T}` and inherits BigFloat working precision (256–2048 bit). `LorenzField` and
  `Puncture` are CPU-only.
- **Data**: the h1 corpus (per-mode `h1-l{l}m{m}.h5` files) is intrinsic — `LorenzField` *writes*
  it (via `Data.write_mode_h5` against the versioned `H1Schema` contract) and `SecondOrder` /
  `Puncture` *re-ingest* it (`read_h1`, `read_hret_mode`). All corpus locations are injected
  (`datadir=`); the legacy `ref-code/h1Lorenz/run_r8/fields` relative-path escape is gone.
- **Tests** run against baked fixtures only — captured one-shot from the compiled g++/GSL
  `r2_oracle.cpp`/`h1.cc` oracles, wolframscript/Mathematica (Puncture), analytic identities, and
  a round-trip corpus write/read (which fires the Data HDF5 extension). The heavy predecessor codes
  are never live test-deps.

## Provenance

| Predecessor code / files | Lands in |
|---|---|
| `H1Lorenz` — `orbit, sources, boundary_conditions, integrator, scaling, static, gauge_fields, coupled_sets, assemble, driver, io_hdf5` | `SelfForce.LorenzField` |
| `Punctures` — `static, h0_coefficients, alpha_integration, puncture_modes*, mode_mixing, regsum, tensor_basis, subtraction*, singular_field, worldline_self_field` | `SelfForce.Puncture` |
| `SecondOrderRicci` — `r2_1..r2_10, r2terms, h1_mass, h1_spin, d2hdr2, modes, assembly, read_h1` | `SelfForce.SecondOrder` |

## What changed from the predecessor codes

- **Radial RHS single-homed.** The a=0 Lorenz radial ODE systems moved to
  `MetricReconstruction.Lorenz.LorenzRadial`; `LorenzField` consumes them instead of owning a
  fork (resolves the H1Lorenz-vs-DiscModes duplicate).
- **`r2_term` Kernels-purity hoist.** The 10 quadratic ``r^{(2)}`` combinations now live in
  the pure-arithmetic `Kernels.r2_term` GPU host, which takes a *precomputed* coupling array
  (per `Kernels.R2_COUPLING_SPEC`) and carries no Wigner machinery. The coupling-evaluating
  **host wrapper lives here** as `SelfForce.SecondOrder.r2_term`: it evaluates the 3-harmonic
  Wigner couplings, then launches the kernel — signature, semantics, and output identical to the
  legacy `SecondOrderRicci.r2_term`.
- **Vendored integrators.** OrdinaryDiffEq is gone; `LorenzField` integrates on the NumericsCore
  vendored DP5+GBS steppers (a CI-gated invariant).
- **Type-generic `{T}`.** Legacy Float64-only paths are generic; BigFloat precision propagates
  end-to-end (only the on-device r2 path is Float64-capped).
- **IO through `Data`.** The HDF5 reader/writer bodies live in `Data`'s weakdep extension behind
  the versioned `H1Schema`; the `SecondOrder.read_h1` corpus reader intentionally keeps its own
  6-column (`∂²_r`-bearing) unpacking, which is a different on-disk contract from the 4-part
  H1Lorenz writer format that `Data.read_h1` serves.
- **Special functions centralised.** Puncture's `_poch`/`_f21_term` → `BigSpecial.pochhammer` /
  `f21_terminating`; digamma tails → `BigSpecial.digamma_big`; Wigner-D fixed-angle rotation →
  `Wigner.wignerD_fixed`; elliptic K/E consumed from `NumericsCore`.

## Public API

```@autodocs
Modules = [SelfForce, SelfForce.LorenzField, SelfForce.SecondOrder, SelfForce.Puncture]
Order = [:module, :type, :constant, :function]
```

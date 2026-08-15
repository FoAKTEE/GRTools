# EMRIWaveforms

Tier-5 capability (design §4.8, §6.8): the **EMRI trajectory + observed-waveform pipeline**,
fusing the two legacy waveform engines **FastEMRIWaveforms** (FEW) and **KerrBackwards** (KB).
This is **HDF5 home #2**: inspiral ODE → flux/amplitude grid models → mode selection → mode
summation → observed-waveform generation.

## Overview

End-to-end extreme-mass-ratio-inspiral waveform generation, in the order a physicist thinks
about it:

1. **Trajectory** — adiabatic inspiral ODE evolution of the orbital elements
   ``(p, e, Y, \Phi_\varphi, \Phi_\theta, \Phi_r)`` driven by the analytic 5PN-``e^{10}`` fluxes
   (vendored NumericsCore DOPR853 integrator, `Spacetime` separatrix-refined stopping,
   Roots/NLsolve parameter inversion). One driver: `run_inspiral(...; backend=:few|:kb)`.
2. **FluxModels** — grid-interpolated flux ODE right-hand sides: closed-form PN normalization
   plus a residual interpolated over the on-disk Teukolsky flux grids loaded through `Data`
   (Schwarzschild-eccentric and Kerr-eccentric-equatorial).
3. **AmplitudeModels** — Teukolsky mode amplitudes ``A_{\ell m n}`` from the ROMAN dense
   neural network (forward pass = `Kernels.roman_forward` GEMM) and the bicubic / clamped-B-spline
   amplitude-grid interpolants. One struct: `RomanAmplitude` +
   `get_amplitudes(...; backend=:few|:kb)`.
4. **ModeSelect** — power-thresholded selection of the contributing ``(\ell, m, n)`` harmonics
   (keep modes until cumulative power reaches ``1-\epsilon`` of the total, ``\pm m`` partner
   pairing).
5. **ModeSum** — waveform assembly
   ``h = \sum_{\ell m n} A_{\ell m n}\, {}_{-2}Y_{\ell m}\, e^{-i\Phi_{mn}}``
   in the time domain (direct / trajectory-interpolated), the frequency domain (SPA), and the
   AAK approximation — all summation kernels live in `Kernels`.
6. **AAK** — augmented-analytic-kludge basis: AK parameter maps, Peters–Mathews Bessel
   coefficients, and the harmonic summation. Single `aak_summation` / `aak_waveform` API.
7. **Mappings** — orbital-parameter coordinate maps (`p_to_y`, Kerr-ecc-eq forward/backward,
   PN constants of motion, ``Y \leftrightarrow x_I``, the ``(\dot E, \dot L_z) \to (\dot p, \dot e)``
   Jacobian), reusing the `NumericsCore` conventions.
8. **Generate** — top-level observed-waveform driver: source frame → SSB/detector frame
   (viewing angles, polarization angle, ``2\psi`` rotation) via `generate_emri_waveform`.

## Submodule map

| Submodule | Contents | Backend structure |
|---|---|---|
| `Trajectory` (+ nested `Pn5Fluxes`) | 5PN-``e^{10}`` flux table, `pn5_rhs`, `ODEModel` interface, DOPR853 inspiral driver | single API, `backend=:few\|:kb` selector |
| `FluxModels` | grid-residual flux RHS | `FluxModels.KerrBackwards` (Schwarz-ecc `GslBicubic` over `.dat` grid) / `FluxModels.FEW` (Kerr-ecc-eq `TricubicNotAKnot` over `KerrEccEqFluxData.h5`) — genuinely distinct physics, kept as qualified forks |
| `AmplitudeModels` (+ nested `KerrAmpSpline`) | `RomanAmplitude` NN, `AmpInterpSchwarzEcc` bicubic, `AmpInterpKerrEqEcc` over the `ZNAmps` B-spline grid | single API, `get_amplitudes(...; backend=:few\|:kb)` |
| `ModeSelect` | power-threshold harmonic selection | `ModeSelect.KerrBackwards` (canonical ``(\ell,m,n)``) / `ModeSelect.FEW` (adds polar index ``k`` + `"all"`/explicit-list branches) |
| `ModeSum` | TD direct/interpolated + FD-SPA reduction | `ModeSum.KerrBackwards` / `ModeSum.FEW` — thin convention adapters over the single `Kernels` kernels |
| `AAK` (+ nested `ParameterMapAAK`) | AK parameter map, `ak_bessel_coeffs`, AAK summation/waveform | single API, `backend=:few\|:kb` (no-op selector; forks were byte-identical) |
| `Mappings` | coordinate maps + flux Jacobian | single |
| `Generate` | frame algebra + `generate_emri_waveform` dispatcher | single (wired: `:FastSchwarzschildEccentricFlux` source-frame, `:Pn5AAKWaveform` detector-frame) |

The remaining `KerrBackwards.*` / `FEW.*` deep namespaces are reachable **only qualified** —
they are separate physics backends (different grids, different index sets), not duplicated
scaffolding, and are excluded from any flat re-export.

## Dependencies and poison-axis declaration

GrTools deps: `NumericsCore`, `Spacetime`, `Harmonics`, `Kernels`, `Data`.
External deps: `HDF5` (**direct** — see below), `OrdinaryDiffEqVerner`, `Roots`, `NLsolve`,
`SpecialFunctions` (`besselj` for `ak_bessel_coeffs`), plus stdlib `Printf` / `LinearAlgebra` /
`Random`.

The binding dep budget is CI-enforced by `test/hygiene.jl` (§7.7 gates 1+2+5), which

- **allows** exactly: `Printf`, `LinearAlgebra`, `Random`, `NumericsCore`, `Spacetime`,
  `Harmonics`, `Kernels`, `Data`, `HDF5`, `OrdinaryDiffEqVerner`, `Roots`, `NLsolve`,
  `SpecialFunctions`;
- **requires** `HDF5` to be declared directly (gate 5): all grid-reader *bodies* live in `Data`'s
  HDF5 weakdep-extension, and Julia extension semantics activate `DataHDF5Ext` only when the
  consumer declares the trigger package in its own `[deps]`;
- **forbids** (may never carry): `Symbolics`, `WignerSymbols`, `JLD2`, `CUDA`, `AMDGPU`,
  `KernelAbstractions` (must arrive only *through* `Kernels`), the monolithic `OrdinaryDiffEq`
  (only the Verner sub-package is sanctioned, §7.5), `SciMLBase`, `BigSpecial`;
- asserts post-load (gate 2) that `Symbolics` / `WignerSymbols` / `JLD2` / `CUDA` / `AMDGPU` /
  `OrdinaryDiffEq` are not in the loaded-module closure.

The capability is **Symbolics-clean**: FEW's four import edges onto
Teukolsky / TeukolskyPN / TeukolskySpinFluxes / SpinningOrbit were grep-verified dead
(0 qualified call sites) and deleted, removing the entire heavy Teukolsky subtree — and with it
the only Symbolics carrier — from the closure.

No heavy live test deps: FastEMRIWaveforms / KerrBackwards are **never** live test oracles;
grounding is by analytic identity and captured fixtures baked from the checksum-verified
`.fewdata` grids.

## GPU / data needs

- **GPU** — Float64-on-device via `Kernels` only: TD mode sum, FD-SPA, AAK summation, ROMAN
  GEMM (one work-item per (mode, sample)). Trajectory ODE and root/separatrix solves stay on
  CPU. CPU fallback is threaded (≤ 20) with the BLAS single-thread guard.
- **Precision** — Float64/GPU-only. This capability *calls* the BigFloat-generic lower tiers
  (Harmonics eigenvalues, Perturbation, MetricReconstruction) at Float64 and imposes no
  precision constraint upward (§9).
- **Data** — all grid IO goes through `Data` (`read_fluxgrid`, `read_ampgrid`, `read_roman_nn`,
  `read_kerrampgrid`, `read_dat`); no EMRIWaveforms source performs direct HDF5 reads. Grids:
  the Schwarzschild-eccentric flux/amplitude grids, the ROMAN NN weights
  (`SchwarzschildEccentricInput.hdf5`), `KerrEccEqFluxData.h5`, and the 4.74 GB Kerr-ecc-eq
  `ZNAmps` amplitude grid (`[FUTURE]`-large, §10).

## Legacy provenance (drawsFrom, design §6.8)

| Legacy source | New home |
|---|---|
| KB `trajectory.jl`, `pn5fluxes.jl`; FEW `trajectory/{emri_inspiral,dopr853,ode_base,pn5_rhs,pn5}.jl` | `EMRIWaveforms.Trajectory` (single `backend=:few\|:kb` API) |
| FEW `trajectory/{flux,kerrflux}.jl`; KB `{schwarzeccflux,fastschwarz_waveform,slowschwarz_waveform}.jl` | `EMRIWaveforms.FluxModels` (KB/FEW forks) |
| FEW `amplitude/{ampinterp2d,romannet}.jl`; KB `{romanamplitude,interp2damplitude}.jl` | `EMRIWaveforms.AmplitudeModels` (single API) |
| FEW `modeselector.jl`; KB `modeselector.jl` | `EMRIWaveforms.ModeSelect` (KB/FEW forks) |
| FEW `summation/{directmodesum,interpolatedmodesum,fdinterpolatedmodesum}.jl`; KB `{summation,fdsummation}.jl` | `EMRIWaveforms.ModeSum` (adapters; kernels → `Kernels`) |
| FEW `summation/aaksummation.jl`; KB `{summation,bessel_ak,parametermap_aak}.jl` | `EMRIWaveforms.AAK` (single API; summation → `Kernels.aak_summation`) |
| FEW `mappings/{kerrecceq,pn,jacobian,schwarzecc}.jl`, `geodesic_elq.jl`; KB `inclination.jl` | `EMRIWaveforms.Mappings` |
| FEW `waveform/generate.jl`; KB `waveform.jl` | `EMRIWaveforms.Generate` |
| FEW `geodesic.jl`, `separatrix.jl`; KB `utility.jl` geometry | → `Spacetime` (migrated out) |
| both forks: interpolants / `GslBicubic` / `ModeIndexTable` / constants / Carlson elliptic | → `NumericsCore` (hoisted, single copy) |
| FEW imports of Teukolsky / TeukolskySpinFluxes / SpinningOrbit / TeukolskyPN | **deleted** (dead, 0 call sites) |

## What changed vs legacy

- **Single-API backend collapse (§4.8/§6.8 milestone, CLOSED for three of six collision
  survivors).** The three fork pairs that were the *same* decision — the `Trajectory` inspiral
  driver, `AmplitudeModels.RomanAmplitude`, and `AAK.aak_summation`/`aak_waveform` — are each
  ONE canonical API with a `backend::Symbol` selector (`:few` default, `:kb`) that reproduces
  both legacy forks bit-for-bit; the redundant fork scaffolding is deleted. `FluxModels`,
  `ModeSelect`, and `ModeSum` retain qualified `KerrBackwards`/`FEW` submodules because those
  are genuinely distinct backends (Schwarz-ecc vs Kerr-ecc-eq flux grids; ``(\ell,m,n)`` vs
  ``(\ell,m,k,n)`` selection; TD vs FD-SPA assembly).
- **FD-SPA moved to `Kernels` (S4-a).** The frequency-domain stationary-phase mode sum
  (`fd_interpolated_mode_sum` / `fd_flip_average`) was hoisted to
  `Kernels.fd_interpolated_mode_sum` (single copy, bit-identical to the pre-hoist local kernel
  per the baked fixture); the FEW names in `ModeSum.FEW` are thin facades.
- **HDF5 reader bodies moved to `Data` (S4-b).** `read_roman_nn`, `read_ampgrid`,
  `read_fluxgrid`, `read_kerrampgrid`, `kerr_mode_table` live in `DataHDF5Ext`;
  EMRIWaveforms' direct `HDF5` dep exists solely as the extension trigger.
- **`ak_bessel_coeffs` ownership (S4-c purity hoist).** The Peters–Mathews Bessel basis moved
  verbatim *out of* `Kernels` into `EMRIWaveforms.AAK` (its canonical home), and is injected
  into `Kernels.aak_summation` as the `bessel_coeffs` provider — `Kernels` carries no
  `SpecialFunctions` machinery (§5.5/§7.7 purity gate); `SpecialFunctions` was promoted to a
  regular dep here.
- **Geometry migrated out; primitives hoisted.** Separatrix/frequencies/ISCO live in
  `Spacetime`; the ``{}_sY_{\ell m}`` carrier in `Harmonics`; cubic/bicubic/tricubic
  interpolants, `ModeIndexTable`, constants, Carlson elliptic, `leaky_relu` in
  `NumericsCore`/`Kernels` (single copies).
- **Known holes.** KB `amplitude.jl` `mode_amplitudes` remains a NotYetPorted stub;
  `SlowSchwarzschildEccentricFlux` is `[FUTURE]`; the full Python end-to-end FEW match is a
  documented `[HOLE]` (every component is grounded individually).

## Public API

### Top-level module

```@autodocs
Modules = [EMRIWaveforms]
Order = [:module, :type, :constant, :function]
```

### Trajectory

```@autodocs
Modules = [EMRIWaveforms.Trajectory, EMRIWaveforms.Trajectory.Pn5Fluxes]
Order = [:module, :type, :constant, :function]
```

### Flux models

```@autodocs
Modules = [EMRIWaveforms.FluxModels, EMRIWaveforms.FluxModels.KerrBackwards, EMRIWaveforms.FluxModels.FEW]
Order = [:module, :type, :constant, :function]
```

### Amplitude models

```@autodocs
Modules = [EMRIWaveforms.AmplitudeModels, EMRIWaveforms.AmplitudeModels.KerrAmpSpline]
Order = [:module, :type, :constant, :function]
```

### Mode selection

```@autodocs
Modules = [EMRIWaveforms.ModeSelect, EMRIWaveforms.ModeSelect.KerrBackwards, EMRIWaveforms.ModeSelect.FEW]
Order = [:module, :type, :constant, :function]
```

### Mode summation

```@autodocs
Modules = [EMRIWaveforms.ModeSum, EMRIWaveforms.ModeSum.KerrBackwards, EMRIWaveforms.ModeSum.FEW]
Order = [:module, :type, :constant, :function]
```

### AAK

```@autodocs
Modules = [EMRIWaveforms.AAK, EMRIWaveforms.AAK.ParameterMapAAK]
Order = [:module, :type, :constant, :function]
```

### Mappings

```@autodocs
Modules = [EMRIWaveforms.Mappings]
Order = [:module, :type, :constant, :function]
```

### Waveform generation

```@autodocs
Modules = [EMRIWaveforms.Generate]
Order = [:module, :type, :constant, :function]
```

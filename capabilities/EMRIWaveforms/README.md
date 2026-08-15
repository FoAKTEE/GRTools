# EMRIWaveforms

GrTools tier-5 capability: the **EMRI trajectory + observed-waveform
pipeline**, fusing the two legacy waveform engines **FastEMRIWaveforms** and **KerrBackwards**.
This is **HDF5 home #2** (inspiral ODE → flux/amplitude grid models → mode selection → mode
summation → generate).

## Physics owned

End-to-end extreme-mass-ratio-inspiral (EMRI) waveform generation:

1. **Trajectory** — adiabatic inspiral ODE evolution of the orbital elements (vendored
   NumericsCore DOPR853/Verner integrator, Spacetime geodesics + separatrix stopping, Roots /
   NLsolve parameter inversion).
2. **FluxModels** — energy/angular-momentum flux RHS interpolated over the `Data` HDF5 flux grids
   (Schwarzschild-eccentric, Kerr-eccentric-equatorial).
3. **AmplitudeModels** — Teukolsky mode amplitudes from the RomanAmplitude neural network
   (`Kernels` amplitude GEMM) and the `Data` amplitude-grid bicubic interpolant.
4. **ModeSelect** — power-thresholded selection of the contributing `(l,m,n)` harmonics.
5. **ModeSum** — waveform assembly `h = Σ A·ₛYₗₘ·e^{-iΦ}` in the time domain, frequency-domain
   (SPA), and AAK approximations (summation kernels in `Kernels`).
6. **Mappings** — orbital-parameter coordinate maps (`p_to_y`, Kerr-ecc-eq, PN, Jacobians),
   reusing the `NumericsCore` conventions.
7. **Generate** — top-level observed-waveform driver (source frame → detector frame).

## Provenance

- the earlier `FastEMRIWaveforms` research code (FEW fork)
- the earlier `KerrBackwards` research code (KB fork)

Geometry collisions are migrated out to `Spacetime`; the spherical-harmonic `ₛYₗₘ` carrier to
`Harmonics`; shared numeric primitives (cubic/bicubic/tricubic interpolants, `ModeIndexTable`,
constants, Carlson, `leaky_relu`) are hoisted to `NumericsCore` / `Kernels`. The **4 dead FEW
import edges** (Teukolsky / TeukolskyPN / TeukolskySpinFluxes / SpinningOrbit — 0 qualified call
sites) are **deleted**, so this capability is **Symbolics-clean**.

## v0.1 is a TRANSITIONAL two-fork state

The collision survivors (`RomanAmplitude`, `SchwarzEccFlux`, `aak_waveform`, the mode-sum /
selector / inspiral driver / frame functions) exist in BOTH forks with differing conventions. In
v0.1 they are preserved as DISTINCT deep namespaces — `EMRIWaveforms.KerrBackwards.*` vs
`EMRIWaveforms.FEW.*` — never flattened, reachable only qualified, and collapsed to a single API
only after a per-submodule **bit-equivalence gate**. So `using EMRIWaveforms` is a documented
transitional surface, not yet one coherent convention.

## Precision / device

Float64 / GPU-only. It **calls** the BigFloat-generic lower tiers (Harmonics eigenvalues,
Perturbation, MetricReconstruction) at Float64 and imposes no precision constraint upward. GPU
work flows only through `Kernels` (Float64-on-device).

## Dependencies

- GrTools: `NumericsCore`, `Spacetime`, `Harmonics`, `Kernels`, `Data`.
- External: `HDF5` (DIRECT — the trigger that activates Data's HDF5 weakdep-extension for grid
  IO), `OrdinaryDiffEqVerner`, `Roots`, `NLsolve`, plus std `Printf`/`LinearAlgebra`/`Random`.

## Grounding / testing

No heavy live test-deps (FastEMRIWaveforms/KerrBackwards are NEVER live test-deps). Each function
is grounded tests-first by analytic identity (`p_to_y` round-trip, separatrix, hand-computed
mode-sum, parallel == serial bit-identical) and CAPTURED fixtures baked from the checksum-verified
`.fewdata` grids and one-shot legacy oracle outputs. The full Python end-to-end FEW match is a
documented `[HOLE]` (not buildable in-sandbox); every component is grounded individually.

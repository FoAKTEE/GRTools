# Kernels

The GPU/parallel compute-kernel spine of GrTools (capability **C4**, design §5.5 + §8). It owns the
hot inner loops that several physics capabilities share, each consolidated into a single canonical
copy — collapsing the FastEMRIWaveforms-vs-KerrBackwards forks of the *same* arithmetic. Everything
else (trajectory ODEs, HDF5 IO, Wigner couplings, Bessel machinery) stays with the calling
capability; Kernels receives it as precomputed data or an injected provider.

## Physics this capability owns

| module        | kernels | physics |
|---------------|---------|---------|
| `modesum.jl`  | [`direct_mode_sum`](@ref), [`interpolated_mode_sum`](@ref), [`fd_interpolated_mode_sum`](@ref), [`spa_func`](@ref), [`fd_flip_average`](@ref) | time-domain + frequency-domain (SPA) waveform mode-sum ``h = \sum_{lmn} A_{lmn}\, {}_sY_{lm}\, e^{-i\Phi_{mn}}`` |
| `aak.jl`      | [`aak_summation`](@ref) (+ internal `_aak_point`, `_aak_rot_coeff`) | Augmented-Analytic-Kludge harmonic summation |
| `r2.jl`       | [`r2_term`](@ref) dispatcher over the work-item kernels `r2_1!..r2_10!`; [`R2_COUPLING_SPEC`](@ref) | the 10 quadratic second-order-Ricci source kernels |
| `amp_gemm.jl` | [`roman_forward`](@ref), [`leaky_relu`](@ref) | ROMAN dense-NN Teukolsky-amplitude surrogate forward pass |
| `backend.jl`  | [`select_backend`](@ref), [`ka_device`](@ref), [`with_kernel_threads`](@ref) | KernelAbstractions backend selection over `NumericsCore.backend_iface`; the shared CPU execution seam |

All five files are `include`d into the single `Kernels` module — there are no nested submodules.
CUDA/AMDGPU support lives in the package extensions `ext/KernelsCUDAExt.jl` /
`ext/KernelsAMDGPUExt.jl`, loaded only when the operator brings the device package.

## Parallelism + precision contract (design §8, BINDING)

- **One work-item per output axis + a serial inner reduction** — encoded once in
  [`with_kernel_threads`](@ref) and shared by every CPU kernel — so the CPU-serial and CPU-threaded
  paths are **bit-identical** (the only bit-equality gate; it is TESTED).
- CPU partitioning is bounded ≤ 20 chunks; any LAPACK inside a work-item is pinned to a single BLAS
  thread via the mandatory `NumericsCore` BLAS guard.
- **Float64 ON DEVICE ONLY** — the BigFloat generation tier never touches the GPU; the type-generic
  CPU paths (e.g. the `_modesum_eltype` promotion in `modesum.jl`) preserve BigFloat inheritance.
- CPU-vs-GPU divergence (`sincos`/`fma`) is absorbed by a `:gpu` ULP tolerance class.

## Dependencies and poison-axis declaration

Post-S4 hard dep budget (design §5.5), enforced by the §7.7 `TestKit.assert_dep_budget` gate in
`core/Kernels/test/hygiene.jl`:

- **May carry (`allow`, closed budget):** `NumericsCore` (backend seam, BLAS guard, spline/mode-index
  primitives), `KernelAbstractions`, plus stdlib `LinearAlgebra` (the `mul!`/GEMM host).
- **Weakdeps only (`weakdeps_allow`):** `CUDA`, `AMDGPU` — realized as package extensions, never
  hard `[deps]`.
- **May NOT carry (`forbid`):** `HDF5`, `JLD2`, `Symbolics`, `OrdinaryDiffEq*`, `WignerSymbols`,
  `Data`, `Roots`, `NLsolve`, `SciMLBase`, `SpecialFunctions`, `HypergeometricFunctions`, `QuadGK`.
  A second closure-load gate (`TestKit.assert_not_loaded`) asserts none of these — nor CUDA/AMDGPU —
  is loaded after `using Kernels` in the CPU-only test env.

Two purity hoists keep the forbidden axes out (§5.5/§7.7.1):

- **Injected-provider AAK.** The Peters–Mathews Bessel builder `ak_bessel_coeffs` lives in
  `EMRIWaveforms.AAK`; [`aak_summation`](@ref) receives it as the `bessel_coeffs` keyword provider,
  so Kernels carries no `SpecialFunctions`. The provider cannot be a precomputed array: the
  coefficients depend on the spline-resampled eccentricity at the dense output nodes, which only
  exists inside the kernel — so the coefficient *builder* is injected (host-side dependency
  inversion) and the per-work-item body stays Bessel-free.
- **`R2_COUPLING_SPEC` contract.** The 3-harmonic angular coupling array `c` consumed by the ten
  `r2` kernels is precomputed by the caller (canonical host wrapper:
  `SelfForce.SecondOrder.r2_term`, which evaluates `Wigner.coupling`). [`R2_COUPLING_SPEC`](@ref)
  is the pure-data description of exactly which spin-weight triples ``(s_3, s_1, s_2)`` each kernel
  needs, so Kernels carries no Wigner/3j machinery.

## GPU / data needs

- **GPU is an optional, CPU-first extension.** The CPU backend is the dependency-free default;
  CUDA/AMDGPU enter only as `[weakdeps]` extensions that add the `gpu=true` method to
  `select_backend` and the device map to `ka_device`. GPU device runs are **[FUTURE]** this round.
- **No data artifacts.** Kernels reads no files: mode amplitudes, splines' input grids, NN weights
  (the HDF5 ROMAN weight loading stays in `EMRIWaveforms.AmplitudeModels`), first-order field
  dictionaries (`read_h1` layout, packed to dense matrices internally) and coupling arrays all
  arrive as in-memory arguments. Test grounding uses baked fixtures under `test/fixtures/`
  (`r2_oracle_ref.jl`, `r2_coupling_ref.jl`, `fdspa_{inputs,ref}.jl`, `aak_bessel_ref.jl`) — no
  live legacy dep.

## Legacy provenance (drawsFrom)

Rewrite-redesign, not copy; the legacy packages remain green as the characterization oracle and are
never live test deps (they pull OrdinaryDiffEq/HDF5/SWSH/Teukolsky closures).

| legacy package / files | landed here as |
|---|---|
| FastEMRIWaveforms `summation/{directmodesum,interpolatedmodesum,fdinterpolatedmodesum}.jl`; KerrBackwards `{summation,fdsummation}.jl` | `modesum.jl` — [`direct_mode_sum`](@ref), [`interpolated_mode_sum`](@ref), [`fd_interpolated_mode_sum`](@ref), [`spa_func`](@ref), [`fd_flip_average`](@ref) (single copy; `EMRIWaveforms.ModeSum` keeps the forked conventions and delegates) |
| FastEMRIWaveforms `summation/aaksummation.jl`; KerrBackwards `{summation,bessel_ak,parametermap_aak}.jl` | `aak.jl` — [`aak_summation`](@ref) (harmonic-sum kernel only; AK parameter maps + `ak_bessel_coeffs` stay in `EMRIWaveforms.AAK`) |
| SecondOrderRicci `r2terms.jl` + `r2_1.jl..r2_10.jl` (transliterations of the upstream `R2_*.cc`) | `r2.jl` — [`r2_term`](@ref) + `r2_1!..r2_10!` (math byte-identical; coupling evaluation hoisted to `SelfForce.SecondOrder`) |
| KerrBackwards `romanamplitude.jl` (`roman_amplitudes`, `leaky_relu`); FastEMRIWaveforms `amplitude/romannet.jl` (`leaky_relu`) | `amp_gemm.jl` — [`roman_forward`](@ref), [`leaky_relu`](@ref) (pure GEMM; weight IO and conj-partner mode selection stay in `EMRIWaveforms.AmplitudeModels`) |

## What changed vs legacy

- **Fork collapse.** The FEW and KerrBackwards copies of the same mode-sum/AAK/ROMAN arithmetic are
  collapsed to one canonical kernel each; the capabilities keep thin, convention-preserving facades.
- **FD-SPA is NEW here (S4-a).** The frequency-domain stationary-phase kernel
  ([`fd_interpolated_mode_sum`](@ref)/[`spa_func`](@ref)/[`fd_flip_average`](@ref)) was hoisted
  from the transitional `EMRIWaveforms.ModeSum.FEW` local copy into Kernels, fulfilling the §5.5
  "TD + FD-SPA + AAK mode-sum live in Kernels" single-copy mandate; `EMRIWaveforms` now delegates.
- **Dependency inversion instead of dependency import.** Legacy kernels evaluated Bessel
  coefficients and Wigner couplings inline; here the AAK Bessel *builder* is injected
  (`bessel_coeffs` provider) and the `r2` couplings arrive precomputed per
  [`R2_COUPLING_SPEC`](@ref), keeping the hard dep budget at {`NumericsCore`,
  `KernelAbstractions`} exactly.
- **Deterministic threading.** Legacy ad-hoc loops are replaced by the shared
  [`with_kernel_threads`](@ref) seam (bounded ≤ 20 chunks, disjoint outputs), making
  `parallel == serial` a tested bit-identity rather than a hope.
- **No Dict hashing in hot loops.** The `(i,l,m)`-keyed first-order field Dicts of the legacy `r2`
  path are packed once into dense `10×N` complex matrices before kernel launch.

## Public API

```@autodocs
Modules = [Kernels]
Order = [:module, :type, :constant, :function]
```

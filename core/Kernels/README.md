# Kernels

The GPU/parallel compute-kernel spine of GrTools (capability **C4**, design §5.5 + §8). It owns the hot
inner loops that several physics capabilities share, each consolidated into a single canonical copy
(collapsing the FastEMRIWaveforms-vs-KerrBackwards forks of the *same* arithmetic — AHA).

## Physics this capability owns

| module        | kernel | physics |
|---------------|--------|---------|
| `modesum.jl`  | `direct_mode_sum` / `interpolated_mode_sum` / `fd_interpolated_mode_sum` / `spa_func` … | time-domain + frequency-domain (SPA) waveform mode-sum `h = Σ_{lmn} A·{}_sY_{lm}·e^{-iΦ}` |
| `aak.jl`      | `aak_summation` / `_aak_point` / `_aak_rot_coeff` … | Augmented-Analytic-Kludge harmonic summation |
| `r2.jl`       | `r2_1`…`r2_10` + `r2_term` | the 10 quadratic second-order-Ricci kernels |
| `amp_gemm.jl` | `roman_forward` / `leaky_relu` | ROMAN dense-NN amplitude surrogate forward pass |
| `backend.jl`  | `select_backend` / `ka_device` | KernelAbstractions backend selection over `NumericsCore.backend_iface` |

## Legacy packages absorbed (rewrite-redesign, not copy)

- **FastEMRIWaveforms** — `summation/{directmodesum,interpolatedmodesum,fdinterpolatedmodesum,aaksummation}`, `amplitude/romannet` (`leaky_relu`).
- **KerrBackwards** — `summation` / `fdsummation` (mode-sum, FD-SPA, AAK), `romanamplitude` (ROMAN GEMM).
- **SecondOrderRicci** — `r2terms` + `r2_1..r2_10` (the 10 quadratic kernels).

These remain green as the characterization oracle; Kernels never edits them and never live-deps them in
test (they pull OrdinaryDiffEq/HDF5/SWSH/Teukolsky closures). Grounding is done WITHOUT a live dep:
hand-computed reference sums + `parallel==serial` bit-identity for modesum/aak/amp_gemm, and captured
seeded-field g++/GSL fixtures for `r2_1..r2_10`.

## Parallelism + precision contract (design §8, BINDING)

- **One work-item per output axis + serial inner reduction** ⇒ CPU-serial == CPU-threaded **bit-identical**
  (the only bit-equality gate — it is TESTED).
- CPU partitioning bounded ≤20; mandatory `NumericsCore` BLAS guard (`with_blas_threads`) around any LAPACK.
- **Float64 ON DEVICE ONLY** — the BigFloat generation tier never touches the GPU; the type-generic CPU
  paths preserve BigFloat inheritance.
- CPU-vs-GPU divergence (`sincos`/`fma`) is absorbed by a `:gpu` ULP tolerance class.

## GPU is an optional, CPU-first extension

The CPU backend is the dependency-free default. CUDA/AMDGPU are `[weakdeps]` realized by
`ext/KernelsCUDAExt.jl` and `ext/KernelsAMDGPUExt.jl`, loaded only when the operator brings the device
package in. GPU device runs are **[FUTURE]** this round.

## Public API (scaffold)

`select_backend`, `CPUBackend`, `AbstractBackend`, `backend_array`, `backend_collect` (re-exported from
NumericsCore), `ka_device`. Kernel entry points (`direct_mode_sum`, `aak_summation`, `r2_term`,
`roman_forward`, …) land tests-first in subsequent rounds.

## Dependencies

`NumericsCore` (backend seam, BLAS guard), `KernelAbstractions`. Weak: `CUDA`, `AMDGPU`.

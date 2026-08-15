"""
    Kernels

GPU/parallel compute-kernel spine of GrTools (capability C4, design §5.5 + §8).

This package owns the *hot inner loops* that several physics capabilities share, consolidated into a
single canonical copy each (AHA — collapse the FEW-vs-KerrBackwards forks of the same arithmetic):

  - `modesum.jl`  — time-domain (`direct_mode_sum`/`interpolated_mode_sum`) + frequency-domain SPA
                    (`fd_interpolated_mode_sum`/`spa_func`/…) waveform mode-sum `h = Σ A·sYlm·e^{-iΦ}`.
  - `aak.jl`      — Augmented-Analytic-Kludge harmonic summation (`aak_summation`/`_aak_point`/…).
  - `r2.jl`       — the 10 quadratic second-order-Ricci kernels `r2_1..r2_10` + `r2_term` dispatcher.
  - `amp_gemm.jl` — ROMAN dense-NN forward pass (`roman_forward`/`leaky_relu`).
  - `backend.jl`  — KernelAbstractions backend selection layered over `NumericsCore.backend_iface`.

Absorbs (rewrite-redesign, not copy) the kernel cores of legacy `FastEMRIWaveforms`, `KerrBackwards`
(mode-sum / FD-SPA / AAK / ROMAN GEMM) and `SecondOrderRicci` (`r2_1..r2_10`).

PARALLELISM CONTRACT (design §8, BINDING): every kernel is **one work-item per output axis + a serial
inner reduction**, so the CPU-serial and CPU-threaded paths are **bit-identical** (the bit-equality
gate). CPU partitioning is bounded ≤20 with a mandatory `NumericsCore.blas_guard` around any LAPACK.
Float64 ON DEVICE ONLY — the BigFloat generation tier never touches the GPU.

GPU is a CPU-FIRST optional extension: CUDA/AMDGPU are `[weakdeps]` (`ext/KernelsCUDAExt.jl`,
`ext/KernelsAMDGPUExt.jl`); the CPU backend is the dependency-free default.

DEP BUDGET (design §5.5, BINDING — the §7.7.1 dep-budget gate): `[deps]` == {`NumericsCore`,
`KernelAbstractions`} EXACTLY (+ stdlib `LinearAlgebra`); CUDA/AMDGPU only as `[weakdeps]`. In
particular NO `Wigner` (the r2 coupling arrays are precomputed by the caller per `R2_COUPLING_SPEC`;
host wrapper: `SelfForce.SecondOrder.r2_term`) and NO `SpecialFunctions` (the AAK Peters–Mathews
Bessel builder `ak_bessel_coeffs` lives in `EMRIWaveforms.AAK` and is injected as `aak_summation`'s
`bessel_coeffs` provider).
"""
module Kernels

using LinearAlgebra
using NumericsCore
using KernelAbstractions

# --- Backend selection (KA over NumericsCore.backend_iface; CPU default) ----------------------------
include("backend.jl")

# --- Compute kernels (modules built tests-first in subsequent rounds) -------------------------------
include("modesum.jl")
include("aak.jl")
include("r2.jl")
include("amp_gemm.jl")

end # module Kernels

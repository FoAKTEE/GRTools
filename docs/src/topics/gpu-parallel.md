# GPU / Parallel

Single-source topic page (design §8 / §12): **this page is the home of the GrTools parallelism and
GPU-gating rules**. Package pages ([Kernels](../core/kernels.md), [NumericsCore](../core/numericscore.md),
[EMRIWaveforms](../capabilities/emriwaveforms.md), [SelfForce](../capabilities/selfforce.md),
[Hybrid](../capabilities/hybrid.md), [AccretionDisc](../capabilities/accretiondisc.md)) link here rather
than restating the contract.

The model in one sentence: **heavy kernels** (waveform mode-sum, FD-SPA, AAK summation, ROMAN GEMM, the
ten second-order-Ricci `r2` kernels) live in the `Kernels` package and may run on an *operator-leased*
GPU; **everything else** is CPU threads/processes over mode/orbit/trajectory axes, bounded to ≤ 20
threads, with a mandatory BLAS guard and a parallel == serial bit-identity gate.

## Heavy-run GPU gating (operator-leased, never free-memory probed)

GPU support is a **CPU-first optional extension** of `Kernels`:

- `CUDA` and `AMDGPU` are `[weakdeps]` of `Kernels` (see `core/Kernels/Project.toml`); the
  package extensions `ext/KernelsCUDAExt.jl` and `ext/KernelsAMDGPUExt.jl` load only when the operator's
  environment brings the GPU package in. The hard `[deps]` are exactly `NumericsCore` +
  `KernelAbstractions` (+ stdlib `LinearAlgebra`) — the CPU path is GPU-dependency-free.
- The dispatch seam lives one tier lower, in `NumericsCore` (`src/backend_iface.jl`):
  `AbstractBackend` / `CPUBackend` / `select_backend` / `backend_array` / `backend_collect`.
  The NumericsCore fallback `select_backend(; gpu=false, kind=:cpu)` **always returns `CPUBackend()`** —
  even when `gpu=true` is requested but no extension is loaded.
- Each extension injects a device tag (`CUDAGPUBackend`, `ROCMGPUBackend`) and the more-specific
  `select_backend` method that honors `gpu=true, kind=:cuda` / `kind=:rocm`. `Kernels.ka_device`
  (`src/backend.jl`) maps tags to concrete KernelAbstractions devices: `CPUBackend → KA.CPU()`,
  and `CUDA.CUDABackend()` / `AMDGPU.ROCBackend()` from the extensions. `backend_array` /
  `backend_collect` allocate/collect on the right side of the boundary.

**Lease model.** GPU occupancy is *operator-signalled*: the device is used only when the operator
(a) loads `CUDA`/`AMDGPU` into the environment and (b) explicitly requests
`select_backend(gpu=true, kind=:cuda)` (or `:rocm`). There is deliberately **no free-memory probing** —
no code path in GrTools queries device memory or `CUDA.functional()`-style availability to auto-select a
backend (grep-verified). When the massive GPU is leased to another heavy run, the operator simply does
not request it, and everything resolves to `CPUBackend()`.

Status: the extensions are declared CPU-first (`[FUTURE]` in the extension headers) — actual on-device
runs are gated to a heavy GPU round. **Float64 on device only** (see the precision boundary below).

## The CPU ≤ 20 fallback

When the GPU is unavailable or leased, heavy kernels fall back to bounded CPU threading through one
shared execution seam, `Kernels.with_kernel_threads(f, n; threaded=false, maxthreads=20)`
(`core/Kernels/src/backend.jl`):

- serial by default (`threaded=false`), so light/test runs pay no task overhead;
- the threaded path partitions `1:n` into at most `maxthreads` (default and policy cap: **20**)
  contiguous chunks dispatched with `Threads.@spawn` — it never spawns one task per index;
- the whole threaded region runs inside `blas_single_thread` (next section).

Every CPU kernel in `Kernels` (`direct_mode_sum`, `interpolated_mode_sum`, `fd_interpolated_mode_sum`,
`aak_summation`, `roman_forward`, `r2_term`) dispatches its output axis through this seam. Capability-level
batch loops follow the same launch policy — e.g. both `int_inspiral_batch` methods in
`Hybrid/src/Drivers/inspiral_driver.jl` carry the "keep Julia nthreads ≤ 20 CPUs" rule at their
`Threads.@threads` sites. Rationale (measured, per project memory): heavy BigFloat/LA workloads scale
*worse* beyond a few threads due to GC + BLAS contention — scale across **processes over modes**, not by
oversubscribing threads.

Note on the DG time-domain solver: `Perturbation.TimeDomain` (and the `NumericsCore.dg` core plus the
`:GSN` radial backend) contains **no threaded region** — each `solve_*` evolution is serial per
``(\ell, m)`` problem (grep-verified: zero `Threads.@threads`/`@spawn`/`Distributed` sites in
`Perturbation/src/TimeDomain/`, `Radial/gsn_*.jl`, `NumericsCore/src/dg.jl`). Parallelism over a TD mode
sweep is caller-side (independent problems, processes or bounded threads under the same ≤ 20 policy);
no bit-identity gate is needed because no in-package threaded path exists.

## The BLAS-guard rule (binding)

Running OpenBLAS from inside a `Threads.@threads` / `@spawn` region oversubscribes and can **deadlock**
on the BLAS internal locks (nested Julia↔OpenBLAS threading), and a multi-thread BLAS reduction is
nondeterministic in the last bits. The rule (§8, binding): *every threaded LAPACK/BLAS region pins
OpenBLAS to one thread, scoped*.

The canonical primitives live in `NumericsCore` (`src/blas_guard.jl`; tests in
`test/test_blas_guard.jl`):

```julia
with_blas_threads(f, n)   # run f() with OpenBLAS pinned to n threads, restored on exit (try/finally)
blas_single_thread(f)     # = with_blas_threads(f, 1) — the deadlock-safe, bit-deterministic setting
```

Real call sites (all `NumericsCore.blas_single_thread`):

| Site | What it guards |
|---|---|
| `Kernels/src/backend.jl` — `with_kernel_threads` | every bounded-threaded kernel sweep (any LAPACK inside a work-item) |
| `SelfForce/src/LorenzField/Solver.jl` — `calculate_scaling_coefficients` | the dense `Φ \ J` scaling-coefficient solve |
| `SelfForce/src/LorenzField/Solver.jl` — `compute_modes` | the threaded `(l,m)` mode sweep (`parallel=:threads`) |
| `Hybrid/src/Drivers/inspiral_driver.jl` — `int_inspiral_batch` (2-var and 7-var-spin methods) | `Threads.@threads` batch of independent inspiral solves |
| `EMRIWaveforms/src/AmplitudeModels.jl` — `AmpInterpKerrEqEcc` constructor | threaded per-mode B-spline coefficient-bank build |
| `AccretionDisc/src/Force/Pipeline.jl` — `build_forcing` (+ the per-slot source assembly) | the 𝒬 4×4 LAPACK inverse inside a threaded sweep |
| `AccretionDisc/src/Force/KerrFrontend.jl` — `kerr_mode` / `kerr_modes` | race-free preallocated per-slot Kerr mode builds |
| `MetricReconstruction/src/Lorenz/Closure.jl` — unit-probe matrix assembly (`probe_parallel=true`) | LAPACK inside the threaded closure probe |
| `Harmonics/src/SpheroidalEigenvalues.jl` — `_swsh_eigenvalue_spectral` / `_swsh_spectral_coeffs` / `_swsh_eigenvalue_dgamma` | `eigvals`/`eigen(Symmetric(M))` so callers may thread over modes safely |

## parallel == serial bit-identity gates

The §8 parallelism contract: every kernel is **one work-item per output index + a serial inner
reduction inside the work-item**, and each work-item writes *only its own output index* (disjoint
writes). Consequence: there is no cross-item reduction whose summation order could change, so the
CPU-serial and CPU-threaded paths are **bit-identical** — and CI asserts it with `===`/`isequal`
element-wise, not with a tolerance. This is the only bit-equality gate; **CPU-vs-GPU is *not*
bit-gated** — device `sincos`/`fma` divergence is absorbed by a `:gpu` ULP tolerance class (design §8).

Where the gates live:

| Test file | Gate |
|---|---|
| `Kernels/test/test_backend.jl` | `with_kernel_threads`: exactly-once index coverage + serial `===` threaded |
| `Kernels/test/test_modesum.jl` | `direct_mode_sum`, `interpolated_mode_sum`, `fd_interpolated_mode_sum` (NaN-aware `isequal`, combined + separate-modes paths) |
| `Kernels/test/test_aak.jl` | `aak_summation` over dense nodes |
| `Kernels/test/test_amp_gemm.jl` | `roman_forward` over output rows |
| `Kernels/test/test_r2.jl` | `r2_term` (all ten `r2_1..r2_10` couplings) |
| `AccretionDisc/test/Force/test_pipeline.jl` | threaded `build_forcing` == serial (+ speedup check) |
| `Hybrid/test/test_drivers.jl` | 2-var inspiral batch serial == threaded |
| `MetricReconstruction/test/lorenz_closure.jl` | `probe_parallel=true` == serial closure assembly |
| `MetricReconstruction/test/radiation_assembly.jl` | threaded `h_coord_modesum` == serial reference sum |

If you add a threaded path anywhere, you inherit this contract: disjoint per-index writes, BLAS guard,
and a serial-vs-threaded bit-identity test *before* claiming any speedup.

## BigFloat-CPU / Float64-GPU boundary

The high-precision **generation tier** — Harmonics spheroidal eigenvalues, Perturbation MST/CF,
the MetricReconstruction Kerr radial functions — runs at 256–2048-bit `BigFloat` and **never touches
the device**. `Kernels` sits on the other side of the boundary: it is **Float64-only on device** and
imposes no precision constraint upward (the BigFloat-generic lower tiers are simply *called* at
Float64 by the waveform tier). The `BigFloat → Float64 → GPU-grid` handoff is gated end-to-end
against the existing waveform oracle tolerances, not bit-identity. For the BigFloat inheritance chain
itself (which types are `{T}`-generic, and the named Float64 ceilings like
`angular_precision_floor()`), see [Precision](precision.md).

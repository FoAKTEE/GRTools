# KernelAbstractions backend selection layered over NumericsCore.backend_iface.
#
# NumericsCore owns the GPU-dependency-free `AbstractBackend`/`CPUBackend`/`select_backend` seam. Here we
# map those tags onto a concrete KernelAbstractions device so kernels launch uniformly across CPU/GPU. The
# CPU device is the default; CUDA/AMDGPU devices are injected by the package extensions, which add the
# `gpu=true` method to `select_backend` (NumericsCore's fallback returns `CPUBackend()`).

import KernelAbstractions as KA

# Re-export the NumericsCore backend seam so downstream capabilities depend on Kernels alone for dispatch.
export AbstractBackend, CPUBackend, select_backend, backend_array, backend_collect
export ka_device, with_kernel_threads

"""
    ka_device(backend::AbstractBackend) -> KernelAbstractions backend

Concrete KernelAbstractions device for a GrTools `AbstractBackend`. The host `CPUBackend` maps to
`KernelAbstractions.CPU()`; GPU backends (from the CUDA/AMDGPU extensions) override with `CUDABackend()`
/ `ROCBackend()`. One place for the tag→device map keeps the kernel launch sites backend-agnostic.
"""
ka_device(::CPUBackend) = KA.CPU()

"""
    with_kernel_threads(f, n; threaded=false, maxthreads=20)

Execution seam every CPU kernel here shares: apply the per-output-index closure `f(i)` for `i = 1 … n`.

This encodes the design's parallelism contract (§8): **one work-item per output index + a serial inner
reduction inside `f`**, so the serial path (`threaded=false`) and the bounded-threaded path are
**bit-identical** — there is no cross-item reduction whose summation order could change. The threaded
path partitions `1:n` into at most `maxthreads` (≤20, bounded) contiguous chunks dispatched with
`Threads.@spawn`; it never spawns one task per index. Any LAPACK inside `f` is pinned to a single BLAS
thread (`NumericsCore.blas_single_thread`) so the threaded region neither oversubscribes nor loses
determinism. `f` must write only its own index `i` (disjoint outputs) for the bit-equality gate to hold.
"""
function with_kernel_threads(f, n::Integer; threaded::Bool = false, maxthreads::Integer = 20)
    if !threaded || n <= 1
        @inbounds for i in 1:n
            f(i)
        end
        return nothing
    end
    nchunks = min(Int(maxthreads), Int(n))
    base, rem = divrem(Int(n), nchunks)
    blas_single_thread() do
        @sync begin
            lo = 1
            for c in 1:nchunks
                len = base + (c <= rem ? 1 : 0)
                hi = lo + len - 1
                rng = lo:hi
                Threads.@spawn @inbounds for i in rng
                    f(i)
                end
                lo = hi + 1
            end
        end
    end
    return nothing
end

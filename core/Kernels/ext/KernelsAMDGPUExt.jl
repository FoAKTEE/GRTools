# AMDGPU (ROCm) device extension for Kernels (design §5.5 / §8). Loads only when the operator brings in
# AMDGPU; keeps the CPU path dependency-free. [FUTURE] — declared CPU-first; actual on-device runs are
# gated to a heavy GPU round. Float64 ON DEVICE ONLY.
module KernelsAMDGPUExt

import AMDGPU
import KernelAbstractions as KA
import NumericsCore
using NumericsCore: AbstractBackend, CPUBackend
import Kernels

"GPU backend tag for AMD ROCm devices (the `gpu=true` resolution of `select_backend`)."
struct ROCMGPUBackend <: AbstractBackend end

NumericsCore.select_backend(; gpu::Bool=false, kind::Symbol=:cpu) =
    (gpu && kind === :rocm) ? ROCMGPUBackend() : CPUBackend()

Kernels.ka_device(::ROCMGPUBackend) = AMDGPU.ROCBackend()
NumericsCore.backend_array(::ROCMGPUBackend, ::Type{T}, dims::Integer...) where {T} = AMDGPU.ROCArray{T}(undef, dims...)
NumericsCore.backend_collect(::ROCMGPUBackend, src::AbstractArray) = Array(src)

end # module KernelsAMDGPUExt

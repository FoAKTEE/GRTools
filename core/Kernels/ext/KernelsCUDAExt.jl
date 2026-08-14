# CUDA device extension for Kernels (design §5.5 / §8). Loads only when the operator brings in CUDA;
# keeps the CPU path dependency-free. [FUTURE] — declared CPU-first; actual on-device runs are gated to a
# heavy GPU round. Float64 ON DEVICE ONLY.
module KernelsCUDAExt

import CUDA
import KernelAbstractions as KA
import NumericsCore
using NumericsCore: AbstractBackend, CPUBackend
import Kernels

"GPU backend tag for CUDA devices (the `gpu=true` resolution of `select_backend`)."
struct CUDAGPUBackend <: AbstractBackend end

# More-specific method honoring gpu=true once CUDA is loaded; NumericsCore's fallback handles the rest.
NumericsCore.select_backend(; gpu::Bool=false, kind::Symbol=:cpu) =
    (gpu && kind === :cuda) ? CUDAGPUBackend() : CPUBackend()

Kernels.ka_device(::CUDAGPUBackend) = CUDA.CUDABackend()
NumericsCore.backend_array(::CUDAGPUBackend, ::Type{T}, dims::Integer...) where {T} = CUDA.CuArray{T}(undef, dims...)
NumericsCore.backend_collect(::CUDAGPUBackend, src::AbstractArray) = Array(src)

end # module KernelsCUDAExt

# CPU/GPU dispatch seam — the `AbstractBackend` abstraction that `GrToolsKernels` extends via a package
# extension. This is the Kernels boundary: it keeps NumericsCore GPU-dependency-free (NO CUDA/AMDGPU/
# KernelAbstractions referenced here). A GPU backend type + the `gpu=true` branch of `select_backend`
# are injected by the Kernels extension; without it everything resolves to the host `CPUBackend`.

"Supertype of compute backends. Concrete subtypes (`CPUBackend` here; `GPUBackend` from the Kernels ext)
select where arrays live and how kernels dispatch."
abstract type AbstractBackend end

"Host backend: plain `Base.Array` storage, serial/Threads execution. The dependency-free default."
struct CPUBackend <: AbstractBackend end

"""
    select_backend(; gpu=false, kind=:cpu) -> AbstractBackend

Resolve a compute backend. Returns `CPUBackend()` by default and whenever no GPU extension is loaded —
including when `gpu=true` is requested but the `GrToolsKernels` GPU extension (which adds the more
specific method that honors `gpu`/`kind`) is absent. Keeping the fallback here means NumericsCore never
hard-depends on a GPU package.
"""
select_backend(; gpu::Bool=false, kind::Symbol=:cpu) = CPUBackend()

"Allocate an uninitialized `T` array of the given dims on `backend`. Host path = `Array{T}(undef, …)`."
backend_array(::CPUBackend, ::Type{T}, dims::Integer...) where {T} = Array{T}(undef, dims...)

"Copy `src` into host memory for `backend`, returning a plain `Array` (identity copy on `CPUBackend`)."
backend_collect(::CPUBackend, src::AbstractArray) = collect(src)

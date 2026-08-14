using Test
using Kernels

@testset "Kernels" begin
    @testset "smoke: loads + backend seam" begin
        # `using Kernels` resolves; the CPU backend seam is re-exported from NumericsCore and maps onto a
        # concrete KernelAbstractions device (the scaffold's only live surface — kernels land tests-first
        # in subsequent rounds).
        b = select_backend()
        @test b isa CPUBackend
        # KernelAbstractions is a transitive dep (not a test dep), so identify the device by name.
        @test nameof(typeof(ka_device(b))) === :CPU
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests
    include("test_backend.jl")
    include("test_modesum.jl")
    include("test_amp_gemm.jl")
    include("test_aak.jl")
    include("test_r2.jl")
end

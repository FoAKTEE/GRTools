# Behavior-contract tests for the CPU/GPU dispatch seam (no legacy oracle — genuinely new).
# AbstractBackend / CPUBackend / select_backend is the boundary GrToolsKernels extends via a package
# extension; NumericsCore itself stays GPU-dependency-free.

using NumericsCore
using NumericsCore: AbstractBackend, CPUBackend, select_backend, backend_array, backend_collect

@testset "backend_iface" begin

    @testset "type relations" begin
        @test CPUBackend <: AbstractBackend
        @test CPUBackend() isa AbstractBackend
        @test isa(CPUBackend(), CPUBackend)
    end

    @testset "select_backend defaults to CPU" begin
        @test select_backend() isa CPUBackend
        @test select_backend(; gpu=false) isa CPUBackend
        # gpu=true WITHOUT the Kernels GPU extension loaded must still resolve to CPUBackend
        # (no hard GPU dependency lives here); a future GPU ext overrides this dispatch.
        @test select_backend(; gpu=true) isa CPUBackend
        @test select_backend(; kind=:cpu) isa CPUBackend
        @test select_backend(; gpu=true, kind=:cuda) isa CPUBackend
    end

    @testset "array allocation / host-copy roundtrip identity on CPUBackend" begin
        b = CPUBackend()
        a = backend_array(b, Float64, 3, 4)
        @test a isa Array{Float64,2}
        @test size(a) == (3, 4)
        src = collect(reshape(1.0:12.0, 3, 4))
        dst = backend_collect(b, src)
        @test dst == src
        @test dst isa Array{Float64,2}
        # BigFloat element type roundtrip (precision inheritance)
        bigsrc = BigFloat[big"1.5", big"2.25", big"3.125"]
        @test backend_collect(b, bigsrc) == bigsrc
        @test backend_array(b, BigFloat, 2) isa Vector{BigFloat}
    end

    @testset "no GPU package imported in src (dependency contract)" begin
        # the contract is that NumericsCore has NO GPU dependency: assert no `using`/`import` line
        # pulls a GPU package (prose mentions in comments are fine).
        codelines = filter(readlines(joinpath(@__DIR__, "..", "src", "backend_iface.jl"))) do ln
            s = strip(ln)
            startswith(s, "using ") || startswith(s, "import ")
        end
        @test isempty(codelines)
        for gpu in ("CUDA", "AMDGPU", "Metal", "oneAPI", "KernelAbstractions", "GPUArrays")
            @test !any(ln -> occursin(gpu, ln), codelines)
        end
        # and the package's declared deps carry no GPU package
        proj = read(joinpath(@__DIR__, "..", "Project.toml"), String)
        for gpu in ("CUDA", "AMDGPU", "Metal", "oneAPI", "GPUArrays")
            @test !occursin(gpu, proj)
        end
    end
end

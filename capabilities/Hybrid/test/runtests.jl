using Test
using Hybrid

@testset "Hybrid" begin
    @testset "smoke" begin
        @test isdefined(Main, :Hybrid)
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module test files (added TDD per round):
    include("test_models_interpolation.jl")
    include("test_binding_energy.jl")
    include("test_inspiral_0PASchwarz.jl")
    include("test_amplitude_0PASchwarz.jl")
    include("test_inspiral_1PAT1.jl")
    include("test_inspiral_1PAT1ea.jl")
    include("test_inspiral_0PAKerr.jl")
    include("test_amplitude_1PAT1.jl")
    include("test_amplitude_1PAT1ea.jl")
    include("test_amplitude_0PAKerr.jl")
    include("test_drivers.jl")
    include("test_registry.jl")
    include("test_binary_waveform.jl")
    include("test_introspection.jl")
    include("test_hybridc_deps.jl")
    include("test_hybridc_codegen.jl")
    include("test_hybridc_inspiral.jl")
    include("test_data_jld2_ext.jl")
end

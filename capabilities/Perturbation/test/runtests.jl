using Test
using Perturbation

@testset "Perturbation" begin
    @testset "smoke" begin
        @test isdefined(Main, :Perturbation)
        @test Perturbation isa Module
    end

    @testset "firewall closure-load gate" begin
        # Design §4.3: Perturbation is the LIGHT-DEP FIREWALL — its runtime closure must carry
        # NONE of the poison axes. This asserts the package was NOT loaded with any of
        # HDF5/JLD2/Symbolics/OrdinaryDiffEq/CUDA/AMDGPU/KernelAbstractions in its dependency
        # closure (they would appear as loaded modules only if a transitive runtime dep pulled
        # them). External-ODE-free / HDF5-free is a structural invariant of this capability.
        poison = ("HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq", "OrdinaryDiffEqVerner",
                  "CUDA", "AMDGPU", "KernelAbstractions")
        loaded = Set(string(m) for m in values(Base.loaded_modules))
        for p in poison
            @test !(p in loaded)
        end
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (TOML-static dep budget + closure-load),
                            # on top of the (untouched) firewall testset above; before heavy tests

    # Per-module grounded test files (test_radial.jl, test_amplitudes_fluxes.jl,
    # test_scattering.jl, ...) are added TDD in the BUILD rounds — grounded vs the live LIGHT
    # oracles Teukolsky + ScatteringPhaseMST at their frozen tolerances + analytic identities.
    # Each is included here as it lands.
    for f in sort(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end

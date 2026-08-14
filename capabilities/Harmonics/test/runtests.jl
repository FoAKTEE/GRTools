using Test
using Harmonics

@testset "Harmonics" begin
    @testset "smoke" begin
        @test isdefined(Main, :Harmonics)
        @test Harmonics isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module grounded test files (test_spheroidal_eigenvalues.jl, test_spheroidal_harmonics.jl,
    # test_spherical_harmonics.jl, test_precision.jl) are added TDD in the BUILD rounds —
    # grounded vs the SpinWeightedSpheroidalHarmonics oracle at its tolerances + analytic
    # identities. Each is included here as it lands.
    for f in sort(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end

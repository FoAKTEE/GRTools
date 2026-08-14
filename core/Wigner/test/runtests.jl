using Test
using Wigner

@testset "Wigner" begin
    @testset "smoke" begin
        @test isdefined(Main, :Wigner)
        @test Wigner isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module grounded test files (test_wigner3j.jl, test_clebsch_gordan.jl, test_coupling.jl,
    # test_wignerD.jl) are added TDD in the BUILD rounds; include each as it lands.
    for f in sort(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end

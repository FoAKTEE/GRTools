using Test
using Wigner

@testset "Wigner" begin
    @testset "smoke" begin
        @test isdefined(Main, :Wigner)
        @test Wigner isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # EXPLICIT test-file list (was `readdir` discovery: an untracked file in a dirty working
    # tree silently entered the suite, and a renamed/deleted file silently left it — repro-audit
    # G19, same fix as EllipticCore/Perturbation). Keep the list ALPHABETICAL: that is the order
    # `sort(readdir(...))` used, so the sequence in which the files run is unchanged. Adding a
    # test file means adding it here.
    TEST_FILES_WIGNER = [
        "test_clebsch_gordan.jl",
        "test_coupling.jl",
        "test_wigner3j.jl",
        "test_wignerD.jl",
    ]

    @testset "test-file registry is complete and exact" begin
        listed = Set(TEST_FILES_WIGNER)
        ondisk = Set(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        @test isempty(setdiff(listed, ondisk))    # every listed file exists ...
        @test isempty(setdiff(ondisk, listed))    # ... and no test_*.jl here is unregistered
    end

    for f in TEST_FILES_WIGNER
        include(joinpath(@__DIR__, f))
    end
end

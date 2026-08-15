using Test
using EllipticCore
using REPL   # Julia >=1.11: docstring lookup (@doc) is a no-op unless REPL is loaded

@testset "EllipticCore" begin
    # Binding hygiene gates (dep budget + closure-load poison) run FIRST, before any other
    # test file (elliptic-ssf TOPOLOGICAL_SCHEDULE.md W1-A completion criteria).
    include("hygiene.jl")

    @testset "smoke: scaffold" begin
        @test EllipticCore isa Module
        @test pkgversion(EllipticCore) == v"0.1.0"
        docstr = string(@doc EllipticCore)
        @test occursin("EllipticCore (GrToolsElliptic)", docstr)   # alias-line convention
        @test occursin("LinearAlgebra, SparseArrays, NumericsCore", docstr)  # dep budget stated
    end

    # Per-module grounded test files (test_refelement.jl, test_mesh2d.jl, …) land TDD with
    # their own DAG nodes; each is picked up automatically as it arrives.
    for f in sort(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end

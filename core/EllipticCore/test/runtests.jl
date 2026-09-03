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

    # EXPLICIT test-file list (was `readdir` discovery: an untracked file in a dirty working
    # tree silently entered the suite, and a renamed/deleted file silently left it). Keep the
    # list ALPHABETICAL — that is the order `sort(readdir(...))` used, so the sequence in which
    # the files run is unchanged. Adding a test file means adding it here.
    TEST_FILES_EC = [
        "test_amr.jl", "test_bottom_solve.jl",
        "test_boundary_hooks.jl", "test_dg_matrices.jl",
        "test_dg_operator.jl", "test_dg_penalty.jl",
        "test_gmres.jl", "test_matrixfree.jl",
        "test_mesh2d.jl", "test_mortar_hooks.jl", "test_mortars.jl", "test_multigrid.jl",
        "test_refelement.jl", "test_schwarz.jl",
        "test_schwarz_shift.jl", "test_solver_stack.jl",
    ]

    @testset "test-file registry is complete and exact" begin
        listed = Set(TEST_FILES_EC)
        ondisk = Set(filter(n -> startswith(n, "test_") && endswith(n, ".jl"), readdir(@__DIR__)))
        # every listed file exists ...
        @test isempty(setdiff(listed, ondisk))
        # ... and no test_*.jl in this directory is unregistered (a new test file must be added
        # to the list above; a stray scratch file must not be left in test/).
        @test isempty(setdiff(ondisk, listed))
    end

    for f in TEST_FILES_EC
        include(joinpath(@__DIR__, f))
    end
end

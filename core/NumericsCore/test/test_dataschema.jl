# Characterization of the dep-free DataSchema contract structs.
#
# Pure structs (no compute, no HDF5) ⇒ grounding is ROUND-TRIP + FIELD-INVARIANT. The H1 layout
# literals (inhom_left cols = r0index+1, inhom_right cols = gs−r0index, rows = 4·nf) are read from
# packages/H1Lorenz/src/io_hdf5.jl:29-56 (NOT live-dep'd — H1Lorenz pulls HDF5). The amplitude/flux
# axis invariants (≥3 strictly-increasing nodes) come from FEW ampinterp2d / kerrflux.

using NumericsCore
using NumericsCore: H1Schema, AmpGridSchema, FluxGridSchema, validate_schema, schema_version, SCHEMA_VERSION

@testset "dataschema" begin
    @testset "H1Schema round-trip + literal layout (io_hdf5.jl:29-56)" begin
        gs = 6; r0i = 3; nf = 2
        grid = collect(2.0:1.0:(2.0 + gs - 1))            # length 6
        L = reshape(collect(1.0:(4nf * (r0i + 1))), 4nf, r0i + 1)   # (8, 4)
        R = reshape(collect(1.0:(4nf * (gs - r0i))), 4nf, gs - r0i) # (8, 3)
        C_inf = collect(1.0:(2nf)); C_horiz = collect(2.0:(2nf + 1))
        s = H1Schema{Float64}(grid, L, R, C_inf, C_horiz, r0i, 0, gs - 1)
        # field round-trip identity
        @test s.grid == grid
        @test s.inhom_left == L && s.inhom_right == R
        @test s.C_inf == C_inf && s.C_horiz == C_horiz
        @test s.r0index == r0i && s.ra_index == 0 && s.rb_index == gs - 1
        # the literal layout invariants hold ⇒ validates
        @test size(s.inhom_left, 2) == r0i + 1
        @test size(s.inhom_right, 2) == gs - r0i
        @test validate_schema(s) == true
        @test schema_version(s) == SCHEMA_VERSION
    end

    @testset "H1Schema validators reject layout violations" begin
        gs = 6; r0i = 3; nf = 2
        grid = collect(2.0:1.0:(2.0 + gs - 1))
        L = reshape(collect(1.0:(4nf * (r0i + 1))), 4nf, r0i + 1)
        R = reshape(collect(1.0:(4nf * (gs - r0i))), 4nf, gs - r0i)
        Cv = collect(1.0:(2nf))
        # r0index out of range
        @test_throws ArgumentError validate_schema(H1Schema{Float64}(grid, L, R, Cv, Cv, gs, 0, gs - 1))
        # wrong inhom_left column count (r0index+1 violated)
        Lbad = reshape(collect(1.0:(4nf * (r0i + 2))), 4nf, r0i + 2)
        @test_throws DimensionMismatch validate_schema(H1Schema{Float64}(grid, Lbad, R, Cv, Cv, r0i, 0, gs - 1))
        # wrong inhom_right column count (gs−r0index violated)
        Rbad = reshape(collect(1.0:(4nf * (gs - r0i + 1))), 4nf, gs - r0i + 1)
        @test_throws DimensionMismatch validate_schema(H1Schema{Float64}(grid, L, Rbad, Cv, Cv, r0i, 0, gs - 1))
        # mismatched amplitude vectors
        @test_throws DimensionMismatch validate_schema(H1Schema{Float64}(grid, L, R, Cv, collect(1.0:(2nf + 1)), r0i, 0, gs - 1))
        # rows not a multiple of 4
        L5 = reshape(collect(1.0:(5 * (r0i + 1))), 5, r0i + 1)
        R5 = reshape(collect(1.0:(5 * (gs - r0i))), 5, gs - r0i)
        @test_throws DimensionMismatch validate_schema(H1Schema{Float64}(grid, L5, R5, Cv, Cv, r0i, 0, gs - 1))
    end

    @testset "AmpGridSchema axis invariants" begin
        w = [0.0, 0.2, 0.4, 0.6]; u = [0.1, 0.3, 0.5]
        modes = [(2, 2, 0), (3, 2, 1)]
        s = AmpGridSchema{Float64}(w, u, modes)
        @test s.w == w && s.u == u && s.modes == modes
        @test validate_schema(s) == true
        @test schema_version(s) == SCHEMA_VERSION
        @test_throws ArgumentError validate_schema(AmpGridSchema{Float64}([0.0, 0.1], u, modes))      # <3 nodes
        @test_throws ArgumentError validate_schema(AmpGridSchema{Float64}([0.0, 0.0, 0.2], u, modes)) # not increasing
        @test_throws ArgumentError validate_schema(AmpGridSchema{Float64}(w, u, NTuple{3,Int}[]))     # empty modes
    end

    @testset "FluxGridSchema region + tricubic axes" begin
        u = [0.0, 0.25, 0.5, 0.75, 1.0]; w = [0.0, 0.5, 1.0]; z = [0.0, 0.5, 1.0]
        s = FluxGridSchema{Float64}(:A, u, w, z)
        @test s.region == :A && s.u == u && s.w == w && s.z == z
        @test validate_schema(s) == true
        @test validate_schema(FluxGridSchema{Float64}(:B, u, w, z)) == true
        @test_throws ArgumentError validate_schema(FluxGridSchema{Float64}(:C, u, w, z))            # bad region
        @test_throws ArgumentError validate_schema(FluxGridSchema{Float64}(:A, [0.0, 1.0], w, z))   # <3 nodes
        @test_throws ArgumentError validate_schema(FluxGridSchema{Float64}(:A, u, [1.0, 0.5, 0.0], z)) # decreasing
    end

    @testset "BigFloat genericity" begin
        u = big.([0.0, 0.5, 1.0]); w = big.([0.0, 0.5, 1.0]); z = big.([0.0, 0.5, 1.0])
        s = FluxGridSchema{BigFloat}(:A, u, w, z)
        @test eltype(s.u) === BigFloat
        @test validate_schema(s) == true
    end
end

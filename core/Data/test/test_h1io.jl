# Module 4: h1io — versioned H1Schema HDF5 reader/writer (de-dups H1Lorenz/SecondOrderRicci/Punctures/
# DiscEMRIForce readers). GROUNDING:
#  * ROUND-TRIP (primary): build a validate_schema-passing H1Schema -> write_mode_h5 (ext) -> read_h1 ->
#    assert field/byte identity. HDF5 is loaded here to fire the weakdep extension.
#  * h1_parity_fields / hret_block: the LITERAL legacy index tables (SecondOrderRicci read_h1.jl L10-18,
#    Punctures subtraction.jl L23-46), transcribed as expected literals (pure closed form).
#  * If the real run_r8 corpus is present (GRTOOLS_H1LORENZ_ROOT), read h1-l2m2.h5 and assert the
#    on-disk layout; when absent the gate records a Broken (TestKit.oracle_or_skip), never a silent pass.

using Test
using Data
using HDF5   # fires DataHDF5Ext — read_h1/write_mode_h5 bodies activate
using TestKit: oracle_root, oracle_or_skip

@testset "h1io" begin
    @testset "h1_parity_fields literal table" begin
        # even l+m
        @test Data.h1_parity_fields(2, 2) == [1, 3, 5, 6, 7, 2, 4]
        @test Data.h1_parity_fields(3, 1) == [1, 3, 5, 6, 7, 2, 4]
        @test Data.h1_parity_fields(1, 1) == [1, 3, 5, 6, 2, 4]   # dipole
        @test Data.h1_parity_fields(0, 0) == [1, 3, 6, 2]         # monopole
        # odd l+m
        @test Data.h1_parity_fields(2, 1) == [9, 10, 8]
        @test Data.h1_parity_fields(3, 2) == [9, 10, 8]
        @test Data.h1_parity_fields(1, 0) == [9, 8]               # odd dipole
    end

    @testset "hret_block literal map" begin
        # even, l>=2
        @test [Data.hret_block(i, 2) for i in (1, 3, 5, 6, 7, 2, 4)] == [1, 2, 3, 4, 5, 6, 7]
        # even dipole l==1
        @test [Data.hret_block(i, 1) for i in (1, 3, 5, 6, 2, 4)] == [1, 2, 3, 4, 5, 6]
        # monopole l==0
        @test [Data.hret_block(i, 0) for i in (1, 3, 6, 2)] == [1, 2, 3, 4]
        # odd, l>=2
        @test [Data.hret_block(i, 2) for i in (9, 10, 8)] == [1, 2, 3]
        # odd dipole l==1
        @test [Data.hret_block(i, 1) for i in (9, 8)] == [1, 2]
        # absent components throw
        @test_throws ArgumentError Data.hret_block(7, 1)   # field-7 absent at dipole
        @test_throws ArgumentError Data.hret_block(5, 0)   # field-5 absent at monopole
        @test_throws ArgumentError Data.hret_block(10, 1)  # field-10 absent at odd dipole
    end

    @testset "write_mode_h5 / read_h1 round-trip" begin
        gs = 5; r0index = 2; nf = 1
        grid = collect(1.0:gs)
        inhom_left = reshape(collect(1.0:(4nf * (r0index + 1))), 4nf, r0index + 1)
        inhom_right = reshape(collect(100.0 .+ (1.0:(4nf * (gs - r0index)))), 4nf, gs - r0index)
        C_inf = [1.5, -2.5]
        C_horiz = [3.5, 4.25]
        s = H1Schema{Float64}(grid, inhom_left, inhom_right, C_inf, C_horiz, r0index, 0, gs - 1)
        @test validate_schema(s)

        path = joinpath(mktempdir(), "h1-l2m2.h5")
        @test Data.write_mode_h5(path, s) == path
        r = Data.read_h1(path)
        @test r isa H1Schema{Float64}
        @test r.grid == grid
        @test r.inhom_left == inhom_left
        @test r.inhom_right == inhom_right
        @test r.C_inf == C_inf
        @test r.C_horiz == C_horiz
        @test r.r0index == r0index
        @test r.ra_index == 0
        @test r.rb_index == gs - 1
        @test validate_schema(r)

        # BigFloat widening: element type follows T, values preserved.
        rb = Data.read_h1(path; T = BigFloat)
        @test rb isa H1Schema{BigFloat}
        @test eltype(rb.grid) == BigFloat
        @test rb.grid == BigFloat.(grid)
        @test rb.inhom_left == BigFloat.(inhom_left)
    end

    @testset "real run_r8 corpus layout (if present)" begin
        f = joinpath(oracle_root("GRTOOLS_H1LORENZ_ROOT",
                                 "/data/haiyangw/claude/emri-disk/ref-code/h1Lorenz"),
                     "run_r8", "fields", "h1-l2m2.h5")
        if oracle_or_skip(isfile(f), "run_r8 h1 corpus on-disk layout"; env = "GRTOOLS_H1LORENZ_ROOT")
            s = Data.read_h1(f)
            @test length(s.grid) == 6196
            @test s.r0index == 1049
            @test s.ra_index == 0
            @test s.rb_index == 6195
            @test size(s.inhom_left) == (28, 1050)    # (4nf, r0index+1), nf=7 even parity
            @test size(s.inhom_right) == (28, 5147)   # (4nf, gs-r0index)
            @test length(s.C_inf) == 10
            @test length(s.C_horiz) == 10
            @test validate_schema(s)
        end
    end
end

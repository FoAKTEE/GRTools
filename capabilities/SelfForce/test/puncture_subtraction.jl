# Regular-field subtraction hR = hret − hS + retarded-mode HDF5 read (Puncture/Subtraction.jl).
#
# Grounding (CAPTURE-FIXTURE + identity):
#  * hret_block parity/dipole/monopole map: full enumeration of the (consumed Data) block literals.
#  * read_hret_mode column layout (the GATE): an H1Schema fixture is MINTED in-repo with the known
#    encoding _enc(k,p,gi) and written by the canonical Data.write_mode_h5; reading it back through
#    read_hret_mode must recover, for each physical i, the block hret_block(i,l) with rows
#    [Re h, Im h, Re ∂_r h, Im ∂_r h] — the (k-1)*4+part layout the corpus stores.
#  * compute_hR/compute_ddhR: subtraction-linearity identity hR+hS===hret (against the oracle-grounded
#    hS_mode) + a known residual recovered exactly.

using SelfForce.Puncture: hret_block, read_hret_mode, compute_hR, compute_ddhR, hS_mode
import Data
import NumericsCore: H1Schema

# encode block k, part p (1..4), joined-grid index gi (1-based) into a distinct real.
_enc(k, p, gi) = k * 1.0 + p * 1.0e-3 + gi * 1.0e-6

# Mint an H1Schema with the _enc encoding: nf fields (4nf rows), gs grid points, worldline at
# 0-based r0index; the joined inhom_left ++ inhom_right carries full[(k-1)*4+p, gi] = _enc(k,p,gi).
function _mint_schema(nf, gs, r0index)
    njoined = gs + 1                              # worldline column doubled
    full = zeros(Float64, 4nf, njoined)
    for k in 1:nf, p in 1:4, gi in 1:njoined
        full[(k - 1) * 4 + p, gi] = _enc(k, p, gi)
    end
    grid = collect(range(8.0, 12.0; length = gs))
    left = full[:, 1:(r0index + 1)]
    right = full[:, (r0index + 2):end]
    return H1Schema{Float64}(grid, left, right, Float64[], Float64[],
                             r0index, 0, gs - 1)
end

@testset "Subtraction (hR, hret read, column map)" begin
    @testset "hret_block map (full enumeration vs consumed Data literals)" begin
        # even-parity general (l>=2): R1,R3,R5,R6,R7 coupled, R2,R4 gauge.
        @test hret_block(1, 2) == 1
        @test hret_block(3, 2) == 2
        @test hret_block(5, 2) == 3
        @test hret_block(6, 2) == 4
        @test hret_block(7, 2) == 5
        @test hret_block(2, 2) == 6
        @test hret_block(4, 2) == 7
        # even-dipole (l=1): no R7; gauge shifts down.
        @test hret_block(2, 1) == 5
        @test hret_block(4, 1) == 6
        @test_throws ArgumentError hret_block(7, 1)
        # monopole (l=0).
        @test hret_block(6, 0) == 3
        @test hret_block(2, 0) == 4
        # odd-parity general (l>=2): R9,R10 coupled, R8 gauge.
        @test hret_block(9, 2) == 1
        @test hret_block(10, 2) == 2
        @test hret_block(8, 2) == 3
        # odd-dipole (l=1): R9 coupled, R8 gauge.
        @test hret_block(9, 1) == 1
        @test hret_block(8, 1) == 2
    end

    @testset "read_hret_mode round-trip vs write_mode_h5 (GATE: column layout)" begin
        mktempdir() do dir
            fe = joinpath(dir, "h1-l2m2.h5")
            Data.write_mode_h5(fe, _mint_schema(7, 8, 3))        # even file (l=2): 7 fields
            for i in (1, 3, 5, 6, 7, 2, 4)
                grid, hret, dhret = read_hret_mode(fe, i, 2)
                k = hret_block(i, 2)
                @test length(grid) == 8
                @test length(hret) == 9                          # joined: gs+1
                @test real(hret[1]) ≈ _enc(k, 1, 1)
                @test imag(hret[1]) ≈ _enc(k, 2, 1)
                @test real(dhret[1]) ≈ _enc(k, 3, 1)
                @test imag(dhret[1]) ≈ _enc(k, 4, 1)
            end
            fo = joinpath(dir, "h1-l3m2.h5")
            Data.write_mode_h5(fo, _mint_schema(3, 8, 3))        # odd file (l=3): 3 fields
            for i in (9, 10, 8)
                _, hret, dhret = read_hret_mode(fo, i, 3)
                k = hret_block(i, 3)
                @test real(hret[1]) ≈ _enc(k, 1, 1)
                @test imag(hret[1]) ≈ _enc(k, 2, 1)
                @test real(dhret[1]) ≈ _enc(k, 3, 1)
                @test imag(dhret[1]) ≈ _enc(k, 4, 1)
            end
        end
    end

    @testset "exact relation hR = hret − hS vs oracle-grounded hS (L485)" begin
        grid = [-1.0, -0.5, 0.5, 1.0]
        for (i, l, m) in ((1, 2, 2), (6, 3, 1), (8, 2, 1))
            hS = hS_mode(i, l, m, 10.0, 1.0, grid)
            residual = ComplexF64[0.1k - 0.2k * im for k in 1:length(grid)]
            hret = hS .+ residual
            hR = compute_hR(hret, hS)
            @test hR ≈ residual
            @test hR .+ hS ≈ hret
        end
        @test_throws DimensionMismatch compute_hR([1.0 + 0im], [1.0 + 0im, 2.0])
    end

    @testset "compute_ddhR exact-relation primitive (L515)" begin
        a = ComplexF64[1, 2 + im, 3 - 2im]
        b = ComplexF64[0.5, 1 + 0.5im, 1 - im]
        @test compute_ddhR(a, b) ≈ a .- b
        @test compute_ddhR(a, b) .+ b ≈ a
        @test_throws DimensionMismatch compute_ddhR(a, b[1:2])
    end
end

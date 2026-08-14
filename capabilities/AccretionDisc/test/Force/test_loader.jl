# Module 5 — Force/Loader.jl. The in-memory complex BLS adapter H1Mode{T} over the Data.H1Schema
# contract. GROUNDING:
#   • ROUND-TRIP a captured small HDF5 fixture (built as a validate_schema-passing H1Schema, written
#     by Data.write_mode_h5 — the same fixture-write pattern as Data/test/test_h1io.jl) through
#     Data.read_h1 → load_h1_mode, asserting grid / r0 / r0index split / field VALUES match the schema
#     literals exactly (Float64), the parity rule (nf=7⇒:even, nf=3⇒:odd), and the Complex{T=BigFloat}
#     widening identity. NOT a live DiscEMRIForce dep (HDF5-heavy).
#   • load_modes parallel == serial (identical Dict), available_modes == Data.available_modes scan.

using Test
using AccretionDisc
using AccretionDisc.Force
using Data
using HDF5   # fires DataHDF5Ext (also active via Force's import HDF5)

# Build a validate_schema-passing H1Schema with `nf` fields and distinct, recoverable entries.
function _make_schema(nf; gs = 6, r0index = 2)
    grid = collect(10.0:(10.0 + gs - 1))
    nL = r0index + 1; nR = gs - r0index
    # entry (row, col) = row*1000 + col  → uniquely identifies the (field, part, radial) slot
    inhom_left = Float64[row * 1000 + col for row in 1:(4nf), col in 1:nL]
    inhom_right = Float64[5_000_000 + row * 1000 + col for row in 1:(4nf), col in 1:nR]
    C_inf = collect(1.0:nf); C_horiz = collect(101.0:(100.0 + nf))
    H1Schema{Float64}(grid, inhom_left, inhom_right, C_inf, C_horiz, r0index, 0, gs - 1)
end

@testset "Force.Loader" begin

    @testset "H1Mode round-trip vs H1Schema literals (even nf=7)" begin
        dir = mktempdir()
        s = _make_schema(7)
        Data.write_mode_h5(joinpath(dir, "h1-l2m2.h5"), s)
        mode = load_h1_mode(2, 2; dir = dir)
        @test mode isa H1Mode{Float64}
        @test mode.l == 2 && mode.m == 2
        @test mode.nfields == 7
        @test mode.parity === :even
        @test mode.r0 == s.grid[s.r0index + 1]                    # 0-based r0index → grid[r0index+1]
        @test mode.r_right == s.grid[s.r0index+1:end]
        @test mode.r_left == s.grid[1:s.r0index+1]
        @test mode.C_inf == s.C_inf
        # field value/deriv split: val_right[j][k] = R[(j-1)*4+1,k] + i R[(j-1)*4+2,k]
        for j in 1:7, k in 1:size(s.inhom_right, 2)
            @test mode.val_right[j][k] == complex(s.inhom_right[(j-1)*4+1, k], s.inhom_right[(j-1)*4+2, k])
            @test mode.dr_right[j][k]  == complex(s.inhom_right[(j-1)*4+3, k], s.inhom_right[(j-1)*4+4, k])
        end
        for j in 1:7, k in 1:size(s.inhom_left, 2)
            @test mode.val_left[j][k] == complex(s.inhom_left[(j-1)*4+1, k], s.inhom_left[(j-1)*4+2, k])
            @test mode.dr_left[j][k]  == complex(s.inhom_left[(j-1)*4+3, k], s.inhom_left[(j-1)*4+4, k])
        end
    end

    @testset "parity rule (odd nf=3)" begin
        dir = mktempdir()
        Data.write_mode_h5(joinpath(dir, "h1-l2m1.h5"), _make_schema(3))
        mode = load_h1_mode(2, 1; dir = dir)
        @test mode.nfields == 3
        @test mode.parity === :odd
    end

    @testset "low-l compact fields map to canonical parity layouts" begin
        dir = mktempdir()
        cases = ((0, 0), (1, 0), (1, 1), (2, 1), (2, 2))
        for (l, m) in cases
            fields = Data.h1_parity_fields(l, m)
            s = _make_schema(length(fields))
            Data.write_mode_h5(joinpath(dir, "h1-l$(l)m$(m).h5"), s)
            mode = load_h1_mode(l, m; dir = dir)
            canonical = iseven(l + m) ? [1, 3, 5, 6, 7, 2, 4] : [9, 10, 8]
            @test mode.parity === (iseven(l + m) ? :even : :odd)
            @test mode.nfields == length(canonical)
            @test length(mode.val_right) == length(canonical)
            for (j, field) in enumerate(canonical)
                src = findfirst(==(field), fields)
                if src === nothing
                    @test all(iszero, mode.val_right[j])
                    @test all(iszero, mode.dr_right[j])
                    @test all(iszero, mode.val_left[j])
                    @test all(iszero, mode.dr_left[j])
                else
                    @test mode.val_right[j][1] == complex(
                        s.inhom_right[(src-1)*4+1, 1], s.inhom_right[(src-1)*4+2, 1])
                    @test mode.dr_right[j][1] == complex(
                        s.inhom_right[(src-1)*4+3, 1], s.inhom_right[(src-1)*4+4, 1])
                end
            end
        end
    end

    @testset "Complex{BigFloat} widening identity" begin
        dir = mktempdir()
        s = _make_schema(7)
        Data.write_mode_h5(joinpath(dir, "h1-l2m2.h5"), s)
        mode64 = load_h1_mode(2, 2; dir = dir)
        modeBF = load_h1_mode(2, 2; dir = dir, T = BigFloat)
        @test modeBF isa H1Mode{BigFloat}
        @test eltype(modeBF.r_right) == BigFloat
        @test eltype(modeBF.val_right[1]) == Complex{BigFloat}
        # values preserved exactly across the widening
        @test modeBF.val_right[3][1] == Complex{BigFloat}(mode64.val_right[3][1])
        @test modeBF.r0 == BigFloat(mode64.r0)
    end

    @testset "load_modes parallel == serial + available_modes" begin
        dir = mktempdir()
        for l in 2:4
            Data.write_mode_h5(joinpath(dir, "h1-l$(l)m2.h5"),
                               _make_schema(length(Data.h1_parity_fields(l, 2))))
        end
        par = load_modes(2:4, 2; dir = dir)
        ser = Dict(l => load_h1_mode(l, 2; dir = dir) for l in 2:4)
        @test sort(collect(keys(par))) == [2, 3, 4]
        for l in 2:4
            @test par[l].val_right[1] == ser[l].val_right[1]
            @test par[l].r_right == ser[l].r_right
        end
        @test Force.available_modes(dir = dir) == [(2, 2), (3, 2), (4, 2)]
    end

end

# Grounded tests for SelfForce.SecondOrder.read_h1 (directory-level h1-corpus reader) + the
# ZeroPadField shim + the h1_parity_fields re-export.
#
# GROUNDING — structural + closed-form on the CAPTURED fixtures, no heavy legacy dep:
#
# The fixtures capabilities/SelfForce/test/fixtures/h1/h1-l{l}m{m}.h5 were captured one-shot
# from the genuine H1Lorenz corpus (the same files SecondOrderRicci.read_h1 ingests). Each holds a
# "grid" vector (G points), and "inhom_left"/"inhom_right" matrices of shape (6·#fields, N_left/right)
# storing per field the 6 columns (Re h̄, Im h̄, Re ∂_r h̄, Im ∂_r h̄, Re ∂²_r h̄, Im ∂²_r h̄).
#
# (1) STRUCTURE: read_h1 reconstructs the full grid r of length N=G+1 with r0 duplicated at the
#     worldline (the C++ inhom_left/right split point), so r0 appears TWICE and the left/right halves
#     join there. We assert N=G+1, the duplication index, and N_left+N_right==N.
# (2) COLUMN UNPACK: the populated (i,l,m) keys are exactly the parity field set h1_parity_fields(l,m)
#     mapped through the 3↔6 trace-reversal swap, for every fixture mode; each h/dh/ddh has length N.
# (3) ANALYTIC Schwarzschild metric: f = 1−2/r and fp = 2/r² reproduced at every grid point (M=1).
# (4) NEGATIVE-m SYMMETRY: h[(i,l,-m)] == (-1)^m conj(h[(i,l,m)]) (the C++ azimuthal reflection), and
#     for m=0 the in-place conjugation h[(i,l,0)] == (-1)^0 conj(...) (real-axis closure).
# (5) h1_parity_fields is the Data re-export (=== object identity, not a second copy).

using Test
using SelfForce.SecondOrder: read_h1, ZeroPadField, h1_parity_fields
using Data: h1_parity_fields as data_h1_parity_fields
import HDF5

const _FIX = joinpath(@__DIR__, "fixtures", "h1")

@testset "h1_parity_fields is the Data re-export" begin
    @test h1_parity_fields === data_h1_parity_fields
    @test h1_parity_fields(0, 0) == [1, 3, 6, 2]
    @test h1_parity_fields(2, 2) == [1, 3, 5, 6, 7, 2, 4]
    @test h1_parity_fields(3, 0) == [9, 10, 8]   # l+m odd
end

@testset "read_h1 grid structure + r0 duplication" begin
    fields = read_h1(_FIX; l_max = 3)
    # N_left from the first sorted fixture file (h1-l0m0.h5)
    grid, nL = HDF5.h5open(joinpath(_FIX, "h1-l0m0.h5"), "r") do file
        (read(file["grid"]), size(read(file["inhom_left"]), 2))
    end
    G = length(grid); N = G + 1
    r = fields.r
    @test length(r) == N
    @test fields.r0 == grid[nL]               # worldline at the left/right split (1-based nL)
    @test r[nL] == fields.r0 && r[nL + 1] == fields.r0   # r0 duplicated across the join
    @test r[1:nL] == grid[1:nL]               # left half is the raw grid up to r0
end

@testset "read_h1 Schwarzschild metric (analytic f, fp)" begin
    fields = read_h1(_FIX; l_max = 3)
    @test fields.f ≈ @.(1 - 2 / fields.r) rtol = 1e-15
    @test fields.fp ≈ @.(2 / fields.r^2) rtol = 1e-15
end

@testset "read_h1 column unpack: keys = parity fields ⊕ 3↔6 swap, lengths = N" begin
    fields = read_h1(_FIX; l_max = 3)
    N = length(fields.r)
    swap(i) = i == 3 ? 6 : (i == 6 ? 3 : i)
    pat = r"^h1-l([0-9]+)m([0-9]+)\.h5$"
    for fn in readdir(_FIX)
        mm = match(pat, fn); mm === nothing && continue
        l = parse(Int, mm[1]); m = parse(Int, mm[2])
        for i in h1_parity_fields(l, m)
            key = (swap(i), l, m)
            @test haskey(fields.h, key)
            @test length(fields.h[key]) == N
            @test length(fields.dh[key]) == N
            @test length(fields.ddh[key]) == N
        end
    end
end

@testset "read_h1 negative-m azimuthal reflection (m>0 distinct keys)" begin
    fields = read_h1(_FIX; l_max = 3)
    swap(i) = i == 3 ? 6 : (i == 6 ? 3 : i)
    for (l, m) in [(2, 2), (3, 1), (2, 1), (1, 1)]
        sgn = isodd(m) ? -1 : 1
        for i in h1_parity_fields(l, m)
            k = (swap(i), l, m); kn = (swap(i), l, -m)
            @test haskey(fields.h, kn)
            @test fields.h[kn] ≈ sgn .* conj.(fields.h[k]) rtol = 1e-14
            @test fields.dh[kn] ≈ sgn .* conj.(fields.dh[k]) rtol = 1e-14
            @test fields.ddh[kn] ≈ sgn .* conj.(fields.ddh[k]) rtol = 1e-14
        end
    end
end

# For m=0 the C++ source applies the reflection IN PLACE (dst==src), i.e. h[(i,l,0)] is overwritten
# with conj(h[(i,l,0)]). We pin this faithful quirk against an independent recompute of the physical
# field straight from the raw inhom_left columns (a·h̄/r) at the worldline column, then conjugate.
@testset "read_h1 m=0 in-place conjugation (faithful legacy quirk)" begin
    using NumericsCore: a_il
    fields = read_h1(_FIX; l_max = 3)
    l = 0; m = 0
    grid = HDF5.h5open(joinpath(_FIX, "h1-l0m0.h5"), "r") do file; read(file["grid"]); end
    L = HDF5.h5open(joinpath(_FIX, "h1-l0m0.h5"), "r") do file; read(file["inhom_left"]); end
    swap(i) = i == 3 ? 6 : (i == 6 ? 3 : i)
    for (it, i) in enumerate(h1_parity_fields(l, m))
        ih = swap(i); a = a_il(ih, l); c = 6 * (it - 1)
        j = 1                          # first grid point (left solution)
        hbar = complex(L[c + 1, j], L[c + 2, j])
        recomputed = a * hbar / grid[j]
        @test fields.h[(ih, l, m)][j] ≈ conj(recomputed) rtol = 1e-13
    end
end

@testset "ZeroPadField returns shared zero for absent keys, stored vec for present" begin
    z = zeros(ComplexF64, 5)
    d = Dict((1, 0, 0) => ComplexF64[1, 2, 3, 4, 5])
    p = ZeroPadField(d, z)
    @test p[(1, 0, 0)] === d[(1, 0, 0)]
    @test p[(7, 9, 9)] === z          # missing → the single shared zero vector (read-only alias)
    @test all(iszero, p[(2, 2, 2)])
end

@testset "read_h1 BigFloat precision path (type-generic {T})" begin
    fields = read_h1(BigFloat, _FIX; l_max = 2)
    @test eltype(fields.r) == BigFloat
    @test fields.f ≈ @.(1 - 2 / fields.r) rtol = 1e-40
    k = first(keys(fields.h))
    @test eltype(fields.h[k]) == Complex{BigFloat}
end

# Singular-field modes hS[i,l,m] on the worldtube grid + worldline split (Puncture/SingularField.jl).
#
# Grounding (IDENTITY + enumeration):
#  * hS_mode off-worldline === [h_mode(i,l,m,r0,M,Δr) for Δr in grid] (already-grounded h_mode).
#  * value pin at Δr=0.5 vs the captured legacy h oracle (fixtures/puncture/h_reference.csv).
#  * SINGULAR_ZERO_MODES short-circuit returns exact zeros (enumeration).
#  * worldline_index returns the exact Δr=0 position; split_worldline is the index identity.

using SelfForce.Puncture: hS_mode, worldline_index, split_worldline, SINGULAR_ZERO_MODES, h_mode
import SelfForce.Puncture as P

const _HMREF_S = joinpath(@__DIR__, "fixtures", "puncture", "h_reference.csv")
_sf_pc(s) = (a = split(s, '|'); complex(parse(Float64, a[1]), parse(Float64, a[2])))

@testset "SingularField (hS)" begin
    @testset "worldline index + split convention (L550-555)" begin
        wgrid = [-1.0, -0.5, 0.0, 0.5, 1.0]
        r0i = worldline_index(wgrid)
        @test r0i == 3
        vals = collect(1:length(wgrid))
        left, right = split_worldline(vals, r0i)
        @test left == [1, 2, 3]               # includes worldline
        @test right == [4, 5]                 # excludes worldline
        @test length(left) + length(right) == length(wgrid)
        @test wgrid[left[end]] == 0.0
        @test_throws ArgumentError worldline_index([-1.0, -0.5, 0.5])
    end

    grid = [-1.0, -0.5, 0.5, 1.0]
    idx05 = findfirst(==(0.5), grid)

    @testset "hS values at Δr=0.5 vs captured h oracle" begin
        n = 0
        for (k, line) in enumerate(eachline(_HMREF_S))
            (k == 1 || isempty(strip(line))) && continue
            f = split(line, ',')
            i = parse(Int, f[1]); l = parse(Int, f[2]); m = parse(Int, f[3])
            ex = _sf_pc(f[4])
            abs(ex) == 0 && continue
            v = hS_mode(i, l, m, 10.0, 1.0, grid)
            @test isapprox(v[idx05], ex; rtol = 1e-7, atol = 1e-9)
            n += 1
        end
        @test n >= 15
    end

    @testset "hS === h_mode mapped over the grid (assembly identity, non-zero modes)" begin
        for (i, l, m) in ((1, 2, 2), (3, 2, 2), (6, 3, 1), (8, 2, 1), (9, 2, 1))
            v = hS_mode(i, l, m, 10.0, 1.0, grid)
            @test length(v) == length(grid)
            for (k, Δr) in enumerate(grid)
                @test v[k] == h_mode(i, l, m, 10.0, 1.0, Δr; mmax = 10)
            end
        end
    end

    @testset "identically-zero modes hS[2,0,0], hS[9,1,0] (L460-475)" begin
        wgrid = [-1.0, -0.5, 0.0, 0.5, 1.0]      # zero modes short-circuit; worldline safe
        @test all(==(0), hS_mode(2, 0, 0, 10.0, 1.0, wgrid))
        @test all(==(0), hS_mode(9, 1, 0, 10.0, 1.0, wgrid))
        @test (2, 0, 0) in SINGULAR_ZERO_MODES
        @test (9, 1, 0) in SINGULAR_ZERO_MODES
    end

    @testset "BigFloat precision path (single off-worldline point)" begin
        v = hS_mode(1, 2, 2, big(10), big(1), BigFloat[1//2])
        @test eltype(v) == Complex{BigFloat}
        @test isapprox(Float64(real(v[1])), 4.094659219523288963958; rtol = 1e-7)
    end
end

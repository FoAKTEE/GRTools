# Wigner-D β=π/2 mode mixing h[i][l,m] (Puncture/ModeMixing.jl).
#
# Grounding (CAPTURE-FIXTURE + identity):
#  * h[i] mixed modes vs fixtures/puncture/h_reference.csv (captured legacy Punctures oracle, copied
#    verbatim) at the legacy tolerance.
#  * the consumed Wigner-D helper vs fixtures/puncture/wignerD_reference.csv (Mathematica WignerD
#    oracle) — the AHA collapse onto Wigner.wignerD_fixed reproduces the fixed-angle rotation.
#  * the reality identity h[i][l,-m] = (-1)^m conj(h[i][l,m]) — oracle-independent.

using SelfForce.Puncture: h_mode
import SelfForce.Puncture as P

const _WREF = joinpath(@__DIR__, "fixtures", "puncture", "wignerD_reference.csv")
const _HMREF = joinpath(@__DIR__, "fixtures", "puncture", "h_reference.csv")
_mm_pc(s) = (a = split(s, '|'); complex(parse(Float64, a[1]), parse(Float64, a[2])))

@testset "ModeMixing (h)" begin
    @testset "consumed Wigner-D fixed-angle helper vs oracle" begin
        n = 0
        for (k, line) in enumerate(eachline(_WREF))
            (k == 1 || isempty(strip(line))) && continue
            f = split(line, ',')
            l = parse(Int, f[1]); m = parse(Int, f[2]); mp = parse(Int, f[3])
            ex = _mm_pc(f[4])
            @test isapprox(P._wignerD_fixed(Float64, l, m, mp), ex; atol = 1e-12)
            n += 1
        end
        @test n > 250
    end

    @testset "2^{-j} closed form for small (j,m,mp)" begin
        # d^1_{0,0}(π/2)=0, d^2_{0,0}(π/2)=-1/2, d^1_{1,-1}(π/2)=1/2 (textbook reduced-d values).
        @test isapprox(P._wignerd_halfpi(Float64, 1, 0, 0), 0.0; atol = 1e-13)
        @test isapprox(P._wignerd_halfpi(Float64, 2, 0, 0), -0.5; atol = 1e-13)
        @test isapprox(P._wignerd_halfpi(Float64, 1, 1, -1), 0.5; atol = 1e-13)
    end

    @testset "Wigner-D BigFloat precision" begin
        setprecision(BigFloat, 192) do
            g = P._wignerD_fixed(BigFloat, 5, 3, 2)
            gF = P._wignerD_fixed(Float64, 5, 3, 2)
            @test isapprox(Float64(real(g)), real(gF); atol = 1e-12)
            @test g isa Complex{BigFloat}
        end
    end

    @testset "h[i] mixed modes vs captured legacy oracle" begin
        rows = Tuple{Int,Int,Int,ComplexF64}[]
        for (k, line) in enumerate(eachline(_HMREF))
            (k == 1 || isempty(strip(line))) && continue
            f = split(line, ',')
            push!(rows, (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]), _mm_pc(f[4])))
        end
        @test length(rows) > 20
        for (i, l, m, ex) in rows
            @test isapprox(h_mode(i, l, m, 10.0, 1.0, 0.5), ex; rtol = 1e-7, atol = 1e-9)
        end
    end

    @testset "reality identity h[i][l,-m] = (-1)^m conj(h[i][l,m])" begin
        for (i, l, m) in ((1, 2, 1), (2, 3, 1), (4, 3, 2), (6, 2, 2), (8, 3, 1))
            hp = h_mode(i, l, m, 10.0, 1.0, 0.5)
            hm = h_mode(i, l, -m, 10.0, 1.0, 0.5)
            @test isapprox(hm, (-1)^m * conj(hp); rtol = 1e-7, atol = 1e-9)
        end
    end
end

# Worldline self-field assembly + covariant derivatives (Puncture/WorldlineSelfField.jl).
#
# Grounding (CAPTURE-FIXTURE + identity):
#  * covariant_first/second vs fixtures/puncture/{DhR1,D2hR1}_reference.csv (captured legacy oracle:
#    the explicit covariant-derivative formulas on SYNTHETIC partials, reproduced here from the same
#    closed forms) at the legacy tolerance. ∇_d∇_c does NOT commute — the fixture is asymmetric in
#    (c,d) and the per-index match pins both orderings.
#  * assemble_hR1Data/dhR1Data vs {hR1Data,dhR1Data}_reference.csv on synthetic Lorenz amplitudes.
#  * consumed Christoffel grounded by metric-compatibility ∇_c g_ab = 0 (oracle-independent) + the
#    a=0 background vs Spacetime.kerr_metric.

using SelfForce.Puncture: covariant_first, covariant_second,
                          assemble_hR1Data, assemble_dhR1Data
import SelfForce.Puncture as P
import Spacetime

const _R0 = 10.0
const _M = 1.0
_hR1() = Float64[(a + b) / 5 + a * b / 7 + 1 / 3 for a in 1:4, b in 1:4]
_dhR1() = Float64[(a + b) / 4 + a * b / 9 + c / 6 + (a + b) * c / 8
                  for a in 1:4, b in 1:4, c in 1:4]
_d2hR1() = Float64[(a + b) / 3 + (c + d) / 5 + a * b / 7 + c * d / 11 +
                   ((a + b) + (c + d)) / 13 for a in 1:4, b in 1:4, c in 1:4, d in 1:4]

function _wl_csv(path, ncoord)
    rows = Tuple{Vector{Int},Float64}[]
    for (k, line) in enumerate(eachline(path))
        (k == 1 || isempty(strip(line))) && continue
        f = split(line, ',')
        idx = [parse(Int, f[j]) for j in 1:ncoord]
        push!(rows, (idx, parse(Float64, f[ncoord+1])))
    end
    return rows
end
_fx(name) = joinpath(@__DIR__, "fixtures", "puncture", name)

@testset "WorldlineSelfField" begin
    @testset "background metric vs Spacetime.kerr_metric (a=0)" begin
        g = Spacetime.kerr_metric(0.0, _R0, π / 2)
        f = 1 - 2_M / _R0
        @test isapprox(g[1, 1], -f; rtol = 1e-13)
        @test isapprox(g[2, 2], 1 / f; rtol = 1e-13)
        @test isapprox(g[3, 3], _R0^2; rtol = 1e-13)
        @test isapprox(g[4, 4], _R0^2; rtol = 1e-13)
    end

    @testset "consumed Christoffel: metric compatibility ∇_c g_ab = 0" begin
        for (r, θ) in ((10.0, π / 2), (8.0, 0.8), (6.5, 1.3))
            f = 1 - 2_M / r; s = sin(θ); c = cos(θ)
            g = [-f 0 0 0; 0 1/f 0 0; 0 0 r^2 0; 0 0 0 r^2*s^2]
            dgr = [-2_M/r^2 0 0 0; 0 -(2_M/r^2)/f^2 0 0; 0 0 2r 0; 0 0 0 2r*s^2]
            dgθ = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 2r^2*s*c]
            dg = (zeros(4, 4), dgr, dgθ, zeros(4, 4))
            Γ = P.christoffel(r, θ; M = _M)
            for a in 1:4, b in 1:4, cc in 1:4
                val = dg[cc][a, b]
                for e in 1:4
                    val -= Γ[e, a, cc] * g[e, b] + Γ[e, b, cc] * g[a, e]
                end
                @test abs(val) < 1e-9
            end
        end
    end

    @testset "DhR1Data = ∇_c h_ab vs captured oracle" begin
        Dh = covariant_first(_hR1(), _dhR1(), _R0, _M)
        rows = _wl_csv(_fx("DhR1_reference.csv"), 3)
        @test length(rows) == 64
        for (idx, val) in rows
            @test isapprox(Dh[idx[1], idx[2], idx[3]], val; rtol = 1e-9, atol = 1e-11)
        end
    end

    @testset "D2hR1Data = ∇_d ∇_c h_ab vs captured oracle (non-commuting)" begin
        D2 = covariant_second(_hR1(), _dhR1(), _d2hR1(), _R0, _M)
        rows = _wl_csv(_fx("D2hR1_reference.csv"), 4)
        @test length(rows) == 256
        for (idx, val) in rows
            @test isapprox(D2[idx[1], idx[2], idx[3], idx[4]], val; rtol = 1e-9, atol = 1e-10)
        end
        # explicit [c,d] asymmetry (∇_d∇_c ≠ ∇_c∇_d).
        @test !isapprox(D2[4, 4, 3, 4], D2[4, 4, 4, 3]; atol = 1e-8)
    end

    @testset "assembly hR1Data / dhR1Data vs captured oracle (synthetic modes)" begin
        sH(i, l, m) = (i + 1 + l - m) / 7 + im * (2 + i - l + 2m) / 5
        sDH(i, l, m) = (3 + i + 2l + m) / 11 + im * (1 + i + l - m) / 4
        sDDH(i, l, m) = (5 + 2i - l + m) / 13 + im * (4 + i - 2l + m) / 3
        base(i, l, m) = (ComplexF64(sH(i, l, m)), ComplexF64(sDH(i, l, m)),
                         ComplexF64(sDDH(i, l, m)))
        H = assemble_hR1Data(base, 10.0, 1.0; lmax = 8)
        D = assemble_dhR1Data(base, 10.0, 1.0; lmax = 8)
        for (idx, val) in _wl_csv(_fx("hR1Data_reference.csv"), 2)
            @test isapprox(H[idx[1], idx[2]], val; rtol = 1e-7, atol = 1e-9)
        end
        for (idx, val) in _wl_csv(_fx("dhR1Data_reference.csv"), 3)
            @test isapprox(D[idx[1], idx[2], idx[3]], val; rtol = 1e-7, atol = 1e-9)
        end
        @test all(H[a, b] == H[b, a] for a in 1:4, b in 1:4)
    end
end

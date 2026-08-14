# Grounding for Lorenz.Operators — geometry/tetrad/Chandrasekhar scalars + jet algebra + CG mixing.
#
# (1) ANALYTIC: geom/Kfn/Qfn closed forms + ginv_tetrad/kerr_metric_inv_eq == Spacetime kerr_metric_inv
#     (the already-grounded LIGHT oracle) at θ=π/2 (rtol 1e-12).
# (2) cg == Wigner.clebsch_gordan (the already-built oracle); frak_* are deterministic CG products.
# (3) ONE-SHOT KLC literals (captured via julia --project=packages/KerrLorenzCirc scratchpad/cap12.jl)
#     for frak_s/frak_c/frak_c2 + jet_ops mul/inv at sample (r0,a,ω,m).

using Test
import MetricReconstruction.Lorenz as L
import Spacetime
import Wigner

@testset "Lorenz.Operators (geometry + jet algebra + CG mixing)" begin

    @testset "geometry/tetrad scalars vs Spacetime closed forms (analytic)" begin
        for (a, r) in ((0.0, 6.0), (0.99, 6.0), (0.5, 8.3))
            θ = π/2
            g = L.geom(a, r, θ)
            @test isapprox(g.Δ, Spacetime.delta_kerr(a, r); rtol=1e-12)
            @test isapprox(g.Σ, Spacetime.sigma_kerr(a, r, θ); rtol=1e-12)
            @test g.Δp == 2*(r - 1)
            @test isapprox(L.Kfn(a, r, 0.4, 2), Spacetime.kfun(a, r, 0.4, 2); rtol=1e-12)
            @test isapprox(L.Qfn(a, θ, 0.4, 2), Spacetime.qfun(a, 0.4, 2, θ); rtol=1e-12)
            # inverse metric from tetrad reconstruction == Spacetime kerr_metric_inv
            gi_t = L.ginv_tetrad(a, r, θ)
            gi_e = L.kerr_metric_inv_eq(a, r)
            gi_ref = Spacetime.kerr_metric_inv(a, r, θ)
            for I in CartesianIndices(gi_ref)
                @test isapprox(real(gi_t[I]), gi_ref[I]; rtol=1e-11, atol=1e-13)
                @test isapprox(gi_e[I], gi_ref[I]; rtol=1e-11, atol=1e-13)
            end
        end
        # a=0 Schwarzschild closed form: g^{tt}=-1/f, g^{rr}=Δ/Σ=f, g^{θθ}=1/r², g^{tφ}=0
        gi = L.kerr_metric_inv_eq(0.0, 6.0); f = 1 - 2/6.0
        @test isapprox(gi[1,1], -1/f; rtol=1e-12)
        @test isapprox(gi[2,2], f; rtol=1e-12)
        @test isapprox(gi[3,3], 1/36.0; rtol=1e-12)
        @test abs(gi[1,4]) < 1e-14
        # a=0.99 g^{rr}=Δ/Σ at the equator
        gi99 = L.kerr_metric_inv_eq(0.99, 6.0)
        @test isapprox(gi99[2,2], (6.0^2 - 12 + 0.99^2)/36.0; rtol=1e-12)
    end

    @testset "cg == Wigner.clebsch_gordan (3j↔CG identity)" begin
        for (j1,m1,j2,m2,j3,m3) in ((2,0,1,0,3,0),(2,-1,1,1,2,0),(4,2,2,0,4,2),(3,1,1,-1,2,0))
            @test L.cg(j1,m1,j2,m2,j3,m3) == Wigner.clebsch_gordan(j1,m1,j2,m2,j3,m3)
        end
    end

    @testset "frak_* CG mixing — KLC one-shot literals (m=2)" begin
        @test isapprox(L.frak_s(1,0,2,3,2),  -0.30860669992418377; rtol=1e-13)
        @test isapprox(L.frak_s(-1,0,3,2,2), -0.4364357804719848;  rtol=1e-13)
        @test isapprox(L.frak_s(2,1,2,2,2),   0.6666666666666666;  rtol=1e-13)
        @test isapprox(L.frak_s(-2,-1,3,3,2), 0.5270462766947299;  rtol=1e-13)
        @test isapprox(L.frak_s(-1,1,2,3,2),  0.5345224838248487;  rtol=1e-13)
        @test isapprox(L.frak_s2(2,0,2,4,2),  0.06734350297014736; rtol=1e-13)
        @test isapprox(L.frak_s2(-2,0,4,2,2), 0.26082026547865045; rtol=1e-13)
        @test isapprox(L.frak_c(2,2,3,2),     0.2817180849095055;  rtol=1e-13)
        @test isapprox(L.frak_c(0,3,4,2),     0.4364357804719848;  rtol=1e-13)
        @test isapprox(L.frak_c(-1,4,5,2),    0.451260859854213;   rtol=1e-13)
        @test isapprox(L.frak_c2(0,2,2,2),    0.14285714285714282; rtol=1e-13)
        @test isapprox(L.frak_c2(0,2,4,2),    0.1649572197684645;  rtol=1e-13)
        @test isapprox(L.frak_c2(2,3,3,2),    0.3333333333333333;  rtol=1e-13)
        # ladder diagonals
        @test isapprox(L.Λ̂(2), L.λ̂(2)*L.λ̂2(2); rtol=1e-14)
        @test isapprox(L.λ̂(2), sqrt(6.0); rtol=1e-14)
    end

    @testset "jet algebra: Leibniz/reciprocal identities + KLC jet_ops literals" begin
        C = ComplexF64
        r = 2.3
        u = C[r^3, 3r^2, 6r, 6, 0]; v = C[r^2, 2r, 2, 0, 0]
        @test maximum(abs.(L.jmul(u, v) .- C[r^5, 5r^4, 20r^3, 60r^2, 120r])) < 1e-10
        iv = L.jinv(v)
        @test maximum(abs.(L.jmul(v, iv) .- C[1,0,0,0,0])) < 1e-12
        @test L.jderiv(u) == C[3r^2, 6r, 6, 0, 0]

        # KLC jmul/jinv one-shot literals
        f = C[1.0+0im, 2.0, -1.0, 0.5, 0,0,0]; g = C[2.0+1im, -1.0, 0.3, 0,0,0,0]
        @test L.jmul(f,g) ≈ C[2.0+1.0im, 3.0+2.0im, -5.7-1.0im, 5.8+0.5im, -3.8+0.0im, 1.5+0.0im, 0.0+0.0im]
        @test L.jinv(f) ≈ C[1.0, -2.0, 9.0, -60.5, 542.0, -6070.0, 81575.0]

        # jet_ops(6.0,0.5,0.4,2;N=6) coefficient jets — KLC literals
        ops = L.jet_ops(6.0, 0.5, 0.4, 2; N=6)
        @test ops.Δjet[1:3] ≈ C[24.25, 10.0, 2.0]
        @test ops.Kjet[1:3] ≈ C[13.5, 4.800000000000001, 0.8]
        @test isapprox(ops.invΔjet[1], 0.041237113402061855; rtol=1e-13)
        @test isapprox(ops.cjet[1], 0.0 + 0.5567010309278351im; rtol=1e-13)
        @test isapprox(ops.cjet[2], -0.0 - 0.03162929110426185im; rtol=1e-12)
    end

    @testset "BigFloat inheritance of jet_ops" begin
        setprecision(BigFloat, 256) do
            ops = L.jet_ops(big"6.0", big"0.99", big"0.1", 2; N=6)
            @test ops.CT === Complex{BigFloat}
            @test ops.Δjet[1] isa Complex{BigFloat}
            @test abs(ops.Δjet[1] - (big"6.0"^2 - 12 + big"0.99"^2)) < big"1e-60"
            @test maximum(abs.(L.jmul(ops.Δjet, ops.invΔjet) .- vcat([one(Complex{BigFloat})], zeros(Complex{BigFloat},6)))) < big"1e-60"
        end
    end
end

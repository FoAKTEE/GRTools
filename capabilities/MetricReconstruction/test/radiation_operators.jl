# Grounding for Radiation Operators — the separable jet algebra + GHP/Teukolsky ladder operators.
#
# ANALYTIC closed-form, NO live dep (mirrors legacy MetricReconstruction verify_operators.jl G-M1 +
# cov_operators_jets.jl). Every expected value is a derivative identity, a Leibniz/inverse identity,
# or a central-FD cross-check — never == whatever the code emits.

using Test
using MetricReconstruction.Radiation
import MetricReconstruction.Radiation as RAD   # for the underscore-prefixed coefficient jets

@testset "Radiation.Operators (jet algebra + ladders)" begin
    C = ComplexF64

    @testset "Taylor jets: _jet_mul (Leibniz), _jet_inv (reciprocal), _jet_deriv, _jet_add" begin
        r = 2.3
        u = C[r^3, 3r^2, 6r, 6, 0]      # r³ jet
        v = C[r^2, 2r, 2, 0, 0]         # r² jet
        # u·v = r⁵ with exact derivs
        w = RAD._jet_mul(u, v)
        @test maximum(abs.(w .- C[r^5, 5r^4, 20r^3, 60r^2, 120r])) < 1e-10
        # 1/v of r²: [r⁻²,−2r⁻³,6r⁻⁴,−24r⁻⁵,120r⁻⁶]
        iv = RAD._jet_inv(v)
        @test maximum(abs.(iv .- C[r^-2, -2r^-3, 6r^-4, -24r^-5, 120r^-6])) < 1e-12
        @test maximum(abs.(RAD._jet_mul(v, iv) .- C[1,0,0,0,0])) < 1e-12   # v·(1/v)=1
        @test RAD._jet_deriv(u) == C[3r^2, 6r, 6, 0, 0]
        @test RAD._jet_add(u, v) == u .+ v
    end

    @testset "analytic coefficient jets _Q_jet, _cot_jet vs finite difference" begin
        a, θ, ω, m = 0.5, 0.8, 0.4, 2
        Qj = RAD._Q_jet(C, a, θ, ω, m)
        Q(x) = -a*ω*sin(x) + m/sin(x)
        @test isapprox(Qj[1], Q(θ); rtol=1e-13)
        h = 1e-4
        @test isapprox(Qj[2], (Q(θ+h)-Q(θ-h))/(2h); rtol=1e-6)
        @test isapprox(Qj[3], (Q(θ+h)-2Q(θ)+Q(θ-h))/h^2; rtol=1e-4)
        cj = RAD._cot_jet(C, θ)
        cotf(x) = cos(x)/sin(x)
        @test isapprox(cj[1], cotf(θ); rtol=1e-13)
        @test isapprox(cj[2], (cotf(θ+h)-cotf(θ-h))/(2h); rtol=1e-6)
        @test isapprox(cj[3], (cotf(θ+h)-2cotf(θ)+cotf(θ-h))/h^2; rtol=1e-4)
        Δj = RAD._Delta_jet(C, a, 6.0)
        @test Δj[1] == 6.0^2 - 12 + a^2 && Δj[2] == 10.0 && Δj[3] == 2
    end

    @testset "radial/angular coefficient jets _radcoef_jet, _angcoef_jet" begin
        a, r, θ, ω, m = 0.5, 6.0, 0.8, 0.4, 2
        n, sgn = 0, -1
        cj = RAD._radcoef_jet(C, sgn, n, a, r, ω, m)
        K = (r^2+a^2)*ω - a*m; Δ = r^2-2r+a^2
        @test isapprox(cj[1], sgn*im*K/Δ + 2n*(r-1)/Δ; rtol=1e-12)
        cjb = RAD._radcoef_jet(C, sgn, 0, a, r, ω, m; bar=true)   # i→−i (real ω ⇒ K unchanged)
        @test isapprox(cjb[1], -sgn*im*K/Δ; rtol=1e-12)
        Qj = RAD._Q_jet(C, a, θ, ω, m); cotj = RAD._cot_jet(C, θ)
        acj = RAD._angcoef_jet(C, +1, 1, a, θ, ω, m)
        @test isapprox(acj[1], Qj[1] + cotj[1]; rtol=1e-12)        # sgn·Q + n cotθ
    end

    @testset "ladders Dn, Ln vs analytic action + finite difference (G-M1)" begin
        a, ω, m, r, θ = 0.1, 0.3, 2, 5.0, 1.1
        rj = C[r^3, 3r^2, 6r, 6, 0]
        aj = C[cos(θ), -sin(θ), -cos(θ), sin(θ), cos(θ)]
        F = sepfield(rj, aj)              # f = r³ ⊗ cosθ
        Δ = r^2 - 2r + a^2; K = (r^2+a^2)*ω - a*m
        # Dn[-1,0]: g = f' − iK/Δ f
        g_an = (3r^2 - im*K/Δ*r^3)*cos(θ)
        g_co = value(Dn(-1, 0, F; a=a, r=r, ω=ω, m=m))
        @test abs(g_co - g_an)/abs(g_an) < 1e-12
        # Ln[+1,1]: g = f_θ + (Q + cotθ) f
        Q = -a*ω*sin(θ) + m/sin(θ); cotθ = cos(θ)/sin(θ)
        gL_an = r^3*(-sin(θ) + (Q + cotθ)*cos(θ))
        gL_co = value(Ln(+1, 1, F; a=a, θ=θ, ω=ω, m=m))
        @test abs(gL_co - gL_an)/abs(gL_an) < 1e-12
        # composed Dn∘Dn vs central-FD of Dn[F] in r
        DF = Dn(-1, 0, F; a=a, r=r, ω=ω, m=m)
        DDF = Dn(-1, 0, DF; a=a, r=r, ω=ω, m=m)
        fval(rr) = value(Dn(-1, 0, sepfield(C[rr^3,3rr^2,6rr,6,0], aj); a=a, r=rr, ω=ω, m=m))
        hfd = 1e-5
        ddf_fd = (fval(r+hfd)-fval(r-hfd))/(2hfd) - im*K/Δ*value(DF)
        @test abs(value(DDF) - ddf_fd)/abs(ddf_fd) < 1e-9
        # a=0 hand calc
        a0, r0 = 0.0, 5.0
        F0 = sepfield(C[r0^4, 4r0^3, 12r0^2, 24r0, 24], C[1,0,0,0,0])
        Δ0 = r0^2-2r0; K0 = r0^2*ω
        @test abs(value(Dn(-1,0,F0; a=a0, r=r0, ω=ω, m=m)) - (4r0^3 - im*K0/Δ0*r0^4)) < 1e-10
    end

    @testset "RTJet algebra: product/quotient/inverse identities vs central-FD" begin
        a, r, θ = 0.3, 4.0, 0.9
        Z  = RAD._rt_zeta(C, a, r, θ)        # ζ = r − i a cosθ
        Zb = RAD._rt_zetabar(C, a, r, θ)
        # rtinv·Z == 1 (2-jet identity)
        one_jet = RAD.rtinv(Z) * Z
        @test isapprox(one_jet.v, 1; atol=1e-13)
        @test all(abs.((one_jet.r, one_jet.θ, one_jet.rr, one_jet.θθ, one_jet.rθ)) .< 1e-12)
        # Σ = ζ ζ̄ (real); compare against _rt_Sigma
        S1 = Z * Zb; S2 = RAD._rt_Sigma(C, a, r, θ)
        @test isapprox(S1.v, S2.v; rtol=1e-13) && isapprox(S1.r, S2.r; rtol=1e-12) &&
              isapprox(S1.θ, S2.θ; rtol=1e-12) && isapprox(S1.rr, S2.rr; rtol=1e-12)
        # RTJet ∂_r of 1/Δ vs central-FD
        invΔ = RAD.rtinv(RAD._rt_Delta(C, a, r))
        Δof(x) = 1/(x^2 - 2x + a^2); h = 1e-5
        @test isapprox(invΔ.r, (Δof(r+h)-Δof(r-h))/(2h); rtol=1e-7)
        # sepfield_rtjet matches direct value/derivs of f=r³⊗cosθ
        F = sepfield(C[r^3,3r^2,6r,6,0], C[cos(θ),-sin(θ),-cos(θ),sin(θ),cos(θ)])
        J = sepfield_rtjet(F)
        @test isapprox(J.v, r^3*cos(θ); rtol=1e-13)
        @test isapprox(J.r, 3r^2*cos(θ); rtol=1e-13)
        @test isapprox(J.θ, -r^3*sin(θ); rtol=1e-13)
        @test isapprox(J.rθ, 3r^2*(-sin(θ)); rtol=1e-13)
    end

    @testset "BigFloat inheritance (Complex{BigFloat})" begin
        setprecision(BigFloat, 256) do
            CB = Complex{BigFloat}
            r = big"5.0"; a = big"0.1"; ω = big"0.3"; m = 2
            rj = CB[r^3, 3r^2, 6r, 6, 0]
            aj = CB[1, 0, 0, 0, 0]
            g = value(Dn(-1, 0, sepfield(rj, aj); a=a, r=r, ω=ω, m=m))
            Δ = r^2 - 2r + a^2; K = (r^2+a^2)*ω - a*m
            @test g isa CB
            @test abs(g - (3r^2 - im*K/Δ*r^3)) < big"1e-60"
        end
    end
end

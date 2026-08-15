# Grounding for Radiation Modes — teuk_lambda / radial_jet / angular_jet / clear_mode_caches!.
#
# Characterized against the already-green GrTools Harmonics/Perturbation surfaces (themselves
# grounded vs SWSH/Teukolsky) plus closed-form/field-equation identities (mirrors legacy
# MetricReconstruction verify_modes G-M2/G-M3):
#   • teuk_lambda  vs Harmonics.SpinWeightedSpheroidalEigenvalue + a=0,s=−2 closed form (l−1)(l+2).
#   • radial_jet[1]/[2] == direct TeukolskyRadial value / teukolsky_radial_deriv; the eqR residual
#     Δ R'' + (s+1)Δ' R' + V R ≈ 0 (R''…R'''' built by collapsing eqR); R'' vs central-FD of R'.
#   • angular_jet[1] == direct SpinWeightedSpheroidalHarmonicS value; the eqS residual ≈ 0.

using Test
using MetricReconstruction.Radiation
using Harmonics: SpinWeightedSpheroidalEigenvalue, SpinWeightedSpheroidalHarmonicS, swsh_derivative
using Perturbation: TeukolskyRadial, teukolsky_radial_deriv
import Perturbation   # Move C: radial_jet's radial/λ caches now live in Perturbation

@testset "Radiation.Modes (teuk_lambda / radial_jet / angular_jet)" begin
    @testset "teuk_lambda vs Harmonics eigenvalue + a=0 closed form" begin
        a, ω = 0.1, 0.3
        for (s, l, m) in ((-2,2,2), (2,2,1), (-2,3,2))
            @test teuk_lambda(s, l, m, a, ω) ≈ SpinWeightedSpheroidalEigenvalue(s, l, m, a*ω) rtol=1e-13
        end
        # a=0, s=−2 ⇒ γ=0, spin-weighted spherical eigenvalue λ = l(l+1) − s(s+1) = (l−1)(l+2)
        for l in (2, 3, 4)
            @test teuk_lambda(-2, l, 0, 0.0, 0.5) ≈ (l-1)*(l+2) rtol=1e-12
        end
    end

    @testset "radial_jet: direct value/deriv + eqR residual + FD (G-M3)" begin
        a, ω, l, m, r = 0.1, 0.3, 2, 2, 6.0
        for s in (-2, 2)
            λ = teuk_lambda(s, l, m, a, ω)
            sols = TeukolskyRadial(s, l, m, a, ω; method=:MST, lambda=λ)
            rj = radial_jet(s, l, m, a, ω, r; bc=:Up)
            # [1],[2] match the GrTools surface directly
            @test rj[1] ≈ sols.Up(r) rtol=1e-12
            @test rj[2] ≈ teukolsky_radial_deriv(sols.Up, r) rtol=1e-12
            # eqR residual (R'' built by collapsing this equation ⇒ ≈ 0 by construction)
            R, Rp, Rpp = rj[1], rj[2], rj[3]
            Δ = r^2-2r+a^2; dΔ = 2r-2; K = (r^2+a^2)*ω - a*m
            V = (K^2 - 2im*s*(r-1)*K)/Δ + 4im*s*ω*r - λ
            @test abs(Δ*Rpp + (s+1)*dΔ*Rp + V*R)/abs(R) < 1e-9
            # independent: R'' from ODE vs central-FD of R' (teukolsky_radial_deriv)
            hf = 1e-4
            Rpp_fd = (teukolsky_radial_deriv(sols.Up, r+hf) - teukolsky_radial_deriv(sols.Up, r-hf))/(2hf)
            @test abs(Rpp - Rpp_fd)/abs(Rpp_fd) < 1e-6
            # R''' vs central-FD of R'' (eqR-derivative collapse)
            rjp = radial_jet(s, l, m, a, ω, r+hf; bc=:Up); rjm = radial_jet(s, l, m, a, ω, r-hf; bc=:Up)
            @test abs(rj[4] - (rjp[3]-rjm[3])/(2hf))/abs(rj[4]) < 1e-5
        end
    end

    @testset "angular_jet: direct value + eqS residual (G-M2)" begin
        a, ω, l, m, θ = 0.1, 0.3, 2, 2, 1.0
        γ = a*ω
        for s in (-2, 2)
            S = SpinWeightedSpheroidalHarmonicS(s, l, m, γ; method=:SphericalExpansion)
            aj = angular_jet(s, l, m, a, ω, θ)
            @test aj[1] ≈ swsh_derivative(S, 0, 0, θ, 0.0) rtol=1e-12
            @test aj[2] ≈ swsh_derivative(S, 1, 0, θ, 0.0) rtol=1e-12
            Sv, Sθ, Sθθ = aj[1], aj[2], aj[3]
            λ = teuk_lambda(s, l, m, a, ω); A = λ + 2*m*γ - γ^2
            sθ = sin(θ); cθ = cos(θ)
            res = (cθ/sθ*Sθ + Sθθ) + (γ^2*cθ^2 - 2*γ*s*cθ - (m+s*cθ)^2/sθ^2 + s + A)*Sv
            @test abs(res) < 1e-9
        end
    end

    @testset "clear_mode_caches! resets the memo tables (local λ + relocated radial_jet caches)" begin
        radial_jet(-2, 2, 2, 0.1, 0.3, 6.0)         # populate the (relocated) Perturbation radial cache
        teuk_lambda(-2, 2, 2, 0.1, 0.3)             # populate the local MR λ cache
        @test clear_mode_caches!() === nothing
        @test isempty(MetricReconstruction.Radiation._LAMBDA_CACHE)
        @test isempty(Perturbation._radial_jet_cache)        # Move C: radial-solution cache moved to Perturbation
        @test isempty(Perturbation._radial_jet_lambda_cache) # Move C: radial_jet's internal λ cache moved too
    end
end

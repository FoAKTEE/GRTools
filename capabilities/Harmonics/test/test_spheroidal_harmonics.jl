# Module 4 — SpheroidalHarmonics: the callable S_{lm}(θ,φ) Float64-VALUE carrier in
# SphericalExpansion (+Leaver) reps with mixed ∂_θ/∂_φ/∂_γ. Depends up the in-batch DAG on
# module 1 (Y seed), module 2 (sign-fixed eigenvector coeffs), module 3 (small-γ d-tables).
# The Leaver coeff norm consumes NumericsCore.kummer_1F1 (no local Kummer copy).
#
# Grounding:
#  (A) CHARACTERIZATION vs the LIGHT legacy oracle at FULL precision — S(θ,φ),
#      expansion_coefficients, swsh_derivative, swsh_gamma_derivative, SWSHSLeaverFunction(θ,φ),
#      SpinWeightedSpheroidalHarmonicSSeries over (s,l,m,γ,θ,φ), both method=:SphericalExpansion
#      and :Leaver.
#  (B) ANALYTIC — γ→0 limit S → {}_sY_{lm} (the module-1 seed); unit-norm spectral coeffs;
#      ∂_φ S = im·m·S consistency.

using Test
using Harmonics
# Legacy oracle SpinWeightedSpheroidalHarmonics {SpinWeightedSpheroidalHarmonicS
# (.eigenvalue/.coeffs/callable, :SphericalExpansion + :Leaver), expansion_coefficients,
# swsh_derivative, swsh_gamma_derivative, SpinWeightedSpheroidalHarmonicSSeries} BAKED into a
# fixture (captured once, same (s,l,m,γ,d1,d2,θ,φ,order) grids as the loops below). Same `==`.
include(joinpath(@__DIR__, "fixtures", "spheroidal_harmonics_ref.jl"))
# ^ FIXSH_Seig, FIXSH_Scoeffs, FIXSH_Sval, FIXSH_expc, FIXSH_dS, FIXSH_gdS,
#   FIXSH_Slcoeffs, FIXSH_Slval, FIXSH_SldS, FIXSH_Sseries

@testset "SpheroidalHarmonics" begin

    _cases = [(-2, 2, 2), (-2, 3, -1), (0, 3, 0), (0, 4, 2), (2, 4, -3), (-2, 5, 5)]
    _angles = [(0.4, 0.0), (1.1, 0.9), (2.6, 1.7)]

    # ---- (A) characterization: SphericalExpansion S(θ,φ) + coeffs ----------------------------
    @testset "SphericalExpansion S(θ,φ) == legacy (full precision)" begin
        for (s, l, m) in _cases, γ in (0.3, 1.2)
            Sgot = SpinWeightedSpheroidalHarmonicS(s, l, m, γ)
            @test Sgot.eigenvalue == FIXSH_Seig[(s, l, m, γ)]
            @test Sgot.coeffs == FIXSH_Scoeffs[(s, l, m, γ)]
            for (θ, φ) in _angles
                @test Sgot(θ, φ) == FIXSH_Sval[(s, l, m, γ, θ, φ)]
            end
        end
    end

    @testset "expansion_coefficients == legacy" begin
        for (s, l, m) in [(-2, 3, 2), (0, 4, -1)], γ in (0.5, 1.5)
            Sgot = SpinWeightedSpheroidalHarmonicS(s, l, m, γ)
            @test expansion_coefficients(Sgot) == FIXSH_expc[(s, l, m, γ)]
        end
    end

    @testset "swsh_derivative (SphericalExpansion) == legacy" begin
        for (s, l, m) in _cases, γ in (0.4, 1.0), d1 in (0, 1, 2), d2 in (0, 1)
            Sgot = SpinWeightedSpheroidalHarmonicS(s, l, m, γ)
            for (θ, φ) in _angles
                @test swsh_derivative(Sgot, d1, d2, θ, φ) == FIXSH_dS[(s, l, m, γ, d1, d2, θ, φ)]
            end
        end
    end

    @testset "swsh_gamma_derivative == legacy" begin
        for (s, l, m) in _cases, γ in (0.4, 1.1), d1 in (0, 1), d2 in (0, 1)
            for (θ, φ) in _angles
                @test swsh_gamma_derivative(s, l, m, γ, d1, d2, θ, φ) ==
                      FIXSH_gdS[(s, l, m, γ, d1, d2, θ, φ)]
            end
        end
    end

    @testset "swsh_gamma_derivative is regular at γ=0" begin
        s, l, m, θ, φ = -2, 3, 2, 0.7, 0.2
        got = swsh_gamma_derivative(s, l, m, 0.0, 0, 0, θ, φ)
        h = 1e-5
        Sp = SpinWeightedSpheroidalHarmonicS(s, l, m, h)
        Sm = SpinWeightedSpheroidalHarmonicS(s, l, m, -h)
        fd = (Sp(θ, φ) - Sm(θ, φ)) / (2h)
        @test isfinite(got)
        @test got ≈ fd rtol = 1e-7 atol = 1e-9
    end

    # ---- (A) characterization: Leaver rep (Kummer-norm rewrite over NumericsCore) -------------
    @testset "Leaver SWSHSLeaverFunction(θ,φ) + coeffs == legacy" begin
        for (s, l, m) in _cases, γ in (0.3, 1.2)
            Sgot = SpinWeightedSpheroidalHarmonicS(s, l, m, γ; method=:Leaver)
            @test Sgot.coeffs == FIXSH_Slcoeffs[(s, l, m, γ)]
            for (θ, φ) in _angles
                @test Sgot(θ, φ) == FIXSH_Slval[(s, l, m, γ, θ, φ)]
                @test swsh_derivative(Sgot, 1, 0, θ, φ) == FIXSH_SldS[(s, l, m, γ, 1, 0, θ, φ)]
                @test swsh_derivative(Sgot, 2, 1, θ, φ) == FIXSH_SldS[(s, l, m, γ, 2, 1, θ, φ)]
            end
        end
    end

    @testset "SpinWeightedSpheroidalHarmonicSSeries == legacy (full precision)" begin
        for (s, l, m) in _cases, order in (4, 6)
            for (θ, φ) in _angles
                ref = FIXSH_Sseries[(s, l, m, order, θ, φ)]
                got = SpinWeightedSpheroidalHarmonicSSeries(s, l, m, θ, φ; order=order)
                @test got == ref
            end
        end
    end

    # ---- (B) analytic identities --------------------------------------------------------------
    @testset "γ→0 limit S → {}_sY_{lm} (module-1 seed)" begin
        for (s, l, m) in _cases
            S0 = SpinWeightedSpheroidalHarmonicS(s, l, m, 0.0)
            for (θ, φ) in _angles
                @test S0(θ, φ) ≈ SpinWeightedSphericalHarmonicY(s, l, m, θ, φ) rtol = 1e-12
            end
        end
    end

    @testset "unit-norm spectral coeffs" begin
        for (s, l, m) in _cases, γ in (0.3, 1.4)
            S = SpinWeightedSpheroidalHarmonicS(s, l, m, γ)
            @test sum(abs2, S.coeffs) ≈ 1.0 rtol = 1e-12
        end
    end

    @testset "∂_φ S = im·m·S" begin
        for (s, l, m) in _cases, γ in (0.4, 1.0)
            S = SpinWeightedSpheroidalHarmonicS(s, l, m, γ)
            for (θ, φ) in _angles
                @test swsh_derivative(S, 0, 1, θ, φ) ≈ im * m * S(θ, φ) rtol = 1e-12
            end
        end
    end
end

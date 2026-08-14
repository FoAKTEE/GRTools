# Characterization of HybridC module 1 (wasabic_deps): the WaSABI-C interpolant bundle +
# dependent-variable closed forms + the BigSpecial-routed complex PolyGamma.
#
# GROUNDING (WaSABI NOT a live dep — baked WL-oracle literals from the legacy oracle16a/16b runs,
# captured verbatim from packages/WaSABI/test/inspiral_WaSABI-C_deps_test.jl / _test.jl which graded
# them vs `wolframscript -file oracle16a.wl`/`oracle16b.wl` over the upstream paclet):
#   * 2D Kerr-grid interpolants (FI1/FH1/FIS2/FHS2/DU) at WL `Interpolation[...][x]` points (rtol 1e-8).
#   * (2,2) amplitude Re22/Im22 + mixed partials vs WL Derivative[…][itp] (rtol 1e-8).
#   * 1D Schwarz fluxes (order-8 / log-exp) vs WL (rtol 1e-8); 3D ωcrit vs WL (rtol 1e-8).
#   * x=ω^(2/3), s, σ analytic + WL paclet probe literals (rtol 1e-14).
#   * intrinsic node recovery (a global Hermite-All interpolant reproduces its grid nodes).
#   * complex PolyGamma[0,3+2iχ/√(1-χ²)] vs WL Im[PolyGamma] literals (rtol 1e-12) — routed to
#     BigSpecial.PolyGamma (the AHA-collapsed canonical complex digamma).

using Test
using Hybrid.HybridC
using Hybrid.Interpolation: partial

# Baked fixture copy of the upstream WaSABI data (provenance: test/fixtures/wasabi/README.md).
const _HC_BASE = joinpath(@__DIR__, "fixtures", "wasabi")
const _HC_INSP = joinpath(_HC_BASE, "InspiralModels", "sf_data")
const _HC_AMP = joinpath(_HC_BASE, "AmplitudeModels", "sf_amp_data")

@testset "HybridC.wasabic_deps" begin
    mdl = wasabic_model(_HC_INSP, _HC_AMP)

    @testset "2D Kerr-grid interpolants vs WL Interpolation[...][x]" begin
        @test isapprox(mdl.FI1SFKerr(-0.7, 0.8), 5.341818848111235; rtol = 1e-8)
        @test isapprox(mdl.FI1SFKerr(-2.0, 0.5), 5.101032214159012; rtol = 1e-8)
        @test isapprox(mdl.FH1SFKerr(-0.7, 0.8), -0.9110084787826039; rtol = 1e-8)
        @test isapprox(mdl.FH1SFKerr(-2.0, 0.5), -3.5470261476190834; rtol = 1e-8)
        @test isapprox(mdl.FIS2Kerr(-0.7, 0.8), -0.40037167727102424; rtol = 1e-8)
        @test isapprox(mdl.FIS2Kerr(-2.0, 0.5), -0.10035992004953022; rtol = 1e-8)
        @test isapprox(mdl.FHS2Kerr(-0.7, 0.8), -0.47281338021615965; rtol = 1e-8)
        @test isapprox(mdl.FHS2Kerr(-2.0, 0.5), -5.996893863080762; rtol = 1e-8)
        @test isapprox(mdl.DU(-0.7, 0.8), -0.20447296014665062; rtol = 1e-8)
        @test isapprox(mdl.DU(-2.0, 0.5), -0.11257564662146778; rtol = 1e-8)
    end

    @testset "(2,2) amplitude interpolants + partials vs WL" begin
        @test isapprox(mdl.Re22(0.6, 0.2), -2.751168837640576; rtol = 1e-8)
        @test isapprox(mdl.Re22(0.5, -0.3), -2.9738731783672163; rtol = 1e-8)
        @test isapprox(mdl.Im22(0.6, 0.2), 0.6550167777151532; rtol = 1e-8)
        @test isapprox(mdl.Im22(0.5, -0.3), 0.2806880593495512; rtol = 1e-8)
        @test isapprox(partial(mdl.Re22, 0.6, 0.2, 0, 1), 0.5515082601534531; rtol = 1e-8)
        @test isapprox(partial(mdl.Re22, 0.5, -0.3, 0, 1), 0.49864166587755904; rtol = 1e-8)
        @test isapprox(partial(mdl.Im22, 0.6, 0.2, 1, 0), -0.6869005555609442; rtol = 1e-8)
        @test isapprox(partial(mdl.Im22, 0.5, -0.3, 1, 0), -0.42275227846807023; rtol = 1e-8)
    end

    @testset "1D Schwarz fluxes (order-8 / log-exp) vs WL" begin
        @test isapprox(mdl.FI1SFSchw(7.0), 0.00039963398939231535; rtol = 1e-8)
        @test isapprox(mdl.FI1SFSchw(12.5), 0.000019756512511570108; rtol = 1e-8)
        @test isapprox(mdl.FI1SFSchw(25.0), 6.157899216259335e-7; rtol = 1e-8)
        @test isapprox(mdl.FI2SFSchw(2.0), -0.00018601426872875133; rtol = 1e-8)
        @test isapprox(mdl.FI2SFSchw(2.5), -6.035202961095554e-6; rtol = 1e-8)
        @test isapprox(mdl.FI2SFSchw(3.0), -2.715512719839147e-7; rtol = 1e-8)
    end

    @testset "3D ωcrit (StopHybridFixedChi) vs WL" begin
        @test isapprox(mdl.ωcrit(0.1, 0.5, 0.3), 0.06536131922368031; rtol = 1e-8)
        @test isapprox(mdl.ωcrit(0.2, 0.8, 0.7), 0.0639397431291707; rtol = 1e-8)
    end

    @testset "interpolant node recovery (intrinsic)" begin
        for itp in (mdl.FI1SFKerr, mdl.FH1SFKerr, mdl.DU, mdl.Re22, mdl.Im22)
            for (i, j) in ((1, 1), (5, 7), (length(itp.x), length(itp.y)))
                @test isapprox(itp(itp.x[i], itp.y[j]), itp.z[i, j]; rtol = 1e-9, atol = 1e-10)
            end
        end
    end

    @testset "analytic dependent vars x, s, σ (no oracle)" begin
        @test wc_x(0.02) == 0.02^(2 / 3)
        ν, χ1, χ2 = 0.1, 0.5, 0.3
        r = sqrt(1 - 4ν)
        @test wc_s(ν, χ1, χ2) ≈ ((1 + r)^2 * χ1 + (-1 + r)^2 * χ2) / 4
        @test wc_sigma(ν, χ1, χ2) ≈ (-(1 + r) * χ1 - (-1 + r) * χ2) / 2
        @test isapprox(wc_s(ν, χ1, χ2), 0.39745966692414836; rtol = 1e-14)
        @test isapprox(wc_sigma(ν, χ1, χ2), -0.4098386676965934; rtol = 1e-14)
    end

    @testset "complex PolyGamma[0,·] (BigSpecial-routed) vs WL Im[PolyGamma]" begin
        # WL: Im[PolyGamma[0, 3 + 2 I χ/Sqrt[1-χ²]]] (oracle16b.wl "DIGAMMA" rows).
        for (χ, wl) in ((0.5, 0.428626696619868), (0.7, 0.6610320482782235),
                        (-0.4, -0.3323527553853915))
            v = imag(PolyGamma(0, 3 + 2im * χ / sqrt(1 - χ^2)))
            @test isapprox(v, wl; rtol = 1e-12)
        end
    end
end

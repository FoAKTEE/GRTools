# Characterization of the 0PAKerr (2,2) strain-amplitude model (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — the literal 0PAKerr.m closed form + a synthetic 2D-grid fixture
# + the light live Spacetime.kerr_isco):
#   * real_amp22 == (1/(2(1-χ1 ω)^(2/3))) M ν ω^(2/3) RealAmp22[c1, c2] with the VERBATIM
#     0PAKerr.m abscissa remap c1 = 0.7988880369434854 + 0.2901432317650783 Log[1-χ1],
#     c2 = -1 + 2 √(ω^(2/3) rISCO[χ1]/(1-χ1 ω)^(2/3))  (constants copied from 0PAKerr.m:41).
#   * mode_keys_0PAKerr == [(2,2)].

using Test
using Hybrid.Models
using Spacetime: kerr_isco

# Literal magic constants, copied verbatim from ref-code/WaSABI/AmplitudeModels/0PAKerr.m:41.
const C1_OFFSET = 0.7988880369434854
const C1_SLOPE = 0.2901432317650783

@testset "Models.amplitude_0PAKerr" begin
    # Synthetic separable RealAmp22 surface over the (c1, c2) abscissa range.
    c1n = collect(range(0.55, 0.85; length = 5))
    c2n = collect(range(-0.9, 0.9; length = 6))
    R(c1, c2) = 3.0 + 0.8 * c1 - 0.4 * c2 + 0.2 * c1 * c2
    grid = [((c1, c2), R(c1, c2)) for c1 in c1n for c2 in c2n]
    model = Amp0PAKerr(interp2_hermite(grid))

    @testset "real_amp22 literal remap + prefactor" begin
        for (ω, χ1, ν, M) in ((0.04, 0.1, 1e-3, 1.0), (0.05, 0.3, 2e-3, 1.0), (0.03, 0.5, 1e-3, 2.0))
            c1 = C1_OFFSET + C1_SLOPE * log(1 - χ1)
            c2 = -1 + 2 * sqrt(ω^(2 / 3) * kerr_isco(χ1) / (1 - χ1 * ω)^(2 / 3))
            exp_amp = (1 / (2 * (1 - χ1 * ω)^(2 / 3))) * M * ν * ω^(2 / 3) * R(c1, c2)
            @test real_amp22(model, ω, χ1, ν, M) ≈ exp_amp rtol = 1e-9
        end
    end

    @testset "mode_keys" begin
        @test mode_keys_0PAKerr() == [(2, 2)]
    end
end

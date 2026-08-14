# Module 3 — SpheroidalEigenvalueSeries: small-γ Taylor (Stewart/Press-Teukolsky recurrence) +
# large-γ Laurent (Casals-Ottewill-Warburton arXiv:1810.00432) eigenvalue coefficient lists.
# Pure recurrences, independent of the spectral/Leaver dispatcher.
#
# Grounding:
#  (A) CHARACTERIZATION vs the LIGHT legacy oracle at FULL precision —
#      SpinWeightedSpheroidalEigenvalueSeries / ...LargeGamma coefficient VECTORS and the
#      private _swsh_alpha/_swsh_beta/_swsh_q over representative (s,l,m,order).
#  (B) ANALYTIC/self-consistency — c_0 = l(l+1)-s(s+1); the small-γ series Σ c_i γ^i reproduces
#      SpinWeightedSpheroidalEigenvalue(s,l,m,γ) at small γ to series tolerance; the large-γ
#      leading entry is 2(q-m); imaginary part negligible (<1e-9).

using Test
using Harmonics
# Legacy oracle SpinWeightedSpheroidalHarmonics {SpinWeightedSpheroidalEigenvalueSeries,
# ...LargeGamma, _swsh_q, _swsh_alpha, _swsh_beta} BAKED into a fixture (captured once, same
# args as the loops below). Same characterization at the SAME `==`.
include(joinpath(@__DIR__, "fixtures", "spheroidal_eigenvalue_series_ref.jl"))
# ^ FIXES_series, FIXES_largeG, FIXES_q, FIXES_alpha, FIXES_beta

@testset "SpheroidalEigenvalueSeries" begin

    _cases = [(-2, 2, 0), (-2, 3, 2), (-2, 4, -1), (0, 2, 1), (0, 5, 0), (2, 3, -2), (2, 4, 4)]

    # ---- (A) characterization: small-γ coefficient vectors ------------------------------------
    @testset "small-γ EigenvalueSeries == legacy (full precision)" begin
        for (s, l, m) in _cases, order in (4, 6, 8)
            ref = FIXES_series[(s, l, m, order)]
            got = SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=order)
            @test got == ref
        end
    end

    @testset "large-γ EigenvalueSeriesLargeGamma == legacy (full precision)" begin
        for (s, l, m) in _cases, order in (2, 4, 5)
            ref = FIXES_largeG[(s, l, m, order)]
            got = SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s, l, m; order=order)
            @test got == ref
        end
    end

    @testset "private α/β/q helpers == legacy (full precision)" begin
        for (s, l, m) in _cases
            @test Harmonics._swsh_q(s, l, m) == FIXES_q[(s, l, m)]
            for ll in (l-1):(l+2)
                @test Harmonics._swsh_alpha(s, ll, m) == FIXES_alpha[(s, ll, m)]
                @test Harmonics._swsh_beta(s, ll, m) == FIXES_beta[(s, ll, m)]
            end
        end
    end

    # ---- (B) analytic / self-consistency ------------------------------------------------------
    @testset "c_0 = l(l+1) - s(s+1)" begin
        for (s, l, m) in _cases
            @test SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=4)[1] ==
                  l * (l + 1) - s * (s + 1)
        end
    end

    @testset "small-γ Σ c_i γ^i reproduces the eigenvalue at small γ" begin
        for (s, l, m) in _cases
            c = SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=8)
            for γ in (0.02, 0.05)
                ser = sum(c[i+1] * γ^i for i in 0:length(c)-1)
                exact = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SphericalExpansion)
                @test ser ≈ exact rtol = 1e-7
            end
        end
    end

    @testset "large-γ leading entry = 2(q-m)" begin
        for (s, l, m) in _cases, order in (2, 4)
            cl = SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s, l, m; order=order)
            @test cl[end] == 2.0 * (Harmonics._swsh_q(s, l, m) - m)   # γ^1 coeff of λ
        end
    end

    @testset "imaginary-part guard (<1e-9): series coeffs are real" begin
        # the guard is internal; a complex coefficient would throw. Exercise it.
        for (s, l, m) in _cases
            @test SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=10) isa Vector{Float64}
        end
        @test_throws ArgumentError SpinWeightedSpheroidalEigenvalueSeries(-2, 3, 2; order=-1)
    end
end

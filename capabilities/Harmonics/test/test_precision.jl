# Module — Precision: the thin type-generic precision-policy primitive for the angular solve
# (design §4.2 'Precision asymmetry', §6 row 375). Surfaces the load-bearing asymmetry KLC lives
# with: the Leaver eigenvalue path is type-generic (γ::BigFloat inherits 256..2048-bit → FLOOR
# eps(T)) while the harmonic VALUE tier (SWSHSFunction Float64 coeffs, downcast Y) is Float64-
# capped (CEILING eps(Float64)).
#
# Grounding — NO oracle call needed, CLOSED-FORM contract + analytic self-convergence:
#  (A) angular_precision_floor() reports eps(Float64) for the value tier (closed form) and
#      eps(BigFloat) tracking the ACTIVE setprecision for the eigenvalue tier; the type-generic
#      angular_precision_floor(T) reports eps(T).
#  (B) precision-DOUBLING self-convergence demonstrating the floor is REAL: the BigFloat Leaver
#      eigenvalue at 256 vs 512 bits AGREES to ~the reported eigenvalue floor (path inherits
#      precision), while the Float64-capped value tier does NOT improve past eps(Float64) —
#      characterizing the named ceiling vs floor.

using Test
using Harmonics
# NOTE: this file is closed-form/self-convergence only — it never called the legacy oracle
# (the import was unused). No fixture needed. (Phase7-B: legacy import dropped.)

@testset "Precision (angular_precision_floor)" begin

    # ---- (A) closed-form policy report --------------------------------------------------------
    @testset "value tier == eps(Float64) ceiling (closed form)" begin
        p = angular_precision_floor()
        @test p.value_tier == eps(Float64)
        @test angular_precision_floor(Float64).value_tier == eps(Float64)
        # value tier is the Float64 ceiling regardless of active BigFloat precision
        setprecision(BigFloat, 512) do
            @test angular_precision_floor().value_tier == eps(Float64)
        end
    end

    @testset "eigenvalue tier == eps(T) floor, tracking active setprecision" begin
        # type-generic: explicit T threads eps(T) (no hard-Float64 on the generic path)
        @test angular_precision_floor(Float64).eigenvalue_tier == eps(Float64)
        setprecision(BigFloat, 256) do
            @test angular_precision_floor().eigenvalue_tier == eps(BigFloat)   # tracks 256-bit
            @test angular_precision_floor(BigFloat).eigenvalue_tier == eps(BigFloat)
        end
        setprecision(BigFloat, 512) do
            @test angular_precision_floor().eigenvalue_tier == eps(BigFloat)   # tracks 512-bit
        end
        # the floor is strictly below the value-tier ceiling once BigFloat is engaged
        setprecision(BigFloat, 256) do
            p = angular_precision_floor()
            @test p.eigenvalue_tier < p.value_tier
            @test p.eigenvalue_tier < 1e-70
        end
    end

    # ---- (B) precision-doubling self-convergence: floor real, ceiling capped -------------------
    @testset "eigenvalue path inherits floor; Float64 value tier capped at ceiling" begin
        s, l, m = -2, 4, 2
        λ512 = setprecision(BigFloat, 512) do
            SpinWeightedSpheroidalEigenvalue(s, l, m, big"0.5")     # :Auto ⇒ Leaver for BigFloat
        end
        setprecision(BigFloat, 256) do
            floor256 = angular_precision_floor().eigenvalue_tier   # == eps(BigFloat) @256
            ceil_f64 = angular_precision_floor().value_tier        # == eps(Float64)
            λ256 = SpinWeightedSpheroidalEigenvalue(s, l, m, big"0.5")
            relerr_big = abs(λ256 - λ512) / abs(λ512)
            # eigenvalue tier INHERITS precision: agrees to ~the reported 256-bit floor
            @test relerr_big < floor256 * 16
            @test relerr_big > floor256 * 1e-6                     # genuinely at the floor, not exact
            # Float64-capped value tier: spectral eigenvalue stuck at ~the eps(Float64) ceiling
            λf64 = SpinWeightedSpheroidalEigenvalue(s, l, m, 0.5; method=:SphericalExpansion)
            relerr_f64 = abs(λf64 - λ512) / abs(λ512)
            @test relerr_f64 < ceil_f64 * 1e3                      # near the ceiling, not below it
            @test relerr_f64 > ceil_f64 * 1e-3                     # does NOT reach the BigFloat floor
            @test relerr_f64 > floor256 * 1e40                     # ceiling ≫ eigenvalue floor
        end
    end
end

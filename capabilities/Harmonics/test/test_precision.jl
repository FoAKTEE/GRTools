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

    # ---- (C) Leaver COEFFICIENT recurrence inherits the caller precision (repro-audit G10) -----
    # `_swsh_leaver_coeffs` used to run its BigFloat recurrence at a HARD-CODED 256 bits, so a
    # BigFloat-512 caller had its intermediate silently truncated — contradicting the "precision
    # inherited end-to-end" contract. The working precision is now `max(256, precision(BigFloat))`.
    @testset "Leaver coeff recurrence: max(256, precision(BigFloat))" begin
        # (C1) the policy itself: 256 is a FLOOR, the ambient precision wins above it
        for bits in (64, 128, 256)
            setprecision(BigFloat, bits) do
                @test Harmonics._leaver_coeff_precision() == 256
            end
        end
        for bits in (512, 1024, 2048)
            setprecision(BigFloat, bits) do
                @test Harmonics._leaver_coeff_precision() == bits
            end
        end
        @test Harmonics._leaver_coeff_precision() == max(256, precision(BigFloat))

        # (C2) no hard-coded 256-bit intermediate survives in the source
        src = read(joinpath(pkgdir(Harmonics), "src", "SpheroidalHarmonics.jl"), String)
        @test !occursin("setprecision(BigFloat, 256)", src)

        # (C3) SELF-CONVERGENCE of the Float64-capped output: raising the recurrence from 256 to
        # 512/1024 bits does not move the returned coefficients or the harmonic value — the
        # 256-bit floor was already converged for these modes, so the fix is result-neutral here
        # while removing the truncation for callers that need more.
        for (s, l, m, γ) in ((-2, 2, 2, 0.3), (-2, 4, 2, 0.5), (0, 3, 1, 0.7), (2, 5, -3, 1.1))
            λ256, c256, k1, k2 = setprecision(BigFloat, 256) do
                Harmonics._swsh_leaver_coeffs(s, l, m, γ)
            end
            for bits in (512, 1024)
                λb, cb, k1b, k2b = setprecision(BigFloat, bits) do
                    Harmonics._swsh_leaver_coeffs(s, l, m, γ)
                end
                @test (k1b, k2b) == (k1, k2)
                @test λb == λ256                                   # Float64 γ ⇒ Float64 eigenvalue
                @test length(cb) == length(c256)
                @test maximum(abs, cb .- c256) <= 1e-15 * maximum(abs, c256)
            end
            # the callable Leaver harmonic agrees with the DEFAULT (Float64-ambient) evaluation
            S_def = SpinWeightedSpheroidalHarmonicS(s, l, m, γ; method = :Leaver)
            S_512 = setprecision(BigFloat, 512) do
                SpinWeightedSpheroidalHarmonicS(s, l, m, γ; method = :Leaver)
            end
            for θ in (0.4, 1.1, 2.3)
                v = S_def(θ, 0.6); w = S_512(θ, 0.6)
                @test abs(w - v) <= 1e-15 * max(abs(v), 1e-300)
            end
        end

        # (C4) the BigFloat EIGENVALUE tier stays at its own floor (unchanged by C1-C3):
        # 512 vs 1024 bits self-converge far below 1e-70.
        λ512 = setprecision(BigFloat, 512) do
            SpinWeightedSpheroidalEigenvalue(-2, 4, 2, big"0.5")
        end
        λ1024 = setprecision(BigFloat, 1024) do
            SpinWeightedSpheroidalEigenvalue(-2, 4, 2, big"0.5")
        end
        @test abs(λ512 - λ1024) / abs(λ1024) < 1e-70
        @test abs(Float64(λ512) - Float64(SpinWeightedSpheroidalEigenvalue(-2, 4, 2, 0.5;
                                                method = :Leaver))) < 1e-15 * abs(Float64(λ512))
    end
end

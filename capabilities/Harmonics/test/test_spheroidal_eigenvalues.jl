# Module 2 — SpheroidalEigenvalues: λ(s,l,m,γ) + ∂_γλ via ONE dispatcher over the
# pentadiagonal-LAPACK spectral (Float64) + type-generic Leaver continued-fraction (BigFloat)
# backends. The Leaver path consumes NumericsCore.cf_finite/cf_convergent (no local CF copy).
#
# Grounding:
#  (A) CHARACTERIZATION vs the LIGHT legacy oracle (LinearAlgebra-only) at FULL precision —
#      SpinWeightedSpheroidalEigenvalue (both :SphericalExpansion and :Leaver), its Deriv, and
#      the private _swsh_eigenvalue_leaver / _swsh_spectral_coeffs over (s∈{-2,0,2}, several l,m,
#      γ small/large). The new CF-primitive rewrite must reproduce the legacy verbatim CF bit-for-bit.
#  (B) ANALYTIC — γ=0 limit λ=l(l+1)-s(s+1) and dλ/dγ|_0 = -2m s²/(l+l²)-2m; spectral↔Leaver
#      cross-agreement; BigFloat γ at setprecision(256) beats the Float64 spectral path
#      (precision inheritance).

using Test
using Harmonics
# Legacy oracle SpinWeightedSpheroidalHarmonics {SpinWeightedSpheroidalEigenvalue(:SphericalExpansion
# /:Leaver), ...Deriv, _swsh_eigenvalue_leaver, _swsh_spectral_coeffs} BAKED into a fixture
# (captured once, same args as the loops below). Same characterization at the SAME `==`.
include(joinpath(@__DIR__, "fixtures", "spheroidal_eigenvalues_ref.jl"))
# ^ FIXEV_spec, FIXEV_deriv, FIXEV_leaver, FIXEV_leaverfront, FIXEV_coeffs

@testset "SpheroidalEigenvalues" begin

    _cases = [(s, l, m) for s in (-2, 0, 2) for l in (max(abs(s), 2):max(abs(s), 2)+2)
              for m in (-l, -1, 0, 1, l) if abs(m) <= l]

    # ---- (A) characterization: spectral eigenvalue + deriv, full precision -------------------
    @testset "spectral eigenvalue == legacy (full precision)" begin
        for (s, l, m) in _cases, γ in (0.1, 0.7, 2.5)
            ref = FIXEV_spec[(s, l, m, γ)]
            got = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SphericalExpansion)
            @test got == ref
        end
    end

    @testset "eigenvalue derivative == legacy (full precision)" begin
        for (s, l, m) in _cases, γ in (0.0, 0.3, 1.4)
            ref = FIXEV_deriv[(s, l, m, γ)]
            got = SpinWeightedSpheroidalEigenvalueDeriv(s, l, m, γ)
            @test got == ref
        end
    end

    # ---- (A) characterization: Leaver CF backend (rewritten over NumericsCore CF) ------------
    @testset "Leaver eigenvalue == legacy (full precision, CF-primitive rewrite)" begin
        for (s, l, m) in _cases, γ in (0.2, 0.9, 3.0)
            ref = FIXEV_leaver[(s, l, m, γ)]
            got = Harmonics._swsh_eigenvalue_leaver(s, l, m, γ)
            @test got == ref          # bit-identical: NumericsCore CF == legacy verbatim CF (baked)
            ref2 = FIXEV_leaverfront[(s, l, m, γ)]
            got2 = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver)
            @test got2 == ref2
        end
    end

    @testset "spectral eigenvector coeffs == legacy (full precision)" begin
        for (s, l, m) in [(-2, 3, 2), (0, 4, -1), (2, 5, 0)], γ in (0.3, 1.2)
            λr, cr, lor, hir = FIXEV_coeffs[(s, l, m, γ)]
            λg, cg, log_, hig = Harmonics._swsh_spectral_coeffs(s, l, m, γ)
            @test λg == λr
            @test cg == cr
            @test (log_, hig) == (lor, hir)
        end
    end

    # ---- (B) analytic identities --------------------------------------------------------------
    @testset "γ=0 limit λ = l(l+1)-s(s+1)" begin
        for (s, l, m) in _cases
            @test SpinWeightedSpheroidalEigenvalue(s, l, m, 0.0) == l * (l + 1) - s * (s + 1)
        end
    end

    @testset "dλ/dγ|_0 = -2m s²/(l+l²) - 2m" begin
        for (s, l, m) in _cases
            @test SpinWeightedSpheroidalEigenvalueDeriv(s, l, m, 0.0) ==
                  -2 * m * s^2 / (l + l^2) - 2 * m
        end

        # The scalar l=0 mode has a removable 0/0 in the unsimplified formula.
        @test SpinWeightedSpheroidalEigenvalueDeriv(0, 0, 0, 0.0f0) === 0.0f0
        zbig = SpinWeightedSpheroidalEigenvalueDeriv(0, 0, 0, zero(BigFloat))
        @test zbig isa BigFloat && iszero(zbig)
        zc = SpinWeightedSpheroidalEigenvalueDeriv(0, 0, 0, ComplexF32(0))
        @test zc isa ComplexF32 && iszero(zc)
    end

    @testset "spectral ↔ Leaver cross-agreement" begin
        for (s, l, m) in [(-2, 3, 2), (0, 4, -2), (2, 5, 1)], γ in (0.4, 1.5)
            a = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SphericalExpansion)
            b = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver)
            @test a ≈ b rtol = 1e-9
        end
    end

    # ---- (B) precision inheritance: BigFloat γ → Leaver, beats Float64 spectral ---------------
    @testset "BigFloat γ Leaver beats Float64 spectral (precision inheritance)" begin
        setprecision(BigFloat, 256) do
            s, l, m = -2, 4, 2
            γ = big"0.5"
            λbig = SpinWeightedSpheroidalEigenvalue(s, l, m, γ)         # :Auto ⇒ Leaver for BigFloat
            @test λbig isa BigFloat
            # reference: high-precision Leaver value is self-consistent under tighter precision
            λref = setprecision(BigFloat, 512) do
                SpinWeightedSpheroidalEigenvalue(s, l, m, big"0.5")
            end
            err_big = abs(λbig - λref)
            err_f64 = abs(SpinWeightedSpheroidalEigenvalue(s, l, m, 0.5; method=:SphericalExpansion) - λref)
            @test err_big < err_f64                                    # BigFloat path is strictly tighter
            @test err_big < 1e-50
        end
    end

    # ---- (C) asymptotic-series routes through the ONE dispatcher (design §4.2: THREE backends) --
    # Tolerances MEASURED on _acases (includes the worst case l=5), not tuned:
    #  * small-γ order-6 vs the Float64 Leaver backend agrees to ≤2e-16 for γ≤1e-2 (series
    #    truncation is below the Float64 floor there) → test at rtol=1e-12.
    #  * large-γ vs the BigFloat(256) Leaver truth decays monotonically in γ, as an asymptotic
    #    series must — order-4 worst: 1.65e-5 @γ=20, 1.22e-6 @γ=30, 2.0e-7 @γ=40; order-6 worst:
    #    1.74e-6 / 5.6e-8 / 5.2e-9. The Float64 Leaver itself carries secant/CF noise up to
    #    8.5e-5 at large γ (measured, non-monotone), so the large-γ reference is BigFloat Leaver.
    _acases = [(-2, 2, 0), (-2, 3, 2), (-2, 4, -1), (0, 2, 1), (0, 5, 0), (2, 3, -2), (2, 4, 4)]

    @testset "dispatcher :SmallGamma route == series eval, ≈ Leaver in-window" begin
        for (s, l, m) in _acases
            c = SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=6)
            for γ in (1e-4, 1e-3, 1e-2)
                λr = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SmallGamma)
                @test λr == evalpoly(γ, c)                     # route wiring is EXACT
                @test λr ≈ SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver) rtol = 1e-12
            end
        end
    end

    @testset "dispatcher :LargeGamma route == series eval, ≈ BigFloat Leaver in-window" begin
        rtol4 = Dict(20.0 => 5e-5, 30.0 => 5e-6, 40.0 => 1e-6)   # measured order-4 window + ~3-5x headroom
        rtol6 = Dict(20.0 => 6e-6, 30.0 => 2e-7, 40.0 => 2e-8)   # measured order-6 window + ~3-4x headroom
        setprecision(BigFloat, 256) do
            for (s, l, m) in _acases
                cl = SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s, l, m; order=4)
                for γ in (20.0, 30.0, 40.0)
                    λref = Float64(SpinWeightedSpheroidalEigenvalue(s, l, m, BigFloat(γ); method=:Leaver))
                    λ4 = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:LargeGamma)
                    @test λ4 == evalpoly(γ, cl) / γ^4          # route wiring is EXACT
                    @test λ4 ≈ λref rtol = rtol4[γ]
                    λ6 = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:LargeGamma, numterms=6)
                    @test λ6 ≈ λref rtol = rtol6[γ]            # higher order ⇒ tighter, still in-window
                end
            end
        end
    end

    @testset "asymptotic routes leave existing routes untouched" begin
        for (s, l, m) in [(-2, 3, 2), (0, 2, 1)]
            # :Auto default for Float64 γ is STILL the spectral backend (bit-identical)
            @test SpinWeightedSpheroidalEigenvalue(s, l, m, 0.3) ==
                  SpinWeightedSpheroidalEigenvalue(s, l, m, 0.3; method=:SphericalExpansion)
            # γ=0 shortcut holds on the new routes too
            @test SpinWeightedSpheroidalEigenvalue(s, l, m, 0.0; method=:SmallGamma) ==
                  l * (l + 1) - s * (s + 1)
            @test SpinWeightedSpheroidalEigenvalue(s, l, m, 0.0; method=:LargeGamma) ==
                  l * (l + 1) - s * (s + 1)
        end
        @test_throws ArgumentError SpinWeightedSpheroidalEigenvalue(-2, 2, 0, 0.1; method=:Bogus)
    end
end

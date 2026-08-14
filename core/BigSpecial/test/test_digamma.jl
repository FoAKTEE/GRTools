# Characterization of the hand-rolled complex digamma (promotes WaSABI `_digamma`/`PolyGamma`;
# SpecialFunctions.digamma has NO Complex{BigFloat} method). Expected values are ANALYTIC
# digamma identities, or SpecialFunctions.digamma(ComplexF64) as the upstream oracle for the
# WaSABI-C complex forcing point. The shift-to-Re≥12 + 7-term Bernoulli series floors at ≈1e-16.

const _γ = Base.MathConstants.eulergamma

@testset "digamma" begin
    # ---- analytic special values (ψ(1)=-γ, ψ(1/2)=-γ-2ln2, ψ(2)=1-γ) -------------------
    @test real(digamma_big(complex(1.0))) ≈ -_γ atol=1e-13
    @test abs(imag(digamma_big(complex(1.0)))) < 1e-13
    @test real(digamma_big(complex(0.5))) ≈ -_γ - 2*log(2) atol=1e-12
    @test real(digamma_big(complex(2.0))) ≈ 1 - _γ atol=1e-13

    # ---- recurrence ψ(z+1) = ψ(z) + 1/z --------------------------------------------------
    let z = ComplexF64(2.4, 1.3)
        @test digamma_big(z + 1) ≈ digamma_big(z) + 1/z atol=1e-12
    end

    # ---- reflection ψ(1-z) - ψ(z) = π cot(πz) -------------------------------------------
    let z = ComplexF64(0.3, 0.7)
        @test digamma_big(complex(1 - z)) - digamma_big(z) ≈ π * cot(π * z) atol=1e-11
    end

    # ---- real-axis cross-check vs SpecialFunctions.digamma(real) -------------------------
    for x in (1.5, 3.7, 8.2)
        @test real(digamma_big(complex(x))) ≈ SpecialFunctions.digamma(x) atol=1e-12
    end

    # ---- WaSABI-C complex forcing point z = 3 + 2iχ/√(1-χ²), χ=0.5: vs ComplexF64 oracle ----
    let χ = 0.5, z = ComplexF64(3.0, 2*χ/sqrt(1-χ^2))
        @test digamma_big(z) ≈ SpecialFunctions.digamma(z) atol=1e-12
    end

    # ---- type-genericity: Complex{BigFloat} flows through (no Float64-literal cap) -------
    setprecision(BigFloat, 256) do
        zb = Complex{BigFloat}(big"1.0", big"0.0")
        d = digamma_big(zb)
        @test d isa Complex{BigFloat}
        @test abs(real(d) + big(_γ)) < 1e-16    # series-truncation floor (Re≥12 + B_14)
        # recurrence holds at BigFloat working precision
        zc = Complex{BigFloat}(big"2.4", big"1.3")
        @test abs(digamma_big(zc + 1) - (digamma_big(zc) + 1/zc)) < 1e-16
    end

    # ---- PolyGamma wrapper: n==0 delegates to digamma_big; else ArgumentError ------------
    @test PolyGamma(0, ComplexF64(3.0, 1.0)) == digamma_big(ComplexF64(3.0, 1.0))
    @test PolyGamma(0, 2.5) == digamma_big(complex(2.5))
    @test_throws ArgumentError PolyGamma(1, ComplexF64(3.0, 1.0))
    @test_throws ArgumentError PolyGamma(2, 1.0)
end

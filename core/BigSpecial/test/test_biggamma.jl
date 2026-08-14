# Characterization of the Complex{BigFloat}-correct Γ/logΓ/factorial stratum (consolidates
# Teukolsky `_gamma_big`/`_g`/`_lg`/`_needs_refl`, Punctures `_fact`, H1Lorenz `_fact`).
# Expected values are ANALYTIC (Γ(1/2)=√π, reflection identity, integer table) or the LIGHT
# live Teukolsky oracle at the legacy 256-bit tolerance (cov_internals: gb_* < 1e-70).

@testset "biggamma" begin
    # ---- needs_refl predicate: real non-positive Complex{BigFloat} (upstream's throw domain) ----
    setprecision(BigFloat, 256) do
        @test needs_refl(Complex{BigFloat}(big"-1.5", big"0.0")) == true
        @test needs_refl(Complex{BigFloat}(big"0.0", big"0.0")) == true   # real(z)<=0
        @test needs_refl(Complex{BigFloat}(big"-1.5", big"0.3")) == false # imag≠0
        @test needs_refl(Complex{BigFloat}(big"2.0", big"0.0")) == false  # real>0
        @test needs_refl(ComplexF64(-1.5, 0.0)) == false                  # not Complex{BigFloat}
        @test needs_refl(big"-1.5") == false                              # not Complex
        # BAKED oracle: Teukolsky._needs_refl(z) captured once (true, false, false)
        @test needs_refl(Complex{BigFloat}(big"-1.5"))            == true
        @test needs_refl(Complex{BigFloat}(big"-1.5", big"0.3")) == false
        @test needs_refl(Complex{BigFloat}(big"2.0"))            == false
    end

    # ---- gamma_big: analytic anchors ----------------------------------------------------
    setprecision(BigFloat, 256) do
        # Γ(1/2) = √π (exact); Re z = 1/2 takes the direct SpecialFunctions branch
        @test abs(gamma_big(Complex{BigFloat}(big"0.5")) - sqrt(big(pi))) < 1e-70
        # Γ(n) = (n-1)! integer table
        @test abs(gamma_big(Complex{BigFloat}(big"5.0")) - 24) < 1e-70
        @test abs(gamma_big(Complex{BigFloat}(big"1.0")) - 1) < 1e-70
        @test abs(gamma_big(Complex{BigFloat}(big"2.0")) - 1) < 1e-70
        # reflection identity Γ(z)Γ(1-z) = π/sin(πz) on a Re z<1/2 arg
        z = Complex{BigFloat}(big"-1.3", big"0.4")
        @test abs(gamma_big(z) * gamma_big(1 - z) - big(pi) / sin(big(pi) * z)) < 1e-68
        # recurrence Γ(z+1) = z Γ(z) (the cov_internals gb_recur check)
        @test abs(gamma_big(z + 1) - z * gamma_big(z)) < 1e-70
        # Re z > 1/2: matches SpecialFunctions.gamma exactly (no reflection there)
        zr = Complex{BigFloat}(big"2.3", big"1.1")
        @test gamma_big(zr) == SpecialFunctions.gamma(zr)
        # edge: real-negative half-integer -3/2 → finite reflection value 4√π/3
        gm = gamma_big(Complex{BigFloat}(big"-1.5", big"0.0"))
        @test isfinite(real(gm)) && isfinite(imag(gm))
        @test abs(gm - 4 * sqrt(big(pi)) / 3) < 1e-68

        # Exact nonpositive integers are poles, not huge finite reflection artifacts. The public
        # convention is direction-independent +Inf for both Γ and logΓ.
        for n in (0, -1, -2, -7)
            gr = gamma_big(BigFloat(n))
            gc = gamma_big(Complex{BigFloat}(BigFloat(n), zero(BigFloat)))
            lg = loggamma_big(Complex{BigFloat}(BigFloat(n), zero(BigFloat)))
            @test gr isa BigFloat && gr == BigFloat(Inf)
            @test gc isa Complex{BigFloat} && real(gc) == BigFloat(Inf) && iszero(imag(gc))
            @test lg isa Complex{BigFloat} && real(lg) == BigFloat(Inf) && iszero(imag(lg))
        end
        @test gamma_big(-1.0) === Inf
        @test gamma_big(ComplexF64(-2.0, 0.0)) === ComplexF64(Inf, 0.0)
        @test loggamma_big(-2.0) === ComplexF64(Inf, 0.0)
        @test isfinite(gamma_big(Complex{BigFloat}(big"-2.0", eps(BigFloat))))
    end

    # ---- gamma_big: BAKED Teukolsky._gamma_big oracle (captured once; bit-identical, same
    #      reflection). RHS literals == the legacy 256-bit values -> characterization preserved.
    setprecision(BigFloat, 256) do
        # Teukolsky._gamma_big(Complex{BigFloat}(big"0.5"))
        @test gamma_big(Complex{BigFloat}(big"0.5")) == Complex{BigFloat}(parse(BigFloat, "1.772453850905516027298167483341145182797549456122387128213807789852911284591025"), parse(BigFloat, "0.0"))
        # Teukolsky._gamma_big(Complex{BigFloat}(big"5.0"))
        @test gamma_big(Complex{BigFloat}(big"5.0")) == Complex{BigFloat}(parse(BigFloat, "24.0"), parse(BigFloat, "0.0"))
        # Teukolsky._gamma_big(Complex{BigFloat}(big"-1.3", big"0.4"))
        @test gamma_big(Complex{BigFloat}(big"-1.3", big"0.4")) == Complex{BigFloat}(parse(BigFloat, "1.088661863120153869523304447696766072433931845698365216462672772677391930866162"), parse(BigFloat, "1.112780331676832146518694134301417749118230313358712934756661414915144425085234"))
        # Teukolsky._gamma_big(Complex{BigFloat}(big"-1.5"))
        @test gamma_big(Complex{BigFloat}(big"-1.5")) == Complex{BigFloat}(parse(BigFloat, "2.363271801207354703064223311121526910396732608163182837618410386470548379454688"), parse(BigFloat, "-0.0"))
        # Teukolsky._gamma_big(Complex{BigFloat}(big"2.3", big"1.1"))
        @test gamma_big(Complex{BigFloat}(big"2.3", big"1.1")) == Complex{BigFloat}(parse(BigFloat, "0.6429386029517632935962313816547722859569678948677798071204279936492533768027916"), parse(BigFloat, "0.562237642153628482290268019682549534924175164318665157864908983422992585861802"))
    end

    # ---- loggamma_big: exp(loggamma_big(z)) == gamma_big(z) (branch-robust grounding) ----
    setprecision(BigFloat, 256) do
        for z in (Complex{BigFloat}(big"2.3", big"1.1"), Complex{BigFloat}(big"0.3", big"0.1"),
                  Complex{BigFloat}(big"-1.3", big"0.4"))
            @test abs(exp(loggamma_big(z)) - gamma_big(z)) < 1e-66 * abs(gamma_big(z))
        end
        # Re z > 1/2: equals SpecialFunctions.loggamma directly
        zr = Complex{BigFloat}(big"3.2", big"0.7")
        @test loggamma_big(zr) == SpecialFunctions.loggamma(zr)
    end

    # ---- fact_big = gamma_big(n+1): analytic factorial table -----------------------------
    @test fact_big(0) ≈ 1 rtol=0 atol=1e-12
    @test fact_big(5) ≈ 120 rtol=0 atol=1e-9
    @test fact_big(6) ≈ 720 rtol=1e-12
    setprecision(BigFloat, 256) do
        @test abs(fact_big(big"5.0") - 120) < 1e-70
        @test abs(fact_big(big"5.0") - gamma_big(Complex{BigFloat}(big"6.0"))) < 1e-66
    end

    # ---- factorial_big(T, n): exact integer product (H1Lorenz overflow-safe), distinct prim ----
    @test factorial_big(BigFloat, 0) == big(1)
    @test factorial_big(BigFloat, 1) == big(1)
    @test factorial_big(BigFloat, 5) == big(120)
    # overflow-safe at high l: exact == BigFloat factorial; Int64 factorial would overflow
    @test factorial_big(BigFloat, 25) == prod(big(2):big(25))
    setprecision(BigFloat, 256) do
        @test abs(factorial_big(BigFloat, 6) - real(gamma_big(Complex{BigFloat}(big"7.0")))) < 1e-66
    end
end

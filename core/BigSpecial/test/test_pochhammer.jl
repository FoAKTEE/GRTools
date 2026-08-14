# Characterization of `pochhammer` (the canonical rising/reciprocal-falling Pochhammer that
# consolidates Teukolsky `_pochhammer` + Punctures `_poch`). Expected values are ANALYTIC
# closed forms (never == new-code); plus a LIGHT live cross-check vs the Teukolsky oracle.

@testset "pochhammer" begin
    # ---- analytic closed forms ---------------------------------------------------------
    # (a)_0 = 1
    @test pochhammer(2.5, 0) == 1.0
    @test pochhammer(7.3, 0) == 1.0

    # (a)_n = ∏_{k=0}^{n-1}(a+k): explicit literal product
    @test pochhammer(2.5, 3) ≈ 2.5 * 3.5 * 4.5 rtol=0 atol=1e-13
    @test pochhammer(3.0, 4) == 3.0 * 4.0 * 5.0 * 6.0

    # (1)_n = n!
    for n in 1:8
        @test pochhammer(1.0, n) ≈ factorial(n) rtol=0 atol=1e-9
    end

    # (1/2)_n = (2n)!/(4^n n!)
    for n in 1:7
        @test pochhammer(0.5, n) ≈ factorial(2n) / (4.0^n * factorial(n)) rtol=1e-12
    end

    # reciprocal-falling: (a)_{-n} = 1/∏_{k=1}^{n}(a-k)
    @test pochhammer(5.0, -2) ≈ 1 / (4.0 * 3.0) rtol=0 atol=1e-14   # 1/((5-1)(5-2))
    @test pochhammer(5.0, -3) ≈ 1 / (4.0 * 3.0 * 2.0) rtol=1e-13

    # ratio identity (a)_{n+1} = (a+n)(a)_n
    let a = 1.7
        for n in 0:6
            @test pochhammer(a, n + 1) ≈ (a + n) * pochhammer(a, n) rtol=1e-12
        end
    end

    # ---- type-genericity: Complex{BigFloat} flows through (precision inherited) ----------
    setprecision(BigFloat, 256) do
        ab = Complex{BigFloat}(big"1.8", big"0.3")
        # (a)_3 explicit product, full precision
        @test pochhammer(ab, 3) ≈ ab * (ab + 1) * (ab + 2)
        @test pochhammer(ab, 0) == one(ab)
        @test eltype(real(pochhammer(ab, 4))) == BigFloat
        # reciprocal branch in BigFloat
        @test pochhammer(ab, -2) ≈ one(ab) / ((ab - 1) * (ab - 2))
    end

    # ---- BAKED FIXTURE oracle: Teukolsky `_pochhammer` (captured once; its closure is
    #      KerrGeodesics+SWSH, light). RHS is now a baked literal == the legacy value,
    #      preserving the bit-identical characterization without a live packages/* dep.
    @test pochhammer(2.5, 3) == 39.375                     # Teukolsky._pochhammer(2.5, 3)
    @test pochhammer(7.3, 0) == 1.0                        # Teukolsky._pochhammer(7.3, 0)
    @test pochhammer(5.0, -2) == 0.08333333333333333       # Teukolsky._pochhammer(5.0, -2)
    setprecision(BigFloat, 256) do
        ab = Complex{BigFloat}(big"1.8", big"0.3")
        # Teukolsky._pochhammer(ab, 5)
        @test pochhammer(ab, 5) == Complex{BigFloat}(parse(BigFloat, "489.0907799999999999999999999999999999999999999999999999999999999999999999999958"), parse(BigFloat, "245.229030000000000000000000000000000000000000000000000000000000000000000000002"))
        # Teukolsky._pochhammer(ab, -3)
        @test pochhammer(ab, -3) == Complex{BigFloat}(parse(BigFloat, "1.694249881195892477117295811896940019421888882001694249881195892477117295811906"), parse(BigFloat, "2.004173639951238661955825533585404657121014897002004173639951238661955825533546"))
    end
end

using Test
using AccretionDisc.Modes: pl00_tauL_inv, pl00_taue_inv, cn06_tauL_inv, cn06_taue_inv

@testset "FR24 v1 literature prescriptions (52)-(55)" begin
    @testset "circular normalisations and physical-e input" begin
        @test pl00_tauL_inv(0, 0.05, 0.5) ≈ 0.3665
        @test pl00_taue_inv(0, 0.05, 0.5) ≈ 85.2
        @test cn06_tauL_inv(0, 0.05, 1.5) ≈ 0.10875
        for (Q, expected) in ((0.1, 156.0), (0.2, 78.0), (1.0, 15.6))
            @test cn06_taue_inv(0, 0.05, Q) ≈ expected
            @test cn06_taue_inv(0.05, 0.05, Q) ≈ expected * 4 / 5
            @test cn06_taue_inv(0.10, 0.05, Q) ≈ expected / 3
        end
        @test pl00_taue_inv(0.05, 0.05, 0.5) ≈ 68.16
        @test pl00_taue_inv(0.10, 0.05, 0.5) ≈ 28.4
        # Independent exact fractions from the printed fits at e/h=2.
        @test pl00_tauL_inv(0.10, 0.05, 0.5) ≈ -425382592947//1124458082000
        @test cn06_tauL_inv(0.10, 0.05, 1.5) ≈ -50488793433//449783232800
        @test pl00_tauL_inv(0, 1, 1) isa Float64
        @test cn06_taue_inv(0, 1, 1) isa Float64
    end

    @testset "reversal, softening, and similarity scalings" begin
        for f in (pl00_tauL_inv, cn06_tauL_inv)
            third = f === pl00_tauL_inv ? 0.5 : 1.5
            @test f(0.054, 0.05, third) > 0
            @test abs(f(0.055, 0.05, third)) < 1e-15
            @test f(0.056, 0.05, third) < 0
        end
        for e in (0.0, 0.025, 0.08, 0.12), b in (0.3, 0.4, 0.5, 0.8)
            @test pl00_tauL_inv(e, 0.05, b) / pl00_tauL_inv(e, 0.05, 0.5) ≈ (2b)^(-7/4)
            @test pl00_taue_inv(e, 0.05, b) / pl00_taue_inv(e, 0.05, 0.5) ≈ (2b)^(-5/2)
            @test pl00_tauL_inv(2e, 0.10, b) ≈ 2pl00_tauL_inv(e, 0.05, b)
            @test pl00_taue_inv(2e, 0.10, b) ≈ pl00_taue_inv(e, 0.05, b)/2
        end
        @test cn06_tauL_inv(0.02, 0.05, 2.0) / cn06_tauL_inv(0.02, 0.05, 0.0) ≈ 49/27
        # Equation (55) is uncapped: FR24 discusses a separate cap, but does not plot it.
        @test cn06_taue_inv(0.01, 0.05, 0.1) > cn06_taue_inv(0.04, 0.05, 0.1)
    end

    @testset "inherited BigFloat precision and mixed inputs" begin
        for bits in (128, 256)
            setprecision(BigFloat, bits) do
                e, h, b = BigFloat(1)/10, BigFloat(1)/20, BigFloat(1)/2
                # At e/h=2 all four expectations below are exact rationals.
                shape = (1 - (big(20)//11)^4) / (1 + (big(20)//13)^5)
                expected = (big(733)//2000*shape, big(142)//5,
                            big(87)//800*shape, big(26)//1)
                actual = (pl00_tauL_inv(e,h,b), pl00_taue_inv(e,h,b),
                          cn06_tauL_inv(e,h,3//2), cn06_taue_inv(e,h,1//5))
                for (got, want) in zip(actual, expected)
                    @test got isa BigFloat
                    @test precision(got) == bits
                    @test isapprox(got, BigFloat(want); rtol=32eps(BigFloat))
                end
                @test isapprox(cn06_taue_inv(e,h,1//10), BigFloat(52); rtol=8eps(BigFloat))
                for f in (pl00_tauL_inv, pl00_taue_inv)
                    mixed = f(0.125, 0.0625, BigFloat(3)/8)
                    full = f(BigFloat(1)/8, BigFloat(1)/16, BigFloat(3)/8)
                    @test mixed == full
                    @test precision(mixed) == bits
                end
                # Independent root-based softening factors at non-binary b=0.3.
                soft = BigFloat(3)/5
                @test isapprox(pl00_tauL_inv(e,h,3//10),
                               BigFloat(expected[1])/sqrt(sqrt(soft))^7; rtol=32eps(BigFloat))
                @test isapprox(pl00_taue_inv(e,h,3//10),
                               BigFloat(expected[2])/sqrt(soft)^5; rtol=32eps(BigFloat))
            end
        end
    end
end

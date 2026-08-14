# Grounding for source_amplitudes (Sourcesamplitudes, Sourceincoordinates.nb SC l.1114).
#
# ORACLE: frozen `wolframscript` replication of the upstream `Sourcesamplitudes` composition —
# same +ω/−ω join, minus-mode conjugation σ=(−1)^(l+kloc), and (mi→−m, ω→−ω) substitution
# transcribed verbatim from the notebook source, wired to the BHPToolkit `Teukolsky` paclet
# `RadHzpotmode` and the extracted analytic `source_kernels` (kernels + RadHZ are independently
# grounded upstream). Agreement is MST-truncation limited (~1e-7 relative).
#
# Reference values (circular a=0 p=10 l=3 m=2 at r=7; spherical a=0 p=12 x=0.5 l=3 m=2 k=1 at r=9):

@testset "source_amplitudes" begin
    @testset "circular oracle (a=0,p=10,l=3,m=2) at r=7" begin
        S = source_amplitudes(3, 2, 0, 0, 0.0, 10.0, 0.0, 1.0)
        @test S.omega ≈ 0.06324555320336758 rtol=1e-10
        @test S.Sl(7.0)    ≈ (-0.030996469282640743 - 0.01156955943612227im)  rtol=1e-7
        @test S.Sn(7.0)    ≈ (0.011070146555085172 - 0.0028453970890759614im) rtol=1e-7
        @test S.Sm(7.0)    ≈ (-0.0077181785535011364 - 0.0005276647998963633im) rtol=1e-7
        @test S.Smdag(7.0) ≈ (0.007718178553501135 + 0.0005276647998963633im) rtol=1e-7
        @test S.Snorm(7.0) == 0
    end

    @testset "spherical oracle (a=0,p=12,x=0.5,l=3,m=2,k=1) at r=9" begin
        S = source_amplitudes(3, 2, 0, 1, 0.0, 12.0, 0.0, 0.5)
        @test S.omega ≈ 0.07216878364870322 rtol=1e-10
        @test S.Sm(9.0)    ≈ (0.01293272199462648 + 6.063064914813834e-5im)  rtol=1e-7
        @test S.Smdag(9.0) ≈ (0.012932721994626481 + 6.063064914813834e-5im) rtol=1e-7
    end

    @testset "l=2 dominant mode runs (|s|>l harmonics vanish, no throw)" begin
        S = source_amplitudes(2, 2, 0, 0, 0.0, 10.0, 0.0, 1.0)
        @test isfinite(abs(S.Sm(7.0)))
        @test isfinite(abs(S.Sl(7.0)))
    end
end

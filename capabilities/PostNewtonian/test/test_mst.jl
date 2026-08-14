# Characterization of PostNewtonian.MST (symbolic Mano-Suzuki-Takasugi engine).
#
# GROUNDING:
#  * Analytic anchor: the eps->0 (Schwarzschild) limit of the renormalized angular momentum is
#    nu -> l, so the eta^0 coefficient of nuMST(s=-2,l=2,m=2,a=0) is exactly l=2.
#  * Captured one-shot fixtures from the LEGACY heavy TeukolskyPN.MSTCoefficients (run ONCE, baked as
#    exact-rational literals here — NOT a live test dep). Reproduced commands in the comments.
#  * ES Laurent round-trip identity: _es_mul(a, _es_inv(a)) == 1 (exact).
#  * imagunit numeric substitution: 𝕚^2 = -1.

using PostNewtonian: MST
using PostNewtonian.MST: ES, MSTCoefficients, imagunit,
    _es_from, _es_mul, _es_inv, _es_coeff, _es_const
using Symbolics
import Symbolics

@testset "MST" begin
    @variables ω

    # --- ES Laurent round-trip: a * (1/a) == 1 (exact rationals) -----------------------------
    @testset "ES round-trip" begin
        ord = 8
        a = _es_from(Complex{Rational{BigInt}}[3, 1, -2, 5, 7], ord)   # 3 + eps - 2eps^2 + ...
        one_es = _es_mul(a, _es_inv(a, ord), ord)
        @test _es_coeff(one_es, 0) == 1
        for k in 1:ord
            @test _es_coeff(one_es, k) == 0
        end
        # a Laurent (negative-lead) series inverts to a positive-lead one
        b = MST._es_norm(-2, Complex{Rational{BigInt}}[4, 0, 1, 3], ord)
        ob = _es_mul(b, _es_inv(b, ord), ord)
        @test _es_coeff(ob, 0) == 1
        @test _es_coeff(ob, 1) == 0
    end

    # --- imagunit: 𝕚^2 = -1 -------------------------------------------------------------------
    @testset "imagunit" begin
        iu = imagunit()
        @test Symbolics.value(Symbolics.substitute(iu^2, Dict(iu => im))) == -1
        @test Symbolics.value(Symbolics.substitute(iu, Dict(iu => im))) == im
    end

    # --- nu -> l Schwarzschild limit + full captured fixture ----------------------------------
    # LEGACY ORACLE (one-shot capture):
    #   MSTCoefficients(-2,2,2, 0//1, :ω, (:η,6))[:νMST]
    #     => PNSeries(:η, [2//1,0,0,0//1,0,0,(-214//105)*ω^2], nmin=0, nmax=7, den=1)
    @testset "nuMST Schwarzschild" begin
        res = MSTCoefficients(-2, 2, 2, 0 // 1, :ω, (:η, 6))
        ν = res[:νMST]
        @test ν.nmin == 0
        @test ν.nmax == 7
        @test ν.den == 1
        # eps->0 limit: leading coeff == l == 2 (analytic anchor)
        @test isequal(ν.coeffs[1], Num(2))
        # captured higher coeffs
        @test isequal(ν.coeffs[2], Num(0))
        @test isequal(ν.coeffs[4], Num(0))
        @test isequal(Symbolics.expand(ν.coeffs[7] - (-214 // 105) * ω^2), Num(0))
    end

    # --- aMST coefficients captured fixture ---------------------------------------------------
    # LEGACY ORACLE (one-shot capture), s=-2,l=2,m=2,a=0,(:η,6):
    #   aMST[1]  => PNSeries(:η, [(5//3)*ω*𝕚, 0, 0, (10//9)*ω^2], nmin=3, nmax=7, den=1)
    #   aMST[-1] => pure O[] (nmin==nmax==7)
    @testset "aMST coefficients" begin
        res = MSTCoefficients(-2, 2, 2, 0 // 1, :ω, (:η, 6))
        a1 = res[:aMST][1]
        @test a1.nmin == 3
        @test a1.nmax == 7
        iu = imagunit()
        # leading: (5/3) ω 𝕚
        @test isequal(Symbolics.expand(a1.coeffs[1] - (5 // 3) * ω * iu), Num(0))
        @test isequal(Symbolics.expand(a1.coeffs[4] - (10 // 9) * ω^2), Num(0))
        am1 = res[:aMST][-1]
        @test am1.nmin == am1.nmax     # pure O[]
    end
end

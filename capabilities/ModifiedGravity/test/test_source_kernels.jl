# Grounding for SourceKernels (Sourceace{l,n,m,m†}2± / Sourcenorm2±, Sourceincoordinates.nb).
#
# ORACLE: each kernel is the closed-form expression from the notebook's saved "Analytical formulas"
# output cells. The reference numbers are a `wolframscript` evaluation of those SAME extracted
# expressions at a fixed point: MB=1, R=8, ω=0.05, mi=2, li=3, rhz=0.1-0.2i, drhz=-0.03+0.04i,
# rhzb=0.2+0.1i, drhzb=0.05-0.02i, yv(s,m)={}_sY_{3,m}(π/2). This is structural/exact agreement of
# the same analytic algebra (frozen literals transcribed from the legacy source_kernels_test.jl),
# not a float fit. `yv` is independently validated against the oracle's SWSH values (~1e-12), so
# the kernel test isolates the kernel algebra.

using Harmonics: SpinWeightedSphericalHarmonicY
import ModifiedGravity: _source_kernel_l_plus, _source_kernel_l_minus,
    _source_kernel_n_plus, _source_kernel_n_minus, _source_kernel_m_plus, _source_kernel_m_minus,
    _source_kernel_mdag_plus, _source_kernel_mdag_minus, _source_kernel_norm_plus,
    _source_kernel_norm_minus

@testset "SourceKernels" begin
    MB, R, ω, mi, li = 1, 8.0, 0.05, 2, 3
    rhz   = 0.1 - 0.2im
    drhz  = -0.03 + 0.04im
    rhzb  = 0.2 + 0.1im
    drhzb = 0.05 - 0.02im
    yv(s, m) = SpinWeightedSphericalHarmonicY(s, li, m, π/2, 0.0)

    @testset "Julia SWSH matches oracle yv (sanity)" begin
        @test yv(2, 2)  ≈ 0.3731763325901154  atol=1e-12
        @test yv(3, 2)  ≈ 0.2285228997322329  atol=1e-12
        @test yv(-2, -2) ≈ 0.3731763325901154 atol=1e-12
    end

    @testset "kernels == wolframscript analytic cells" begin
        @test _source_kernel_l_plus(MB,R,ω,mi,li,rhz,drhz,yv) ≈
            0.9328014679884595 + 0.5497657688511733im rtol=1e-10
        @test _source_kernel_l_minus(MB,R,ω,mi,li,rhzb,drhzb,yv) ≈
            0.3729089161824997 - 1.0570820835929584im rtol=1e-10
        @test _source_kernel_n_plus(MB,R,ω,mi,li,rhz,drhz,yv) ≈
            -0.2981529460652004 - 0.22797805172240812im rtol=1e-10
        @test _source_kernel_n_minus(MB,R,ω,mi,li,rhzb,drhzb,yv) ≈
            0.07198150633636399 + 0.37102979073123143im rtol=1e-10
        @test _source_kernel_m_plus(MB,R,ω,mi,li,rhz,drhz,yv) ≈
            -0.4733110404149359 + 0.7317164815703933im rtol=1e-10
        @test _source_kernel_m_minus(MB,R,ω,mi,li,rhzb,drhzb,yv) ≈
            -0.8246109854774907 - 0.37704512699972786im rtol=1e-10
        @test _source_kernel_mdag_plus(MB,R,ω,mi,li,rhz,drhz,yv) ≈
            -0.4082484066901961 + 0.8170245644110976im rtol=1e-10
        @test _source_kernel_mdag_minus(MB,R,ω,mi,li,rhzb,drhzb,yv) ≈
            -0.8010384109883697 + 0.09597258076761492im rtol=1e-10
        @test _source_kernel_norm_plus(MB,R,ω,mi,li,rhz,drhz,yv) ≈
            -0.28916222302036626 - 1.5796764162952457im rtol=1e-10
        @test _source_kernel_norm_minus(MB,R,ω,mi,li,rhzb,drhzb,yv) ≈
            1.0150396246451137 + 4.874808565357167im rtol=1e-10
    end
end

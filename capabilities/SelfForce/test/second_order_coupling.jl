# Grounded tests for SelfForce.SecondOrder.coupling (3-harmonic spin-weighted angular coupling).
#
# GROUNDING — pure delegation + analytic, no heavy legacy dep:
#
# (1) SecondOrder.coupling is the FRONT DOOR onto Wigner.coupling (the canonical home, already
#     GSL-grounded vs gsl_sf_coupling_3j in the Wigner core). We assert OBJECT IDENTITY (=== the
#     core function) so there is no second implementation to drift, and re-check the VALUES against
#     the core on a representative (lc,mc,sc;la,ma,sa;lb,mb,sb) grid at Float64 + BigFloat.
#
# (2) Analytic 3j selection-rule zeros (the legacy GSL contract): coupling → 0 whenever a 3j slot is
#     out of range (|m|>j or |s|>l), the m-sum -mc+ma+mb ≠ 0, the s-sum sc-sa-sb ≠ 0, or the
#     (lc,la,lb) triangle inequality is violated.
#
# (3) Closed-form anchor: with all spins 0 and m's 0, coupling reduces to the Gaunt-type product
#     √((2la+1)(2lb+1)(2lc+1)/4π) · 3j(lc la lb;0 0 0)^2, which we cross-check against an explicit
#     value computed from Wigner.wigner3j.

using Test
using SelfForce.SecondOrder: coupling
using Wigner: coupling as core_coupling, wigner3j

@testset "coupling is the Wigner re-export (identity + values)" begin
    @test coupling === core_coupling
    # representative non-trivial couplings (selection rules satisfied)
    cases = [
        (2, 0, 0, 2, 0, 0, 2, 0, 0),
        (2, 1, -2, 2, 1, 2, 2, 0, 0),
        (3, 1, 0, 2, 1, -2, 1, 0, 2),
        (4, 2, 0, 2, 2, 2, 2, 0, -2),
        (2, 0, -2, 1, 1, 0, 3, -1, 2),
    ]
    for c in cases
        @test coupling(c...) == core_coupling(c...)
        @test coupling(Float64, c...) == core_coupling(Float64, c...)
    end
end

@testset "coupling BigFloat precision path" begin
    setprecision(BigFloat, 256) do
        v = coupling(BigFloat, 2, 1, -2, 2, 1, 2, 2, 0, 0)
        @test v isa BigFloat
        @test v == core_coupling(BigFloat, 2, 1, -2, 2, 1, 2, 2, 0, 0)
        # BigFloat carries more precision than the Float64 path (not bit-identical downcast)
        @test Float64(v) ≈ coupling(Float64, 2, 1, -2, 2, 1, 2, 2, 0, 0) rtol = 1e-12
    end
end

@testset "coupling selection-rule zeros (GSL contract)" begin
    # |m| > j slot
    @test coupling(2, 3, 0, 2, 0, 0, 2, 0, 0) == 0.0
    # |s| > l slot
    @test coupling(2, 0, 3, 2, 0, 0, 2, 0, 0) == 0.0
    # m-sum -mc+ma+mb ≠ 0  (here -1+0+0)
    @test coupling(2, 1, 0, 2, 0, 0, 2, 0, 0) == 0.0
    # s-sum sc-sa-sb ≠ 0   (here 1-0-0)
    @test coupling(2, 0, 1, 2, 0, 0, 2, 0, 0) == 0.0
    # triangle violation: lc=5 > la+lb = 2+2
    @test coupling(5, 0, 0, 2, 0, 0, 2, 0, 0) == 0.0
end

@testset "coupling Gaunt closed-form anchor (all spins/m zero)" begin
    for (lc, la, lb) in ((2, 2, 2), (4, 2, 2), (0, 1, 1))
        norm = sqrt((2la + 1) * (2lb + 1) * (2lc + 1) / (4 * pi))
        w = wigner3j(Float64, lc, la, lb, 0, 0, 0)
        @test coupling(lc, 0, 0, la, 0, 0, lb, 0, 0) ≈ norm * w * w rtol = 1e-13
    end
end

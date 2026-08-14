# Grounded tests for Wigner.clebsch_gordan (Racah CG, promoted from KerrLorenzCirc's `cg`).
#
# Grounding sources:
#  (1) KerrLorenzCirc literal cg values baked from packages/KerrLorenzCirc/test/cov_internal.jl:92-96
#      (rtol=1e-12). KLC is HEAVY (Teukolsky/SWSH/QuadGK) so we do NOT live-dep it — these are
#      dead-simple closed forms copied verbatim.
#  (2) 3j↔CG relation CG(j1m1,j2m2|j3m3) = (-1)^(j1-j2+m3)·√(2j3+1)·3j(j1 j2 j3; m1 m2 -m3),
#      cross-checked against the in-package wigner3j.
#  (3) WignerSymbols.clebschgordan over a grid (light TEST-ONLY dep).

using Wigner, Test
using WignerSymbols: clebschgordan as wsCG

@testset "clebsch_gordan KLC baked literals (cov_internal.jl:92-96)" begin
    @test clebsch_gordan(1, 0, 1, 0, 0, 0) ≈ -1 / sqrt(3) rtol=1e-12   # ⟨10;10|00⟩
    @test clebsch_gordan(1, 1, 1, -1, 0, 0) ≈ 1 / sqrt(3) rtol=1e-12
    @test clebsch_gordan(2, 0, 2, 0, 0, 0) ≈ sqrt(1 / 5) rtol=1e-12    # ⟨20;20|00⟩=1/√5
    @test clebsch_gordan(1, 1, 1, 1, 2, 2) ≈ 1.0 rtol=1e-12            # stretched state
    @test clebsch_gordan(1, 0, 1, 0, 1, 0) == 0.0                      # parity selection rule
end

@testset "clebsch_gordan ↔ wigner3j relation" begin
    js = 0:4
    for j1 in js, j2 in js, j3 in abs(j1 - j2):(j1 + j2)
        for m1 in -j1:j1, m2 in -j2:j2
            m3 = m1 + m2
            abs(m3) > j3 && continue
            sgn = iseven(j1 - j2 + m3) ? 1.0 : -1.0
            rel = sgn * sqrt(2j3 + 1) * wigner3j(Float64, j1, j2, j3, m1, m2, -m3)
            @test clebsch_gordan(j1, m1, j2, m2, j3, m3) ≈ rel atol=1e-13
        end
    end
end

@testset "clebsch_gordan vs WignerSymbols oracle" begin
    js = 0:4
    for j1 in js, j2 in js, j3 in abs(j1 - j2):(j1 + j2)
        for m1 in -j1:j1, m2 in -j2:j2
            m3 = m1 + m2
            abs(m3) > j3 && continue
            @test clebsch_gordan(j1, m1, j2, m2, j3, m3) ≈ Float64(wsCG(j1, m1, j2, m2, j3, m3)) atol=1e-13
        end
    end
end

@testset "clebsch_gordan selection-rule zeros" begin
    @test clebsch_gordan(1, 1, 1, 0, 0, 0) == 0.0    # m1+m2 ≠ m3
    @test clebsch_gordan(1, 0, 1, 0, 3, 0) == 0.0    # triangle violated
    @test clebsch_gordan(1, 2, 1, 0, 1, 0) == 0.0    # |m1|>j1
end

@testset "clebsch_gordan BigFloat precision inherited" begin
    setprecision(BigFloat, 256) do
        gb = clebsch_gordan(BigFloat, 1, 0, 1, 0, 0, 0)
        @test gb isa BigFloat
        @test gb ≈ -1 / sqrt(big(3)) rtol=1e-70
    end
end

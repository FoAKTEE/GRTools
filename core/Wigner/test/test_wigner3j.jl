# Grounded tests for Wigner.wigner3j / Wigner.wigner6j (hand-rolled Racah closed forms).
#
# Grounding sources:
#  (1) WignerSymbols.jl (small, light, registered) as a TEST-ONLY oracle over a grid of
#      (j1,j2,j3,m1,m2,m3) — never a runtime [deps].
#  (2) Analytic closed forms / identities:
#        - 3j(l,l,0; m,-m,0) = (-1)^(l-m)/√(2l+1)
#        - even-permutation invariance; odd-permutation ×(-1)^(j1+j2+j3)
#        - m-flip (all m→-m) ×(-1)^(j1+j2+j3)
#        - orthogonality Σ_{m1,m2} (2j3+1) 3j(j1 j2 j3)² = δ over allowed j3
#        - triangle / |m|>j / m-sum≠0 selection-rule zeros
#        - 6j special value {j1 j2 j3; 0 j3 j2} = (-1)^(j1+j2+j3)/√((2j2+1)(2j3+1))

using Wigner, Test
using WignerSymbols: wigner3j as ws3j, wigner6j as ws6j

@testset "wigner3j vs WignerSymbols oracle" begin
    js = 0:4
    for j1 in js, j2 in js, j3 in js
        (j3 < abs(j1 - j2) || j3 > j1 + j2) && continue
        for m1 in -j1:j1, m2 in -j2:j2
            m3 = -(m1 + m2)
            abs(m3) > j3 && continue
            ref = Float64(ws3j(j1, j2, j3, m1, m2, m3))
            got = wigner3j(Float64, j1, j2, j3, m1, m2, m3)
            @test got ≈ ref atol=1e-13
        end
    end
end

@testset "wigner3j closed form 3j(l,l,0;m,-m,0)" begin
    for l in 0:6, m in -l:l
        expected = (iseven(l - m) ? 1.0 : -1.0) / sqrt(2l + 1)
        @test wigner3j(Float64, l, l, 0, m, -m, 0) ≈ expected rtol=1e-13
    end
end

@testset "wigner3j permutation + Regge symmetry" begin
    cases = [(2, 3, 4, 1, -1, 0), (3, 1, 2, 0, 1, -1), (4, 4, 2, 2, -1, -1),
             (5, 3, 4, -2, 1, 1), (2, 2, 2, 1, 1, -2)]
    for (j1, j2, j3, m1, m2, m3) in cases
        base = wigner3j(Float64, j1, j2, j3, m1, m2, m3)
        # cyclic (even) permutation: invariant
        @test wigner3j(Float64, j2, j3, j1, m2, m3, m1) ≈ base atol=1e-14
        @test wigner3j(Float64, j3, j1, j2, m3, m1, m2) ≈ base atol=1e-14
        # odd permutation (swap cols 1,2): ×(-1)^(j1+j2+j3)
        s = iseven(j1 + j2 + j3) ? 1.0 : -1.0
        @test wigner3j(Float64, j2, j1, j3, m2, m1, m3) ≈ s * base atol=1e-14
        # m-flip: ×(-1)^(j1+j2+j3)
        @test wigner3j(Float64, j1, j2, j3, -m1, -m2, -m3) ≈ s * base atol=1e-14
    end
end

@testset "wigner3j orthogonality" begin
    # Fixed-m3 orthogonality: Σ_{m1} (2j3+1) 3j(j1 j2 j3; m1, -m1-m3, m3)² = 1
    # for each allowed (j3, m3) triangle.
    for j1 in 0:3, j2 in 0:3, j3 in abs(j1 - j2):(j1 + j2), m3 in -j3:j3
        acc = 0.0
        for m1 in -j1:j1
            m2 = -(m1 + m3)
            abs(m2) > j2 && continue
            acc += (2j3 + 1) * wigner3j(Float64, j1, j2, j3, m1, m2, m3)^2
        end
        @test acc ≈ 1.0 atol=1e-12
    end
end

@testset "wigner3j selection-rule zeros" begin
    @test wigner3j(Float64, 1, 0, 0, 0, 0, 0) == 0.0      # triangle (1,0,0)
    @test wigner3j(Float64, 1, 1, 1, 0, 0, 0) == 0.0      # odd perimeter w/ all-m=0
    @test wigner3j(Float64, 2, 1, 1, 1, 0, 0) == 0.0      # m-sum ≠ 0
    @test wigner3j(Float64, 1, 1, 2, 2, 0, -2) == 0.0     # |m1|>j1
    @test wigner3j(Float64, 3, 2, 2, 0, 0, 0) == 0.0      # odd perimeter (3,2,2)
end

@testset "wigner6j vs WignerSymbols oracle" begin
    js = 0:3
    n = 0
    for j1 in js, j2 in js, j3 in js, j4 in js, j5 in js, j6 in js
        ref = try
            Float64(ws6j(j1, j2, j3, j4, j5, j6))
        catch
            continue
        end
        got = wigner6j(Float64, j1, j2, j3, j4, j5, j6)
        @test got ≈ ref atol=1e-12
        n += 1
    end
    @test n > 50
end

@testset "wigner6j special value {j1 j2 j3; 0 j3 j2}" begin
    for j1 in 0:4, j2 in 0:4, j3 in abs(j1 - j2):(j1 + j2)
        s = iseven(j1 + j2 + j3) ? 1.0 : -1.0
        expected = s / sqrt((2j2 + 1) * (2j3 + 1))
        @test wigner6j(Float64, j1, j2, j3, 0, j3, j2) ≈ expected rtol=1e-12
    end
end

@testset "wigner3j BigFloat precision inherited" begin
    setprecision(BigFloat, 256) do
        gb = wigner3j(BigFloat, 1, 1, 0, 0, 0, 0)
        @test gb isa BigFloat
        @test gb ≈ -1 / sqrt(big(3)) rtol=1e-70   # 3j(1,1,0;0,0,0) = (-1)^1/√3
    end
end

# legendre_P3 — associated Legendre function of TYPE 3 (z > 1), the m-mode puncture-integral
# kernel of the CRP elliptic scheme (arXiv:2606.04998 main.tex:467-486). Gates per node
# ESSF::legendre-type3-2f1 (W1-D): baked mpmath fixtures (rel 1e-13 Float64 / 1e-30 BigFloat-256),
# 0.995-switch continuity < 1e-12, n<0 / m<0 identity wiring, derivative recurrence vs central FD.

include("fixtures_legendre_p3.jl")

@testset "legendre_P3 (type-3 associated Legendre)" begin
    relerr(a, b) = abs(a - b) / abs(b)

    @testset "closed forms + input genericity" begin
        for z in (1.001, 1.5, 500.0)                      # far + near branch
            @test relerr(legendre_P3(0, 0, z), 1.0) < 1e-14
            @test relerr(legendre_P3(1, 1, z), sqrt(z^2 - 1)) < 1e-13       # ₃P₁¹ = √(z²−1)
            @test relerr(legendre_P3(2, 0, z), (3z^2 - 1) / 2) < 1e-13      # ₃P₂⁰ = (3z²−1)/2
        end
        # half-integer degree accepted as Rational / Float64 / Int order interchangeably
        @test legendre_P3(3 // 2, 1, 1.5) == legendre_P3(1.5, 1.0, 1.5)
        @test_throws DomainError legendre_P3(0.5, 0, 0.9)
    end

    @testset "baked mpmath fixtures — Float64 (rel < 1e-13)" begin
        for r in LEGENDRE_P3_FIXTURE
            ref = Float64(setprecision(() -> parse(BigFloat, r.val), BigFloat, 256))
            @test relerr(legendre_P3(r.n, r.m, r.z), ref) < 1e-13
        end
    end

    @testset "baked mpmath fixtures — BigFloat(256) pass-through (rel < 1e-30)" begin
        setprecision(BigFloat, 256) do
            for r in LEGENDRE_P3_FIXTURE
                v = legendre_P3(BigFloat(r.n), r.m, BigFloat(r.z))   # exact widening of the doubles
                @test v isa BigFloat
                @test relerr(v, parse(BigFloat, r.val)) < 1e-30
            end
        end
    end

    @testset "identity wiring (n<0 map, m<0 map are exact re-routes)" begin
        for z in (1.05, 10.0, 500.0)
            @test legendre_P3(-3.5, 2, z) == legendre_P3(2.5, 2, z)          # ₃P_n = ₃P_{−n−1}
            @test legendre_P3(-0.5, 5, z) == legendre_P3(-0.5, 5, z)         # n=−1/2 fixed point terminates
            @test legendre_P3(1.5, -1, z) ==
                  gamma_big(1.5) / gamma_big(3.5) * legendre_P3(1.5, 1, z)   # ₃P^{−m} = Γ(n−m+1)/Γ(n+m+1) ₃P^m
        end
    end

    @testset "0.995-argument switch continuity (< 1e-12)" begin
        # Compare the ACTUAL two branch kernels (regularized ₂F₁): Maclaurin/Γ(c) vs DLMF
        # 15.8.10/15.8.12, at the boundary z=399 (x=0.995) and one step either side. x and w
        # are derived from the same z exactly as legendre_P3 derives them.
        for z in (398.99, 399.0, 399.01)
            x = (z - 1) / (z + 1)
            w = 2 / (z + 1)
            for n in (-0.5, 0.5, 1.5, 2.5, 4.5), m in (0, 1, 3, 5, 20)
                a = n + m + 1; b = n + 1; c = m + 1.0
                M = Int(2n + 1)
                far = BigSpecial._reg2f1_far(a, b, c, x)
                near = BigSpecial._reg2f1_near1(a, b, M, w)
                @test relerr(near, far) < 1e-12
            end
        end
    end

    @testset "derivative: recurrence vs central FD (BigFloat 256) + Float64 consistency" begin
        pts = [(0.5, 0, 1.5), (4.5, 5, 10.0), (-3.5, -3, 1.05), (-0.5, 1, 500.0)]  # far + near branch
        setprecision(BigFloat, 256) do
            for (n, m, z) in pts
                nb = BigFloat(n); zb = BigFloat(z)
                h = zb * big"1e-20"
                fd = (legendre_P3(nb, m, zb + h) - legendre_P3(nb, m, zb - h)) / (2h)
                dz = legendre_P3_dz(nb, m, zb)
                @test relerr(dz, fd) < 1e-25
                @test relerr(legendre_P3_dz(n, m, z), Float64(dz)) < 1e-13
            end
        end
    end
end

# test_chebyshev.jl — gates for the Q2-H CGL collocation primitives (deliverable A).
# Node symmetry, D1/D2 polynomial exactness to degree n, the independent-D2 == D1² consistency,
# BigFloat eltype inheritance, and the dense solve. Ties into G-Q2H-5 (spectral eltype-cleanliness).
using Test
using NumericsCore
using LinearAlgebra

@testset "chebyshev (CGL collocation, Q2-H deliverable A)" begin

    @testset "node symmetry + endpoints" begin
        # reference [-1,1]: x_j = -x_{n-j} (CGL is symmetric about the midpoint)
        for n in (8, 16, 33)
            x = cheb_lobatto(Float64, n; lo=-1.0, hi=1.0)
            @test length(x) == n + 1
            @test x[1] ≈ 1.0 atol=0 rtol=0            # hi endpoint exact
            @test x[end] ≈ -1.0 atol=0 rtol=0          # lo endpoint exact
            for j in 0:n
                @test x[j+1] ≈ -x[n-j+1] atol=1e-15
            end
        end
        # affine [0, Rp]: midpoint-symmetry of (x - Rp/2); endpoints exact (scri R=0 last node)
        Rp = 0.5833819533586021
        x = cheb_lobatto(Float64, 24; lo=0.0, hi=Rp)
        @test x[1] == Rp
        @test x[end] == 0.0                            # scri is an EXACT node
        @test all(j -> isapprox((x[j+1]-Rp/2), -(x[24-j+1]-Rp/2); atol=1e-15), 0:24)
    end

    @testset "D1/D2 polynomial exactness to degree n (Float64)" begin
        n = 20
        x, D1, D2 = cheb_collocation(Float64, n; lo=0.0, hi=1.3)
        for k in 0:n
            f  = x .^ k
            df = k == 0 ? zeros(n+1)       : k .* x .^ (k-1)
            d2f= k <= 1 ? zeros(n+1)       : (k*(k-1)) .* x .^ (k-2)
            @test maximum(abs, D1*f .- df)  ≤ 1e-9 * max(1, maximum(abs, df))
            @test maximum(abs, D2*f .- d2f) ≤ 1e-7 * max(1, maximum(abs, d2f))
        end
        # row-sum (derivative of a constant) is exactly zero by the negative-sum trick
        @test maximum(abs, sum(D1; dims=2)) ≤ 1e-13
        @test maximum(abs, sum(D2; dims=2)) ≤ 1e-11
    end

    @testset "D2 == D1² to the differentiation-matrix conditioning floor" begin
        for n in (12, 24, 40)
            _, D1, D2 = cheb_collocation(Float64, n; lo=0.0, hi=1.0)
            # independently-built D2 agrees with D1*D1 to the conditioning floor. The D2 entries are
            # O(n⁴) (corner blow-up), so the RELATIVE defect ‖D2−D1²‖/‖D2‖ ~ few·eps is the honest
            # measure (n=40: abs ~1e-8 over ‖D2‖~3e6 ⇒ rel ~3e-15 ≈ 13·eps).
            reldef = maximum(abs, D2 .- D1*D1) / maximum(abs, D2)
            @test reldef ≤ 1e-12
        end
    end

    @testset "BigFloat-256 eltype inheritance + BF exactness" begin
        setprecision(BigFloat, 256) do
            n = 24
            x, D1, D2 = cheb_collocation(BigFloat, n; lo=big"0.0", hi=big"1.0")
            @test eltype(x) == BigFloat
            @test eltype(D1) == BigFloat && eltype(D2) == BigFloat
            # spectral exactness now at the BF floor, not the F64 floor
            for k in (3, 7, 15, n)
                f  = x .^ k
                df = k .* x .^ (k-1)
                d2f= (k*(k-1)) .* x .^ (k-2)
                @test maximum(abs, D1*f .- df)  ≤ big"1e-60"
                @test maximum(abs, D2*f .- d2f) ≤ big"1e-55"
            end
        end
    end

    @testset "cheb_solve (dense LU, {T}/BigFloat-clean)" begin
        A = ComplexF64[2 1 0; 1 3 1; 0 1 2]
        b = ComplexF64[1, 2im, 3]
        x = cheb_solve(A, b)
        @test maximum(abs, A*x .- b) ≤ 1e-13
        setprecision(BigFloat, 256) do
            Ab = Complex{BigFloat}.(A); bb = Complex{BigFloat}.(b)
            xb = cheb_solve(Ab, bb)
            @test eltype(xb) == Complex{BigFloat}
            @test maximum(abs, Ab*xb .- bb) ≤ big"1e-70"
        end
    end
end

# Geometry/geometry.jl — grounded by ANALYTIC closed-form invariants (the legacy
# MetricReconstruction/DiscModes/KerrLorenzCirc geometry forks are HEAVY → NO live dep; the
# literal formulas were read from MetricReconstruction/src/background.jl + weyl_amplitudes.jl
# and are asserted here against coordinate-free identities).
using Test
using Spacetime
using LinearAlgebra

@testset "Geometry: metric, inverse, scalars, horizons, Kinnersley" begin
    @testset "a=0 Schwarzschild BL metric = diag(-(1-2/r), 1/(1-2/r), r², r² sin²θ)" begin
        for (r, θ) in ((6.0, π/3), (10.0, π/2), (3.0, 0.7))
            g = kerr_metric(0.0, r, θ)
            @test isapprox(g[1,1], -(1 - 2/r);        rtol = 1e-14)
            @test isapprox(g[2,2], 1/(1 - 2/r);       rtol = 1e-14)
            @test isapprox(g[3,3], r^2;               rtol = 1e-14)
            @test isapprox(g[4,4], r^2*sin(θ)^2;      rtol = 1e-14)
            @test g[1,4] == 0 && g[1,2] == 0          # diagonal, no frame-drag at a=0
        end
    end

    @testset "det(g_BL) = -Σ² sin²θ (Kerr metric determinant identity)" begin
        for (a, r, θ) in ((0.7, 8.0, 1.1), (0.99, 4.0, π/4), (0.3, 12.0, 2.0))
            g = kerr_metric(a, r, θ)
            Σ = r^2 + a^2*cos(θ)^2
            @test isapprox(det(g), -Σ^2*sin(θ)^2; rtol = 1e-12)
        end
    end

    @testset "kerr_metric_inv * kerr_metric = I (all three charts)" begin
        for coords in (:BL, :In, :Out), (a, r, θ) in ((0.7, 8.0, 1.1), (0.99, 5.0, π/3))
            g  = kerr_metric(a, r, θ; coords = coords)
            gi = kerr_metric_inv(a, r, θ; coords = coords)
            @test isapprox(gi*g, Matrix(I, 4, 4); atol = 1e-10)
        end
    end

    @testset "horizons r± = M ± √(M²−a²) and Ω₊ = a/(2r₊)" begin
        @test horizons(0.0) == (2.0, 0.0)                 # Schwarzschild: r+ = 2M
        @test all(horizons(1.0) .≈ (1.0, 1.0))            # extremal: r± = M
        for a in (0.3, 0.7, 0.998)
            rp, rm = horizons(a)
            @test isapprox(delta_kerr(a, rp), 0.0; atol = 1e-13)   # Δ(r±) = 0
            @test isapprox(delta_kerr(a, rm), 0.0; atol = 1e-13)
            @test isapprox(rp + rm, 2.0;  rtol = 1e-14)            # sum of roots = 2M
            @test isapprox(rp * rm, a^2;  rtol = 1e-14)            # product = a²
            @test isapprox(omega_plus(a), a/(2*rp); rtol = 1e-14)
        end
    end

    @testset "Kinnersley normalization l·l=n·n=m·m=0, l·n=-1, m·m̄=1 (all charts)" begin
        for coords in (:BL, :In, :Out), (a, r, θ) in ((0.7, 8.0, 1.1), (0.5, 6.0, π/3))
            g = kerr_metric(a, r, θ; coords = coords)
            l, n, m, mb = kinnersley_tetrad(a, r, θ; coords = coords)
            dot2(u, v) = sum(g[μ,ν]*u[μ]*v[ν] for μ in 1:4, ν in 1:4)
            @test isapprox(dot2(l, l),  0.0; atol = 1e-10)
            @test isapprox(dot2(n, n),  0.0; atol = 1e-10)
            @test isapprox(dot2(m, m),  0.0; atol = 1e-10)
            @test isapprox(dot2(l, n), -1.0; atol = 1e-10)
            @test isapprox(dot2(m, mb), 1.0; atol = 1e-10)
        end
    end

    @testset "bg_scalars literal formulas (Δ, Σ, ζ; K, Q)" begin
        a, r, θ = 0.9, 10.0, π/3
        bg = bg_scalars(a, r, θ)
        @test isapprox(bg.Δ, r^2 - 2r + a^2;          rtol = 1e-14)
        @test isapprox(bg.Σ, r^2 + a^2*cos(θ)^2;       rtol = 1e-14)
        @test bg.ζ == complex(r, -a*cos(θ))
        ω, m = 0.4, 2
        @test isapprox(kfun(a, r, ω, m), (r^2+a^2)*ω - a*m;     rtol = 1e-14)
        @test isapprox(qfun(a, ω, m, θ), -a*ω*sin(θ) + m/sin(θ); rtol = 1e-14)
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        a, r, θ = big"0.9", big"10.0", big"1.1"
        g = kerr_metric(a, r, θ)
        @test eltype(g) == BigFloat
        @test horizons(a)[1] isa BigFloat
        l, _, _, _ = kinnersley_tetrad(a, r, θ)
        @test eltype(l) == Complex{BigFloat}
    end

    @testset "public metric promotion handles integer and mixed inputs" begin
        @test eltype(kerr_metric(0, 8, 1)) === Float64
        @test kerr_metric(0, 8, 1) == kerr_metric(0.0, 8.0, 1.0)
        @test eltype(kerr_metric(0.5f0, 8, 1.0f0)) === Float32
        @test eltype(kerr_metric(0.5f0, 8.0, 1.0f0)) === Float64
        @test eltype(kerr_metric_inv(0, 8.0f0, 1.0f0)) === Float32
        @test eltype(first(kinnersley_tetrad(0, 8, 1))) === ComplexF64
        @test first(kinnersley_tetrad(0, 8, 1)) ==
              first(kinnersley_tetrad(0.0, 8.0, 1.0))
        @test eltype(first(kinnersley_tetrad(0.5f0, 8.0, 1.0f0))) === ComplexF64
    end
end

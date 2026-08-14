# Geometry/christoffel.jl — grounded ANALYTIC (legacy MetricReconstruction/Punctures are HEAVY →
# NO live dep): the M-only Schwarzschild closed forms are baked verbatim from
# Punctures/worldline_self_field.jl (lines 93-130); the Kerr Γ is checked against a
# finite-difference of geometry.jl's kerr_metric (½g^{λσ}(∂g+∂g−∂g) self-consistency); the
# Schwarzschild-Γ r/θ derivatives are checked against finite differences of `christoffel`.
using Test
using Spacetime

# Central finite difference of a 4×4×4 Γ-producing closure.
function _fd(f, x; h = 1e-6)
    (f(x + h) .- f(x - h)) ./ (2h)
end

@testset "Christoffel symbols (Kerr BL + Schwarzschild M-only)" begin
    @testset "Schwarzschild M-only literals (baked from Punctures source)" begin
        for (r, θ, M) in ((8.0, 1.1, 1.0), (6.0, π/3, 2.0), (10.0, 0.7, 1.5))
            Γ = christoffel(r, θ; M = M)
            f = 1 - 2M/r; s = sin(θ); c = cos(θ)
            @test Γ[1,1,2] ≈ M/(r^2*f)        # Γ^t_{tr} = M/(r(r-2M))
            @test Γ[2,1,1] ≈ M*f/r^2          # Γ^r_{tt} = M(r-2M)/r³
            @test Γ[2,2,2] ≈ -M/(r^2*f)       # Γ^r_{rr}
            @test Γ[2,3,3] ≈ -(r - 2M)        # Γ^r_{θθ}
            @test Γ[2,4,4] ≈ -(r - 2M)*s^2    # Γ^r_{φφ}
            @test Γ[3,2,3] ≈ 1/r              # Γ^θ_{rθ}
            @test Γ[3,4,4] ≈ -s*c             # Γ^θ_{φφ}
            @test Γ[4,2,4] ≈ 1/r              # Γ^φ_{rφ}
            @test Γ[4,3,4] ≈ c/s              # Γ^φ_{θφ}
        end
    end

    @testset "lower-index symmetry Γ^λ_{μν} = Γ^λ_{νμ}" begin
        for (a, r, θ) in ((0.7, 8.0, 1.1), (0.99, 5.0, π/3))
            Γ = christoffel(a, r, θ)
            for λ in 1:4, μ in 1:4, ν in 1:4
                @test Γ[λ,μ,ν] == Γ[λ,ν,μ]
            end
        end
        Γs = christoffel(8.0, 1.1; M = 1.0)
        for λ in 1:4, μ in 1:4, ν in 1:4
            @test Γs[λ,μ,ν] == Γs[λ,ν,μ]
        end
    end

    @testset "Kerr Γ = ½ g^{λσ}(∂_μ g_σν + ∂_ν g_σμ − ∂_σ g_μν) via FD of kerr_metric" begin
        for (a, r, θ) in ((0.7, 8.0, 1.1), (0.5, 12.0, π/4), (0.9, 6.0, 2.0))
            Γ = christoffel(a, r, θ)
            h = 1e-6
            gr = (kerr_metric(a, r+h, θ) .- kerr_metric(a, r-h, θ)) ./ (2h)
            gθ = (kerr_metric(a, r, θ+h) .- kerr_metric(a, r, θ-h)) ./ (2h)
            gi = kerr_metric_inv(a, r, θ)
            dg = zeros(4,4,4); dg[2,:,:] = gr; dg[3,:,:] = gθ
            Γfd = zeros(4,4,4)
            for λ in 1:4, μ in 1:4, ν in 1:4
                Γfd[λ,μ,ν] = 0.5*sum(gi[λ,σ]*(dg[μ,σ,ν] + dg[ν,σ,μ] - dg[σ,μ,ν]) for σ in 1:4)
            end
            @test isapprox(Γ, Γfd; atol = 1e-5)
        end
    end

    @testset "a→0 Kerr reduces to Schwarzschild (M=1) connection" begin
        for (r, θ) in ((8.0, 1.1), (6.0, π/3), (12.0, 0.7))
            @test isapprox(christoffel(0.0, r, θ), christoffel(r, θ; M = 1.0); atol = 1e-12)
        end
    end

    @testset "christoffel_dr / christoffel_dtheta == FD of christoffel" begin
        for (r, θ, M) in ((8.0, 1.1, 1.0), (6.0, π/3, 1.5))
            Γr_fd = _fd(ρ -> christoffel(ρ, θ; M = M), r)
            Γθ_fd = _fd(ϑ -> christoffel(r, ϑ; M = M), θ)
            @test isapprox(christoffel_dr(r, θ; M = M),     Γr_fd; atol = 1e-5)
            @test isapprox(christoffel_dtheta(r, θ; M = M), Γθ_fd; atol = 1e-5)
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        Γ = christoffel(big"0.9", big"8.0", big"1.1")
        @test eltype(Γ) == BigFloat
        Γs = christoffel(big"8.0", big"1.1"; M = big"1.0")
        @test eltype(Γs) == BigFloat
    end


    @testset "full-Kerr public promotion handles integer and mixed inputs" begin
        @test eltype(christoffel(0, 8, 1)) === Float64
        @test christoffel(0, 8, 1) == christoffel(0.0, 8.0, 1.0)
        @test eltype(christoffel(0.5f0, 8, 1.0f0)) === Float32
        @test eltype(christoffel(0.5f0, 8.0, 1.0f0)) === Float64
    end
end

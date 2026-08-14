using Test
using MetricReconstruction.GHZ

const GHZ_FP = UInt128(0x7d62a643d19e4c6b9ab21750a8829e31)

mode(::Val{N}, S; T=Float64, l=2, m=1, omega=0.13 - 0.07im,
     fp=GHZ_FP, normalization=EinsteinRHS()) where {N} =
    TetradSourceMode{N}(S; M=T(1), l, m, omega=Complex{T}(omega),
                        source_fingerprint=fp, normalization)

@testset "GHZ — outgoing-Schwarzschild IRG corrector head" begin
    @test GHZIRG !== HertzIRG

    @testset "zero source and typed order" begin
        for N in (1, 2)
            src = mode(Val(N), r -> 0.0)
            bc = InnerRayIVP(src; rb=2.0, x_mmb=2 - 3im, dr_x_mmb=0.25 + 0.5im)
            sol = solve_irg_corrector_head(src, bc)
            @test sol isa IRGCorrectorHeadMode
            @test perturbative_order(src) == N
            @test perturbative_order(bc) == N
            @test perturbative_order(sol) == N
            @test sol.gauge isa GHZIRG
            for r in (2.0, 3.5, 8.0)
                expected = bc.x0 + bc.rb^2 * bc.dx0 * (1 / bc.rb - 1 / r)
                @test x_mmb(sol, r) ≈ expected atol=2e-14
                @test dr_x_mmb(sol, r) ≈ bc.rb^2 * bc.dx0 / r^2 atol=2e-14
            end
        end
    end

    @testset "constant complex source fixes the overall sign" begin
        C = 0.7 - 0.2im
        rb = 2.5
        src = mode(Val(2), r -> C)
        sol = solve_irg_corrector_head(src, InnerRayIVP(src; rb))
        for r in (rb, 3.0, 7.25)
            expected_x = -C * ((r^2 - rb^2) / 2 - (r^3 - rb^3) / (3r))
            expected_dx = -C * (r^3 - rb^3) / (3r^2)
            expected_d2x = -C - 2expected_dx / r
            @test x_mmb(sol, r) ≈ expected_x rtol=2e-13 atol=2e-14
            @test dr_x_mmb(sol, r) ≈ expected_dx rtol=2e-13 atol=2e-14
            @test irg_head_residual(r, expected_dx, expected_d2x, C) ≈ 0 atol=2e-15
            # The opposite overall sign is decisively nonzero away from the boundary.
            r > rb && @test abs(irg_head_residual(r, -expected_dx, -expected_d2x, C)) > 1
        end
    end

    @testset "linearity" begin
        A(r) = (0.2 + 0.1im) * r
        B(r) = (-0.4 + 0.3im) / r
        a, b = 1.7 - 0.2im, -0.6 + 0.8im
        srcA = mode(Val(2), A; fp=GHZ_FP + 1)
        srcB = mode(Val(2), B; fp=GHZ_FP + 2)
        srcAB = mode(Val(2), r -> a * A(r) + b * B(r); fp=GHZ_FP + 3)
        solA = solve_irg_corrector_head(srcA, InnerRayIVP(srcA; rb=2.2))
        solB = solve_irg_corrector_head(srcB, InnerRayIVP(srcB; rb=2.2))
        solAB = solve_irg_corrector_head(srcAB, InnerRayIVP(srcAB; rb=2.2))
        for r in (2.2, 3.1, 5.6)
            @test x_mmb(solAB, r) ≈ a * x_mmb(solA, r) + b * x_mmb(solB, r) rtol=3e-13 atol=2e-14
            @test dr_x_mmb(solAB, r) ≈ a * dr_x_mmb(solA, r) + b * dr_x_mmb(solB, r) rtol=3e-13 atol=2e-14
        end
    end

    @testset "BigFloat precision is preserved" begin
        setprecision(BigFloat, 256) do
            C = Complex{BigFloat}(big"0.125", big"-0.375")
            src = mode(Val(2), r -> C; T=BigFloat, omega=Complex{BigFloat}(big"0.03", big"0.01"))
            rb, r = big"2.0", big"4.75"
            sol = solve_irg_corrector_head(src, InnerRayIVP(src; rb))
            got = x_mmb(sol, r; rtol=big"1e-60")
            expected = -C * ((r^2 - rb^2) / 2 - (r^3 - rb^3) / (3r))
            @test got isa Complex{BigFloat}
            @test dr_x_mmb(sol, r; rtol=big"1e-60") isa Complex{BigFloat}
            @test abs(got - expected) < big"1e-55"
        end
    end

    @testset "normalization, metadata, and boundary rejection" begin
        @test_throws ArgumentError TetradSourceMode{0}(r -> 0.0; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_FP)
        @test_throws ArgumentError mode(Val(2), r -> 0.0; l=-1)
        @test_throws ArgumentError mode(Val(2), r -> 0.0; l=2, m=3)
        @test_throws ArgumentError TetradSourceMode{2}(r -> 0.0; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_FP, tetrad=:wrong)
        @test_throws ArgumentError TetradSourceMode{2}(r -> 0.0; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_FP, slicing=:wrong)

        physical = mode(Val(2), r -> 1.0; normalization=PhysicalStressEnergy())
        @test_throws ArgumentError solve_irg_corrector_head(
            physical, InnerRayIVP(physical; rb=2.0))

        src = mode(Val(2), r -> 1.0)
        @test_throws ArgumentError InnerRayIVP(src; rb=1.99)
        @test_throws ArgumentError InnerRayIVP(src; rb=2.0, x_mmb=Inf)
        sol = solve_irg_corrector_head(src, InnerRayIVP(src; rb=2.0))
        @test_throws ArgumentError x_mmb(sol, 1.99)
        @test_throws ArgumentError dr_x_mmb(sol, NaN)

        # Boundary copied from a different mode/fingerprint cannot be reused.
        other_l = mode(Val(2), r -> 1.0; l=3, m=1, fp=GHZ_FP + 9)
        @test_throws ArgumentError solve_irg_corrector_head(src, InnerRayIVP(other_l; rb=2.0))
        other_N = mode(Val(1), r -> 1.0; fp=GHZ_FP)
        @test_throws ArgumentError solve_irg_corrector_head(src, InnerRayIVP(other_N; rb=2.0))
    end

    @testset "checked harmonic mode indices" begin
        indexed_source(l, m) = TetradSourceMode{2}(r -> 0.0; M=1.0, l, m,
            omega=0.0, source_fingerprint=GHZ_FP + 60)

        @test_throws ArgumentError indexed_source(true, 0)
        @test_throws ArgumentError indexed_source(2, true)
        @test_throws ArgumentError indexed_source(2, typemin(Int))
        @test_throws ArgumentError indexed_source(typemax(Int), typemin(Int))

        too_big = BigInt(typemax(Int)) + 1
        too_uint = UInt(typemax(Int)) + one(UInt)
        for invalid in (too_big, too_uint)
            @test_throws ArgumentError indexed_source(invalid, 0)
            @test_throws ArgumentError indexed_source(2, invalid)
        end

        for (l, m, expected) in ((BigInt(2), BigInt(-2), (2, -2)),
                (UInt(2), UInt(2), (2, 2)), (Int128(2), Int128(-2), (2, -2)))
            src = indexed_source(l, m)
            @test (src.l, src.m) == expected
        end

        l0 = indexed_source(BigInt(0), UInt(0))
        @test (l0.l, l0.m) == (0, 0)
        for m in (-typemax(Int), typemax(Int))
            extreme = indexed_source(typemax(Int), m)
            @test (extreme.l, extreme.m) == (typemax(Int), m)
        end
    end

    @testset "Boolean scalar rejection" begin
        radial = r -> 0.0
        @test_throws ArgumentError TetradSourceMode{2}(radial; M=true, l=2, m=0,
            omega=0.0, source_fingerprint=GHZ_FP)
        @test_throws ArgumentError TetradSourceMode{2}(radial; M=1.0, l=2, m=0,
            omega=true, source_fingerprint=GHZ_FP)
        @test_throws ArgumentError TetradSourceMode{2}(radial; M=1.0, l=2, m=0,
            omega=Complex{Bool}(true, false), source_fingerprint=GHZ_FP)
        @test_throws ArgumentError TetradSourceMode{2}(radial; M=1.0, l=2, m=0,
            omega=0.0, source_fingerprint=true)

        low_mass_src = TetradSourceMode{2}(radial; M=0.25, l=2, m=0, omega=0.0,
            source_fingerprint=GHZ_FP + 42)
        for bad_radius in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError InnerRayIVP(low_mass_src; rb=bad_radius)
        end
        low_mass_sol = solve_irg_corrector_head(
            low_mass_src, InnerRayIVP(low_mass_src; rb=0.5))
        for bad_radius in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError x_mmb(low_mass_sol, bad_radius)
            @test_throws ArgumentError dr_x_mmb(low_mass_sol, bad_radius)
        end

        src = mode(Val(2), radial)
        @test_throws ArgumentError InnerRayIVP(src; rb=2.0, x_mmb=true)
        @test_throws ArgumentError InnerRayIVP(src; rb=2.0, dr_x_mmb=true)
        @test_throws ArgumentError InnerRayIVP(src; rb=2.0,
            x_mmb=Complex{Bool}(true, false))

        boundary_bool = mode(Val(2), r -> true; fp=GHZ_FP + 40)
        @test_throws ArgumentError solve_irg_corrector_head(
            boundary_bool, InnerRayIVP(boundary_bool; rb=2.0))
        late_bool = mode(Val(2), r -> r > 2.1 ? true : 0.0; fp=GHZ_FP + 41)
        late_sol = solve_irg_corrector_head(late_bool, InnerRayIVP(late_bool; rb=2.0))
        @test_throws ArgumentError x_mmb(late_sol, 3.0)

        sol = solve_irg_corrector_head(src, InnerRayIVP(src; rb=2.0))
        @test_throws ArgumentError x_mmb(sol, 2.0; atol=true)
        @test_throws ArgumentError dr_x_mmb(sol, 2.0; rtol=true)
    end

    @testset "residual scalar validation and result types" begin
        valid = (2.5, 0.25 - 0.5im, -0.1 + 0.2im, 0.3 + 0.4im)
        expected = -(valid[3] + 2valid[2] / valid[1]) - valid[4]
        got = irg_head_residual(valid...)
        @test got isa ComplexF64
        @test got == expected

        for slot in eachindex(valid), invalid in (true, Complex{Bool}(true, false))
            args = Base.setindex(valid, invalid, slot)
            @test_throws ArgumentError irg_head_residual(args...)
        end
        for slot in eachindex(valid), invalid in (Inf, NaN, complex(1.0, Inf))
            args = Base.setindex(valid, invalid, slot)
            @test_throws ArgumentError irg_head_residual(args...)
        end

        f16 = irg_head_residual(Float16(3), Float16(0.2), Float16(-0.1),
                                Float16(0.3))
        @test f16 isa Float16
        @test isfinite(f16)
        setprecision(BigFloat, 256) do
            args = (big"3.0", Complex{BigFloat}(big"0.2", big"-0.1"),
                    Complex{BigFloat}(big"-0.3", big"0.4"),
                    Complex{BigFloat}(big"0.5", big"0.6"))
            value = irg_head_residual(args...)
            @test value isa Complex{BigFloat}
            @test value == -(args[3] + 2args[2] / args[1]) - args[4]
        end
    end

    @testset "scale-safe finite inputs and tolerance rejection" begin
        M = 9e153
        rb = 2M
        src = TetradSourceMode{2}(r -> 0.0; M, l=0, m=0,
                                  source_fingerprint=GHZ_FP + 20)
        bc = InnerRayIVP(src; rb, x_mmb=2.0, dr_x_mmb=1e-154)
        sol = solve_irg_corrector_head(src, bc)
        @test x_mmb(sol, rb) == 2.0 + 0im
        @test dr_x_mmb(sol, rb) == 1e-154 + 0im

        R = 1e155
        C = inv(R) * inv(R)
        scaled = mode(Val(2), r -> C; l=0, m=0, fp=GHZ_FP + 21)
        scaled_sol = solve_irg_corrector_head(scaled, InnerRayIVP(scaled; rb=2.0))
        @test isapprox(real(x_mmb(scaled_sol, R)), -1 / 6; rtol=2e-12)
        @test isapprox(real(dr_x_mmb(scaled_sol, R)), -inv(R) / 3; rtol=2e-12)

        for kwargs in ((; atol=-1.0), (; rtol=-1.0), (; atol=NaN), (; rtol=Inf),
                       (; atol=0.0, rtol=0.0))
            @test_throws ArgumentError x_mmb(sol, rb; kwargs...)
            @test_throws ArgumentError dr_x_mmb(sol, rb; kwargs...)
        end
    end
end

using Test
using MetricReconstruction.GHZ

const GHZ_MIDDLE_FP = UInt128(0x81a994e4293548549ac6d6057b229d11)

function middle_mode(::Val{N}, Sll; Slm=nothing, T=Float64, l=2, m=1,
        omega=0.11 - 0.04im, fp=GHZ_MIDDLE_FP) where {N}
    common = (; M=T(1), l, m, omega=Complex{T}(omega), source_fingerprint=fp)
    isnothing(Slm) ? TetradSourceMode{N}(Sll; common...) :
                     TetradSourceMode{N}(Sll; Slm, common...)
end

@testset "GHZ — outgoing-Schwarzschild IRG corrector middle" begin
    @testset "zero force, homogeneous data, and typed order" begin
        for N in (1, 2)
            src = middle_mode(Val(N), r -> zero(r))
            head = solve_irg_corrector_head(src, InnerRayIVP(src; rb=2.0))
            bc = MiddleRayIVP(head; x_nm=1.2 - 0.7im, dr_x_nm=-0.3 + 0.2im)
            sol = solve_irg_corrector_middle(head, bc)
            @test sol isa IRGCorrectorMiddleMode
            @test perturbative_order(bc) == N
            @test perturbative_order(sol) == N
            @test sol.gauge isa GHZIRG
            for r in (2.0, 3.25, 7.0)
                A = bc.dx0 + 2bc.x0 / bc.rb
                z = bc.rb^2 * bc.x0 + A * (r^3 - bc.rb^3) / 3
                @test x_nm(sol, r) ≈ z / r^2 atol=2e-14
                @test dr_x_nm(sol, r) ≈ A - 2z / r^3 atol=2e-14
            end
        end
        src = middle_mode(Val(2), r -> zero(r); T=BigFloat)
        z = src.Slm(big"3.25")
        @test z == 0
        @test z isa BigFloat
    end

    @testset "constant complex S_lm fixes the sign" begin
        C, rb = 0.35 - 0.12im, 2.4
        src = middle_mode(Val(2), r -> 0.0; Slm=r -> C)
        head = solve_irg_corrector_head(src, InnerRayIVP(src; rb))
        sol = solve_irg_corrector_middle(head, MiddleRayIVP(head))
        F = -2C
        for r in (rb, 3.1, 6.75)
            z = F / 3 * (r^3 * (r - rb) - (r^4 - rb^4) / 4)
            expected_x = z / r^2
            expected_dx = F / 3 * ((r - rb) + (r^4 - rb^4) / (2r^3))
            expected_d2x = F - 2expected_dx / r + 2expected_x / r^2
            @test x_nm(sol, r) ≈ expected_x rtol=3e-13 atol=2e-14
            @test dr_x_nm(sol, r) ≈ expected_dx rtol=3e-13 atol=2e-14
            @test irg_middle_residual(r, expected_x, expected_dx, expected_d2x,
                                      C, 0.0, src.l) ≈ 0 atol=3e-15
            r > rb && @test abs(irg_middle_residual(
                r, -expected_x, -expected_dx, -expected_d2x, C, 0.0, src.l)) > 0.5
        end
    end

    @testset "manufactured head coupling and linearity" begin
        rb, l = 2.25, 2
        D = 0.4 - 0.15im
        src = middle_mode(Val(2), r -> 0.0; l)
        head = solve_irg_corrector_head(
            src, InnerRayIVP(src; rb, dr_x_mmb=D))
        sol = solve_irg_corrector_middle(head, MiddleRayIVP(head))
        K = sqrt(l * (l + 1) / 2) * D * rb^2
        for r in (rb, 3.0, 5.5)
            expected_x = K / 3 * (r / (2rb^2) - 3 / (2r) + rb / r^2)
            expected_dx = K / 3 * (1 / (2rb^2) + 3 / (2r^2) - 2rb / r^3)
            expected_d2x = K / 3 * (-3 / r^3 + 6rb / r^4)
            dxhead = D * rb^2 / r^2
            @test x_nm(sol, r) ≈ expected_x rtol=4e-13 atol=2e-14
            @test dr_x_nm(sol, r) ≈ expected_dx rtol=4e-13 atol=2e-14
            @test irg_middle_residual(r, expected_x, expected_dx, expected_d2x,
                                      0.0, dxhead, l) ≈ 0 atol=4e-15
        end

        LA(r) = (0.1 + 0.2im) * r
        LB(r) = (-0.3 + 0.1im) / r
        a, b = 1.3 - 0.2im, -0.4 + 0.6im
        srcA = middle_mode(Val(2), r -> 0.0; Slm=LA, fp=GHZ_MIDDLE_FP + 1)
        srcB = middle_mode(Val(2), r -> 0.0; Slm=LB, fp=GHZ_MIDDLE_FP + 2)
        srcAB = middle_mode(Val(2), r -> 0.0;
            Slm=r -> a * LA(r) + b * LB(r), fp=GHZ_MIDDLE_FP + 3)
        make_sol(src, d) = begin
            h = solve_irg_corrector_head(src, InnerRayIVP(src; rb=2.2, dr_x_mmb=d))
            solve_irg_corrector_middle(h, MiddleRayIVP(h))
        end
        solA, solB = make_sol(srcA, 0.2 - 0.1im), make_sol(srcB, -0.15 + 0.3im)
        solAB = make_sol(srcAB, a * (0.2 - 0.1im) + b * (-0.15 + 0.3im))
        for r in (2.2, 3.4, 5.0)
            @test x_nm(solAB, r) ≈ a * x_nm(solA, r) + b * x_nm(solB, r) rtol=5e-12 atol=4e-14
            @test dr_x_nm(solAB, r) ≈ a * dr_x_nm(solA, r) + b * dr_x_nm(solB, r) rtol=5e-12 atol=4e-14
        end
    end

    @testset "BigFloat precision" begin
        setprecision(BigFloat, 256) do
            C = Complex{BigFloat}(big"0.0625", big"-0.1875")
            rb, r = big"2.0", big"4.25"
            src = middle_mode(Val(2), q -> zero(q); Slm=q -> C, T=BigFloat,
                              omega=Complex{BigFloat}(big"0.02", big"-0.01"))
            head = solve_irg_corrector_head(src, InnerRayIVP(src; rb))
            sol = solve_irg_corrector_middle(head, MiddleRayIVP(head))
            F = -2C
            expected = F / 3 * (r^3 * (r - rb) - (r^4 - rb^4) / 4) / r^2
            got = x_nm(sol, r; rtol=big"1e-60")
            @test got isa Complex{BigFloat}
            @test dr_x_nm(sol, r; rtol=big"1e-60") isa Complex{BigFloat}
            @test abs(got - expected) < big"1e-55"
        end
    end

    @testset "middle spin domain excludes l=0 without narrowing the head" begin
        src0 = middle_mode(Val(2), r -> 0.0; l=0, m=0, fp=GHZ_MIDDLE_FP + 50)
        head0 = solve_irg_corrector_head(src0, InnerRayIVP(src0; rb=2.0))
        @test x_mmb(head0, 3.0) == 0.0 + 0.0im
        @test dr_x_mmb(head0, 3.0) == 0.0 + 0.0im
        @test_throws ArgumentError MiddleRayIVP(src0; rb=2.0)
        @test_throws ArgumentError MiddleRayIVP(head0)

        forged = MiddleRayIVP{2,Float64}(src0.M, src0.l, src0.m, src0.omega,
            src0.source_fingerprint, 2.0, 0.0 + 0.0im, 0.0 + 0.0im)
        @test_throws ArgumentError solve_irg_corrector_middle(head0, forged)
        raw_mode = IRGCorrectorMiddleMode(head0, forged, GHZIRG())
        @test_throws ArgumentError x_nm(raw_mode, forged.rb)
        @test_throws ArgumentError dr_x_nm(raw_mode, forged.rb)
        @test_throws ArgumentError irg_middle_residual(2.0, 0, 0, 0, 0, 0, 0)
        @test_throws ArgumentError GHZ._middle_coupling(Float64, 0)

        too_big = BigInt(typemax(Int)) + 1
        too_uint = UInt(typemax(Int)) + one(UInt)
        for invalid in (true, too_big, too_uint)
            @test_throws ArgumentError irg_middle_residual(2.0, 0, 0, 0, 0, 0, invalid)
            @test_throws ArgumentError GHZ._middle_coupling(Float64, invalid)
        end
        for valid in (BigInt(1), UInt(1), Int128(1))
            @test irg_middle_residual(2.0, 0, 0, 0, 0, 0, valid) == 0.0
            @test GHZ._middle_coupling(Float64, valid) ≈ 1.0
        end

        src1 = middle_mode(Val(2), r -> 0.0; l=1, m=0, fp=GHZ_MIDDLE_FP + 51)
        head1 = solve_irg_corrector_head(src1, InnerRayIVP(src1; rb=2.0))
        sol1 = solve_irg_corrector_middle(head1, MiddleRayIVP(head1))
        @test sol1 isa IRGCorrectorMiddleMode
        @test x_nm(sol1, 3.0) == 0.0 + 0.0im
        @test irg_middle_residual(2.0, 0, 0, 0, 0, 0, 1) == 0.0
    end

    @testset "source, boundary, and ray rejection" begin
        src = middle_mode(Val(2), r -> 0.0)
        head = solve_irg_corrector_head(src, InnerRayIVP(src; rb=2.0))
        bc = MiddleRayIVP(head)
        sol = solve_irg_corrector_middle(head, bc)
        @test sol.head.source === src
        @test (bc.M, bc.l, bc.m, bc.omega, bc.source_fingerprint) ==
              (src.M, src.l, src.m, src.omega, src.source_fingerprint)
        @test_throws ArgumentError MiddleRayIVP(src; rb=1.99)
        @test_throws ArgumentError MiddleRayIVP(src; rb=2.0, x_nm=Inf)
        @test_throws ArgumentError x_nm(sol, 1.9)
        @test_throws ArgumentError dr_x_nm(sol, NaN)
        @test_throws ArgumentError irg_middle_residual(0.0, 0, 0, 0, 0, 0, 2)
        @test_throws ArgumentError irg_middle_residual(2.0, 0, 0, 0, 0, 0, -1)

        different = middle_mode(Val(2), r -> 0.0; l=3, fp=GHZ_MIDDLE_FP + 9)
        @test_throws ArgumentError solve_irg_corrector_middle(
            head, MiddleRayIVP(different; rb=2.0))
        different_N = middle_mode(Val(1), r -> 0.0; fp=GHZ_MIDDLE_FP)
        @test_throws ArgumentError solve_irg_corrector_middle(
            head, MiddleRayIVP(different_N; rb=2.0))
        @test_throws ArgumentError solve_irg_corrector_middle(
            head, MiddleRayIVP(src; rb=2.5))

        bad = middle_mode(Val(2), r -> 0.0; Slm=r -> Inf, fp=GHZ_MIDDLE_FP + 20)
        bad_head = solve_irg_corrector_head(bad, InnerRayIVP(bad; rb=2.0))
        @test_throws ArgumentError solve_irg_corrector_middle(bad_head, MiddleRayIVP(bad_head))
        late = middle_mode(Val(2), r -> 0.0;
            Slm=r -> r > 2.1 ? NaN : 0.0, fp=GHZ_MIDDLE_FP + 21)
        late_head = solve_irg_corrector_head(late, InnerRayIVP(late; rb=2.0))
        late_sol = solve_irg_corrector_middle(late_head, MiddleRayIVP(late_head))
        @test_throws ArgumentError x_nm(late_sol, 3.0)
    end

    @testset "Boolean scalar rejection" begin
        low_mass_src = TetradSourceMode{2}(r -> 0.0; Slm=r -> 0.0, M=0.25,
            l=2, m=0, omega=0.0, source_fingerprint=GHZ_MIDDLE_FP + 43)
        low_mass_head = solve_irg_corrector_head(
            low_mass_src, InnerRayIVP(low_mass_src; rb=0.5))
        for bad_radius in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError MiddleRayIVP(low_mass_head; rb=bad_radius)
        end
        low_mass_sol = solve_irg_corrector_middle(
            low_mass_head, MiddleRayIVP(low_mass_head))
        for bad_radius in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError x_nm(low_mass_sol, bad_radius)
            @test_throws ArgumentError dr_x_nm(low_mass_sol, bad_radius)
        end

        src = middle_mode(Val(2), r -> 0.0; fp=GHZ_MIDDLE_FP + 40)
        head = solve_irg_corrector_head(src, InnerRayIVP(src; rb=2.0))
        @test_throws ArgumentError MiddleRayIVP(head; x_nm=true)
        @test_throws ArgumentError MiddleRayIVP(head; dr_x_nm=true)
        @test_throws ArgumentError MiddleRayIVP(head;
            x_nm=Complex{Bool}(true, false))

        boundary_bool = middle_mode(Val(2), r -> 0.0; Slm=r -> true,
            fp=GHZ_MIDDLE_FP + 41)
        boundary_head = solve_irg_corrector_head(
            boundary_bool, InnerRayIVP(boundary_bool; rb=2.0))
        @test_throws ArgumentError solve_irg_corrector_middle(
            boundary_head, MiddleRayIVP(boundary_head))

        late_bool = middle_mode(Val(2), r -> 0.0;
            Slm=r -> r > 2.1 ? true : 0.0, fp=GHZ_MIDDLE_FP + 42)
        late_head = solve_irg_corrector_head(late_bool, InnerRayIVP(late_bool; rb=2.0))
        late_sol = solve_irg_corrector_middle(late_head, MiddleRayIVP(late_head))
        @test_throws ArgumentError x_nm(late_sol, 3.0)

        sol = solve_irg_corrector_middle(head, MiddleRayIVP(head))
        @test_throws ArgumentError x_nm(sol, 2.0; atol=true)
        @test_throws ArgumentError dr_x_nm(sol, 2.0; rtol=true)
    end

    @testset "residual scalar validation and result types" begin
        valid = (2.5, 0.2 - 0.1im, 0.3 + 0.4im, -0.5 + 0.6im,
                 0.7 - 0.8im, -0.9 + 0.1im, 2)
        c = sqrt(3.0)
        expected = valid[4] + 2valid[3] / valid[1] - 2valid[2] / valid[1]^2 +
                   2valid[5] - c * valid[6] / valid[1]
        got = irg_middle_residual(valid...)
        @test got isa ComplexF64
        @test got == expected

        for slot in eachindex(valid), invalid in (true, Complex{Bool}(true, false))
            args = Base.setindex(valid, invalid, slot)
            @test_throws ArgumentError irg_middle_residual(args...)
        end
        for slot in eachindex(valid), invalid in (Inf, NaN, complex(1.0, Inf))
            args = Base.setindex(valid, invalid, slot)
            @test_throws ArgumentError irg_middle_residual(args...)
        end

        f16 = irg_middle_residual(Float16(3), Float16(0.2), Float16(0.1),
                                  Float16(-0.05), Float16(0.1), Float16(0), 256)
        @test f16 isa Float16
        @test isfinite(f16)
        setprecision(BigFloat, 256) do
            args = (big"3.0", Complex{BigFloat}(big"0.2", big"-0.1"),
                    Complex{BigFloat}(big"0.3", big"0.4"),
                    Complex{BigFloat}(big"-0.5", big"0.6"),
                    Complex{BigFloat}(big"0.7", big"-0.8"),
                    Complex{BigFloat}(big"-0.9", big"0.1"), 2)
            value = irg_middle_residual(args...)
            @test value isa Complex{BigFloat}
            @test value == args[4] + 2args[3] / args[1] - 2args[2] / args[1]^2 +
                           2args[5] - sqrt(big"3.0") * args[6] / args[1]
        end

        for callable in (irg_head_residual, irg_middle_residual)
            ms = collect(methods(callable))
            @test all(!Base.isambiguous(ms[i], ms[j])
                for i in eachindex(ms) for j in (i + 1):length(ms))
        end
    end


    @testset "Float16 angular coupling avoids intermediate overflow" begin
        src = middle_mode(Val(2), r -> Float16(0); Slm=r -> Float16(0.1),
                          T=Float16, l=256, m=0, omega=Float16(0),
                          fp=GHZ_MIDDLE_FP + 30)
        head = solve_irg_corrector_head(src, InnerRayIVP(src; rb=Float16(2)))
        sol = solve_irg_corrector_middle(head, MiddleRayIVP(head))
        @test isfinite(x_nm(sol, Float16(3)))
        @test isfinite(dr_x_nm(sol, Float16(3)))
        @test isfinite(irg_middle_residual(Float16(3), Float16(0.2), Float16(0.1),
                                           Float16(-0.05), Float16(0.1), Float16(0), 256))
    end
end

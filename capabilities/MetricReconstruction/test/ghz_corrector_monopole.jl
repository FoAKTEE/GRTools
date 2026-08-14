using Test
using MetricReconstruction
using MetricReconstruction.GHZ

const GHZ_MONO = MetricReconstruction.GHZ
const GHZ_MONO_SOURCE_FP = UInt128(0x6b9129c415e64ebcac0ba3b1d4e03b71)
const GHZ_MONO_BOUNDARY_FP = UInt128(0xd9e3145b6fd24fcd89a325850ac623a1)

function monopole_head(::Val{N}, ::Type{T}=Float64; M=one(T), l=0, m=0,
        omega=zero(Complex{T}), rb=T(2) * M, head_x0=zero(Complex{T}),
        head_dx0=zero(Complex{T}), source_fingerprint=GHZ_MONO_SOURCE_FP,
        Sll=q -> zero(q), Slm=GHZ_MONO.ZeroRadialSource(),
        Sln=GHZ_MONO.ZeroRadialSource()) where {N,T<:AbstractFloat}
    src = TetradSourceMode{N}(Sll; Slm, Sln, M, l, m, omega,
        source_fingerprint)
    solve_irg_corrector_head(src,
        InnerRayIVP(src; rb, x_mmb=head_x0, dr_x_mmb=head_dx0))
end

function monopole_solution(head; x0=0,
        boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
        policy=FiniteRayFreeIVP())
    bc = GHZ_MONO.MonopoleTailRayIVP(head; x_nn=x0,
        boundary_fingerprint, policy)
    GHZ_MONO.solve_irg_corrector_monopole_tail(head, bc)
end

function raw_monopole_boundary(bc::GHZ_MONO.MonopoleTailRayIVP{N,T};
        M=bc.M, l=bc.l, m=bc.m, omega=bc.omega,
        source_fingerprint=bc.source_fingerprint, head_x0=bc.head_x0,
        head_dx0=bc.head_dx0, rb=bc.rb, x0=bc.x0,
        boundary_fingerprint=bc.boundary_fingerprint,
        policy=bc.policy) where {N,T}
    GHZ_MONO.MonopoleTailRayIVP{N,T}(M, l, m, omega,
        source_fingerprint, head_x0, head_dx0, rb, x0,
        boundary_fingerprint, policy)
end

function raw_head_labels(head::IRGCorrectorHeadMode{N,T}, l, m) where {N,T}
    src, bc = head.source, head.boundary
    raw_src = typeof(src)(src.Sll, src.Slm, src.Sln, src.M, l, m,
        src.omega, src.source_fingerprint, src.normalization, src.tetrad,
        src.slicing, src.phase)
    raw_bc = InnerRayIVP{N,T}(bc.M, l, m, bc.omega,
        bc.source_fingerprint, bc.rb, bc.x0, bc.dx0)
    IRGCorrectorHeadMode(raw_src, raw_bc, GHZIRG())
end

@testset "GHZ — head-only outgoing-Schwarzschild monopole tail" begin
    @testset "finite-ray policy and homogeneous C/r" begin
        head = monopole_head(Val(2); omega=0.3 - 0.2im)
        mode = monopole_solution(head; x0=1.25 - 0.5im)
        @test perturbative_order(mode) == 2
        @test mode.boundary.policy isa FiniteRayFreeIVP
        for r in (2.0, 3.5, 8.0)
            want = mode.boundary.rb * mode.boundary.x0 / r
            @test x_nn(mode, r) ≈ want atol=2e-15
            @test dr_x_nn(mode, r) ≈ -want / r atol=2e-15
            @test irg_tail_residual(r, want, -want / r, 0.0) == 0
        end
        @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(head; x_nn=0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP, policy=:retarded)
        @test_throws UndefKeywordError GHZ_MONO.MonopoleTailRayIVP(head; x_nn=0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP)

        F = floatmax(Float64)
        scaled = monopole_solution(head; x0=0.75F,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 1)
        @test x_nn(scaled, 4.0) == 0.375F + 0im
    end

    @testset "dynamic A=1 oracle and sign" begin
        omega = 0.2 - 0.1im
        head = monopole_head(Val(2); omega, head_x0=1.0)
        x0 = 0.4 + 0.3im
        mode = monopole_solution(head; x0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 2)
        rb = mode.boundary.rb
        for r in (rb, 2.75, 4.5)
            Q0 = 2im * omega / r + inv(r^2)
            want = (rb * x0 + im * omega * (r^2 - rb^2) + (r - rb)) / r
            dwant = r * Q0 - want / r
            @test x_nn(mode, r) ≈ want rtol=3e-13 atol=3e-14
            @test dr_x_nn(mode, r) ≈ dwant rtol=3e-13 atol=3e-14
            @test irg_tail_residual(r, want, dwant, Q0) ≈ 0 atol=3e-15
            r > rb && @test abs(irg_tail_residual(r, -want, -dwant, Q0)) > 0.1
        end

        coefficients = GHZ_MONO.corrector_monopole_coefficients(mode, 3.0)
        @test coefficients isa GHZ_MONO.IRGCorrectorMonopoleCoefficients{2,Float64}
        @test coefficients.A == 1.0 + 0im
        @test coefficients.N ≈ x_nn(mode, 3.0) rtol=3e-13
        @test coefficients.spin_one_components isa
            GHZ_MONO.MonopoleSpinOneComponentsAbsent
        names = fieldnames(typeof(coefficients))
        @test :A in names && :N in names
        @test all(name -> !(name in names),
            (:B, :C, :partner, :pair_fingerprint, :partner_source_fingerprint))
        @test (coefficients.l, coefficients.m) == (0, 0)
        @test coefficients.source_fingerprint == head.source.source_fingerprint
        @test coefficients.tail_boundary_fingerprint ==
            mode.boundary.boundary_fingerprint
        @test coefficients.policy isa FiniteRayFreeIVP
        @test coefficients.gauge isa GHZIRG
    end

    @testset "manufactured A=r full-chain oracle" begin
        C = 0.35 - 0.2im
        omega = 0.15 + 0.05im
        head = monopole_head(Val(2); omega, Sll=r -> -2 / r,
            Sln=r -> C, head_x0=2.0, head_dx0=1.0,
            source_fingerprint=GHZ_MONO_SOURCE_FP + 20)
        N0 = 0.4 - 0.15im
        mode = monopole_solution(head; x0=N0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 20)
        rb = mode.boundary.rb
        for r in (rb, 2.75, 4.5)
            Q0 = C + 3im * omega + 3 / r - 3 / r^2
            want = (rb * N0 + (C + 3im * omega) * (r^3 - rb^3) / 3 +
                3 * (r^2 - rb^2) / 2 - 3 * (r - rb)) / r
            dwant = r * Q0 - want / r
            @test x_mmb(head, r) ≈ r rtol=3e-13 atol=3e-14
            @test dr_x_mmb(head, r) ≈ 1 rtol=3e-13 atol=3e-14
            @test x_nn(mode, r) ≈ want rtol=8e-12 atol=5e-13
            @test dr_x_nn(mode, r) ≈ dwant rtol=8e-12 atol=5e-13
            @test irg_tail_residual(r, want, dwant, Q0) ≈ 0 atol=5e-15
            r > rb && @test abs(
                irg_tail_residual(r, -want, -dwant, Q0)) > 0.1
        end
    end

    @testset "Float16 long-ray scale gate" begin
        T = Float16
        head = monopole_head(Val(2), T; omega=zero(Complex{T}),
            Sln=r -> inv(r), source_fingerprint=GHZ_MONO_SOURCE_FP + 21)
        mode = monopole_solution(head; x0=zero(T),
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 21)
        @test x_nn(mode, T(400)) == Complex{T}(T(200), zero(T))
    end

    @testset "Float16 dynamic-coefficient scale gate" begin
        T = Float16
        omega = Complex{T}(T(20000), zero(T))
        head = monopole_head(Val(2), T; omega, head_x0=one(T),
            source_fingerprint=GHZ_MONO_SOURCE_FP + 22)
        mode = monopole_solution(head; x0=zero(T),
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 22)
        @test dr_x_nn(mode, T(2)) == Complex{T}(T(0.5), T(40000))
    end

    @testset "subnormal-radius scale gate" begin
        M = 1.0e-310
        head = monopole_head(Val(2), Float64; M, rb=2M,
            source_fingerprint=GHZ_MONO_SOURCE_FP + 23)
        mode = monopole_solution(head; x0=0.0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 23)
        @test x_nn(mode, 2M) == 0.0 + 0.0im
        @test dr_x_nn(mode, 2M) == 0.0 + 0.0im

        M2 = 5.0e-309
        head2 = monopole_head(Val(2), Float64; M=M2, rb=2M2,
            head_dx0=2.6, source_fingerprint=GHZ_MONO_SOURCE_FP + 24)
        mode2 = monopole_solution(head2; x0=0.0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 24)
        @test dr_x_nn(mode2, 2M2) ≈ 1.3 rtol=2eps(Float64)

        r3 = floatmin(Float64)
        omega3 = ComplexF64(0.0, inv(r3) / 2)
        head3 = monopole_head(Val(2), Float64; M=r3 / 2, rb=r3,
            omega=omega3, head_x0=1.0,
            source_fingerprint=GHZ_MONO_SOURCE_FP + 25)
        mode3 = monopole_solution(head3; x0=0.0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 25)
        @test dr_x_nn(mode3, r3) == 0.0 + 0.0im

        head4 = monopole_head(Val(2), Float64; M=r3 / 2, rb=r3,
            omega=omega3, head_dx0=8.0,
            source_fingerprint=GHZ_MONO_SOURCE_FP + 26)
        mode4 = monopole_solution(head4; x0=0.0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 26)
        @test dr_x_nn(mode4, r3) == 0.0 + 0.0im
    end

    @testset "complex BigFloat-256 oracle" begin
        setprecision(BigFloat, 256) do
            omega = Complex{BigFloat}(big"0.1875", big"-0.0625")
            x0 = Complex{BigFloat}(big"0.25", big"-0.125")
            head = monopole_head(Val(2), BigFloat; omega,
                head_x0=one(BigFloat))
            mode = monopole_solution(head; x0,
                boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 3)
            r, rb = big"3.75", mode.boundary.rb
            want = (rb * x0 + Complex{BigFloat}(0, 1) * omega *
                (r^2 - rb^2) + (r - rb)) / r
            got = x_nn(mode, r; rtol=big"1e-55")
            dgot = dr_x_nn(mode, r; rtol=big"1e-55")
            Q0 = Complex{BigFloat}(0, 2) * omega / r + inv(r^2)
            @test got isa Complex{BigFloat}
            @test dgot isa Complex{BigFloat}
            @test abs(got - want) < big"1e-50"
            @test abs(dgot - (r * Q0 - want / r)) < big"1e-50"
            @test abs(irg_tail_residual(r, got, dgot, Q0)) < big"1e-65"
        end
    end

    @testset "S_lm is structurally absent and never called" begin
        arbitrary = monopole_head(Val(2); Slm=q -> zero(q))
        @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(arbitrary;
            x_nn=0, boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
            policy=FiniteRayFreeIVP())

        calls = Ref(0)
        poison(q) = (calls[] += 1; error("S_lm must never be called"))
        poisoned = monopole_head(Val(2); Slm=poison)
        error = try
            GHZ_MONO.MonopoleTailRayIVP(poisoned; x_nn=0,
                boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
                policy=FiniteRayFreeIVP())
            nothing
        catch exception
            exception
        end
        @test error isa ArgumentError
        @test calls[] == 0

        raw_bc = GHZ_MONO.MonopoleTailRayIVP{2,Float64}(
            poisoned.source.M, 0, 0, poisoned.source.omega,
            poisoned.source.source_fingerprint, poisoned.boundary.x0,
            poisoned.boundary.dx0, poisoned.boundary.rb, 0.0 + 0im,
            GHZ_MONO_BOUNDARY_FP + 10, FiniteRayFreeIVP())
        raw_mode = GHZ_MONO.IRGCorrectorMonopoleTailMode(
            poisoned, raw_bc, GHZIRG())
        @test_throws ArgumentError x_nn(raw_mode, raw_bc.rb)
        @test calls[] == 0
    end

    @testset "mode, boundary, scalar, and raw rejection" begin
        for (l, m) in ((1, 0), (1, 1))
            bad_head = monopole_head(Val(2); l, m)
            @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(bad_head;
                x_nn=0, boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
                policy=FiniteRayFreeIVP())
        end
        head = monopole_head(Val(2))
        raw_m = raw_head_labels(head, 0, 1)
        @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(raw_m;
            x_nn=0, boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
            policy=FiniteRayFreeIVP())

        for bad in (true, Inf, Complex{Bool}(true, false))
            @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(head;
                x_nn=bad, boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
                policy=FiniteRayFreeIVP())
        end
        for fp in (0, true, 2.0)
            @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(head;
                x_nn=0, boundary_fingerprint=fp, policy=FiniteRayFreeIVP())
        end
        for rb in (true, Inf, 1.9, 2.5)
            @test_throws ArgumentError GHZ_MONO.MonopoleTailRayIVP(head; rb,
                x_nn=0, boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
                policy=FiniteRayFreeIVP())
        end

        bc = GHZ_MONO.MonopoleTailRayIVP(head; x_nn=0.2 - 0.1im,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP,
            policy=FiniteRayFreeIVP())
        mode = GHZ_MONO.solve_irg_corrector_monopole_tail(head, bc)
        for bad_r in (true, Complex{Bool}(true, false), Inf, 1.9)
            @test_throws ArgumentError x_nn(mode, bad_r)
            @test_throws ArgumentError dr_x_nn(mode, bad_r)
            @test_throws ArgumentError GHZ_MONO.corrector_monopole_coefficients(
                mode, bad_r)
        end
        for kwargs in ((; atol=true), (; rtol=true), (; atol=NaN),
                (; rtol=Inf), (; atol=0.0, rtol=0.0))
            @test_throws ArgumentError x_nn(mode, 2.0; kwargs...)
            @test_throws ArgumentError dr_x_nn(mode, 2.0; kwargs...)
            @test_throws ArgumentError GHZ_MONO.corrector_monopole_coefficients(
                mode, 2.0; kwargs...)
        end

        invalid_boundaries = (
            raw_monopole_boundary(bc; M=2.0),
            raw_monopole_boundary(bc; l=1),
            raw_monopole_boundary(bc; m=1),
            raw_monopole_boundary(bc; omega=1.0 + 0im),
            raw_monopole_boundary(bc;
                source_fingerprint=bc.source_fingerprint + 1),
            raw_monopole_boundary(bc; head_x0=bc.head_x0 + 1),
            raw_monopole_boundary(bc; head_dx0=bc.head_dx0 + 1),
            raw_monopole_boundary(bc; rb=2.5),
            raw_monopole_boundary(bc; x0=Inf + 0im),
            raw_monopole_boundary(bc; boundary_fingerprint=0),
        )
        for invalid in invalid_boundaries
            @test_throws ArgumentError GHZ_MONO.solve_irg_corrector_monopole_tail(
                head, invalid)
            raw = GHZ_MONO.IRGCorrectorMonopoleTailMode(head, invalid, GHZIRG())
            @test_throws ArgumentError x_nn(raw, max(2.0, invalid.rb))
            @test_throws ArgumentError dr_x_nn(raw, max(2.0, invalid.rb))
            @test_throws ArgumentError GHZ_MONO.corrector_monopole_coefficients(
                raw, max(2.0, invalid.rb))
        end

        head1 = monopole_head(Val(1); source_fingerprint=GHZ_MONO_SOURCE_FP + 1)
        bc1 = GHZ_MONO.MonopoleTailRayIVP(head1; x_nn=0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 1,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError GHZ_MONO.solve_irg_corrector_monopole_tail(
            head, bc1)
        big_head = monopole_head(Val(2), BigFloat;
            source_fingerprint=GHZ_MONO_SOURCE_FP + 4)
        big_bc = GHZ_MONO.MonopoleTailRayIVP(big_head; x_nn=0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 4,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError GHZ_MONO.solve_irg_corrector_monopole_tail(
            head, big_bc)

        bad_sln = monopole_head(Val(2); Sln=q -> true,
            source_fingerprint=GHZ_MONO_SOURCE_FP + 2)
        bad_sln_bc = GHZ_MONO.MonopoleTailRayIVP(bad_sln; x_nn=0,
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 2,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError GHZ_MONO.solve_irg_corrector_monopole_tail(
            bad_sln, bad_sln_bc)
        raw_bad_sln = GHZ_MONO.IRGCorrectorMonopoleTailMode(
            bad_sln, bad_sln_bc, GHZIRG())
        @test_throws ArgumentError x_nn(raw_bad_sln, bad_sln_bc.rb)
        @test_throws ArgumentError dr_x_nn(raw_bad_sln, bad_sln_bc.rb)

        late = monopole_solution(monopole_head(Val(2);
            Sln=q -> q > 2.2 ? Inf : 0.0,
            source_fingerprint=GHZ_MONO_SOURCE_FP + 3);
            boundary_fingerprint=GHZ_MONO_BOUNDARY_FP + 3)
        @test_throws ArgumentError x_nn(late, 3.0)
        @test_throws ArgumentError dr_x_nn(late, 3.0)
    end

    @testset "dispatch has no internal ambiguities" begin
        for callable in (GHZ_MONO.MonopoleTailRayIVP,
                GHZ_MONO.solve_irg_corrector_monopole_tail, x_nn, dr_x_nn,
                GHZ_MONO.corrector_monopole_coefficients)
            ms = collect(methods(callable))
            @test all(!Base.isambiguous(ms[i], ms[j])
                for i in eachindex(ms) for j in (i + 1):length(ms))
        end
    end
end

using Test
using MetricReconstruction
using MetricReconstruction.GHZ

const GHZ_ADJ_PRIMARY_FP = UInt128(0x6c55c27f89dd49f0bb3497e3f26359b1)
const GHZ_ADJ_PARTNER_FP = UInt128(0x52ac8f24946d4d178a82e78da51d47e3)
const GHZ_ADJ_PAIR_FP = UInt128(0x027e1ab8326d4c4092d2797f4bc8e6ab)
const GHZ_ADJ_SOURCE_FP = UInt128(0x9e1340b9cfb54652a9ae6d55bc35f826)
const GHZ_ADJ_BOUNDARY_FP = UInt128(0xf023706d8b6243cc8066c40f5981200d)

function adjoint_middle_mode(::Val{N}, ::Type{T}=Float64; M=one(T), l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), rb=T(2) * M,
        source_fingerprint=GHZ_ADJ_PRIMARY_FP, head_x0=zero(Complex{T}),
        head_dx0=zero(Complex{T}), middle_x0=zero(Complex{T}),
        middle_dx0=zero(Complex{T}), Sll=nothing, Slm=nothing, Sln=nothing) where
        {N,T<:AbstractFloat}
    sll = isnothing(Sll) ? (q -> zero(q)) : Sll
    slm = isnothing(Slm) ? (q -> zero(q)) : Slm
    sln = isnothing(Sln) ? (q -> zero(q)) : Sln
    src = TetradSourceMode{N}(sll; Slm=slm, Sln=sln, M, l, m, omega,
        source_fingerprint)
    head = solve_irg_corrector_head(src,
        InnerRayIVP(src; rb, x_mmb=head_x0, dr_x_mmb=head_dx0))
    solve_irg_corrector_middle(head,
        MiddleRayIVP(head; x_nm=middle_x0, dr_x_nm=middle_dx0))
end

function adjoint_pair(::Val{N}=Val(2), ::Type{T}=Float64; l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), primary_kwargs=NamedTuple(),
        partner_kwargs=NamedTuple(), pair_fingerprint=GHZ_ADJ_PAIR_FP) where
        {N,T<:AbstractFloat}
    primary = adjoint_middle_mode(Val(N), T; l, m, omega,
        source_fingerprint=GHZ_ADJ_PRIMARY_FP, primary_kwargs...)
    partner = adjoint_middle_mode(Val(N), T; l, m=-m, omega=-conj(omega),
        source_fingerprint=GHZ_ADJ_PARTNER_FP, partner_kwargs...)
    PhysicalCorrectorMiddlePair(primary, partner; pair_fingerprint)
end

function adjoint_source(pair, Smbmb; adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP,
        kwargs...)
    ReducedAdjointSourceMode(pair, Smbmb; adjusted_source_fingerprint, kwargs...)
end

function eta_solution(pair, Smbmb; eta0=0, deta0=0,
        adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP,
        boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP, rb=pair.primary.boundary.rb,
        policy=FiniteRayFreeIVP())
    src = adjoint_source(pair, Smbmb; adjusted_source_fingerprint)
    bc = EtaRayIVP(src; rb, eta=eta0, dr_eta=deta0, boundary_fingerprint, policy)
    solve_irg_eta(src, bc)
end

function raw_adjoint_source(src::ReducedAdjointSourceMode{N,T};
        Smbmb=src.Smbmb, pair=src.pair, M=src.M, l=src.l, m=src.m,
        omega=src.omega, primary_source_fingerprint=src.primary_source_fingerprint,
        partner_source_fingerprint=src.partner_source_fingerprint,
        pair_fingerprint=src.pair_fingerprint, primary_head_x0=src.primary_head_x0,
        primary_head_dx0=src.primary_head_dx0,
        primary_middle_x0=src.primary_middle_x0,
        primary_middle_dx0=src.primary_middle_dx0,
        partner_head_x0=src.partner_head_x0,
        partner_head_dx0=src.partner_head_dx0,
        partner_middle_x0=src.partner_middle_x0,
        partner_middle_dx0=src.partner_middle_dx0, pair_rb=src.pair_rb,
        adjusted_source_fingerprint=src.adjusted_source_fingerprint,
        normalization=src.normalization, gauge=src.gauge, tetrad=src.tetrad,
        slicing=src.slicing, phase=src.phase) where {N,T}
    P, F = typeof(pair), typeof(Smbmb)
    ReducedAdjointSourceMode{N,T,P,F}(Smbmb, pair, M, l, m, omega,
        primary_source_fingerprint, partner_source_fingerprint, pair_fingerprint,
        primary_head_x0, primary_head_dx0, primary_middle_x0, primary_middle_dx0,
        partner_head_x0, partner_head_dx0, partner_middle_x0, partner_middle_dx0,
        pair_rb, adjusted_source_fingerprint, normalization, gauge, tetrad, slicing, phase)
end

function raw_eta_boundary(bc::EtaRayIVP{N,T}; M=bc.M, l=bc.l, m=bc.m,
        omega=bc.omega, primary_source_fingerprint=bc.primary_source_fingerprint,
        partner_source_fingerprint=bc.partner_source_fingerprint,
        pair_fingerprint=bc.pair_fingerprint,
        adjusted_source_fingerprint=bc.adjusted_source_fingerprint, rb=bc.rb,
        eta0=bc.eta0, deta0=bc.deta0,
        boundary_fingerprint=bc.boundary_fingerprint, policy=bc.policy) where {N,T}
    EtaRayIVP{N,T}(M, l, m, omega, primary_source_fingerprint,
        partner_source_fingerprint, pair_fingerprint, adjusted_source_fingerprint,
        rb, eta0, deta0, boundary_fingerprint, policy)
end

@testset "GHZ — supplied reduced adjoint source to eta" begin
    @testset "exact sign and factor oracle: eta=r^2, S_mbar_mbar=-3" begin
        pair = adjoint_pair()
        mode = eta_solution(pair, r -> -3.0; eta0=4.0, deta0=4.0)
        @test mode.source.normalization isa TwoRealT0Dagger
        @test mode.boundary.policy isa FiniteRayFreeIVP
        @test perturbative_order(mode.source) == 2
        @test perturbative_order(mode.boundary) == 2
        @test perturbative_order(mode) == 2
        for r in (2.0, 2.5, 3.0, 7.0)
            @test eta(mode, r) ≈ r^2 rtol=2e-13 atol=3e-14
            @test dr_eta(mode, r) ≈ 2r rtol=2e-13 atol=3e-14
            @test d2r_eta(mode, r) ≈ 2.0 rtol=2e-13 atol=3e-14
            @test irg_eta_residual(r, 2r, 2.0, -3.0) == 0.0
            @test irg_eta_residual(r, 2r, 2.0, 3.0) == 12.0
        end
    end

    @testset "free homogeneous data at the shared finite rb" begin
        rb, a, b = 2.75, 1.2 - 0.4im, -0.3 + 0.7im
        pair = adjoint_pair(primary_kwargs=(; rb), partner_kwargs=(; rb),
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 1)
        mode = eta_solution(pair, r -> 0.0; rb, eta0=a, deta0=b,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 1,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 1)
        @test mode.boundary.rb == rb
        for r in (rb, 3.5, 8.0)
            want = a + rb * b * (1 - rb / r)
            dwant = b * (rb / r)^2
            d2want = -2b * rb^2 / r^3
            @test eta(mode, r) ≈ want rtol=2e-13 atol=3e-14
            @test dr_eta(mode, r) ≈ dwant rtol=2e-13 atol=3e-14
            @test d2r_eta(mode, r) ≈ d2want rtol=2e-13 atol=3e-14
            @test irg_eta_residual(r, dwant, d2want, 0) ≈ 0 atol=4e-16
        end
        @test_throws ArgumentError EtaRayIVP(mode.source; eta=0, dr_eta=0,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP, policy=:retarded)
        @test_throws UndefKeywordError EtaRayIVP(mode.source; eta=0, dr_eta=0,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP)
        @test_throws ArgumentError EtaRayIVP(mode.source; rb=rb + 0.25, eta=0,
            dr_eta=0, boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP,
            policy=FiniteRayFreeIVP())
    end

    @testset "complex source, BigFloat-256, and linearity" begin
        C = 0.25 - 0.375im
        pair = adjoint_pair(pair_fingerprint=GHZ_ADJ_PAIR_FP + 2)
        mode = eta_solution(pair, r -> C;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 2,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 2)
        rb = mode.boundary.rb
        for r in (rb, 3.0, 5.5)
            want = C * (-r^2 / 3 + rb^2 - 2rb^3 / (3r))
            dwant = -2C * (r - rb^3 / r^2) / 3
            d2want = -2C - 2dwant / r
            @test eta(mode, r) ≈ want rtol=4e-13 atol=4e-14
            @test dr_eta(mode, r) ≈ dwant rtol=4e-13 atol=4e-14
            @test d2r_eta(mode, r) ≈ d2want rtol=4e-13 atol=4e-14
        end

        C1, C2 = 0.2 - 0.1im, -0.35 + 0.15im
        e1, e2 = 0.4 + 0.2im, -0.1 + 0.3im
        d1, d2 = -0.2 + 0.5im, 0.25 - 0.15im
        alpha, beta = 1.1 - 0.4im, -0.3 + 0.2im
        function linear_mode(C, e, d, offset)
            p = adjoint_pair(pair_fingerprint=GHZ_ADJ_PAIR_FP + offset)
            eta_solution(p, r -> C; eta0=e, deta0=d,
                adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + offset,
                boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + offset)
        end
        u1, u2 = linear_mode(C1, e1, d1, 10), linear_mode(C2, e2, d2, 11)
        u12 = linear_mode(alpha * C1 + beta * C2, alpha * e1 + beta * e2,
            alpha * d1 + beta * d2, 12)
        for r in (2.0, 3.25, 6.0)
            @test eta(u12, r) ≈ alpha * eta(u1, r) + beta * eta(u2, r) rtol=7e-13
            @test dr_eta(u12, r) ≈
                alpha * dr_eta(u1, r) + beta * dr_eta(u2, r) rtol=7e-13
        end

        setprecision(BigFloat, 256) do
            pairb = adjoint_pair(Val(2), BigFloat;
                omega=Complex{BigFloat}(big"0.125", big"-0.0625"),
                pair_fingerprint=GHZ_ADJ_PAIR_FP + 20)
            modeb = eta_solution(pairb, r -> -big"3"; eta0=big"4", deta0=big"4",
                adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 20,
                boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 20)
            r = big"3.75"
            got = eta(modeb, r; rtol=big"1e-55")
            @test got isa Complex{BigFloat}
            @test dr_eta(modeb, r; rtol=big"1e-55") isa Complex{BigFloat}
            @test d2r_eta(modeb, r; rtol=big"1e-55") isa Complex{BigFloat}
            @test abs(got - r^2) < big"1e-50"
        end
    end

    @testset "scale-safe finite Green assembly" begin
        F = floatmax(Float64)
        minsub = nextfloat(0.0)
        tiniest_boundary = 2 * minsub
        @test MetricReconstruction.GHZ._eta_scaled_ray_gap(
            tiniest_boundary, F) == tiniest_boundary
        @test MetricReconstruction.GHZ._eta_scaled_ray_gap(minsub, F) == minsub

        # Select the algebraic ordering whose divided factor is at least one
        # half.  Here q/F is a nonzero subnormal, but using that rounded ratio
        # first loses every factor discriminator and returns the same value.
        unit_q = F * minsub
        for factor in (0.51, 0.6, 0.75)
            q = factor * unit_q
            rounded = setprecision(BigFloat, 8192) do
                value = BigFloat(q) * (BigFloat(F) - BigFloat(q)) / BigFloat(F)
                Float64(value)
            end
            @test !iszero(q / F)
            @test issubnormal(q / F)
            got = MetricReconstruction.GHZ._eta_scaled_ray_gap(q, F)
            @test !iszero(got)
            @test got == rounded
        end

        # Multiplying by two before dividing by r must not reject a finite
        # homogeneous second derivative.
        pair_d2 = adjoint_pair(pair_fingerprint=GHZ_ADJ_PAIR_FP + 70)
        mode_d2 = eta_solution(pair_d2, r -> 0.0; deta0=0.75F,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 70,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 70)
        @test d2r_eta(mode_d2, 2.0) == ComplexF64(-0.75F, 0.0)

        # Exercise both exact source/derivative cancellation and the fallback
        # in which deta/r itself overflows although the combined equation is
        # representable.
        small_M, small_rb = 0.25, 0.5
        pair_small = adjoint_pair(
            primary_kwargs=(; M=small_M, rb=small_rb),
            partner_kwargs=(; M=small_M, rb=small_rb),
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 71)
        S_cancel = 0.75F
        mode_cancel = eta_solution(pair_small, r -> S_cancel;
            deta0=-small_rb * S_cancel,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 71,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 71)
        @test d2r_eta(mode_cancel, small_rb) == 0.0 + 0.0im

        fallback_derivative = exp2(1023.0)
        mode_fallback = eta_solution(pair_small, r -> -F;
            deta0=fallback_derivative,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 72,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 72)
        @test d2r_eta(mode_fallback, small_rb) ==
            ComplexF64(-exp2(972.0), 0.0)

        # Rounding S+deta/r before the final factor of two is one ulp wrong in
        # this cancellation-generated subnormal.  Widen the complete original
        # expression, whose correctly rounded magnitude ends in hex 0x1555.
        hex_source = reinterpret(Float64, UInt64(0x0010000000000000))
        hex_derivative = -reinterpret(Float64, UInt64(0x0027fffffffff000))
        hex_magnitude = reinterpret(Float64, UInt64(0x0000000000001555))
        pair_hex = adjoint_pair(
            primary_kwargs=(; rb=3.0), partner_kwargs=(; rb=3.0),
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 79)
        mode_hex = eta_solution(pair_hex, r -> hex_source; deta0=hex_derivative,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 79,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 79)
        hex_d2 = d2r_eta(mode_hex, 3.0)
        @test hex_d2 == ComplexF64(-hex_magnitude, 0.0)
        @test reinterpret(UInt64, abs(real(hex_d2))) == UInt64(0x1555)

        # Keep the residual in the same scale-safe equation grouping as the
        # mode evaluator.  The first case must not form `2*deta` eagerly; the
        # second must recover a finite cancellation after `deta/r` overflows.
        @test irg_eta_residual(2.0, 0.75F, 0.0, 0.0) == 0.75F
        @test irg_eta_residual(small_rb, fallback_derivative,
            -exp2(972.0), -F) == 0.0
        hex_residual = irg_eta_residual(3.0, hex_derivative, 0.0, hex_source)
        @test hex_residual == hex_magnitude
        @test reinterpret(UInt64, abs(hex_residual)) == UInt64(0x1555)
        @test_throws ArgumentError irg_eta_residual(1.0, F, F, F)

        residual16 = irg_eta_residual(Float16(2), Float16(0.25),
            Float16(-0.5), Float16(0.125))
        @test residual16 isa Float16
        @test residual16 == zero(Float16)
        setprecision(BigFloat, 256) do
            residual_big = irg_eta_residual(big"2", big"0.25", big"-0.5",
                big"0.125")
            @test residual_big isa BigFloat
            @test iszero(residual_big)
        end

        # Each ratio multiplication must be applied before the next one.  A
        # precomputed ratio^2 underflows here even though the final result is
        # normal and exactly representable.
        long_ratio = exp2(-600.0)
        large_value = ComplexF64(exp2(1000.0), 0.0)
        @test MetricReconstruction.GHZ._eta_ratio_square_scale(
            long_ratio, large_value) == ComplexF64(exp2(-200.0), 0.0)
        mode_long_ratio = eta_solution(pair_d2, r -> 0.0; deta0=exp2(1000.0),
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 78,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 78)
        @test dr_eta(mode_long_ratio, exp2(601.0)) ==
            ComplexF64(exp2(-200.0), 0.0)

        # All three eta terms are individually large enough that the naive
        # homogeneous sum overflows.  The quadrature contract has already
        # rounded the source moment to F/2, so compose those stored Float64
        # values rather than imposing the nominal correlation between the
        # separately rounded `.75F` and `.09375F` inputs.
        pair_three_term = adjoint_pair(
            primary_kwargs=(; rb=4.0), partner_kwargs=(; rb=4.0),
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 77)
        mode_three_term = eta_solution(pair_three_term, r -> 0.09375F;
            deta0=0.75F,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 77,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 77)
        three_term_expected = setprecision(BigFloat, 8192) do
            Float64(BigFloat(2) * BigFloat(0.75F) -
                BigFloat(2) * BigFloat(0.5F))
        end
        @test three_term_expected == prevfloat(0.5F)
        @test eta(mode_three_term, 8.0) == ComplexF64(three_term_expected, 0.0)

        # The homogeneous increment must retain a one-ulp cell at a very large
        # shared boundary instead of subtracting nearly equal unit-scale ratios.
        adjacent_M, adjacent_rb = 5.0e307, 1.0e308
        pair_adjacent = adjoint_pair(
            primary_kwargs=(; M=adjacent_M, rb=adjacent_rb),
            partner_kwargs=(; M=adjacent_M, rb=adjacent_rb),
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 73)
        mode_adjacent = eta_solution(pair_adjacent, r -> 0.0; deta0=1.0,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 73,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 73)
        adjacent_r = nextfloat(adjacent_rb)
        adjacent_want = (adjacent_rb / adjacent_r) * (adjacent_r - adjacent_rb)
        @test eta(mode_adjacent, adjacent_r) == ComplexF64(adjacent_want, 0.0)

        tiny_source = nextfloat(0.0)
        source_cell_r = nextfloat(adjacent_rb, 1024)
        mode_source_cell = eta_solution(pair_adjacent, r -> tiny_source;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 76,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 76)
        source_cell_want = setprecision(BigFloat, 256) do
            rb, r = BigFloat(adjacent_rb), BigFloat(source_cell_r)
            delta = r - rb
            kernel = delta^2 / 2 - delta^3 / (3r)
            Float64(-BigFloat(2) * BigFloat(tiny_source) * kernel)
        end
        @test eta(mode_source_cell, source_cell_r) ≈
            ComplexF64(source_cell_want, 0.0) rtol=4eps(Float64)

        # Large but finite source moments may cancel finite homogeneous data.
        # Do not materialize twice the moment before that cancellation.
        large_M, large_rb = 5000.0, 10000.0
        pair_source = adjoint_pair(
            primary_kwargs=(; M=large_M, rb=large_rb),
            partner_kwargs=(; M=large_M, rb=large_rb),
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 74)
        source_value = 0.006F

        eta_r = 10015.0
        mode_eta_source = eta_solution(pair_source, r -> source_value; eta0=F,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 74,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 74)
        eta_want = setprecision(BigFloat, 256) do
            rb, r = BigFloat(large_rb), BigFloat(eta_r)
            kernel = r^2 / 6 - rb^2 / 2 + rb^3 / (3r)
            Float64(BigFloat(F) - BigFloat(2) * BigFloat(source_value) * kernel)
        end
        @test eta(mode_eta_source, eta_r) ≈ ComplexF64(eta_want, 0.0) rtol=5e-13

        derivative_r = 10100.0
        mode_dr_source = eta_solution(pair_source, r -> source_value; deta0=F,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 75,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 75)
        derivative_want = setprecision(BigFloat, 256) do
            rb, r = BigFloat(large_rb), BigFloat(derivative_r)
            kernel = (r^3 - rb^3) / (3r^2)
            Float64((rb / r)^2 * BigFloat(F) -
                BigFloat(2) * BigFloat(source_value) * kernel)
        end
        @test dr_eta(mode_dr_source, derivative_r) ≈
            ComplexF64(derivative_want, 0.0) rtol=5e-13
    end

    @testset "source is supplied and independent of corrector-tail data" begin
        poison(q) = error("S_ln must not participate in the reduced eta solve")
        pair = adjoint_pair(primary_kwargs=(; Sln=poison),
            partner_kwargs=(; Sln=poison), pair_fingerprint=GHZ_ADJ_PAIR_FP + 30)
        mode = eta_solution(pair, r -> -3.0; eta0=4, deta0=4,
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 30,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 30)
        @test eta(mode, 3.0) ≈ 9.0 atol=4e-14
        @test dr_eta(mode, 3.0) ≈ 6.0 atol=4e-14
        @test all(!occursin("tail", lowercase(String(name))) &&
            !occursin("x_nn", lowercase(String(name)))
            for name in fieldnames(typeof(mode.source)))
        @test all(!occursin("tail", lowercase(String(name))) &&
            !occursin("x_nn", lowercase(String(name)))
            for name in fieldnames(typeof(mode.boundary)))
    end

    @testset "constructor, convention, radius, and tolerance rejection" begin
        pair = adjoint_pair(pair_fingerprint=GHZ_ADJ_PAIR_FP + 40)
        for bad in (true, Inf, Complex{Bool}(true, false))
            @test_throws ArgumentError adjoint_source(pair, r -> bad;
                adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 40)
        end
        for fp in (0, true, 2.0, pair.primary.head.source.source_fingerprint,
                pair.partner.head.source.source_fingerprint, pair.pair_fingerprint)
            @test_throws ArgumentError adjoint_source(pair, r -> 0;
                adjusted_source_fingerprint=fp)
        end
        for kwargs in ((; normalization=EinsteinRHS()), (; gauge=HertzIRG()),
                (; tetrad=:wrong), (; slicing=:wrong), (; phase=:wrong))
            @test_throws ArgumentError ReducedAdjointSourceMode(pair, r -> 0;
                adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 40, kwargs...)
        end
        pair_l1 = adjoint_pair(l=1, m=1, pair_fingerprint=GHZ_ADJ_PAIR_FP + 41)
        @test_throws ArgumentError adjoint_source(pair_l1, r -> 0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 41)

        src = adjoint_source(pair, r -> 0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 42)
        for bad in (true, Inf, Complex{Bool}(true, false))
            @test_throws ArgumentError EtaRayIVP(src; eta=bad, dr_eta=0,
                boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 40,
                policy=FiniteRayFreeIVP())
            @test_throws ArgumentError EtaRayIVP(src; eta=0, dr_eta=bad,
                boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 40,
                policy=FiniteRayFreeIVP())
        end
        for fp in (0, true, 2.0)
            @test_throws ArgumentError EtaRayIVP(src; eta=0, dr_eta=0,
                boundary_fingerprint=fp, policy=FiniteRayFreeIVP())
        end
        for bad_rb in (true, Inf, Complex{Bool}(true, false), 1.9)
            @test_throws ArgumentError EtaRayIVP(src; rb=bad_rb, eta=0, dr_eta=0,
                boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 40,
                policy=FiniteRayFreeIVP())
        end

        bc = EtaRayIVP(src; eta=0, dr_eta=0,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 42,
            policy=FiniteRayFreeIVP())
        mode = solve_irg_eta(src, bc)
        for bad_r in (true, Inf, Complex{Bool}(true, false), 1.9)
            @test_throws ArgumentError eta(mode, bad_r)
            @test_throws ArgumentError dr_eta(mode, bad_r)
            @test_throws ArgumentError d2r_eta(mode, bad_r)
        end
        for kwargs in ((; atol=true), (; rtol=true), (; atol=NaN), (; rtol=Inf),
                (; atol=0.0, rtol=0.0))
            @test_throws ArgumentError eta(mode, 2.0; kwargs...)
            @test_throws ArgumentError dr_eta(mode, 2.0; kwargs...)
            @test_throws ArgumentError d2r_eta(mode, 2.0; kwargs...)
        end
        @test_throws ArgumentError irg_eta_residual(true, 0, 0, 0)
        @test_throws ArgumentError irg_eta_residual(2.0, true, 0, 0)
        @test_throws ArgumentError irg_eta_residual(2.0, 0, Inf, 0)
        @test_throws ArgumentError irg_eta_residual(2.0, 0, 0, Complex{Bool}(true, false))
    end

    @testset "late invalid callable and complete raw revalidation" begin
        pair = adjoint_pair(pair_fingerprint=GHZ_ADJ_PAIR_FP + 50)
        late = eta_solution(pair, r -> r > 2.2 ? Inf : 0.0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 50,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 50)
        @test_throws ArgumentError eta(late, 3.0)
        @test_throws ArgumentError dr_eta(late, 3.0)
        @test_throws ArgumentError d2r_eta(late, 3.0)
        late_bool = eta_solution(pair, r -> r > 2.2 ? true : 0.0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 51,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 51)
        @test_throws ArgumentError eta(late_bool, 3.0)
        endpoint_only = eta_solution(pair, r -> r == 3.0 ? Inf : 0.0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 53,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 53)
        @test_throws ArgumentError eta(endpoint_only, 3.0)
        @test_throws ArgumentError dr_eta(endpoint_only, 3.0)
        @test_throws ArgumentError d2r_eta(endpoint_only, 3.0)

        src = adjoint_source(pair, r -> 0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 52)
        bad_sources = (
            raw_adjoint_source(src; M=2.0),
            raw_adjoint_source(src; l=3),
            raw_adjoint_source(src;
                primary_source_fingerprint=src.primary_source_fingerprint + 1),
            raw_adjoint_source(src; pair_fingerprint=src.pair_fingerprint + 1),
            raw_adjoint_source(src; primary_middle_x0=src.primary_middle_x0 + 1),
            raw_adjoint_source(src; partner_head_dx0=src.partner_head_dx0 + 1),
            raw_adjoint_source(src; pair_rb=src.pair_rb + 1),
            raw_adjoint_source(src; adjusted_source_fingerprint=0),
            raw_adjoint_source(src;
                adjusted_source_fingerprint=src.primary_source_fingerprint),
        )
        for bad in bad_sources
            @test_throws ArgumentError EtaRayIVP(bad; eta=0, dr_eta=0,
                boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 52,
                policy=FiniteRayFreeIVP())
        end

        bc = EtaRayIVP(src; eta=0.2 - 0.1im, dr_eta=-0.3 + 0.4im,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 52,
            policy=FiniteRayFreeIVP())
        bad_boundaries = (
            raw_eta_boundary(bc; M=2.0),
            raw_eta_boundary(bc; l=3),
            raw_eta_boundary(bc;
                partner_source_fingerprint=bc.partner_source_fingerprint + 1),
            raw_eta_boundary(bc; pair_fingerprint=bc.pair_fingerprint + 1),
            raw_eta_boundary(bc;
                adjusted_source_fingerprint=bc.adjusted_source_fingerprint + 1),
            raw_eta_boundary(bc; rb=1.5),
            raw_eta_boundary(bc; eta0=Complex(Inf, 0)),
            raw_eta_boundary(bc; deta0=Complex(NaN, 0)),
            raw_eta_boundary(bc; boundary_fingerprint=0),
        )
        for bad in bad_boundaries
            @test_throws ArgumentError solve_irg_eta(src, bad)
            raw_mode = IRGEtaMode(src, bad, GHZIRG())
            @test_throws ArgumentError eta(raw_mode, max(2.0, bad.rb))
            @test_throws ArgumentError dr_eta(raw_mode, max(2.0, bad.rb))
            @test_throws ArgumentError d2r_eta(raw_mode, max(2.0, bad.rb))
        end
    end

    @testset "order and precision separation" begin
        pair2 = adjoint_pair(pair_fingerprint=GHZ_ADJ_PAIR_FP + 60)
        src2 = adjoint_source(pair2, r -> 0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 60)
        bc2 = EtaRayIVP(src2; eta=0, dr_eta=0,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 60,
            policy=FiniteRayFreeIVP())

        pair1 = adjoint_pair(Val(1); pair_fingerprint=GHZ_ADJ_PAIR_FP + 61)
        src1 = adjoint_source(pair1, r -> 0;
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 61)
        bc1 = EtaRayIVP(src1; eta=0, dr_eta=0,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 61,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_eta(src2, bc1)

        pairb = adjoint_pair(Val(2), BigFloat;
            pair_fingerprint=GHZ_ADJ_PAIR_FP + 62)
        srcb = adjoint_source(pairb, r -> zero(r);
            adjusted_source_fingerprint=GHZ_ADJ_SOURCE_FP + 62)
        bcb = EtaRayIVP(srcb; eta=0, dr_eta=0,
            boundary_fingerprint=GHZ_ADJ_BOUNDARY_FP + 62,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_eta(src2, bcb)
        @test solve_irg_eta(src2, bc2) isa IRGEtaMode{2,Float64}
    end

    @testset "dispatch has no internal ambiguities" begin
        for callable in (ReducedAdjointSourceMode, EtaRayIVP, solve_irg_eta,
                eta, dr_eta, d2r_eta)
            ms = collect(methods(callable))
            @test all(!Base.isambiguous(ms[i], ms[j])
                for i in eachindex(ms) for j in (i + 1):length(ms))
        end
    end
end

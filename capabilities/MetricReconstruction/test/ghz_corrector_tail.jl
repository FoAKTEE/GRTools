using Test
using MetricReconstruction.GHZ

const GHZ_TAIL_PRIMARY_FP = UInt128(0x4698be73c4ad45539e46a6d627dc0c51)
const GHZ_TAIL_PARTNER_FP = UInt128(0xced57f063cc244b793bb8dc2c9770e12)
const GHZ_TAIL_PAIR_FP = UInt128(0x358372b4db58418583b88db26a558d94)
const GHZ_TAIL_BOUNDARY_FP = UInt128(0xa6b5d31112894571a77f1044af82ab28)

function tail_middle_mode(::Val{N}, ::Type{T}=Float64; M=one(T), l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), rb=T(2) * M,
        source_fingerprint=GHZ_TAIL_PRIMARY_FP, head_x0=zero(Complex{T}),
        head_dx0=zero(Complex{T}), middle_x0=zero(Complex{T}),
        middle_dx0=zero(Complex{T}), Sll=nothing, Slm=nothing, Sln=nothing) where
        {N,T<:AbstractFloat}
    sll = isnothing(Sll) ? (q -> zero(q)) : Sll
    slm = isnothing(Slm) ? (q -> zero(q)) : Slm
    common = (; M, l, m, omega, source_fingerprint)
    src = isnothing(Sln) ? TetradSourceMode{N}(sll; Slm=slm, common...) :
                           TetradSourceMode{N}(sll; Slm=slm, Sln, common...)
    head = solve_irg_corrector_head(src,
        InnerRayIVP(src; rb, x_mmb=head_x0, dr_x_mmb=head_dx0))
    solve_irg_corrector_middle(head,
        MiddleRayIVP(head; x_nm=middle_x0, dr_x_nm=middle_dx0))
end

function tail_pair(primary, partner; pair_fingerprint=GHZ_TAIL_PAIR_FP)
    PhysicalCorrectorMiddlePair(primary, partner; pair_fingerprint)
end

function tail_solution(primary, partner; x0=0, pair_fingerprint=GHZ_TAIL_PAIR_FP,
        boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP, policy=FiniteRayFreeIVP())
    pair = tail_pair(primary, partner; pair_fingerprint)
    bc = TailRayIVP(pair; x_nn=x0, boundary_fingerprint, policy)
    solve_irg_corrector_tail(pair, bc)
end

function zero_tail_inputs(::Val{N}=Val(2), ::Type{T}=Float64; l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), primary_kwargs=NamedTuple(),
        partner_kwargs=NamedTuple()) where {N,T<:AbstractFloat}
    primary = tail_middle_mode(Val(N), T; l, m, omega,
        source_fingerprint=GHZ_TAIL_PRIMARY_FP, primary_kwargs...)
    partner = tail_middle_mode(Val(N), T; l, m=-m, omega=-conj(omega),
        source_fingerprint=GHZ_TAIL_PARTNER_FP, partner_kwargs...)
    primary, partner
end

function raw_tail_boundary(bc::TailRayIVP{N,T}; M=bc.M, l=bc.l, m=bc.m,
        omega=bc.omega, primary_source_fingerprint=bc.primary_source_fingerprint,
        partner_source_fingerprint=bc.partner_source_fingerprint,
        pair_fingerprint=bc.pair_fingerprint, primary_head_x0=bc.primary_head_x0,
        primary_head_dx0=bc.primary_head_dx0,
        primary_middle_x0=bc.primary_middle_x0,
        primary_middle_dx0=bc.primary_middle_dx0,
        partner_head_x0=bc.partner_head_x0,
        partner_head_dx0=bc.partner_head_dx0,
        partner_middle_x0=bc.partner_middle_x0,
        partner_middle_dx0=bc.partner_middle_dx0, rb=bc.rb, x0=bc.x0,
        boundary_fingerprint=bc.boundary_fingerprint, policy=bc.policy) where {N,T}
    TailRayIVP{N,T}(M, l, m, omega, primary_source_fingerprint,
        partner_source_fingerprint, pair_fingerprint, primary_head_x0,
        primary_head_dx0, primary_middle_x0, primary_middle_dx0, partner_head_x0,
        partner_head_dx0, partner_middle_x0, partner_middle_dx0, rb, x0,
        boundary_fingerprint, policy)
end

function manufactured_tail()
    root3 = sqrt(3.0)
    omega = 0.2 + 0.0im
    primary = tail_middle_mode(Val(2); l=2, m=2, omega,
        Sll=r -> -6 + 8 / r,
        Slm=r -> root3 * (r - 2) / r - 2 / (root3 * r^2), Sln=r -> 1.0,
        middle_dx0=inv(root3), source_fingerprint=GHZ_TAIL_PRIMARY_FP + 10)
    partner = tail_middle_mode(Val(2); l=2, m=-2, omega=-conj(omega),
        Slm=r -> -4 / (root3 * r^2), middle_dx0=2 / root3,
        source_fingerprint=GHZ_TAIL_PARTNER_FP + 10)
    tail_solution(primary, partner; x0=1.5 - 0.25im,
        pair_fingerprint=GHZ_TAIL_PAIR_FP + 10,
        boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 10)
end

@testset "GHZ — bounded outgoing-Schwarzschild x_nn tail" begin
    @testset "finite-ray policy and homogeneous C/r" begin
        primary, partner = zero_tail_inputs()
        mode = tail_solution(primary, partner; x0=1.25 - 0.5im)
        @test mode.boundary.policy isa FiniteRayFreeIVP
        @test perturbative_order(mode) == 2
        for r in (2.0, 3.5, 8.0)
            want = mode.boundary.rb * mode.boundary.x0 / r
            @test x_nn(mode, r) ≈ want atol=2e-14
            @test dr_x_nn(mode, r) ≈ -want / r atol=2e-14
            @test irg_tail_residual(r, want, -want / r, 0.0) ≈ 0 atol=2e-15
        end
        @test_throws ArgumentError TailRayIVP(mode.pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP, policy=:retarded)
        @test_throws UndefKeywordError TailRayIVP(mode.pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP)

        l0 = TetradSourceMode{2}(r -> 0.0; M=1.0, l=0, m=0,
            source_fingerprint=GHZ_TAIL_PRIMARY_FP + 1)
        head0 = solve_irg_corrector_head(l0, InnerRayIVP(l0; rb=2.0))
        @test x_mmb(head0, 3.0) == 0.0 + 0.0im
        @test_throws ArgumentError MiddleRayIVP(head0)
    end

    @testset "constant complex S_ln sign and BigFloat" begin
        C = 0.3 - 0.2im
        primary, partner = zero_tail_inputs(primary_kwargs=(; Sln=r -> C))
        mode = tail_solution(primary, partner)
        for r in (2.0, 3.25, 6.0)
            want = C * (r^3 - 2.0^3) / (3r)
            dwant = r * C - want / r
            @test x_nn(mode, r) ≈ want rtol=3e-13 atol=2e-14
            @test dr_x_nn(mode, r) ≈ dwant rtol=3e-13 atol=2e-14
            @test irg_tail_residual(r, want, dwant, C) ≈ 0 atol=3e-15
            r > 2 && @test abs(irg_tail_residual(r, -want, -dwant, C)) > 0.5
        end

        setprecision(BigFloat, 256) do
            Cb = Complex{BigFloat}(big"0.1875", big"-0.0625")
            primary, partner = zero_tail_inputs(Val(2), BigFloat;
                omega=Complex{BigFloat}(big"0.05", big"-0.02"),
                primary_kwargs=(; Sln=r -> Cb))
            mode = tail_solution(primary, partner; x0=zero(Cb),
                pair_fingerprint=GHZ_TAIL_PAIR_FP + 2,
                boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 2)
            r = big"3.5"
            want = Cb * (r^3 - big"8") / (3r)
            got = x_nn(mode, r; rtol=big"1e-55")
            @test got isa Complex{BigFloat}
            @test dr_x_nn(mode, r; rtol=big"1e-55") isa Complex{BigFloat}
            @test abs(got - want) < big"1e-50"
        end
    end

    @testset "linearity" begin
        C1, C2 = 0.12 - 0.05im, -0.08 + 0.11im
        n1, n2 = 0.4 + 0.2im, -0.3 + 0.1im
        a, b = 1.2 - 0.3im, -0.4 + 0.5im
        function make(C, n0, offset)
            primary, partner = zero_tail_inputs(primary_kwargs=(; Sln=r -> C))
            tail_solution(primary, partner; x0=n0,
                pair_fingerprint=GHZ_TAIL_PAIR_FP + offset,
                boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + offset)
        end
        mode1, mode2 = make(C1, n1, 20), make(C2, n2, 21)
        mode12 = make(a * C1 + b * C2, a * n1 + b * n2, 22)
        for r in (2.0, 3.0, 5.0)
            @test x_nn(mode12, r) ≈ a * x_nn(mode1, r) + b * x_nn(mode2, r) rtol=5e-12
            @test dr_x_nn(mode12, r) ≈
                a * dr_x_nn(mode1, r) + b * dr_x_nn(mode2, r) rtol=5e-12
        end
    end

    @testset "full-chain polynomial oracle" begin
        mode = manufactured_tail()
        r = 3.0
        A, B, C = x_mmb(mode.pair.primary.head, r), x_nm(mode.pair.primary, r),
            x_nmb(mode.pair, r)
        @test A ≈ 1.0 atol=3e-13
        @test B ≈ inv(sqrt(3.0)) atol=3e-13
        @test C ≈ -2 / sqrt(3.0) atol=4e-13
        Q = 25 / 9 + 8im / 15
        want = 55 / 9 + 13im / 30
        dwant = 170 / 27 + 131im / 90
        @test x_nn(mode, r) ≈ want rtol=8e-11 atol=5e-13
        @test dr_x_nn(mode, r) ≈ dwant rtol=8e-11 atol=5e-13
        @test irg_tail_residual(r, want, dwant, Q) ≈ 0 atol=3e-15
    end

    @testset "partner dependence and tetrad assembly" begin
        primary, partner = zero_tail_inputs(primary_kwargs=(;
            middle_x0=0.4 - 0.2im, middle_dx0=-0.1 + 0.3im))
        partner2 = tail_middle_mode(Val(2); l=2, m=-1,
            omega=-conj(primary.head.source.omega), middle_x0=1.1 + 0.7im,
            middle_dx0=-0.6 + 0.2im, source_fingerprint=GHZ_TAIL_PARTNER_FP + 30)
        mode1 = tail_solution(primary, partner; pair_fingerprint=GHZ_TAIL_PAIR_FP + 30,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 30)
        mode2 = tail_solution(primary, partner2; pair_fingerprint=GHZ_TAIL_PAIR_FP + 31,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 31)
        @test x_nmb(mode1.pair, 3.0) != x_nmb(mode2.pair, 3.0)
        @test x_nn(mode1, 3.0) != x_nn(mode2, 3.0)

        mode = manufactured_tail()
        components = corrector_tetrad_components(mode, 3.0)
        @test components isa IRGCorrectorTetradComponents{2,Float64}
        assembled = (components.A, components.B, components.C, components.N)
        direct = (x_mmb(mode.pair.primary.head, 3.0), x_nm(mode.pair.primary, 3.0),
            x_nmb(mode.pair, 3.0), x_nn(mode, 3.0))
        @test all(isapprox.(assembled, direct; rtol=3e-13, atol=3e-14))
        @test components.policy isa FiniteRayFreeIVP
        @test components.gauge isa GHZIRG
        @test components.phase isa ExpMinusIomegaUPlusImPhi
        psrc, qsrc = mode.pair.primary.head.source, mode.pair.partner.head.source
        @test (components.M, components.l, components.m, components.omega,
            components.rb, components.r) ==
            (psrc.M, psrc.l, psrc.m, psrc.omega, mode.boundary.rb, 3.0)
        @test (components.primary_source_fingerprint,
            components.partner_source_fingerprint, components.pair_fingerprint,
            components.tail_boundary_fingerprint) ==
            (psrc.source_fingerprint, qsrc.source_fingerprint,
             mode.pair.pair_fingerprint, mode.boundary.boundary_fingerprint)
        @test components.normalization isa EinsteinRHS
        @test components.tetrad isa KinnersleyMinusPlus
        @test components.slicing isa SchwarzschildOutgoing
    end

    @testset "scale-safe Green composition and angular precision" begin
        F = floatmax(Float64)
        primary, partner = zero_tail_inputs()
        mode = tail_solution(primary, partner; x0=0.75F,
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 40,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 40)
        @test x_nn(mode, 4.0) == 0.375F + 0im

        T = Float16
        primary, partner = zero_tail_inputs(Val(2), T; l=256, m=0,
            omega=zero(Complex{T}))
        mode16 = tail_solution(primary, partner; x0=zero(T),
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 41,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 41)
        @test x_nn(mode16, T(2)) == zero(Complex{T})

        primary, partner = zero_tail_inputs(Val(2), T; l=2, m=1,
            omega=zero(Complex{T}), primary_kwargs=(; Sln=q -> inv(q)))
        scaled16 = tail_solution(primary, partner; x0=zero(T),
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 43,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 43)
        @test x_nn(scaled16, T(400)) == Complex{T}(T(200), zero(T))

        huge = typemax(Int)
        primary, partner = zero_tail_inputs(Val(2), Float64; l=huge, m=huge,
            omega=zero(ComplexF64))
        pair = tail_pair(primary, partner; pair_fingerprint=GHZ_TAIL_PAIR_FP + 42)
        bc = TailRayIVP(pair; x_nn=0.0, boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 42,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_corrector_tail(pair, bc)
    end

    @testset "source, boundary, radius, tolerance, and raw rejection" begin
        primary, partner = zero_tail_inputs()
        pair = tail_pair(primary, partner)
        for bad in (true, Inf, Complex{Bool}(true, false))
            @test_throws ArgumentError TailRayIVP(pair; x_nn=bad,
                boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP, policy=FiniteRayFreeIVP())
        end
        for fp in (0, true, 2.0)
            @test_throws ArgumentError TailRayIVP(pair; x_nn=0,
                boundary_fingerprint=fp, policy=FiniteRayFreeIVP())
        end
        @test_throws ArgumentError TailRayIVP(pair; rb=2.5, x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP, policy=FiniteRayFreeIVP())
        for bad_rb in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError TailRayIVP(pair; rb=bad_rb, x_nn=0,
                boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP, policy=FiniteRayFreeIVP())
        end

        bc = TailRayIVP(pair; x_nn=0.2 - 0.1im,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP, policy=FiniteRayFreeIVP())
        mode = solve_irg_corrector_tail(pair, bc)
        for bad_r in (true, Complex{Bool}(true, false), Inf, 1.9)
            @test_throws ArgumentError x_nn(mode, bad_r)
            @test_throws ArgumentError dr_x_nn(mode, bad_r)
        end
        for kwargs in ((; atol=true), (; rtol=true), (; atol=NaN), (; rtol=Inf),
                (; atol=0.0, rtol=0.0))
            @test_throws ArgumentError x_nn(mode, 2.0; kwargs...)
            @test_throws ArgumentError dr_x_nn(mode, 2.0; kwargs...)
            @test_throws ArgumentError corrector_tetrad_components(mode, 2.0; kwargs...)
        end
        @test_throws ArgumentError irg_tail_residual(true, 0, 0, 0)
        @test_throws ArgumentError irg_tail_residual(2.0, true, 0, 0)
        @test_throws ArgumentError irg_tail_residual(2.0, 0, Inf, 0)

        invalid_boundaries = (
            raw_tail_boundary(bc; l=3),
            raw_tail_boundary(bc; pair_fingerprint=bc.pair_fingerprint + 1),
            raw_tail_boundary(bc;
                primary_source_fingerprint=bc.primary_source_fingerprint + 1),
            raw_tail_boundary(bc; primary_middle_x0=bc.primary_middle_x0 + 1),
            raw_tail_boundary(bc; rb=2.5),
            raw_tail_boundary(bc; boundary_fingerprint=0),
        )
        for invalid in invalid_boundaries
            @test_throws ArgumentError solve_irg_corrector_tail(pair, invalid)
            raw = IRGCorrectorTailMode(pair, invalid, GHZIRG())
            @test_throws ArgumentError x_nn(raw, invalid.rb)
            @test_throws ArgumentError dr_x_nn(raw, invalid.rb)
            @test_throws ArgumentError corrector_tetrad_components(raw, invalid.rb)
        end

        bad_primary, bad_partner = zero_tail_inputs(primary_kwargs=(; Sln=r -> true))
        bad_pair = tail_pair(bad_primary, bad_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 50)
        bad_bc = TailRayIVP(bad_pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 50, policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_corrector_tail(bad_pair, bad_bc)
        raw_bad = IRGCorrectorTailMode(bad_pair, bad_bc, GHZIRG())
        @test_throws ArgumentError x_nn(raw_bad, 2.0)
        @test_throws ArgumentError corrector_tetrad_components(raw_bad, 2.0)

        inf_primary, inf_partner = zero_tail_inputs(primary_kwargs=(; Sln=r -> Inf))
        inf_pair = tail_pair(inf_primary, inf_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 53)
        inf_bc = TailRayIVP(inf_pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 53, policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_corrector_tail(inf_pair, inf_bc)

        late_primary, late_partner = zero_tail_inputs(primary_kwargs=(;
            Sln=r -> r > 2.2 ? Inf : 0.0))
        late = tail_solution(late_primary, late_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 51,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 51)
        @test_throws ArgumentError x_nn(late, 3.0)
        @test_throws ArgumentError dr_x_nn(late, 3.0)

        late_bool_primary, late_bool_partner = zero_tail_inputs(primary_kwargs=(;
            Sln=r -> r > 2.2 ? true : 0.0))
        late_bool = tail_solution(late_bool_primary, late_bool_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 54,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 54)
        @test_throws ArgumentError x_nn(late_bool, 3.0)

        poison(q) = error("partner S_ln must remain inert")
        inert_primary, inert_partner = zero_tail_inputs(partner_kwargs=(; Sln=poison))
        inert = tail_solution(inert_primary, inert_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 52,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 52)
        @test x_nn(inert, 3.0) == 0.0 + 0.0im

        overflow_primary, overflow_partner = zero_tail_inputs(
            omega=1.0 + 0.0im, primary_kwargs=(; head_x0=floatmax(Float64)))
        overflow_pair = tail_pair(overflow_primary, overflow_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 55)
        overflow_bc = TailRayIVP(overflow_pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 55,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_corrector_tail(overflow_pair, overflow_bc)
        overflow_raw = IRGCorrectorTailMode(overflow_pair, overflow_bc, GHZIRG())
        @test_throws ArgumentError x_nn(overflow_raw, 2.0)
        @test_throws ArgumentError corrector_tetrad_components(overflow_raw, 2.0)

        order_primary, order_partner = zero_tail_inputs(Val(1))
        order_pair = tail_pair(order_primary, order_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 56)
        order_bc = TailRayIVP(order_pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 56,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_corrector_tail(pair, order_bc)
        big_primary, big_partner = zero_tail_inputs(Val(2), BigFloat)
        big_pair = tail_pair(big_primary, big_partner;
            pair_fingerprint=GHZ_TAIL_PAIR_FP + 57)
        big_bc = TailRayIVP(big_pair; x_nn=0,
            boundary_fingerprint=GHZ_TAIL_BOUNDARY_FP + 57,
            policy=FiniteRayFreeIVP())
        @test_throws ArgumentError solve_irg_corrector_tail(pair, big_bc)
    end

    @testset "dispatch has no internal ambiguities" begin
        for callable in (TailRayIVP, solve_irg_corrector_tail, x_nn, dr_x_nn,
                corrector_tetrad_components)
            ms = collect(methods(callable))
            @test all(!Base.isambiguous(ms[i], ms[j])
                for i in eachindex(ms) for j in (i + 1):length(ms))
        end
    end
end

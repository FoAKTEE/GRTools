using Test
using MetricReconstruction.GHZ

const GHZ_PAIR_PRIMARY_FP = UInt128(0x1d82f9d1457b49e6a2b7185f44f02131)
const GHZ_PAIR_PARTNER_FP = UInt128(0x8a43cf2251ac4cedb817fe3d5f922e47)
const GHZ_PAIR_FP = UInt128(0xf32e6b1ce94b47d7af2d08c7319185a4)

function corrector_middle_mode(::Val{N}, ::Type{T}=Float64; M=one(T), l=2, m=1,
        omega=Complex{T}(T(0.2), T(-0.1)), rb=T(2) * M,
        source_fingerprint=GHZ_PAIR_PRIMARY_FP, head_x0=zero(Complex{T}),
        head_dx0=zero(Complex{T}), middle_x0=zero(Complex{T}),
        middle_dx0=zero(Complex{T}), Sll=nothing, Slm=nothing, Sln=nothing,
        phase=ExpMinusIomegaUPlusImPhi()) where {N,T<:AbstractFloat}
    sll = isnothing(Sll) ? (q -> zero(q)) : Sll
    slm = isnothing(Slm) ? (q -> zero(q)) : Slm
    common = (; M, l, m, omega, source_fingerprint, phase)
    src = isnothing(Sln) ? TetradSourceMode{N}(sll; Slm=slm, common...) :
                           TetradSourceMode{N}(sll; Slm=slm, Sln, common...)
    head = solve_irg_corrector_head(src,
        InnerRayIVP(src; rb, x_mmb=head_x0, dr_x_mmb=head_dx0))
    solve_irg_corrector_middle(head,
        MiddleRayIVP(head; x_nm=middle_x0, dr_x_nm=middle_dx0))
end

function raw_corrector_sources(mode::IRGCorrectorMiddleMode{N,T};
        Sll=mode.head.source.Sll, Slm=mode.head.source.Slm,
        Sln=mode.head.source.Sln) where {N,T}
    src = mode.head.source
    raw_src = TetradSourceMode{N,T,typeof(src.normalization),typeof(Sll),typeof(Slm),
        typeof(Sln)}(Sll, Slm, Sln, src.M, src.l, src.m, src.omega,
        src.source_fingerprint, src.normalization, src.tetrad, src.slicing, src.phase)
    raw_head = IRGCorrectorHeadMode(raw_src, mode.head.boundary, GHZIRG())
    IRGCorrectorMiddleMode(raw_head, mode.boundary, GHZIRG())
end

@testset "GHZ — physical corrector middle pair prerequisites" begin
    @testset "phase pin and inert typed S_ln" begin
        mode = corrector_middle_mode(Val(2))
        @test mode.head.source.phase isa ExpMinusIomegaUPlusImPhi
        explicit = corrector_middle_mode(Val(2); phase=ExpMinusIomegaUPlusImPhi(),
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 1)
        @test explicit.head.source.phase isa ExpMinusIomegaUPlusImPhi
        @test_throws ArgumentError corrector_middle_mode(Val(2); phase=:wrong,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 2)

        setprecision(BigFloat, 256) do
            default_sln = corrector_middle_mode(Val(2), BigFloat;
                source_fingerprint=GHZ_PAIR_PRIMARY_FP + 3)
            value = default_sln.head.source.Sln(big"3.25")
            @test value == 0
            @test value isa BigFloat

            sln(q) = Complex{BigFloat}(q, -q / 2)
            custom_sln = corrector_middle_mode(Val(2), BigFloat; Sln=sln,
                source_fingerprint=GHZ_PAIR_PRIMARY_FP + 4)
            @test custom_sln.head.source.Sln === sln
            @test custom_sln.head.source.Sln(big"2.0") ==
                Complex{BigFloat}(big"2.0", big"-1.0")
        end

        poison_sln(q) = error("S_ln is not part of the head or middle equations")
        inert = corrector_middle_mode(Val(2); Sln=poison_sln,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 5)
        @test x_mmb(inert.head, 3.0) == 0.0 + 0.0im
        @test x_nm(inert, 3.0) == 0.0 + 0.0im
    end

    @testset "explicit partner parity and provenance" begin
        B, dB = 1.25 - 0.75im, -0.4 + 0.6im
        for (offset, m) in enumerate((0, 1, -1, 2, -2))
            l = max(2, abs(m))
            omega = 0.17 - 0.03im
            primary = corrector_middle_mode(Val(2); l, m, omega,
                middle_x0=10 + offset * im, middle_dx0=-7 + offset * im,
                source_fingerprint=GHZ_PAIR_PRIMARY_FP + offset)
            partner = corrector_middle_mode(Val(2); l, m=-m, omega=-conj(omega),
                middle_x0=B, middle_dx0=dB,
                source_fingerprint=GHZ_PAIR_PARTNER_FP + offset)
            pair = PhysicalCorrectorMiddlePair(primary, partner;
                pair_fingerprint=GHZ_PAIR_FP + offset)
            parity = isodd(m) ? 1 : -1
            @test perturbative_order(pair) == 2
            @test x_nmb(pair, 2.0) == parity * conj(B)
            @test dr_x_nmb(pair, 2.0) == parity * conj(dB)
            @test x_nmb(pair, 3.0) ≈ parity * conj(x_nm(partner, 3.0))
            @test dr_x_nmb(pair, 3.0) ≈ parity * conj(dr_x_nm(partner, 3.0))
        end

        omega = 0.21 - 0.08im
        primary = corrector_middle_mode(Val(2); m=1, omega, middle_x0=99 - 2im,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 20,
            Sll=q -> zero(q), Slm=q -> zero(q))
        partner = corrector_middle_mode(Val(2); m=-1, omega=-conj(omega),
            middle_x0=B, middle_dx0=dB, source_fingerprint=GHZ_PAIR_PARTNER_FP + 20,
            Sll=q -> zero(q) + zero(q), Slm=q -> zero(q) - zero(q))
        alternate_primary = corrector_middle_mode(Val(2); m=1, omega,
            middle_x0=-300 + 40im, middle_dx0=81 - 9im,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 21)
        pair = PhysicalCorrectorMiddlePair(primary, partner;
            pair_fingerprint=GHZ_PAIR_FP + 20)
        alternate_pair = PhysicalCorrectorMiddlePair(alternate_primary, partner;
            pair_fingerprint=GHZ_PAIR_FP + 21)
        @test x_nmb(pair, 3.0) == x_nmb(alternate_pair, 3.0)
        @test dr_x_nmb(pair, 3.0) == dr_x_nmb(alternate_pair, 3.0)
        reverse_pair = PhysicalCorrectorMiddlePair(partner, primary;
            pair_fingerprint=GHZ_PAIR_FP + 22)
        @test x_nmb(reverse_pair, 2.0) == conj(primary.boundary.x0)
    end

    @testset "BigFloat and overflow-safe extreme parity" begin
        setprecision(BigFloat, 256) do
            omega = Complex{BigFloat}(big"0.125", big"-0.0625")
            B = Complex{BigFloat}(big"1.75", big"-0.375")
            dB = Complex{BigFloat}(big"-0.25", big"0.875")
            primary = corrector_middle_mode(Val(2), BigFloat; m=1, omega,
                source_fingerprint=GHZ_PAIR_PRIMARY_FP + 30)
            partner = corrector_middle_mode(Val(2), BigFloat; m=-1,
                omega=-conj(omega), middle_x0=B, middle_dx0=dB,
                source_fingerprint=GHZ_PAIR_PARTNER_FP + 30)
            pair = PhysicalCorrectorMiddlePair(primary, partner;
                pair_fingerprint=GHZ_PAIR_FP + 30)
            @test x_nmb(pair, big"2.0") == conj(B)
            @test dr_x_nmb(pair, big"2.0") == conj(dB)
            @test x_nmb(pair, big"3.0") isa Complex{BigFloat}
            @test dr_x_nmb(pair, big"3.0") isa Complex{BigFloat}
        end

        m = typemax(Int)
        B = 0.75 - 0.5im
        primary = corrector_middle_mode(Val(2); l=m, m, omega=0.0 + 0.0im,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 31)
        partner = corrector_middle_mode(Val(2); l=m, m=-m, omega=-0.0 + 0.0im,
            middle_x0=B, source_fingerprint=GHZ_PAIR_PARTNER_FP + 31)
        pair = PhysicalCorrectorMiddlePair(primary, partner;
            pair_fingerprint=GHZ_PAIR_FP + 31)
        @test x_nmb(pair, 2.0) == conj(B)
    end

    @testset "key, precision, fingerprint, and explicit-partner rejection" begin
        omega = 0.2 - 0.1im
        primary = corrector_middle_mode(Val(2); m=1, omega,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 40)
        partner = corrector_middle_mode(Val(2); m=-1, omega=-conj(omega),
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 40)
        @test PhysicalCorrectorMiddlePair(primary, partner;
            pair_fingerprint=GHZ_PAIR_FP + 40) isa PhysicalCorrectorMiddlePair

        wrong_m = corrector_middle_mode(Val(2); m=0, omega=-conj(omega),
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 41)
        wrong_omega = corrector_middle_mode(Val(2); m=-1, omega=0.3 + 0.2im,
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 42)
        wrong_M = corrector_middle_mode(Val(2); M=2.0, rb=4.0, m=-1,
            omega=-conj(omega), source_fingerprint=GHZ_PAIR_PARTNER_FP + 43)
        wrong_l = corrector_middle_mode(Val(2); l=3, m=-1, omega=-conj(omega),
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 44)
        wrong_rb = corrector_middle_mode(Val(2); rb=2.5, m=-1,
            omega=-conj(omega), source_fingerprint=GHZ_PAIR_PARTNER_FP + 45)
        wrong_order = corrector_middle_mode(Val(1); m=-1, omega=-conj(omega),
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 46)
        wrong_precision = corrector_middle_mode(Val(2), BigFloat; m=-1,
            omega=Complex{BigFloat}(-conj(omega)),
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 47)
        for invalid in (wrong_m, wrong_omega, wrong_M, wrong_l, wrong_rb,
                wrong_order, wrong_precision)
            @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, invalid;
                pair_fingerprint=GHZ_PAIR_FP + 41)
        end

        @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, partner)
        for invalid_fp in (0, true, 2.0)
            @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, partner;
                pair_fingerprint=invalid_fp)
        end
        @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary;
            pair_fingerprint=GHZ_PAIR_FP + 42)
        @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, nothing;
            pair_fingerprint=GHZ_PAIR_FP + 42)

        static_mode = corrector_middle_mode(Val(2); m=0, omega=0.0 + 0.0im,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 48)
        @test_throws ArgumentError PhysicalCorrectorMiddlePair(static_mode;
            pair_fingerprint=GHZ_PAIR_FP + 48)
    end

    @testset "raw constituent and pair revalidation" begin
        omega = 0.2 - 0.1im
        primary = corrector_middle_mode(Val(2); m=1, omega,
            source_fingerprint=GHZ_PAIR_PRIMARY_FP + 50)
        partner = corrector_middle_mode(Val(2); m=-1, omega=-conj(omega),
            source_fingerprint=GHZ_PAIR_PARTNER_FP + 50)
        bc = partner.boundary
        wrong_key_bc = MiddleRayIVP{2,Float64}(bc.M, bc.l, bc.m, bc.omega,
            bc.source_fingerprint + 1, bc.rb, bc.x0, bc.dx0)
        wrong_key_mode = IRGCorrectorMiddleMode(partner.head, wrong_key_bc, GHZIRG())
        @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, wrong_key_mode;
            pair_fingerprint=GHZ_PAIR_FP + 50)

        wrong_rb_bc = MiddleRayIVP{2,Float64}(bc.M, bc.l, bc.m, bc.omega,
            bc.source_fingerprint, 2.5, bc.x0, bc.dx0)
        wrong_rb_mode = IRGCorrectorMiddleMode(partner.head, wrong_rb_bc, GHZIRG())
        P, Q = typeof(primary), typeof(wrong_rb_mode)
        raw_constituent_pair = PhysicalCorrectorMiddlePair{2,Float64,P,Q}(
            primary, wrong_rb_mode, GHZ_PAIR_FP + 51, Val(:validated))
        @test_throws ArgumentError x_nmb(raw_constituent_pair, 2.5)
        @test_throws ArgumentError dr_x_nmb(raw_constituent_pair, 2.5)

        P, Q = typeof(primary), typeof(partner)
        raw_zero_pair = PhysicalCorrectorMiddlePair{2,Float64,P,Q}(
            primary, partner, UInt128(0), Val(:validated))
        @test_throws ArgumentError x_nmb(raw_zero_pair, 2.0)
        @test_throws ArgumentError dr_x_nmb(raw_zero_pair, 2.0)

        bool_src = raw_corrector_sources(partner; Sll=q -> true)
        @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, bool_src;
            pair_fingerprint=GHZ_PAIR_FP + 52)
        nonfinite_src = raw_corrector_sources(partner; Slm=q -> Inf)
        @test_throws ArgumentError PhysicalCorrectorMiddlePair(primary, nonfinite_src;
            pair_fingerprint=GHZ_PAIR_FP + 53)
    end
end

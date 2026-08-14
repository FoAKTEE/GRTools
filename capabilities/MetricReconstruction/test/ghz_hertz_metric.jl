using Test
using MetricReconstruction.GHZ: ConjugateHertzRadialJet, PhysicalHertzPair,
    IRGHertzMetricMode, QuarterThorn4, HertzIRG, GHZIRG, KinnersleyMinusPlus,
    SchwarzschildOutgoing, ExpMinusIomegaUPlusImPhi,
    FixedOutgoingURadialDerivative, FiniteHertzRadialIVP,
    LinearizedPsi0, ReducedWeylMode, HertzRadialBoundary,
    invert_hertz, conjugate_hertz_radial_jet, perturbative_order, physical_hertz_jet,
    hertz_metric, static_hertz_metric

const GHZ_METRIC_SFP = UInt128(0xc67aa38dbe7d4f7188e513509abc6a10)
const GHZ_METRIC_BFP = UInt128(0xeb3870c691044e27a882fe4658dfa9e6)
const GHZ_METRIC_PFP = UInt128(0x246f111a0e6e44aa94d695320034d865)

function static_bar_jet(::Val{N}, ::Type{T}=Float64; M=one(T), l=2, m=0,
        omega=zero(Complex{T}), r=T(4), values=nothing, sfp=GHZ_METRIC_SFP,
        bfp=GHZ_METRIC_BFP, phase=ExpMinusIomegaUPlusImPhi(),
        radial_derivative=FixedOutgoingURadialDerivative(),
        policy=FiniteHertzRadialIVP(), kwargs...) where {N,T<:AbstractFloat}
    data = isnothing(values) ?
        (r^4, T(4) * r^3, T(12) * r^2, T(24) * r, T(24)) : values
    ConjugateHertzRadialJet{N}(data; M=T(M), l, m, omega=Complex{T}(omega), r,
        source_fingerprint=sfp, boundary_fingerprint=bfp, phase, radial_derivative,
        policy, kwargs...)
end

function paired_static_jets(::Val{N}, ::Type{T}=Float64; l=2, m=0, r=T(4),
        omega=zero(Complex{T}), values=nothing, partner_values=values,
        sfp=GHZ_METRIC_SFP, bfp=GHZ_METRIC_BFP,
        partner_sfp=GHZ_METRIC_SFP + 1, partner_bfp=GHZ_METRIC_BFP + 1,
        kwargs...) where {N,T<:AbstractFloat}
    omegaT = Complex{T}(omega)
    mode = static_bar_jet(Val(N), T; l, m, omega=omegaT, r, values, sfp, bfp,
        kwargs...)
    partner = static_bar_jet(Val(N), T; l, m=-m, r, values=partner_values,
        omega=-conj(omegaT), sfp=partner_sfp, bfp=partner_bfp, kwargs...)
    mode, partner
end

hertz_coefficients(metric) = (metric.k_nn, metric.k_nbar, metric.k_nm,
    metric.k_barbar, metric.k_mm)

function same_hertz_metric(left, right)
    typeof(left) === typeof(right) && all(name ->
        isequal(getfield(left, name), getfield(right, name)), fieldnames(typeof(left)))
end

@testset "GHZ — solver-independent outgoing-u IRG Hertz metric injection" begin
    @testset "exact l=2, m=0 r^4 oracle and convention factors" begin
        for N in (1, 2)
            r = 5.0
            mode, partner = paired_static_jets(Val(N); r)
            pair = PhysicalHertzPair(mode, partner; pair_fingerprint=GHZ_METRIC_PFP)
            metric = static_hertz_metric(pair)
            @test same_hertz_metric(metric, hertz_metric(pair))
            want = (-sqrt(6.0) * r^2, sqrt(2.0) * r^2, -sqrt(2.0) * r^2,
                -2r^2, -2r^2)

            @test metric isa IRGHertzMetricMode{N,Float64}
            @test perturbative_order(mode) == N
            @test perturbative_order(pair) == N
            @test perturbative_order(metric) == N
            got = (metric.k_nn, metric.k_nbar, metric.k_nm,
                metric.k_barbar, metric.k_mm)
            @test all(isapprox.(got, want; rtol=5e-15))
            @test metric.k_nn != sqrt(6.0) * r^2       # mostly-plus sign
            @test metric.k_nn != -2sqrt(6.0) * r^2    # real-field factor is not doubled
            @test metric.k_nbar == -metric.k_nm
            @test metric.gauge isa HertzIRG
            @test metric.tetrad isa KinnersleyMinusPlus
            @test metric.slicing isa SchwarzschildOutgoing
            @test metric.normalization isa QuarterThorn4
        end
    end

    @testset "adapter from a typed conjugate-Hertz solve" begin
        source(m, off; omega=0.0) = ReducedWeylMode{2}(r -> zero(ComplexF64);
            M=1.0, l=2, m, omega, source_fingerprint=GHZ_METRIC_SFP + off,
            content=LinearizedPsi0(), phase=ExpMinusIomegaUPlusImPhi())
        solve(m, off; omega=0.0) = begin
            src = source(m, off; omega)
            bc = HertzRadialBoundary(src; rb=4.0, b0=4.0^4, b1=4 * 4.0^3,
                b2=12 * 4.0^2, b3=24 * 4.0,
                boundary_fingerprint=GHZ_METRIC_BFP + off,
                radial_derivative=FixedOutgoingURadialDerivative(),
                policy=FiniteHertzRadialIVP())
            invert_hertz(src, bc)
        end
        solved_mode, solved_partner = solve(0, 0), solve(0, 1)
        mode = conjugate_hertz_radial_jet(solved_mode, 4.0)
        partner = conjugate_hertz_radial_jet(solved_partner, 4.0)
        pair = PhysicalHertzPair(mode, partner; pair_fingerprint=GHZ_METRIC_PFP)
        metric = static_hertz_metric(pair)
        @test same_hertz_metric(metric, hertz_metric(pair))
        @test mode isa ConjugateHertzRadialJet{2,Float64}
        @test mode.derivatives == (256, 256, 192, 96, 0)
        @test mode.phase === solved_mode.phase
        @test mode.radial_derivative === solved_mode.radial_derivative
        @test mode.policy === solved_mode.policy
        @test metric.k_nn ≈ -sqrt(6.0) * 16
        @test metric.k_nbar ≈ sqrt(2.0) * 16
        @test metric.k_nm ≈ -sqrt(2.0) * 16
        @test metric.mode_source_fingerprint == GHZ_METRIC_SFP
        @test metric.partner_source_fingerprint == GHZ_METRIC_SFP + 1
        @test metric.mode_boundary_fingerprint == GHZ_METRIC_BFP
        @test metric.partner_boundary_fingerprint == GHZ_METRIC_BFP + 1
        @test metric.pair_fingerprint == GHZ_METRIC_PFP

        dynamic_solved = solve(1, 20; omega=0.2 + 0.1im)
        dynamic_jet = conjugate_hertz_radial_jet(dynamic_solved, 4.0)
        @test dynamic_jet.omega == 0.2 + 0.1im
        @test dynamic_jet.phase === dynamic_solved.source.phase
        @test dynamic_jet.phase === dynamic_solved.boundary.phase
        @test dynamic_jet.radial_derivative === dynamic_solved.boundary.radial_derivative
        @test dynamic_jet.policy === dynamic_solved.boundary.policy
    end

    @testset "explicit caller-asserted partner and +/-m phase" begin
        r = 4.0
        B = (2 + 3im, -1 + 4im, 5 - 2im, 7im, -3 + im)
        Q = (-4 + 2im, 3 - im, 6 + 5im, -2im, 8 - im)
        mode, partner = paired_static_jets(Val(2); l=3, m=1, r, values=B,
            partner_values=Q)
        pair = PhysicalHertzPair(mode, partner; pair_fingerprint=GHZ_METRIC_PFP)
        P = ntuple(i -> -conj(Q[i]), 5)
        @test physical_hertz_jet(pair) == P
        @test physical_hertz_jet(pair) != ntuple(i -> -conj(B[i]), 5)

        lambda2 = sqrt((3 - 1) * 3 * (3 + 1) * (3 + 2))
        c = sqrt((3 - 1) * (3 + 2)) / (2sqrt(2.0))
        metric = static_hertz_metric(pair)
        @test same_hertz_metric(metric, hertz_metric(pair))
        @test metric.k_nn ≈ -lambda2 * (P[1] + B[1]) / (4r^2)
        @test metric.k_nbar ≈ c * (P[2] / r - 2P[1] / r^2)
        @test metric.k_nm ≈ -c * (B[2] / r - 2B[1] / r^2)
        @test metric.k_barbar ≈ -(P[3] - 2P[2] / r) / 2
        @test metric.k_mm ≈ -(B[3] - 2B[2] / r) / 2

        reverse = PhysicalHertzPair(partner, mode; pair_fingerprint=GHZ_METRIC_PFP + 1)
        @test physical_hertz_jet(reverse) == ntuple(i -> -conj(B[i]), 5)
        @test_throws ArgumentError PhysicalHertzPair(mode)
        @test_throws ArgumentError PhysicalHertzPair(mode, nothing)
    end

    @testset "finite-frequency fixed-outgoing-u kernel" begin
        r = 4.75
        B = (2 + 3im, -1 + 4im, 5 - 2im, 7im, -3 + im)
        Q = (-4 + 2im, 3 - im, 6 + 5im, -2im, 8 - im)
        static_mode, static_partner = paired_static_jets(Val(2); l=3, m=1, r,
            values=B, partner_values=Q)
        static_pair = PhysicalHertzPair(static_mode, static_partner;
            pair_fingerprint=GHZ_METRIC_PFP)
        static_metric = hertz_metric(static_pair)

        for (j, omega) in enumerate((0.25 + 0.0im, 0.25 + 0.125im))
            mode, partner = paired_static_jets(Val(2); l=3, m=1, omega, r,
                values=B, partner_values=Q,
                sfp=GHZ_METRIC_SFP + 10j, bfp=GHZ_METRIC_BFP + 10j,
                partner_sfp=GHZ_METRIC_SFP + 10j + 1,
                partner_bfp=GHZ_METRIC_BFP + 10j + 1)
            pair = PhysicalHertzPair(mode, partner;
                pair_fingerprint=GHZ_METRIC_PFP + j)
            metric = hertz_metric(pair)

            @test mode.omega == omega
            @test partner.omega == -conj(omega)
            @test metric.omega == omega
            @test metric.omega != static_metric.omega
            @test hertz_coefficients(metric) == hertz_coefficients(static_metric)
            @test metric.mode_source_fingerprint != static_metric.mode_source_fingerprint
            @test metric.phase isa ExpMinusIomegaUPlusImPhi
            @test metric.radial_derivative isa FixedOutgoingURadialDerivative
            @test metric.policy isa FiniteHertzRadialIVP
            @test metric.normalization isa QuarterThorn4
            @test_throws ArgumentError static_hertz_metric(pair)
            if !iszero(imag(omega))
                wrong_partner = static_bar_jet(Val(2); l=3, m=-1, omega=-omega, r,
                    values=Q, sfp=GHZ_METRIC_SFP + 100,
                    bfp=GHZ_METRIC_BFP + 100)
                @test_throws ArgumentError PhysicalHertzPair(mode, wrong_partner;
                    pair_fingerprint=GHZ_METRIC_PFP + 100)
            end
        end
    end

    @testset "independent fixed-t to fixed-u derivative oracle" begin
        # At M=1, r=4, f=1/2, f'=1/8, and omega=1/4, the outgoing-u radial
        # operator on a BL-separated amplitude is L=d/dr-i*omega/f.  For a locally
        # constant BL amplitude, L(1)=-i/2 and L^2(1)=-1/4+i/8.
        M, r, omega = 1.0, 4.0, 0.25 + 0.0im
        f = 1 - 2M / r
        df = 2M / r^2
        fixed_t = (1.0 + 0.0im, 0.0im, 0.0im, 0.0im, 0.0im)
        alpha = -im * omega / f
        dalpha = im * omega * df / f^2
        fixed_u = (fixed_t[1], alpha, alpha^2 + dalpha, 0.0im, 0.0im)
        @test fixed_u[2] == -0.5im
        @test fixed_u[3] == -0.25 + 0.125im

        mode, partner = paired_static_jets(Val(2); M, l=2, m=1, r, omega,
            values=fixed_u, partner_values=ntuple(_ -> 0.0im, 5),
            sfp=GHZ_METRIC_SFP + 40, bfp=GHZ_METRIC_BFP + 40,
            partner_sfp=GHZ_METRIC_SFP + 41, partner_bfp=GHZ_METRIC_BFP + 41)
        metric = hertz_metric(PhysicalHertzPair(mode, partner;
            pair_fingerprint=GHZ_METRIC_PFP + 40))
        expected = (-sqrt(24.0) / 64, 0.0im, sqrt(2.0) * (1 + im) / 16,
            0.0im, 1 / 8 - 3im / 16)
        @test all(isapprox.(hertz_coefficients(metric), expected; rtol=4e-15, atol=2e-16))

        naive_mode, naive_partner = paired_static_jets(Val(2); M, l=2, m=1, r, omega,
            values=fixed_t, partner_values=ntuple(_ -> 0.0im, 5),
            sfp=GHZ_METRIC_SFP + 42, bfp=GHZ_METRIC_BFP + 42,
            partner_sfp=GHZ_METRIC_SFP + 43, partner_bfp=GHZ_METRIC_BFP + 43)
        naive = hertz_metric(PhysicalHertzPair(naive_mode, naive_partner;
            pair_fingerprint=GHZ_METRIC_PFP + 41))
        @test metric.k_nm != naive.k_nm
        @test metric.k_mm != naive.k_mm
    end

    @testset "BigFloat precision" begin
        setprecision(BigFloat, 256) do
            r = big"7.25"
            mode, partner = paired_static_jets(Val(2), BigFloat; r)
            metric = static_hertz_metric(PhysicalHertzPair(mode, partner;
                pair_fingerprint=GHZ_METRIC_PFP))
            @test same_hertz_metric(metric, hertz_metric(PhysicalHertzPair(mode, partner;
                pair_fingerprint=GHZ_METRIC_PFP)))
            @test metric isa IRGHertzMetricMode{2,BigFloat}
            @test all(z -> z isa Complex{BigFloat},
                (metric.k_nn, metric.k_nbar, metric.k_nm,
                    metric.k_barbar, metric.k_mm))
            @test abs(metric.k_nn + sqrt(big"6") * r^2) < big"1e-70"
            @test abs(metric.k_nbar - sqrt(big"2") * r^2) < big"1e-70"
            @test abs(metric.k_nm + sqrt(big"2") * r^2) < big"1e-70"
            @test metric.k_barbar == -2r^2
            @test metric.k_mm == -2r^2

            omega = Complex{BigFloat}(big"0.125", big"0.0625")
            dynamic_mode, dynamic_partner = paired_static_jets(Val(2), BigFloat;
                r, omega)
            dynamic_pair = PhysicalHertzPair(dynamic_mode, dynamic_partner;
                pair_fingerprint=GHZ_METRIC_PFP + 1)
            dynamic_metric = hertz_metric(dynamic_pair)
            @test dynamic_metric isa IRGHertzMetricMode{2,BigFloat}
            @test dynamic_metric.omega == omega
            @test hertz_coefficients(dynamic_metric) == hertz_coefficients(metric)
            @test_throws ArgumentError static_hertz_metric(dynamic_pair)
        end
    end

    @testset "partner key and metadata rejection" begin
        mode, partner = paired_static_jets(Val(2); l=3, m=1)
        @test PhysicalHertzPair(mode, partner; pair_fingerprint=GHZ_METRIC_PFP) isa
            PhysicalHertzPair{2,Float64}

        bad(; kwargs...) = static_bar_jet(Val(2); l=3, m=-1, kwargs...)
        pair(candidate) = PhysicalHertzPair(mode, candidate;
            pair_fingerprint=GHZ_METRIC_PFP)
        @test_throws ArgumentError pair(bad(l=4, m=-1))
        @test_throws ArgumentError pair(bad(m=0))
        @test_throws ArgumentError pair(bad(omega=0.2 + 0im))
        @test_throws ArgumentError pair(bad(M=2.0))
        @test_throws ArgumentError pair(bad(r=5.0))
        @test_throws ArgumentError PhysicalHertzPair(mode,
            static_bar_jet(Val(1); l=3, m=-1); pair_fingerprint=GHZ_METRIC_PFP)
        @test_throws ArgumentError PhysicalHertzPair(mode,
            static_bar_jet(Val(2), Float32; l=3, m=-1);
            pair_fingerprint=GHZ_METRIC_PFP)

        distinct = bad(sfp=GHZ_METRIC_SFP + 91, bfp=GHZ_METRIC_BFP + 92)
        distinct_pair = pair(distinct)
        distinct_metric = static_hertz_metric(distinct_pair)
        @test same_hertz_metric(distinct_metric, hertz_metric(distinct_pair))
        @test distinct_pair.partner.source_fingerprint == GHZ_METRIC_SFP + 91
        @test distinct_pair.partner.boundary_fingerprint == GHZ_METRIC_BFP + 92
        @test distinct_metric.mode_source_fingerprint == mode.source_fingerprint
        @test distinct_metric.partner_source_fingerprint == distinct.source_fingerprint
        @test distinct_metric.mode_boundary_fingerprint == mode.boundary_fingerprint
        @test distinct_metric.partner_boundary_fingerprint == distinct.boundary_fingerprint
        @test distinct_metric.pair_fingerprint == GHZ_METRIC_PFP
        @test_throws ArgumentError PhysicalHertzPair(mode, partner)
        @test_throws ArgumentError PhysicalHertzPair(mode, partner; pair_fingerprint=0)

        forged_zero_fp = PhysicalHertzPair{2,Float64}(
            mode, partner, UInt128(0), Val(:validated))
        wrong_key_partner = static_bar_jet(Val(2); l=3, m=1,
            sfp=GHZ_METRIC_SFP + 101, bfp=GHZ_METRIC_BFP + 101)
        forged_wrong_key = PhysicalHertzPair{2,Float64}(
            mode, wrong_key_partner, GHZ_METRIC_PFP, Val(:validated))
        for forged in (forged_zero_fp, forged_wrong_key)
            @test_throws ArgumentError physical_hertz_jet(forged)
            @test_throws ArgumentError hertz_metric(forged)
            @test_throws ArgumentError static_hertz_metric(forged)
        end


        forged_mode = ConjugateHertzRadialJet{2,Float64}(
            mode.M, mode.l, mode.m, mode.omega, UInt128(0), mode.boundary_fingerprint,
            mode.r, mode.derivatives, mode.normalization, mode.gauge, mode.tetrad,
            mode.slicing, mode.phase, mode.radial_derivative, mode.policy)
        forged_jet_pair = PhysicalHertzPair{2,Float64}(
            forged_mode, partner, GHZ_METRIC_PFP, Val(:validated))
        @test_throws ArgumentError physical_hertz_jet(forged_jet_pair)
        @test_throws ArgumentError hertz_metric(forged_jet_pair)
        @test_throws ArgumentError static_hertz_metric(forged_jet_pair)

        data = ntuple(_ -> 0.0, 5)
        common = (; M=1.0, l=2, m=0, omega=0.0, r=4.0,
            source_fingerprint=GHZ_METRIC_SFP, boundary_fingerprint=GHZ_METRIC_BFP,
            phase=ExpMinusIomegaUPlusImPhi(),
            radial_derivative=FixedOutgoingURadialDerivative(),
            policy=FiniteHertzRadialIVP())
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            normalization=:wrong)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., gauge=GHZIRG())
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., tetrad=:wrong)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., slicing=:wrong)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., phase=:wrong)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            radial_derivative=:fixed_t)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            policy=:retarded)
        @test_throws UndefKeywordError ConjugateHertzRadialJet{2}(data; M=1.0, l=2,
            m=0, omega=0.0, r=4.0, source_fingerprint=GHZ_METRIC_SFP,
            boundary_fingerprint=GHZ_METRIC_BFP)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data[1:4]; common...)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}((data[1:4]..., Inf); common...)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            source_fingerprint=0)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            boundary_fingerprint=0)
        uint_mode = ConjugateHertzRadialJet{2}(data; common..., l=UInt(2), m=UInt(0))
        @test (uint_mode.l, uint_mode.m) == (2, 0)
        int128_mode = ConjugateHertzRadialJet{2}(data; common...,
            l=Int128(2), m=Int128(-2))
        @test (int128_mode.l, int128_mode.m) == (2, -2)
        bigint_mode = ConjugateHertzRadialJet{2}(data; common...,
            l=BigInt(2), m=BigInt(-2))
        @test (bigint_mode.l, bigint_mode.m) == (2, -2)
        extreme_mode = ConjugateHertzRadialJet{2}(data; common...,
            l=typemax(Int), m=-typemax(Int))
        @test (extreme_mode.l, extreme_mode.m) == (typemax(Int), -typemax(Int))
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., l=true)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., m=true)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            l=2, m=typemin(Int))
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            l=typemax(Int), m=typemin(Int))
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            l=Int128(typemax(Int)) + 1)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
            m=UInt(typemax(Int)) + one(UInt))
        for invalid_fp in (true, 2.0, 2 // 1)
            @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
                source_fingerprint=invalid_fp)
            @test_throws ArgumentError PhysicalHertzPair(mode, partner;
                pair_fingerprint=invalid_fp)
        end
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., M=true)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., r=true)
        @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common..., omega=true)
        for invalid_omega in (Inf, NaN, complex(0.0, Inf), complex(NaN, 0.0))
            @test_throws ArgumentError ConjugateHertzRadialJet{2}(data; common...,
                omega=invalid_omega)
        end
    end

    @testset "static-only boundary" begin
        complex_frequency = 0.0 + 0.2im
        mode = static_bar_jet(Val(2); l=2, m=1, omega=complex_frequency)
        partner = static_bar_jet(Val(2); l=2, m=-1,
            omega=-conj(complex_frequency))
        pair = PhysicalHertzPair(mode, partner; pair_fingerprint=GHZ_METRIC_PFP)
        @test_throws ArgumentError static_hertz_metric(pair)

        mode_real = static_bar_jet(Val(2); l=2, m=1, omega=0.2)
        partner_real = static_bar_jet(Val(2); l=2, m=-1, omega=-0.2)
        @test_throws ArgumentError static_hertz_metric(
            PhysicalHertzPair(mode_real, partner_real; pair_fingerprint=GHZ_METRIC_PFP))
    end

    @testset "subnormal coefficient ordering" begin
        tiny = nextfloat(0.0)
        zero_jet = ntuple(_ -> 0.0 + 0.0im, 5)
        function subnormal_metric(B, Q, offset)
            mode, partner = paired_static_jets(Val(2); M=0.5, l=2, m=0, r=1.0,
                values=B, partner_values=Q,
                sfp=GHZ_METRIC_SFP + offset, bfp=GHZ_METRIC_BFP + offset,
                partner_sfp=GHZ_METRIC_SFP + offset + 1,
                partner_bfp=GHZ_METRIC_BFP + offset + 1)
            hertz_metric(PhysicalHertzPair(mode, partner;
                pair_fingerprint=GHZ_METRIC_PFP + offset))
        end

        p0 = subnormal_metric(zero_jet,
            (tiny + 0.0im, 0.0im, 0.0im, 0.0im, 0.0im), 200)
        p1 = subnormal_metric(zero_jet,
            (0.0im, tiny + 0.0im, 0.0im, 0.0im, 0.0im), 202)
        b1 = subnormal_metric(
            (0.0im, tiny + 0.0im, 0.0im, 0.0im, 0.0im), zero_jet, 204)
        want_nn, want_nbar, want_nm = setprecision(BigFloat, 256) do
            xb = BigFloat(tiny)
            (Float64(-sqrt(BigFloat(24)) * xb / 4),
             Float64(sqrt(BigFloat(2)) * xb / 2),
             Float64(-sqrt(BigFloat(2)) * xb / 2))
        end
        @test tiny / 4 == 0.0
        @test tiny / 2 == 0.0
        @test p0.k_nn == ComplexF64(want_nn, 0.0)
        @test p1.k_nbar == ComplexF64(want_nbar, 0.0)
        @test b1.k_nm == ComplexF64(want_nm, 0.0)
        @test (want_nn, want_nbar, want_nm) == (-tiny, tiny, -tiny)

        cancellation_P0 = nextfloat(floatmin())
        cancellation_P1 = reinterpret(Float64, 9007199254740979)
        cancellation_jet = (cancellation_P0 + 0.0im,
            cancellation_P1 + 0.0im, 0.0im, 0.0im, 0.0im)
        cancellation_nbar = subnormal_metric(zero_jet, cancellation_jet, 206)
        cancellation_nm = subnormal_metric(cancellation_jet, zero_jet, 208)
        want_cancellation = setprecision(BigFloat, 8192) do
            Float64((sqrt(BigFloat(2)) / 2) *
                (BigFloat(cancellation_P1) - 2BigFloat(cancellation_P0)))
        end
        @test want_cancellation == -11tiny
        @test cancellation_nbar.k_nbar == ComplexF64(want_cancellation, 0.0)
        @test cancellation_nm.k_nm == ComplexF64(-want_cancellation, 0.0)
    end


    @testset "scale-safe huge finite data and angular-label overflow" begin
        r = 1.0e155
        B = (8.0e307, 1.0e155, 0.5, 0.0, 0.0)
        Q = (7.0e307, 2.0e155, 0.75, 0.0, 0.0)
        mode, partner = paired_static_jets(Val(2); r, values=B, partner_values=Q)
        metric = static_hertz_metric(PhysicalHertzPair(mode, partner;
            pair_fingerprint=GHZ_METRIC_PFP))
        @test same_hertz_metric(metric, hertz_metric(PhysicalHertzPair(mode, partner;
            pair_fingerprint=GHZ_METRIC_PFP)))
        R = BigFloat(r)
        PB, BB = BigFloat.(Q), BigFloat.(B)
        lambda2, c = sqrt(BigFloat(24)), inv(sqrt(BigFloat(2)))
        want = (
            -(lambda2 / 4) * (PB[1] / R^2 + BB[1] / R^2),
            c * (PB[2] / R - 2PB[1] / R^2),
            -c * (BB[2] / R - 2BB[1] / R^2),
            -PB[3] / 2 + PB[2] / R,
            -BB[3] / 2 + BB[2] / R,
        )
        got = (metric.k_nn, metric.k_nbar, metric.k_nm,
            metric.k_barbar, metric.k_mm)
        @test all(isfinite, got)
        @test all(isapprox.(got, want; rtol=8e-15, atol=0))

        huge_l_mode, huge_l_partner = paired_static_jets(Val(2); l=typemax(Int))
        huge_l_pair = PhysicalHertzPair(huge_l_mode, huge_l_partner;
            pair_fingerprint=GHZ_METRIC_PFP)
        @test_throws ArgumentError static_hertz_metric(huge_l_pair)
    end

    @testset "small-radius near-maximum cancellation" begin
        F = floatmax(Float64)
        M, r = 0.25, 0.75
        B = (0.0, 0.0, 0.0, 0.0, 0.0)
        Q = (0.0, 0.9F, F, 0.0, 0.0)
        mode, partner = paired_static_jets(Val(2); r, values=B, partner_values=Q, M)
        metric = static_hertz_metric(PhysicalHertzPair(mode, partner;
            pair_fingerprint=GHZ_METRIC_PFP))
        @test same_hertz_metric(metric, hertz_metric(PhysicalHertzPair(mode, partner;
            pair_fingerprint=GHZ_METRIC_PFP)))

        R = BigFloat(r)
        PB = BigFloat.(Q)
        d = sqrt(BigFloat(2))
        want = (
            zero(BigFloat),
            d * (PB[2] / (2R)),
            zero(BigFloat),
            PB[2] / R - PB[3] / 2,
            zero(BigFloat),
        )
        got = (metric.k_nn, metric.k_nbar, metric.k_nm,
            metric.k_barbar, metric.k_mm)
        @test all(isfinite, got)
        @test all(isapprox.(got, want; rtol=5e-15, atol=0))
        @test real(metric.k_nbar) > 0.84F
        @test real(metric.k_barbar) > 0.69F
    end

    @testset "dispatch has no internal ambiguities" begin
        for callable in (ConjugateHertzRadialJet, PhysicalHertzPair,
                hertz_metric, static_hertz_metric)
            ms = collect(methods(callable))
            @test all(!Base.isambiguous(ms[i], ms[j])
                for i in eachindex(ms) for j in (i + 1):length(ms))
        end
    end
end

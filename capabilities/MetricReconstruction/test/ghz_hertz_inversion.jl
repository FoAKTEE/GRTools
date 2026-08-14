using Test
using MetricReconstruction.GHZ

const GHZ_HERTZ_FP = UInt128(0xd2736ba0d54a47b29f396b65362baef1)
const GHZ_HERTZ_BFP = UInt128(0xb6ee7131655e4da687558a68bcd55b7c)

function hertz_mode(::Val{N}, psi0; T=Float64, M=one(T), l=2, m=1,
        omega=0.17 - 0.03im, fp=GHZ_HERTZ_FP, content=LinearizedPsi0(),
        phase=ExpMinusIomegaUPlusImPhi()) where {N}
    ReducedWeylMode{N}(psi0; M=T(M), l, m, omega=Complex{T}(omega),
        source_fingerprint=fp, content, phase)
end

function hertz_boundary(src; radial_derivative=FixedOutgoingURadialDerivative(),
        policy=FiniteHertzRadialIVP(), kwargs...)
    HertzRadialBoundary(src; radial_derivative, policy, kwargs...)
end

function raw_hertz_boundary(bc::HertzRadialBoundary{N,T};
        source_fingerprint=bc.source_fingerprint,
        boundary_fingerprint=bc.boundary_fingerprint, rb=bc.rb, b0=bc.b0,
        b1=bc.b1, b2=bc.b2, b3=bc.b3) where {N,T}
    HertzRadialBoundary{N,T}(bc.M, bc.l, bc.m, bc.omega, source_fingerprint,
        bc.tetrad, bc.slicing, bc.phase, bc.radial_derivative, bc.policy, rb,
        b0, b1, b2, b3, boundary_fingerprint, bc.normalization, bc.gauge)
end

@testset "GHZ — outgoing-Schwarzschild IRG Hertz inversion" begin
    @testset "analytic polynomial and the positive quarter-thorn4 factor" begin
        rb = 3.0
        src = hertz_mode(Val(2), r -> r)
        bc = hertz_boundary(src; rb, boundary_fingerprint=GHZ_HERTZ_BFP)
        mode = invert_hertz(src, bc)
        @test mode isa ConjugateHertzMode
        @test mode.normalization isa QuarterThorn4
        @test mode.gauge isa HertzIRG
        @test mode.source.content isa LinearizedPsi0
        @test mode.source.phase isa ExpMinusIomegaUPlusImPhi
        @test mode.boundary.phase isa ExpMinusIomegaUPlusImPhi
        @test mode.phase isa ExpMinusIomegaUPlusImPhi
        @test mode.boundary.radial_derivative isa FixedOutgoingURadialDerivative
        @test mode.radial_derivative isa FixedOutgoingURadialDerivative
        @test mode.boundary.policy isa FiniteHertzRadialIVP
        @test mode.policy isa FiniteHertzRadialIVP

        expected(delta) = (
            delta^4 / 2 + delta^5 / 30,
            2delta^3 + delta^4 / 6,
            6delta^2 + 2delta^3 / 3,
            12delta + 2delta^2,
            12 + 4delta,
        )
        for r in (rb, 3.5, 5.0, 8.25)
            want = expected(r - rb)
            for n in 0:4
                @test bar_hertz_derivative(mode, r, n) ≈ want[n + 1] rtol=3e-13 atol=3e-14
            end
            jet = bar_hertz_jet(mode, r)
            @test jet isa NTuple{5,ComplexF64}
            @test all(isapprox.(jet, want; rtol=3e-13, atol=3e-14))

            # Canonical mostly-plus IRG normalization: +1/4 d_r^4 bar(Phi) = psi0L.
            d4 = bar_hertz_derivative(mode, r, 4)
            @test d4 / 4 == ComplexF64(r)
            @test d4 / 2 != ComplexF64(r)
            @test -d4 / 4 != ComplexF64(r)
        end
    end

    @testset "homogeneous constants, complex data, and linearity" begin
        rb = 3.0
        b = (1.2 - 0.4im, -0.3 + 0.8im, 0.5 - 0.2im, -0.7 + 0.1im)
        zero_src = hertz_mode(Val(2), r -> zero(ComplexF64); fp=GHZ_HERTZ_FP + 1)
        bc = hertz_boundary(zero_src; rb, b0=b[1], b1=b[2], b2=b[3], b3=b[4],
            boundary_fingerprint=GHZ_HERTZ_BFP + 1)
        homogeneous = invert_hertz(zero_src, bc)
        for r in (rb, 4.25, 7.0), n in 0:3
            delta = r - rb
            want = sum(b[j + 1] * delta^(j - n) / factorial(j - n) for j in n:3)
            @test bar_hertz_derivative(homogeneous, r, n) ≈ want atol=2e-14
        end
        @test bar_hertz_derivative(homogeneous, 5.0, 4) == 0

        A(r) = (0.2 + 0.3im) * r
        B(r) = (-0.4 + 0.1im) / r
        a, c = 1.3 - 0.2im, -0.6 + 0.7im
        srcA = hertz_mode(Val(2), A; fp=GHZ_HERTZ_FP + 2)
        srcB = hertz_mode(Val(2), B; fp=GHZ_HERTZ_FP + 3)
        srcL = hertz_mode(Val(2), r -> a * A(r) + c * B(r); fp=GHZ_HERTZ_FP + 4)
        solve(src, off) = invert_hertz(src, hertz_boundary(src; rb,
            boundary_fingerprint=GHZ_HERTZ_BFP + off))
        solA, solB, solL = solve(srcA, 2), solve(srcB, 3), solve(srcL, 4)
        for r in (rb, 4.0, 6.5), n in 0:4
            @test bar_hertz_derivative(solL, r, n) ≈
                a * bar_hertz_derivative(solA, r, n) +
                c * bar_hertz_derivative(solB, r, n) rtol=5e-12 atol=4e-14
        end
    end

    @testset "perturbative order one and two" begin
        for N in (1, 2)
            src = hertz_mode(Val(N), r -> 0.25r; fp=GHZ_HERTZ_FP + 10 + N)
            bc = hertz_boundary(src; rb=3.0,
                boundary_fingerprint=GHZ_HERTZ_BFP + 10 + N)
            mode = invert_hertz(src, bc)
            @test perturbative_order(src) == N
            @test perturbative_order(bc) == N
            @test perturbative_order(mode) == N
            @test bar_hertz_derivative(mode, 4.0, 4) == 4.0 + 0im
        end
    end

    @testset "BigFloat precision" begin
        setprecision(BigFloat, 256) do
            C = Complex{BigFloat}(big"0.125", big"-0.375")
            rb, r = big"3.0", big"5.25"
            src = hertz_mode(Val(2), q -> C * q; T=BigFloat,
                omega=Complex{BigFloat}(big"0.02", big"-0.01"), fp=GHZ_HERTZ_FP + 20)
            bc = hertz_boundary(src; rb, boundary_fingerprint=GHZ_HERTZ_BFP + 20)
            mode = invert_hertz(src, bc)
            delta = r - rb
            want = C * (delta^4 / 2 + delta^5 / 30)
            got = bar_hertz_derivative(mode, r, 0; rtol=big"1e-60")
            @test got isa Complex{BigFloat}
            @test all(x -> x isa Complex{BigFloat},
                bar_hertz_jet(mode, r; rtol=big"1e-60"))
            @test abs(got - want) < big"1e-55"
            @test bar_hertz_derivative(mode, r, 4; rtol=big"1e-60") == 4C * r
        end
    end

    @testset "metadata, boundary, and full-Weyl rejection" begin
        @test_throws ArgumentError hertz_mode(Val(0), r -> r)
        @test_throws ArgumentError hertz_mode(Val(2), r -> r; l=1)
        @test_throws ArgumentError hertz_mode(Val(2), r -> r; l=2, m=3)
        @test_throws ArgumentError hertz_mode(Val(2), r -> r; M=0.0)
        @test_throws ArgumentError hertz_mode(Val(2), r -> r; omega=Inf + 0im)
        @test_throws ArgumentError hertz_mode(Val(2), r -> r; fp=0)
        @test_throws ArgumentError ReducedWeylMode{2}(r -> r; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_HERTZ_FP, tetrad=:wrong,
            phase=ExpMinusIomegaUPlusImPhi())
        @test_throws ArgumentError ReducedWeylMode{2}(r -> r; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_HERTZ_FP, slicing=:wrong,
            phase=ExpMinusIomegaUPlusImPhi())
        @test_throws ArgumentError ReducedWeylMode{2}(r -> r; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_HERTZ_FP, phase=:wrong)
        @test_throws UndefKeywordError ReducedWeylMode{2}(r -> r; M=1.0, l=2, m=0,
            source_fingerprint=GHZ_HERTZ_FP)

        src = hertz_mode(Val(2), r -> r)
        @test_throws ArgumentError hertz_boundary(src; rb=1.99,
            boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0, b2=Inf,
            boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0,
            boundary_fingerprint=0)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP, normalization=:wrong)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP, gauge=GHZIRG())
        @test_throws ArgumentError HertzRadialBoundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP,
            radial_derivative=:fixed_t, policy=FiniteHertzRadialIVP())
        @test_throws ArgumentError HertzRadialBoundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP,
            radial_derivative=FixedOutgoingURadialDerivative(), policy=:retarded)
        @test_throws UndefKeywordError HertzRadialBoundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP)

        bc = hertz_boundary(src; rb=3.0, boundary_fingerprint=GHZ_HERTZ_BFP)
        mode = invert_hertz(src, bc)
        @test bc.boundary_fingerprint == GHZ_HERTZ_BFP
        @test (bc.M, bc.l, bc.m, bc.omega, bc.source_fingerprint) ==
            (src.M, src.l, src.m, src.omega, src.source_fingerprint)

        raw_zero_boundary_fp = raw_hertz_boundary(bc; boundary_fingerprint=UInt128(0))
        @test_throws ArgumentError invert_hertz(src, raw_zero_boundary_fp)
        @test_throws ArgumentError bar_hertz_jet(
            ConjugateHertzMode{2,Float64,typeof(src.psi0)}(src, raw_zero_boundary_fp,
                mode.normalization, mode.gauge, mode.phase, mode.radial_derivative,
                mode.policy), 3.0)
        raw_zero_source = ReducedWeylMode{2,Float64,LinearizedPsi0,typeof(src.psi0)}(
            src.psi0, src.M, src.l, src.m, src.omega, UInt128(0), src.content,
            src.tetrad, src.slicing, src.phase)
        raw_zero_source_bc = raw_hertz_boundary(bc; source_fingerprint=UInt128(0))
        @test_throws ArgumentError invert_hertz(raw_zero_source, raw_zero_source_bc)

        other = hertz_mode(Val(2), r -> r; l=3, fp=GHZ_HERTZ_FP + 31)
        @test_throws ArgumentError invert_hertz(src, hertz_boundary(other; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 31))
        other_N = hertz_mode(Val(1), r -> r; fp=GHZ_HERTZ_FP)
        @test_throws ArgumentError invert_hertz(src, hertz_boundary(other_N; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 32))

        full = hertz_mode(Val(2), r -> r; content=FullSecondOrderWeyl(),
            fp=GHZ_HERTZ_FP + 33)
        full_bc = hertz_boundary(full; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 33)
        @test full.content isa FullSecondOrderWeyl
        @test_throws ArgumentError invert_hertz(full, full_bc)

        @test_throws ArgumentError bar_hertz_derivative(mode, 2.9, 0)
        @test_throws ArgumentError bar_hertz_derivative(mode, Inf, 0)
        @test_throws ArgumentError bar_hertz_derivative(mode, 3.0, -1)
        @test_throws ArgumentError bar_hertz_derivative(mode, 3.0, 5)
        @test_throws ArgumentError bar_hertz_derivative(mode, 3.0, true)
        @test_throws ArgumentError bar_hertz_derivative(mode, 3.0, 1.0)
    end

    @testset "checked harmonic mode indices" begin
        indexed_weyl(l, m) = ReducedWeylMode{2}(r -> r; M=1.0, l, m, omega=0.0,
            source_fingerprint=GHZ_HERTZ_FP + 120,
            phase=ExpMinusIomegaUPlusImPhi())

        @test_throws ArgumentError indexed_weyl(true, 0)
        @test_throws ArgumentError indexed_weyl(2, true)
        @test_throws ArgumentError indexed_weyl(2, typemin(Int))
        @test_throws ArgumentError indexed_weyl(typemax(Int), typemin(Int))

        too_big = BigInt(typemax(Int)) + 1
        too_uint = UInt(typemax(Int)) + one(UInt)
        for invalid in (too_big, too_uint)
            @test_throws ArgumentError indexed_weyl(invalid, 0)
            @test_throws ArgumentError indexed_weyl(2, invalid)
        end

        for (l, m, expected) in ((BigInt(2), BigInt(-2), (2, -2)),
                (UInt(2), UInt(2), (2, 2)), (Int128(2), Int128(-2), (2, -2)))
            src = indexed_weyl(l, m)
            @test (src.l, src.m) == expected
        end

        for m in (-typemax(Int), typemax(Int))
            extreme = indexed_weyl(typemax(Int), m)
            @test (extreme.l, extreme.m) == (typemax(Int), m)
        end
    end

    @testset "Boolean scalar rejection" begin
        radial = r -> r
        @test_throws ArgumentError ReducedWeylMode{2}(radial; M=true, l=2, m=0,
            omega=0.0, source_fingerprint=GHZ_HERTZ_FP,
            phase=ExpMinusIomegaUPlusImPhi())
        @test_throws ArgumentError ReducedWeylMode{2}(radial; M=1.0, l=2, m=0,
            omega=true, source_fingerprint=GHZ_HERTZ_FP,
            phase=ExpMinusIomegaUPlusImPhi())
        @test_throws ArgumentError ReducedWeylMode{2}(radial; M=1.0, l=2, m=0,
            omega=Complex{Bool}(true, false), source_fingerprint=GHZ_HERTZ_FP,
            phase=ExpMinusIomegaUPlusImPhi())
        @test_throws ArgumentError ReducedWeylMode{2}(radial; M=1.0, l=2, m=0,
            omega=0.0, source_fingerprint=true,
            phase=ExpMinusIomegaUPlusImPhi())

        low_mass_src = ReducedWeylMode{2}(radial; M=0.25, l=2, m=0, omega=0.0,
            source_fingerprint=GHZ_HERTZ_FP + 113,
            phase=ExpMinusIomegaUPlusImPhi())
        for bad_radius in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError hertz_boundary(low_mass_src; rb=bad_radius,
                boundary_fingerprint=GHZ_HERTZ_BFP + 113)
        end
        low_mass_mode = invert_hertz(low_mass_src,
            hertz_boundary(low_mass_src; rb=0.5,
                boundary_fingerprint=GHZ_HERTZ_BFP + 113))
        for bad_radius in (true, Complex{Bool}(true, false))
            @test_throws ArgumentError bar_hertz_derivative(low_mass_mode, bad_radius, 0)
            @test_throws ArgumentError bar_hertz_jet(low_mass_mode, bad_radius)
        end

        src = hertz_mode(Val(2), radial; fp=GHZ_HERTZ_FP + 110)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0, b0=true,
            boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0, b1=true,
            boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0, b2=true,
            boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0, b3=true,
            boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0,
            b0=Complex{Bool}(true, false), boundary_fingerprint=GHZ_HERTZ_BFP)
        @test_throws ArgumentError hertz_boundary(src; rb=3.0,
            boundary_fingerprint=true)

        boundary_bool = hertz_mode(Val(2), r -> true; fp=GHZ_HERTZ_FP + 111)
        @test_throws ArgumentError invert_hertz(boundary_bool,
            hertz_boundary(boundary_bool; rb=3.0,
                boundary_fingerprint=GHZ_HERTZ_BFP + 111))
        late_bool = hertz_mode(Val(2), r -> r > 3.2 ? true : r;
            fp=GHZ_HERTZ_FP + 112)
        late_mode = invert_hertz(late_bool, hertz_boundary(late_bool; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 112))
        @test_throws ArgumentError bar_hertz_derivative(late_mode, 4.0, 4)

        mode = invert_hertz(src, hertz_boundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 110))
        @test_throws ArgumentError bar_hertz_derivative(mode, 3.0, 4; atol=true)
        @test_throws ArgumentError bar_hertz_jet(mode, 3.0; rtol=true)
    end

    @testset "nonfinite source and tolerance rejection" begin
        bad = hertz_mode(Val(2), r -> NaN; fp=GHZ_HERTZ_FP + 40)
        bad_bc = hertz_boundary(bad; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 40)
        @test_throws ArgumentError invert_hertz(bad, bad_bc)

        late = hertz_mode(Val(2), r -> r > 3.2 ? Inf : r; fp=GHZ_HERTZ_FP + 41)
        late_mode = invert_hertz(late, hertz_boundary(late; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 41))
        @test_throws ArgumentError bar_hertz_derivative(late_mode, 4.0, 0)
        @test_throws ArgumentError bar_hertz_jet(late_mode, 4.0)

        src = hertz_mode(Val(2), r -> r; fp=GHZ_HERTZ_FP + 42)
        mode = invert_hertz(src, hertz_boundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 42))
        for kwargs in ((; atol=-1.0), (; rtol=-1.0), (; atol=NaN), (; rtol=Inf),
                       (; atol=0.0, rtol=0.0))
            @test_throws ArgumentError bar_hertz_derivative(mode, 3.0, 4; kwargs...)
            @test_throws ArgumentError bar_hertz_jet(mode, 3.0; kwargs...)
        end
    end

    @testset "finite large-ray scaling" begin
        # Forming delta^2 first overflows Float16, despite the weighted term being O(1).
        T16 = Float16
        R16 = T16(400)
        src16 = hertz_mode(Val(2), r -> zero(T16); T=T16, omega=zero(T16),
            fp=GHZ_HERTZ_FP + 50)
        b2 = inv(R16) * inv(R16)
        mode16 = invert_hertz(src16, hertz_boundary(src16; rb=T16(2), b2,
            boundary_fingerprint=GHZ_HERTZ_BFP + 50))
        want16 = b2 * (R16 - T16(2)) * (R16 - T16(2)) / T16(2)
        @test isfinite(bar_hertz_derivative(mode16, R16, 0))
        @test bar_hertz_derivative(mode16, R16, 0) ≈ want16 rtol=T16(0.01)

        # Float64 discriminators for quadratic and cubic homogeneous weights.
        R2 = 1e155
        src2 = hertz_mode(Val(2), r -> 0.0; fp=GHZ_HERTZ_FP + 51)
        b2large = inv(R2) * inv(R2)
        mode2 = invert_hertz(src2, hertz_boundary(src2; rb=3.0, b2=b2large,
            boundary_fingerprint=GHZ_HERTZ_BFP + 51))
        want2 = (b2large * (R2 - 3.0)) * ((R2 - 3.0) / 2)
        @test bar_hertz_derivative(mode2, R2, 0) ≈ want2 rtol=3e-15

        R3 = 1e103
        b3large = inv(R3)^3
        mode3 = invert_hertz(src2, hertz_boundary(src2; rb=3.0, b3=b3large,
            boundary_fingerprint=GHZ_HERTZ_BFP + 52))
        want3 = ((b3large * (R3 - 3.0)) * ((R3 - 3.0) / 2)) *
            ((R3 - 3.0) / 3)
        @test bar_hertz_derivative(mode3, R3, 0) ≈ want3 rtol=4e-15

        # Source-first recurrence prevents an overflowing cubic kernel and preserves an exact
        # zero source: neither path may create Inf*0 inside QuadGK.
        C = inv(R3)^3
        sourced = hertz_mode(Val(2), r -> C; fp=GHZ_HERTZ_FP + 53)
        sourced_mode = invert_hertz(sourced, hertz_boundary(sourced; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 53))
        delta = R3 - 3.0
        want_source = (((C * delta) * (delta / 2)) * (delta / 3)) * delta
        got_source = bar_hertz_derivative(sourced_mode, R3, 0; rtol=1e-12)
        @test isfinite(got_source)
        @test got_source ≈ want_source rtol=3e-12
        @test isfinite(bar_hertz_derivative(mode3, R3, 0))
    end

    @testset "large homogeneous cancellation" begin
        src = hertz_mode(Val(2), r -> 0.0; fp=GHZ_HERTZ_FP + 60)

        # Direct refuter discriminator: each power-expanded term overflows, but their Taylor
        # polynomial evaluated from the actual Float64 inputs is finite.
        R0 = 1e155
        b3_0 = -3 / R0
        mode0 = invert_hertz(src, hertz_boundary(src; rb=3.0, b2=1.0, b3=b3_0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 60))
        want0 = setprecision(BigFloat, 256) do
            delta = BigFloat(R0) - BigFloat(3.0)
            Float64(delta^2 / 2 + BigFloat(b3_0) * delta^3 / 6)
        end
        got0 = bar_hertz_derivative(mode0, R0, 0)
        @test isfinite(got0)
        @test got0 ≈ want0 rtol=3e-15

        # Scaling the Horner offset by its small integer denominator may itself overflow; the
        # alternate fused form must preserve this finite large-coefficient polynomial.
        large_offset = invert_hertz(src, hertz_boundary(src; rb=3.0, b2=1e308,
            boundary_fingerprint=GHZ_HERTZ_BFP + 64))
        @test bar_hertz_derivative(large_offset, 4.0, 0) ≈ 5e307 rtol=2e-15

        # A denominator-scaled numerator would overflow to 1.8e308 before division by two.
        finite_n1 = invert_hertz(src, hertz_boundary(src; rb=3.0, b2=7e307,
            b3=4e307, boundary_fingerprint=GHZ_HERTZ_BFP + 65))
        @test bar_hertz_derivative(finite_n1, 4.0, 1) ≈ 9e307 rtol=2e-15

        # Exact binary cancellations exercise the analogous n=1 and n=2 truncations without
        # relying on decimal-rounding accidents.  The separate n=1 terms are both infinite.
        R = ldexp(1.0, 800)
        b2_1 = ldexp(1.0, 400)
        b3_1 = -ldexp(1.0, -399)
        mode1 = invert_hertz(src, hertz_boundary(src; rb=3.0, b1=7.0,
            b2=b2_1, b3=b3_1, boundary_fingerprint=GHZ_HERTZ_BFP + 61))
        @test !isfinite(b2_1 * R)
        @test bar_hertz_derivative(mode1, R, 1) == 7.0 + 0im

        b2_2 = ldexp(1.0, 900)
        b3_2 = -ldexp(1.0, 100)
        mode2 = invert_hertz(src, hertz_boundary(src; rb=3.0, b2=b2_2, b3=b3_2,
            boundary_fingerprint=GHZ_HERTZ_BFP + 62))
        @test bar_hertz_derivative(mode2, R, 2) == 0.0 + 0im
        @test bar_hertz_derivative(mode2, R, 3) == ComplexF64(b3_2)

        # Fuse the particular factor of four with the homogeneous contribution: `4I` alone
        # overflows here, while the combined derivative is finite.
        C = 0.375 * floatmax(Float64)
        cancel_src = hertz_mode(Val(2), r -> C; fp=GHZ_HERTZ_FP + 63)
        cancel_mode = invert_hertz(cancel_src, hertz_boundary(cancel_src; rb=3.0,
            b3=-floatmax(Float64), boundary_fingerprint=GHZ_HERTZ_BFP + 63))
        @test !isfinite(4C)
        got_cancel = bar_hertz_derivative(cancel_mode, 4.0, 3)
        @test isfinite(got_cancel)
        @test got_cancel ≈ 0.5floatmax(Float64) rtol=2e-14
    end

    @testset "widened finite-result fallback" begin
        function oracle(T, data, delta, n, integral=zero(Complex{T}))
            setprecision(BigFloat, 256) do
                db = BigFloat(delta)
                total = Complex{BigFloat}(BigFloat(real(integral)), BigFloat(imag(integral))) * 4
                for j in n:3
                    k = j - n
                    bj = Complex{BigFloat}(BigFloat(real(data[j + 1])),
                                           BigFloat(imag(data[j + 1])))
                    total += bj * db^k / BigFloat(factorial(k))
                end
                Complex{T}(T(real(total)), T(imag(total)))
            end
        end

        for (offset, T) in enumerate((Float16, Float32, Float64))
            F = floatmax(T)
            data = Complex{T}.((T(0.1) * F, -T(0.5) * F, T(0.9) * F, T(0.8) * F))
            src = hertz_mode(Val(2), r -> zero(T); T, omega=zero(T),
                fp=GHZ_HERTZ_FP + 70 + offset)
            mode = invert_hertz(src, hertz_boundary(src; rb=T(2),
                b0=data[1], b1=data[2], b2=data[3], b3=data[4],
                boundary_fingerprint=GHZ_HERTZ_BFP + 70 + offset))
            delta = T(0.5)
            got = ntuple(i -> bar_hertz_derivative(mode, T(2.5), i - 1), 4)
            want = ntuple(i -> oracle(T, data, delta, i - 1), 4)
            @test got[1] == want[1]
            @test got[2] == want[2]
            @test !isfinite(got[3]) && !isfinite(want[3])
            @test got[4] == want[4]
        end

        # The homogeneous n=2 polynomial is 1.6F and the particular piece is -F; only their
        # complete widened composition recovers the finite 0.6F result.
        F = floatmax(Float64)
        C = -0.5F
        src = hertz_mode(Val(2), r -> C; fp=GHZ_HERTZ_FP + 80)
        data = ComplexF64.((0.0, 0.0, 0.8F, 0.8F))
        mode = invert_hertz(src, hertz_boundary(src; rb=3.0, b2=data[3], b3=data[4],
            boundary_fingerprint=GHZ_HERTZ_BFP + 80))
        integral = ComplexF64(C / 2)
        want = oracle(Float64, data, 1.0, 2, integral)
        got = bar_hertz_derivative(mode, 4.0, 2)
        @test isfinite(got)
        @test got ≈ want rtol=3e-15

        # The ordinary T-valued integral is 1.2F and overflows before composition.  Risky
        # boundary data must route directly to BigFloat accumulation so -5F + 4.8F stays finite.
        for (offset, T) in enumerate((Float16, Float32, Float64))
            F = floatmax(T)
            C = T(0.15) * F
            delta = T(4)
            data = Complex{T}.((zero(T), zero(T), -F, -F))
            src = hertz_mode(Val(2), r -> C; T, omega=zero(T),
                fp=GHZ_HERTZ_FP + 90 + offset)
            mode = invert_hertz(src, hertz_boundary(src; rb=T(2), b2=-F, b3=-F,
                boundary_fingerprint=GHZ_HERTZ_BFP + 90 + offset))
            want = setprecision(BigFloat, 256) do
                db, cb, fb = BigFloat(delta), BigFloat(C), BigFloat(F)
                Complex{T}(T(-fb - fb * db + 4cb * db^2 / 2), zero(T))
            end
            got = bar_hertz_derivative(mode, T(6), 2)
            @test isfinite(got)
            @test got ≈ want rtol=8eps(T)

            # With no cancelling boundary the widened retry should return a true overflow,
            # rather than leak QuadGK's nonfinite-integrand DomainError.
            overflow_mode = invert_hertz(src, hertz_boundary(src; rb=T(2),
                boundary_fingerprint=GHZ_HERTZ_BFP + 93 + offset))
            @test !isfinite(bar_hertz_derivative(overflow_mode, T(6), 2))
        end

        # A nonconstant source sees only Float64 ray coordinates on both paths, and forced-wide
        # accumulation agrees with the ordinary path after adding the known homogeneous seven.
        typed_ray = Ref(true)
        radial(q) = (typed_ray[] &= q isa Float64; (0.01 - 0.02im) * q)
        src = hertz_mode(Val(2), radial; fp=GHZ_HERTZ_FP + 100)
        normal = invert_hertz(src, hertz_boundary(src; rb=3.0,
            boundary_fingerprint=GHZ_HERTZ_BFP + 100))
        forced_wide = invert_hertz(src, hertz_boundary(src; rb=3.0,
            b0=7.0, b2=1.0, b3=-3.0, boundary_fingerprint=GHZ_HERTZ_BFP + 101))
        ordinary_value = bar_hertz_derivative(normal, 4.0, 0)
        typed_ray[] = true
        wide_value = bar_hertz_derivative(forced_wide, 4.0, 0)
        @test typed_ray[]
        @test wide_value ≈ ordinary_value + 7 rtol=3e-13 atol=3e-14
    end
end

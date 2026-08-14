# Grounding for Radiation WeylAmplitudes — Teukolsky–Starobinsky constants + CCK reconstruction presets.
#
# Reference values are the UPSTREAM Mathematica (wolframscript) oracle from ExampleUsage.nb cell-17/19 at
# the cell-21 config (a=1/10, M=1, l=3, m=0, ω=0.9, αin=1, αout=0), baked verbatim from the legacy
# MetricReconstruction verify_weyl_amplitudes.jl reference block (20-digit). teuk_lambda comes from the
# already-grounded GrTool Harmonics surface; _horizons/_Omega_plus are delegated to Spacetime.
#
# [HOLE] crosswalk H10: the :psi0_IRG / :psi4_ORG presets are implemented and their numeric B-values match
# the .nb oracle below, BUT the preset CHOICE is not exercised by any upstream vacuum-Einstein gate — the
# physics-gate validation of those two presets is the residual open hole (the calD/calC/Dhat core and the
# self-contained :psi0_ORG/:psi4_IRG presets are fully grounded).

using Test
using MetricReconstruction.Radiation
using Spacetime: horizons, omega_plus

@testset "Radiation.WeylAmplitudes (Teukolsky–Starobinsky + CCK presets)" begin
    a = 1/10; M = 1; l = 3; m = 0; ω = 0.9

    @testset "TS constants vs cell-21 wolframscript oracle" begin
        @test teuk_lambda(2, l, m, a, ω) ≈ 6.005400441683527 rtol=1e-12
        @test calC(ω, l, m, a)      ≈ 14516.639997533973 rtol=1e-12
        @test calD(ω, l, m, a)      ≈ 14399.999997533976 rtol=1e-12
        @test Dhat(ω, l, m, a)      ≈ 113.01215056379583 rtol=1e-12
        @test DhatPrime(ω, l, m, a) ≈ 127.41992719982011 rtol=1e-12
    end

    @testset "|m|≥2 TS factors: continuation, identity, and numeric type" begin
        # Literal paper-oracle factors F(|m|), not values recomputed through the
        # implementation under test.  For m≤-2 this is the finite continuation
        # of (m+1)!/(m-3)!; evaluating those factorials separately is singular.
        factors = (2 => 24, 3 => 120, 4 => 360)

        for T in (Float32, Float64)
            ωT = T(0.9); aT = T(0.1)
            for (q, factor) in factors, mm in (-q, q)
                d = Dhat(ωT, q, mm, aT)
                dp = DhatPrime(ωT, q, mm, aT)
                @test typeof(d) === T
                @test typeof(dp) === T
                @test d * dp ≈ calD(ωT, q, mm, aT) rtol=50eps(T)
                mm < 0 && @test d == T(factor)
                mm > 0 && @test dp == T(factor)

                p0o = weyl_amplitudes(:psi0_ORG, ωT, q, mm, aT)
                p4i = weyl_amplitudes(:psi4_IRG, ωT, q, mm, aT)
                @test typeof(p0o.B) === T
                @test typeof(p4i.B) === T
                @test all(isfinite, (p0o.A, p0o.Abar, p0o.B, p0o.Bbar))
                @test all(isfinite, (p4i.A, p4i.Abar, p4i.B, p4i.Bbar))
            end

            ωc = complex(ωT, T(0.03))
            for mm in (-2, 2)
                d = Dhat(ωc, 2, mm, aT)
                dp = DhatPrime(ωc, 2, mm, aT)
                @test typeof(d) === Complex{T}
                @test typeof(dp) === Complex{T}
                @test d * dp ≈ calD(ωc, 2, mm, aT) rtol=50eps(T)
                @test typeof(weyl_amplitudes(:psi0_ORG, ωc, 2, mm, aT).B) === Complex{T}
                @test typeof(weyl_amplitudes(:psi4_IRG, ωc, 2, mm, aT).B) === Complex{T}
            end
        end

        # Float inputs must be completely independent of ambient MPFR precision.
        function float_snapshot(bits, T)
            return setprecision(BigFloat, bits) do
                ωT = T(0.9); aT = T(0.1); ll = 2; mm = -2
                p0o = weyl_amplitudes(:psi0_ORG, ωT, ll, mm, aT)
                p4i = weyl_amplitudes(:psi4_IRG, ωT, ll, mm, aT)
                (Dhat(ωT, ll, mm, aT), DhatPrime(ωT, ll, mm, aT), p0o.B, p4i.B)
            end
        end
        for T in (Float32, Float64)
            lo = float_snapshot(32, T)
            hi = float_snapshot(512, T)
            @test all(x -> typeof(x) === T, lo)
            @test all(x -> typeof(x) === T, hi)
            @test isequal(lo, hi)
        end

        setprecision(BigFloat, 256) do
            ωb = big"0.9"; ab = big"0.1"
            for mm in (-2, 2)
                d = Dhat(ωb, 2, mm, ab)
                dp = DhatPrime(ωb, 2, mm, ab)
                @test d isa BigFloat
                @test dp isa BigFloat
                @test precision(d) == 256
                @test precision(dp) == 256
                @test d * dp ≈ calD(ωb, 2, mm, ab) rtol=eps(BigFloat) * 50
            end
        end

        # Constant-factor branches still preserve the harmonic-domain contract
        # formerly enforced as a side effect of an unnecessary eigensolve.
        @test_throws ArgumentError Dhat(0.9, 1, -2, 0.1)
        @test_throws ArgumentError DhatPrime(0.9, 1, 2, 0.1)
        @test_throws ArgumentError Dhat(0.9, 2, typemin(Int), 0.1)
        @test_throws MethodError Dhat(0.9, 2, -2, 0.1 + 0im)
        @test_throws MethodError DhatPrime(0.9, 2, 2, 0.1 + 0im)
    end

    @testset "radial constants Γ/Γ̃ + in/out (cell 17)" begin
        Γ_ref  = 2864.804014718387 + 1587.5588718000337im
        Γt_ref = 2864.804014718387 - 1587.5588718000337im
        @test radial_in(ω, l, m, a)       ≈ Γ_ref  rtol=1e-12
        @test radial_outPrime(ω, l, m, a) ≈ Γt_ref rtol=1e-12
        @test radial_out(ω, l, m, a)      ≈ 3.876722529441985 + 2.148323310600947im rtol=1e-12
        @test radial_inPrime(ω, l, m, a)  ≈ 3.876722529441985 - 2.148323310600947im rtol=1e-12
        # _horizons/_Omega_plus delegate to Spacetime (already grounded)
        @test Radiation._horizons(a) == horizons(a)
        @test Radiation._Omega_plus(a) == omega_plus(a)
    end

    @testset "self-contained presets ψ0-ORG / ψ4-IRG (cell 19)" begin
        A_org = 192*im*ω*M / calC(ω, l, m, a)
        p0o = weyl_amplitudes(:psi0_ORG, ω, l, m, a; M=M)
        @test p0o.A ≈ A_org rtol=1e-12
        @test p0o.B ≈ 16*conj(Dhat(ω, l, m, a))/conj(calC(ω, l, m, a)) rtol=1e-12
        @test p0o.Abar == conj(p0o.A) && p0o.Bbar == conj(p0o.B)
        p4i = weyl_amplitudes(:psi4_IRG, ω, l, m, a; M=M)
        @test p4i.A ≈ -A_org rtol=1e-12
        @test p4i.B ≈ 16*conj(DhatPrime(ω, l, m, a))/conj(calC(ω, l, m, a)) rtol=1e-12
    end

    # [HOLE] H10: presets implemented + value-grounded vs the .nb oracle; physics-gate validation open.
    @testset "[HOLE H10] radial-asymptotic presets ψ0-IRG / ψ4-ORG vs cell-21 oracle" begin
        p0i = weyl_amplitudes(:psi0_IRG, ω, l, m, a; M=M, αin=1, αout=0)
        @test p0i.A == 0
        @test p0i.B ≈ 0.0010682148293546018 - 0.00059196158641831599im rtol=1e-10
        @test p0i.Bbar == conj(p0i.B)
        p4o = weyl_amplitudes(:psi4_ORG, ω, l, m, a; M=M, αin=1, αout=0)
        @test p4o.A == 0
        @test p4o.B ≈ 12.630158009919863 - 6.9991243023497299im rtol=1e-10
        @test p4o.Bbar == conj(p4o.B)
        # B identity at a non-trivial (αin,αout)
        αin, αout = 0.3, 0.7
        p0i2 = weyl_amplitudes(:psi0_IRG, ω, l, m, a; M=M, αin=αin, αout=αout)
        @test p0i2.B ≈ 4/(αin*radial_in(ω,l,m,a) + αout*radial_out(ω,l,m,a)) rtol=1e-12
        p4o2 = weyl_amplitudes(:psi4_ORG, ω, l, m, a; M=M, αin=αin, αout=αout)
        @test p4o2.B ≈ 64/conj(αin*radial_inPrime(ω,l,m,a) + αout*radial_outPrime(ω,l,m,a)) rtol=1e-12
    end

    @testset "BigFloat inheritance" begin
        setprecision(BigFloat, 256) do
            cC = calC(big"0.9", l, m, big"0.1")
            @test cC isa Union{BigFloat,Complex{BigFloat}}
            @test Float64(real(cC)) ≈ 14516.639997533973 rtol=1e-10
        end
    end
end

using Test
using AccretionDisc

@testset "Shared forcing vocabulary" begin
    @testset "Total frequency and Newtonian adapter" begin
        for (m, l, n) in ((2, 3, 2.0), (3, 3, 1), (1, 0, 1.0), (2, -1, 0.5),
                          (0, 2, 0.25), (0, 0, 1.0))
            mode = NewtonianMode(m, l, n)
            @test mode isa HarmonicMode{Float64}
            @test (mode.m, mode.k, mode.omega) == (m, l - m, l * n)
        end
        # Kepler limit: Omega_phi = Omega_r = n; l = m + k.
        for (m, k) in ((3, 0), (2, 1), (2, -3), (1, -1), (0, 2))
            gr = HarmonicMode(m, k, m * 0.25 + k * 0.25)
            newtonian = NewtonianMode(m, m + k, 0.25)
            @test gr == newtonian
        end
        @test HarmonicMode(Int32(2), Int32(-1), 1).omega === 1.0
        @test HarmonicMode{Rational{Int}}(2, -1, 1//3).omega == 1//3
    end

    @testset "Signed gravitational power into the disc" begin
        # omega = 6, pattern speed = 3, so a torque of 4 gives power 12.
        @test mode_power(NewtonianMode(2, 3, 2.0), 4) == 12.0
        @test mode_power(NewtonianMode(3, 3, 2.0), -4) == -8.0
        @test mode_power(NewtonianMode(2, -1, 0.5), 4) == -1.0
        @test mode_power(NewtonianMode(1, 0, 1.0), 4) == 0.0
        # GR apsidal mode: Omega_phi - Omega_r is nonzero.
        @test mode_power(HarmonicMode(1, -1, 0.125 - 0.0625), 4) == 0.25
        # Axisymmetric power needs its own integrand, even when torque is zero.
        @test_throws DomainError mode_power(HarmonicMode(0, 2, 0.5), 0.0)
        @test_throws DomainError mode_power(HarmonicMode(0, 0, 0.0), 0.0)
    end

    @testset "Constructor validation" begin
        @test_throws DomainError HarmonicMode(-1, 0, 1.0)
        @test_throws DomainError HarmonicMode{BigFloat}(-1, 0, 1)
        @test_throws DomainError NewtonianMode(-1, 2, 1)
        for bad in (NaN, Inf, -Inf)
            @test_throws DomainError HarmonicMode(1, 0, bad)
            @test_throws DomainError HarmonicMode{BigFloat}(1, 0, bad)
            @test_throws DomainError NewtonianMode(1, 2, bad)
        end
        for bad in (0, -1, -0.0)
            @test_throws DomainError NewtonianMode(1, 0, bad)
        end
        # Reject overflow after conversion, rather than storing a nonfinite frequency.
        @test_throws DomainError HarmonicMode{Float32}(1, 0, 1e100)
        @test_throws DomainError NewtonianMode(1, 2, floatmax(Float64))
        @test_throws OverflowError NewtonianMode(1, typemin(Int), 1.0)
        @test_throws InexactError HarmonicMode(big(typemax(Int)) + 1, 0, 1.0)
        @test_throws MethodError HarmonicMode(1.5, 0, 1.0)
        @test_throws MethodError NewtonianMode(1, 1.5, 1.0)
    end

    @testset "Softening policies are distinct, validated inputs" begin
        @test FixedSoftening(0).eps === 0.0
        @test SemimajorScaleHeight(0.3).b === 0.3
        @test InstantaneousScaleHeight(0.3).b === 0.3
        @test LocalScaleHeight(0.3).b === 0.3
        # FD25's constant length must never collapse into FR24's moving b*H(r_p).
        @test length(Set(typeof.((FixedSoftening(0.3), SemimajorScaleHeight(0.3),
                                 InstantaneousScaleHeight(0.3), LocalScaleHeight(0.3))))) == 4
        # FD25: pass the constant length once; no time/field radius is stored.
        @test FixedSoftening((3//10) * 0.1 * 2).eps ≈ 0.06
        for policy in (FixedSoftening, SemimajorScaleHeight,
                       InstantaneousScaleHeight, LocalScaleHeight)
            @test policy(1) isa SofteningPolicy{Float64}
            @test policy(0) isa SofteningPolicy{Float64}
            @test policy{BigFloat}(1) isa SofteningPolicy{BigFloat}
            @test policy{Rational{Int}}(1//3) isa SofteningPolicy{Rational{Int}}
            for bad in (-1, NaN, Inf, -Inf)
                @test_throws DomainError policy(bad)
                @test_throws DomainError policy{BigFloat}(bad)
            end
            @test_throws DomainError policy{Float32}(1e100)
        end
    end

    @testset "Closure and paper-mask tags" begin
        @test Adiabatic isa ClosureTag
        @test LocallyIsothermalOptionA isa ClosureTag
        @test Cooling isa ClosureTag
        @test length(Set((Adiabatic, LocallyIsothermalOptionA, Cooling))) == 3
        @test all(tag -> tag isa ModeMaskTag, (FR22Mask, FR24Mask, FD25Mask))
        @test length(Set((FR22Mask, FR24Mask, FD25Mask))) == 3
    end

    @testset "Value semantics for mode sets and cache keys" begin
        setprecision(BigFloat, 256) do
            a, b = HarmonicMode(2, 1, big"0.75"), HarmonicMode(2, 1, big"0.75")
            @test a == b
            @test isequal(a, b) && hash(a) == hash(b)
            @test haskey(Dict(a => 1), b)
            @test length(unique([a, b])) == 1
            @test a != HarmonicMode(2, 0, big"0.75") && a != HarmonicMode(1, 1, big"0.75")
            @test a != HarmonicMode(2, 1, big"0.5")
            @test HarmonicMode(2, 1, 0.75) == HarmonicMode{BigFloat}(2, 1, 0.75)
            @test hash(HarmonicMode(2, 1, 0.75)) == hash(HarmonicMode{BigFloat}(2, 1, 0.75))
            for policy in (FixedSoftening, SemimajorScaleHeight,
                           InstantaneousScaleHeight, LocalScaleHeight)
                p, q = policy(big"0.3"), policy(big"0.3")
                @test p == q
                @test isequal(p, q) && hash(p) == hash(q)
                @test haskey(Dict(p => 1), q)
                @test p != policy(big"0.5")
            end
        end
        @test HarmonicMode(1, -1, 0.0) == HarmonicMode(1, -1, -0.0)
        # Equal scale, different policy: distinct keys.
        @test length(unique((FixedSoftening(0.3), SemimajorScaleHeight(0.3),
                             InstantaneousScaleHeight(0.3), LocalScaleHeight(0.3)))) == 4
    end

    @testset "Precision inherited at 128 and 256 bits" begin
        for bits in (128, 256)
            setprecision(BigFloat, bits) do
                n = 1 + eps(BigFloat)
                mode = NewtonianMode(2, 3, n)
                @test mode isa HarmonicMode{BigFloat}
                @test mode.omega == 3n
                @test mode.omega != BigFloat(Float64(mode.omega))
                @test mode_power(mode, 4) == 6n
                @test mode_power(HarmonicMode(2, 0, 3.0), BigFloat(4)) isa BigFloat
                @test mode_power(HarmonicMode(3, 0, 1.0), BigFloat(1)) == BigFloat(1) / 3
                @test precision(mode.omega) == bits
                @test precision(mode_power(mode, 4)) == bits
                @test HarmonicMode{BigFloat}(2, 1, mode.omega).omega == mode.omega
                b = BigFloat(3) / 10
                @test FixedSoftening(b * n).eps == b * n
                for policy in (SemimajorScaleHeight, InstantaneousScaleHeight, LocalScaleHeight)
                    @test policy(b) isa SofteningPolicy{BigFloat}
                    @test policy(b).b == b
                    @test precision(policy(b).b) == bits
                end
                @test precision(BigFloat) == bits
            end
        end
    end
end

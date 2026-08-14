# Separatrix/separatrix.jl — grounded ANALYTIC + BAKED (KerrBackwards/FEW/WaSABI are HEAVY → NO
# live dep): the Stein–Warburton polynomials are re-stated here verbatim from
# KerrBackwards/src/utility.jl and the canonical `separatrix` is asserted to return their
# bracketed root (residual ≈ 0). Plus closed-form ISCO/Schwarzschild invariants.
using Test
using Spacetime

# Stein–Warburton 2020 polynomials, transcribed verbatim from KerrBackwards/src/utility.jl
# (lines 37–66) as the BAKED root oracle. p_sep must zero the relevant polynomial.
_oracle_equat(p, a, e) = a^4*(-3 - 2e + e^2)^2 + p^2*(-6 - 2e + p)^2 -
    2*a^2*(1 + e)*p*(14 + 2*e^2 + 3*p - e*p)
_oracle_polar(p, a, e) = a^6*(-1 + e)^2*(1 + e)^4 + p^5*(-6 - 2e + p) +
    a^2*p^3*(-4*(-1 + e)*(1 + e)^2 + (3 + e*(2 + 3e))*p) -
    a^4*(1 + e)^2*p*(6 + 2*e^3 + 2*e*(-1 + p) - 3*p - 3*e^2*(2 + p))

@testset "Separatrix + ISCO" begin
    @testset "a=0 Schwarzschild: separatrix == 6 + 2e exactly (all e)" begin
        for e in (0.0, 0.1, 0.3, 0.5, 0.7)
            for x in (1.0, -1.0, 0.0, 0.5)
                @test separatrix(0.0, e, x) == 6 + 2e
            end
        end
    end

    @testset "kerr_isco Bardeen limits: 6 / 1 / 9" begin
        @test kerr_isco(0.0) == 6
        @test isapprox(kerr_isco(1.0),  1.0; atol = 1e-12)   # prograde extremal
        @test isapprox(kerr_isco(-1.0), 9.0; atol = 1e-12)   # retrograde extremal
        @test kerr_isco(0.5) < 6 && kerr_isco(0.5) > 1        # prograde < Schwarzschild
        @test kerr_isco(-0.5) > 6 && kerr_isco(-0.5) < 9      # retrograde > Schwarzschild
    end

    @testset "separatrix(a,0,±1) == kerr_isco(±a): ISCO–separatrix consistency" begin
        for a in (0.3, 0.7, 0.9, 0.998)
            @test separatrix(a, 0.0, 1.0)  == kerr_isco(a)
            @test separatrix(a, 0.0, -1.0) == kerr_isco(-a)
        end
    end

    @testset "equatorial separatrix is a root of the baked Stein–Warburton polynomial" begin
        for a in (0.3, 0.7, 0.9, 0.998), e in (0.05, 0.2, 0.4)
            ps = separatrix(a, e, 1.0)
            @test isapprox(_oracle_equat(ps, a, e), 0.0; atol = 1e-7)
            @test 1 + e <= ps <= 6 + 2e                      # correct bracket
            psr = separatrix(a, e, -1.0)                      # retrograde equatorial
            @test isapprox(_oracle_equat(psr, a, e), 0.0; atol = 1e-7)
        end
    end

    @testset "polar separatrix is a root of the baked polar polynomial" begin
        for a in (0.3, 0.7, 0.9), e in (0.05, 0.2, 0.4)
            ps = separatrix(a, e, 0.0)
            @test isapprox(_oracle_polar(ps, a, e), 0.0; atol = 1e-6)
        end
    end

    @testset "generic separatrix bracketed by equatorial ≤ p_sep ≤ polar (prograde)" begin
        for a in (0.5, 0.9), e in (0.1, 0.3), x in (0.3, 0.6, 0.9)
            equat = separatrix(a, e, 1.0)
            polar = separatrix(a, e, 0.0)
            ps    = separatrix(a, e, x)
            @test equat <= ps <= polar
        end
    end

    @testset "get_separatrix is a thin alias of separatrix (FEW convention)" begin
        for (a, e, x) in ((0.9, 0.2, 0.7), (0.5, 0.0, 1.0), (0.0, 0.3, 0.5))
            @test get_separatrix(a, e, x) == separatrix(a, e, x)
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        @test separatrix(big"0.0", big"0.3", big"1.0") isa BigFloat
        @test separatrix(big"0.9", big"0.0", big"1.0") isa BigFloat   # ISCO branch
        ps = separatrix(big"0.9", big"0.2", big"1.0")                 # A42 root branch
        @test ps isa BigFloat
        @test abs(_oracle_equat(ps, big"0.9", big"0.2")) < 1e-20      # BigFloat-precision root
    end
end

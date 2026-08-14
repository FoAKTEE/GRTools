# Characterization of the vendored 1D root-finders.
#
# Grounding is ANALYTIC ONLY (no light legacy oracle — the separatrix/_brentq solvers live in heavy
# FEW). Expected roots are independent closed forms / known constants, never the solver's own output:
#   - cos(x)=x        → Dottie number 0.7390851332151607 (known transcendental constant)
#   - x²−2=0          → sqrt(2)   (hardware/library sqrt, independent of the root-finders)
#   - x³−2=0          → cbrt(2)   (independent builtin)
#   - eˣ−2=0          → log(2)    (independent builtin)
#   - x²−6x+8=0       → exact integer roots 2, 4
# Invariants asserted: bracketing methods require f(a)·f(b)<0 (throw otherwise); secant converges
# superlinearly; A42 reaches the root in few iterations; BigFloat secant resolves √2 to 2^-200.

using NumericsCore
using NumericsCore: secant, bisect, brent, a42

const DOTTIE = 0.7390851332151607   # cos(x)=x

@testset "rootfind" begin
    @testset "analytic roots — every solver hits the closed-form value" begin
        # (f, bracket lo, bracket hi, two secant seeds, expected root)
        cases = [
            (x -> cos(x) - x,    0.0, 1.0, 0.0, 1.0, DOTTIE),
            (x -> x^2 - 2,       1.0, 2.0, 1.0, 2.0, sqrt(2.0)),
            (x -> x^3 - 2,       1.0, 2.0, 1.0, 2.0, cbrt(2.0)),
            (x -> exp(x) - 2,    0.0, 1.0, 0.0, 1.0, log(2.0)),
            (x -> x^2 - 6x + 8,  3.0, 5.0, 3.5, 5.0, 4.0),
            (x -> x^2 - 6x + 8,  1.0, 3.0, 1.0, 2.5, 2.0),
        ]
        for (f, lo, hi, s0, s1, root) in cases
            @test secant(f, s0, s1) ≈ root atol=1e-12
            @test bisect(f, lo, hi; atol=1e-13) ≈ root atol=1e-12
            @test brent(f, lo, hi) ≈ root atol=1e-12
            @test a42(f, lo, hi) ≈ root atol=1e-12
            # residual sanity
            @test abs(f(brent(f, lo, hi))) < 1e-11
            @test abs(f(a42(f, lo, hi))) < 1e-11

            # Bracketing front doors normalize reversed bounds to the same root.
            @test bisect(f, hi, lo; atol=1e-13) ≈ root atol=1e-12
            @test brent(f, hi, lo) ≈ root atol=1e-12
            @test a42(f, hi, lo) ≈ root atol=1e-12
        end
    end

    @testset "bracketing invariant — non-bracketing input is rejected" begin
        g = x -> x^2 - 2          # positive on both ends of [2,3]
        @test_throws ArgumentError bisect(g, 2.0, 3.0)
        @test_throws ArgumentError brent(g, 2.0, 3.0)
        @test_throws ArgumentError a42(g, 2.0, 3.0)
    end

    @testset "exhausted iteration budgets fail rather than return non-roots" begin
        f = x -> x^2 - 2
        @test_throws ErrorException secant(f, 1.0, 2.0; rtol=0.0, atol=0.0, maxiter=1)
        @test_throws ErrorException bisect(f, 1.0, 2.0; atol=0.0, maxiter=1)
        @test_throws ErrorException brent(f, 1.0, 2.0; rtol=0.0, atol=0.0, maxiter=1)
        @test_throws ErrorException a42(f, 1.0, 2.0; rtol=0.0, atol=0.0, maxiter=1)

        # Exact seeds/endpoints remain successful even with a zero iteration budget.
        g = x -> x - 1
        @test secant(g, 1.0, 2.0; maxiter=0) == 1.0
        @test bisect(g, 1.0, 2.0; maxiter=0) == 1.0
        @test brent(g, 2.0, 1.0; maxiter=0) == 1.0
        @test a42(g, 2.0, 1.0; maxiter=0) == 1.0
    end

    @testset "A42 high-order: converges in few iterations" begin
        f = x -> x^3 - 2
        r = a42(f, 1.0, 2.0; maxiter=8)
        @test r ≈ cbrt(2.0) atol=1e-13
    end

    @testset "BigFloat precision inheritance — √2 to 2^-200" begin
        setprecision(BigFloat, 256) do
            ref = sqrt(big(2))
            f = x -> x^2 - 2
            rs = secant(f, big(1.0), big(2.0); atol=big(2.0)^-220, rtol=big(2.0)^-220, maxiter=200)
            @test rs isa BigFloat
            @test abs(rs - ref) < big(2.0)^-200
            rb = a42(f, big(1.0), big(2.0); maxiter=200)
            @test abs(rb - ref) < big(2.0)^-200
            rbr = brent(f, big(1.0), big(2.0); maxiter=200)
            @test abs(rbr - ref) < big(2.0)^-200
        end
    end
end

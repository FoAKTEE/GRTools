# Characterization of the shared continued-fraction primitives.
#
# Grounding is ANALYTIC closed-form Gauss-CF identities (computed to BigFloat precision), plus a
# captured-literal regression of the Teukolsky K-form behavior (read from
# packages/Teukolsky/src/continued_fraction.jl — dep-free body, but Teukolsky pulls HypGeo/SWSH so it
# is NOT live-dep'd; the no-leading-term + A/B-seed contract is baked as an exact rational instead):
#   - golden ratio:  a≡1, b≡1, n0=1  →  1/(1+1/(1+⋯)) = (√5−1)/2 = φ−1
#   - √2−1:          a≡1, b≡2, n0=1  →  1/(2+1/(2+⋯)) = √2−1
#   - tan z:         a(1)=z, a(n≥2)=−z², b(n)=2n−1     →  tan z   (Gauss CF, transcendental)
#   - tanh z:        a(1)=z, a(n≥2)=+z², b(n)=2n−1     →  tanh z  (hyperbolic Gauss CF)
#   - cf_finite no-leading-term: nmax=0 ⇒ f(n)/g(n); exact rational cf_finite(1,2,1,3)=12/29
# cf_lentz is asserted to reproduce cf_convergent (same K-form, independent iteration).

using NumericsCore
using NumericsCore: cf_finite, cf_convergent, cf_lentz

# Gauss-CF coefficient closures
golden_a(_) = 1.0;  golden_b(_) = 1.0
sqrt2_a(_) = 1.0;   sqrt2_b(_) = 2.0
tan_a(z)  = n -> n == 1 ? z : -z^2
tanh_a(z) = n -> n == 1 ? z :  z^2
odd_b(n) = 2n - 1.0

@testset "continued_fraction" begin
    @testset "analytic Gauss-CF identities (cf_convergent)" begin
        φm1 = (sqrt(5.0) - 1) / 2
        @test cf_convergent(golden_a, golden_b, 1) ≈ φm1 rtol=1e-14
        @test cf_convergent(sqrt2_a, sqrt2_b, 1) ≈ sqrt(2.0) - 1 rtol=1e-14
        for z in (0.3, 0.5, 0.9, 1.2)
            @test cf_convergent(tan_a(z), odd_b, 1) ≈ tan(z) rtol=1e-13
            @test cf_convergent(tanh_a(z), odd_b, 1) ≈ tanh(z) rtol=1e-13
        end
    end

    @testset "modified Lentz reproduces the A/B recurrence (same K-form)" begin
        @test cf_lentz(golden_a, golden_b, 1) ≈ cf_convergent(golden_a, golden_b, 1) rtol=1e-13
        @test cf_lentz(sqrt2_a, sqrt2_b, 1) ≈ cf_convergent(sqrt2_a, sqrt2_b, 1) rtol=1e-13
        for z in (0.3, 0.5, 0.9, 1.2)
            @test cf_lentz(tan_a(z), odd_b, 1) ≈ tan(z) rtol=1e-13
            @test cf_lentz(tanh_a(z), odd_b, 1) ≈ tanh(z) rtol=1e-13
        end
    end

    @testset "cf_finite K-form contract (Teukolsky source, captured literal)" begin
        # no leading term: nmax=0 ⇒ f(n)/g(n)
        f = i -> Float64(i^2); g = i -> Float64(2i + 1)
        @test cf_finite(f, g, 5, 0) == f(5) / g(5)
        # exact rational: f≡1, g≡2, n=1, nmax=3 → 12/29 (hand-evaluated K-form, no leading term)
        @test cf_finite(i -> 1.0, i -> 2.0, 1, 3) ≈ 12 / 29 rtol=1e-15
        # finite truncation converges to the convergent CF (√2−1) as nmax grows
        @test cf_finite(i -> 1.0, i -> 2.0, 1, 60) ≈ sqrt(2.0) - 1 rtol=1e-13
        # finite tan-CF truncation approaches tan z
        @test cf_finite(tan_a(0.5), odd_b, 1, 40) ≈ tan(0.5) rtol=1e-12
    end

    @testset "BigFloat precision inheritance" begin
        setprecision(BigFloat, 256) do
            φm1 = (sqrt(big(5)) - 1) / 2
            got = cf_convergent(n -> big(1.0), n -> big(1.0), 1; rtol=big(2.0)^-240)
            @test got isa BigFloat
            @test abs(got - φm1) < big(2.0)^-200
            s2 = cf_lentz(n -> big(1.0), n -> big(2.0), 1; rtol=big(2.0)^-240)
            @test abs(s2 - (sqrt(big(2)) - 1)) < big(2.0)^-200
        end
    end
end

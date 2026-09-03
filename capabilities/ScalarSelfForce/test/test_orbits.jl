# test_orbits.jl — ESSF node `orbits-adapters` (schedule W2-C), written BEFORE the
# implementation (TDD). Gates (GLOBAL_DAG #24 oracle column + node contract):
#   * characterization vs BAKED Spacetime fixtures (test/fixtures/orbits_circular.jl,
#     generation recipe in its header), rtol 1e-12, grid a ∈ {0, ±0.5, ±0.9, ±0.998} ×
#     r0 ∈ {isco, 1.001·isco, 4, 6, 10};
#   * puncture-scale constants z_c² and ϱ²-prefactors vs main.tex:445-455 literal formulas;
#   * ISCO closed-form cross-gates (main.tex:207-209 condition) at ±{0, 0.5, 0.9, 0.998};
#   * behavior AT the exact ISCO pinned (upstream Spacetime.kerr_frequencies circ-eq
#     DomainErrors there — the adapter must avoid that path; error-ledger row 2026-07-28,
#     node ESSF::fix-kerr-frequencies-isco);
#   * circularity dV_eff/dr = 0 at r0 to 1e-12 (normalized). NOTE: the V_eff printed in
#     main.tex:190 has a LAST-TERM SIGN TYPO (+M(L−aE)²/r_p³); its derivative does NOT vanish
#     on the paper's own E, L solutions (measured ≈0.5 normalized). The correct condition,
#     from (u^r)² = E²−1 + 2M/r − (L²+a²(1−E²))/r² + 2M(L−aE)²/r³, is
#         M/r² − (L²+a²(1−E²))/r³ + 3M(L−aE)²/r⁴ = 0,
#     which the captured grid satisfies to ≤2.2e-15 (see the capture recipe's diagnostics);
#   * explicit retrograde (a < 0) sanity; type-genericity at BigFloat(256).
using Test
using ScalarSelfForce
using Spacetime: kerr_isco

include(joinpath(@__DIR__, "fixtures", "orbits_circular.jl"))

# corrected circularity residual (see header note), normalized by the term-magnitude sum
function dveff_dr_normalized(a, r0, E, L)
    t1 = 1 / r0^2
    t2 = -(L^2 + a^2 * (1 - E^2)) / r0^3
    t3 = 3 * (L - a * E)^2 / r0^4
    abs(t1 + t2 + t3) / (abs(t1) + abs(t2) + abs(t3))
end

const ORBIT_SPINS = [0.0, 0.5, -0.5, 0.9, -0.9, 0.998, -0.998]

@testset "Orbits: CircularOrbit{T} adapter over Spacetime (orbits-adapters)" begin

    @testset "characterization vs baked Spacetime fixtures (rtol 1e-12)" begin
        for fx in ORBITS_CIRCULAR_FIXTURES
            o = CircularOrbit(fx.a, fx.r0)
            @test o isa CircularOrbit{Float64}
            @test o.M == 1.0
            @test o.a == fx.a
            @test o.r0 == fx.r0
            @test o.prograde == (fx.a >= 0)
            @test isapprox(o.E, fx.E; rtol = 1e-12)
            @test isapprox(o.Lz, fx.Lz; rtol = 1e-12)
            @test isapprox(o.ut, fx.ut; rtol = 1e-12)
            @test isapprox(o.uphi, fx.uphi; rtol = 1e-12)
            @test isapprox(o.Omega, fx.Omega; rtol = 1e-12)
            @test isapprox(o.rp, fx.rp; rtol = 1e-12)
            @test isapprox(o.rm, fx.rm; rtol = 1e-12)   # rm = 0 exactly at a = 0 on both sides
            @test isapprox(o.r_isco, fx.risco; rtol = 1e-12)
        end
    end

    @testset "puncture-scale constants vs main.tex:445-455 literals" begin
        for fx in ORBITS_CIRCULAR_FIXTURES
            o = CircularOrbit(fx.a, fx.r0)
            a, r0, L = fx.a, fx.r0, fx.Lz
            # main.tex:452  ϱ² = r_p²/(r_p² − 2M r_p + a²)·Δr² + r_p²·Δθ²   (M = 1)
            @test isapprox(o.rho2_dr2, r0^2 / (r0^2 - 2 * r0 + a^2); rtol = 1e-12)
            @test o.rho2_dtheta2 == r0^2
            # main.tex:453  z_c² = 4(L² + r_p² + a² + 2M a²/r_p)
            @test isapprox(o.zc2, 4 * (L^2 + r0^2 + a^2 + 2 * a^2 / r0); rtol = 1e-12)
        end
    end

    @testset "ISCO closed-form cross-gates (main.tex:207-209)" begin
        for a in ORBIT_SPINS
            r = kerr_isco(a)
            # paper ISCO condition: 0 = 1 − 6M/r + 8a√(M/r³) − 3a²/r²
            @test abs(1 - 6 / r + 8 * a * sqrt(1 / r^3) - 3 * a^2 / r^2) < 1e-12
            # the struct caches exactly Spacetime's ISCO
            @test CircularOrbit(a, 10.0).r_isco == r
        end
        # Bardeen limits (exact in Float64)
        @test isapprox(kerr_isco(0.0), 6.0; rtol = 4 * eps())
        @test isapprox(kerr_isco(1.0), 1.0; rtol = 4 * eps())
        @test isapprox(kerr_isco(-1.0), 9.0; rtol = 4 * eps())
    end

    @testset "pinned AT the exact ISCO" begin
        # The ISCO is a regular point of the upstream frequencies now: Spacetime's circ-eq Υ_r
        # radicand is the Bardeen ISCO polynomial, so Υ_r ≡ 0 there for prograde orbits and the
        # old DomainError at r0 = kerr_isco(a) is gone (node ESSF::fix-kerr-frequencies-isco;
        # the 2026-07-28 error-ledger trial is discharged). The adapter still takes Ω from the
        # closed form Ω = u^φ/u^t = 1/(a + √(r0³/M)) (main.tex:196) by design, and must
        # construct fine at the exact ISCO.
        for a in ORBIT_SPINS
            r0 = kerr_isco(a)
            o = CircularOrbit(a, r0)                       # must not throw
            @test isfinite(o.E) && isfinite(o.Lz) && isfinite(o.Omega)
            @test isapprox(o.Omega, 1 / (a + sqrt(r0^3)); rtol = 1e-12)
            @test isapprox(o.Omega, o.uphi / o.ut; rtol = 4 * eps())
            @test o.E < 1                                  # ISCO orbit is bound
        end
    end

    @testset "circular-orbit physics identities (main.tex:176-196)" begin
        for fx in ORBITS_CIRCULAR_FIXTURES
            o = CircularOrbit(fx.a, fx.r0)
            a, r0 = fx.a, fx.r0
            # E = (1 − 2M/r)u^t + (2Ma/r)u^φ ;  L = (r² + a² + 2Ma²/r)u^φ − (2Ma/r)u^t
            @test isapprox(o.E, (1 - 2 / r0) * o.ut + (2 * a / r0) * o.uphi; rtol = 1e-10)
            @test isapprox(o.Lz, (r0^2 + a^2 + 2 * a^2 / r0) * o.uphi - (2 * a / r0) * o.ut;
                           rtol = 1e-10)
            # circularity: dV_eff/dr = 0 at r0 (corrected sign — see file header)
            @test dveff_dr_normalized(a, r0, o.E, o.Lz) < 1e-12
            # Ω u^t = u^φ exactly by construction
            @test isapprox(o.Omega * o.ut, o.uphi; rtol = 4 * eps())
            # stable branch (r0 ≥ isco): bound orbit
            if fx.r0 >= fx.risco
                @test o.E < 1
            end
        end
    end

    @testset "retrograde (a < 0) explicit" begin
        for apos in (0.5, 0.9, 0.998)
            op = CircularOrbit(apos, 10.0)
            or = CircularOrbit(-apos, 10.0)
            @test op.prograde && !or.prograde
            # sign convention: motion is ALWAYS +φ; the spin axis is mirrored
            @test or.uphi > 0 && or.Omega > 0 && or.Lz > 0
            # frame dragging: at fixed r0, Ω = 1/(a + r0^{3/2}) is LARGER retrograde
            @test or.Omega > op.Omega
            # retrograde orbits are less bound and carry more angular momentum
            @test or.E > op.E
            @test or.Lz > op.Lz
            # retrograde ISCO lies further out
            @test or.r_isco > op.r_isco
        end
        @test CircularOrbit(0.0, 6.0).prograde            # a = 0 counts as prograde
    end

    @testset "type-genericity: BigFloat(256) inheritance" begin
        setprecision(BigFloat, 256) do
            a = BigFloat(1) / 2
            r0 = BigFloat(10)
            o = CircularOrbit(a, r0)
            @test o isa CircularOrbit{BigFloat}
            for f in (:M, :a, :r0, :E, :Lz, :Omega, :ut, :uphi, :rp, :rm, :r_isco,
                      :rho2_dr2, :rho2_dtheta2, :zc2)
                v = getfield(o, f)
                @test v isa BigFloat
                @test precision(v) == 256
            end
            # Float64 fidelity vs the baked fixture row (a = 0.5, r0 = 10)
            fx = only(filter(r -> r.a == 0.5 && r.r0 == 10.0, ORBITS_CIRCULAR_FIXTURES))
            @test isapprox(Float64(o.E), fx.E; rtol = 1e-13)
            @test isapprox(Float64(o.Lz), fx.Lz; rtol = 1e-13)
            @test isapprox(Float64(o.Omega), fx.Omega; rtol = 1e-13)
            @test isapprox(Float64(o.ut), fx.ut; rtol = 1e-13)
            # genuine precision inheritance: identity residuals at the BigFloat floor,
            # far below any Float64-contaminated path (eps(BigFloat, 256) ≈ 1.7e-77)
            @test abs(o.E - ((1 - 2 / r0) * o.ut + (2 * a / r0) * o.uphi)) < big"1e-70"
            @test abs(o.Lz - ((r0^2 + a^2 + 2 * a^2 / r0) * o.uphi - (2 * a / r0) * o.ut)) < big"1e-70"
            @test dveff_dr_normalized(a, r0, o.E, o.Lz) < big"1e-70"
            # mixed-precision arguments promote to BigFloat
            @test CircularOrbit(0.5, BigFloat(10)) isa CircularOrbit{BigFloat}
        end
    end
end

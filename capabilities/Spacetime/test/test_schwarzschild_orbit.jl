# Geometry/schwarzschild_orbit.jl — PURE ANALYTIC grounding (H1Lorenz/DiscModes are HEAVY —
# never a live test-dep; source closed forms baked as expected literals) plus the in-package
# Geodesics consistency identity (a=0 reduction of constants_of_motion / kerr_frequencies).
using Test
using Spacetime

@testset "Geometry: Schwarzschild orbit / tortoise / potential" begin
    @testset "lapse f = 1 − 2/r and f' = 2/r² (closed form + FD)" begin
        for r in (3.0, 5.0, 10.0, 50.0)
            @test schwarzschild_lapse(r) ≈ 1 - 2 / r rtol = 1e-15
            @test schwarzschild_dlapse(r) ≈ 2 / r^2  rtol = 1e-15
            # finite-difference cross-check of the analytic derivative
            h = 1e-6
            fd = (schwarzschild_lapse(r + h) - schwarzschild_lapse(r - h)) / (2h)
            @test schwarzschild_dlapse(r) ≈ fd rtol = 1e-7
        end
        @test schwarzschild_lapse(2.0) == 0.0          # horizon
    end

    @testset "BLS potential V_l = ¼ f (f'/r + l(l+1)/r²)" begin
        for (r, l) in ((6.0, 2), (10.0, 3), (20.0, 0))
            f = 1 - 2 / r; df = 2 / r^2
            @test bls_potential(r, l) ≈ 0.25 * f * (df / r + l * (l + 1) / r^2) rtol = 1e-14
        end
    end

    @testset "tortoise r*(r) = r + 2 ln(r/2 − 1) and its inverse" begin
        # Closed-form Lambert-branch reference: pick y = r/2−1 > 0 ⇒ r = 2(y+1),
        # r*(r) = 2y + 2 + 2 ln y  (independent of the implementation's Newton solver).
        for y in (0.05, 0.5, 1.0, 4.0, 20.0)
            r  = 2 * (y + 1)
            rs = 2y + 2 + 2 * log(y)
            @test tortoise_schw(r) ≈ rs rtol = 1e-13
            @test tortoise_schw_inverse(rs) ≈ r rtol = 1e-12
        end
        # round-trip
        for r in (2.001, 2.5, 4.0, 10.0, 100.0)
            @test tortoise_schw_inverse(tortoise_schw(r)) ≈ r rtol = 1e-11
        end
        # The old rs<-50 shortcut returned the horizon itself, where the forward map is undefined.
        for rs in (-60.0, -100.0)
            r = tortoise_schw_inverse(rs)
            @test r > 2.0
            @test isfinite(tortoise_schw(r))
        end
        @test tortoise_schw(tortoise_schw_inverse(-60.0)) ≈ -60.0 atol = 0.02
    end

    @testset "circular-orbit constants: closed form + a=0 Geodesics reduction" begin
        for p in (6.5, 8.0, 10.0, 20.0)
            @test schw_orbit_energy(p)  ≈ (1 - 2 / p) / sqrt(1 - 3 / p) rtol = 1e-14
            @test schw_orbit_ang_mom(p) ≈ sqrt(p) / sqrt(1 - 3 / p)     rtol = 1e-14
            @test schw_orbit_phi_freq(p) ≈ p^(-1.5)                     rtol = 1e-14
            # consistency identity: reduce to the in-package generic Geodesics at a=0,e=0,x=1.
            @test schw_orbit_energy(p)   ≈ geodesic_energy(0.0, p, 0.0, 1.0)           rtol = 1e-13
            @test schw_orbit_ang_mom(p)  ≈ geodesic_angular_momentum(0.0, p, 0.0, 1.0) rtol = 1e-13
            @test schw_orbit_phi_freq(p) ≈ kerr_frequencies(0.0, p, 0.0, 1.0).Omegaph  rtol = 1e-13
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        pb = big"10.0"
        @test schw_orbit_energy(pb) isa BigFloat
        @test schw_orbit_phi_freq(pb) ≈ pb^(big(-3) // 2) rtol = big"1e-60"
        rb = big"10.0"
        @test tortoise_schw_inverse(tortoise_schw(rb)) ≈ rb rtol = big"1e-40"
        rdeep = tortoise_schw_inverse(big"-100.0")
        @test rdeep > 2
        # Forming r=2(1+y) at y≈e^-51 loses ~22 decimal digits to the final addition; the
        # cancellation-free solve still retains ~56 reliable decimal digits in the round trip.
        @test tortoise_schw(rdeep) ≈ big"-100.0" rtol = big"1e-56"
    end
end

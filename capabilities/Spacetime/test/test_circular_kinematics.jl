# Geodesics/circular_kinematics.jl — W1-C (elliptic-SSF, ESSF::spacetime-circular-kinematics).
# Grounded by LITERAL transcription of the CRP closed forms (ref-paper/arxiv-2606.04998/src/
# main.tex:192-196 u^t/u^φ/Ω; :223-226 φ̄−φ chart offset) plus coordinate-free gates: the
# metric normalization g_{μν}u^μu^ν = −1 and the radial geodesic equation Γ^r_{μν}u^μu^ν = 0
# built from Spacetime's OWN kerr_metric/christoffel — both valid for prograde AND retrograde
# (a < 0) spins, so the a-sign convention is pinned by physics, not by the same formula twice.
using Test
using Spacetime

@testset "CircularKinematics: kerr_circular_ut_uphi + kerr_varphi_offset(_dr)" begin
    # paper convention (main.tex:196-198): sign of a selects prograde/retrograde, u^φ > 0 always
    spins = (0.0, 0.5, -0.5, 0.998, -0.998)
    radii(a) = (kerr_isco(a), 1.001 * kerr_isco(a), 10.0, 30.0)   # incl. exact + near-ISCO

    @testset "closed forms vs main.tex:192-196 literal transcription (M=1)" begin
        for a in spins, r0 in radii(a)
            ut, uphi = kerr_circular_ut_uphi(a, r0)
            v = sqrt(r0^3)                                        # √(r0³/M), M = 1
            den = sqrt(r0^3 - 3 * r0^2 + 2 * a * v)
            @test isapprox(ut,   (a + v) / den; rtol = 1e-14)
            @test isapprox(uphi, 1 / den;       rtol = 1e-14)
            @test uphi > 0                                        # always +φ motion
            @test ut > 0
            @test isapprox(uphi / ut, 1 / (a + v); rtol = 1e-13)  # Ω = 1/(a + √(r0³/M))
        end
    end

    @testset "normalization g_{μν}u^μu^ν = -1 (Spacetime's own Kerr metric, incl. a<0, ISCO)" begin
        for a in spins, r0 in radii(a)
            ut, uphi = kerr_circular_ut_uphi(a, r0)
            g = kerr_metric(a, r0, π / 2)
            guu = g[1, 1] * ut^2 + 2 * g[1, 4] * ut * uphi + g[4, 4] * uphi^2
            @test isapprox(guu, -1.0; atol = 1e-13)
        end
    end

    @testset "radial geodesic equation Γ^r_{μν}u^μu^ν = 0 (christoffel cross-check, incl. a<0)" begin
        for a in spins, r0 in radii(a)
            ut, uphi = kerr_circular_ut_uphi(a, r0)
            Γ = christoffel(a, r0, π / 2)                         # Γ^λ_{μν}, λ first
            terms = (Γ[2, 1, 1] * ut^2, 2 * Γ[2, 1, 4] * ut * uphi, Γ[2, 4, 4] * uphi^2)
            scale = sum(abs, terms)
            @test abs(sum(terms)) <= 1e-12 * scale
        end
    end

    @testset "Ω and E/Lz reuse cross-checks (kerr_frequencies, constants_of_motion)" begin
        for a in (0.0, 0.5, 0.998), r0 in radii(a)
            ut, uphi = kerr_circular_ut_uphi(a, r0)
            if r0 != kerr_isco(a)   # kerr_frequencies circ-eq Υr radicand → 0⁻ by roundoff AT isco
                f = kerr_frequencies(a, r0, 0.0, 1.0)             # circ-eq closed-form path
                @test isapprox(uphi / ut, f.Omegaph; rtol = 1e-12)
            end
            E, Lz, _ = constants_of_motion(a, r0, 0.0, 1.0)
            g = kerr_metric(a, r0, π / 2)
            @test isapprox(g[1, 1] * ut + g[1, 4] * uphi, -E; rtol = 1e-12)  # u_t = -E
            @test isapprox(g[4, 1] * ut + g[4, 4] * uphi,  Lz; rtol = 1e-12) # u_φ = Lz
        end
    end

    @testset "φ̄-φ offset: main.tex:223-226 literal + dφ̄/dr = a/Δ vs central FD" begin
        for a in (0.5, -0.5, 0.9, 0.998, -0.998), r in (4.0, 6.0, 10.0)
            rp, rm = horizons(a)
            @test isapprox(kerr_varphi_offset(a, r),
                           a / (rp - rm) * log((r - rp) / (r - rm)); rtol = 1e-14)
            # offset has the opposite sign of a on the exterior (log argument < 1)
            @test sign(kerr_varphi_offset(a, r)) == -sign(a)
            # exported derivative = a/Δ, and matches the offset's central finite difference
            # (4th-order 5-point stencil: the a/(r₊−r₋) prefactor amplifies log roundoff, putting
            # the 2nd-order-stencil floor at ~8e-10 relative near r=4, a=±0.998 — measured)
            @test isapprox(kerr_varphi_offset_dr(a, r), a / delta_kerr(a, r); rtol = 1e-14)
            h = 1e-3
            f(x) = kerr_varphi_offset(a, x)
            fd = (f(r - 2h) - 8 * f(r - h) + 8 * f(r + h) - f(r + 2h)) / (12h)
            @test isapprox(fd, kerr_varphi_offset_dr(a, r); rtol = 1e-9)
        end
    end

    @testset "a=0 degeneracy: zero offset, Schwarzschild u^t/u^φ forms" begin
        for r0 in (6.0, 10.0, 30.0)
            @test kerr_varphi_offset(0.0, r0) == 0
            @test kerr_varphi_offset_dr(0.0, r0) == 0
            ut, uphi = kerr_circular_ut_uphi(0.0, r0)
            @test isapprox(ut, 1 / sqrt(1 - 3 / r0);                 rtol = 1e-14)
            @test isapprox(uphi, r0^(-3 // 2) / sqrt(1 - 3 / r0);    rtol = 1e-14)
            @test isapprox(uphi / ut, schw_orbit_phi_freq(r0);       rtol = 1e-14)
            g = kerr_metric(0.0, r0, π / 2)                          # u_t=-E, u_φ=L (Schw)
            @test isapprox(g[1, 1] * ut, -schw_orbit_energy(r0);     rtol = 1e-13)
            @test isapprox(g[4, 4] * uphi, schw_orbit_ang_mom(r0);   rtol = 1e-13)
        end
    end

    @testset "type-genericity: BigFloat(256) pass-through, Int promotion" begin
        setprecision(BigFloat, 256) do
            a, r0 = big"0.5", big"10.0"
            ut, uphi = kerr_circular_ut_uphi(a, r0)
            @test ut isa BigFloat && precision(ut) == 256
            g = kerr_metric(a, r0, BigFloat(π) / 2)
            guu = g[1, 1] * ut^2 + 2 * g[1, 4] * ut * uphi + g[4, 4] * uphi^2
            @test abs(guu + 1) < big"1e-70"                       # impossible at Float64
            off = kerr_varphi_offset(big"0.7", big"8.0")
            @test off isa BigFloat && precision(off) == 256
            h = big"1e-20"                                        # FD at BigFloat precision
            fd = (kerr_varphi_offset(big"0.7", big"8.0" + h) -
                  kerr_varphi_offset(big"0.7", big"8.0" - h)) / (2h)
            @test abs(fd - kerr_varphi_offset_dr(big"0.7", big"8.0")) < big"1e-35"
        end
        @test kerr_circular_ut_uphi(0, 10).ut isa Float64          # Int inputs promote
        @test kerr_varphi_offset(0, 10) == 0
    end
end

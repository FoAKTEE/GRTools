# Geometry/equatorial_potentials.jl — EGD::spacetime-kerr-equatorial-potentials (arXiv:2608.17003).
# The Page–Thorne-form potentials (ν, ψ, μ, ω̃) of the Kerr equatorial line element
#     ds² = −e^{2ν}dt² + e^{2ψ}(dφ − ω̃ dt)² + e^{2μ}dr² + dz²          (main.tex:148, eq (4))
# are grounded THREE independent ways: (i) coordinate-free identities against Spacetime's OWN
# `kerr_metric(a, r, π/2)` (e^{2ψ}=g_φφ, ω̃=−g_tφ/g_φφ, e^{2ν}=−g_tt+g_tφ²/g_φφ, e^{2μ}=g_rr);
# (ii) the LITERAL appendix transcription (appendix.tex:7-21, eqs A1a–e, M=1, χ=a); (iii) the
# a=0 Schwarzschild limits. The analytic r-derivatives are gated by central finite differences at
# Float64 AND BigFloat(256); the circular dt/dτ (main.tex:673) is pinned to the existing
# `kerr_circular_ut_uphi` (REUSED, not re-derived) and to the paper's e^{-ν}/√(1−v²) decomposition
# (main.tex:172-176), which ties the potentials and the four-velocity together.
using Test
using Spacetime

@testset "KerrEquatorialPotentials: (ν, ψ, μ, ω̃), d/dr, circular dt/dτ" begin
    spins = (-0.99, -0.5, 0.0, 0.5, 0.9, 0.998)
    radii(a) = (kerr_isco(a), 1.001 * kerr_isco(a), 8.0, 20.0, 50.0)

    @testset "(i) metric identities vs kerr_metric(a, r, π/2) — Float64 1e-14" begin
        for a in spins, r in radii(a)
            p = kerr_equatorial_potentials(a, r)
            @test keys(p) == (:ν, :ψ, :μ, :ω̃)
            g = kerr_metric(a, r, π / 2)
            @test isapprox(exp(2p.ψ), g[4, 4];                          rtol = 1e-14)
            @test isapprox(p.ω̃,       -g[1, 4] / g[4, 4];               rtol = 1e-14, atol = 1e-300)
            @test isapprox(exp(2p.ν), -g[1, 1] + g[1, 4]^2 / g[4, 4];   rtol = 1e-14)
            @test isapprox(exp(2p.μ), g[2, 2];                          rtol = 1e-14)
            # e^{2ν} is the lapse² = −1/g^{tt} (positive everywhere outside the horizon, even
            # inside the equatorial ergosphere r < 2)
            gi = kerr_metric_inv(a, r, π / 2)
            @test isapprox(exp(2p.ν), -1 / gi[1, 1]; rtol = 1e-13)
        end
        # inside the equatorial ergoregion (r₊ < r < 2) the lapse stays positive: ν is finite
        for a in (0.9, 0.998, -0.998)
            rp = horizons(a)[1]
            r = 0.5 * (rp + 2)
            p = kerr_equatorial_potentials(a, r)
            @test isfinite(p.ν) && exp(2p.ν) > 0
            g = kerr_metric(a, r, π / 2)
            @test isapprox(exp(2p.ν), -g[1, 1] + g[1, 4]^2 / g[4, 4]; rtol = 1e-12)
            @test g[1, 1] > 0      # confirms we really are inside the ergoregion
        end
    end

    @testset "(ii) literal appendix.tex:7-21 (A1a-e) transcription, M=1, χ=a" begin
        for a in spins, r in radii(a)
            χ = a; M = 1.0
            Δ  = r^2 - 2M * r + χ^2 * M^2                              # (A1e)
            ω̃  = 2M^2 * χ / (2M^3 * χ^2 + 2M * r^2 + r * Δ)              # (A1a)
            e2ν = 1 - 2M / r + 2χ * M^2 / r * ω̃                         # (A1b)
            e2ψ = 2M^3 * χ^2 / r + 2M * r + Δ                           # (A1c)
            e2μ = r^2 / Δ                                               # (A1d)
            p = kerr_equatorial_potentials(a, r)
            @test isapprox(p.ω̃,   ω̃;            rtol = 1e-14, atol = 1e-300)
            @test isapprox(p.ν,   log(e2ν) / 2;  rtol = 1e-13, atol = 1e-15)
            @test isapprox(p.ψ,   log(e2ψ) / 2;  rtol = 1e-14)
            @test isapprox(p.μ,   log(e2μ) / 2;  rtol = 1e-13, atol = 1e-15)
            @test Δ == delta_kerr(a, r)                                 # (A1e) IS delta_kerr
        end
    end

    @testset "(iii) a=0 Schwarzschild limits" begin
        for r in (6.0, 10.0, 30.0, 50.0)
            p = kerr_equatorial_potentials(0.0, r)
            @test isapprox(p.ν,  log(1 - 2 / r) / 2;   rtol = 1e-14)
            @test isapprox(p.ψ,  log(r);               rtol = 1e-14)
            @test isapprox(p.μ, -log(1 - 2 / r) / 2;   rtol = 1e-14)
            @test p.ω̃ == 0
            @test isapprox(p.ν, log(schwarzschild_lapse(r)) / 2; rtol = 1e-14)   # f = e^{2ν}: lapse² reuse
            d = kerr_equatorial_potentials_dr(0.0, r)
            @test isapprox(d.dν,  1 / (r * (r - 2));  rtol = 1e-14)
            @test isapprox(d.dψ,  1 / r;              rtol = 1e-14)
            @test isapprox(d.dμ, -1 / (r * (r - 2));  rtol = 1e-14)
            @test d.dω̃ == 0
        end
    end

    @testset "d/dr: analytic derivatives vs 5-point central FD (Float64)" begin
        for a in spins, r in (1.001 * kerr_isco(a), 8.0, 20.0, 50.0)
            d = kerr_equatorial_potentials_dr(a, r)
            @test keys(d) == (:dν, :dψ, :dμ, :dω̃)
            h = 1e-3 * (r - horizons(a)[1])      # the potentials vary on the scale r − r₊
            f(x) = kerr_equatorial_potentials(a, x)
            fm2, fm1, fp1, fp2 = f(r - 2h), f(r - h), f(r + h), f(r + 2h)
            for (k, dk) in ((:ν, :dν), (:ψ, :dψ), (:μ, :dμ), (:ω̃, :dω̃))
                fd = (fm2[k] - 8fm1[k] + 8fp1[k] - fp2[k]) / (12h)
                if k === :ω̃ && a == 0
                    @test d[dk] == 0 && fd == 0
                else
                    # 4th-order stencil: truncation ~ (h/(r−r₊))⁴ ≈ 1e-12 relative
                    @test isapprox(fd, d[dk]; rtol = 2e-9)
                end
            end
        end
    end

    @testset "d/dr: cancellation-free single-fraction forms == logarithmic sums; large-r accuracy" begin
        # (a) the evaluated single-fraction forms equal the logarithmic-sum forms in BigFloat(256)
        setprecision(BigFloat, 256) do
            for a in (big"-0.9", big"0.0", big"0.5", big"0.998"), r in (big"2.5", big"12.0", big"1e3", big"1e6")
                d = kerr_equatorial_potentials_dr(a, r)
                Δ = delta_kerr(a, r); D = 2a^2 + 2r^2 + r * Δ; dΔ = 2(r - 1); dD = 3r^2 + a^2
                @test abs(d.dν - (1 / r + dΔ / Δ - dD / D) / 2) < big"1e-70" * abs(d.dν)
                @test abs(d.dψ - (dD / D - 1 / r) / 2) < big"1e-70" * abs(d.dψ)
                @test abs(d.dμ - (1 / r - (r - 1) / Δ)) < big"1e-70" * abs(d.dμ)
                @test abs(d.dω̃ + (2a / D) * dD / D) ≤ big"1e-70" * max(abs(d.dω̃), big"1e-300")
            end
        end
        # (b) Float64 relative accuracy stays at a few eps out to r = 1e6 (the logarithmic sums lose
        #     ~r·eps: 1e-13 at r = 1e3, 1e-10 at r = 1e6)
        for a in (-0.9, 0.0, 0.5, 0.998), r in (1e2, 1e3, 1e4, 1e6)
            d = kerr_equatorial_potentials_dr(a, r)
            ref = setprecision(() -> kerr_equatorial_potentials_dr(big(a), big(r)), BigFloat, 256)
            @test abs(d.dν - ref.dν) ≤ 8 * eps(Float64) * abs(ref.dν)
            @test abs(d.dψ - ref.dψ) ≤ 8 * eps(Float64) * abs(ref.dψ)
            @test abs(d.dμ - ref.dμ) ≤ 8 * eps(Float64) * abs(ref.dμ)
            a == 0 || @test abs(d.dω̃ - ref.dω̃) ≤ 8 * eps(Float64) * abs(ref.dω̃)
        end
    end

    @testset "d/dr vs BigFloat(256) central FD (h = 1e-25, tol 1e-40 abs)" begin
        setprecision(BigFloat, 256) do
            for a in (big"-0.9", big"0.0", big"0.5", big"0.998"), r in (big"6.0", big"12.0", big"40.0")
                d = kerr_equatorial_potentials_dr(a, r)
                h = big"1e-25"
                fp, fm = kerr_equatorial_potentials(a, r + h), kerr_equatorial_potentials(a, r - h)
                for (k, dk) in ((:ν, :dν), (:ψ, :dψ), (:μ, :dμ), (:ω̃, :dω̃))
                    fd = (fp[k] - fm[k]) / (2h)
                    @test abs(fd - d[dk]) < big"1e-40"
                    @test d[dk] isa BigFloat && precision(d[dk]) == 256
                end
            end
        end
    end

    @testset "(i) metric identities at BigFloat(256) — 1e-60" begin
        setprecision(BigFloat, 256) do
            θ = BigFloat(π) / 2
            for a in (big"-0.99", big"-0.5", big"0.0", big"0.5", big"0.9", big"0.998"),
                r in (kerr_isco(a), big"8.0", big"50.0")
                p = kerr_equatorial_potentials(a, r)
                @test p.ν isa BigFloat && precision(p.ν) == 256
                g = kerr_metric(a, r, θ)
                @test abs(exp(2p.ψ) - g[4, 4]) < big"1e-60" * abs(g[4, 4])
                @test abs(p.ω̃ + g[1, 4] / g[4, 4]) < big"1e-60" * (1 + abs(p.ω̃))
                @test abs(exp(2p.ν) - (-g[1, 1] + g[1, 4]^2 / g[4, 4])) < big"1e-60" * abs(exp(2p.ν))
                @test abs(exp(2p.μ) - g[2, 2]) < big"1e-60" * abs(g[2, 2])
            end
        end
    end

    @testset "circular dt/dτ (main.tex:673) — REUSES kerr_circular_ut_uphi" begin
        for a in spins, r in radii(a)
            ut = kerr_circular_dtdtau(a, r)
            @test ut == kerr_circular_ut_uphi(a, r).ut                      # same object, no re-derivation
            # literal main.tex:673 transcription: a⋆ = r/M (M=1), χ = a
            as, χ = r, a
            @test isapprox(ut, (as^(3 / 2) + χ) / (as^(3 / 4) * sqrt(sqrt(as) * (as - 3) + 2χ)); rtol = 1e-14)
            # main.tex:172-176 decomposition u^t = e^{-ν}/√(1−v²), v = (Ω−ω̃)e^{ψ−ν}, Ω = 1/(a + r^{3/2})
            p = kerr_equatorial_potentials(a, r)
            Ω = 1 / (a + r^(3 / 2))
            v = (Ω - p.ω̃) * exp(p.ψ - p.ν)
            @test isapprox(ut, exp(-p.ν) / sqrt(1 - v^2); rtol = 1e-13)
            # u_φ = e^ψ v/√(1−v²) must equal the geodesic Lz (constants_of_motion, x=+1)
            _, Lz, _ = constants_of_motion(a, r, 0.0, 1.0)
            @test isapprox(exp(p.ψ) * v / sqrt(1 - v^2), Lz; rtol = 1e-12)
        end
        @test kerr_circular_dtdtau(0, 10) isa Float64                        # Int inputs promote
        @test isapprox(kerr_circular_dtdtau(0.0, 10.0), 1 / sqrt(1 - 3 / 10); rtol = 1e-14)
        setprecision(BigFloat, 256) do
            ut = kerr_circular_dtdtau(big"0.5", big"10.0")
            @test ut isa BigFloat && precision(ut) == 256
        end
    end

    @testset "type-genericity + domain" begin
        @test kerr_equatorial_potentials(0, 10).ν isa Float64                 # Int promotion
        @test kerr_equatorial_potentials_dr(0, 10).dν isa Float64
        @test kerr_equatorial_potentials(0.5f0, 10.0f0).ψ isa Float32
        # across the horizon Δ < 0: e^{2μ} < 0 and e^{2ν} < 0 ⇒ log of a negative → DomainError
        @test_throws DomainError kerr_equatorial_potentials(0.5, 1.5)
    end
end

# Geometry/geometry.jl — tortoise_kerr + tortoise_kerr_inverse (DG-P2). Grounded by:
#   (i)  the EXACT a=0 reduction onto the existing tortoise_schw (delegation ⇒ bit equality),
#   (ii) the defining ODE d(r*)/dr = (r²+a²)/Δ via central finite differences (incl. near-horizon),
#   (iii) round-trip identities r*(r*⁻¹(rs)) = rs at 1e-12 for a ∈ {0, 0.5, 0.9, 0.998},
#   (iv) BigFloat(256) round-trip beyond Float64 precision (type-generic {T} inheritance),
#   (v)  characterization values baked ONE-SHOT from the READ-ONLY DG reference
#        /data/haiyangw/claude/DG-Teulolsky/src/Backgrounds.jl `tortoise(r; M=1, a)`
#        (captured 2026-07-02 via a single include+print run — never a live test-dep).
using Test
using Spacetime

@testset "Geometry: Kerr tortoise coordinate + inverse" begin
    @testset "a = 0 reduction: tortoise_kerr(r, 0) == tortoise_schw(r) exactly" begin
        for r in (2.001, 2.5, 3.0, 6.0, 10.0, 50.0, 300.0)
            @test tortoise_kerr(r, 0.0) == tortoise_schw(r)   # exact: a=0 delegates
        end
    end

    @testset "defining ODE d(r*)/dr = (r²+a²)/Δ (central FD, incl. near-horizon)" begin
        # measured worst-case FD-vs-analytic rel. error 2.0e-10 over these points with
        # h = 1e-5·min(r−r₊, 1) (DG reference arithmetic); rtol 1e-7 ⇒ ~500× headroom.
        for (a, dr) in ((0.0, 1.0), (0.5, 8.0), (0.9, 4.5), (0.9, 0.05), (0.998, 0.05), (0.998, 0.44))
            rp = horizons(a)[1]
            r  = rp + dr
            h  = 1e-5 * min(dr, 1.0)
            fd = (tortoise_kerr(r + h, a) - tortoise_kerr(r - h, a)) / (2h)
            @test fd ≈ (r^2 + a^2) / delta_kerr(a, r) rtol = 1e-7
        end
    end

    @testset "round-trip r*(r*⁻¹(rs)) = rs to 1e-12, a ∈ {0, 0.5, 0.9, 0.998}" begin
        # rs-grid spans deep near-horizon (rs = −12 ⇒ r−r₊ ~ 1e-3 at a=0) to far field (1e3).
        # Float64 floor is the forward-map ulp quantization (measured ≤ 5e-13 on this grid).
        for a in (0.0, 0.5, 0.9, 0.998), rs in (-12.0, -6.0, -2.0, 0.0, 3.0, 10.0, 50.0, 200.0, 1000.0)
            r = tortoise_kerr_inverse(rs, a)
            @test r > horizons(a)[1]                            # exterior root only
            @test isapprox(tortoise_kerr(r, a), rs; atol = 1e-12, rtol = 1e-12)
        end
        # r-space direction (inverse∘forward), well-conditioned outside the horizon
        for a in (0.0, 0.5, 0.9, 0.998), dr in (1e-3, 0.1, 1.0, 8.0, 98.0)
            r = horizons(a)[1] + dr
            @test tortoise_kerr_inverse(tortoise_kerr(r, a), a) ≈ r rtol = 1e-12
        end
    end

    @testset "BigFloat(256) round-trip beyond Float64 precision" begin
        setprecision(BigFloat, 256) do
            for a in (big"0.5", big"0.9", big"0.998")
                for rs in (big"-8.0", big"0.0", big"25.0")
                    r = tortoise_kerr_inverse(rs, a)
                    @test r isa BigFloat                        # no Float64 down-cast
                    @test isapprox(tortoise_kerr(r, a), rs; atol = big"1e-40", rtol = big"1e-40")
                end
                rb = big"6.0"
                @test tortoise_kerr_inverse(tortoise_kerr(rb, a), a) ≈ rb rtol = big"1e-40"
            end
        end
    end

    @testset "characterization vs DG-Teukolsky reference (baked one-shot)" begin
        # Provenance: captured 2026-07-02 from the READ-ONLY reference repo by
        #   include("/data/haiyangw/claude/DG-Teulolsky/src/Backgrounds.jl");
        #   Backgrounds.tortoise(r; M=1.0, a=a)   (M = 1, matching Spacetime's convention)
        # printed at %.17g. NOT a live dependency.
        for (r, a, rs_ref) in ((10.0, 0.5,   12.775940471833115),
                               (6.0,  0.9,   7.423934457216939),
                               (4.0,  0.998, 4.143720757992674),
                               (3.0,  0.7,   1.8175272315687372))
            @test tortoise_kerr(r, a) ≈ rs_ref rtol = 1e-14
        end
    end
end

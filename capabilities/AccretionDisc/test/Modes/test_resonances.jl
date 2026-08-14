# Module 2 — Modes/Resonances.jl. Analytic Wronskian (eq:wronskain) + Lindblad/corotation
# resonance structure (RT-A4). GROUNDING:
#   • lindblad_locations(a=0) — CHARACTERIZED vs the PUBLISHED Hirata tab:ilrS literals (the
#     2010-link, verbatim from DiscModes/test/verify_rel_disc.jl): r_p=20,m=2 → ILR 14.50 /
#     OLR 25.47 ; r_p=20,m=1 → 10.14 / 30.65 ; r_p=8,m=2 → 6.97 / 9.55. Published values are
#     2-decimal ⇒ assert at rtol 1e-3.
#   • closed-form identity: lindblad_residual == 0 ⟺ m(Ω−Ω_p)=±√κ² at the returned radii.
#   • wronskian — CANONICAL cold-residual form (F6_QRES collapsed): characterized vs captured
#     DiscModes fixtures (full Float64) + the analytic sign/zero structure (negative in the
#     evanescent [ILR,OLR] band, positive outside).

using Test
using AccretionDisc.Modes

@testset "Modes.Resonances" begin

    @testset "lindblad_locations(a=0) vs PUBLISHED Hirata tab:ilrS" begin
        d0 = RelativisticDisc(a = 0.0)
        # (r_p, m, ILR_pub, OLR_pub) — 2-decimal Hirata 2010 literals
        for (rp, m, ilr, olr) in [(20.0, 2, 14.50, 25.47), (20.0, 1, 10.14, 30.65), (8.0, 2, 6.97, 9.55)]
            rILR, rOLR = lindblad_locations(d0, rp, m)
            @test rILR ≈ ilr rtol = 1e-3
            @test rOLR ≈ olr rtol = 1e-3
        end
    end

    @testset "closed-form identity: residual=0 ⟺ m(Ω−Ω_p)=±√κ² at the locations" begin
        for (a, rp, m) in [(0.0, 20.0, 2), (0.0, 20.0, 1), (0.0, 8.0, 2), (0.6, 15.0, 2)]
            d = RelativisticDisc(a = a)
            rILR, rOLR = lindblad_locations(d, rp, m)
            for rr in (rILR, rOLR)
                isfinite(rr) || continue
                @test abs(lindblad_residual(d, rr, rp, m)) < 1e-12
                Ω = disc_OmegaK(d, rr); Ωp = disc_OmegaK(d, rp); κ = sqrt(kappar2(d, rr))
                @test minimum(abs.((m * (Ω - Ωp) - κ, m * (Ω - Ωp) + κ))) < 1e-12
            end
        end
    end

    @testset "corotation_residual zero at r=r_p" begin
        d = RelativisticDisc(a = 0.6)
        @test abs(corotation_residual(d, 15.0, 15.0)) < 1e-14
        @test corotation_residual(d, 10.0, 15.0) > 0      # Ω decreasing ⇒ inner faster
        @test corotation_residual(d, 30.0, 15.0) < 0
    end

    @testset "wronskian: captured DiscModes fixtures (cold canonical) + sign structure" begin
        d = RelativisticDisc(a = 0.6)
        # captured from packages/DiscModes wronskian(d,r,15,2) (F6_QRES off = cold canonical)
        @test wronskian(d, 10.0, 15.0, 2) ≈ 1.583963758894032 rtol = 1e-13
        @test wronskian(d, 15.0, 15.0, 2) ≈ -0.9710178730797355 rtol = 1e-13
        @test wronskian(d, 30.0, 15.0, 2) ≈ 1.8483961350930593 rtol = 1e-13
        # analytic sign/zero structure: negative inside the [ILR,OLR] evanescent band (resid<0),
        # positive in the propagating regions outside it
        d0 = RelativisticDisc(a = 0.0)
        rILR, rOLR = lindblad_locations(d0, 20.0, 2)
        @test wronskian(d0, (rILR + rOLR) / 2, 20.0, 2) < 0
        @test wronskian(d0, rILR - 1.0, 20.0, 2) > 0
        @test wronskian(d0, rOLR + 1.0, 20.0, 2) > 0
        # the Wronskian numerator IS the Lindblad residual (cold canonical form)
        @test wronskian(d0, rILR, 20.0, 2) ≈ 0.0 atol = 1e-12
    end

    @testset "BigFloat inheritance" begin
        d = RelativisticDisc(a = big"0.0", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
        rILR, rOLR = lindblad_locations(d, big"20.0", 2)
        @test rILR isa BigFloat
        @test abs(rILR - 14.498141476285593) < 1e-6
        @test wronskian(d, big"10.0", big"15.0", 2) isa BigFloat
    end

end

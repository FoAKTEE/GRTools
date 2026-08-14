# Mappings characterization.
#
# Grounding:
#  (A) CAPTURED FIXTURES — legacy FastEMRIWaveforms run ONCE (heavy pkg, never a live dep), outputs
#      baked verbatim:  FEW._PN_E/_PN_L/_PN_C, FEW.xI_to_Y, FEW.kerrecceq_forward_map (flux+amplitude),
#      FEW.ELdot_to_PEdot_Jacobian (Kerr + Schwarzschild). See capture command in the migration log.
#  (B) ANALYTIC identities — pn_carter = Lz²(1/Y²−1) and pn_carter(Y=1)=0 (closed form); kerrecceq
#      forward∘backward = id and backward∘forward = id on each chart; Y_to_xI∘xI_to_Y = id with x=±1⇒Y=±1;
#      p_to_y closed form log(p − 2e − 21/10).

using Test
using EMRIWaveforms.Mappings

@testset "Mappings" begin
    # ---- (A) PN constants of motion: captured FEW literals ----
    @testset "pn constants vs captured FEW" begin
        # (a, p, e, Y) => (E, L, C) from FastEMRIWaveforms._PN_E/_PN_L/_PN_C
        cases = (
            ((0.5, 10.0, 0.3, 0.8), (0.9579979289332909, 2.913037914834036, 4.773256814959102)),
            ((0.9, 12.0, 0.1, 0.6), (0.9610667749796716, 2.3012227829211453, 9.414446749573939)),
            ((0.0, 8.0,  0.2, 1.0), (0.9503819266229829, 3.592106040535496, 0.0)),
        )
        for ((a, p, e, Y), (E, L, C)) in cases
            @test pn_energy(a, p, e, Y)            ≈ E rtol=1e-13
            @test pn_angular_momentum(a, p, e, Y)  ≈ L rtol=1e-13
            if C == 0.0
                @test pn_carter(a, p, e, Y) == 0.0
            else
                @test pn_carter(a, p, e, Y)        ≈ C rtol=1e-13
            end
        end
    end

    # ---- (B) PN Carter identity Q = Lz²(1/Y² − 1) ----
    @testset "pn_carter identity" begin
        for (a, p, e, Y) in ((0.5, 10.0, 0.3, 0.8), (0.9, 12.0, 0.1, 0.6), (-0.4, 9.0, 0.25, 0.7))
            L = pn_angular_momentum(a, p, e, Y)
            @test pn_carter(a, p, e, Y) ≈ L^2 * (1 / Y^2 - 1) rtol=1e-13
        end
        @test pn_carter(0.5, 10.0, 0.3, 1.0) == 0.0     # equatorial closed form
    end

    # ---- (A) inclination: captured FEW literal + analytic inverse ----
    @testset "inclination maps" begin
        @test xI_to_Y(0.5, 10.0, 0.3, 0.9) ≈ 0.899865541707073 rtol=1e-13
        # x = ±1 ⇒ Y = ±1 (equatorial: Q = 0 ⇒ Y = sign(Lz))
        @test xI_to_Y(0.5, 10.0, 0.3, 1.0) ≈ 1.0 atol=1e-12
        # Y_to_xI ∘ xI_to_Y = id (KB Brent is the numeric inverse)
        for (a, p, e, x) in ((0.5, 10.0, 0.3, 0.7), (0.9, 12.0, 0.1, 0.6), (-0.3, 9.0, 0.2, 0.5))
            Y = xI_to_Y(a, p, e, x)
            @test Y_to_xI(a, p, e, Y) ≈ x rtol=1e-8
        end
        @test Y_to_xI(0.5, 10.0, 0.3, 1.0) == 1.0       # |Y|>0.998 ⇒ equatorial passthrough
        @test Y_to_xI(0.5, 10.0, 0.3, -1.0) == -1.0
    end

    # ---- (A) kerrecceq forward: captured FEW literals ----
    @testset "kerrecceq forward vs captured FEW" begin
        u, w, y, z = kerrecceq_forward_map(0.5, 10.0, 0.3, 1.0)
        @test all((u, w, y, z) .≈ (0.818546607784113, 0.35072382359109744, 1.0, 0.5981667063674488))
        ua, wa, ya, za = kerrecceq_forward_map(0.5, 10.0, 0.3, 1.0; kind = "amplitude")
        @test all((ua, wa, ya, za) .≈ (0.8750420873861291, 0.35072382359109744, 1.0, 0.5981667063674488))
    end

    # ---- (B) kerrecceq round-trips (both charts, both kinds) ----
    @testset "kerrecceq round-trips" begin
        for kind in ("flux", "amplitude")
            # region A apex (p close to separatrix)
            for (a, p, e) in ((0.5, 10.0, 0.3), (0.9, 8.0, 0.1), (-0.4, 9.0, 0.2))
                u, w, y, z, near = kerrecceq_forward_map(a, p, e, 1.0; kind = kind, return_mask = true)
                @test near
                a2, p2, e2, x2 = kerrecceq_backward_map(u, w, y, z; regionA = true, kind = kind)
                @test all((a2, p2, e2, x2) .≈ (a, p, e, 1.0))
            end
            # region B apex (p far from separatrix)
            for (a, p, e) in ((0.5, 50.0, 0.3), (0.0, 80.0, 0.1))
                u, w, y, z, near = kerrecceq_forward_map(a, p, e, 1.0; kind = kind, return_mask = true)
                @test !near
                a2, p2, e2, x2 = kerrecceq_backward_map(u, w, y, z; regionA = false, kind = kind)
                @test all((a2, p2, e2, x2) .≈ (a, p, e, 1.0))
            end
        end
        @test_throws ArgumentError kerrecceq_forward_map(0.5, 10.0, 0.3, 0.9)   # xI≠1 unsupported
    end

    # ---- (A) flux Jacobian: captured FEW literals ----
    @testset "ELdot_to_PEdot_Jacobian vs captured FEW" begin
        pd, ed = ELdot_to_PEdot_Jacobian(0.5, 10.0, 0.3, 1.0, -1e-4, -2e-4)
        @test pd ≈ 0.0018819194508862339 rtol=1e-12
        @test ed ≈ -0.0041218974158078675 rtol=1e-12
        pds, eds = ELdot_to_PEdot_Jacobian(0.0, 10.0, 0.3, 1.0, -1e-4, -2e-4)  # Schwarzschild branch
        @test pds ≈ 0.006059754073153514 rtol=1e-12
        @test eds ≈ -0.005068992312549547 rtol=1e-12
    end

    # ---- (B) p_to_y re-export closed form ----
    @testset "p_to_y re-export" begin
        @test p_to_y(10.0, 0.3) ≈ log(10.0 - 2 * 0.3 - 21 / 10) rtol=1e-15
    end
end

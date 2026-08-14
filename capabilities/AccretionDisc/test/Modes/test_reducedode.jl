# Module 3 — Modes/ReducedODE.jl. RT-A1 reduced radial-ODE coefficients 𝔥''+A_m𝔥'+B_m𝔥=𝔉_m.
# GROUNDING:
#   • corotation_exponent(a=0.6,r_p=15,m=2) — vs the PUBLISHED RT-A6 value −0.0818 (the
#     NON-CIRCULAR validation), plus the Frobenius indicial-equation identity s²−s−C=0 with
#     C=4κ²/(9m²r_pΔΩ⁴) (closed form).
#   • A_m — vs an INDEPENDENT central-FD of the closed-form wronskian (analytic −d/dr log𝒲),
#     plus captured DiscModes A_m fixtures (default FD step).
#   • B_m_wkb — k² (=lindblad_residual/c_s²) turning points coincide with Resonances
#     .lindblad_locations (cross-check residual=0 there), plus captured DiscModes fixtures.
#   • B_m_corotation / B_m_full — captured DiscModes fixtures + the B_m_full=B_m_wkb+B_m_corotation
#     identity. [HOLE]: the FULL height-averaged B_m/𝔉_m reduction — B_m_full is an honest
#     physics-motivated construction (WKB dispersion + analytic A_m + corotation gate).

using Test
using AccretionDisc.Modes

@testset "Modes.ReducedODE" begin

    @testset "corotation_exponent vs PUBLISHED RT-A6 + indicial identity" begin
        d = RelativisticDisc(a = 0.6)
        s = corotation_exponent(d, 15.0, 2)
        @test s ≈ -0.0818 atol = 5e-4                       # PUBLISHED RT-A6 non-circular value
        # captured DiscModes fixtures (the other RT-A6 spot-checks)
        @test corotation_exponent(RelativisticDisc(a = 0.0), 20.0, 2) ≈ -0.08001702827280832 rtol = 1e-12
        @test corotation_exponent(RelativisticDisc(a = 0.9), 10.0, 5) ≈ -0.01384867565400838 rtol = 1e-12
        @test corotation_exponent(RelativisticDisc(a = 0.6), 15.0, 1) ≈ -0.27728060353005235 rtol = 1e-12
        # Frobenius indicial identity s²−s−C=0, C=4κ²/(9m²r_pΔΩ⁴)
        for (a, rp, m) in [(0.6, 15.0, 2), (0.0, 20.0, 2), (0.9, 10.0, 5)]
            dd = RelativisticDisc(a = a)
            C = 4 * kappar2(dd, rp) / (9 * m^2 * rp * delta_r(dd, rp) * disc_OmegaK(dd, rp)^4)
            ss = corotation_exponent(dd, rp, m)
            @test abs(ss^2 - ss - C) < 1e-12
        end
    end

    @testset "A_m = −d/dr log𝒲 (independent FD of closed-form wronskian) + fixtures" begin
        d0 = RelativisticDisc(a = 0.0); rp = 20.0; m = 2
        for r in (10.0, 14.0, 30.0)
            h = 1e-6
            Wp = wronskian(d0, r + h, rp, m); Wm = wronskian(d0, r - h, rp, m); W0 = wronskian(d0, r, rp, m)
            ref = -(Wp - Wm) / (2h) / W0                     # independent analytic −dlog𝒲
            @test A_m(d0, r, rp, m; h = h) ≈ ref rtol = 1e-12
        end
        # captured DiscModes A_m fixtures (default FD step)
        @test A_m(d0, 10.0, rp, m) ≈ 0.6307316329801105 rtol = 1e-9
        @test A_m(d0, 14.0, rp, m) ≈ 2.3155706246669197 rtol = 1e-9
        @test A_m(d0, 30.0, rp, m) ≈ -0.19769353241268975 rtol = 1e-9
    end

    @testset "B_m_wkb turning points at lindblad_locations + fixtures" begin
        d0 = RelativisticDisc(a = 0.0); rp = 20.0; m = 2
        rILR, rOLR = lindblad_locations(d0, rp, m)
        # the WKB wavenumber k² = lindblad_residual/c_s² vanishes at the Lindblad resonances
        for rr in (rILR, rOLR)
            @test abs(lindblad_residual(d0, rr, rp, m)) / csound(d0, rr)^2 < 1e-10
        end
        # captured DiscModes B_m_wkb fixtures (cold canonical residual)
        @test B_m_wkb(d0, 10.0, rp, m) ≈ 0.9828028665369648 rtol = 1e-9
        @test B_m_wkb(d0, 14.0, rp, m) ≈ 3.407215965157347 rtol = 1e-9
        @test B_m_wkb(d0, 30.0, rp, m) ≈ 0.23490293674624146 rtol = 1e-9
    end

    @testset "B_m_corotation / B_m_full fixtures + composition identity" begin
        d0 = RelativisticDisc(a = 0.0); rp = 20.0; m = 2
        @test B_m_corotation(d0, 10.0, rp, m) ≈ -0.0008641975308641976 rtol = 1e-12
        @test B_m_corotation(d0, 14.0, rp, m) ≈ -0.0024005486968449933 rtol = 1e-12
        @test B_m_full(d0, 10.0, rp, m) ≈ 0.9819386690061006 rtol = 1e-9
        @test B_m_full(d0, 30.0, rp, m) ≈ 0.23403873921537727 rtol = 1e-9
        for r in (10.0, 14.0, 30.0)
            @test B_m_full(d0, r, rp, m) ≈ B_m_wkb(d0, r, rp, m) + B_m_corotation(d0, r, rp, m) rtol = 1e-14
        end
    end

    @testset "BigFloat inheritance" begin
        d = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
        @test corotation_exponent(d, big"15.0", 2) isa BigFloat
        @test A_m(d, big"12.0", big"15.0", 2) isa BigFloat
    end

    # ── Lin-Shu / Goldreich-Tremaine density-wave dispersion D = κ²−ω̃² (PRIMARY GOAL) ──
    @testset "dispersion D=κ²−ω̃² ≡ −lindblad_residual + turning points" begin
        d = RelativisticDisc(a = 0.6, h0 = 0.1); rp = 15.0; m = 2
        # EXACT algebraic identity: D(d,r,r_p,m) == −lindblad_residual(d,r,r_p,m) (rtol 1e-14).
        for r in (10.0, 12.0, 18.0, 25.0)
            @test dispersion(d, r, rp, m) ≈ -lindblad_residual(d, r, rp, m) rtol = 1e-14
        end
        # D=0 ⇔ the Lindblad turning points (the dispersion vanishes at lindblad_locations).
        rILR, rOLR = lindblad_locations(d, rp, m)
        @test abs(dispersion(d, rILR, rp, m)) < 1e-10
        @test abs(dispersion(d, rOLR, rp, m)) < 1e-10
        # BigFloat inheritance.
        db = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
        @test dispersion(db, big"12.0", big"15.0", 2) isa BigFloat
    end

    # ── reduced_operator: the full assembled operator 𝔥''+A_m𝔥'+B_m𝔥=𝔉_m (composition identity) ──
    @testset "reduced_operator composition + cross-checks" begin
        d = RelativisticDisc(a = 0.6, h0 = 0.1); rp = 15.0; m = 2
        for r in (12.0, 18.0, 25.0)
            op = reduced_operator(d, r, rp, m)
            # principal symbol normalised to 1.
            @test op.principal == 1
            # fields == the already-validated primitives BIT-FOR-BIT.
            @test op.A_m == A_m(d, r, rp, m)
            @test op.B_m == B_m_full(d, r, rp, m)
            @test op.forcing == forcing_Fm(d, r, rp, m)
            @test op.dispersion == dispersion(d, r, rp, m)
            @test op.corotation == corotation_exponent(d, rp, m)
            @test op.lindblad == lindblad_locations(d, rp, m)
            # master-PDE cross-check threaded: sign + order-of-magnitude vs the exact Abel A_m.
            @test sign(real(op.master.A_m)) == sign(A_m(d, r, rp, m))
            let ratio = abs(real(op.master.A_m) / A_m(d, r, rp, m))
                @test max(ratio, 1 / ratio) < 10
            end
        end
        # corotation field == the published RT-A6 −0.0818.
        @test reduced_operator(d, 18.0, rp, m).corotation ≈ -0.0818 atol = 5e-4
        # BigFloat inheritance: every numeric field flows BigFloat.
        db = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
        opb = reduced_operator(db, big"18.0", big"15.0", 2)
        @test opb.A_m isa BigFloat
        @test opb.dispersion isa BigFloat
        @test opb.corotation isa BigFloat
    end

end

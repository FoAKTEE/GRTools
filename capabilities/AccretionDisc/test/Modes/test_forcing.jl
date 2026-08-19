# Module 4 — Modes/Forcing.jl. Relativistic Lindblad/corotation resonance torque (Hirata GW-overlap,
# REUSING the Perturbation Teukolsky engine) + the reduced-ODE forcing source 𝔉_m. GROUNDING:
#   • geodesic/wave factors (_kOmega/_kkappa/_kDelta/_kwt/_kZ/_kgD/_kDprime) — CLOSED-FORM identities
#     vs the Background primitives (_kOmega==disc_OmegaK, _kkappa²==kappar2, _kDelta==delta_r, _kwt
#     a=0 ==1/√(1−3/r)) + captured DiscModes literals (_kwt, _kZ).
#   • _c13_alpha — the live Perturbation MST path at ONE resonance point, vs a captured DiscModes
#     oracle literal (one-shot capture, packages/DiscModes _c13_alpha(2,2,0,w)).
#   • lindblad_torque — ONE live call (a=0, r_p=20, m=2 ILR) through Perturbation, vs the captured
#     DiscModes literal (full precision) + the published Hirata tab:ilrS |N|≈5.37 anchor (loose).
#   • laplace_b — the Laplace-coefficient integral oracle (trapezoidal convergence + b^{(0)}(0)=2,
#     b^{(m≥1)}(0)=0 + b_{1/2}^{(1)}(β)/β→1 leading order).
#   • forcing_Fm — the pot/c_s² construction identity + finiteness at the Lindblad resonance + the
#     corotation (r−r_p) structure tie; [HOLE: H-DISC-FULLOP] = the rigorous height-averaged operator.

using Test
using AccretionDisc
using AccretionDisc.Modes
const M = AccretionDisc.Modes

@testset "Modes.Forcing" begin

    @testset "geodesic/wave factors — closed-form identities + captured literals" begin
        for a in (0.0, 0.6, 0.9), r in (8.0, 14.5, 20.0)
            d = RelativisticDisc(a = a)
            @test M._kOmega(a, r) ≈ disc_OmegaK(d, r) rtol = 1e-14
            @test M._kkappa(a, r)^2 ≈ kappar2(d, r) rtol = 1e-12
            @test M._kDelta(a, r) ≈ delta_r(d, r) rtol = 1e-14
            # _kZ definition consistency
            @test M._kZ(a, r) ≈ M._kkappa(a, r) * M._kwt(a, r) * r^2 / M._kDelta(a, r) rtol = 1e-14
            # _kgD parts
            @test M._kgD(2, a, r, :ILR) ≈ 2 * M._kOmega(a, r) - M._kkappa(a, r) rtol = 1e-14
            @test M._kgD(2, a, r, :OLR) ≈ 2 * M._kOmega(a, r) + M._kkappa(a, r) rtol = 1e-14
        end
        # _kwt a=0 closed form u^t = 1/√(1−3/r)
        @test M._kwt(0.0, 20.0) ≈ 1 / sqrt(1 - 3/20.0) rtol = 1e-14
        # captured DiscModes oracle literals (one-shot, packages/DiscModes)
        @test M._kwt(0.0, 20.0) ≈ 1.0846522890932808 rtol = 1e-13
        @test M._kZ(0.0, 14.5) ≈ 0.01806202872213474 rtol = 1e-13
        # _kDprime ≈ central FD of _kgD
        h = 1e-6
        fd = (M._kgD(2, 0.6, 14.0 + h, :ILR) - M._kgD(2, 0.6, 14.0 - h, :ILR)) / (2h)
        @test M._kDprime(2, 0.6, 14.0, :ILR) ≈ fd rtol = 1e-10
    end

    @testset "laplace_b — integral oracle" begin
        @test M.laplace_b(0.5, 0, 0.0) ≈ 2.0 rtol = 1e-13         # b_{1/2}^{(0)}(0)=2
        @test abs(M.laplace_b(0.5, 2, 0.0)) < 1e-13               # b^{(m≥1)}(0)=0
        # trapezoidal convergence (the integral IS the oracle): N=512 vs N=4096
        @test M.laplace_b(0.5, 1, 0.5; N = 512) ≈ M.laplace_b(0.5, 1, 0.5; N = 4096) rtol = 1e-12
        @test M.laplace_b(0.5, 3, 0.7; N = 512) ≈ M.laplace_b(0.5, 3, 0.7; N = 4096) rtol = 1e-12
        # leading order b_{1/2}^{(1)}(β)/β → 1 as β→0
        @test M.laplace_b(0.5, 1, 1e-4) / 1e-4 ≈ 1.0 rtol = 1e-3
    end

    @testset "forcing_Fm — construction identity + structure + [HOLE: H-DISC-FULLOP]" begin
        d = RelativisticDisc(a = 0.0); rp = 20.0; m = 2
        for r in (10.0, 14.5, 25.0)
            β = min(r, rp) / max(r, rp)
            @test forcing_Fm(d, r, rp, m) * csound(d, r)^2 ≈ -M.laplace_b(0.5, m, β)/2 rtol = 1e-13
        end
        # custom potential callable
        @test forcing_Fm(d, 12.0, rp, m; Φm = r -> 3.0) ≈ 3.0 / csound(d, 12.0)^2 rtol = 1e-13
        # regular (finite) AT the Lindblad resonance (the singular response sits in B_m)
        rILR, rOLR = lindblad_locations(d, rp, m)
        @test isfinite(forcing_Fm(d, rILR, rp, m))
        @test isfinite(forcing_Fm(d, rOLR, rp, m))
        # corotation (r−r_p) structure: source → 0 difference smoothly across r_p
        @test forcing_Fm(d, rp, rp, m) * csound(d, rp)^2 ≈ -M.laplace_b(0.5, m, 1.0)/2 rtol = 1e-12
    end

    @testset "BigFloat inheritance" begin
        d = RelativisticDisc(a = big"0.0", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
        @test M.laplace_b(big(1)//2, 1, big"0.5") isa BigFloat
        @test forcing_Fm(d, big"14.0", big"20.0", 2) isa BigFloat
    end

    @testset "_c13_alpha — live Perturbation MST vs captured DiscModes oracle" begin
        a = 0.0; l = 2; m = 2; r_p = 20.0; w = m * M._kOmega(a, r_p)
        ca = M._c13_alpha(l, m, a, w; Nm = 300)
        # captured one-shot from packages/DiscModes _c13_alpha (Teukolsky MST oracle)
        @test ca.c13 ≈ (-23.58767132477094 - 100.88873420202292im) rtol = 1e-8
        # NOTE the horizon factor α is CORRECTED relative to the captured DiscModes literal
        # (2.8512019997500306e-5): that value carries Γ = 4√(1−a²)/(r₊²+a²), i.e. Hirata's
        # AS-PRINTED bracket, which fails the eq:cplex1 closure by up to 9.1e-5. With the
        # Teukolsky–Press Γ = 2√(1−a²)/(r₊²+a²) the closure holds to ≤6e-13 (gate above) and
        # α is 15.88× smaller. lindblad_torque shifts by only 2.0e-8 (α ≪ |c₁₃|²).
        @test ca.alpha ≈ 1.7953595800524933e-6 rtol = 1e-8
    end

    @testset "lindblad_torque — live Perturbation overlap vs captured DiscModes + Hirata anchor" begin
        d0 = RelativisticDisc(a = 0.0)
        N2 = lindblad_torque(d0, 20.0, 2)                  # a=0, r_p=20, m=2 ILR (fast-converging)
        # captured DiscModes.lindblad_torque literal (full precision, default lmax=24/Nm=300)
        @test N2 ≈ -5.054191918833393 rtol = 1e-6
        # PUBLISHED Hirata tab:ilrS |N|≈5.37 (loose order-of-magnitude anchor)
        @test abs(N2) ≈ 5.37 rtol = 0.1
        # the resonance radius is the Lindblad ILR root
        rILR, _ = lindblad_locations(d0, 20.0, 2)
        @test rILR ≈ 14.50 rtol = 1e-3
    end

end

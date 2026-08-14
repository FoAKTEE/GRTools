# Trajectory characterization.
#
# Grounding:
#  (A) CAPTURED FIXTURES — legacy FastEMRIWaveforms run ONCE (heavy pkg, never a live dep): the 5PN-e10
#      fluxes dEdt8H_5PNe10/dLdt8H_5PNe10, pn5_rhs vector, and a full run_inspiral trajectory (the vendored
#      NumericsCore DOPR853 is a verbatim FEW port ⇒ node-for-node reproduction), all baked verbatim.
#  (B) ANALYTIC — pn5_rhs generic-Y branch: KB Ω_θ sign convention, Ẏ=0 at Y=±1, Peters-Mathews
#      ṗ<0/ė<0 shrinkage; separatrix-refined termination: an inspiral that reaches the wall ends at
#      p ≈ p_sep + buffer (cross-checked vs Spacetime.get_separatrix).

using Test
import EMRIWaveforms.Trajectory
using EMRIWaveforms.Trajectory: dEdt8H_5PNe10, dLdt8H_5PNe10, pn5_rhs,
                                 PN5FluxODE, EMRIInspiral, run_inspiral
using Spacetime: get_separatrix

@testset "Trajectory" begin
    # ---- (A) 5PN flux table vs captured FEW ----
    @testset "5PN fluxes vs captured FEW" begin
        @test dEdt8H_5PNe10(0.5, 10.0, 0.3, 0.8, 10, 10) ≈ -6.591129139798274e-5 rtol=1e-13
        @test dLdt8H_5PNe10(0.5, 10.0, 0.3, 0.8, 10, 10) ≈ -0.0014439794071542892 rtol=1e-13
        @test dEdt8H_5PNe10(0.9, 12.0, 0.1, 1.0, 10, 10) ≈ -2.136136463981765e-5 rtol=1e-13
        @test dLdt8H_5PNe10(0.9, 12.0, 0.1, 1.0, 10, 10) ≈ -0.0008886388978593625 rtol=1e-13
        @test dEdt8H_5PNe10(0.3, 15.0, 0.2, 1.0, 10, 10) ≈ -8.171356841815017e-6 rtol=1e-13
        @test dLdt8H_5PNe10(0.3, 15.0, 0.2, 1.0, 10, 10) ≈ -0.0004400908934649867 rtol=1e-13
    end

    # ---- (A) pn5_rhs (equatorial) vs captured FEW ----
    @testset "pn5_rhs vs captured FEW" begin
        rhs = pn5_rhs(0.5, 10.0, 0.3, 1.0)
        expected = [-0.012680531753428606, -0.0006064344520426035, 0.0,
                    0.027881254842826192, 0.02707964320821458, 0.02004827076804142]
        @test rhs ≈ expected rtol=1e-11
    end

    # ---- (A) pn5_rhs generic-inclination (Y ≠ 1) vs captured legacy KerrBackwards ----
    # Reference = legacy KerrBackwards pn5_deriv assembled ONCE (epsilon=1 => FEW raw-flux convention)
    # at the same inclined points; KB's kerr_frequencies/Y_to_xI derive from the same BHPToolkit forms
    # as Spacetime's, so the assembled 6-vector matches to rtol=1e-11 (parity with the equatorial fixture).
    # Capture command (run once):
    #   julia --project=packages/KerrBackwards -e 'using KerrBackwards;
    #     for (a,p,e,Y) in [(0.5,10.0,0.3,0.8),(0.9,12.0,0.1,0.5),(0.3,15.0,0.2,-0.7)]
    #       println(collect(KerrBackwards.pn5_deriv(1.0,a,p,e,Y))) end'
    @testset "pn5_rhs generic-Y vs captured KerrBackwards" begin
        @test pn5_rhs(0.5, 10.0, 0.3, 0.8) ≈
            [-0.01335412425235242, -0.000621294951185097, -8.563704736552314e-6,
             0.02815484984499537, 0.027330301838549165, 0.019726954422231435] rtol=1e-11
        @test pn5_rhs(0.9, 12.0, 0.1, 0.5) ≈
            [-0.007350731314923963, -9.760066875729245e-5, -1.0069597974678027e-5,
             0.023879473333014746, 0.022956694740743463, 0.018027687475459353] rtol=1e-11
        @test pn5_rhs(0.3, 15.0, 0.2, -0.7) ≈
            [-0.004407662969456472, -9.155154340032214e-5, -7.917488059499781e-7,
             -0.016375807853005523, -0.016552797280482173, 0.012409199342657475] rtol=1e-11
    end

    # ---- (B) analytic checks on the generic-Y branch ----
    @testset "pn5_rhs generic-Y analytic limits" begin
        # (a) KB Ω_θ sign: prograde (xI>0) ⇒ Ω_θ>0; retrograde (xI<0) ⇒ Ω_θ<0.
        @test pn5_rhs(0.5, 10.0, 0.3, 0.8)[5]  > 0
        @test pn5_rhs(0.3, 15.0, 0.2, -0.7)[5] < 0
        # (b) Ẏ = 0 recovered at the equatorial/anti-equatorial fixed points Y = ±1.
        @test pn5_rhs(0.5, 10.0, 0.3, -1.0)[3] == 0.0    # generic branch, dY/dt carries (1-Y^2)
        @test pn5_rhs(0.5, 10.0, 0.3, 1.0)[3]  == 0.0    # equatorial branch
        # (c) Peters-Mathews shrinkage: ṗ<0, ė<0.
        for (a, p, e, Y) in [(0.5, 10.0, 0.3, 0.8), (0.3, 15.0, 0.2, -0.7)]
            r = pn5_rhs(a, p, e, Y)
            @test r[1] < 0 && r[2] < 0
        end
        # Signed-spin equatorial mapping is continuous with the generic-inclination
        # branch and must remain an inspiral for a retrograde (a < 0, Y = 1) orbit.
        req = pn5_rhs(-0.5, 12.0, 0.3, 1.0)
        rnear = pn5_rhs(-0.5, 12.0, 0.3, prevfloat(1.0))
        @test req[1] < 0 && req[2] < 0
        @test isapprox(req[1], rnear[1]; rtol = 0.1)
    end

    # ---- (A) full inspiral trajectory vs captured FEW ----
    # The vendored DOPR853 is the verbatim FEW integrator, so node count and times match to ~1e-9; (p,e)
    # match to ~1e-9 and the year-long accumulated phases (~1.4e5 rad) to ~1e-7 — the residual is the
    # Spacetime.kerr_frequencies vs FEW get_fundamental_frequencies implementation difference in the RHS.
    @testset "run_inspiral vs captured FEW" begin
        insp = EMRIInspiral(PN5FluxODE(0.5))
        tr = run_inspiral(insp, 1.0e6, 10.0, 0.5, [12.0, 0.3, 1.0, 0.0, 0.0, 0.0]; T = 1.0, dt = 10.0)
        @test length(tr.t) == 11
        @test tr.t[1] == 0.0
        @test tr.t[end] ≈ 3.1558149763545595e7 rtol=1e-12
        @test tr.p[1] == 12.0
        @test tr.p[end] ≈ 11.514210324545212  rtol=1e-8
        @test tr.e[end] ≈ 0.28117926286695355 rtol=1e-8
        @test tr.Phi_phi[end] ≈ 140778.0280730805 rtol=1e-7
        @test tr.Phi_r[end]   ≈ 107420.94164073259 rtol=1e-7
        # midpoint node
        @test tr.t[6] ≈ 15550.000000000004 rtol=1e-9
        @test tr.p[6] ≈ 11.999776141807663 rtol=1e-9
        @test tr.e[6] ≈ 0.29999125811273736 rtol=1e-9
        @test tr.Phi_phi[6] ≈ 66.78771727966434 rtol=1e-7
    end

    # ---- fork-collapse: backend selector is a no-op over the shared DOPR853 driver ----
    # AHA: KB Verner and FEW DOPR853 are the SAME adaptive-integration decision on the identical
    # pn5_rhs ODE, so :few and :kb canonicalize on ONE driver and are bit-identical to each other and
    # to the captured trajectory fixture. Default backend is :few.
    @testset "run_inspiral backend selector" begin
        insp = EMRIInspiral(PN5FluxODE(0.5))
        y0 = [12.0, 0.3, 1.0, 0.0, 0.0, 0.0]
        tr_default = run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; T = 1.0, dt = 10.0)
        tr_few = run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; T = 1.0, dt = 10.0, backend = :few)
        tr_kb  = run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; T = 1.0, dt = 10.0, backend = :kb)
        # default == :few, and :few == :kb bit-identical (same canonical driver)
        for f in (:t, :p, :e, :x, :Phi_phi, :Phi_theta, :Phi_r)
            @test getfield(tr_default, f) == getfield(tr_few, f)
            @test getfield(tr_few, f) == getfield(tr_kb, f)
        end
        # both reproduce the captured FEW fixture terminal node
        @test tr_kb.p[end] ≈ 11.514210324545212 rtol = 1e-8
        # unknown backend rejected
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; backend = :nope)
    end

    @testset "run_inspiral validates spin and iteration completion" begin
        insp = EMRIInspiral(PN5FluxODE(0.5))
        y0 = [12.0, 0.3, 1.0, 0.0, 0.0, 0.0]
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.4, y0; T = 0.01)
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; max_iter = 0)
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; max_iter = -1)
        @test_throws ErrorException run_inspiral(insp, 1.0e6, 10.0, 0.5, y0;
                                                  T = 1.0, dt = 10.0, max_iter = 1)
        @test_throws ArgumentError EMRIInspiral(PN5FluxODE(0.5); err = 0.0)
        @test_throws ArgumentError EMRIInspiral(PN5FluxODE(0.5), 0.0, false)
        @test_throws ArgumentError run_inspiral(insp, 0.0, 10.0, 0.5, y0)
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; T = -1.0)
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.5, y0; dt = 0.0)
        @test_throws DimensionMismatch run_inspiral(insp, 1.0e6, 10.0, 0.5, y0[1:5])
        @test_throws ArgumentError run_inspiral(insp, 1.0e6, 10.0, 0.5,
                                                [12.0, 0.3, 1.0, NaN, 0.0, 0.0])
        @test_throws DomainError pn5_rhs(0.0, 12.0, 0.3, 1.0)
        @test_throws ArgumentError run_inspiral(insp, floatmax(Float64), floatmax(Float64),
                                                0.5, y0)
        @test_throws ArgumentError run_inspiral(insp, 1e308, 1e-308, 0.5, y0)
        retro = EMRIInspiral(PN5FluxODE(-0.5))
        @test_throws DomainError run_inspiral(
            retro, 1.0e6, 10.0, -0.5, [6.0, 0.3, 1.0, 0.0, 0.0, 0.0]; T = 0.0)
    end

    # ---- (B) separatrix-refined termination ----
    @testset "separatrix-refined termination" begin
        # Long T so the inspiral reaches the wall; terminal p sits on p_sep + buffer.
        insp = EMRIInspiral(PN5FluxODE(0.5))
        tr = run_inspiral(insp, 1.0e6, 30.0, 0.5, [9.0, 0.2, 1.0, 0.0, 0.0, 0.0]; T = 5.0, dt = 10.0)
        p_end, e_end, x_end = tr.p[end], tr.e[end], tr.x[end]
        p_sep = get_separatrix(0.5, e_end, x_end)
        @test p_end ≈ p_sep + 0.1 atol=1e-6     # buffer = separatrix_buffer_dist = 0.1
    end
end

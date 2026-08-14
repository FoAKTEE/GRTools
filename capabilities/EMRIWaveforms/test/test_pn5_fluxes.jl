# Pn5Fluxes generic-inclination direct fluxes (dp/dt, de/dt, dY/dt) characterization.
#
# Grounding:
#  (A) ONE-SHOT CAPTURED FIXTURE (no heavy live dep). The legacy KerrBackwards pn5fluxes.jl is a
#      self-contained standalone (zero module/using/import/include — pure Base arithmetic), so it was
#      `include`d directly in a scratch Julia session and the three generic-Y direct fluxes evaluated
#      at 4 representative inclined (a,p,e,Y) points at the FEW truncation orders. Baked as literals,
#      asserted rtol=1e-13 (same pattern/tol as the dEdt/dLdt fixture in test_trajectory.jl).
#      Capture command (run once at repo root):
#        julia -e 'include("packages/KerrBackwards/src/pn5fluxes.jl");
#          for (q,p,e,Y) in [(0.5,10.0,0.3,0.8),(0.9,12.0,0.1,0.5),(0.3,15.0,0.2,-0.7),(0.5,10.0,0.3,1.0)]
#            println(dpdt8H_5PNe10(q,p,e,Y,10,10)," ",dedt8H_5PNe10(q,p,e,Y,10,8)," ",dYdt8H_5PNe10(q,p,e,Y,7,10)) end'
#  (B) ANALYTIC LIMITS — dY/dt → 0 at Y = ±1 (equatorial / anti-equatorial fixed points: no
#      inclination change; the flux carries a (1 - Y^2) factor); dp/dt < 0 and de/dt < 0
#      (Peters-Mathews leading-order radiation-reaction: the orbit shrinks and circularizes).

using Test
using EMRIWaveforms.Trajectory: dpdt8H_5PNe10, dedt8H_5PNe10, dYdt8H_5PNe10

@testset "Pn5Fluxes generic-inclination" begin
    # ---- (A) captured legacy KerrBackwards fixtures ----
    @testset "dp/dt de/dt dY/dt vs captured KerrBackwards" begin
        @test dpdt8H_5PNe10(0.5, 10.0, 0.3, 0.8, 10, 10) ≈ -0.01335412425235242    rtol=1e-13
        @test dedt8H_5PNe10(0.5, 10.0, 0.3, 0.8, 10, 8)  ≈ -0.000621294951185097   rtol=1e-13
        @test dYdt8H_5PNe10(0.5, 10.0, 0.3, 0.8, 7, 10)  ≈ -8.563704736552314e-6   rtol=1e-13

        @test dpdt8H_5PNe10(0.9, 12.0, 0.1, 0.5, 10, 10) ≈ -0.007350731314923963   rtol=1e-13
        @test dedt8H_5PNe10(0.9, 12.0, 0.1, 0.5, 10, 8)  ≈ -9.760066875729245e-5   rtol=1e-13
        @test dYdt8H_5PNe10(0.9, 12.0, 0.1, 0.5, 7, 10)  ≈ -1.0069597974678027e-5  rtol=1e-13

        @test dpdt8H_5PNe10(0.3, 15.0, 0.2, -0.7, 10, 10) ≈ -0.004407662969456472  rtol=1e-13
        @test dedt8H_5PNe10(0.3, 15.0, 0.2, -0.7, 10, 8)  ≈ -9.155154340032214e-5  rtol=1e-13
        @test dYdt8H_5PNe10(0.3, 15.0, 0.2, -0.7, 7, 10)  ≈ -7.917488059499781e-7  rtol=1e-13

        @test dpdt8H_5PNe10(0.5, 10.0, 0.3, 1.0, 10, 10) ≈ -0.012775225530141323   rtol=1e-13
        @test dedt8H_5PNe10(0.5, 10.0, 0.3, 1.0, 10, 8)  ≈ -0.0005953619272209951  rtol=1e-13
    end

    # ---- (B) analytic limits ----
    @testset "dY/dt = 0 at Y = ±1 (equatorial fixed points)" begin
        @test dYdt8H_5PNe10(0.5, 10.0, 0.3, 1.0, 7, 10)  == 0.0
        @test dYdt8H_5PNe10(0.5, 10.0, 0.3, -1.0, 7, 10) == 0.0
        @test dYdt8H_5PNe10(0.9, 12.0, 0.1, 1.0, 7, 10)  == 0.0
    end

    @testset "Peters-Mathews shrinkage: dp/dt < 0, de/dt < 0" begin
        for (q, p, e, Y) in [(0.5, 10.0, 0.3, 0.8), (0.9, 12.0, 0.1, 0.5), (0.3, 15.0, 0.2, -0.7)]
            @test dpdt8H_5PNe10(q, p, e, Y, 10, 10) < 0
            @test dedt8H_5PNe10(q, p, e, Y, 10, 8) < 0
        end
    end
end

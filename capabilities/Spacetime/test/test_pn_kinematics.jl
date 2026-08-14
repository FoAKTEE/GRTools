# PNKinematics — PURE ANALYTIC grounding (no oracle: WaSABI is HEAVY, never a live test-dep).
# Every RHS is an independent closed form, mirroring packages/WaSABI/test/pn_kinematics_test.jl.
# Upstream relations baked as expected literals:
#   Ω == √(1/r0³)   (ref-code/WaSABI/InspiralModels/0PASchwarz.m:13)
#   r0 == Ω^(-2/3)  (ref-code/WaSABI/AmplitudeModels/0PASchwarz.m:11)
using Test
using Spacetime

@testset "PNKinematics (analytic)" begin
    # --- round-trip to machine precision ---------------------------------
    for r in (6.0, 7.5, 10.0, 20.0, 30.0, 100.0)
        @test r0_of_omega(kepler_omega(r)) ≈ r rtol = 1e-13
    end

    # --- closed-form spot value (independent of the implementation) ------
    # 6^3 = 216, so Ω(6) = 1/√216; also equals 6^(-3/2).
    @test kepler_omega(6.0) ≈ 1 / sqrt(216.0) rtol = 1e-15
    @test kepler_omega(6.0) ≈ 6.0^(-3 / 2)    rtol = 1e-15

    # --- x = Ω^(2/3) consistency -----------------------------------------
    # pn_x and r0_of_omega are reciprocals: x(Ω)·r0(Ω) = Ω^(2/3)·Ω^(-2/3) = 1.
    for Ω in (0.001, 0.01, 0.05, 0.1)
        @test pn_x(Ω) * r0_of_omega(Ω) ≈ 1.0 rtol = 1e-14
    end
    # On a Kepler orbit, x = 1/r0.
    for r in (6.0, 10.0, 30.0)
        @test pn_x(kepler_omega(r)) ≈ 1 / r rtol = 1e-13
    end

    # --- Unicode aliases agree with the ASCII names ----------------------
    @test kepler_Ω(6.0) == kepler_omega(6.0)
    @test r0_of_Ω(0.05) == r0_of_omega(0.05)

    # --- BigFloat inheritance (type-generic {T}) -------------------------
    rb = big"10.0"
    @test kepler_omega(rb) isa BigFloat
    @test r0_of_omega(kepler_omega(rb)) ≈ rb rtol = big"1e-50"
end

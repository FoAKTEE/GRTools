# Grounding for RingDynamics (Urhoamplitudes + PertGeodesiceqsolution, densitysubstrule,
# LHsnormaftersubst, BLtoNPruleU, ResultingfunctionsRing; Sourceincoordinates.nb).
#
# Oracle: baked `wolframscript` literals of the SAME closed forms transcribed verbatim from the
# notebook (scratchpad/ora_urho.wl, oracle/ring_functions_oracle.wl) — the upstream Mathematica
# source IS the reference; the real Mathematica engine evaluates the published rational functions
# and the Julia port must reproduce them. These cells are pure closed-form algebra (no xAct), so
# the capture is exact to ~16 digits. Literals + tolerances transcribed verbatim from the legacy
# ring_geodesic_test.jl / ring_functions_test.jl. No value taken from new code output.

using ModifiedGravity, Test

@testset "RingDynamics" begin

    # (a) ring_urho_amplitudes on row 5 (m=5) of i1.0_p12.0_e0_r20.0.dat (circular Schwarzschild,
    #     R=20, ringmass=1e-4), read back with read_source_table. wolframscript Urhoamplitudes.
    @testset "ring_urho_amplitudes" begin
        datafile = joinpath(@__DIR__, "data", "i1.0_p12.0_e0_r20.0.dat")
        modes = read_source_table(datafile; radius=20.0, ringmass=1e-4)
        md5 = modes[5]
        @test md5.m == 5
        r = ring_urho_amplitudes(md5)

        @testset "four-velocity vs wolframscript Urhoamplitudes" begin
            @test r.Fourvelocity.U1 ≈ (-0.0000646935949898721 - 0.0031247183360925933im) rtol=1e-12
            @test r.Fourvelocity.U2 ≈ (-0.0010347904424896693 + 0.008068562793766174im)  rtol=1e-12
            @test r.Fourvelocity.U3 ≈ (0.00010911080942188466 + 0.004040975914824931im)  rtol=1e-12
            @test r.Fourvelocity.U4 ≈ (-0.00010911080942188458 - 0.004040975914824931im) rtol=1e-12
        end

        @testset "density response vs wolframscript" begin
            @test r.density.drho ≈ (-1.0942485604835497e-9 - 3.721935153550344e-11im) rtol=1e-11
            @test r.density.dtheta ≈ 0.0 atol=1e-14    # ≈0 (equatorial circular: U4=-U3, cancellation noise)
            @test r.density.dr ≈ (0.09674381701683 + 0.005742016629832808im) rtol=1e-12
        end

        @testset "continuity residual vs wolframscript" begin
            @test r.normalisation.normSource == 0          # Snorm ≡ 0 (Sourcesamplitudes)
            @test r.normalisation.normSourcecheck ≈ (-0.0016215974561984044 - 0.0010231114037686365im) rtol=1e-11
        end

        @testset "_bl_to_np_U fixed-linear-map sanity (Schwarzschild M=1)" begin
            # BLtoNPruleU is a fixed linear map; check the documented combination on a unit input.
            bl = ModifiedGravity._bl_to_np_U(1.0, 20.0, 1.0+0im, 0.0+0im, 0.0+0im, 0.0+0im)
            @test bl.Ut ≈ -(20.0/(2*1.0 - 20.0))
            @test bl.Ur ≈ 1.0
            @test bl.Uθ == 0 && bl.Uφ == 0
        end
    end

    # (b) ring_functions at (t,φ)=(13,0.7) for the inclined i0.5 r=20 (631-mode) and r=400
    #     (211-mode) tables. wolframscript ring_functions_oracle.wl; observables real (Σ + c.c.).
    @testset "ring_functions" begin
        tt, pp = 13.0, 0.7

        @testset "reject heterogeneous ring metadata" begin
            table = read_source_table(joinpath(@__DIR__, "data", "i1.0_p12.0_e0_r20.0.dat");
                                      radius=20.0, ringmass=1e-4)
            a = table[1]
            different_radius = SourceMode(a.m, a.n, a.k, a.ω, 21.0, a.ringmass,
                                          a.Sl, a.Sn, a.Sm, a.Smdag, a.Snorm)
            different_mass = SourceMode(a.m, a.n, a.k, a.ω, a.radius, 2e-4,
                                        a.Sl, a.Sn, a.Sm, a.Smdag, a.Snorm)
            @test_throws ArgumentError ring_functions([a, different_radius])
            @test_throws ArgumentError ring_functions([a, different_mass])
        end

        @testset "r=20 inclined table (631 modes) at (t,φ)=(13,0.7)" begin
            table = read_source_table(joinpath(@__DIR__, "data", "i0.5_p12.0_e0_r20.0.dat");
                                      radius=20.0, ringmass=1e-4)
            @test length(table) == 631
            rf = ring_functions(table)
            @test rf.radius == 20.0 && rf.ringmass == 1e-4

            @test rf.FourvelocityNP.Ul(tt, pp) ≈ (0.27144831492254984 + 0im) rtol=1e-10
            @test rf.FourvelocityNP.Um(tt, pp) ≈ (0.08729626511437064 + 0.19144387094788323im) rtol=1e-10
            @test rf.density.drho(tt, pp)   ≈ (4.6772068443100086e-8 + 0im) rtol=1e-9
            @test rf.density.dtheta(tt, pp) ≈ (1.3425449285001188 + 0im) rtol=1e-10
            @test rf.density.dr(tt, pp)     ≈ (1.7780471154150486 + 0im) rtol=1e-10
            @test rf.Fourvelocity.Ut(tt, pp) ≈ (-0.0783390922521745 + 0im) rtol=1e-10
            @test rf.Fourvelocity.Ur(tt, pp) ≈ (0.6134018128720566 + 0im) rtol=1e-10
            @test rf.Fourvelocity.Uθ(tt, pp) ≈ (0.006172778103463011 + 0im) rtol=1e-10
            @test rf.Fourvelocity.Uφ(tt, pp) ≈ (-0.013537125936385052 + 0im) rtol=1e-10
            @test rf.normalisation(tt, pp)   ≈ (0.021617804319191827 + 0im) rtol=1e-10
            # physical observables are real
            @test abs(imag(rf.density.dr(tt, pp))) < 1e-12
        end

        @testset "r=400 inclined table (211 modes) at (t,φ)=(13,0.7)" begin
            table = read_source_table(joinpath(@__DIR__, "data", "i0.5_p12.0_e0_r400.0.dat");
                                      radius=400.0, ringmass=1e-4)
            @test length(table) == 211
            rf = ring_functions(table)
            @test rf.radius == 400.0

            @test rf.FourvelocityNP.Um(tt, pp) ≈ (-3.002439583373886e-6 + 1.5002042743604788e-5im) rtol=1e-9
            @test rf.density.dr(tt, pp)     ≈ (1.5569274847743824e-5 + 0im) rtol=1e-9
            @test rf.density.dtheta(tt, pp) ≈ (4.027657396298147e-7 + 0im) rtol=1e-9
            @test rf.Fourvelocity.Ut(tt, pp) ≈ (-8.359668881887025e-7 + 0im) rtol=1e-9
            @test rf.Fourvelocity.Uφ(tt, pp) ≈ (-5.304023077826691e-8 + 0im) rtol=1e-9
            @test rf.normalisation(tt, pp)   ≈ (-4.5976247781164425e-7 + 0im) rtol=1e-9
        end
    end
end

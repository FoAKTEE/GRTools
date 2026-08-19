# SpinningSecondary spin-curvature force terms (S5-1) — BAKED characterization vs the
# SpinningOrbit oracle (values captured ONCE by scratchpad/capture_spinforce_s51.jl from the live
# legacy functions at the exact evaluation points of the oracle's own Mathematica-verified test
# refs; see fixtures/spacetime_spinforce_ref.jl provenance header). Tolerances mirror the oracle
# suite (packages/SpinningOrbit/test/runtests.jl): rtol 1e-12 on the closed-form force terms
# (pure algebra), 1e-8/1e-9 on the Darwin-angle maps (trajectory-mediated), 1e-7/1e-6 on
# DFT/δu/solve/assembly (DFT-truncation-limited gates), and rtol 1e-5 on the UrS/UzS pieces
# (finite-difference χ-map derivative — the oracle's loosest comparable gate, its spherical
# UzS at runtests L357). The 8,8 golden reproduces the AUTHORS' published tests.nb Out[]
# literals (vskoupy-KerrSpinningFluxes tests.nb, re-verified in the capture run).
using Test
using Spacetime
const SSF = Spacetime.SpinningSecondary
include(joinpath(@__DIR__, "fixtures", "spacetime_spinforce_ref.jl"))  # SPINFORCE_*

@testset "SpinningSecondary: spin-curvature force terms + correction engine (S5-1)" begin

    @testset "worldline deviations δt3/δr3/δz3/δφ3 (baked, rtol 1e-12)" begin
        for t in SPINFORCE_DEV3
            @test isapprox(SSF._δt3(t.r, t.z, t.Ur, t.Uz, t.K, t.a), t.dt3;   rtol = 1e-12)
            @test isapprox(SSF._δr3(t.r, t.z, t.En, t.Jz, t.K, t.a), t.dr3;   rtol = 1e-12)
            @test isapprox(SSF._δz3(t.r, t.z, t.En, t.Jz, t.K, t.a), t.dz3;   rtol = 1e-12)
            @test isapprox(SSF._δφ3(t.r, t.z, t.Ur, t.Uz, t.K, t.a), t.dphi3; rtol = 1e-12)
        end
    end

    @testset "main force/coefficient layer (28 outputs × 3 tuples, baked, rtol 1e-12)" begin
        for t in SPINFORCE_MAIN
            @test isapprox(SSF._Rt_force(t.r, t.z, t.dr, t.dz, t.K, t.a),   t.Rt;   rtol = 1e-12)
            @test isapprox(SSF._Rphi_force(t.r, t.z, t.dr, t.dz, t.K, t.a), t.Rphi; rtol = 1e-12)
            @test isapprox(SSF._I2_coef(t.r, t.z, t.En, t.Lz, t.K, t.a),    t.I2;   rtol = 1e-12)
            @test isapprox(SSF._I3_coef(t.r, t.z, t.En, t.Lz, t.K, t.a),    t.I3;   rtol = 1e-12)
            @test isapprox(SSF._J_source(t.r, t.z, t.dut, t.dphi, t.En, t.Lz, t.K, t.a), t.J; rtol = 1e-12)
            @test isapprox(SSF._d2r_dλ2(t.r, t.En, t.Lz, t.K, t.a), t.d2rf; rtol = 1e-12)
            @test isapprox(SSF._d2z_dλ2(t.z, t.En, t.Lz, t.Q, t.a), t.d2zf; rtol = 1e-12)
            @test isapprox(SSF._GR_force(t.r, t.z, t.dr, t.a),  t.GR;  rtol = 1e-12)
            @test isapprox(SSF._GZ_force(t.r, t.z, t.dz, t.a),  t.GZ;  rtol = 1e-12)
            @test isapprox(SSF._I1r_coef(t.r, t.z, t.dr, t.d2r, t.Ur, t.a), t.I1r;  rtol = 1e-12)
            @test isapprox(SSF._I1θ_coef(t.r, t.z, t.dz, t.Uth, t.a),       t.I1th; rtol = 1e-12)
            @test isapprox(SSF._KR_force(t.r, t.z, t.dr, t.a),  t.KR;  rtol = 1e-12)
            @test isapprox(SSF._KZ_force(t.r, t.z, t.dz, t.a),  t.KZ;  rtol = 1e-12)
            @test isapprox(SSF._HR_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a), t.HR; rtol = 1e-12)
            @test isapprox(SSF._HZ_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a), t.HZ; rtol = 1e-12)
            @test isapprox(SSF._N1r_coef(t.r, t.z, t.dr, t.Ur, t.a),   t.N1r;  rtol = 1e-12)
            @test isapprox(SSF._N1th_coef(t.r, t.z, t.dz, t.Uth, t.a), t.N1th; rtol = 1e-12)
            @test isapprox(SSF._N2_coef(t.r, t.z, t.En, t.Lz, t.a),    t.N2;   rtol = 1e-12)
            @test isapprox(SSF._N3_coef(t.r, t.z, t.En, t.Lz, t.a),    t.N3;   rtol = 1e-12)
            @test isapprox(SSF._U2_coef(t.r, t.z, t.En, t.Lz, t.a),    t.U2;   rtol = 1e-12)
            @test isapprox(SSF._U3_coef(t.r, t.z, t.En, t.Lz, t.a),    t.U3;   rtol = 1e-12)
            @test isapprox(SSF._QZ_force(t.z),                         t.QZ;   rtol = 1e-12)
            @test isapprox(SSF._MZ_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a), t.MZ; rtol = 1e-12)
            @test isapprox(SSF._MR_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a), t.MR; rtol = 1e-12)
            @test isapprox(SSF._Vgeo(t.r, t.z, t.En, t.Lz, t.K, t.a),  t.Vgeo; rtol = 1e-12)
            @test isapprox(SSF._V_source(t.r, t.z, t.dut, t.dphi, t.En, t.Lz, t.K, t.a), t.V; rtol = 1e-12)
            @test isapprox(SSF._P_source(t.r, t.z, t.dut, t.dphi, t.En, t.Lz, t.a),      t.P; rtol = 1e-12)
        end
    end

    @testset "χ/sin(I)-dependent force layer (11 × 3 tuples, baked, rtol 1e-12)" begin
        for t in SPINFORCE_CHI4
            @test isapprox(SSF._Fr_force(t.chir, t.e, t.p), t.Fr; rtol = 1e-12)
            @test isapprox(SSF._Gr_force(t.r, t.z, t.dr, t.chir, t.a, t.e, t.p), t.Gr; rtol = 1e-12)
            @test isapprox(SSF._Gth_force(t.r, t.z, t.dz, t.chith, t.incl, t.a), t.Gth; rtol = 1e-12)
            @test isapprox(SSF._Kr_force(t.r, t.z, t.dr, t.chir, t.a, t.e, t.p), t.Kr; rtol = 1e-12)
            @test isapprox(SSF._Kth_force(t.r, t.z, t.dz, t.chith, t.incl, t.a), t.Kth; rtol = 1e-12)
            @test isapprox(SSF._Qth_force(t.z, t.chith, t.incl), t.Qth; rtol = 1e-12)
            @test isapprox(SSF._Sr_force(t.r, t.z, t.dr, t.chir, t.a, t.e, t.p), t.Sr; rtol = 1e-12)
            @test isapprox(SSF._SR_force(t.r, t.z, t.dr, t.a), t.SR; rtol = 1e-12)
            @test isapprox(SSF._Sth_force(t.r, t.z, t.dz, t.chith, t.incl, t.a), t.Sth; rtol = 1e-12)
            @test isapprox(SSF._SZ_force(t.r, t.z, t.dz, t.a), t.SZ; rtol = 1e-12)
            @test isapprox(SSF._U1r_coef(t.r, t.z, t.dr, t.Ur, t.a), t.U1r; rtol = 1e-12)
        end
    end

    @testset "heaviest 2nd-order χ forces (9 × 2 tuples, baked, rtol 1e-12)" begin
        for t in SPINFORCE_CHI5
            @test isapprox(SSF._Hr_force(t.r, t.z, t.dr, t.dz, t.d2r, t.chir, t.En, t.Lz, t.a, t.e, t.p), t.Hr; rtol = 1e-12)
            @test isapprox(SSF._Hth_force(t.r, t.z, t.dr, t.dz, t.chith, t.incl, t.En, t.Lz, t.a), t.Hth; rtol = 1e-12)
            @test isapprox(SSF._Mr_force(t.r, t.z, t.dr, t.dz, t.chir, t.En, t.Lz, t.a, t.e, t.p), t.Mr; rtol = 1e-12)
            @test isapprox(SSF._Mth_force(t.r, t.z, t.dr, t.dz, t.chith, t.incl, t.En, t.Lz, t.a), t.Mth; rtol = 1e-12)
            @test isapprox(SSF._Tr_force(t.r, t.z, t.dr, t.dz, t.chir, t.En, t.Lz, t.a, t.e, t.p), t.Tr; rtol = 1e-12)
            @test isapprox(SSF._TR_force(t.r, t.z, t.dr, t.dz, t.En, t.Lz, t.a), t.TR; rtol = 1e-12)
            @test isapprox(SSF._Tth_force(t.r, t.z, t.dr, t.dz, t.d2z, t.chith, t.incl, t.En, t.Lz, t.a), t.Tth; rtol = 1e-12)
            @test isapprox(SSF._TZ_force(t.r, t.z, t.dr, t.dz, t.d2z, t.incl, t.En, t.Lz, t.a), t.TZ; rtol = 1e-12)
            @test isapprox(SSF._U1th_coef(t.r, t.z, t.dz, t.d2z, t.chith, t.En, t.Lz, t.a, t.Uth), t.U1th; rtol = 1e-12)
        end
    end

    @testset "Darwin-angle maps χrhat/χθhat (baked, rtol 1e-8 atol 1e-9)" begin
        for t in SPINFORCE_CHIHAT
            @test isapprox(SSF._chir_hat(t.a, t.p, t.e, t.x, t.lam),  t.chir;  rtol = 1e-8, atol = 1e-9)
            @test isapprox(SSF._chith_hat(t.a, t.p, t.e, t.x, t.lam), t.chith; rtol = 1e-8, atol = 1e-9)
        end
    end

    @testset "2-torus force DFT Im[𝓡t/𝓡φ]four (baked, rtol 1e-7 atol 1e-12)" begin
        d = SPINFORCE_DFT
        Rtf, Rpf = SSF._numeric_force_dft(d.a, d.p, d.e, d.x, d.nmax, d.kmax)
        for n in -d.nmax:d.nmax, k in -d.kmax:d.kmax
            @test isapprox(imag(Rtf[n + d.nmax + 1, k + d.kmax + 1]), d.RtIm[n + d.nmax + 1, k + d.kmax + 1];     rtol = 1e-7, atol = 1e-12)
            @test isapprox(imag(Rpf[n + d.nmax + 1, k + d.kmax + 1]), d.RphiIm[n + d.nmax + 1, k + d.kmax + 1];   rtol = 1e-7, atol = 1e-12)
        end
    end

    @testset "spin-shift potentials δu_{ds,t/φ} (baked, rtol 1e-6 atol 1e-10)" begin
        d = SPINFORCE_DFT
        Rtf, Rpf = SSF._numeric_force_dft(d.a, d.p, d.e, d.x, d.nmax, d.kmax)
        fr = kerr_frequencies(d.a, d.p, d.e, d.x; Time = "Mino")
        for t in SPINFORCE_DELU
            @test isapprox(SSF._delta_u_st(Rtf, fr.Ur, fr.Uth, d.nmax, d.kmax, t.lr, t.lt), t.delut;  rtol = 1e-6, atol = 1e-10)
            @test isapprox(SSF._delta_u_st(Rpf, fr.Ur, fr.Uth, d.nmax, d.kmax, t.lr, t.lt), t.deluph; rtol = 1e-6, atol = 1e-10)
        end
    end

    @testset "FINAL numeric-correction M-matrix solve (baked, rtol 1e-6)" begin
        s = SPINFORCE_SOLVE
        c = SSF._numeric_correction(s.a, s.p, s.e, s.x, s.nmax, s.kmax)
        @test isapprox(c.ΥrS,   s.YrS;     rtol = 1e-6)
        @test isapprox(c.ΥθS,   s.YthS;    rtol = 1e-6)
        @test isapprox(c.udt0S, s.udt0S;   rtol = 1e-6)
        @test isapprox(c.udφ0S, s.udphi0S; rtol = 1e-6)
        @test isapprox(c.error, s.err;     rtol = 1e-6)   # LSQ residual is a reproducible DFT-truncation property
    end

    @testset "orbit-correction helpers h1/h2/h3/Σg** (baked, rtol 1e-12)" begin
        for t in SPINFORCE_OC_HELPERS
            @test isapprox(SSF._oc_h1(t.r, t.z, t.K, t.a),             t.h1;      rtol = 1e-12)
            @test isapprox(SSF._oc_h2(t.r, t.z, t.K, t.En, t.Lz, t.a), t.h2;      rtol = 1e-12)
            @test isapprox(SSF._oc_h3(t.r, t.z, t.K, t.a),             t.h3;      rtol = 1e-12)
            @test isapprox(SSF._oc_Sgtt(t.r, t.z, t.a),                t.Sgtt;    rtol = 1e-12)
            @test isapprox(SSF._oc_Sgphph(t.r, t.z, t.a),              t.Sgphph;  rtol = 1e-12)
            @test isapprox(SSF._oc_Sgtph(t.r, t.a),                    t.Sgtph;   rtol = 1e-12)
            @test isapprox(SSF._oc_dSgttdr(t.r, t.a),                  t.dSgttdr; rtol = 1e-12)
            @test isapprox(SSF._oc_dSgttdz(t.z, t.a),                  t.dSgttdz; rtol = 1e-12)
            @test isapprox(SSF._oc_dSgtphdr(t.r, t.a),                 t.dSgtphdr; rtol = 1e-12)
            @test isapprox(SSF._oc_dSgphphdr(t.r, t.a),                t.dSgphphdr; rtol = 1e-12)
            @test isapprox(SSF._oc_dSgphphdz(t.z, t.a),                t.dSgphphdz; rtol = 1e-12)
        end
    end

    @testset "END-TO-END spin_orbit_correction + spin_orbit_numeric (baked, rtol 1e-6)" begin
        r = SPINFORCE_ORBCORR
        corr = SSF.spin_orbit_correction(r.a, r.p, r.e, r.x, r.nmax, r.kmax)
        @test isapprox(corr.ScriptI, r.ScriptI; rtol = 1e-12)
        @test isapprox(corr.Ehat,  r.Ehat;  rtol = 1e-9)
        @test isapprox(corr.Lzhat, r.Lzhat; rtol = 1e-9)
        @test isapprox(corr.Khat,  r.Khat;  rtol = 1e-9)
        @test isapprox(corr.MinoFrequenciesGeo[1], r.GammaHat; rtol = 1e-8)
        @test isapprox(corr.MinoFrequenciesGeo[2], r.UrHat;    rtol = 1e-8)
        @test isapprox(corr.MinoFrequenciesGeo[3], r.UthHat;   rtol = 1e-8)
        @test isapprox(corr.MinoFrequenciesGeo[4], r.UphHat;   rtol = 1e-8)
        # HEADLINE: linear-in-spin Mino-frequency corrections {ΓS,ΥrS,ΥθS,ΥφS}
        @test isapprox(corr.MinoFrequenciesCorrection[1], r.GammaS; rtol = 1e-6)
        @test isapprox(corr.MinoFrequenciesCorrection[2], r.UrS;    rtol = 1e-6)
        @test isapprox(corr.MinoFrequenciesCorrection[3], r.UthS;   rtol = 1e-6)
        @test isapprox(corr.MinoFrequenciesCorrection[4], r.UphS;   rtol = 1e-6)
        @test isapprox(corr.ES,    r.ES;  rtol = 1e-6)
        @test isapprox(corr.LS,    r.LS;  rtol = 1e-6)
        @test isapprox(corr.error, r.err; rtol = 1e-6)
        for (i, (wr, wz)) in enumerate(zip(r.sampWr, r.sampWz))
            @test isapprox(corr.TrajectoryCorrection[1](wr, wz), r.DtS[i];   rtol = 1e-6, atol = 1e-12)
            @test isapprox(corr.TrajectoryCorrection[2](wr, wz), r.rS[i];    rtol = 1e-6, atol = 1e-12)
            @test isapprox(corr.TrajectoryCorrection[3](wr, wz), r.zS[i];    rtol = 1e-6, atol = 1e-12)
            @test isapprox(corr.TrajectoryCorrection[4](wr, wz), r.DphiS[i]; rtol = 1e-6, atol = 1e-12)
            oc = corr.OrbitCorrection(wr, wz)
            @test isapprox(oc.ΔtS, r.ocDtS[i];   rtol = 1e-6, atol = 1e-12)
            @test isapprox(oc.rS,  r.ocrS[i];    rtol = 1e-6, atol = 1e-12)
            @test isapprox(oc.zS,  r.oczS[i];    rtol = 1e-6, atol = 1e-12)
            @test isapprox(oc.ΔφS, r.ocDphiS[i]; rtol = 1e-6, atol = 1e-12)
            # UrS/UzS go through the finite-difference χ-map derivative (h = 1e-6) — the FD noise
            # floor across independent trajectory implementations is ~1e-6 absolute, so the gate is
            # the oracle's loosest comparable one (its spherical UzS gate, runtests L357: rtol 1e-5).
            @test isapprox(oc.UrS, r.ocUrS[i]; rtol = 1e-5, atol = 1e-8)
            @test isapprox(oc.UzS, r.ocUzS[i]; rtol = 1e-5, atol = 1e-8)
        end
        num = SSF.spin_orbit_numeric(corr, r.spar)
        @test isapprox(num.MinoFrequencies[1], r.numGamma; rtol = 1e-6)
        @test isapprox(num.MinoFrequencies[2], r.numUr;    rtol = 1e-6)
        @test isapprox(num.MinoFrequencies[3], r.numUth;   rtol = 1e-6)
        @test isapprox(num.MinoFrequencies[4], r.numUph;   rtol = 1e-6)
        @test isapprox(num.BLFrequencies[1],   r.numBLr;   rtol = 1e-6)
        @test isapprox(num.BLFrequencies[2],   r.numBLth;  rtol = 1e-6)
        @test isapprox(num.BLFrequencies[3],   r.numBLph;  rtol = 1e-6)
        @test isapprox(num.ES, r.numES; rtol = 1e-6)
        @test isapprox(num.LS, r.numLS; rtol = 1e-6)
        @test isapprox(num.Ehat,  r.numEhat;  rtol = 1e-9)
        @test isapprox(num.Lzhat, r.numLzhat; rtol = 1e-9)
        @test isapprox(num.Khat,  r.numKhat;  rtol = 1e-9)
        @test num.s == r.spar
    end

    @testset "PUBLISHED tests.nb golden: spin_orbit_correction(0.95,10,0.2,0.8,8,8)" begin
        # The authors' own vskoupy-KerrSpinningFluxes tests.nb Out[] literals (also the oracle's
        # runtests L444-459 golden); the capture run re-verified the oracle reproduces them.
        g = SPINFORCE_GOLDEN88
        ocg = SSF.spin_orbit_correction(g.a, g.p, g.e, g.x, g.nmax, g.kmax)
        @test isapprox(ocg.Ehat, 0.9544746770348815028245737919294451657838972804860443023072; rtol = 1e-12)
        @test isapprox(ocg.MinoFrequenciesGeo[1], 124.7629736005981287023759422670878097995901511945922024627606; rtol = 1e-12)
        @test isapprox(ocg.MinoFrequenciesCorrection[1], -0.30438197634945600312974191602274806466; rtol = 1e-10)
        @test isapprox(ocg.BLFrequenciesCorrection[1], 0.00128402591293930362965445382249630736; rtol = 1e-10)
        # and the same values as captured from the Julia oracle (round-trip consistency)
        @test isapprox(ocg.Ehat, g.Ehat; rtol = 1e-12)
        @test isapprox(ocg.MinoFrequenciesGeo[1], g.GammaHat; rtol = 1e-12)
        @test isapprox(ocg.MinoFrequenciesCorrection[1], g.GammaS; rtol = 1e-10)
        @test isapprox(ocg.BLFrequenciesCorrection[1], g.BLc1; rtol = 1e-10)
    end

    @testset "shifted_geodesic_correction: libration zeroed, secular part untouched" begin
        # The SG member of arXiv:2603.12189: rS = zS = UrS = UzS = 0 with ΔtS/ΔφS and every other
        # field (frequencies, constants, ES/LS) preserved. Promoted out of the figure scripts.
        oc = SSF.spin_orbit_correction(0.9, 12.0, 0.2, sqrt(3) / 2, 4, 4)
        sg = SSF.shifted_geodesic_correction(oc)
        for wr in (0.0, 0.7, 2.4, 5.9), wz in (0.0, 1.3, 4.1)
            c = oc.OrbitCorrection(wr, wz)
            g = sg.OrbitCorrection(wr, wz)
            @test g.rS == 0.0 && g.zS == 0.0 && g.UrS == 0.0 && g.UzS == 0.0
            @test g.ΔtS === c.ΔtS
            @test g.ΔφS === c.ΔφS
        end
        # the exact correction really does librate at this orbit, so the wrapper is not a no-op
        @test maximum(abs(oc.OrbitCorrection(wr, wz).rS) for wr in 0:0.5:6, wz in 0:0.5:6) > 1e-6
        # every non-closure field is carried through untouched
        for k in keys(oc)
            k === :OrbitCorrection && continue
            @test getfield(sg, k) === getfield(oc, k)
        end
        # and the SG object is still a valid input to the downstream frequency assembly
        n0 = SSF.spin_orbit_numeric(oc, 1e-3)
        ns = SSF.spin_orbit_numeric(sg, 1e-3)
        @test all(ns.BLFrequencies .=== n0.BLFrequencies)
    end

    @testset "guards: kmax even, nmax > 0" begin
        @test_throws ArgumentError SSF.spin_orbit_correction(0.9, 12.0, 0.2, 0.5, 2, 3)
        @test_throws ArgumentError SSF.spin_orbit_correction(0.9, 12.0, 0.2, 0.5, 0, 2)
    end

    @testset "BigFloat inheritance of the algebraic force terms (type-generic {T})" begin
        # the force/deviation/helper terms are pure arithmetic — BigFloat in ⇒ BigFloat out,
        # agreeing with the Float64 value; (the DFT/pinv engine is Float64-tier by design).
        t = SPINFORCE_MAIN[1]
        vb = SSF._Rt_force(big(t.r), big(t.z), big(t.dr), big(t.dz), big(t.K), big(t.a))
        @test vb isa BigFloat
        @test isapprox(Float64(vb), t.Rt; rtol = 1e-12)
        db = SSF._δt3(big"10.2", big"0.3", big"0.5", big"-0.4", big"13.0", big"0.9")
        @test db isa BigFloat
        @test isapprox(Float64(db), SPINFORCE_DEV3[1].dt3; rtol = 1e-12)
        hb = SSF._oc_h2(big(t.r), big(t.z), big(t.K), big(t.En), big(t.Lz), big(t.a))
        @test hb isa BigFloat
        mb = SSF._MR_force(big(t.r), big(t.z), big(t.dr), big(t.dz), big(t.En), big(t.Lz), big(t.a))
        @test mb isa BigFloat
        @test isapprox(Float64(mb), t.MR; rtol = 1e-12)
    end
end

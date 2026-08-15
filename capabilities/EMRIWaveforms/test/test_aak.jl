# AAK characterization.
#
# Grounding:
#  (A) ak_bessel_coeffs (canonical home HERE since the §5.5 Kernels-purity hoist; Kernels only takes
#      the injected provider) vs SpecialFunctions.besselj closed form, incl. J_{−1}=−J_1 at n=1,
#      plus the captured legacy KerrBackwards literal (tests moved here from Kernels with the fn).
#  (B) AK parameter-map OmegaPhi → Newtonian Kepler limit v³/M as v→0 (closed form).
#  (C) CAPTURED FIXTURES — legacy KerrBackwards run ONCE (heavy pkg, never a live dep): the AK rate
#      polynomials drdm/dtdm/dthetadm/dphidm, OmegaPhi, _ak_mapped_freqs, and aak_summation output samples,
#      baked verbatim.
#  (D) aak_waveform is a documented [HOLE] (needs Trajectory.run_inspiral) → errors.

using Test
using SpecialFunctions: besselj
import EMRIWaveforms.AAK
using EMRIWaveforms.AAK: ak_bessel_coeffs
using EMRIWaveforms.Trajectory: EMRIInspiral, PN5FluxODE, run_inspiral
using Spacetime: kerr_frequencies
using EMRIWaveforms.Mappings: xI_to_Y
using NumericsCore: MTSUN_SI
const PM = AAK.ParameterMapAAK

@testset "AAK" begin
    # ---- (A) ak_bessel_coeffs (canonical home; moved from Kernels with its tests) ----
    @testset "ak_bessel_coeffs" begin
        e = 0.3
        @test ak_bessel_coeffs(1, e)[1] == -besselj(1, 1 * e)     # J_{-1} = -J_1
        J = ak_bessel_coeffs(3, e)
        @test all(J .≈ (besselj(1, 3e), besselj(2, 3e), besselj(3, 3e), besselj(4, 3e), besselj(5, 3e)))
        # moved-from-Kernels grounding (pre-hoist test body, verbatim): exact besselj identity at
        # (n,e) representatives + the J_{-1}(x) = -J_1(x) identity + captured legacy literal.
        for (n, ee) in ((1, 0.3), (2, 0.5), (3, 0.7), (5, 0.2))
            ne = n * ee
            Jn = ak_bessel_coeffs(n, ee)
            expect = n == 1 ?
                (-besselj(1, ne), besselj(0, ne), besselj(1, ne), besselj(2, ne), besselj(3, ne)) :
                (besselj(n - 2, ne), besselj(n - 1, ne), besselj(n, ne), besselj(n + 1, ne), besselj(n + 2, ne))
            @test all(Jn[k] == expect[k] for k in 1:5)
        end
        J1 = ak_bessel_coeffs(1, 0.3)
        @test J1[1] == -J1[3]
        # captured legacy literal (KerrBackwards.ak_bessel_coeffs(1,0.3)):
        @test J1[1] ≈ -0.148318816273104 rtol = 1e-13
        # harmonic guard (n ≥ 1):
        @test_throws ArgumentError ak_bessel_coeffs(0, 0.3)
    end

    # ---- (B) OmegaPhi Newtonian Kepler limit ----
    @testset "OmegaPhi Kepler limit v³/M" begin
        M = 1.0e6
        for v in (1e-3, 1e-4, 1e-5)
            r = PM.OmegaPhi(v, 0.0, 1.0, 0.0, M) * M / v^3
            @test r ≈ 1.0 atol = 50 * v^2     # PN correction is O(v²)
        end
    end

    # ---- (C) captured legacy KerrBackwards parametermap literals ----
    @testset "AK rate polynomials vs captured KB" begin
        @test PM.drdm(0.2, 0.3, 0.8, 0.5)     ≈ 4.711844989819999 rtol=1e-13
        @test PM.dtdm(0.2, 0.3, 0.8, 0.5)     ≈ 759.881870590556  rtol=1e-13
        @test PM.dthetadm(0.2, 0.3, 0.8, 0.5) ≈ 5.287870681820001 rtol=1e-13
        @test PM.dphidm(0.2, 0.3, 0.8, 0.5)   ≈ 5.327964281820001 rtol=1e-13
        @test PM.OmegaPhi(0.2, 0.3, 0.8, 0.5, 1.0e6) ≈ 7.011569150451341e-9 rtol=1e-13
        Or, Ot, Op = PM._ak_mapped_freqs(0.2, 1.0e6, 0.5, 0.3, acos(0.8), 1.0e6)
        @test Or ≈ 0.000986881580198545  rtol=1e-13
        @test Ot ≈ 0.0011075284067355187 rtol=1e-13
        @test Op ≈ 0.0011159258891249722 rtol=1e-13
    end

    # ---- (C) captured legacy KerrBackwards.aak_summation samples ----
    @testset "aak_summation vs captured KB" begin
        t = [0.0, 100.0, 200.0, 300.0, 400.0]
        e = [0.3, 0.29, 0.28, 0.27, 0.26]
        Phivec = [0.0, 1.0, 2.0, 3.0, 4.0]
        gimvec = [0.0, 0.5, 1.0, 1.5, 2.0]
        alpvec = [0.0, 0.2, 0.4, 0.6, 0.8]
        nuvec  = [2.0e-3, 2.01e-3, 2.02e-3, 2.03e-3, 2.04e-3]
        gimdotvec = fill(1.0e-5, 5)
        OmegaPhi = [1.3e-2, 1.31e-2, 1.32e-2, 1.33e-2, 1.34e-2]
        lam = fill(0.6, 5)
        expected = ComplexF64[
            -1.3451585570718374e-22 + 2.214593342976386e-22im,
             2.336833580096155e-22  - 3.7101972843177205e-24im,
            -1.322283799822049e-22  - 2.1129462320733045e-23im,
             1.3822786266437883e-22 - 4.2884571354236225e-23im,
            -1.1103177453632862e-22 + 6.884537892123118e-23im,
        ]
        # fork-collapse: backend=:few|:kb are no-op selectors over the SAME Kernels kernel; default=:few.
        hdef = AAK.aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                                 1.0e6, 0.5, 10.0, 0.7, 0.8, 0.9, 1.0, 1.5, 4; dt = 100.0, out_len = 5)
        # PRE-HOIST BEHAVIOR-IDENTITY GATE (§5.5 Kernels purity, S4-c): these literals are the EXACT
        # (repr round-trip) output of the pre-hoist committed GrTools pipeline (Kernels @ 51efa9f,
        # ak_bessel_coeffs internal) on these inputs, baked 2026-07-02 by
        # scratchpad/s4c_capture_fixtures.jl BEFORE the provider-injection refactor. Public behavior
        # must be IDENTICAL post-hoist (far tighter than the KB rtol=1e-10 gate below).
        pre_hoist = ComplexF64[
            -1.3451585570718374e-22 + 2.214593342976386e-22im,
            2.336833580096155e-22 - 3.7101972843177205e-24im,
            -1.322283799822049e-22 - 2.1129462320733045e-23im,
            1.3822786266437883e-22 - 4.2884571354236225e-23im,
            -1.1103177453632862e-22 + 6.884537892123118e-23im,
        ]
        for i in 1:5
            @test hdef[i] ≈ pre_hoist[i] rtol = 1e-13
        end
        for backend in (:few, :kb)
            h = AAK.aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec, OmegaPhi, lam,
                                  1.0e6, 0.5, 10.0, 0.7, 0.8, 0.9, 1.0, 1.5, 4;
                                  backend = backend, dt = 100.0, out_len = 5)
            @test length(h) == 5
            @test h == hdef                 # backend is a no-op; default == :few
            for i in 1:5
                @test h[i] ≈ expected[i] rtol=1e-10
            end
        end
        @test_throws ArgumentError AAK.aak_summation(t, e, Phivec, gimvec, alpvec, nuvec, gimdotvec,
                                                     OmegaPhi, lam, 1.0e6, 0.5, 10.0, 0.7, 0.8, 0.9,
                                                     1.0, 1.5, 4; backend = :nope, dt = 100.0, out_len = 5)
    end

    # ---- (D) aak_waveform: self-consistency + captured KB fixture ----
    # Grounding:
    #  · SELF-CONSISTENCY (bit-exact): aak_waveform == its own components composed
    #    (run_inspiral → AK basis map → Kernels.aak_summation).
    #  · CAPTURED FIXTURE: one-shot legacy KerrBackwards aak_waveform (Y0=1, T=0.02) early samples.
    #    GrTools run_inspiral differs from the legacy KB integrator at the documented ~1e-7 phase level,
    #    so the early samples (tiny accumulated phase) are matched at the trajectory-limited tolerance.
    @testset "aak_waveform self-consistency" begin
        M, mu, a, p0, e0, Y0 = 1.0e6, 10.0, 0.5, 12.0, 0.3, 1.0
        qS, phiS, qK, phiK, dist = 0.7, 0.8, 0.9, 1.0, 1.5
        T, dt = 0.02, 10.0
        @test_throws ArgumentError AAK.aak_waveform(
            M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist; T = 0.0, dt = dt)
        h = AAK.aak_waveform(M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist; T = T, dt = dt)

        # rebuild the components explicitly (must be bit-identical to the convenience driver).
        insp = EMRIInspiral(PN5FluxODE(a); err = 1e-10)
        tr = run_inspiral(insp, M, mu, a, [p0, e0, Y0, 0.0, 0.0, 0.0]; T = T, dt = dt)
        nn = length(tr.t); Msec = M * MTSUN_SI
        OmP = similar(tr.t); OmR = similar(tr.t); OmT = similar(tr.t); Yv = similar(tr.t)
        for i in 1:nn
            f = kerr_frequencies(a, tr.p[i], tr.e[i], tr.x[i])
            OmP[i] = f.Omegaph / Msec
            OmT[i] = copysign(f.Omegath, tr.x[i]) / Msec
            OmR[i] = f.Omegar / Msec
            Yv[i] = xI_to_Y(a, tr.p[i], tr.e[i], tr.x[i])
        end
        lam = acos.(Yv)
        AAK._aak_lam_fixup!(lam)
        href = AAK.aak_summation(tr.t, tr.e, tr.Phi_r, tr.Phi_theta .- tr.Phi_r,
                                    tr.Phi_phi .- tr.Phi_theta, OmR ./ (2π), OmT .- OmR, OmP, lam,
                                    M, a, mu, qS, phiS, qK, phiK, dist, max(4, floor(Int, 30 * e0));
                                    mich = false, dt = dt)
        @test h == href     # bit-identical composition

        # User-supplied initial phases seed the trajectory state (in the integrator's
        # mass-ratio-scaled convention) and therefore change the waveform.
        hph = AAK.aak_waveform(M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist;
                               T = T, dt = dt, Phi_phi0 = 0.7,
                               Phi_theta0 = 0.2, Phi_r0 = 1.1)
        @test hph != h

        # captured legacy KB aak_waveform early samples (one-shot, baked literals)
        kb = ComplexF64[
            -6.85543261727384609e-23 + 1.38932374246283525e-22im,
            -8.60639156136692719e-23 + 1.28624358066796753e-22im,
            -1.01916237638620770e-22 + 1.15988264666315179e-22im,
            -1.15776562582832936e-22 + 1.01321515606991222e-22im,
            -1.27369578831885728e-22 + 8.49730988432968663e-23im,
            -1.36487692539926568e-22 + 6.73305739724157221e-23im,
        ]
        for i in 1:length(kb)
            @test h[i] ≈ kb[i] rtol = 1e-4    # trajectory-implementation-limited cross-check
        end
    end
end

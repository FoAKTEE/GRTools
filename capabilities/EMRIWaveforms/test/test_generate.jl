# Generate characterization — SSB source/detector frame algebra + the top-level generator dispatcher.
#
# Grounding (NO heavy live dep):
#  (A) ANALYTIC closed forms — get_viewing_angles / ssb_polarization_angle / to_SSB_frame are pure
#      transcriptions of KerrBackwards/waveform.jl; expected values computed INDEPENDENTLY in python
#      and baked as literals, PLUS structural invariants:
#        · to_SSB_frame is a rotation by 2ψ ⇒ h₊²+h_×² is invariant (random angles/strains);
#        · aligned spin qS=qK, phiS=phiK ⇒ R·S=1 ⇒ θ=arccos(−1)=π, φ=−π/2;
#        · ssb_polarization_angle on the dw=0 axis ⇒ ψ=π/2.
#  (B) generate_emri_waveform dispatch — chains the already-grounded submodules; the x0→Y₀ map is
#      Mappings.xI_to_Y (grounded). The real Pn5AAK runtime is a [HOLE] (AAK.aak_waveform needs the
#      spine inspiral driver), so the generator backend is an injectable seam: a STUB generator grounds
#      the dispatch/frame composition/return-shape; the default backend errors as the documented [HOLE].
#      The full few-Python end-to-end array match remains [HOLE].

using Test
using EMRIWaveforms.Generate
using EMRIWaveforms.Generate: get_viewing_angles, ssb_polarization_angle, to_SSB_frame,
                              generate_emri_waveform, fast_schwarzschild_eccentric_flux
using EMRIWaveforms.Mappings: xI_to_Y
using EMRIWaveforms.AmplitudeModels: RomanAmplitude
import EMRIWaveforms.AAK
using NumericsCore: Gpc_in_m, MRSUN_SI
using Random
import Data
using HDF5
using TestKit: oracle_root, oracle_or_skip
# Local FEW data corpus (.fewdata): override GRTOOLS_FEWDATA_DIR at sites with their
# own mirror; when absent, the Data artifact/override chain stays authoritative.
const _FEWDATA_GEN = oracle_root("GRTOOLS_FEWDATA_DIR",
    "/data/haiyangw/claude/emri-disk/ref-code/FastEMRIWaveforms/.fewdata")
if oracle_or_skip(isdir(_FEWDATA_GEN), "FEW data root for the datadir-backed Generate gates";
                  env = "GRTOOLS_FEWDATA_DIR")
    Data.set_datadir!(_FEWDATA_GEN)
end

@testset "Generate" begin

    # =========================================================================
    # (A) Frame algebra
    # =========================================================================
    @testset "get_viewing_angles" begin
        θ, φ = get_viewing_angles(1.2, 0.7, 0.9, 2.1)
        @test θ ≈ 1.9276594629913906 rtol = 1e-13      # python oracle
        @test φ == -π / 2
        # aligned-spin limit: R·S = 1 ⇒ θ = arccos(−1) = π, φ = −π/2.
        # (acos is sqrt-conditioned at the endpoint, so θ = π − O(1e-8); atol, not rtol.)
        θa, φa = get_viewing_angles(0.8, 1.3, 0.8, 1.3)
        @test θa ≈ π atol = 1e-7
        @test φa == -π / 2
        # Roundoff can make a valid self-dot exceed one by a few ulps.
        q, ph = 1.0295521883268661, 1.6555512199244948
        @test get_viewing_angles(q, ph, q, ph)[1] == Float64(π)
    end

    @testset "ssb_polarization_angle" begin
        @test ssb_polarization_angle(1.2, 0.7, 0.9, 2.1) ≈ 2.538936689282923 rtol = 1e-13  # python
        # dw = sin(qK) sin(phiS−phiK) = 0 on the aligned-azimuth axis ⇒ ψ = π/2
        @test ssb_polarization_angle(1.2, 0.7, 0.9, 0.7) == 0.5 * π
    end

    @testset "to_SSB_frame rotation" begin
        # python literal at (hp,hc)=(0.3,−0.4), angles (1.2,0.7,0.9,2.1)
        hp, hc = to_SSB_frame(0.3, -0.4, 1.2, 0.7, 0.9, 2.1)
        @test hp ≈ -0.2663597762663558 rtol = 1e-12
        @test hc ≈ -0.42314591997009365 rtol = 1e-12

        # rotation invariant h₊²+h_×² for random strains/angles
        rng = MersenneTwister(20260630)
        for _ in 1:200
            a, b = randn(rng), randn(rng)
            qS, phiS = π * rand(rng), 2π * rand(rng)
            qK, phiK = π * rand(rng), 2π * rand(rng)
            ap, bp = to_SSB_frame(a, b, qS, phiS, qK, phiK)
            @test ap^2 + bp^2 ≈ a^2 + b^2 rtol = 1e-12
        end

        # elementwise over arrays (vectorized rotation)
        hpv, hcv = to_SSB_frame([0.3, 0.1], [-0.4, 0.2], 1.2, 0.7, 0.9, 2.1)
        @test hpv[1] ≈ -0.2663597762663558 rtol = 1e-12
        @test length(hpv) == 2 && length(hcv) == 2
    end

    # =========================================================================
    # (B) generate_emri_waveform dispatch (injectable AAK backend)
    # =========================================================================
    @testset "generate_emri_waveform dispatch" begin
        # input validation
        @test_throws ArgumentError generate_emri_waveform(:NotAClass, 1e6, 10.0, 0.5, 12.0, 0.3,
                                                          1.0, 1.0, 1.2, 0.7, 0.9, 2.1)
        @test_throws ArgumentError generate_emri_waveform(:Pn5AAKWaveform, 1e6, 10.0, 0.5, 12.0, 0.3,
                                                          1.0, 1.0, 1.2, 0.7, 0.9, 2.1; frame = :nope)

        # STUB backend: a known complex strain, so the dispatch + frame composition is exercised.
        M, mu, a, p0, e0, x0, dist = 1e6, 10.0, 0.5, 12.0, 0.3, 0.95, 1.0
        qS, phiS, qK, phiK = 1.2, 0.7, 0.9, 2.1
        dt = 10.0
        hbase = ComplexF64[0.3 - 0.4im, 0.1 + 0.2im, -0.25 + 0.05im]
        captured = Ref{Any}(nothing)
        stub = function (M_, mu_, a_, p0_, e0_, Y0_, qS_, phiS_, qK_, phiK_, dist_; kwargs...)
            captured[] = (; M_, mu_, a_, p0_, e0_, Y0_, dist_)
            return copy(hbase)
        end

        # detector frame (default): x0 → Y0 = xI_to_Y, distance scaling, then SSB rotation.
        res = generate_emri_waveform(:Pn5AAKWaveform, M, mu, a, p0, e0, x0, dist,
                                     qS, phiS, qK, phiK; dt = dt, aak_generator = stub)
        # x0 → Y0 identity passed to the generator
        @test captured[].Y0_ ≈ xI_to_Y(a, p0, e0, x0) rtol = 1e-12
        # time grid is dt·(i−1)
        @test res.t == [0.0, dt, 2dt]
        # detector-frame: h scaled by 1/dist_dimensionless = 1 (dist factor only applies to :source AAK)
        # and rotated to SSB. Reconstruct expected from the grounded frame primitives.
        hp0, hc0 = real.(hbase), -imag.(hbase)
        ehp, ehc = to_SSB_frame(hp0, hc0, qS, phiS, qK, phiK)
        @test res.h ≈ ehp .- im .* ehc rtol = 1e-12

        # return_list shape
        resl = generate_emri_waveform(:Pn5AAKWaveform, M, mu, a, p0, e0, x0, dist,
                                      qS, phiS, qK, phiK; dt = dt, return_list = true,
                                      aak_generator = stub)
        @test resl.hp ≈ ehp rtol = 1e-12
        @test resl.hc ≈ ehc rtol = 1e-12

        # source frame: no SSB rotation; h scaled UP by the inverse distance factor.
        ress = generate_emri_waveform(:Pn5AAKWaveform, M, mu, a, p0, e0, x0, dist,
                                      qS, phiS, qK, phiK; dt = dt, frame = :source,
                                      aak_generator = stub)
        dd = (dist * Gpc_in_m) / (mu * MRSUN_SI)
        hs = hbase .* dd
        @test ress.h ≈ real.(hs) .- im .* (-imag.(hs)) rtol = 1e-10
    end

    # =========================================================================
    # (C) FastSchwarzschildEccentricFlux source-frame generator
    # =========================================================================
    # Grounding:
    #  · CAPTURED FIXTURE — one-shot legacy KB fast_schwarzschild_eccentric_flux(M,mu,p0,e0,θ,φ).
    #    The t=0 sample is the implementation-independent IC anchor (Φ=0, p=p0,e0) → tight rtol;
    #    later samples differ only by the trajectory-integrator implementation → loose cross-check.
    #  · SELF-CONSISTENCY — generate_emri_waveform(:FastSchwarzschildEccentricFlux) reconstructs exactly
    #    from fast_schwarzschild_eccentric_flux composed with the grounded frame primitives.
    #  · FRAME-ROTATION INVARIANT — the SSB rotation preserves h₊²+h_×².
    @testset "FastSchwarzschildEccentricFlux" begin
        M, mu, p0, e0 = 1.0e6, 10.0, 12.0, 0.3
        theta, phi = 1.0, 0.5
        T, dt = 0.05, 10.0
        @test_throws ArgumentError fast_schwarzschild_eccentric_flux(
            M, mu, p0, e0, theta, phi; mode_selection = :bogus)
        @test_throws ArgumentError fast_schwarzschild_eccentric_flux(
            M, mu, p0, e0, theta, phi; T = 0.0)
        res = fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi; T = T, dt = dt)
        @test res.t[1] == 0.0 && res.t[2] == dt
        # IC anchor: implementation-independent
        @test res.h[1] ≈ 2.81078116907542086e-01 + 7.14190914278783062e-02im rtol = 1e-10
        # later samples: trajectory-implementation-limited cross-check vs captured KB
        kb = ComplexF64[
            2.81078116907542086e-01 + 7.14190914278783062e-02im,
            2.85216909613525027e-01 + 3.78304150799643399e-02im,
            2.84164672078787262e-01 + 4.43502684645670844e-03im,
            2.77445457624437619e-01 - 2.82202403119960477e-02im,
            2.64989086828170184e-01 - 5.90196044491937064e-02im,
        ]
        for i in 1:length(kb)
            @test res.h[i] ≈ kb[i] rtol = 2e-4
        end
        shifted = fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi;
                                                     T = 0.001, dt = dt,
                                                     Phi_phi0 = 0.7,
                                                     Phi_theta0 = 0.2,
                                                     Phi_r0 = 1.1)
        unshifted = fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi;
                                                       T = 0.001, dt = dt)
        @test shifted.h != unshifted.h
        custom_roman = RomanAmplitude(lmax = 2, nmax = 1)
        custom = fast_schwarzschild_eccentric_flux(M, mu, p0, e0, theta, phi;
                                                    T = 0.001, dt = dt,
                                                    roman = custom_roman)
        @test !isempty(custom.h) && all(isfinite, custom.h)

        # dispatch self-consistency: frame=:source ⇒ gen.h == fast_schwarz(θ,φ).h · (−1).
        a2, x02, dist = 0.0, 1.0, 1.0
        qS, phiS, qK, phiK = 1.2, 0.7, 0.9, 2.1
        th, ph = get_viewing_angles(qS, phiS, qK, phiK)
        base = fast_schwarzschild_eccentric_flux(M, mu, p0, e0, th, ph; T = 0.02, dt = dt)
        gens = generate_emri_waveform(:FastSchwarzschildEccentricFlux, M, mu, a2, p0, e0, x02, dist,
                                      qS, phiS, qK, phiK; frame = :source, T = 0.02, dt = dt)
        @test gens.t == base.t
        @test gens.h ≈ base.h .* (-1) rtol = 1e-12

        # detector frame: distance-scaled, π-rotated, SSB-rotated; rotation preserves h₊²+h_×².
        gend = generate_emri_waveform(:FastSchwarzschildEccentricFlux, M, mu, a2, p0, e0, x02, dist,
                                      qS, phiS, qK, phiK; frame = :detector, return_list = true,
                                      T = 0.02, dt = dt)
        ddf = (dist * Gpc_in_m) / (mu * MRSUN_SI)
        hpre = (base.h ./ ddf) .* (-1)            # pre-SSB-rotation source strain
        hp0, hc0 = real.(hpre), -imag.(hpre)
        @test gend.hp .^ 2 .+ gend.hc .^ 2 ≈ hp0 .^ 2 .+ hc0 .^ 2 rtol = 1e-10
    end

    # =========================================================================
    # (D) Pn5AAKWaveform default backend wired to AAK.aak_waveform (Module 4)
    # =========================================================================
    # Grounding: the default detector-frame generator now IS AAK.aak_waveform — assert the full
    # dispatch reconstructs exactly from aak_waveform composed with x0→Y0 + SSB frame rotation
    # (bit-exact self-consistency, NO Python). Equatorial Y0=1 (the only Trajectory-supported flux).
    @testset "Pn5AAKWaveform default backend (AAK.aak_waveform)" begin
        M, mu, a, p0, e0, x0, dist = 1.0e6, 10.0, 0.5, 12.0, 0.3, 1.0, 1.5
        qS, phiS, qK, phiK = 0.7, 0.8, 0.9, 1.0
        T, dt = 0.01, 10.0
        Y0 = xI_to_Y(a, p0, e0, x0)
        # raw AAK waveform (detector-frame generator output, h₊ − i h_×, distance applied inside)
        hraw = AAK.aak_waveform(M, mu, a, p0, e0, Y0, qS, phiS, qK, phiK, dist; T = T, dt = dt)
        ehp, ehc = to_SSB_frame(real.(hraw), -imag.(hraw), qS, phiS, qK, phiK)

        # default backend (no aak_generator kwarg) ⇒ uses AAK.aak_waveform.
        res = generate_emri_waveform(:Pn5AAKWaveform, M, mu, a, p0, e0, x0, dist,
                                     qS, phiS, qK, phiK; T = T, dt = dt)
        @test length(res.t) == length(hraw)
        @test res.t == [dt * (i - 1) for i in 1:length(hraw)]
        @test res.h ≈ ehp .- im .* ehc rtol = 1e-12       # bit-exact composition

        resl = generate_emri_waveform(:Pn5AAKWaveform, M, mu, a, p0, e0, x0, dist,
                                      qS, phiS, qK, phiK; T = T, dt = dt, return_list = true)
        @test resl.hp ≈ ehp rtol = 1e-12
        @test resl.hc ≈ ehc rtol = 1e-12
        # SSB rotation preserves h₊²+h_×².
        @test resl.hp .^ 2 .+ resl.hc .^ 2 ≈ real.(hraw) .^ 2 .+ imag.(hraw) .^ 2 rtol = 1e-10
    end
end

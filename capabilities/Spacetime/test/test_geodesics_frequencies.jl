# Geodesics/frequencies.jl — grounded by BAKED characterization vs the KerrGeodesics oracle
# (values captured once from KerrGeoRadialRoots/KerrGeoFrequencies/KerrGeoProperFrequencyFactor,
# Phase7-B) at the legacy runtests tolerance, confirming NumericsCore's ellipK/E/Pi reproduce the
# KerrGeodesics elliptic.jl path. Plus the analytic Schwarzschild Kepler limit.
using Test
using Spacetime
include(joinpath(@__DIR__, "fixtures", "spacetime_frequencies_ref.jl"))  # const FREQUENCIES_REF

@testset "Geodesics: radial roots + frequencies" begin
    cases = [
        (0.0,   10.0, 0.2, 1.0),  # Schwarzschild eccentric equatorial
        (0.9,   10.0, 0.0, 1.0),  # circular equatorial prograde
        (0.9,   12.0, 0.1, 0.5),  # generic inclined eccentric
        (0.5,    8.0, 0.3, 0.7),  # generic
        (-0.7,  11.0, 0.2, 0.6),  # retrograde generic
        (0.998,  7.0, 0.15, 0.8), # near-extremal generic
    ]

    @testset "radial_roots == KerrGeoRadialRoots (baked): (a,p,e,x)=$c" for (idx, c) in enumerate(cases)
        a, p, e, x = c
        r = radial_roots(a, p, e, x)
        ro = FREQUENCIES_REF[idx].roots
        for i in 1:4
            @test isapprox(real(r[i]), real(ro[i]); rtol = 1e-10, atol = 1e-12)
            @test isapprox(imag(r[i]), imag(ro[i]); rtol = 1e-10, atol = 1e-12)
        end
    end

    @testset "kerr_frequencies Mino + BL == KerrGeoFrequencies (baked): (a,p,e,x)=$c" for (idx, c) in enumerate(cases)
        a, p, e, x = c
        ref = FREQUENCIES_REF[idx]
        fm  = kerr_frequencies(a, p, e, x; Time = "Mino")
        fmo = ref.mino
        @test isapprox(fm.Ur,  fmo.Ur;  rtol = 1e-9, atol = 1e-11)
        @test isapprox(fm.Uth, fmo.Uth; rtol = 1e-9, atol = 1e-11)
        @test isapprox(fm.Uph, fmo.Uph; rtol = 1e-9, atol = 1e-11)
        @test isapprox(fm.Ut,  fmo.Ut;  rtol = 1e-9, atol = 1e-11)
        fb  = kerr_frequencies(a, p, e, x)                       # default BoyerLindquist
        fbo = ref.bl
        @test isapprox(fb.Omegar,  fbo.Omegar;  rtol = 1e-9, atol = 1e-11)
        @test isapprox(fb.Omegath, fbo.Omegath; rtol = 1e-9, atol = 1e-11)
        @test isapprox(fb.Omegaph, fbo.Omegaph; rtol = 1e-9, atol = 1e-11)
        @test isapprox(proper_frequency_factor(a, p, e, x),
                       ref.pff; rtol = 1e-9, atol = 1e-11)
    end

    @testset "analytic Schwarzschild circular Kepler: Ω_φ = p^{-3/2}" begin
        for p in (8.0, 10.0, 20.0)
            f = kerr_frequencies(0.0, p, 0.0, 1.0)
            @test isapprox(f.Omegaph, p^(-1.5); rtol = 1e-12)
        end
    end

    @testset "Ω_r → 0 approaching the Schwarzschild separatrix p → 6+2e" begin
        # At a=0 the radial frequency vanishes at the separatrix p = 6+2e (marginally bound).
        for e in (0.0, 0.2)
            psep = 6 + 2e
            Ωr_near = kerr_frequencies(0.0, psep + 1e-4, e, 1.0).Omegar
            Ωr_far  = kerr_frequencies(0.0, psep + 1e-2, e, 1.0).Omegar
            @test Ωr_near > 0
            @test Ωr_near < Ωr_far               # monotone descent toward the separatrix
            @test Ωr_near < 1e-2                  # vanishing (∝ 1/K(m→1), logarithmic for e>0)
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        a, p, e, x = big"0.9", big"12.0", big"0.1", big"0.5"
        f = kerr_frequencies(a, p, e, x; Time = "Mino")
        @test f.Ur isa BigFloat
        @test proper_frequency_factor(a, p, e, x) isa BigFloat
    end
end

# ---------------------------------------------------------------------------------------------
# ESSF::fix-kerr-frequencies-isco — kerr_frequencies AT the exact Kerr ISCO.
#
# Physics: on the circular-equatorial branch Υ_r vanishes identically at the ISCO, because the
# Υ_r radicand factors EXACTLY through the Bardeen–Press–Teukolsky ISCO polynomial
#     P_isco(p) = p² − 6p + 8a√p − 3a²,   Υ_r² = p·P_isco(p)/(2a√p + (p−3)p)
# (arXiv:2606.04998 main.tex:176-184, 207-209). Before the fix the Float64 evaluation of that
# radicand at p = kerr_isco(a) landed on the wrong side of zero (−2.2e−16 at a = 0.9, −1.7e−14
# at a = 0.998) and `sqrt` raised a DomainError; at a = 0.5 it landed at +1.3e−15 and returned a
# spurious Ω_r = 1.2e−9.
#
# The pins below were baked from the PRE-FIX code (GRTools 0252b52,
# evidence/fix-kerr-frequencies-isco/bake_regression_pin.jl) and are asserted BITWISE (`===`):
# the fix must not move any value away from the ISCO.
# ---------------------------------------------------------------------------------------------
@testset "kerr_frequencies at the exact ISCO" begin
    # (a, p, e, x) + PRE-FIX Mino/BL outputs, 12 circular-equatorial prograde points (the branch
    # that is edited) + 6 points routed through the other branches (retrograde flip, eccentric,
    # inclined, Schwarzschild). Baked with `repr` ⇒ bit-exact Float64 literals.
    freq_pin = [
    (a = 0.13724, p = 8.4712, e = 0.0, x = 1.0,
     Ur = 2.0797603922188186, Uth = 3.5520412034431446, Uph = 3.5908115783113126, Ut = 89.0267894360038,
     Omegar = 0.023361062500336913, Omegath = 0.03989856565586364, Omegaph = 0.0403340567604377),
    (a = 0.31859, p = 5.2037, e = 0.0, x = 1.0,
     Ur = 0.7416675518480136, Uth = 3.1396383935949976, Uph = 3.302341014694168, Ut = 40.25247444173742,
     Omegar = 0.01842539029299986, Omegath = 0.07799864324215397, Omegaph = 0.08204069589495847),
    (a = 0.47231, p = 22.9184, e = 0.0, x = 1.0,
     Ur = 4.487936504245138, Uth = 5.069045860303883, Uph = 5.1099486283251165, Ut = 563.0649645327285,
     Omegar = 0.007970548314916996, Omegath = 0.009002595046045067, Omegaph = 0.009075238116734392),
    (a = 0.58406, p = 4.1793, e = 0.0, x = 1.0,
     Ur = 0.7247122572613944, Uth = 2.798819777008711, Uph = 3.158627582859784, Ut = 28.831774269653675,
     Omegar = 0.025135888290585582, Omegath = 0.09707414295188051, Omegaph = 0.10955370118114223),
    (a = 0.66172, p = 11.3355, e = 0.0, x = 1.0,
     Ur = 2.9699255570180934, Uth = 3.7216316563801386, Uph = 3.8367969490427427, Ut = 148.96880961982322,
     Omegar = 0.019936559636862982, Omegath = 0.024982623314759325, Omegaph = 0.025755706572633992),
    (a = 0.72948, p = 3.6021, e = 0.0, x = 1.0,
     Ur = 0.7837572596606709, Uth = 2.567084836527457, Uph = 3.0765664530674766, Ut = 23.277232014332803,
     Omegar = 0.033670552374014125, Omegath = 0.11028307983297977, Omegaph = 0.13217063142099975),
    (a = 0.81095, p = 47.6218, e = 0.0, x = 1.0,
     Ur = 6.718738399056085, Uth = 7.078287527718077, Uph = 7.110358903867743, Ut = 2342.451480641355,
     Omegar = 0.0028682508280668925, Omegath = 0.0030217434965953134, Omegaph = 0.0030354348692511455),
    (a = 0.88263, p = 2.9074, e = 0.0, x = 1.0,
     Ur = 0.8686611137443341, Uth = 2.2494950254375876, Uph = 2.994494249617334, Ut = 17.48804620571598,
     Omegar = 0.049671707378061, Omegath = 0.12863043698399762, Omegaph = 0.17123092050377714),
    (a = 0.92517, p = 6.7439, e = 0.0, x = 1.0,
     Ur = 2.2051199999783058, Uth = 2.936876876155293, Uph = 3.1946122384672213, Ut = 58.90359388691754,
     Omegar = 0.037436085889965734, Omegath = 0.049859043945492974, Omegaph = 0.05423458956681323),
    (a = 0.95681, p = 2.1586, e = 0.0, x = 1.0,
     Ur = 0.6709975341252422, Uth = 1.966459333217405, Uph = 3.1789674136327295, Ut = 13.123613074199985,
     Omegar = 0.05112902447911786, Omegath = 0.14984130681841823, Omegaph = 0.24223263789164395),
    (a = 0.97904, p = 15.0832, e = 0.0, x = 1.0,
     Ur = 3.615570633880379, Uth = 4.134495010448294, Uph = 4.251340173494934, Ut = 253.20050238836316,
     Omegar = 0.014279476540432594, Omegath = 0.016328936836415658, Omegaph = 0.01679040970848532),
    (a = 0.99312, p = 1.7245, e = 0.0, x = 1.0,
     Ur = 0.6540044796645852, Uth = 1.738182424424621, Uph = 3.5422066810929755, Ut = 11.539578411415638,
     Omegar = 0.05667490235324412, Omegath = 0.15062789665739498, Omegaph = 0.3069615331517496),
    (a = -0.62145, p = 9.8317, e = 0.0, x = 1.0,
     Ur = 1.803208443451471, Uth = 4.051152837090362, Uph = 3.875643814648731, Ut = 117.06916208525787,
     Omegar = 0.015402932858938973, Omegath = 0.03460478203593899, Omegaph = 0.03310559113617144),
    (a = -0.91038, p = 12.4402, e = 0.0, x = 1.0,
     Ur = 2.4125999439958963, Uth = 4.36574867694318, Uph = 4.164359266034248, Ut = 178.9301746342473,
     Omegar = 0.013483471689040223, Omegath = 0.024399175185891613, Omegaph = 0.023273655628776146),
    (a = 0.77219, p = 10.6538, e = 0.23, x = 1.0,
     Ur = 2.8961861786480054, Uth = 3.6048977027943403, Uph = 3.745362857313548, Ut = 142.33100071167135,
     Omegar = 0.02034824573822106, Omegath = 0.025327565216076875, Omegaph = 0.026314456011594826),
    (a = 0.44107, p = 13.7726, e = 0.0, x = 0.61,
     Ur = 3.246099609690221, Uth = 4.125251720379909, Uph = 4.1935693129176705, Ut = 214.22286253987153,
     Omegar = 0.015152909316978488, Omegath = 0.019256822878146864, Omegaph = 0.019575731848589017),
    (a = -0.53964, p = 11.2891, e = 0.17, x = 0.42,
     Ur = 2.576540031473297, Uth = 4.006015810588084, Uph = 3.887761350368091, Ut = 154.2201459993553,
     Omegar = 0.01670689659108524, Omegath = 0.025975956543348305, Omegaph = 0.02520916657921166),
    (a = 0.0, p = 9.3158, e = 0.31, x = 1.0,
     Ur = 2.2137387587582986, Uth = 3.735385575461579, Uph = 3.735385575461579, Ut = 117.32688934065737,
     Omegar = 0.018868127938947838, Omegath = 0.03183742104178631, Omegaph = 0.03183742104178631),
    ]
    # r0 = 1.001*kerr_isco(a): the ISCO neighbourhood must stay bitwise continuous with pre-fix.
    near_isco_pin = [
    (a = 0.5, r0 = 4.237235532060357,
     Ur = 0.09195761108374322, Uth = 2.909641706788739, Uph = 3.228013344692251, Ut = 29.769246410505552,
     Omegar = 0.0030890137363803545, Omegath = 0.0977398509410009, Omegaph = 0.10843449982506399),
    (a = 0.9, r0 = 2.3232039248036487,
     Ur = 0.06807201058016861, Uth = 2.154477674746914, Uph = 3.271966780175446, Ut = 14.530939789871574,
     Omegar = 0.0046846254656988185, Omegath = 0.14826829550616116, Omegaph = 0.22517241331191035),
    (a = 0.998, r0 = 1.2382076258303598,
     Ur = 0.04960608212214662, Uth = 1.5728809517179516, Uph = 6.925534650100786, Ut = 16.453773792450868,
     Omegar = 0.003014875660008546, Omegath = 0.09559393313402685, Omegaph = 0.4209085853166573),
    (a = -0.5, r0 = 7.562139299226869,
     Ur = 0.12286272438718854, Uth = 3.8870532990184095, Uph = 3.690607763757769, Ut = 74.90231329956978,
     Omegar = 0.0016403061397556892, Omegath = 0.051894969965376694, Omegaph = 0.04927228013635951),
    (a = -0.9, r0 = 8.726069631886094,
     Ur = 0.13198232463322543, Uth = 4.175490381949952, Uph = 3.857648442263387, Ut = 95.96572196639099,
     Omegar = 0.0013753069526163534, Omegath = 0.043510227364436306, Omegaph = 0.04019819122096959),
    (a = -0.998, r0 = 9.003368829258372,
     Ur = 0.1340635445336869, Uth = 4.241316378737122, Uph = 3.896810161307225, Ut = 101.38393793846873,
     Omegar = 0.0013223351475561333, Omegath = 0.04183420436194966, Omegaph = 0.03843616888971359),
    ]
    spins = (0.5, 0.9, 0.998, -0.5, -0.9, -0.998)

    @testset "no throw + Ω_r ≡ 0 + Ω_φ vs sealed u^φ/u^t: a=$a" for a in spins
        r0 = kerr_isco(a)
        fb = kerr_frequencies(a, r0, 0.0, 1.0)                       # must not throw
        fm = kerr_frequencies(a, r0, 0.0, 1.0; Time = "Mino")
        fp = kerr_frequencies(a, r0, 0.0, 1.0; Time = "Proper")
        @test all(isfinite, (fb.Omegar, fb.Omegath, fb.Omegaph))
        @test all(isfinite, (fm.Ur, fm.Uth, fm.Uph, fm.Ut))
        @test all(isfinite, (fp.Omegar, fp.Omegath, fp.Omegaph))
        if a > 0
            # prograde ⇒ the fixed `_mino_freqs_circeq` branch: Υ_r and Ω_r are EXACTLY zero.
            @test fm.Ur == 0
            @test fb.Omegar == 0
            @test fp.Omegar == 0
        else
            # a < 0 flips to (−a, x = −1) and routes through the GENERIC branch, which is not
            # touched by this fix: there r3 → r1 = r2 and `sqrt(complex(...))` leaves a residual
            # instead of 0. Recorded, not fixed (out of node scope). Measured: BL Ω_r =
            # 1.53/2.44/2.31e−9 and Proper Ω_r = 2.00/3.07/2.89e−9 over the three spins.
            @test abs(fb.Omegar) <= 1e-8
            @test abs(fp.Omegar) <= 1e-8
        end
        # Ω_φ against the independently sealed closed form (ESSF::spacetime-circular-kinematics).
        ut, uphi = kerr_circular_ut_uphi(a, r0)
        @test isapprox(fb.Omegaph, uphi/ut; rtol = 1e-12)
        @test isapprox(fm.Uph/fm.Ut, uphi/ut; rtol = 1e-12)
    end

    @testset "Schwarzschild ISCO p = 6 stays exact" begin
        f = kerr_frequencies(0.0, kerr_isco(0.0), 0.0, 1.0)
        @test kerr_isco(0.0) == 6
        @test f.Omegar == 0
        @test isapprox(f.Omegaph, 6.0^(-1.5); rtol = 1e-14)
    end

    @testset "prograde ISCO spin scan: Ω_r ≡ 0 for 201 spins" begin
        nonzero = 0
        for a in range(0.001, 0.9995; length = 201)
            kerr_frequencies(a, kerr_isco(a), 0.0, 1.0).Omegar == 0 || (nonzero += 1)
        end
        @test nonzero == 0
    end

    @testset "BigFloat(256) inherits the exact ISCO zero" begin
        setprecision(BigFloat, 256) do
            for a in (big"0.5", big"0.9", big"0.998")
                r0 = kerr_isco(a)
                fb = kerr_frequencies(a, r0, big"0.0", big"1.0")
                fm = kerr_frequencies(a, r0, big"0.0", big"1.0"; Time = "Mino")
                @test fb.Omegar isa BigFloat
                @test fb.Omegaph isa BigFloat
                @test fm.Ur isa BigFloat
                @test fb.Omegar == 0                       # exact, not merely ≤ 1e-60
                ut, uphi = kerr_circular_ut_uphi(a, r0)
                @test abs(fb.Omegaph - uphi/ut) <= big"1e-60" * abs(uphi/ut)
            end
            # retrograde: generic branch, untouched — residual measured at 2.9e-40 (a = −0.5),
            # exactly 0 at a = −0.9 and a = −0.998
            for a in (big"-0.5", big"-0.9", big"-0.998")
                fb = kerr_frequencies(a, kerr_isco(a), big"0.0", big"1.0")
                @test fb.Omegar isa BigFloat
                @test abs(fb.Omegar) <= big"1e-35"
            end
        end
    end

    @testset "Υ_r radicand == p·P_isco(p)/(2a√p+(p−3)p) — Bardeen factorization" begin
        setprecision(BigFloat, 256) do
            for a in (big"0.3", big"0.7", big"0.95"), fac in (big"1.001", big"1.5", big"4.0")
                p = fac * kerr_isco(a)
                Ur = kerr_frequencies(a, p, big"0.0", big"1.0"; Time = "Mino").Ur
                Pisco = p^2 - 6*p + 8*a*sqrt(p) - 3*a^2
                @test isapprox(Ur^2, p*Pisco/(2*a*sqrt(p) + (p - 3)*p); rtol = big"1e-70")
            end
        end
    end

    @testset "inside the ISCO still raises DomainError (no blanket clamp)" for a in (0.5, 0.9, 0.998)
        for frac in (0.999, 0.99, 0.95)
            @test_throws DomainError kerr_frequencies(a, frac*kerr_isco(a), 0.0, 1.0)
        end
    end

    @testset "r0 = 1.001*ISCO bitwise unchanged: a=$(row.a)" for row in near_isco_pin
        fm = kerr_frequencies(row.a, row.r0, 0.0, 1.0; Time = "Mino")
        fb = kerr_frequencies(row.a, row.r0, 0.0, 1.0)
        @test fm.Ur  === row.Ur
        @test fm.Uth === row.Uth
        @test fm.Uph === row.Uph
        @test fm.Ut  === row.Ut
        @test fb.Omegar  === row.Omegar
        @test fb.Omegath === row.Omegath
        @test fb.Omegaph === row.Omegaph
    end

    @testset "away-from-ISCO regression pin, bitwise: (a,p,e,x)=($(row.a),$(row.p),$(row.e),$(row.x))" for row in freq_pin
        fm = kerr_frequencies(row.a, row.p, row.e, row.x; Time = "Mino")
        fb = kerr_frequencies(row.a, row.p, row.e, row.x)
        @test fm.Ur  === row.Ur
        @test fm.Uth === row.Uth
        @test fm.Uph === row.Uph
        @test fm.Ut  === row.Ut
        @test fb.Omegar  === row.Omegar
        @test fb.Omegath === row.Omegath
        @test fb.Omegaph === row.Omegaph
    end

    @testset "branch routing: a<0 and e≠0/|x|≠1 bypass the circular-equatorial closed forms" begin
        # `_mino_frequencies` sends a < 0 to (−a, −x) and only (a>0, e==0, x==1) reaches
        # `_mino_freqs_circeq`, so the retrograde and generic branches see no edit at all
        # (their values are pinned bitwise above). Assert the flip identity that routes them.
        a, p = 0.7, 11.0
        fneg  = Spacetime._mino_frequencies(-a, p, 0.0, 1.0)
        fflip = Spacetime._mino_frequencies(a, p, 0.0, -1.0)
        @test fneg.Ur  === fflip.Ur
        @test fneg.Uth === fflip.Uth
        @test fneg.Ut  === fflip.Ut
        @test fneg.Uph === -fflip.Uph
        # x = −1 is NOT the circ-eq branch: it differs from the prograde closed form at the same p
        @test fneg.Uph != real(Spacetime._mino_freqs_circeq(a, p).Uph)
    end
end

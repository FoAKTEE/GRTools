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

# Module 2 (AmplitudesFluxes): energy_flux / angular_momentum_flux from Z^∞/Z^H.
# GROUNDING: (a) LIVE Teukolsky oracle — same closed-form on the same (s,l,m,a,ω,λ,Z), matched at
# the legacy tolerance (Teukolsky/test/runtests.jl: ℐ rtol 1e-7, ℋ rtol 1e-3 atol 1e-10); (b) the
# analytic identity JFlux/EFlux = m/ω exact (rtol 1e-12, legacy line 411).
# Legacy oracle RHS BAKED once into fixtures/leg_flux.jl (Phase 7-B); tolerances unchanged.
include("fixtures/leg_flux.jl")  # const REF_FLUX_LEG

@testset "AmplitudesFluxes :: energy/angular_momentum_flux vs Teukolsky oracle (baked from Teukolsky.*_flux, captured once)" begin
    # sample (s,l,m,a,ω,λ) + arbitrary complex Z^∞/Z^H (the flux is a pure function of these).
    samples = [(-2, 2, 2, 0.0, 0.1, 4.0, 1.2 + 0.3im, 0.04 - 0.02im),
               (-2, 2, 1, 0.5, 0.15, 3.8, -0.7 + 1.1im, 0.01 + 0.005im),
               (-1, 2, 2, 0.6, 0.2, 5.1, 0.9 - 0.4im, 0.02 - 0.01im),
               (0, 2, 2, 0.5, 0.2, 6.0, 2.0 + 0.0im, 0.1 + 0.0im),
               (1, 2, 1, 0.9, 0.2, 4.7, 0.3 + 0.6im, 0.008 - 0.003im),
               (2, 3, 3, 0.7, 0.25, 9.0, 0.5 - 0.2im, 0.003 + 0.001im)]
    for (s, l, m, a, ω, λ, ZI, ZH) in samples
        leg = REF_FLUX_LEG[(s, l, m, a, ω)]
        efn = energy_flux(s, l, m, a, ω, λ, ZI, ZH)
        @test isapprox(efn.I, leg.eI; rtol=1e-7)
        @test isapprox(efn.H, leg.eH; rtol=1e-3, atol=1e-10)
        jfn = angular_momentum_flux(s, l, m, a, ω, λ, ZI, ZH)
        @test isapprox(jfn.I, leg.jI; rtol=1e-7)
        @test isapprox(jfn.H, leg.jH; rtol=1e-3, atol=1e-10)
        # analytic identity: J/E = m/ω (exact in the flux formula)
        @test isapprox(jfn.I, efn.I * m / ω; rtol=1e-12)
        @test isapprox(jfn.H, efn.H * m / ω; rtol=1e-12)
    end
end

@testset "AmplitudesFluxes :: flux ω=0 short-circuit" begin
    @test energy_flux(-2, 2, 2, 0.5, 0.0, 4.0, 1.0 + 0.0im, 0.0im) == (I=0.0, H=0.0)
    @test angular_momentum_flux(-2, 2, 2, 0.5, 0.0, 4.0, 1.0 + 0.0im, 0.0im) == (I=0.0, H=0.0)
end

@testset "AmplitudesFluxes :: static Wronskian uses the static infinity amplitude" begin
    # BHPToolkit Teukolsky/Kernel/ConvolveSource.m: static a*m==0 uses I, otherwise I+.
    for (a, m) in ((0.0, 1), (0.5, 0), (0.5, 1))
        Rin = TeukolskyRadial(-2, 2, m, a, 0.0).In
        amp = iszero(a * m) ? Rin.amplitudes.I : Rin.amplitudes.Ip
        W = Perturbation._teuk_wronskian(Rin)
        @test W == -(2 * Rin.l + 1) * amp
        @test isfinite(W)
        @test !iszero(W)
    end
end

# Module 3 (AmplitudesFluxes end-to-end): the public TeukolskyPointParticleMode* entrypoints.
# GROUNDING: LIVE Teukolsky oracle end-to-end from (a,p[,e][,x]) — Zinf rtol 1e-7, fluxes ℐ rtol
# 1e-7, orbit scalars/ω/λ tight, ℋ fluxes rtol 1e-3 atol 1e-10 (Teukolsky/test/runtests.jl lines
# 332-478). PLUS the e→0 circular-limit identity, the Newtonian-quadrupole 32/5 limit, the
# JFlux/EFlux=m/ω identity, and parallel==serial determinism.
# Legacy oracle Mode fields BAKED once into fixtures/leg_pointparticle.jl (Phase 7-B); tolerances unchanged.
include("fixtures/leg_pointparticle.jl")  # REF_PP_CIRC / REF_PP_ECC / REF_PP_SPH / REF_PP_GEN

@testset "PointParticle :: Circular vs Teukolsky oracle (grav + non-grav spins, baked, captured once)" begin
    for (s, l, m, a, p) in [(-2, 2, 2, 0.0, 8.0), (-2, 2, 1, 0.5, 8.0), (-2, 3, 3, 0.9, 6.0),
                            (0, 2, 2, 0.5, 8.0), (-1, 2, 2, 0.6, 8.0), (2, 2, 2, 0.5, 8.0)]
        mode = TeukolskyPointParticleModeCircular(s, l, m, a, p)
        ref = REF_PP_CIRC[(s, l, m, a, p)]
        @test isapprox(mode.En, ref.En; rtol=1e-12)
        @test isapprox(mode.Lz, ref.Lz; rtol=1e-12)
        @test isapprox(mode.Upst, ref.Upst; rtol=1e-10)
        @test isapprox(mode.omega, ref.omega; rtol=1e-10)
        @test isapprox(mode.lambda, ref.lambda; rtol=1e-8)
        @test isapprox(mode.Zinf, ref.Zinf; rtol=1e-7, atol=1e-12)
        @test isapprox(mode.EnergyFlux.I, ref.efI; rtol=1e-7, atol=1e-18)
        @test isapprox(mode.EnergyFlux.H, ref.efH; rtol=1e-3, atol=1e-10)
        @test isapprox(mode.AngularMomentumFlux.I, ref.jfI; rtol=1e-7, atol=1e-18)
        # analytic identity J/E = m/ω
        @test isapprox(mode.AngularMomentumFlux.I, mode.EnergyFlux.I * m / mode.omega; rtol=1e-12)
    end
end

@testset "PointParticle :: Eccentric vs Teukolsky oracle (s=-2, +2, 0)" begin
    for (s, l, m, n, a, p, e) in [(-2, 2, 2, 0, 0.5, 8.0, 0.1), (-2, 2, 1, 1, 0.5, 8.0, 0.2),
                                  (2, 2, 2, 0, 0.5, 8.0, 0.1), (0, 2, 2, -1, 0.5, 8.0, 0.2)]
        md = TeukolskyPointParticleModeEccentric(l, m, n, a, p, e; s=s)
        ref = REF_PP_ECC[(s, l, m, n, a, p, e)]
        @test isapprox(md.omega, ref.omega; rtol=1e-10)
        @test isapprox(md.lambda, ref.lambda; rtol=1e-8, atol=1e-10)
        @test isapprox(md.Zinf, ref.Zinf; rtol=1e-7, atol=1e-12)
        @test isapprox(md.EnergyFlux.I, ref.efI; rtol=1e-7, atol=1e-20)
        @test isapprox(md.AngularMomentumFlux.I, ref.jfI; rtol=1e-7, atol=1e-20)
    end
end

@testset "PointParticle :: Spherical + Generic vs Teukolsky oracle (s=-2)" begin
    sph = TeukolskyPointParticleModeSpherical(2, 2, 0, 0.5, 8.0, 0.9)
    rsph = REF_PP_SPH
    @test isapprox(sph.omega, rsph.omega; rtol=1e-10)
    @test isapprox(sph.Zinf, rsph.Zinf; rtol=1e-7, atol=1e-12)
    @test isapprox(sph.EnergyFlux.I, rsph.efI; rtol=1e-7, atol=1e-20)

    gen = TeukolskyPointParticleModeGeneric(2, 2, 0, 0, 0.5, 8.0, 0.1, 0.9)
    rgen = REF_PP_GEN
    @test isapprox(gen.omega, rgen.omega; rtol=1e-10)
    @test isapprox(gen.Zinf, rgen.Zinf; rtol=1e-7, atol=1e-12)
    @test isapprox(gen.EnergyFlux.I, rgen.efI; rtol=1e-7, atol=1e-20)
end

@testset "PointParticle :: eccentric e→0 reduces to circular n=0" begin
    a, p, l, m = 0.5, 8.0, 2, 2
    circ = TeukolskyPointParticleModeCircular(-2, l, m, a, p)
    ecc = TeukolskyPointParticleModeEccentric(l, m, 0, a, p, 0.002)
    @test isapprox(ecc.Zinf, circ.Zinf; rtol=1e-3)
    @test isapprox(ecc.EnergyFlux.I, circ.EnergyFlux.I; rtol=2e-3)
end

@testset "PointParticle :: public static m=0 mode has finite amplitudes and zero flux" begin
    mode = TeukolskyPointParticleModeCircular(-2, 2, 0, 0.5, 8.0)
    @test iszero(mode.omega)
    @test isfinite(mode.Zinf)
    @test isfinite(mode.ZH)
    @test mode.EnergyFlux == (I=0.0, H=0.0)
    @test mode.AngularMomentumFlux == (I=0.0, H=0.0)
end

@testset "PointParticle :: Newtonian-quadrupole 32/5 limit (Schwarzschild circular l=m=2)" begin
    # Leading-PN total Edot of a circular Schwarzschild orbit → (32/5)(M/p)^5 (μ=M=1). The l=2,m=±2
    # quadrupole pair dominates; at large p the residual is the O(M/p) PN + higher-l correction.
    p = 80.0
    e22 = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.0, p).EnergyFlux.I
    e2m2 = TeukolskyPointParticleModeCircular(-2, 2, -2, 0.0, p).EnergyFlux.I
    newt = (32 / 5) * p^(-5)
    @test isapprox(e22 + e2m2, newt; rtol=5e-2)
end

@testset "PointParticle :: parallel mode-sum == serial (determinism)" begin
    a, p = 0.5, 8.0
    modes = [(l, m) for l in 2:3 for m in vcat(-l:-1, 1:l)]
    n = length(modes)
    efI = zeros(n)
    Threads.@threads for i in 1:n
        l, m = modes[i]
        efI[i] = TeukolskyPointParticleModeCircular(-2, l, m, a, p).EnergyFlux.I
    end
    efI_ser = [TeukolskyPointParticleModeCircular(-2, l, m, a, p).EnergyFlux.I for (l, m) in modes]
    @test efI == efI_ser
end

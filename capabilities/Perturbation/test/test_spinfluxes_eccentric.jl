# SpinFluxes eccentric-equatorial linear-in-spin s=-2 flux. GROUNDING:
# (1) PRIMARY characterization vs REF_EQ_SPIN_FLUX (ref_eq_spin_flux.jl) — the BHPToolkit
#     Mathematica/JSON UPSTREAM oracle — at the FROZEN legacy tolerances copied verbatim from
#     TeukolskySpinFluxes/test/runtests.jl: geodesic Energy.I/AngularMomentum.I rtol=1e-6,
#     FluxesCorrection dE1I/dL1I rtol=2e-3, horizon dE1H/dL1H rtol=5e-3 atol=1e-10.
# (2) IN-PACKAGE reduction gate: equatorial_spin_flux(l,m,0,a,p,1e-5) == circular_spin_flux(l,m,a,p)
#     (dE1I/dL1I/CPlus0/CMinus0 rtol=1e-4) — the e→0 circular limit of the wr Mino sum.
# (3) PARALLEL determinism: Threads.@threads n-sum of dE1I bit-identical to the serial vector.
# (4) ANALYTIC: single-mode GW flux identity L̇ = (m/ω) Ė on the geodesic I-flux.
using Perturbation
using Perturbation: circular_spin_flux, equatorial_spin_flux

include("ref_eq_spin_flux.jl")

@testset "SpinFluxes eccentric :: TeukolskySpinModeEquatorialCorrectionAnalyticalPhases vs MMA oracle (6 modes)" begin
    for r in REF_EQ_SPIN_FLUX
        md = TeukolskySpinModeEquatorialCorrectionAnalyticalPhases(r.l, r.m, r.n, r.a, r.p, r.e)
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=1e-6)            # geodesic Edot ∞ (exact wr sum)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=1e-6)   # geodesic Ldot ∞
        # linear-in-spin: FD-limited; analytic ω-derivatives → ~1e-7 is [FUTURE]
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=2e-3)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=2e-3)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=5e-3, atol=1e-10)
        @test isapprox(md.FluxesCorrection.AngularMomentum.H, r.dL1H; rtol=5e-3, atol=1e-10)
    end
end

@testset "SpinFluxes eccentric :: e→0 reduction to circular_spin_flux" begin
    # Validates the eccentric wr-sum machinery (2 velocity branches, phases, FabEq+∂_ωF+G_a, W/W1)
    # by its reduction to the oracle-verified circular flux as e→0.
    for (l, m, a, p) in ((2,2,0.5,10.0), (2,2,0.0,10.0), (3,3,0.5,10.0), (2,1,0.9,8.0))
        c  = circular_spin_flux(l, m, a, p)
        ec = equatorial_spin_flux(l, m, 0, a, p, 1e-5)
        @test isapprox(ec.dE1I, c.dE1I; rtol=1e-4)
        @test isapprox(ec.dL1I, c.dL1I; rtol=1e-4)
        @test isapprox(ec.CPlus0, c.CPlus0; rtol=1e-4)
        @test isapprox(ec.CMinus0, c.CMinus0; rtol=1e-4)
    end
end

@testset "SpinFluxes eccentric :: PARALLEL n-sum deterministic (Threads.@threads == serial)" begin
    # The physical eccentric flux Σ_n dE¹_{lmn}: parallel-runnable (per-n on its own thread into a
    # preallocated slot) and bit-identical to the serial vector. n=0 cross-checked vs the oracle.
    l, m, a, p, e = 2, 2, 0.5, 10.0, 0.2
    nrange = collect(0:4)
    dE = zeros(length(nrange))
    Threads.@threads for i in eachindex(nrange)
        dE[i] = equatorial_spin_flux(l, m, nrange[i], a, p, e).dE1I
    end
    dE_ser = [equatorial_spin_flux(l, m, n, a, p, e).dE1I for n in nrange]
    @test dE == dE_ser                                       # bit-identical determinism
    @test isapprox(dE[1], 1.877227537378918e-5; rtol=2e-3)   # n=0 vs BHPToolkit oracle
    @test dE[1] > sum(dE[2:end])                             # n=0 fundamental dominates
end

@testset "SpinFluxes eccentric :: analytic single-mode flux identity L̇=(m/ω)Ė" begin
    for (l, m, n, a, p, e) in ((2,2,0,0.5,10.0,0.2), (3,3,0,0.5,10.0,0.2), (2,2,1,0.5,10.0,0.2))
        md = TeukolskySpinModeEquatorialCorrectionAnalyticalPhases(l, m, n, a, p, e)
        @test isapprox(md.Fluxes.AngularMomentum.I, (m / md.omega) * md.Fluxes.Energy.I; rtol=1e-10)
        @test isapprox(md.Fluxes.AngularMomentum.H, (m / md.omega) * md.Fluxes.Energy.H; rtol=1e-10)
    end
end

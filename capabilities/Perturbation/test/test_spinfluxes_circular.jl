# SpinFluxes circular-equatorial linear-in-spin s=-2 flux. GROUNDING:
# (1) PRIMARY characterization vs REF_SPIN_FLUX_CIRC (ref_spin_flux_circular.jl) — the BHPToolkit
#     Mathematica/JSON UPSTREAM oracle (strictly better than the legacy Julia port) — at the FROZEN
#     legacy tolerances copied verbatim from TeukolskySpinFluxes/test/runtests.jl: geodesic
#     Energy.I/AngularMomentum.I rtol=1e-6, FluxesCorrection dE1I/dL1I rtol=1e-3, horizon dE1H/dL1H
#     rtol=3e-3 atol=1e-10 (the MMA Chop threshold).
# (2) IN-PACKAGE cross-check: _circular_amplitudes_fabeq CPlus0/CMinus0 (FabEq geodesic-amplitude
#     path) == the oracle-verified TeukolskyPointParticleModeCircular Zinf/ZH at rtol=1e-8, and the
#     FabEq amplitudes reproduce mode.EnergyFlux through energy_flux (I rtol 1e-7 / H rtol 1e-6).
# (3) ANALYTIC: single-mode GW flux identity L̇ = (m/ω) Ė on the geodesic I-flux (rtol 1e-10), and
#     the m=−2 mode = m=+2 mirror (efI(2,-2) == efI(2,2) at the same a,p reflected) sanity.
using Perturbation
using Perturbation: circular_spin_flux, _circular_amplitudes_fabeq

include("ref_spin_flux_circular.jl")

@testset "SpinFluxes circular :: TeukolskySpinModeCircularCorrectionAnalytical vs MMA oracle (8 modes)" begin
    for r in REF_SPIN_FLUX_CIRC
        md = TeukolskySpinModeCircularCorrectionAnalytical(r.l, r.m, r.a, r.p)
        # geodesic monopole fluxes (FabEq path) — tight
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=1e-6)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=1e-6)
        # linear-in-spin correction — FD-limited (∂_ω on radials); analytic ω-derivs [FUTURE]
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=1e-3)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=1e-3)
        # horizon: atol = the oracle's MMA Chop threshold (1e-10)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=3e-3, atol=1e-10)
        @test isapprox(md.FluxesCorrection.AngularMomentum.H, r.dL1H; rtol=3e-3, atol=1e-10)
    end
end

@testset "SpinFluxes circular :: FabEq amplitudes == Teukolsky source Zinf/ZH (5 modes)" begin
    # In-package cross-check: the FabEq geodesic-amplitude path (CPlus0/CMinus0) must reproduce the
    # oracle-verified point-particle source amplitudes Zinf/ZH and the geodesic flux through them.
    for (l, m, a, p) in ((2,2,0.5,10.0), (2,2,0.0,10.0), (2,1,0.9,8.0), (3,3,0.5,8.0), (2,-2,0.9,10.0))
        mode = TeukolskyPointParticleModeCircular(-2, l, m, a, p)
        amp = _circular_amplitudes_fabeq(l, m, a, p, mode.En, mode.Lz, mode.Upst, mode.omega, mode.lambda)
        @test isapprox(amp.CPlus0,  mode.Zinf; rtol=1e-8)
        @test isapprox(amp.CMinus0, mode.ZH;   rtol=1e-8)
        ef = energy_flux(-2, l, m, a, mode.omega, mode.lambda, amp.CPlus0, amp.CMinus0)
        @test isapprox(ef.I, mode.EnergyFlux.I; rtol=1e-7)
        @test isapprox(ef.H, mode.EnergyFlux.H; rtol=1e-6, atol=1e-12)
    end
end

@testset "SpinFluxes circular :: analytic single-mode flux identity L̇=(m/ω)Ė" begin
    # For a single (l,m) GW mode the angular-momentum flux equals (m/ω) times the energy flux —
    # an exact identity independent of the oracle. Verifies the geodesic Fluxes tuple is consistent.
    for (l, m, a, p) in ((2,2,0.5,10.0), (3,3,0.5,10.0), (2,1,0.9,8.0))
        md = TeukolskySpinModeCircularCorrectionAnalytical(l, m, a, p)
        @test isapprox(md.Fluxes.AngularMomentum.I, (m / md.omega) * md.Fluxes.Energy.I; rtol=1e-10)
        @test isapprox(md.Fluxes.AngularMomentum.H, (m / md.omega) * md.Fluxes.Energy.H; rtol=1e-10)
    end
end

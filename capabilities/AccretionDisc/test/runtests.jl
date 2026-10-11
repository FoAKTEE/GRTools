using Test
using AccretionDisc

@testset "AccretionDisc" begin
    @testset "smoke" begin
        @test isdefined(AccretionDisc, :Modes)
    end

    include("hygiene.jl")
    include("test_forcing_modes.jl")
    include("test_eccentric_newtonian.jl")       # ECC::fr22-kepler-orbit, ECC::fr22-zz22-overlay
    include("test_newtonian_projector.jl")       # ECC::shared-newtonian-projector: all-mode potential bank (oracle sets: GRTOOLS_SLOW_TESTS=1)
    include("Modes/test_fr22_potential_harmonics.jl") # ECC::fr22-potential-harmonics: Phi_lm gates (cache audit: FR22_STAGE=cache)
    include("test_literature_prescriptions.jl")  # ECC::fr24-literature-prescriptions: FR24 v1 (52)-(55)

    include("Modes/test_background.jl")
    include("Modes/test_velocity.jl")
    include("Modes/test_offequator.jl")
    include("Modes/test_masterpde.jl")
    include("Modes/test_resonances.jl")
    include("Modes/test_reducedode.jl")
    include("Modes/test_forcing.jl")
    include("Modes/test_negative_spin_torque.jl")  # ECC::st-negative-spin-constants: a<0 ILR torque consumer
    include("Modes/test_pressureless_lindblad.jl")  # EGD::analytic-lindblad-hegade: arXiv:2510.03564 (18a,18b,19,20,23,24,25,C1-C2)
    include("Modes/test_newtonian_planet.jl")  # MR19/MR20 Newtonian planet-forced disc modes (arXiv:1811.09628)
    include("Modes/test_phase_gradient.jl")  # ECC::mr19-fig15: MR Appendix A secant step
    include("Modes/test_fr22_phase_gradient.jl")  # ECC::fr22-phase-gradient: opt-in descent safeguard of the secant step
    include("Modes/test_solver_mode_generalisation.jl")  # ECC::fr22-solver-mode-generalisation: frequency-aware (m,l) MR solver
    include("test_mr_core_certification.jl")   # ECC::mr19-core-certification: thermal endpoints
    include("Modes/test_thindisc.jl")   # NT/SS thin-disc structure + total Lindblad torque (2510.02433)
    include("Force/test_loader.jl")
    include("Force/test_assemble.jl")
    include("Force/test_source_fluid.jl")
    include("Force/test_pipeline.jl")
    include("Force/test_source_grid.jl")
    include("Force/test_bigfloat.jl")
    include("Force/test_kerr_frontend_numerics.jl")
    include("Force/test_kerr.jl")
    include("SpiralWaves/test_scaffold.jl")   # arXiv:2608.17003 spiral-wave submodule: load + dep budget
    include("SpiralWaves/test_background.jl") # EGD::sw-background: disc equilibrium (Sec. II.A)
    include("SpiralWaves/test_grid.jl")       # EGD::sw-grid: adaptive corotation/WKB grid (App. C)
    include("SpiralWaves/test_master.jl")     # EGD::sw-master-matrices: master system C, S (47), (50), Landau (62-63)
    include("SpiralWaves/test_source.jl")     # EGD::sw-singular-source: singular field (52)-(60), (B1), Fourier modes
    include("SpiralWaves/test_r2_generated.jl") # EGD::sw-mode-solver: route-R2 generated δû_φ, AMF (72a), torque density (72b)
    include("SpiralWaves/test_derived_equivalence.jl") # EGD::jd-equivalence-gate: route-R3 Julia derivation == the R1/R2 CAS oracles
    include("SpiralWaves/test_solver.jl")     # EGD::sw-mode-solver: implicit midpoint (64) + outgoing-wave BCs (66)-(68)
    include("SpiralWaves/test_boundary_response.jl") # EGD::fig10: WKB amplitude split, first-order BC, edge overrides
    include("SpiralWaves/test_diagnostics.jl") # EGD::sw-diagnostics: AMF (72a), torque (72b), (74), (75), m-cutoffs
    include("SpiralWaves/test_survey.jl")     # EGD::sw-survey-parallel: threaded m-sweep, configuration sweep
    include("SpiralWaves/test_image_synthesis.jl") # EGD::sw-image-synthesis: mode-sum wake reconstruction (Figs. 1-2)
    include("SpiralWaves/test_newtonian.jl")  # EGD::sw-newtonian-limit: (50) on a Newtonian background (Fig. 5 T_Newt)
    include("test_eccentric_legacy.jl")      # ECC::gr-disc-model: pre-context circular bit golden (96 Float64 values)
    include("test_eccentric_disc_model.jl")  # ECC::gr-disc-model: DiscModel units, hydrostatic balance, closure/evolution
    include("test_eccentric_context.jl")     # ECC::gr-frequency-context: arbitrary-ω complex contexts, e → 0 parity
    include("test_eccentric_grid.jl")        # ECC::gr-corotation-grid: physical grids, root/no-root/edge, Landau layer
    include("test_eccentric_li.jl")          # ECC::gr-li-runtime: OPTION A runtime, closure-bound cache and solutions
    include("test_eccentric_boundary.jl")    # ECC::gr-boundary-certificate: arbitrary-frequency radiation rows (66)-(68)
    include("test_eccentric_li_boundary.jl") # ECC::gr-li-boundary-certificate: OPTION A edge rows, thermal transport, GR transparency
    include("test_eccentric_source_bank.jl") # ECC::gr-azimuth-projector, ECC::gr-source-harmonics: eccentric source bank
    include("test_eccentric_modes.jl")       # ECC::gr-mode-solver: end-to-end adiabatic modes, e=0 gate, circular isolation
    include("test_eccentric_li_modes.jl")    # ECC::gr-li-mode-solver: end-to-end OPTION A modes, LI rows, closure law, isolation
    include("test_li_thermostat.jl")         # ECC::gr-li-thermostat-sensitivity: thermostat law, bath work, primary-dispatch separation
    include("test_eccentric_diagnostics.jl") # ECC::gr-{adiabatic,li}-balance-diagnostics, gr-{,li-}axisymmetric-power: direct work, LI Reynolds flux
    include("test_eccentric_convergence.jl") # ECC::gr-quadrature-convergence (PREREG_v2 audit), gr-radial-convergence tooling: norms, gates, tiny ladders
end

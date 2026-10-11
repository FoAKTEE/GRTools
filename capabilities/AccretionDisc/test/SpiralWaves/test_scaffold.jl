# SpiralWaves/test_scaffold.jl — EGD::sw-scaffold (arXiv:2608.17003 spiral-density-wave pipeline).
# Hygiene for the scaffold: the submodule loads and is reachable as `AccretionDisc.SpiralWaves`,
# ships NO physics stubs ahead of their tests (every exported name belongs to a LANDED node — the
# landed set is asserted below and grows with the DAG; planned-but-unlanded entry points must NOT
# exist), and the SparseArrays stdlib it will use for the sparse implicit solver is a DECLARED hard
# dep of the package (gate-5 style must_declare through TestKit.assert_dep_budget) and resolvable
# from the package environment. Physics nodes add their own test files as they land.
using Test
using AccretionDisc
using TestKit

@testset "SpiralWaves scaffold: load + dep budget" begin
    @test isdefined(AccretionDisc, :SpiralWaves)
    @test AccretionDisc.SpiralWaves isa Module
    @test :SpiralWaves in names(AccretionDisc)                   # re-exported module name
    SW = AccretionDisc.SpiralWaves
    # every exported symbol belongs to a LANDED node (no stubs): the landed export set is exact
    landed_background = [:DiscBackground, :rotation_profile, :Omega, :omega_d, :Omega_emri, :eps_profile,
                         :pressure, :gamma1, :dlogp_dr_over_enthalpy, :kappa2, :r_isco, :softening_length,
                         :F_J0, :F_J1, :disc_OmegaK, :kappar2, :csound, :eden0]     # EGD::sw-background
    landed_grid = [:SpiralGrid, :build_grid, :corotation_radius, :corotation_block, :smoothmin,
                   :wkb_wavenumber, :wkb_wavenumber2, :wkb_wavenumber_floored, :WKBRegion, :wkb_region,
                   :cell_count,                                                      # EGD::sw-grid
                   :wkb_local, :GRID_EXTENTS]                             # EGD::sw-grid-gr-refinement
    landed_master = [:MasterInputs, :master_inputs, :master_matrices, :master_matrices_newtonian,
                     :delta_u_phi_newtonian, :landau_shift, :dOmega_dr,
                     :delta_u_phi, :r2_context, :dv_dr, :rotation_derivs]        # EGD::sw-master-matrices (+R2)
    landed_source = [:SingularFieldCoeffs, :singular_field_coeffs, :worldline_velocity, :s0_squared,
                     :s0_squared_dr, :singular_field, :singular_field_mode, :source_modes, :source_modes!,
                     :raise_modes, :strip_halfwidth, :trapezoid_nodes]            # EGD::sw-singular-source
    landed_solver = [:ModeSolution, :solve_mode, :midpoint_system, :solve_midpoint, :midpoint_residual,
                     :eigen2, :outgoing_branch, :wkb_wavenumber_gr, :wkb_wavenumber2_gr,
                     :wkb_split, :wkb_amplitudes, :mode_amplitudes]                 # EGD::sw-mode-solver
    landed_diag = [:sqrt_minus_g, :amf, :amf_mode, :torque_density, :torque_density_mode,
                   :cumulative_torque, :total_torque, :lindblad_torque, :torque_units, :m_cut,
                   :balance_residual, :corotation_split, :corotation_window]         # EGD::sw-diagnostics
    landed_survey = [:ModeTorque, :TorqueSurvey, :torque_survey, :survey_configs, :mode_torques,
                     :mode_lindblad_torques, :compute_budget]                  # EGD::sw-survey-parallel
    landed_image = [:ModeInterp, :mode_interp, :synthesize, :synthesize_cartesian,
                    :polar_to_cartesian, :image_grid, :arm_phase,
                    :wake_modes]                                              # EGD::sw-image-synthesis
    landed_newtonian = [:newtonian_Omega_emri, :newtonian_csound, :newtonian_pressure, :newtonian_gamma1,
                        :newtonian_rotation, :newtonian_omega_d, :newtonian_corotation_radius,
                        :newtonian_potential_coeffs, :newtonian_potential_mode, :newtonian_potential_modes,
                        :newtonian_master_matrices, :newtonian_delta_u_phi, :newtonian_delta_sigma,
                        :newtonian_grid, :newtonian_solve_mode, :newtonian_amf_mode,
                        :newtonian_torque_density_mode, :newtonian_amf, :newtonian_torque_density,
                        :newtonian_balance_residual]                          # EGD::sw-newtonian-limit
    landed_boundary = [:radiation_impedance_row]                  # ECC::gr-boundary-certificate (C-66)
    @test Set(names(SW)) == Set(vcat(:SpiralWaves, landed_background, landed_grid, landed_master,
                                     landed_source, landed_solver, landed_diag, landed_survey,
                                     landed_image, landed_newtonian, landed_boundary))
    @test SW.PAPER == "arXiv:2608.17003"
    @test SW.sw_version() isa VersionNumber
    @test SW.sw_version() == v"0.12.0"
    # no placeholder physics: none of the UNLANDED planned node entry points exist yet
    for sym in (:emri_source, :torque_balance, :spiral_morphology, :pressureless_comparison)
        @test !isdefined(SW, sym)
    end
    # SparseArrays is a declared hard dep (and resolvable from the package env)
    @test Base.identify_package(AccretionDisc, "SparseArrays") !== nothing
    @test isdefined(SW, :SparseArrays) && SW.SparseArrays isa Module
    @test isdefined(SW, :Spacetime)                             # Kerr geometry ONLY via Spacetime
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "..", "Project.toml");
        must_declare=["SparseArrays", "Spacetime"],
        forbid=["SelfForce", "Kernels", "Symbolics", "JLD2", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Roots", "NLsolve", "SciMLBase"])
end

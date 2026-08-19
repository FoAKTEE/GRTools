"""
    AccretionDisc

Tier-6 EMRI-accretion-disc capability (design §4.10) — the TOP of the GrTools capability DAG.
The Lindblad-torque / spiral-density-wave disc front-end, fusing the legacy `DiscModes`
(relativistic disc-mode solver) and `DiscEMRIForce` (circular-equatorial Kerr thin-disc torque)
into two namespaced submodules:

  - [`Modes`](@ref) — the relativistic disc-mode solver: the radial epicyclic frequency `κ_r²`
    from the Spacetime geodesic SHEAR (physics, independent of the metric perturbation 𝒲), the
    height-averaged κ²-dispersion reduced operator `D = κ² − ω̃²` (Lin-Shu / Goldreich-Tremaine
    density-wave reduction), the `B_m` mode functions and the `𝔉_m` forcing source, the
    Lindblad/corotation resonance torque.
  - [`Force`](@ref) — the circular-equatorial Kerr thin-disc torque front-end, CONSUMING
    `MetricReconstruction.Lorenz.closure_kerr` / `reconstruct_kerr` (HDF5-backed captured data).
  - [`SpiralWaves`](@ref) — the arXiv:2608.17003 relativistic spiral-density-wave pipeline on
    the Page–Thorne equatorial line element (Lagrangian perturbation theory, two first-order
    master equations, corotation regularization, adaptive grid, torque balance); built node by
    node tests-first on `Spacetime.kerr_equatorial_potentials`.

Every primitive is type-generic `{T}` so BigFloat precision is inherited.
"""
module AccretionDisc

include("Modes/Modes.jl")
using .Modes
# Re-export the Modes background/resonance primitives at the package top level.
using .Modes: RelativisticDisc, disc_OmegaK, kappar2, delta_r, csound, eden0, risco_disc
using .Modes: lindblad_residual, corotation_residual, wronskian, lindblad_locations
using .Modes: A_m, B_m_corotation, B_m_wkb, B_m_full, corotation_exponent
using .Modes: dispersion, reduced_operator
using .Modes: lindblad_torque, forcing_Fm, laplace_b
using .Modes: SurfaceDensityProfile, PowerLawSigma
using .Modes: surface_density, surface_density_dp, surface_density_slope
using .Modes: hkr_epicyclic_poly, hkr_binding_poly, hkr_delta_poly
using .Modes: hkr_shear_argument, hkr_b0, hkr_b1, hkr_b2
using .Modes: hkr_lindblad_power_quadratic, hkr_lindblad_power_linear, hkr_lindblad_torque_rate
using .Modes: hkr_lindblad_torque_mode, hkr_lindblad_torque
using .Modes: hkr_mode_cutoff, hkr_mode_cutoff_powerlaw, HKR_M45_COEFFICIENT
using .Modes: hkr_torque_reversal_criterion, hkr_torque_reversal_radius
using .Modes: qinv_operator, q_operator
using .Modes: disc_OmegaK_z, alpha_z, Uvec_z, pressure_z, eden_z
using .Modes: qinv_operator_z, q_operator_z, heightavg
using .Modes: master_radial_coeffs
export RelativisticDisc, disc_OmegaK, kappar2, delta_r, csound, eden0, risco_disc
export lindblad_residual, corotation_residual, wronskian, lindblad_locations
export A_m, B_m_corotation, B_m_wkb, B_m_full, corotation_exponent
export dispersion, reduced_operator
export lindblad_torque, forcing_Fm, laplace_b
export SurfaceDensityProfile, PowerLawSigma
export surface_density, surface_density_dp, surface_density_slope
export hkr_epicyclic_poly, hkr_binding_poly, hkr_delta_poly
export hkr_shear_argument, hkr_b0, hkr_b1, hkr_b2
export hkr_lindblad_power_quadratic, hkr_lindblad_power_linear, hkr_lindblad_torque_rate
export hkr_lindblad_torque_mode, hkr_lindblad_torque
export hkr_mode_cutoff, hkr_mode_cutoff_powerlaw, HKR_M45_COEFFICIENT
export hkr_torque_reversal_criterion, hkr_torque_reversal_radius
export qinv_operator, q_operator
export disc_OmegaK_z, alpha_z, Uvec_z, pressure_z, eden_z
export qinv_operator_z, q_operator_z, heightavg
export master_radial_coeffs

# Force submodule (HDF5-backed circular-equatorial Kerr thin-disc torque front-end).
include("Force/Force.jl")
using .Force
using .Force: H1Mode, load_h1_mode, load_modes
using .Force: kerr_parity, tetmodes_from_bls, hmunu, hmunu_eq
using .Force: disc_a0, disc, force_mode, Vvec, source_F, source_F_grid
using .Force: FluidVars, reconstruct, Jadv
using .Force: ForcingInterp, build_forcing, forced_solve
using .Force: f8_torque_balance, f8_torque_balance_ilr, f8_flux_crosscheck, f5_source
using .Force: kerr_mode, kerr_modes
export H1Mode, load_h1_mode, load_modes
export kerr_parity, tetmodes_from_bls, hmunu, hmunu_eq
export disc_a0, disc, force_mode, Vvec, source_F, source_F_grid
export FluidVars, reconstruct, Jadv
export ForcingInterp, build_forcing, forced_solve
export f8_torque_balance, f8_torque_balance_ilr, f8_flux_crosscheck, f5_source
export kerr_mode, kerr_modes

# SpiralWaves submodule (arXiv:2608.17003 spiral-density-wave pipeline; exports are added per
# landed DAG node — the Modes generics disc_OmegaK/kappar2/csound/eden0 gain DiscBackground methods).
include("SpiralWaves/SpiralWaves.jl")
using .SpiralWaves
using .SpiralWaves: DiscBackground, rotation_profile, Omega, omega_d, Omega_emri, eps_profile, pressure
using .SpiralWaves: gamma1, dlogp_dr_over_enthalpy, kappa2, r_isco, softening_length, F_J0, F_J1
export DiscBackground, rotation_profile, Omega, omega_d, Omega_emri, eps_profile, pressure
export gamma1, dlogp_dr_over_enthalpy, kappa2, r_isco, softening_length, F_J0, F_J1
using .SpiralWaves: SpiralGrid, build_grid, corotation_radius, corotation_block, smoothmin
using .SpiralWaves: wkb_wavenumber, wkb_wavenumber2, wkb_wavenumber_floored, WKBRegion, wkb_region, cell_count
export SpiralGrid, build_grid, corotation_radius, corotation_block, smoothmin
export wkb_wavenumber, wkb_wavenumber2, wkb_wavenumber_floored, WKBRegion, wkb_region, cell_count
using .SpiralWaves: MasterInputs, master_inputs, master_matrices, master_matrices_newtonian, delta_u_phi_newtonian
using .SpiralWaves: landau_shift, dOmega_dr
export MasterInputs, master_inputs, master_matrices, master_matrices_newtonian, delta_u_phi_newtonian
export landau_shift, dOmega_dr
using .SpiralWaves: SingularFieldCoeffs, singular_field_coeffs, worldline_velocity, s0_squared, s0_squared_dr
using .SpiralWaves: singular_field, singular_field_mode, source_modes, source_modes!, raise_modes
using .SpiralWaves: strip_halfwidth, trapezoid_nodes
export SingularFieldCoeffs, singular_field_coeffs, worldline_velocity, s0_squared, s0_squared_dr
export singular_field, singular_field_mode, source_modes, source_modes!, raise_modes
export strip_halfwidth, trapezoid_nodes

export Modes, Force, SpiralWaves

end # module AccretionDisc

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
using .Modes: qinv_operator, q_operator
using .Modes: disc_OmegaK_z, alpha_z, Uvec_z, pressure_z, eden_z
using .Modes: qinv_operator_z, q_operator_z, heightavg
using .Modes: master_radial_coeffs
export RelativisticDisc, disc_OmegaK, kappar2, delta_r, csound, eden0, risco_disc
export lindblad_residual, corotation_residual, wronskian, lindblad_locations
export A_m, B_m_corotation, B_m_wkb, B_m_full, corotation_exponent
export dispersion, reduced_operator
export lindblad_torque, forcing_Fm, laplace_b
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

export Modes, Force

end # module AccretionDisc

"""
    AccretionDisc.Modes

The relativistic disc-mode solver (from legacy `DiscModes`). Leaves built tests-first, grounded
WITHOUT heavy live deps:

  - `kappar2`  — the radial epicyclic frequency κ_r² from the Spacetime geodesic SHEAR (closed
                 form, independent of the perturbation 𝒲).
  - the height-averaged κ²-dispersion reduced operator `D = κ² − ω̃²`.
  - `B_m`      — the disc-mode functions.
  - `𝔉_m`      — the forcing source.
  - the Lindblad / corotation resonance torque (Hirata GW-overlap kernel), and the CLOSED-FORM
    pressureless Kerr Lindblad torque of HegadeKR-Gammie-Yunes 2025 (arXiv:2510.03564) that
    arXiv:2608.17003 Fig. 9 compares its linear theory against.

Type-generic `{T}`; BigFloat-inherited.
"""
module Modes

include("Background.jl")
export RelativisticDisc, disc_OmegaK, kappar2, delta_r, csound, eden0, risco_disc

include("Velocity.jl")
export qinv_operator, q_operator

include("OffEquator.jl")
export disc_OmegaK_z, alpha_z, Uvec_z, pressure_z, eden_z
export qinv_operator_z, q_operator_z, heightavg

include("MasterPDE.jl")
export master_radial_coeffs

include("Resonances.jl")
export lindblad_residual, corotation_residual, wronskian, lindblad_locations, lindblad_locations_pressure

include("ReducedODE.jl")
export A_m, B_m_corotation, B_m_wkb, B_m_full, B_m_geometric, corotation_exponent
export dispersion, reduced_operator

include("Forcing.jl")
export lindblad_torque, forcing_Fm, laplace_b, laplace_b_deriv, c13_alpha
export epicyclic_impedance, lindblad_detuning, lindblad_detuning_deriv
export lindblad_strength, lindblad_strength_newtonian

include("ThinDisc.jl")
export NovikovThorneDisc, nt_isco, nt_factors, nt_surface_density, nt_scale_height, nt_aspect_ratio
export ss_surface_density, ss_aspect_ratio, pressure_cutoff, total_lindblad_torque, powerlaw_fit

include("NewtonianMaster.jl")
include("NewtonianPlanet.jl")
include("NewtonianSolver.jl")
include("LindbladCrossing.jl")
include("PhaseGradient.jl")
export phase_gradient_error, phase_error_diagnostic, phase_gradient_refine, mr19_phase_refinement
# MR19/MR20 Newtonian planet-forced disc: the model profiles stay UNEXPORTED (`Modes.Omega`,
# `Modes.Sigma`, … belong to this background only); the solver API and the softened Laplace family
# are exported.
export NewtonianDisc, master_equation_cool, loggrid, solve_m_mode, solve_modes, vars_from_enthalpy
export solve_lm_mode  # ECC::fr22-solver-mode-generalisation: frequency-aware (m,l) entry point
export PotentialCache, passive_mode, passive_sigma
export laplace_b_soft, laplace_b_soft_deriv, laplace_b_soft_2deriv

# DESIGN v3 §3 owners of the FR22/FR24 eccentric Newtonian planet.
include("EccentricPlanet.jl")   # FR22 (24)-(25) Kepler orbit; ZZ22 (4)-(9) wavelets
export EccentricOrbit, kepler_E, planet_orbit
export zz22_wavelet, zz22_track
export NewtonianPotentialBank, potential_bank, harmonic_cache
include("FR24Diagnostics.jl")  # FR24 per-harmonic flux, torque and power (unexported)
include("FR24Evolution.jl")    # Shared Keplerian orbital rates from modal diagnostics
include("FR22Potential.jl")     # ECC::fr22-potential-harmonics: factorised Gauss-E oracle (unexported)
include("EccentricTorques.jl")  # FR24 v1 (52)-(55) PL00/CN06 literature prescriptions
export pl00_tauL_inv, pl00_taue_inv, cn06_tauL_inv, cn06_taue_inv

include("PressurelessLindblad.jl")
export SurfaceDensityProfile, PowerLawSigma
export surface_density, surface_density_dp, surface_density_slope
export hkr_epicyclic_poly, hkr_binding_poly, hkr_delta_poly
export hkr_shear_argument, hkr_b0, hkr_b1, hkr_b2
export hkr_lindblad_power_quadratic, hkr_lindblad_power_linear, hkr_lindblad_torque_rate
export hkr_lindblad_torque_mode, hkr_lindblad_torque
export hkr_mode_cutoff, hkr_mode_cutoff_powerlaw, HKR_M45_COEFFICIENT
export hkr_torque_reversal_criterion, hkr_torque_reversal_radius

end # module Modes

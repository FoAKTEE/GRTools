"""
    AccretionDisc.Modes

The relativistic disc-mode solver (from legacy `DiscModes`). Leaves built tests-first, grounded
WITHOUT heavy live deps:

  - `kappar2`  — the radial epicyclic frequency κ_r² from the Spacetime geodesic SHEAR (closed
                 form, independent of the perturbation 𝒲).
  - the height-averaged κ²-dispersion reduced operator `D = κ² − ω̃²`.
  - `B_m`      — the disc-mode functions.
  - `𝔉_m`      — the forcing source.
  - the Lindblad / corotation resonance torque.

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
export lindblad_residual, corotation_residual, wronskian, lindblad_locations

include("ReducedODE.jl")
export A_m, B_m_corotation, B_m_wkb, B_m_full, corotation_exponent
export dispersion, reduced_operator

include("Forcing.jl")
export lindblad_torque, forcing_Fm, laplace_b

end # module Modes

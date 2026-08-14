"""
    MetricReconstruction.Radiation

Radiation-gauge CCK (Chrzanowski–Cohen–Kegeles) metric reconstruction: the Hertz potential
from Teukolsky amplitudes, the six tetrad H-components (IRG s=−2 / ORG s=+2), the
tetrad→coordinate projection, the regenerated gauge-invariant ψ0/ψ4, and the
linearized-Einstein vacuum residual. Type-generic `{T}` (BigFloat-inherited).
"""
module Radiation

using Harmonics
using Perturbation

include("Radiation/Operators.jl")         # SepField/RTJet jet algebra + Dn/Ln ladder operators
include("Radiation/Modes.jl")             # per-mode radial/angular jets + separation constant
include("Radiation/TetradComponents.jl")  # the six Hertz/CCK tetrad components (+ jets)
include("Radiation/WeylAmplitudes.jl")    # Teukolsky–Starobinsky constants + amplitude presets
include("Radiation/PlusMinus.jl")         # ±m tetrad-component value assembly
include("Radiation/Projection.jl")        # tetrad→coordinate projection (scalar + jet)
include("Radiation/Einstein.jl")          # linearized-Einstein vacuum residual (gate G-M4)
include("Radiation/WeylScalars.jl")       # regenerated ψ0/ψ4 from h (gate H10)
include("Radiation/Assembly.jl")          # parallel mode-sum + residual grid
include("Radiation/API.jl")               # user-facing single-mode surface

export SepField, sepfield, value, RTJet, rtconst, rtinv, sepfield_rtjet,
       Dn, Ln, mul_rjet, mul_ajet
export teuk_lambda, radial_jet, angular_jet, clear_mode_caches!
export HIRG_nn, HIRG_nm, HIRG_mm, HORG_ll, HORG_lm, HORG_mm
export HIRG_nn_jet, HIRG_nm_jet, HIRG_mm_jet, HORG_ll_jet, HORG_lm_jet, HORG_mm_jet
export calD, calC, Dhat, DhatPrime, weyl_amplitudes,
       radial_in, radial_out, radial_inPrime, radial_outPrime
export component_value
export tetrad_matrix, covariant_legs, h_coordinate, h_coordinate_jet
export lin_einstein_residual
export weyl_scalar, lin_riemann_lowered
export h_coord_modesum, einstein_grid
export metric_tetrad, h_tetrad, h_coord, reconstruct, psi

end # module Radiation

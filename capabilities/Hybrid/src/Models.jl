"""
    Hybrid.Models

WaSABI inspiral/amplitude model leaves: the per-model adiabatic forcing / amplitude tables
built over the [`Hybrid.Interpolation`](@ref) foundation. Routes the on-disk `.m`/`.wl` SF
tables through `Data`, the Kerr ISCO abscissa remap through the canonical `Spacetime.kerr_isco`,
and the first-law binding-energy channels through [`Hybrid.BindingEnergy`](@ref).
"""
module Models

using Spacetime: kerr_isco
using Data: read_wl_table, read_wl_table_complex, read_wl_data, read_wl_grid2d, read_wl_grid3d
using ..Interpolation
using ..BindingEnergy: E_FLx_deriv, E_dm, dEdΩ

# Re-export the interpolation foundation so the historical `Hybrid.Models.<interp>` surface holds.
export Interp1, interp1_order8, derivative, derivative2, interp_logexp,
    Interp2Hermite, interp2_hermite, partial, Interp3Lagrange, interp3_lagrange,
    kerr_grid_coords

include("Models/inspiral_0PASchwarz.jl")
export Forcing0PASchwarz, forcing_0PASchwarz, flux_1SF, forcing_F0, omega_dot, stop_r0,
    inspiral_rhs, inspiral_rhs!

include("Models/inspiral_1PAT1.jl")
export Forcing1PAT1, forcing_1PAT1, flux_1SF_deriv, flux_2SF, forcing_F1

include("Models/inspiral_1PAT1ea.jl")
export Forcing1PAT1ea, forcing_1PAT1ea, flux_LH, flux_sigma, flux_a, forcing_F1s1, forcing_F1s2,
    forcing_F1extra, forcing_F1dm, Model1PAT1ea, Model1PAT1R, model_1PAT1ea, model_1PAT1R,
    omega_crit, _spin_eom!

include("Models/inspiral_0PAKerr.jl")
export Forcing0PAKerr, forcing_0PAKerr, forcing_Fω, kerr_inspiral_rhs!

include("Models/amplitude_0PASchwarz.jl")
export Amp0PASchwarz, amplitudes_0PASchwarz, Z_I, Z_I_deriv, h0PA_amp, h0PA_ampΩ,
    mode_keys_0PASchwarz

include("Models/amplitude_1PAT1.jl")
export Amp1PAT1, amplitudes_1PAT1, Z2_I, h0PA_amp_deriv, h1PA_amp, h1PA_ampν, h1PA_ampΩ,
    mode_keys_1PAT1

include("Models/amplitude_1PAT1ea.jl")
export Amp1PAT1ea, amplitudes_1PAT1ea, h1PA_ampΩ_ea, mode_keys_1PAT1ea

include("Models/amplitude_0PAKerr.jl")
export Amp0PAKerr, amplitudes_0PAKerr, real_amp22, mode_keys_0PAKerr

end # module Models

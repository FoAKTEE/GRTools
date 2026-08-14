"""
    Hybrid.BindingEnergy

First-law binding-energy / flux closed forms over the redshift
[`Hybrid.Interpolation.Interp1`](@ref): the EFLx flux function, the secondary-mass correction
EΔm, and the per-frequency flux `dEdΩ`. Pure-algebra leaf consumed by the 1PA forcing.
"""
module BindingEnergy

using ..Interpolation: Interp1, derivative, derivative2

include("BindingEnergy/binding_energy.jl")
export E_FLx, E_FLx_deriv, E_dm, E_dm_deriv, dEdΩ

end # module BindingEnergy

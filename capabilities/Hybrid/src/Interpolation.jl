"""
    Hybrid.Interpolation

WL-faithful interpolation foundation (design §4.9): the moving-window order-8 `Interp1`, the
tensor-product global-Lagrange `Interp2Hermite` / `Interp3Lagrange`, their analytic
derivatives / mixed partials, and the Kerr grid-abscissa remap. This is the bottom of the
Hybrid DAG — `BindingEnergy` and the `Models` forcing/amplitude leaves both build on it, so it
lives in its own submodule to keep that layering acyclic.
"""
module Interpolation

using Spacetime: kerr_isco
using NumericsCore: Interp1WL, interp1_movingwindow, Interp2Hermite, interp2_hermite,
    partial, Interp3Lagrange, interp3_lagrange

include("Models/interpolation.jl")

export Interp1, interp1_order8, derivative, derivative2, interp_logexp,
    Interp2Hermite, interp2_hermite, partial, Interp3Lagrange, interp3_lagrange,
    kerr_grid_coords

end # module Interpolation

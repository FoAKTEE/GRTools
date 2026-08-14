"""
    MetricReconstruction.Lorenz

Dolan-2024 Lorenz-gauge metric reconstruction for a circular Kerr (a≠0) orbit. Owns the
angular/radial closure machinery, component assembly, App-H jump conditions, and the
Lorenz-gauge `weyl_scalar`/`to_BLS` variants. Also hosts [`Lorenz.LorenzRadial`](@ref),
the canonical a=0 Lorenz-gauge radial field-equation RHS consumed downstream.
Type-generic `{T}` (BigFloat-inherited).
"""
module Lorenz

# includes in dependency order; one flat namespace (matching the legacy module layout).

include("Lorenz/Operators.jl")            # geometry/tetrad/Chandrasekhar scalars + jets + CG mixing
export geom, lpm, mpm, ginv_tetrad, kerr_metric_inv_eq, Kfn, Qfn, Dcoef
export jet_ops, jinv, constjet, cg, frak_s, frak_s2, frak_c, frak_c2

include("Lorenz/Background.jl")           # circular-equatorial Kerr orbit kinematics
export OrbitBackground, orbit_background

include("Lorenz/Jumps.jl")                # App-H spin-general source jumps
export jumps_pm, hj1, jumps_eigenvalue

include("Lorenz/Closure.jl")              # spin-general banded jump closure (+ a=0 path)
export closure_kerr, closure_a0

include("Lorenz/KerrRadial.jl")           # a≠0 full radial reconstruction (P_{-2},P_{-1},h,κ)
export reconstruct_kerr, HomFieldK, ScaledHomFieldK, KappaFieldK
export valK, derK, der2K, κvalK, κderK, radial_gates_kerr, jump_data_kerr

include("Lorenz/SchwarzschildRadial.jl")  # a=0 radial reconstruction — surface + a→0 oracle
export reconstruct, HomField, KappaField, val, der, der2, κval, κder, jump_data_radial, radial_gates

include("Lorenz/SchwarzschildComponents.jl")  # a=0 ten tetrad modes → BLS assembly + a→0 oracle
export Recon, build_recon, tetrad_modes, to_BLS

include("Lorenz/KerrSMatrices.jl")        # a≠0 closure S-matrices + radial chains + lplp/mpmp
export smatrices_kerr, radial_chain_kerr, tetrad_lplp_kerr, tetrad_mpmp_kerr

include("Lorenz/KerrComponents.jl")       # a≠0 ten-component tetrad surface + BLS mapping
export tetrad_modes_kerr, to_BLS_kerr, tetrad_ρlpmp_kerr, tetrad_ΣΔlplm_kerr

include("Lorenz/LorenzWeyl.jl")           # gauge-invariant ψ4 cross-check
export weyl_psi4_lorenz, weyl_psi4_lorenz_amp

include("Lorenz/LorenzRadial.jl")         # canonical a=0 Lorenz radial ODE RHS (nested module)

end # module Lorenz

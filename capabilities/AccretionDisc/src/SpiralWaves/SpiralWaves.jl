"""
    AccretionDisc.SpiralWaves

The relativistic spiral-density-wave pipeline of arXiv:2608.17003 (Hegade K R, Yunes et al.:
"Relativistic spiral density waves and disk-EMRI torques" — linear Lagrangian perturbation theory
of a razor-thin, stationary, axisymmetric, purely rotating disc on the Page–Thorne equatorial
line element ds² = −e^{2ν}dt² + e^{2ψ}(dφ − ω̃ dt)² + e^{2μ}dr² + dz², driven by the simplified
gravitational potential of a circular equatorial EMRI).

Design (one file per DAG node, each landed tests-first; the include list below is the planned
dependency ORDER and is uncommented as each node lands — no stub is shipped ahead of its tests):

    # include("background.jl")   # EGD::sw-background LANDED (v0.2.0): disc equilibrium on (ν,ψ,μ,ω̃) —
    #                            #      Ω(r), v, u^t, u_t, u_φ from the exact (11), power laws (12), Ω• (13),
    #                            #      Γ1 (22), κ², r_ISCO, F_{J,0/1} (76); over Spacetime.kerr_equatorial_potentials
    # include("master.jl")       # EGD::sw-master-matrices LANDED (v0.4.0): the 2×2 first-order master
    #                            #      system dy/dr = C y + S, y = (δÛ, −iδû^r) (47) — relativistic C, S
    #                            #      GENERATED from the supplementary notebook (master_generated.jl, route
    #                            #      R1), Newtonian (50a-g), Landau shift (62-63), dΩ/dr; over background + grid
    # include("source.jl")       # EGD::sw-singular-source LANDED (v0.5.0): the singular field h_μν (52) of
    #                            #      the circular EMRI worldline (53)-(55), s0² (58)/(B1), softening (59),
    #                            #      cosine Fourier modes (60) + analytic ∂_r by the periodic trapezoid
    #                            #      rule (node count from the analytic strip), covariant + contravariant,
    #                            #      grid-vectorised `source_modes`; over background
    # include("grid.jl")         # EGD::sw-grid LANDED (v0.3.0): adaptive radial grid (App. C) — r_c (61),
    #                            #      sinh corotation block (C4-5), WKB k (C2,C6), Δr_WKB/smoothMin (C7-10),
    #                            #      cubic-spline node placement, per-mode domain (C1-3); NumericsCore
    # include("solver.jl")       # EGD: implicit two-first-order-equation integrator, sparse
    #                            #      (SparseArrays) banded assembly, boundary conditions
    # include("diagnostics.jl")  # EGD: energy / angular-momentum balance laws, torque per m,
    #                            #      spiral morphology diagnostics
    # include("survey.jl")       # EGD: spin / radius / aspect-ratio / profile surveys, m-cutoff
    #                            #      ~3/h convergence, pressureless-disc comparison

Lives ABOVE `Modes` (the legacy Lin–Shu reduced-operator disc-mode solver) and `Force` (the
Lorenz-closure torque front-end) and reaches Kerr geometry ONLY through `Spacetime`
(`kerr_equatorial_potentials`, `kerr_equatorial_potentials_dr`, `kerr_circular_dtdtau`,
`kerr_isco`, `horizons`). Every primitive is type-generic `{T}`; BigFloat precision is inherited.
Exports are added per landed node; the scaffold exports nothing.
"""
module SpiralWaves

using Spacetime
using SparseArrays
using ..Modes

"The arXiv identifier of the paper this submodule reproduces."
const PAPER = "arXiv:2608.17003"

"""
    sw_version() -> VersionNumber

Version of the `SpiralWaves` submodule API (bumped as DAG nodes land; `v"0.1.0"` = scaffold,
`v"0.2.0"` = background, `v"0.3.0"` = grid, `v"0.4.0"` = master matrices, `v"0.5.0"` = singular source,
`v"0.6.0"` = mode solver).
"""
sw_version() = v"0.6.0"

# ---- EGD::sw-background — disc equilibrium (Sec. II.A) -------------------------------------------
include("background.jl")
export DiscBackground, rotation_profile, Omega, omega_d, Omega_emri, eps_profile, pressure, gamma1
export dlogp_dr_over_enthalpy, kappa2, r_isco, softening_length, F_J0, F_J1
# the EXISTING Modes generics gain DiscBackground methods (one name = one concept); re-exported here
export disc_OmegaK, kappar2, csound, eden0

# ---- EGD::sw-grid — adaptive corotation/WKB radial grid (App. C) ---------------------------------
include("grid.jl")
export SpiralGrid, build_grid, corotation_radius, corotation_block, smoothmin
export wkb_wavenumber, wkb_wavenumber2, wkb_wavenumber_floored, WKBRegion, wkb_region, cell_count

# ---- EGD::sw-master-matrices — relativistic master system (Sec. II.C) -----------------------------
# (included AFTER grid.jl: the Landau shift uses grid.jl's `corotation_radius`; master_generated.jl is
#  the notebook transcription emitted by the emri-disk generator codes/cas/R1_to_julia.wls)
include("master.jl")
export MasterInputs, master_inputs, master_matrices, master_matrices_newtonian, delta_u_phi_newtonian
export delta_u_phi, r2_context, dv_dr, rotation_derivs
export landau_shift, dOmega_dr

# ---- EGD::sw-singular-source — the EMRI singular field and its Fourier modes (Sec. II.D, App. B) -----
# (depends on background.jl only: softening_length, Omega_emri; the modes feed master_inputs(; hmodes))
include("source.jl")
export SingularFieldCoeffs, singular_field_coeffs, worldline_velocity, s0_squared, s0_squared_dr
export singular_field, singular_field_mode, source_modes, source_modes!, raise_modes
export strip_halfwidth, trapezoid_nodes

# ---- EGD::sw-mode-solver — implicit-midpoint global solve + outgoing-wave BCs (Sec. II.E) ---------
# (needs master.jl for C, S and delta_u_phi, source.jl for the modes at the midpoints, grid.jl for r)
include("solver.jl")
export ModeSolution, solve_mode, midpoint_system, solve_midpoint, midpoint_residual
export eigen2, outgoing_branch, wkb_wavenumber_gr, wkb_wavenumber2_gr

end # module SpiralWaves

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
    # include("master.jl")       # EGD: the two first-order master equations for (δh, δu^r) per
    #                            #      azimuthal m (main.tex Sec. master-eqs) + WKB/Lindblad forms
    # include("source.jl")       # EGD: the EMRI potential Φ_m source — s₀² approximation
    #                            #      (appendix.tex Sec. ssq-expr), softening ε, worldline v^μ
    # include("grid.jl")         # EGD: adaptive radial grid — corotation layer + WKB oscillation
    #                            #      resolution, Landau-type corotation regularization
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
`v"0.2.0"` = background).
"""
sw_version() = v"0.2.0"

# ---- EGD::sw-background — disc equilibrium (Sec. II.A) -------------------------------------------
include("background.jl")
export DiscBackground, rotation_profile, Omega, omega_d, Omega_emri, eps_profile, pressure, gamma1
export dlogp_dr_over_enthalpy, kappa2, r_isco, softening_length, F_J0, F_J1
# the EXISTING Modes generics gain DiscBackground methods (one name = one concept); re-exported here
export disc_OmegaK, kappar2, csound, eden0

end # module SpiralWaves

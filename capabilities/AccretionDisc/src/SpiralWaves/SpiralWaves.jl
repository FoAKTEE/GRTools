"""
    AccretionDisc.SpiralWaves

The relativistic spiral-density-wave pipeline of arXiv:2608.17003 (Hegade K R, Yunes et al.:
"Relativistic spiral density waves and disk-EMRI torques" — linear Lagrangian perturbation theory
of a razor-thin, stationary, axisymmetric, purely rotating disc on the Page–Thorne equatorial
line element ds² = −e^{2ν}dt² + e^{2ψ}(dφ − ω̃ dt)² + e^{2μ}dr² + dz², driven by the simplified
gravitational potential of a circular equatorial EMRI).

Design (one file per DAG node, each landed tests-first; the include list below is the planned
dependency ORDER and is uncommented as each node lands — no stub is shipped ahead of its tests):

    # include("background.jl")   # EGD: disc equilibrium on (ν,ψ,μ,ω̃) — Ω(r), v, u^t, enthalpy,
    #                            #      sound speed / aspect ratio, κ², the dpdr equilibrium
    #                            #      (main.tex:172-200); consumes Spacetime.kerr_equatorial_potentials
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

"The arXiv identifier of the paper this submodule reproduces."
const PAPER = "arXiv:2608.17003"

"""
    sw_version() -> VersionNumber

Version of the `SpiralWaves` submodule API (bumped as DAG nodes land; `v"0.1.0"` = scaffold).
"""
sw_version() = v"0.1.0"

end # module SpiralWaves

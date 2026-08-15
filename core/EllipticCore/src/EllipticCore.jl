"""
    EllipticCore (GrToolsElliptic)

Generic complex-parametric 2D hp-DG elliptic solver — Tier-1 core stratum, **zero
physics** (elliptic-ssf DESIGN.md §1). Any capability with a 2D elliptic problem
(the scalar-self-force m-mode operator today; Lorenz-gauge gravitational fields,
modified-gravity elliptic problems later) consumes this package; physics coefficients
enter only through caller-supplied `fluxes!`/`sources!` closures and boundary hooks.

Owns (DESIGN.md §2.1, one DAG node = one commit per module):
  - `refelement.jl` — tensor-product Legendre-**Gauss** reference elements over
    `NumericsCore.gl_nodes_weights` (nodes/weights, barycentric basis, D, face interp)
  - `mesh2d.jl` — multi-block `Mesh2D` (linear + Inverse 1/r maps), refinement-tree
    `ElementId`, 2:1 balance, per-element Jacobians
  - `dg_matrices.jl` — per-element M, MD_i, ML, MLD^T with cached 1D constituents
  - `dg_operator.jl` — matrix-free strong-form first-order-flux DG residual
  - `dg_penalty.jl` — internal-penalty numerical flux (σ per Vu:2024, C=1.5 default)
  - `mortars.jl` — 1D h/p mortar projection + mass-conservative restriction
  - `boundary_hooks.jl` — `ModifyBoundaryData` hooks + external-BC ghost seam
  - `matrixfree.jl` — full `apply!` composition + `build_matrix` column probes
  - `gmres.jl` — complex-generic GMRES (Arnoldi/MGS + complex Givens)
  - `schwarz.jl` — element-centered additive Schwarz with cached (complex-shifted)
    subdomain factorizations
  - `multigrid.jl` — geometric multigrid V-cycle with L2 transfers from mortars
  - `bottom_solve.jl` — coarsest-level cached sparse LU (UMFPACK), dense-LU fallback
  - `solver_stack.jl` — right-preconditioned GMRES composition API `solve!`
  - `amr.jl` — h/p mesh transformations with 2:1 enforcement + solution carry-over

Dep budget (binding, enforced by `test/hygiene.jl`): {LinearAlgebra, SparseArrays,
NumericsCore}. No external Krylov/AMG package — GMRES/multigrid/Schwarz are in-tree
per the GrTools vendoring precedent; no GPU, no AD, no IO, no GrTools capability.

Conventions (DESIGN.md §4): every struct is parametric in the scalar type `T<:Number`
with geometry in `RT = real(T)` — ComplexF64 is first-class, BigFloat-compatible where
cheap (precision inherited from NumericsCore; the UMFPACK bottom solve is
Float64/ComplexF64-only with a dense generic-LU fallback). Field storage is
`(npoints_total, nfields)` with `nfields = 1` for the scalar problem — the 10-DOF
Lorenz-gauge future is a parameter, not a fork. Complex discipline: inner products
conjugate on the left (`LinearAlgebra.dot`); GMRES Givens keeps the cosine real and
the sine complex.
"""
module EllipticCore

using LinearAlgebra
using SparseArrays
using NumericsCore

# Modules land TDD (failing tests first, grounded per GLOBAL_DAG oracle gates) in the
# scheduled waves. Comment-only headers now; each module's exports arrive with its node.
include("refelement.jl")     # DAG: mesh-basis-quadrature
include("mesh2d.jl")         # DAG: mesh2d-blocks
include("dg_matrices.jl")    # DAG: dg-operator-matrices
include("dg_operator.jl")    # DAG: dg-strong-residuals
include("dg_penalty.jl")     # DAG: dg-penalty-flux
include("mortars.jl")        # DAG: dg-mortars
include("boundary_hooks.jl") # DAG: boundary-data-hooks
include("matrixfree.jl")     # DAG: matrix-free-operator
include("gmres.jl")          # DAG: complex-gmres
include("schwarz.jl")        # DAG: schwarz-smoother, schwarz-complex-shift
include("multigrid.jl")      # DAG: multigrid-hierarchy
include("bottom_solve.jl")   # DAG: multigrid-bottom-solver
include("solver_stack.jl")   # DAG: solver-stack-core
include("amr.jl")            # DAG: amr-h-refinement, amr-p-refinement

end # module EllipticCore

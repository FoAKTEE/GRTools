# EllipticCore — hp-DG elliptic solver core

Generic complex-parametric 2D hp-DG elliptic solver — a Tier-1 core stratum with
**zero physics**. Any capability with a 2D elliptic problem (the scalar-self-force
m-mode operator today; Lorenz-gauge gravitational fields and modified-gravity
elliptic problems later) consumes this package; physics coefficients enter only
through caller-supplied `fluxes!`/`sources!` closures and boundary hooks.

## What this package owns

| module | one line |
|---|---|
| `refelement.jl` | tensor-product Legendre-Gauss reference elements over `NumericsCore.gl_nodes_weights` (nodes/weights, barycentric basis, D, face interpolation) |
| `mesh2d.jl` | multi-block `Mesh2D` (linear + inverse-radial 1/r maps), refinement-tree `ElementId`, 2:1 balance, per-element Jacobians |
| `dg_matrices.jl` | per-element M, MD_i, ML, MLDᵀ with cached 1D constituents |
| `dg_operator.jl` | matrix-free strong-form first-order-flux DG residual |
| `dg_penalty.jl` | internal-penalty numerical flux (σ per Vu 2024, C = 1.5 default) |
| `mortars.jl` | 1D h/p mortar projection + mass-conservative restriction |
| `boundary_hooks.jl` | `ModifyBoundaryData` hooks + external-BC ghost seam |
| `matrixfree.jl` | full `apply!` composition + `build_matrix` column probes |
| `gmres.jl` | complex-generic GMRES (Arnoldi/MGS + complex Givens) |
| `schwarz.jl` | element-centered additive Schwarz with cached subdomain LUs |
| `schwarz_shift.jl` | complex-shift policy `helmholtz_shift(m, Ω; slicing)` (`(mΩ)²` in t-slicing, 0 under hyperboloidal/null) + the `schwarz_preconditioner` right-preconditioner adapter for `gmres` |
| `multigrid.jl` | geometric multigrid V-cycle with L2 transfers from mortars |
| `bottom_solve.jl` | coarsest-level cached sparse LU (UMFPACK), dense-LU fallback |
| `solver_stack.jl` | right-preconditioned GMRES composition API `solve!` |
| `amr.jl` | h/p mesh transformations with 2:1 enforcement + solution carry-over |

## Dependencies

Binding dependency budget, enforced by `test/hygiene.jl`:
`{LinearAlgebra, SparseArrays, NumericsCore}`. No external Krylov/AMG package —
GMRES, multigrid, and Schwarz are in-tree per the GrTools vendoring precedent;
no GPU, no AD, no IO, no GrTools capability package.

## Conventions

Every struct is parametric in the scalar type `T <: Number` with geometry in
`RT = real(T)` — ComplexF64 is first-class, BigFloat-compatible where cheap
(precision inherited from NumericsCore; the UMFPACK bottom solve is
Float64/ComplexF64-only with a dense generic-LU fallback). Field storage is
`(npoints_total, nfields)` with `nfields = 1` for the scalar problem — the
10-DOF Lorenz-gauge future is a parameter, not a fork. Complex discipline:
inner products conjugate on the left (`LinearAlgebra.dot`); GMRES Givens keeps
the cosine real and the sine complex.

## Testing

`Pkg.test` standalone. Grounding oracles: exact polynomial/manufactured
solutions, sparse-matrix cross-checks of the matrix-free operator, and
convergence-order gates; the dependency budget is itself a test
(`test/hygiene.jl`).

# EllipticCore

Generic complex-parametric 2D hp-DG elliptic solver — a Tier-1 core stratum with
**zero physics**. Any capability with a 2D elliptic problem (the scalar-self-force
m-mode operator today; Lorenz-gauge gravitational fields and modified-gravity
elliptic problems later) consumes this package; physics coefficients enter only
through caller-supplied `fluxes!`/`sources!` closures and boundary hooks.

## Capability overview

The solver stack, bottom to top:

- **Reference elements** — tensor-product Legendre-Gauss elements over
  `NumericsCore.gl_nodes_weights`: nodes/weights, barycentric basis,
  differentiation matrix, face interpolation ([`RefElement`](@ref),
  [`refelement`](@ref)).
- **Meshes** — multi-block [`Mesh2D`](@ref) with linear and inverse-radial
  (1/r) coordinate maps, refinement-tree `ElementId`s, 2:1 balance, and
  per-element Jacobians.
- **DG operator** — per-element mass/stiffness/lifting matrices with cached 1D
  constituents; a matrix-free strong-form first-order-flux residual; the
  internal-penalty numerical flux (σ per Vu 2024, C = 1.5 default); 1D h/p
  mortar projection with mass-conservative restriction; boundary-data hooks
  and the external-BC ghost seam; the full matrix-free `apply!` composition
  plus `build_matrix` column probes.
- **Linear solvers** — complex-generic GMRES (Arnoldi/MGS with complex
  Givens); element-centered additive Schwarz smoothing with cached
  (complex-shifted) subdomain factorizations; geometric multigrid V-cycles
  with L2 transfers from the mortars; a coarsest-level cached sparse LU
  (UMFPACK, dense-LU fallback); and the right-preconditioned composition API
  `solve!`.
- **AMR** — h/p mesh transformations with 2:1 enforcement and solution
  carry-over.

## Conventions

Every struct is parametric in the scalar type `T <: Number` with geometry in
`RT = real(T)` — ComplexF64 is first-class, BigFloat-compatible where cheap
(precision inherited from NumericsCore; the UMFPACK bottom solve is
Float64/ComplexF64-only with a dense generic-LU fallback). Field storage is
`(npoints_total, nfields)` with `nfields = 1` for the scalar problem — the
10-DOF Lorenz-gauge future is a parameter, not a fork. Complex discipline:
inner products conjugate on the left (`LinearAlgebra.dot`); GMRES Givens keeps
the cosine real and the sine complex.

## Dependencies and poison-axis declaration

Binding dependency budget, enforced by `test/hygiene.jl`:
`{LinearAlgebra, SparseArrays, NumericsCore}`. No external Krylov/AMG package —
GMRES, multigrid, and Schwarz are in-tree per the GrTools vendoring precedent;
no GPU, no AD, no IO, no GrTools capability package.

## Public API

```@autodocs
Modules = [EllipticCore]
Order = [:module, :type, :constant, :function]
```

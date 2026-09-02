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
  (complex-shifted) subdomain factorizations; the complex-shift **policy**
  [`helmholtz_shift`](@ref) (`s = (mΩ)²` in t-slicing, `0` under
  hyperboloidal/null slicing) with its mass-weighted application route and the
  [`schwarz_preconditioner`](@ref) adapter that feeds a Schwarz smoother to
  `gmres` as a right preconditioner. (The right-preconditioned composition API
  `solve!` arrives with the solver-stack node.)
- **Geometric multigrid** — [`MultigridHierarchy`](@ref) h-coarsens a
  [`Mesh2D`](@ref) down to the irreducible block structure, one
  [`coarsen_mesh`](@ref) step at a time (per-dimension level capping: tiling,
  2:1 balance and cross-block compatibility are preserved, and polynomial
  orders p-coarsen by the minimum rule), rediscretizing the caller's
  coefficient closures on every level through an operator factory and giving
  each level its own additive-Schwarz smoother with a private LU cache. The
  transfers are the sealed mortar projections: [`mg_prolong!`](@ref) is the
  nodal Lagrange interpolation onto the children, [`mg_restrict!`](@ref) is
  either the Galerkin restriction `Pᵀ` for the **massive** residual (the
  V-cycle's transfer, exact for any block map) or the mass-conservative nodal
  L2 restriction for fields, and [`mg_transfer_check`](@ref) reports the three
  transfer diagnostics (prolongation exactness, `R∘P = I`, mass-weighted
  adjointness — all ≤ 1.6e-15 at Float64, ≤ 1.1e-76 at BigFloat(256)).
  [`MultigridV`](@ref) is the V-cycle preconditioner itself — pre-smoothing on,
  no post-smoothing at the bottom, one V-cycle per application, two Schwarz
  passes per level (the CRP yaml configuration) — and is directly the `Minv!`
  callable of [`gmres`](@ref); [`vcycle!`](@ref) applies it to an existing
  iterate. Its coarsest-level `bottom` seam defaults to
  [`dense_bottom_solver`](@ref) (one `lu` of `build_matrix`, factored once);
  the cached sparse UMFPACK factorization lands with the bottom-solver node.
  Measured contraction 0.171 at 3, 4 and 5 levels (level-independent), and 7
  GMRES iterations where the Schwarz smoother alone needs 18-32.
- **AMR** — mesh transformations for the CRP hp-AMR strategy, with
  caller-supplied refinement criteria: [`h_split`](@ref) halves the flagged
  elements along their flagged dimensions, [`balance_two_to_one`](@ref)
  restores the two-to-one balance by splitting the offending elements (within
  a block and transversally across a `BlockConnection`),
  [`h_refine`](@ref) composes both into one pure `Mesh2D` → `Mesh2D` step,
  [`p_refine`](@ref) raises the polynomial order of every element *except* the
  flagged ones (whose order is held — the CRP h-refines the puncture elements
  instead of p-refining them, so one flag set drives both halves),
  [`hp_refine`](@ref) is the full AMR step from that one criterion, and
  [`transfer_solution`](@ref) carries the previous level's nodal solution over
  by Lagrange interpolation (exact on every polynomial the source element
  represents, across an h step, a p step, or the two combined).

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

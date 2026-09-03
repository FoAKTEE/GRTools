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
- **Per-mortar boundary data** — a [`MeshOperator`](@ref) `hooks` entry binds
  either one [`BoundaryDataHook`](@ref) for a whole face or **one hook per
  mortar**, in the receiver's mortar order (spectre keys `ModifyBoundaryData`
  by `mortar_id`; the CRP's 2:1 worldtube faces carry different
  ``\Psi^{\mathcal P}`` on their two halves). A single hook is replicated, so
  the per-face spelling stays bitwise unchanged; a wrong-length collection is
  an `ArgumentError` naming the face. [`mortar_points`](@ref) publishes the
  mortar quadrature geometry the operator itself places — the face abscissa and
  the per-mortar transverse points — so a caller building per-mortar data lands
  on the operator's own grid; [`Mesh2D`](@ref)'s `faces` field and
  [`MeshInterface`](@ref) are the public-read neighbour/mortar query behind it.
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
  [`coarsen_mesh`](@ref) step at a time: each step merges every **complete
  sibling group** by one refinement level (largest dimension mask first, each
  group consumed as a whole) and restores the two-to-one balance with the
  sealed anisotropic [`balance_two_to_one`](@ref), so the rule is valid on the
  anisotropically h-refined meshes the CRP's AMR policy produces (a 300-mesh
  census over AMR-reachable anisotropic meshes builds 300/300 valid, balanced,
  nested, exactly-tiling hierarchies) and reduces to spectre's uniform
  decrement on a uniform mesh (64 → 16 → 4 → 1). Polynomial orders p-coarsen by
  the minimum rule. Every level rediscretizes the caller's coefficient closures
  through an operator factory and carries its own additive-Schwarz smoother;
  `caches` shares the subdomain LU factorizations across hierarchy rebuilds and
  AMR steps (rebuilding on identical meshes costs zero new factorizations), and
  `shrink` trades levels for coarse-space quality by composing the transfers of
  steps that do not shrink the problem enough.
  The transfers are the sealed mortar projections:
  [`mg_prolong!`](@ref) is the nodal Lagrange interpolation onto the children;
  [`mg_restrict_residual!`](@ref) is the Galerkin restriction `Pᵀ` of the
  **massive** residual (the V-cycle's transfer, exact for any block map);
  [`mg_restrict_field!`](@ref) is the mass-conservative nodal L2 restriction for
  fields; and [`mg_transfer_check`](@ref) reports the three transfer
  diagnostics — prolongation exactness, `R∘P = I`, and adjointness in the
  **logical (dyadic-measure) Gauss mass** (the physical-L2 form of that identity
  is affine-only, and its O(1) failure on the 1/r compactified block is pinned
  by the test suite) — all ≤ 1.6e-15 at Float64 and ≤ 1.1e-76 at BigFloat(256).
  [`MultigridV`](@ref) is the V-cycle preconditioner itself — pre-smoothing on,
  no post-smoothing at the bottom, one V-cycle per application, two Schwarz
  passes per level (the CRP yaml configuration) — and is directly the `Minv!`
  callable of [`gmres`](@ref); [`vcycle!`](@ref) applies it to an existing
  iterate. Its coarsest-level `bottom` seam defaults to
  [`dense_bottom_solver`](@ref) (one `lu` of `build_matrix`, factored once).
  Measured contraction 0.171 at 3, 4 and 5 levels (level-independent) and on
  the graded puncture cascade, and 7 GMRES iterations where the Schwarz
  smoother alone needs 18-32.
- **Bottom solver** — [`sparse_bottom_solver`](@ref) is the production filling
  of that seam and the paper's own route (main.tex:674): the coarsest matrix
  assembled by the sealed column probes, converted to `SparseArrays.sparse`,
  factored with the stdlib sparse `lu` (UMFPACK; no new dependency), and
  applied allocation-free (measured 0 bytes per solve, relative residual
  ≤ 6.3e-16, agreeing with the dense route to 1.1e-15 and keeping the sealed
  GMRES iteration counts). A caller-owned [`BottomSolverCache`](@ref) keeps the
  factorization across every V-cycle **and** across AMR rebuilds whose coarsest
  grid is unchanged — the CRP claim that the bottom LU survives AMR iterations,
  measured as one build and seven hits over an eight-step AMR run. Its key is
  the scalar type, the mesh identity, the penalty, the caller's `epoch`, and an
  operator fingerprint `hash(A·v)` costing one extra operator apply (a cache
  hit is 1.5 % of a build), so a coefficient, penalty or boundary-condition
  change rebuilds even at an unchanged epoch. `BigFloat`, `Complex{BigFloat}`
  and the 32-bit floats — everything UMFPACK does not take exactly — fall back
  to the generic dense LU in the caller's own precision (4.9e-77 at
  BigFloat(256)).
- **AMR** — mesh transformations for the CRP hp-AMR strategy, with
  caller-supplied refinement criteria: [`h_split`](@ref) halves the flagged
  elements along their flagged dimensions, [`balance_two_to_one`](@ref)
  restores the two-to-one balance by splitting the offending elements (within
  a block and transversally across a `BlockConnection`),
  [`h_refine`](@ref) composes both into one pure `Mesh2D` → `Mesh2D` step,
  [`p_refine`](@ref)`(mesh; hold = …)` raises the polynomial order of every
  element *except* those named by `hold`, whose order is held — the CRP
  h-refines the puncture elements instead of p-refining them, so one criterion
  drives both halves of a step; `hold` is a keyword, and named for what it
  does, because its polarity is the opposite of `h_refine`'s "split these".
  [`hp_refine`](@ref) is the full AMR step from that one criterion — literally
  `h_refine(p_refine(mesh; hold = flags), flags)` — and
  [`transfer_solution`](@ref) carries the previous level's nodal solution over
  by Lagrange interpolation (exact on every polynomial the source element
  represents, across an h step, a p step, or the two combined; an order
  *decrease* is refused, since nodal interpolation onto fewer points is not the
  L2 projection).

  The 2:1 closure follows the CRP's own solver input
  (`anc/ScalarSelfForce.yaml:145`, `Isotropy: Anisotropic`): only the offending
  dimension of a violating pair is split, and the dimensions marked in one sweep
  are applied together, matching spectre's `amr::update_amr_decision` fixpoint
  over pre-split decisions followed by a single domain adjustment. That is not
  the *minimal* closure — a one-split-at-a-time alternative yields up to 16.7 %
  fewer elements — but it is the algorithm the paper ran.

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

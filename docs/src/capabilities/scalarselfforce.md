# ScalarSelfForce

The m-mode **elliptic** Kerr scalar self-force capability: the GrTools
reproduction of arXiv:2606.04998 (Vu et al., *Self-force calculations with
numerical relativity methods*). It owns the physics of the elliptic
formulation and consumes the physics-free hp-DG solver stack from
[EllipticCore](../core/ellipticcore.md).

## Capability overview

- **Orbits** — `CircularOrbit{T}` thin adapters over Spacetime: E, L, Ω, uᵗ,
  u^φ, r_ISCO, tortoise r_*, and the horizon-regular azimuthal offset φ̄.
- **MModeOperator** — the operator pipeline: mode ansatz → t-sliced
  ``\tilde\triangle_m`` → hyperboloidal slicing → vtu junctions →
  horizon-penetrating first-order flux form (``F^i``, ``\beta_m``,
  ``\gamma^i_m``) → equatorial-symmetry ``z^2`` variant → Bayliss–Turkel
  outer radiation boundary condition.
- **Puncture** — 4th-order Detweiler–Whiting ``\Phi^P`` (``A_{ijk}`` tensor
  basis), m-mode integrals over `BigSpecial.legendre_P3`, and analytic
  derivatives of the mode function: first (`mode_integral_Rn_drho2`,
  `puncture_dphi_m`, `puncture_dpsi_m`) and second
  (`mode_integral_Rn_d2rho2`, `puncture_d2phi_m` — the fused
  ``(I_m, \partial I_m, \partial^2 I_m)`` bundle over one ``J_n`` ladder; the
  flux form is diagonal in ``(r, z)``, so no mixed second derivative exists).
- **EffectiveSource** — the direct pointwise m-mode effective source
  ``S^{\rm eff}_m = -\triangle_m \Psi^{\mathcal P}_m`` (`effective_source_m`),
  i.e. the fixed source of ``\triangle_m \Psi^{\mathcal R}_m = S^{\rm eff}_m``;
  `puncture_psi_chart`, the chart-complete ``\Psi^{\mathcal P}_m`` bundle
  carrying the whole m-mode convention transformation chain (the field-point
  rotation ``e^{-im\chi(r)}``, the ``\sin^{|m|}\theta/r`` factoring and the
  ``\Sigma/(r^2+a^2)`` projection scale); and `effective_source_bundle`, which
  adds the worldtube's
  ``(\Psi^{\mathcal P}, \partial\Psi^{\mathcal P})`` data plus the
  near-worldline roundoff-loss diagnostic (`term_scale`, `cancellation`).
  Alongside it, the **near-worldline series** and the route **switch**: the
  direct assembly cancels four orders of singularity, so its relative accuracy
  degrades like ``\varepsilon\,\lambda^{-2}`` towards the worldline, and the
  paper's power series — accurate through ``\mathcal O(\lambda^2)`` with the
  ``\log\varrho`` terms and the ``H_{-m-1/2}`` m-dependence, whose eleven
  coefficients `effective_source_series_coefficients` returns and
  `effective_source_series_m` evaluates — takes over there (it is finite even
  ON the worldline circle). `effective_source_switched_bundle` /
  `effective_source_m_switched` are the combined evaluator: they compare the
  series' predicted truncation with the direct route's predicted roundoff and
  pick the smaller, so the crossover moves inward like
  ``\varepsilon^{1/5}`` as the working precision rises. The truncation estimate
  carries an unconditional geometric floor ``((1+|m|)\varrho/r_0)^3``: the
  retained ``\lambda^2`` bracket is a signed sum whose zeros would otherwise
  make the series look exact arbitrarily far from the worldline. The DG
  fixed-source assembly calls `effective_source_m_switched`.
- **Worldtube** — the worldtube formulation of the regularized region
  (``\triangle_m \Psi^{\mathcal R}_m = S^{\rm eff}_m`` inside,
  ``\triangle_m \Psi_m = 0`` outside,
  ``\Psi^{\mathcal R}_m = \Psi_m - \Psi^{\mathcal P}_m`` across the
  interface) as pure data for the DG layer: `worldtube_data(dom, orbit, config)`
  assembles a `WorldtubeData` in one pass over the sealed effective-source
  bundle. It carries the region partition (`inside_worldtube`,
  `worldtube_elements`, `worldtube_boundary_faces` — the interior faces between
  an inside and an outside element, with their mortar geometry from
  `worldtube_face_points`); the jump data
  ``\pm\Psi^{\mathcal P}_m`` / ``\mp n_i\mathcal{A}^{ij}\partial_j\Psi^{\mathcal P}_m``
  as one `FaceState` per mortar (`worldtube_jump_data`, signed by
  `worldtube_jump_signs` — the single normative statement of the sender-normal
  sign table — and packable as EllipticCore `FixedDataHook`s through
  `worldtube_hooks` (per FACE; a 2:1 face with differing per-mortar data raises
  — `Solve`'s `jump_hooks` binds them per mortar instead), so they
  reach the right-hand side via `operator_on_zero`);
  and the fixed source ``+S^{\rm eff}_m`` inside and exactly zero outside, in
  the nodal (`worldtube_source`) and the massive ``M\cdot S``
  (`worldtube_massive_source`) form. `worldtube_jump_values` is the pointwise
  kernel and the home of the ``z^2`` chain rule
  ``\mathcal{F}^{z^2} = A^{z^2}\partial_w\Psi``, ``\partial_w = \partial_z/(2z)``.
- **Domain** — `RunConfig` runtime configuration (`crp_paper_config`,
  `CRP_PAPER_MMODES`) and the EllipticCore `Mesh2D` block lattice:
  `block_domain(orbit, config)` instantiates the `BlockDomain` — 5 radial × 2
  angular irreducible blocks in the ``(r, z^2)`` production chart (5 × 3 in the
  full ``(r, z)`` chart), linear coordinate maps everywhere except the outer
  wavezone column's ``1/r`` compactification — together with
  `block_lattice_bounds` (the yaml-derived bounds and the angular chart),
  `block_lattice_index` (the yaml block numbering ``i + n_r j``),
  `block_region` (the vtu ``v``/``t``/``u`` tag of a block),
  `worldtube_blocks` (the regularized-region block), `boundary_face_classes`,
  `boundary_class` and `boundary_class_rule` (`:horizon`, `:outer`,
  `:equator`, `:pole` — the classification the ghost table of `Solve` shares
  with the lattice, so the two cannot drift), and the `block_coefficients` /
  `domain_coefficients` bundles the solver hands to `MeshOperator`.
- **Solve** — the two glue tables a `MeshOperator` over the block domain needs.
  Boundary conditions (the `bcs` table):
  `boundary_conditions(dom, orbit, config; mesh = dom.mesh)` returns the
  `MeshOperator` `bcs` table of the block domain, covering exactly
  `boundary_faces(mesh)`. The `:outer` faces get the sealed
  `outer_boundary_ghost` — a `BaylissTurkelGhost` whose Robin data is evaluated
  at that face's own quadrature points — and every other face gets `NoneGhost`,
  the port of spectre's `None` boundary condition: the convention-N7
  transmissive exterior ``(\psi_{\rm int}, -\varphi_{\rm int})``, through which
  the numerical flux reduces bitwise to the interior flux so the face imposes
  nothing. That is what the vanishing perpendicular flux coefficient at the
  horizon (``\mathcal{A}^{rr}(r_+) = 0``) and at the angular boundaries
  (``\mathcal{A}^{z^2}(0) = \mathcal{A}^{z^2}(1) = 0``) demands, and it is also
  the ``m = 0`` outer branch of the production input file. `mesh` may be any
  coarser mesh of the same lattice, which is what a `MultigridHierarchy`
  `build_op` factory needs.

  Jump conditions (the `hooks` table): `jump_hooks(dom, wt, orbit, config, m;
  mesh = dom.mesh)` binds BOTH interior jump conditions of main.tex:596-608 onto
  the EllipticCore `ModifyBoundaryData` seam, one hook **per mortar** in the
  receiver's mortar order — the granularity the production lattice's ``2^L`` 2:1
  worldtube faces need. The sealed `worldtube_jump_data` become `FixedDataHook`s:
  trial-field independent, so they leave the linearized matrix untouched
  (`build_matrix` is bitwise unchanged) and reach the right-hand side through
  `operator_on_zero`. `vtu_junction_faces` tags every interior face that crosses
  the vtu slicing — spectre's `NullSlicingBlocks` rule over the sealed
  `block_region`, cross-checked against the transition radii ``r_v``, ``r_u`` —
  with its sealed `VtuJunction`, and `vtu_junction_hook` turns
  ``[\partial_{r_*}\Psi_m] = -im\Omega\Psi_m`` into the `LinearizedHook`
  ``(n\cdot\mathcal{F})^{\rm ext} \to (n\cdot\mathcal{F})^{\rm ext} -
  im\Omega\langle\!\langle\Psi_m\rangle\!\rangle``, which IS in the GMRES
  matrix. No chart conversion is applied: ``A^{rr} = \Delta/(r^2+a^2)`` is
  ``dr/dr_*``, so ``n\cdot\mathcal{F}`` already is the tortoise derivative. A
  face carrying both composes as `ComposedHook(linearized, fixed)` in the sealed
  W5-B order (it cannot happen on the CRP lattice — the tube spans
  ``[5.2, 14.8]`` and the transitions sit at 2.8 and 24.4 — and is built
  generically anyway). Callable per level mesh like `boundary_conditions`; the
  mesh-bound `WorldtubeData` must be rebuilt by the caller for that mesh.
- **Extraction** — Heffernan m-mode regularization parameters of the radial
  self-force for circular equatorial Kerr orbits (arXiv:1211.6446), in the
  folded ``m \ge 0`` ×2 real-part convention:
  - `reg_param_Fr(order, orbit)` — the m-independent circular-Kerr
    coefficients ``F_{r[2]}``, ``F_{r[4]}`` (`HEFFERNAN_CIRCULAR_FR_ORDERS`),
    purely algebraic in ``(a, r_0, M)``; the printed ``a \to 0`` Schwarzschild
    limits are given in the docstring.
  - `mmode_reg_factor_over_pi(order, m)` — the printed m-dependence
    ``\pi F^m_{r[2j]} / F_{r[2j]}`` as an exact `Rational{BigInt}`, for every
    published order `HEFFERNAN_MMODE_ORDERS = (2, 4, 6, 8)`;
    `mmode_reg_factor([T,] order, m)` is the same divided by ``\pi``.
  - `mmode_reg_Fr(order, orbit, m)` — the assembled per-mode subtraction term
    ``F^m_{r[order]}``.
  - `mmode_reg_factor_tail_over_pi(order, mstart)` /
    `mmode_reg_factor_tail([T,] order, mstart)` /
    `mmode_reg_tail_Fr(order, orbit, mstart)` — the closed-form analytic tail
    ``\sum_{m \ge m_{\rm start}} F^m_{r[order]}`` that completes a truncated
    mode sum.

Planned, each landing written test-first once it is grounded: the rest of
`Solve` (single-mode driver, hp-AMR loop, process-parallel m fan-out), and the
rest of `Extraction` (``F_r^m`` extraction from a solved mode, and the mode sum
itself).

## Conventions (locked package-wide)

- Geometric units ``G = c = 1``, background mass ``M = 1``; Kerr spin
  ``a \in (-1, 1)``, retrograde orbits encoded via ``a < 0`` with ``m \ge 0``
  and the ×2 real-part mode-sum convention.
- Mode ansatz
  ``\Phi = \sum_m \Psi_m(r, z)\, \sin^{|m|}\theta / r\, e^{i m (\bar\varphi - \Omega t)}``
  with ``z = \cos\theta`` and the horizon-regular azimuthal chart
  ``\bar\varphi = \varphi + \frac{a}{r_+ - r_-}\ln\frac{r - r_+}{r - r_-}``.
- Radial chart: horizon-penetrating Boyer–Lindquist ``r`` via
  ``\partial_{r_*} = \frac{\Delta}{r^2 + a^2}\partial_r``; vtu slicing
  (v→t→u null/Boyer slabs) with piecewise boost ``H \in \{-1, 0, +1\}`` and
  junction conditions ``[\partial_{r_*}\Psi] = -i m \Omega \Psi``.
- Equatorial symmetry via the ``z^2`` variant on ``z \in [0, 1]``.

## Dependencies and poison-axis declaration

Binding dependency budget, enforced by `test/hygiene.jl`: exactly
`{LinearAlgebra, Printf, NumericsCore, BigSpecial, Spacetime, EllipticCore,
Data}`. Special functions come only through BigSpecial (no direct
SpecialFunctions / HypergeometricFunctions); the EffectiveSource C oracles are
baked test fixtures, never live deps (no GSL/ccall channel); GPU stacks are
banned.

## Type policy

Type-generic `{T}` throughout so BigFloat precision is inherited from
NumericsCore/BigSpecial (never hard-Float64'd); the production path is
ComplexF64.

## Public API

```@autodocs
Modules = [ScalarSelfForce]
Order = [:module, :type, :constant, :function]
```

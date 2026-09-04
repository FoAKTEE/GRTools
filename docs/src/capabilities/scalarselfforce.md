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
  bundle;  `mesh = dom.mesh` may be any mesh of the same lattice, so a
  multigrid or AMR level rebuilds the mesh-bound data in one call. It carries the region partition (`inside_worldtube`,
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
- **Solve** — the two glue tables a `MeshOperator` over the block domain needs,
  and the assembly that composes them.
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

  Assembly: `assemble_operator(T, dom, wt, orbit, config, m; mesh = dom.mesh)`
  is the whole CRP operator in one call — `MeshOperator` over
  `domain_coefficients`, `jump_hooks`, `boundary_conditions` and
  `config.dg.penalty_parameter` — written once and gated BITWISE
  (`build_matrix`, `operator_on_zero`, `apply!`) against the hand-composed
  spelling, so the figure drivers, the single-mode driver and the AMR driver
  cannot drift apart. `operator_builder(T, dom, orbit, config, m; worldtube)`
  is the `MultigridHierarchy` `build_op` factory over a caller-supplied
  per-level `WorldtubeData` factory (rebuilding that mesh-bound bundle stays
  `worldtube_data`'s own `mesh` keyword — it is not duplicated here), and
  `coarsest_operator(...; worldtube, maxlevels = config.multigrid.initial_coarse_levels + 1)`
  walks `coarsen_mesh` to the multigrid COARSEST level and re-assembles there.
  That default is the yaml twin of `InitialCoarseLevels: L` (which counts only
  the coarser levels), so at the production ``L = 0`` the coarsest level IS the
  initial 20-element / 320-dof lattice that Fig. 2 of arXiv:2606.04998 plots;
  `maxlevels = nothing` coarsens instead to the irreducible block structure
  (10 elements / 160 dofs on the 5×2 lattice, 15 / 240 on the 5×3).
  `linearized_rhs(op, wt)` folds the worldtube fixed source with the operator's
  own fixed-data lift, ``b = M\cdot S - A(0)``, extending EllipticCore's
  general `linearized_rhs(op, b_phys)` method. Its result is ALREADY folded, so
  hand it to `solve_linearized!`, never to `solve!` (which folds for you) —
  double-folding is an O(1) wrong answer that *converges*, so the solver reports
  success on the wrong problem.

  Single-mode solver-stack driver: `solve_mode(orbit, config, m)` is the whole
  m-mode solve — block domain, mesh-bound worldtube bundle, operator,
  `SolverStack` (complex GMRES ∘ one multigrid V-cycle ∘ two additive-Schwarz
  passes per level ∘ a cached direct bottom solve) over `operator_builder`'s
  per-level factory, the physical right-hand side, and the solve — returning a
  `ModeSolution` whose `u` is ``\Psi^{\mathcal R}_m`` inside the worldtube and
  ``\Psi_m`` outside. `solver_options(config)` is the field-by-field
  `SolverStackOptions` twin of the sealed `RunConfig` (the yaml solver block,
  lines 178–205), exposing the two multigrid level-policy knobs `coarse_levels`
  and `coarse_orders` that an AMR sweep must choose; the complex shift it passes
  is `helmholtz_shift(m, Ω; slicing = :vtu)`, exactly zero under the paper's
  null slicing. `mode_fixed_source` builds the VOLUME source through the generic
  `worldtube_massive_source` seam with the **switched** evaluator
  `effective_source_m_switched` inside the tube and exactly zero outside, while
  the worldtube **mortar** data stay on the direct bundle — the one normative
  statement of that rule is its docstring. Read-out: `mode_residual` (the DG
  residual against the FULL operator table), `reconstruct_full_field`
  (``\Psi_m`` everywhere, ``\Psi^{\mathcal R}_m + \Psi^{\mathcal P}_m``
  inside), `field_at` (tensor-product Lagrange interpolation on the element
  containing a chart point, with an `element` keyword for a one-sided trace)
  and `field_grad_at`, its chart gradient
  ``(\partial_{x_1} f, \partial_{x_2} f)`` from the same reference element's
  differentiation matrix — the ``\partial_r \Psi^{\mathcal R}`` at the
  particle that self-force extraction needs.
  `rebuild_for_mesh!` is the one-step per-level rebuild an hp-AMR sweep calls —
  worldtube bundle, operator, factory hand-over and right-hand side, with no
  epoch bump because a pure mesh change is already in both cache keys.
  `resolve_max_overlap(config)` maps `MaxOverlap: Auto` (yaml:198) onto the
  smoother's `overlap` using spectre's own semantics — `Auto` means "no
  restriction, the overlap covers the full neighbour"
  (`Schwarz/Tags.hpp:29-35`) and resolves numerically as
  ``\min(\text{volume extent}, \text{max\_overlap}) =`` the neighbour's
  extent (`Schwarz/OverlapHelpers.cpp:21-24`), i.e. **4** at the production
  ``P = 4``. Any field that changes the solve and cannot be represented is
  *refused*: a non-zero `Multigrid.ComplexShift` (yaml:194) and an unmappable
  `MaxOverlap` both raise rather than being silently dropped.
  Measured at the paper's own configuration (``a = 0.5``, ``r_0 = 10M``,
  ``m = 2``, ``L = 0``, ``P = 4``): converged to the yaml's absolute residual
  ``10^{-12}`` in **1** preconditioned GMRES iteration at a true residual of
  ``1.2\times10^{-15}`` — the paper's own first-AMR-level behaviour, since
  `InitialCoarseLevels: L = 0` makes the V-cycle a direct bottom solve; with a
  genuine two-level hierarchy on the same mesh, 9 iterations at the resolved
  overlap (10 at the package default 2).

  Under the paper's own `RefineAtPuncture` criterion (yaml:132,
  main.tex:552-560 — split only the elements *containing* the particle) over
  five AMR levels (320 → 593 → 950 → 1411 → 1988 grid points), the counts are
  9, 8, 7, 8, 6 with the min coarsening rule (coarsest grid
  160 → 236 → 334 → 450 → 584) and 9, 8, 9, 9, 9 with spectre parity
  (`coarse_orders = P`), which pins the bottom grid at 160 points and builds its
  LU once. (These are the numbers after the EllipticCore `hp_refine`
  two-to-one-closure fix; this page is the one maintained measurement table.) **main.tex:720's "less than 10 iterations per AMR level"
  holds, and the count decreases with refinement depth.** Refining the *whole*
  worldtube instead — which is not the paper's criterion — makes the counts
  grow (9, 8, 9, 10, 11 and 9, 9, 11, 12, 18), so the refinement criterion, not
  the solver, is what decides.

  hp-AMR driver: `solve_mode_amr(orbit, config, m)` runs the paper's
  `Amr: Iterations` levels of "``h``-refine every element that contains the
  puncture, ``p``-refine everything else" (main.tex:552-560) over ONE
  `ModeSolution` — one stack, one epoch, one pair of caches, each level
  warm-started from the transferred previous solution — and returns an
  `AmrHistory`: per level the `Mesh2D`, the puncture elements, the `SolveRecord`
  (iterations, residual staircase, `build_time`, cache counters) and the
  particle diagnostics ``\Psi^{\mathcal R}_m(x_p)``,
  ``\partial_r\Psi^{\mathcal R}_m(x_p)``, plus the final `ModeSolution`.
  `amr_flags` is the criterion: spectre's `RefineAtPuncture` — map the puncture
  into each block, snap ``|\xi_d| < 10^{-10}`` to the block centre or face, and
  keep every element whose CLOSED dyadic interval contains it, so the puncture
  sitting on an element corner flags every element that touches it (two of them
  in the ``z^2`` chart, where the equator is the domain boundary) — merged with
  `IncreaseResolution` exactly as spectre merges criteria, `Split` outranking
  `IncreaseResolution`, which is the sealed `hp_refine` one-container contract.
  The yaml `Limits` (`Auto` = 15 refinement levels, 20 grid points) are enforced
  with `ErrorBeyondLimits: True` semantics: an `ArgumentError`, never a silent
  clip. `amr_level_options(config, step)` is the multigrid level policy — the
  yaml's `InitialCoarseLevels` at level 0, so the first level converges
  immediately on the direct bottom solve (main.tex:701), then coarsening to the
  irreducible block structure, a grid that does not move as AMR adds finer
  levels, so the bottom LU is factorized once and hit at every later level
  (main.tex:674). `amr_step!` is one level, `puncture_chart_position` is
  ``(r_0, 0)`` in the discretization chart. The per-level particle diagnostics
  are read through the sealed `field_at` / `field_grad_at` above, so the driver
  adds no read-out surface of its own.
  Measured over the 8-level sweep at ``a = 0.5``, ``r_0 = 10M``, ``m = 2``:
  20 → 69 elements and 320 → 5696 dofs, GMRES **1 then 8–9** iterations per level
  (main.tex:720's "less than 10" holds at this configuration; 10 is measured at
  one level of ``a = 0`` and at three levels of ``m = 5``), bottom LU built twice
  and hit seven times, and successive differences of
  ``\partial_r\Psi^{\mathcal R}(x_p)`` falling geometrically
  (``8.5\times10^{-5} \to 3.7\times10^{-10}``, ratios 0.06-0.33) — the whole
  sweep in 18.9 s of warm single-threaded compute, which needs
  `OPENBLAS_NUM_THREADS=1` (the dense subdomain solves are tiny, and the host's
  default thread pool costs a factor 8.9). `resolution_error(hist)` and
  `relative_resolution_error(hist)` are the paper's own resolution-error
  estimate (main.tex:730), the difference between the highest and
  second-highest AMR level, and `AmrLevel.series_nodes` records how many in-tube
  nodes select the near-worldline series route — zero throughout the paper's
  eight iterations, first nonzero at level 12.
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
  - `mmode_reg_term_Fr(order, orbit, m)` — THE one helper for the
    folded/unfolded bookkeeping: the term to subtract from an
    Osburn–Nishimura Eq. (3.10) force mode, i.e. ``F^m_{r[order]}`` for
    ``m \ge 1`` and half of it at ``m = 0`` (Eq. (3.10)'s ``m = 0`` branch is
    unfolded while the parameter carries the ×2 at every ``m``). Gated against
    the vanishing of the doubly-infinite parameter sum.
  - `mode_self_force(orbit, m, psi, dpsi_r, dpsi_z)` /
    `mode_self_force(sol::ModeSolution)` /
    `mode_self_force(hist::AmrHistory; level, reduction)` — ``F_r^m``
    extraction at the particle, Osburn–Nishimura Eq. (3.10) (arXiv:2206.07031)
    in the locked ansatz conventions, in the normalized ``(M^2/q^2)F_\alpha``
    of main.tex:428-430 (``q = M = 1``): unfolded at ``m = 0``, and for
    ``m > 0`` the ×2 real part of ``e^{+im(\bar\varphi_p - \Omega t)}
    [\partial_r\Psi^{\mathcal R}_m/r_0 - \Psi^{\mathcal R}_m/r_0^2
    + i m (a/\Delta)\Psi^{\mathcal R}_m/r_0]`` — the last term is the
    dephasing derivative the horizon-regular azimuth contributes to
    ``\partial_r``. The radial chart is horizon-penetrating Boyer–Lindquist
    ``r``, so spectre's ``\alpha = \Delta/(r^2+a^2)`` tortoise factor is
    deliberately absent and a configuration declaring
    `PenetratingHorizon = false` throws. Beyond spectre's radial observable,
    ``F_t^m = -\Omega F_\varphi^m`` mode by mode and ``F_\theta^m`` follows
    from ``\partial_\theta = -\sin\theta\,\partial_z`` (a bitwise zero in the
    production ``z^2`` chart). Ingredients are read through the sealed
    `field_at` / `field_grad_at`. The `ModeSolution` method reduces by the
    **mean over every element containing the puncture** — spectre parity
    (`ObserveSelfForce.hpp:89-107`) — and `element = e` selects one one-sided
    trace instead; the `AmrHistory` method uses `:trace` at every level so a
    convergence sequence is homogeneous.
  - `mode_self_force_history(hist)` / `mode_force_convergence(hist)` — the
    per-AMR-level ``F_r^m`` sequence and its successive-difference table,
    including `resolution_error`, the error budget's numerical-resolution
    floor (main.tex:722), with the same definition as the `Solve`
    `resolution_error(hist; field)` helper. Observed at ``a = 0.5``,
    ``r_0 = 10M``, ``m = 1`` over the paper's 8 AMR levels: ``F_r^m``
    converges to ``-1.3652\times10^{-4} + 9.1381\times10^{-5}\,i`` with a
    resolution error of ``4.6\times10^{-11}`` (``3.4\times10^{-7}``
    relative), reproducing spectre's own regression row
    ``[2, 975, 4, -1.365\times10^{-4}, 9.138\times10^{-5}, 0, 0]`` to
    ``1.7\times10^{-8}``, inside its ``2\times10^{-7}`` tolerance
    *[PRELIMINARY: mesh-dependent, re-measured after 6df9f51; per-level values
    move with the hp-AMR order map, so gate on tolerances and on the shape of
    the sequence, never on a golden value]*.
  - `regularized_mode_force(x, m; order = 4)` — ``F_r^m - F^m_{r[order]}``
    through `mmode_reg_term_Fr`, the per-mode term the mode sum adds up.
    Observed decay of the residual at ``a = 0.5``, ``r_0 = 10M``, 6 AMR
    levels: log-log slope ``-6.007 \pm 0.002`` over the resolved pair
    ``m \in \{2, 5\}`` *[PRELIMINARY]*; the ``m = 10`` residual is only
    ``2.8\times`` its own level-to-level difference at that depth, so fits
    that include it (``-6.32 \pm 0.52`` over ``\{5, 10\}``) are
    resolution-limited. A discharge of the ``m^{-6}`` claim is *[FUTURE]*:
    it needs the process fan-out's 8 levels with ``m`` up to 20 and the same
    propagated uncertainty.

- **Solve — the process-parallel m fan-out** (`solve_modes_parallel`,
  `solve_modes_serial`, `mode_pool`, `mode_config`, `mode_pool_size`,
  `system_load`). The paper's ``m = 0\ldots20`` sweep (`CRP_PAPER_MMODES`,
  main.tex:716 "Each m-mode is independent") run as **one independent process
  per mode**: each worker rebuilds its own `mode_config` — only
  `MModeNumber` and the ``m = 0`` no-outer-BC yaml branch move — runs the
  sealed `solve_mode_amr` with its own stack and its own two caches, pins
  BLAS explicitly, and ships back a `ModeResult` (every `AmrLevel` plus the
  final `sol.u`; the live `ModeSolution` is process-local by the row-42 rule
  and deliberately does not travel — 597 KB per mode against 2.5 GB of
  process memory). The fan-out is over **processes, never threads**: the
  solver stack and its caches are neither thread-safe nor reusable across
  modes. `solve_modes_serial` is the same loop in-process and is the
  **bitwise oracle** of the fan-out — measured bit-identical at the same BLAS
  and Julia thread counts, and different at roundoff only (``3.2\times
  10^{-12}`` relative) when the BLAS thread count changes, which is why the
  pin is applied inside every worker rather than inherited.
  The owner's compute policy is enforced mechanically and stated once, in the
  `solve_modes_parallel` docstring: `mode_pool_size` keeps
  ``n_{\rm workers} \times \max(n_{\rm threads}, n_{\rm BLAS}) \le``
  `PARALLEL_CPU_BUDGET` = 24 for every input and reports every clamp (a CPU
  spent inside OpenBLAS is a CPU the owner lent us), and `system_load` reads
  the live whole-system usage — the worst of `PARALLEL_LOAD_PROBES` = 3
  ``/``proc``/``stat windows — so the pre-launch guard can refuse with the
  numbers when `live + planned` would reach `PARALLEL_LOAD_GUARD` = 0.90.
  **The guard governs the launch decision only.** A monitor task samples the
  machine into `run.samples` (and an optional log file) for the whole run, and
  `run.max_load` is a *record* of what the shared box then did: a crossing is
  warned when the run ends, never silently gated, because nothing in this
  package can stop another user of a shared host once we have launched.
  `mode_pool` is the reusable pool: a fresh worker pays a measured 166 s of
  JIT for its first solve and 45 s more for its first AMR step, against
  15.9 s for a warm 8-level sweep, so several sweeps must share one pool.
  The fan-out is **JIT-bound, not solve-bound**, so its cost has three
  regimes and only the third one is a win: *cold* (a fresh pool per sweep),
  *warm-serial* (one process that has already compiled), and *reused pool*.
  The measured numbers are being re-established on this tree by the 21-mode
  artefact of `evidence/mmode-parallel-driver/` and are quoted there.
  A mode that throws is reported by id in `set.failures`, warned with its
  remote backtrace, and never swallowed; the other modes still return. So is a
  worker that **dies** — `pmap`'s `on_error` does not catch a
  `ProcessExitedException` (measured: it escapes and takes the whole sweep
  down), so the fan-out schedules its own modes, attributes a death to its
  mode and its worker id, and prunes the corpse from the pool; `isopen(pool)`
  asks the cluster, not a flag. `timeout_s` bounds a hung sweep, cancelling
  the outstanding modes the only way a remote call can be cancelled — by
  removing their workers.

Planned, each landing written test-first once it is grounded: the rest of
`Extraction` (the mode sum itself).

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
`{LinearAlgebra, Printf, Distributed, NumericsCore, BigSpecial, Spacetime,
EllipticCore, Data}`. `Distributed` is a **stdlib**, widened into the budget by
the wave-12 orchestrator decision as the IPC of the row-47 m-mode process
fan-out (`src/Solve/parallel.jl`) — stdlib only, no registry package, the same
rule that admitted SparseArrays/SuiteSparse into EllipticCore.
Special functions come only through BigSpecial (no direct
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

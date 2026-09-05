# Changelog

All notable changes to GrTools are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/). Every package in the suite shares
the suite version.

## [Unreleased]

### Added

- **`EllipticCore` — the solver stack is complete: geometric multigrid, the cached
  sparse bottom solve, `SolverStack`, and hp-AMR.** `MultigridHierarchy` h-coarsens a
  `Mesh2D` to the irreducible block structure one `coarsen_mesh` step at a time
  (anisotropic-safe: 300/300 valid, balanced, nested, exactly-tiling hierarchies over an
  AMR-reachable mesh census) or takes an **explicit chain** of coarse meshes,
  `MultigridHierarchy(op, build_op, meshes)`, which is the form arXiv:2606.04998 runs —
  its executable keeps the coarse grids, so the multigrid levels *are* the AMR history.
  `MultigridV` is the V-cycle preconditioner (measured contraction 0.171, level-independent);
  `sparse_bottom_solver` fills its bottom seam with a stdlib sparse `lu` (UMFPACK, no new
  dependency) cached across V-cycles *and* AMR rebuilds — one build and seven hits over an
  eight-step run, which is the paper's own claim. `SolverStack` composes the whole linear
  solver (right-preconditioned complex GMRES ∘ one V-cycle ∘ two additive-Schwarz passes ∘ a
  cached direct bottom solve) from the CRP's own solver input, with `rebuild!` carrying it
  across an AMR step and ONE caller-owned `epoch` governing both caches. AMR adds `h_split`,
  `balance_two_to_one`, `balanced_flags`, `h_refine`, `p_refine`, `hp_refine` and
  `transfer_solution` — spectre's hp-AMR reproduced exactly, including the 2:1 refinement
  decision closed in **every** dimension and the rule that a balance-split element holds its
  polynomial order in the split dimension.

- **`ScalarSelfForce` — the m-mode reproduction of arXiv:2606.04998, end to end.**
  The Detweiler–Whiting effective source (`effective_source_m`, its near-worldline power
  series `effective_source_series_m`, and the precision-aware `effective_source_m_switched`
  route switch); the worldtube formulation as pure DG data (`worldtube_data`,
  `worldtube_jump_data`, `worldtube_massive_source`); the CRP block lattice
  (`block_domain`, `crp_paper_config`, `radially_compressed`); the glue tables and assembly
  (`boundary_conditions`, `jump_hooks`, `assemble_operator`, `operator_builder`,
  `linearized_rhs`); the drivers (`solve_mode`, `solve_mode_amr` with its `AmrHistory`,
  `rebuild_for_mesh!`, `field_at` / `field_grad_at`) and the process-parallel m fan-out
  (`solve_modes_parallel`, `solve_modes_serial`, `mode_pool` — one process per mode,
  bit-identical to the in-process oracle over the paper's whole 21-mode sweep); the
  Heffernan m-mode regularization parameters and their analytic tail; the per-mode
  extraction `mode_self_force`; and the mode sum itself, `total_self_force` /
  `mode_sum_tail` / `mode_decay_fit`. Also `horizon_branch_index` /
  `horizon_axis_profile` / `horizon_order_hold`, the instrument and the opt-in cure for the
  near-horizon non-analytic branch (default off, for reference parity).

- **`reproduce/arXiv-2606.04998/` — six of the paper's seven figures.** Thin drivers and
  JSON goldens for `fig:domain`, `fig:matrix`, `fig:convergence`, `fig:modes`, `fig:linsolv`
  and `fig:ssf_error_budget`, registered in `reproduce/MANIFEST.toml` (all `heavy` tier) and
  re-checked by
  `julia --project=reproduce reproduce/run_all.jl --check --paper 2606.04998`. The hp-AMR
  sweep reproduces the paper's published grid sequence exactly, level by level — 480, 975,
  1604, 2427, 3468, 4751, 6300 points in 30, 47, 59, 71, 83, 95, 107 elements — with no fitted
  constant, and the summed `(M²/q²)F_r` at `r₀ = 10M` matches the printed Warburton–Barack
  Table III at `a/M = 0, ±0.5, ±0.9` to a worst relative deviation of 6.4e-6 (their Table II
  `F_t`, which needs no regularization and no tail, to 3.1e-8). `fig:precond`, the
  preconditioner schematic, is not claimed.

- **`fig:ssf_error_budget` (Fig. 7) — the `F_r` error budget over spin, landed
  `preliminary` with its failures pinned.** Reading the published figure's own vector paths
  established that it samples nine spins rather than the five it labels, and that its
  reference curve carries only the seven with `|a| ≤ 0.9` — Warburton–Barack's tables stop
  there at every radius. Reproduced: the published truncation curve to 0.3–3.6 % at 17 of 18
  configurations, the reference curve to within 23 % at all 8 printed comparison points,
  `|F_r|` to ≤ 0.19 % across the strip, and `F_t` against the printed nine-figure Table II to
  2.91e-8. Not reproduced, and pinned configuration by configuration in `gate_failure_ids`
  rather than scoped away: the relative truncation error at `a = +0.998` at its ISCO (121× the
  published curve) and the relative resolution error at four configurations (more than a decade
  above it). `gates_pass` is a REGISTERED value, so a gate that starts failing — or a known
  failure that disappears — is a golden mismatch that fails `--check`. The figure also narrows
  the error-budget ordering finding: on the paper's own full-z lattice the resolution error
  falls below truncation at 15 of 18 configurations (published 18/18) against 0 of 9 on the
  shipped half domain, with the 1.69× dof difference disclosed as a confound rather than
  argued away.

- **The GMRES count does NOT grow with the mode number** — a `[SOLID]` claim of an earlier
  wave, retracted after re-measurement. The whole `m = 0…20` sweep was re-run at the paper's
  own nine AMR grids on the current tree: at the finest grid the counts are 3, 3, 3, 3, 3, 3,
  2, 2, 2, 2, 2, 2, 3, 3, 4, 3, 3, 3, 3, 4, 4 on the production half domain and 2 … 3 on the
  full-z lattice, and **the maximum over every level and every mode is 7 on both**. So
  main.tex:720's "less than 10 preconditioned iterations per AMR level" holds for *every* mode
  of the sweep, not merely for the paper's one exemplary mode, and the count falls with
  refinement depth at every `m`. The earlier 9 … 20 growth was measured before the
  all-dimension 2:1 closure moved the iteration counts along with the meshes.

- **Two source errata located in arXiv:2606.04998, both diagnosed rather than absorbed.**
  (1) The last term of the printed radial effective potential (main.tex:190) has the wrong
  sign: `dV_eff/dr` of the printed form does not vanish on the paper's own circular-orbit
  `(E, L)` (normalized residual ~0.5 at all 35 fixture points) while the corrected circularity
  condition holds to 2.1e-15 at Float64 and 8.7e-78 at BigFloat-256; every other closed form
  the paper prints is consistent with the correction, so it is a print typo with no downstream
  consequence. (2) The shipped input file sets `ImposeEquatorialSymmetry: True`, but the
  published figures were produced **without** it, on the full-z 30-element / 480-point lattice
  — established by vector extraction of two published figures (`fig:convergence`'s abscissa
  starts at N = 480 for all 21 curves; `fig:linsolv`'s level-0 residual 1.394 is the full-z
  `‖b‖₀`, ours 1.391, against the half domain's 1.744). The package default follows the paper's
  file; `crp_paper_config(; equatorial_symmetry = false)` builds the figure lattice.

- **Documentation: `EllipticCore` and `ScalarSelfForce` are wired into the site**, with a
  `Reproducing arXiv:2606.04998` topic page indexing the drivers, the goldens, the
  `run_all.jl --check` recipe, the two errata and the open holes.

- **`AccretionDisc/derivation` — route R3, the Julia-native derivation of the
  arXiv:2608.17003 master system.** A standalone `Symbolics.jl` project (own
  `Project.toml`, **not** a dependency of `AccretionDisc`) that derives the
  first-order master system (47), the eliminated velocities and the
  angular-momentum diagnostics (72a)/(72b) from the paper's own equations —
  line element (4), fluid equations (3), linearisation (15)-(22), Fourier
  ansatz (45)-(46), the reduction of M:604-608 — certifies every step
  symbolically or by exact `Rational{BigInt}` evaluation (residual exactly
  `0//1`), and emits `src/SpiralWaves/master_derived.jl` and `r2_derived.jl`
  with a deterministic printer whose output is re-parsed and re-checked
  exactly. `scripts/derive_all.jl` runs the derivation, `scripts/emit_sources.jl`
  the emission.

- **`EllipticCore.helmholtz_shift` / `EllipticCore.schwarz_preconditioner`**
  (`core/EllipticCore/src/schwarz_shift.jl`) — the Schwarz complex-shift
  *policy*: `helmholtz_shift(m, Ω; slicing)` returns `s = (mΩ)²` under
  `t`-slicing and `0` under `:hyperboloidal` / `:vtu`, and
  `schwarz_preconditioner(sm; passes = 2)` adapts a `SchwarzSmoother` to the
  `Minv!` right-preconditioner hook of `gmres`. On this package's *massive*
  Gauss operator the shift is applied mass-weighted — build the smoother on the
  damped coefficient `β − i·s` so its subdomain matrices are `A_ss − i·s·M_ss`
  — not as the unweighted `−i·s·I` of the literal `shift =` seam, which at
  operator scale destroys the preconditioner. Gated at a subdomain resonance,
  where the unshifted smoother burns a full 150-iteration budget 2.7e10× above
  tolerance and the policy shift converges in 47.

### Changed — BREAKING: the arXiv:2608.17003 master system is Julia-derived

- **`AccretionDisc.SpiralWaves` no longer consumes the supplementary
  *Mathematica* notebook of arXiv:2608.17003 (nor `wolframscript`) anywhere in
  its build path.** `master.jl` includes `master_derived.jl` and `r2_derived.jl`
  — emitted by the route-R3 Julia derivation above — instead of
  `master_generated.jl` (the notebook transcription) and `r2_generated.jl` (an
  independent Mathematica derivation). Those two files still exist, tracked, as
  `test/SpiralWaves/fixtures/oracle/{master,r2}_generated.jl`, and
  `test/SpiralWaves/test_derived_equivalence.jl` gates the derived sources
  against them forever. Function names, signatures and every public API are
  unchanged; results move by at most Float64 roundoff and are in fact better
  conditioned (worst relative deviation on the 50-point physical Mathematica
  fixture: 2.6e-12 derived vs 7.0e-10 for the notebook form; the diagnostics
  improve from 1e-11/1e-9 to 1e-15/1e-16). **On-shell only:** the three routes
  agree where the background satisfies `W = 1/√(1−v²)` and the hydrostatic
  balance, and carry different off-shell continuations — the off-shell 40-point
  fixture now gates the R1 oracle rather than the build source.

### Removed — BREAKING: the notebook-derived oracle is out of the repository

- **The route-R1 oracle of arXiv:2608.17003 — the transcription of the authors'
  supplementary *Mathematica* notebook — and its off-shell fixture are deleted**
  (owner decision, 2026-08-22, `EGD::jd-notebook-purge`):
  `test/SpiralWaves/fixtures/oracle/master_generated.jl`,
  `test/SpiralWaves/fixtures/master_r1_fixture.jl` (40 points evaluated from the
  notebook expressions), `test/SpiralWaves/fixtures/master_onshell_fixture.jl`
  (baked from the `"R1"` columns of the gate fixture) and
  `fixtures/capture/bake_master_onshell_fixture.py`. This supersedes the note in
  the entry above that those files "still exist, tracked".

- **Equivalence gating was re-sourced onto our own lineage, not weakened.** New
  `test/SpiralWaves/fixtures/oracle/r2_master_generated.jl` provides
  `_r2_master_entries(inp, omt, vp)` — the master entries `C`, `S` of **route R2**,
  our independent Mathematica derivation — emitted by the new emri-disk generator
  `codes/cas/R2_master_to_julia.wls` through the same shared tree walker as
  `r2_generated.jl`, with an exact `Simplify[parse(julia) − original] == 0`
  round-trip on all six entries and byte-identical re-emission. New
  `test/SpiralWaves/fixtures/master_r2_fixture.jl` (+ `bake_master_r2_fixture.py`)
  re-bakes the same 50 physical on-shell points from the `"R2"` columns; its
  numeric payload is byte-for-byte the previous file's, because
  `EGD::R1-vs-R2-gate` had proved all 300 entries equal. Measured on the same 200
  random on-shell points and the same seed: the new oracle reproduces the removed
  one to **3.0e-12** (Float64) / **5.7e-75** (BigFloat 256), and the route-R3 gate
  deviation moved from 3.0e-12 / 1.9e-73 to **1.7e-12 / 1.9e-73** at unchanged
  tolerances.

- **Test-count delta (measured):** the `AccretionDisc` suite goes from 14485 to
  **13914** assertions (13912 pass, 2 broken, green). The −571 nets out as −578
  from the two `test_master.jl` testsets that gated the R1 oracle against the
  off-shell R1 fixture ((a): 1 + 40 points × 6 entries × 2 + 1 inference = 482;
  (a′): 8 × 6 × 2 = 96), −1 from `test_derived_equivalence.jl` (e), and +8 new
  assertions (a permanent 4-way guard in (0) that no removed R1 artefact can
  return, and 4 in (f) covering the new oracle's allocation, inference and
  BigFloat-precision inheritance). **None of the removed assertions gated a
  shipped number** — their inputs are off-shell by construction, where the routes
  provably differ. `test_derived_equivalence.jl` (e) now makes the same FINDING on
  the R2 lineage: R2 eliminated `W` and `ν′`, so it is *exactly* blind to them
  (measured response 0.0) while route R3 moves by 29.5 relative when the same
  points are pushed off the constraint surface, and the two agree to 2.6e-74 on it.

- **Nothing published moved:** `reproduce/run_all.jl --check --paper 2608.17003`
  is green at declared tolerances and `check_printed_numbers` still passes 18/18.
  The R1 lineage, the notebook mirror and the full R1-vs-R2-vs-R3 equivalence
  evidence remain archived in the mission repository (emri-disk
  `ref-code/arXiv-2608.17003-anc/`, `results/emri-gasdisk/paper_EMRI-GASDISK/`
  `evidence/{R1-vs-R2-gate,jd-equivalence-gate,jd-notebook-purge}`), which is the
  canonical archive.

### Changed — de-duplication (one definition per primitive)

- **`NumericsCore.bary_row`** is the single barycentric Lagrange *row* kernel; the private
  copies in `EllipticCore` (`refelement.jl`) and `Perturbation.TimeDomain` (`common.jl`) and
  the inline body of `NumericsCore.lagrange_interp` are gone. `bary_eval`/`bary_row` take a
  `nodetol` node-collision radius (default `0` = the exact `xq == x_j` test).
- **`NumericsCore.secant`** is complex-capable (real tolerances on a complex iterate) and takes
  `throw_on_fail`; the private complex copies `Harmonics._swsh_secant` and
  `Perturbation.QNM._qnm_secant` are gone.
- **`NumericsCore.factorial_int`** is the single integer-product factorial;
  `BigSpecial.factorial_big` forwards to it and `SelfForce.LorenzField`'s `_factT` is gone.
  `SelfForce.LorenzField` consumes `NumericsCore.bc_tol`/`rk_tol` instead of private copies.
- **`Spacetime.SpinningSecondary`** is the single owner of the linear-in-spin orbit-correction
  engine: `Perturbation`'s 662-line `SpinFluxes/spinning_orbit_correction_local.jl` is deleted and
  the BHPToolkit-facing names bind to it.
- **`BigSpecial.{gamma_big,loggamma_big}`** are the single reflection-safe Γ/logΓ bodies
  (`Perturbation`'s private `_lg`/`_g`/`_needs_refl` are gone), and
  **`Spacetime.Geometry`** the single Kerr Δ / r± / Ω_H / K(r) bodies (the `:GSN` backend's
  `_gsn_*` one-liners forward). `Perturbation.Radial`'s `_sn_hermite` uses
  `NumericsCore.cubic_hermite`.

### Changed — BREAKING: export-name collisions removed

- `PostNewtonian.MMAInputForm.read_mma_inputform` → **`read_mma_inputform_symbolic`**
  (`Data.read_mma_inputform`, the numeric `.wl` reader, keeps the plain name).
- `PostNewtonian.MST.imagunit` → **`imagunit_sym`** (`NumericsCore.imagunit`, the numeric
  imaginary unit, keeps the plain name).
- `SelfForce.SecondOrder.read_h1` → **`read_h1_corpus`** (`Data.read_h1`, the one-mode
  `H1Schema` reader, keeps the plain name).
- `SelfForce.SecondOrder.validate_schema` is now a set of METHODS on
  `NumericsCore.validate_schema` rather than a second generic of the same name.

## [0.1.0] — 2026-08-15

First public release: the next-generation, capability-grouped successor to a
family of single-purpose black-hole-perturbation research codes.

### Added

- **Umbrella package `GrTools`** — `using GrTools` loads the whole suite and
  exposes each package as `GrTools.<Name>`.
- **Core tier (7 packages)** — `NumericsCore` (LinearAlgebra-only numerical
  spine), `BigSpecial` (Complex{BigFloat} Γ/ψ/₂F₁/₁F₁), `Wigner` (BigInt-exact
  3j/6j/CG/Wigner-d), `Data` (data root, artifacts, format readers), `Kernels`
  (second-order Ricci kernels, GPU-capable mode-sum), `EllipticCore` (2D hp-DG
  elliptic solver stack), `TestKit` (test-only tolerance oracle).
- **Capability tier (11 packages)** — `Spacetime`, `Harmonics`, `Perturbation`,
  `MetricReconstruction`, `SelfForce`, `ScalarSelfForce`, `PostNewtonian`,
  `ModifiedGravity`, `EMRIWaveforms`, `Hybrid`, `AccretionDisc`.
- End-to-end BigFloat precision inheritance on every numerical path.
- Characterization test suites (~97 000 assertions) grounded on baked
  predecessor-code fixtures, published reference values, and exact identities;
  dependency budgets enforced as tests.
- Documenter site (page per package + topic pages), seven runnable examples,
  CI matrix over all 19 test targets.

### Known limitations

- The `ScalarSelfForce` effective-source / worldtube / solve / extraction
  stages are in development and not part of this release. (`EllipticCore`'s
  Schwarz complex-shift policy layer was listed here too; it has since landed
  — see [Unreleased].)
- Fixture-regeneration tests (EffectiveSource C oracle, FEW `.fewdata` corpus,
  QQNM source tables) require external mirrors and skip when those are absent;
  the committed fixtures remain authoritative.
- GPU paths in `Kernels` are weak-dependency extensions and are not exercised
  in CI.

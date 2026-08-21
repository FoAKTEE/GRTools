# Changelog

All notable changes to GrTools are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/). Every package in the suite shares
the suite version.

## [Unreleased]

### Added

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

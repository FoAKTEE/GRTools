# Changelog

All notable changes to GrTools are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/). Every package in the suite shares
the suite version.

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

- `EllipticCore`'s Schwarz complex-shift policy layer and the
  `ScalarSelfForce` effective-source / worldtube / solve / extraction stages
  are in development and not part of this release.
- Fixture-regeneration tests (EffectiveSource C oracle, FEW `.fewdata` corpus,
  QQNM source tables) require external mirrors and skip when those are absent;
  the committed fixtures remain authoritative.
- GPU paths in `Kernels` are weak-dependency extensions and are not exercised
  in CI.

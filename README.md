# GrTools.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-informational.svg)](CHANGELOG.md)

**A next-generation Julia toolkit for black-hole perturbation theory and
extreme-mass-ratio-inspiral (EMRI) science.**

GrTools consolidates the calculations that this field has historically spread
across many single-purpose research codes — Teukolsky solvers, harmonic
libraries, geodesic integrators, self-force pipelines, waveform generators —
into one coherently engineered suite: 11 physics **capability** packages
layered over 7 physics-free **core** packages, strictly acyclic. It covers
Kerr geometry and geodesics, spin-weighted spherical/spheroidal harmonics, the
Teukolsky equation (radial solutions, point-particle fluxes, scattering,
quasinormal and quadratic-quasinormal modes, time domain), metric
reconstruction, gravitational and scalar self-force (first and second order),
post-Newtonian expansions, EMRI trajectories and waveforms, hybrid inspiral
models, and relativistic accretion-disc modes.

## What makes it next-generation

- **One `using` per concept.** Packages are organized by physics capability,
  not by historical codebase. "I need the Teukolsky radial solution" means
  `using Perturbation` — nothing else to remember. Every function from the
  legacy generation of codes has exactly one home here (see the Migration
  appendix in the docs).
- **Precision is inherited, never truncated.** Every numerical path is
  type-generic in the element type: BigFloat inputs inherit the active
  `setprecision(BigFloat, bits)` (256–2048 bits in production use) end-to-end.
  Deliberate Float64 fast tiers are documented and route BigFloat callers to a
  reference path rather than silently downcasting.
- **Correctness is characterization-grounded.** The suite was built
  test-first against the implementations it supersedes: every numerical
  module carries characterization gates versus baked legacy-oracle fixtures
  (with recorded capture provenance), published reference tables, or exact
  identities — over 100,000 assertions across the 18 package suites.
- **Dependency hygiene is enforced, not aspirational.** Core packages declare
  binding dependency budgets (NumericsCore: LinearAlgebra only) and several
  packages gate their allowed dependency set in CI (`test/hygiene.jl`). Heavy
  stacks (Symbolics, HDF5, GPU) are quarantined to the packages that need
  them or to weak extensions.
- **Modern solver technology.** In-tree hp-adaptive discontinuous-Galerkin
  elliptic stack (matrix-free operators, complex-generic GMRES,
  Schwarz/multigrid preconditioning), MST/Sasaki–Nakamura/GSN radial
  backends, Leaver and fast QNM continued-fraction solvers, and a
  quadratic-QNM (second-order ringdown) pipeline.

## Installation

Clone and instantiate — packages reference each other with relative-path
`[sources]` entries, so no registry setup is needed:

```bash
git clone <this-repo> GrTools.jl
cd GrTools.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then either load the whole suite through the umbrella package:

```julia
using GrTools
GrTools.Spacetime.kerr_isco(0.9)
```

or work in a single package's own project (each stands alone):

```bash
julia --project=capabilities/Perturbation -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Quickstart

```julia
using GrTools
using GrTools.Spacetime, GrTools.Perturbation, GrTools.Harmonics

# Kerr bound geodesic: constants of motion and fundamental frequencies
com  = constants_of_motion(0.9, 8.0, 0.3, 0.8)      # (a, p, e, x) -> E, Lz, Q
freq = kerr_frequencies(0.9, 8.0, 0.3, 0.8)         # Omega_r, Omega_th, Omega_ph

# Gravitational-wave flux from a circular orbit (Teukolsky, s = -2, l = m = 2)
mode = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.5, 8.0)
mode.EnergyFlux.I                                    # energy flux to infinity

# Kerr quasinormal-mode frequency (Leaver continued fraction)
qnm_frequency(-2, 2, 2, 0, 0.7)                      # complex M*omega at chi = 0.7

# Spheroidal eigenvalue at 256-bit precision — same code path, higher T
setprecision(BigFloat, 256) do
    SpinWeightedSpheroidalEigenvalue(-2, 2, 2, big"0.35")
end
```

Runnable, self-checking scripts for all of this live in [`examples/`](examples/).

## The packages

**Capabilities** (`capabilities/`, physics):

| package | scope |
|---|---|
| `Spacetime` | Kerr background geometry and orbits: metric, tetrads, geodesics, frequencies, separatrix/ISCO, spinning secondary |
| `Harmonics` | spin-weighted spherical + spheroidal harmonics (real and complex spheroidicity), sphere collocation |
| `Perturbation` | the Teukolsky engine: radial solutions (MST/Sasaki–Nakamura/GSN), point-particle fluxes, scattering, spin corrections, QNM/QQNM, time-domain DG |
| `MetricReconstruction` | h_μν from Weyl scalars: CCK radiation gauges, Kerr–Lorenz circular, Li–Yunes transport |
| `SelfForce` | first- + second-order gravitational self-force: Lorenz-gauge solver, punctures/effective source, second-order Ricci source |
| `ScalarSelfForce` | m-mode elliptic Kerr scalar self-force (arXiv:2606.04998 reproduction) over the EllipticCore stack |
| `PostNewtonian` | PN-series engine (Symbolics quarantine): PN algebra, symbolic MST, PN self-force, flux catalogues |
| `ModifiedGravity` | beyond-GR (modified Teukolsky) source machinery |
| `EMRIWaveforms` | EMRI trajectory + observed waveform: inspiral ODE, flux/amplitude models, AAK, mode-sum |
| `Hybrid` | self-force–PN hybrid inspiral/waveform (WaSABI) |
| `AccretionDisc` | relativistic disc modes + disc→EMRI back-reaction |

**Core** (`core/`, physics-free infrastructure):

| package | scope |
|---|---|
| `NumericsCore` | the numerical spine (depends on LinearAlgebra only): special functions, elliptic integrals, interpolation, root-finding, vendored ODE integrators, DG spectral-element core |
| `BigSpecial` | `Complex{BigFloat}`-correct Γ/ψ/₂F₁/₁F₁ |
| `Wigner` | 3j/6j/Clebsch–Gordan/Wigner-d, BigInt-exact core |
| `Data` | data root + Pkg artifacts + format readers (HDF5/JLD2 as weak extensions) |
| `Kernels` | second-order Ricci quadratic kernels + GPU-capable mode-sum (CUDA/AMDGPU weak deps) |
| `EllipticCore` | complex-parametric 2D hp-DG elliptic solver stack: reference elements, multi-block meshes, matrix-free DG, GMRES/Schwarz/multigrid, AMR |
| `TestKit` | test-only tolerance oracle + characterization gates |

The dependency DAG is strictly layered and acyclic — see the
[docs index](docs/src/index.md) for the full tier diagram.

## Documentation

One Documenter site covers the whole suite (page per package + cross-cutting
topic pages on conventions, precision, GPU/parallel, and data artifacts):

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl        # output in docs/build/
```

## Testing

Every package carries its own suite and runs standalone:

```bash
julia --project=capabilities/Perturbation -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Tests are grounded: numerics are characterization-gated against stated oracles
(published values, baked upstream fixtures with capture provenance, exact
identities), tolerances go through TestKit's frozen tolerance classes, and
several packages enforce their dependency budget as a test (`test/hygiene.jl`).
CI runs the full 19-suite matrix (18 packages + umbrella) — see
`.github/workflows/CI.yml`.

## Conventions

Geometric units `G = c = 1`, black-hole mass `M = 1` unless stated.
Boyer–Lindquist coordinates; Kinnersley tetrad; fields `∝ e^{i(mφ − ωt)}`.
Spin weight `s ∈ {0, ±1, ±2}`, default −2. Spheroidicity `γ = aω`; Teukolsky
eigenvalue convention `λ = E − 2mγ + γ²` with `λ(0) = l(l+1) − s(s+1)`. The
full statement (including the `m < 0` parity surface and ε↔ω maps) is in
`docs/src/topics/conventions.md`.

## Status

GrTools is at **v0.1.0** — the first public release. The public API of every
package is exercised by its test suite and documented on its docs page; the
`Known limitations` section of [CHANGELOG.md](CHANGELOG.md) lists the
components still in development. Pre-1.0, minor versions may include breaking
changes, which will always be called out in the changelog.

## Contributing and citing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the layout, design rules, and PR
checklist; [CHANGELOG.md](CHANGELOG.md) for release history; and
[CITATION.cff](CITATION.cff) to cite the suite. Licensed under
the [MIT License](LICENSE).

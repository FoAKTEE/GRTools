# GrTools

**GrTools** is a next-generation Julia toolkit for black-hole perturbation
theory and EMRI (extreme-mass-ratio inspiral) science. It consolidates the
calculations this field has historically spread across many single-purpose
research codes into one coherently engineered, capability-grouped monorepo:
11 physics **capability** packages layered over 7 **core/infra** packages,
strictly acyclic, living under `core/` + `capabilities/`. The umbrella package
`GrTools` (repo root) loads the whole suite: `using GrTools` exposes each
package as `GrTools.<Name>`.

## The thesis: one `using` per mental concept

Each package is one concept in the physicist's mental map. If you are thinking
"I need the Teukolsky radial solution", you type `using Perturbation` — not the
name of whichever historical repository the function happened to live in. The
20 predecessor research codes were organized by *provenance*; GrTools is
organized by *capability*. Every function from the predecessor codes has
exactly one new home — see the [Migration](migration.md) appendix to find
yours.

## The 11 capabilities

| # | Capability | Scope (one line) |
|---|---|---|
| 1 | **[Spacetime](capabilities/spacetime.md)** | Kerr background-geometry + orbit foundation (spine root): geometry, geodesics, separatrix/ISCO, spinning secondary, numeric PN kinematics |
| 2 | **[Harmonics](capabilities/harmonics.md)** | Angular-eigenfunction foundation (2nd spine root): SWSH eigenvalues + S + spherical Y + parity unifier; BigFloat eigenvalue inheritance |
| 3 | **[Perturbation](capabilities/perturbation.md)** | Vacuum/test perturbed-field solver (Teukolsky radial FD backends incl. :GSN + the TimeDomain DG solver) + amplitudes/fluxes/scattering/spin-fluxes — the light-dep firewall |
| 4 | **[MetricReconstruction](capabilities/metricreconstruction.md)** | First-order metric reconstruction from the curvature scalar (HDF5-free): CCK/IRG-ORG radiation gauge + Dolan-2024 Kerr-Lorenz circular |
| 5 | **[SelfForce](capabilities/selfforce.md)** | 1st+2nd-order GSF corpus — HDF5 home #1: Lorenz field solver+writer, punctures/effective-source, 2nd-order Ricci source |
| 6 | **[PostNewtonian](capabilities/postnewtonian.md)** | PN-series engine — the Symbolics quarantine: PN-series algebra, symbolic MST, PN self-force, MMA catalog |
| 7 | **[ModifiedGravity](capabilities/modifiedgravity.md)** | Beyond-GR (modified Teukolsky) source machinery — thin research surface |
| 8 | **[EMRIWaveforms](capabilities/emriwaveforms.md)** | EMRI trajectory + observed waveform — HDF5 home #2: inspiral ODE, flux/amp models, AAK, mode-sum, generate |
| 9 | **[Hybrid](capabilities/hybrid.md)** | SF-PN hybrid inspiral/waveform (WaSABI), outside the Symbolics quarantine |
| 10 | **[AccretionDisc](capabilities/accretiondisc.md)** | Disc modes + disc→EMRI back-reaction — top of DAG: Modes (HDF5-free) + Force (HDF5) |
| 11 | **[ScalarSelfForce](capabilities/scalarselfforce.md)** | m-mode elliptic Kerr scalar self-force (arXiv:2606.04998 reproduction): mode ansatz → flux form → puncture → worldtube over the EllipticCore hp-DG stack |

## The 7 core/infra packages

| Unit | Package | One line |
|---|---|---|
| C0 | **[NumericsCore](core/numericscore.md)** | Dimensionless math sink (LinearAlgebra only): Carlson/elliptic/Jacobi, interpolation lib, Legendre seed, mode-index, precision policy, BLAS-guard, vendored DP5/GBS/DOPR853, DG spectral-element core, continued-fraction, root-finders, dep-free DataSchema structs |
| C1 | **[BigSpecial](core/bigspecial.md)** | `Complex{BigFloat}` Γ/ψ/polygamma/pochhammer + self-contained 2F1/1F1 |
| C2 | **[Wigner](core/wigner.md)** | Hand-rolled 3j/6j/Racah-CG, Wigner-D (general + β=π/2), 3-harmonic coupling |
| C3 | **[Data](core/data.md)** | Versioned h1schema contract + deduped grid loaders + MMA/.dat readers; bodies in weakdep-ext, schema structs dep-free |
| C4 | **[Kernels](core/kernels.md)** | TD+FD-SPA+AAK mode-sum, r2_1..r2_10, ROMAN GEMM (KernelAbstractions; CUDA/AMDGPU weakdep) |
| C5 | **[TestKit](core/testkit.md)** | `@approx_oracle` frozen-tolerance harness + characterization gates (test-only) |
| C6 | **[EllipticCore](core/ellipticcore.md)** | Complex-parametric 2D hp-DG elliptic solver stack (zero physics): reference elements, multi-block meshes, matrix-free DG, GMRES/Schwarz/multigrid, AMR |

## Dependency DAG (tiers; deps point DOWN, strictly layered, acyclic)

```
TIER 0  [0] NumericsCore (LinearAlgebra only — unique sink)

TIER 1  [1] BigSpecial → 0           [2] Wigner → 0
        [3] Kernels → 0              [4] Data → 0          [5] TestKit → 4 (test-only)
        [16] EllipticCore → 0 (+SparseArrays)

TIER 2  [6] Spacetime → 0            [7] Harmonics → 0
        (spine roots: zero heavy deps)

TIER 3  [8] Perturbation → 6,7,1,0
        (vacuum engine; HDF5-free, Symbolics-free, external-ODE-free)

TIER 4  [9]  MetricReconstruction → 8,7,2,6,0 (+QuadGK)
        [10] PostNewtonian → 7,0,4 + Symbolics (quarantined; nothing depends UP)
        [11] ModifiedGravity → 8,7,0 (+Printf, Data opt)

TIER 5  [12] SelfForce → 9,7,1,2,3,4,0 + HDF5
        [13] EMRIWaveforms → 6,7,3,4,0 + HDF5, OrdinaryDiffEqVerner, Roots, NLsolve
        [14] Hybrid → 6,1,4,0 + OrdinaryDiffEqVerner, SciMLBase (no HDF5; JLD2 direct)

TIER 6  [15] AccretionDisc → 9,8,7,6,4,0 + HDF5 (Force submodule only)

TIER 3' [17] ScalarSelfForce → 16,6,1,4,0 (+Printf) — elliptic m-mode branch
```

## How to navigate this site

- **Capabilities** — one page per physics capability: submodules, public API,
  dependencies (with the heavy-dep declaration), GPU/data needs, and
  provenance.
- **Core** — reference pages for the 7 infra packages every capability builds on.
- **Topics** — cross-cutting, single-source pages:
  [Conventions](topics/conventions.md) (the `m<0` parity surface, ε↔ω, units),
  [Precision](topics/precision.md) (BigFloat inheritance +
  `angular_precision_floor()`), [GPU/Parallel](topics/gpu-parallel.md)
  (heavy-run gating, lease model, parallel==serial), and
  [Data artifacts](topics/data-artifacts.md) (h1schema contract, datadir injection).
- **[Migration](migration.md)** — the lookup table from the predecessor codes
  to GrTools: look up the old package/function, find its new capability home.

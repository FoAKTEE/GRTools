# Harmonics — capability (Tier-2 spine root)

The angular-eigenfunction **foundation** of GrTools (design §4.2). It owns every angular
harmonic object of the Kerr perturbation problem and exposes them through one consolidated,
type-generic surface:

- spin-weighted **spheroidal eigenvalues** λ(s,l,m,γ) + ∂_γλ,
- the spheroidal harmonic **S_lm(θ)** carrier and its θ/φ/γ derivatives,
- the spin-weighted **spherical-harmonic Y** seed (γ=0), and
- the **unified m<0 parity surface** that collapses the per-package forks into one tested flag.

It is the second spine root: heavy-dep-FREE by construction, so the light numeric path
(`Perturbation`, `MetricReconstruction`, …) inherits no heavy dependency through it.

## Physics this capability owns

- **`SpheroidalEigenvalues`** — λ(s,l,m,γ) + ∂_γλ behind ONE dispatcher over three backends:
  pentadiagonal-LAPACK spectral Float64 fast path · type-generic Leaver continued-fraction
  BigFloat path · small-γ Taylor + large-γ Laurent asymptotics. This is the **BigFloat
  eigenvalue home** (H-CORE-STRAT).
- **`SpheroidalHarmonics`** — the callable `S_lm(θ)` carrier (spherical-expansion + Leaver
  representations) with mixed ∂_θ/∂_φ/∂_γ; the Float64-VALUE tier.
- **`SphericalHarmonics`** — the γ=0 canonical spin-weighted Y seed (via
  `NumericsCore.legendre`), a half-angle BigFloat high-l path, exact θ/θθ derivatives, and ONE
  `SpinWeightedSphericalHarmonicY(s,l,m,θ,φ; parity=:none|:few)` — `:none` by default, `:few`
  applies the `(-1)^l` flip that FEW/KerrBackwards need. One tested flag, never a silent
  per-call fork.
- **`Precision`** — `angular_precision_floor()`, surfacing the Float64 value-ceiling vs the
  BigFloat eigenvalue-floor.

All primitives are type-generic `{T}`, so BigFloat precision (256/512/1024/2048 bits) is
inherited. The spectral LAPACK eigensolve is wrapped in the BLAS-guard.

## Legacy packages absorbed (via REWRITE/promotion, not copy)

| Legacy source | Function(s) | New home |
|---|---|---|
| `SpinWeightedSpheroidalHarmonics/spectral.jl` | `SpinWeightedSpheroidalEigenvalue(+Deriv)`, `_swsh_eigenvalue_spectral` | `Harmonics.SpheroidalEigenvalues` |
| `SpinWeightedSpheroidalHarmonics/leaver.jl` | `_swsh_eigenvalue_leaver`, `_swsh_leaver_coeffs` | `Harmonics.SpheroidalEigenvalues` |
| `SpinWeightedSpheroidalHarmonics/series_coeffs.jl` | `EigenvalueSeries(+LargeGamma)`, `_swsh_alpha` | `Harmonics.SpheroidalEigenvalues` |
| `SpinWeightedSpheroidalHarmonics/harmonic.jl,leaver.jl` | `SpinWeightedSpheroidalHarmonicS`, `SWSHSFunction`, `SWSHSLeaverFunction`, `swsh_derivative` | `Harmonics.SpheroidalHarmonics` |
| `SpinWeightedSpheroidalHarmonics/spherical_Y.jl` | `SpinWeightedSphericalHarmonicY`, `_y_terms`, `_y_dtheta` | `Harmonics.SphericalHarmonics` |
| `FastEMRIWaveforms/harmonics.jl` | `ylm` | `Harmonics.SphericalHarmonics` (`parity=:few` shim) |
| `KerrBackwards/ylm.jl` | `ylm_few`, `ylms_few` | `Harmonics.SphericalHarmonics` (`parity=:few` shim) |
| `DiscModes/metric_recon.jl` | `Ylm_real`, `Ylm_eq` | `Harmonics.SphericalHarmonics` (real-Y path) |
| `KerrLorenzCirc/precision.jl` | `swsh_eigenvalue_big` | `Harmonics.SpheroidalEigenvalues` (fold candidate — confirm signature) |
| NEW | `spin_weighted_spherical_Y`, `angular_parity`, `angular_precision_floor` | `Harmonics` |

The legacy packages stay green as the characterization oracle; they are never edited.

## Public API (built TDD in the BUILD rounds)

Spheroidal eigenvalues + derivative, `S_lm(θ)` carrier + `swsh_derivative`,
`SpinWeightedSphericalHarmonicY(s,l,m,θ,φ; parity)`, `angular_parity`,
`angular_precision_floor`.

## Dependencies

- Runtime: `LinearAlgebra` (spectral eigensolve) + `NumericsCore` (precision, mode-index,
  Legendre primitive, continued-fraction primitive, secant, self-contained Kummer 1F1) +
  `BigSpecial` + `Wigner`. **No heavy deps** (no HDF5/Symbolics/ODE/WignerSymbols package);
  no local Kummer/CF copy.
- Test-only oracle: `SpinWeightedSpheroidalHarmonics` (LIGHT, deps = LinearAlgebra only) in
  `[extras]`/`[targets]`, grounding the eigenvalue solver, `S_lm(θ)` + θ-derivative, and the
  spin-weighted-Y seed at the legacy tolerances.

## Grounding

Test-only characterization against the `SpinWeightedSpheroidalHarmonics` oracle at its
tolerances, plus analytic identities: Y orthonormality, l=|m| closed forms, the γ→0
spheroidal→spherical limit, the eigenvalue `l(l+1) − s(s+1)` at γ=0, and the
Condon-Shortley/Goldberg m↔−m relations.

# Harmonics

The angular-eigenfunction **foundation** of GrTools (design §4.2) and the second Tier-2
spine root. It owns every angular harmonic object of the Kerr perturbation problem behind
one consolidated, type-generic surface:

- spin-weighted **spheroidal eigenvalues** ``\lambda(s,\ell,m,\gamma)`` and
  ``\partial_\gamma\lambda``,
- the callable spheroidal harmonic ``{}_sS_{\ell m}(\theta,\varphi)`` carrier with mixed
  ``\partial_\theta/\partial_\varphi/\partial_\gamma`` derivatives,
- the spin-weighted **spherical-harmonic** ``{}_sY_{\ell m}`` seed (the ``\gamma=0`` limit), and
- the **unified `m<0` parity surface** that collapses the per-package sign forks
  (FEW, KerrBackwards, DiscModes real-Y) into one tested flag.

Every primitive is type-generic `{T}`, so BigFloat precision (256/512/1024/2048 bits) is
inherited by the eigenvalue/Leaver paths; the spectral LAPACK eigensolve is wrapped in the
BLAS-guard. The harmonic *value* tier is deliberately Float64-capped (see
[Precision policy](@ref harmonics-precision-policy) below).

## Submodule map

`Harmonics` is a **single Julia module** — the units below are logical (one source file
each), not nested modules:

| Unit (file) | Owns |
|---|---|
| `SphericalHarmonics.jl` | ``\gamma=0`` spin-weighted ``Y`` front door (`SpinWeightedSphericalHarmonicY`, `spin_weighted_spherical_Y`), the parity unifier `angular_parity`, and the legacy convention shims `ylm`/`ylm_few`/`ylms_few`. Carries **no** Legendre/half-angle copy — the seed (`spin_weighted_Y`, `Ylm_real`, `Ylm_eq`) lives in `NumericsCore` and is re-exported here so callers depend on Harmonics, not NumericsCore directly. |
| `SpheroidalEigenvalues.jl` | ``\lambda(s,\ell,m,\gamma)`` + ``\partial_\gamma\lambda`` behind **one dispatcher** over the backends (`SpinWeightedSpheroidalEigenvalue`, `SpinWeightedSpheroidalEigenvalueDeriv`). The BigFloat eigenvalue home (H-CORE-STRAT). |
| `SpheroidalEigenvalueSeries.jl` | Asymptotic-series backend: small-``\gamma`` Taylor coefficients (`SpinWeightedSpheroidalEigenvalueSeries`, Stewart/Press–Teukolsky) and large-``\gamma`` Laurent coefficients (`SpinWeightedSpheroidalEigenvalueSeriesLargeGamma`, Casals–Ottewill–Warburton, arXiv:1810.00432). |
| `SpheroidalHarmonics.jl` | Callable ``{}_sS_{\ell m}`` carriers in SphericalExpansion (`SWSHSFunction`) and Leaver (`SWSHSLeaverFunction`) representations, `SpinWeightedSpheroidalHarmonicS`, the small-``\gamma`` `SpinWeightedSpheroidalHarmonicSSeries` carrier, `expansion_coefficients`, and the derivative surface `swsh_derivative`/`swsh_gamma_derivative`. The Float64-**value** tier. |
| `precision.jl` | `angular_precision_floor` — the precision policy in code. |

### The eigenvalue dispatcher

`SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Auto, ...)` routes:

- `:Auto` — a `BigFloat` ``\gamma`` goes to the type-generic **Leaver** continued-fraction
  backend; anything else to the pentadiagonal-LAPACK **spectral** Float64 fast path
  (`:SphericalExpansion`).
- `:Leaver` — force the continued-fraction path (the BigFloat route; consumes the
  `NumericsCore` continued-fraction and self-contained Kummer ``{}_1F_1`` primitives).
- `:SmallGamma` — evaluate the small-``\gamma`` Taylor series ``\sum_i c_i\,\gamma^i``
  (coefficients from `SpinWeightedSpheroidalEigenvalueSeries`).
- `:LargeGamma` — evaluate the large-``\gamma`` Laurent series
  ``\mathrm{evalpoly}(\gamma, c^{L})/\gamma^{\text{order}}``
  (coefficients from `SpinWeightedSpheroidalEigenvalueSeriesLargeGamma`); the asymptotic
  error decays as ``\gamma\to\infty``.

Any other symbol throws `ArgumentError`. The asymptotic routes are first-class dispatcher
citizens, not a separate API.

## Dependencies and poison-axis declaration

Runtime `[deps]` (Project.toml): **`LinearAlgebra`** (the spectral eigensolve) and
**`NumericsCore`** (precision policy, mode-index, Legendre seed, continued-fraction
primitive, secant, the promoted self-contained Kummer ``{}_1F_1``). That is the *entire*
runtime closure — the §3.2 spine-root pin — so the light numeric path (`Perturbation`,
`MetricReconstruction`, …) inherits **no heavy dependency** through Harmonics.

This budget is CI-enforced by `test/hygiene.jl` (§7.7 gates 1+2):

- **Gate 1 (dep budget):** `allow = [LinearAlgebra, NumericsCore]`;
  `forbid = [Roots, HDF5, JLD2, Symbolics, OrdinaryDiffEq*, CUDA, AMDGPU,
  KernelAbstractions, WignerSymbols, BigSpecial, Wigner, Kernels, Data, SpecialFunctions,
  HypergeometricFunctions, NLsolve, SciMLBase, QuadGK]`.
- **Gate 2 (closure load):** after `using Harmonics`, none of the poison packages may
  appear in the loaded-module set.

In particular Harmonics carries **no local Kummer/continued-fraction copy** (H-CORE-STRAT:
the LinearAlgebra-only primitives were promoted into NumericsCore's light tier;
`BigSpecial`'s HypGeo-backed `kummer_hyp` is a distinct, separately-named function that
never enters this closure), and **Wigner is deliberately not folded in** — co-locating it
would leak angular-coupling deps onto the in-degree spine (§4.2).

## GPU and data needs

None. Harmonics is CPU-only; parallelism is threads/processes over the ``(\ell,m)`` mode
axis. The BigFloat generation tier (eigenvalues, Leaver) never touches a device — GPU work
lives in `Kernels`/`EMRIWaveforms`, which call these type-generic primitives at Float64.
There are no runtime data artifacts; the only files under `test/fixtures/` are baked
characterization references (Phase 7-B), so the test suite runs without any legacy
`packages/*` checkout.

## Legacy provenance (drawsFrom)

Functions landed here by REWRITE/promotion (design §6.2 crosswalk), grounded by
characterization tests against the legacy `SpinWeightedSpheroidalHarmonics` oracle
(now baked into `test/fixtures/*_ref.jl`):

| Legacy source | Function(s) | New home |
|---|---|---|
| `SpinWeightedSpheroidalHarmonics/spectral.jl` | `SpinWeightedSpheroidalEigenvalue(+Deriv)`, `_swsh_eigenvalue_spectral` | `SpheroidalEigenvalues.jl` |
| `SpinWeightedSpheroidalHarmonics/leaver.jl` | `_swsh_eigenvalue_leaver`, `_swsh_leaver_coeffs` | `SpheroidalEigenvalues.jl` |
| `SpinWeightedSpheroidalHarmonics/series_coeffs.jl` | `EigenvalueSeries(+LargeGamma)`, `_swsh_alpha` | `SpheroidalEigenvalueSeries.jl` |
| `SpinWeightedSpheroidalHarmonics/harmonic.jl, leaver.jl` | `SpinWeightedSpheroidalHarmonicS`, `SWSHSFunction`, `SWSHSLeaverFunction`, `swsh_derivative` | `SpheroidalHarmonics.jl` |
| `SpinWeightedSpheroidalHarmonics/spherical_Y.jl` | `SpinWeightedSphericalHarmonicY`, `_y_terms`, `_y_dtheta` | `SphericalHarmonics.jl` (term machinery promoted to `NumericsCore`) |
| `FastEMRIWaveforms/harmonics.jl` | `ylm` | `SphericalHarmonics.jl` (`convention=:few` shim) |
| `KerrBackwards/ylm.jl` | `ylm_few`, `ylms_few` | `SphericalHarmonics.jl` (`convention=:few` shim) |
| `DiscModes/metric_recon.jl` | `Ylm_real`, `Ylm_eq` | re-exported from `NumericsCore` (real-Y path) |
| `KerrLorenzCirc/precision.jl` | `swsh_eigenvalue_big` | fold candidate into `SpheroidalEigenvalues.jl` (signature-gated) |
| NEW | `spin_weighted_spherical_Y`, `angular_parity`, `angular_precision_floor` | this capability |

## [Precision policy](@id harmonics-precision-policy)

`angular_precision_floor(T)` returns the `(value_tier, eigenvalue_tier)` pair that codifies
the load-bearing asymmetry: the eigenvalue/Leaver paths inherit `BigFloat` end-to-end, but
the harmonic **value** carriers (`SWSHSFunction`, `SWSHSLeaverFunction`) store Float64
coefficients — a Float64 value-ceiling above a BigFloat eigenvalue-floor.

## What changed vs legacy

- **One parity flag, not four silent forks.** Legacy code disagreed on the negative-``m``
  sign: FEW `GetYlms` and KerrBackwards `ylm_few` apply ``(-1)^\ell``; the DiscModes real-Y
  path uses ``(-1)^m \cdot \mathrm{conj}`` reality relations; the SWSH package applies
  nothing. These are now one primitive,
  `spin_weighted_spherical_Y(s, l, m, θ, φ; convention=:none|:few|:reality)`, with the sign
  isolated in `angular_parity(l, m; convention)`. `parity=:none` is bit-identical to the
  legacy SWSH implementation; the legacy names survive as thin documented shims.
- **BigFloat eigenvalues.** The Leaver continued-fraction backend is type-generic, removing
  the legacy Float64 caps; `:Auto` routes on the type of ``\gamma``.
- **Asymptotic series are dispatcher routes.** The small-``\gamma`` Taylor and
  large-``\gamma`` Laurent expansions, previously coefficient tables on the side, are
  reachable as `method=:SmallGamma` / `:LargeGamma` on the same public entry point.
- **No local special-function copies.** The Kummer ``{}_1F_1``, continued-fraction,
  Legendre/half-angle-Y, and term-list evaluation machinery were promoted to `NumericsCore`
  and are consumed from there; the restored runtime closure is exactly
  `{LinearAlgebra, NumericsCore}`.
- **Oracle baked.** The `SpinWeightedSpheroidalHarmonics` test oracle was captured once into
  fixtures; the runtime package and its tests carry no legacy dependency.

## Public API

```@autodocs
Modules = [Harmonics]
Order = [:module, :type, :constant, :function]
```

!!! note
    `Ylm_real` and `Ylm_eq` are defined in `NumericsCore` and re-exported by `Harmonics`
    (the DiscModes real-``Y`` path); their docstrings render on the NumericsCore page.

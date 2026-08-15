# BigSpecial (`GrToolsBigSpecial`)

Tier-1 core stratum: the `Complex{BigFloat}` special-functions layer of GrTools.
It exists as a **separate** package so that
`HypergeometricFunctions`/`SpecialFunctions` never leak onto the
`NumericsCore`-only spine (`KerrGeodesics` empty-deps,
`SpinWeightedSpheroidalHarmonics` LinearAlgebra-only). Any caller needing
high-precision ``\Gamma`` / ``\psi`` / ``{}_2F_1`` / ``{}_1F_1`` depends on
BigSpecial *explicitly*; the heavy hypergeometric closure reaches users only
through this package. Downstream consumers: Perturbation (MST/Leaver),
PostNewtonian, Hybrid (WaSABI digamma).

Everything is type-generic `{T}`, so `Complex{BigFloat}` precision is
*inherited*, never hard-`Float64`'d — the property the MST/Leaver expansions
at 256–2048 bits rely on.

## Physicist's map

Three families, one file each:

- **Gamma family** (`src/biggamma.jl`, `src/digamma.jl`, `src/pochhammer.jl`) —
  `gamma_big` / `loggamma_big` are valid for **all** ``z``: for
  ``\operatorname{Re} z < 1/2`` they route through the Euler reflection
  ``\Gamma(z) = \pi / (\sin(\pi z)\,\Gamma(1-z))`` so the inner ``\Gamma`` only
  ever sees the half-plane where `SpecialFunctions` is trustworthy at BigFloat
  precision (`needs_refl` is the documentary predicate for the failing argument
  set). `fact_big(n) = Γ(n+1)` and the exact integer-product `factorial_big`;
  `digamma_big` is the asymptotic-series ``\psi(z)`` (recurrence-shifted to
  ``\operatorname{Re} z \ge 12``), and `PolyGamma(n, z)` is the WaSABI-named
  entry point (currently `n = 0` only — higher orders throw `ArgumentError`,
  preserving the legacy contract). `pochhammer(a, n)` is the dependency-free
  rising factorial ``(a)_n`` (exact products, no ``\Gamma`` — valid for
  negative `n` too).
- **Gauss ``{}_2F_1``** (`src/hyp2f1_primitive.jl`) — self-contained:
  `f21_maclaurin` (incremental-term series, ``|z|<1``), `f21_zinf` (the
  ``z\to\infty`` connection, DLMF 15.8.2 — two series in ``1/z``),
  `f21_terminating` (the polynomial case ``a = -l``), and `safe_2f1` — the
  robust dispatcher that calls `HypergeometricFunctions._₂F₁` (byte-identical
  for `Float64`/`ComplexF64`) and, when upstream throws on a
  `Complex{BigFloat}` argument, falls back to the Pfaff transform
  ``z \to w = z/(z-1)`` plus the ``\Gamma``/``\psi``-free Maclaurin series.
- **Confluent ``{}_1F_1``** (`src/kummer.jl`) — `kummer_hyp(a, b, z)`, the
  HypergeometricFunctions-**delegated** Kummer ``M`` (via `pFq`). This is a
  deliberately **distinct named function** from the series-self-contained
  `NumericsCore.kummer_1F1`: HypGeo-delegated vs. vendored-series are different
  decisions, so they are two named primitives, not one unacknowledged fork.

## Module map

Single flat module — no submodules:

| File | Exports |
|---|---|
| `src/pochhammer.jl` | `pochhammer` |
| `src/biggamma.jl` | `gamma_big`, `loggamma_big`, `needs_refl`, `fact_big`, `factorial_big` |
| `src/digamma.jl` | `digamma_big`, `PolyGamma` |
| `src/hyp2f1_primitive.jl` | `f21_maclaurin`, `f21_zinf`, `safe_2f1`, `f21_terminating` |
| `src/kummer.jl` | `kummer_hyp` |

## Dependencies and poison-axis budget

Hard deps (`Project.toml`): `NumericsCore` (precision policy, conventions) plus
the two sanctioned external backends `HypergeometricFunctions` and
`SpecialFunctions`. BigSpecial is **the only place** the heavy
special-function stack is allowed below tier 3.

The CI dependency budget (`test/hygiene.jl`) is binding:

- **allow**: `NumericsCore`, `HypergeometricFunctions`, `SpecialFunctions`;
- **forbid** (dep-budget gate): `HDF5`, `JLD2`, `Symbolics`,
  `OrdinaryDiffEq*`, `CUDA`, `AMDGPU`, `KernelAbstractions`, `WignerSymbols`,
  `Kernels`, `Data`, `Roots`, `NLsolve`, `SciMLBase`, `QuadGK`;
- **closure-load gate**: after `using BigSpecial`, none of the forbidden
  poison axes may appear in `Base.loaded_modules`
  (`SpecialFunctions`/`HypergeometricFunctions` are legitimately loaded and
  excluded from that list).

## GPU / data needs

None. Pure CPU, no data artifacts, no I/O: the BigFloat generation tier never
touches the device, and the budget above forbids every GPU and
I/O axis. Test-oracle values are **baked fixture literals** in `test/*.jl`
(captured once at 256-bit, bit-identical `==` preserved), so the test suite
needs no dependency on the predecessor codes and no data files.

## Provenance

Consolidates four predecessor forks of the same decision into one canonical
type-generic primitive per function:

| Legacy (pkg / file / fn) | → BigSpecial |
|---|---|
| Teukolsky / `hypergeom.jl`, `amplitudes.jl` / `_f21_maclaurin`, `_f21_zinf`, `_gamma_big`, `_pochhammer` | `f21_maclaurin`, `f21_zinf`, `gamma_big`, `pochhammer` |
| H1Lorenz / `math_funcs.jl` / `_fact` | `fact_big` / `factorial_big` |
| Punctures / `alpha_integration.jl`, `regsum.jl` / `_poch`, `_f21_term`, `_regsum_tail` | `pochhammer`, `f21_terminating`, `digamma_big` |
| WaSABI / inspiral deps / hand-rolled `_digamma`, `PolyGamma` | `digamma_big`, `PolyGamma` |

## What changed from the predecessor codes

- **Fork collapse.** Four per-package copies of ``\Gamma``/pochhammer/``{}_2F_1``
  machinery became one type-generic canonical primitive each; every move was
  characterization-gated against the legacy oracle at the legacy tolerance.
- **AHA exception kept local.** Teukolsky's `safe_2f1` ``z\to\infty``
  branch-*selection* heuristic is a per-caller decision and was **not** pulled
  here — only the dispatcher + Pfaff fallback live in BigSpecial; the fast
  large-``|z|`` route selection stays with the caller.
- **Two Kummers, on purpose.** `kummer_hyp` (HypGeo-delegated, here) and
  `NumericsCore.kummer_1F1` (self-contained series, LinearAlgebra-only tier)
  are distinct named functions, not a fork — the split is what keeps the
  Harmonics spine heavy-dep-free.
- **`PolyGamma` narrowed honestly.** Only `n = 0` (digamma) is implemented;
  other orders throw rather than silently degrade — the legacy WaSABI call
  sites only used `n = 0`.
- **Tests de-legacied.** Characterization oracle values (Teukolsky MST/Leaver
  ``{}_2F_1`` + `gamma_big`/`pochhammer`) are baked into fixture literals; the
  suite runs with zero live legacy deps.

## Public API

```@autodocs
Modules = [BigSpecial]
Order = [:module, :type, :constant, :function]
```

# BigSpecial

Tier-1 core stratum: the `Complex{BigFloat}` special-functions layer for GrTools.
Kept SEPARATE from `NumericsCore` so that `HypergeometricFunctions`
and `SpecialFunctions` never leak onto the NumericsCore-only spine
(`KerrGeodesics` empty-deps, `SpinWeightedSpheroidalHarmonics` LinearAlgebra-only).
Any caller needing high-precision Γ / ψ / polygamma / ₂F₁ / ₁F₁ depends on
BigSpecial *explicitly*; the heavy hypergeometric closure reaches users only here.

## Physics / numerics owned

- **`biggamma.jl`** — `gamma_big`, `loggamma_big`, `digamma_big`, `polygamma_big`,
  `pochhammer`, `factorial_big` / `fact_big`. One type-generic canonical primitive
  consolidating four legacy forks of the same decision. Exact nonpositive-integer Γ/logΓ poles
  use a documented direction-independent `+Inf` (`+Inf + 0im` for complex/logΓ) convention.
- **`hyp2f1_primitive.jl`** — self-contained Gauss ₂F₁: Maclaurin series, the z→∞
  connection (DLMF 15.8.2), a safe dispatcher, and the terminating (polynomial) case.
- **`kummer.jl`** — `kummer_hyp`, the HypGeo-delegated confluent ₁F₁. A **distinct**
  named function from `NumericsCore.kummer_1f1` (series-self-contained): HypGeo-delegated
  vs series-self-contained are different decisions (AHA), so two named functions —
  not one unacknowledged fork.

## Provenance

| Legacy (pkg / file / fn) | → BigSpecial |
|---|---|
| Teukolsky / hypergeom.jl, amplitudes.jl / `_f21_maclaurin`, `_f21_zinf`, `_gamma_big`, `_pochhammer` | `hyp2f1`, `gamma_big`, `pochhammer` |
| H1Lorenz / math_funcs.jl / `_fact` | `fact_big` / `factorial_big` |
| Punctures / alpha_integration.jl, regsum.jl / `_poch`, `_f21_term`, `_regsum_tail` | `pochhammer`, `hyp2f1`, `digamma_big` |
| WaSABI / inspiral deps / hand-rolled `_digamma`, `PolyGamma` | `digamma_big`, `polygamma_big` |

**AHA exception:** Teukolsky's `safe_2f1` z→∞ branch-*selection* heuristic is a
per-caller decision and stays LOCAL to Teukolsky — it is not pulled here.

## Dependencies

- `NumericsCore` (path `../NumericsCore`) — precision policy, conventions.
- `HypergeometricFunctions`, `SpecialFunctions` — external hypergeometric/Γ backends.

Characterization oracle: the predecessor `Teukolsky` code served as the oracle during
construction, grounding the MST/Leaver ₂F₁ + Kummer + `gamma_big`/`pochhammer` primitives at
its tolerances; its outputs are baked as literal fixtures in the test files, so it is not a
test dependency.

## Status

Scaffold. Module loads; each primitive is ported tests-first,
with expected values grounded by analytic identity or characterization vs the predecessor
oracle, never by the new code's own output.

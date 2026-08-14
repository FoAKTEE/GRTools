"""
    BigSpecial (GrToolBigSpecial)

The `Complex{BigFloat}` special-functions stratum for GrTool (design §5.2). It is a
SEPARATE tier so that `HypergeometricFunctions`/`SpecialFunctions` never leak onto the
`NumericsCore`-only spine (KerrGeodesics empty-deps, SWSH LinearAlgebra-only); anything
needing high-precision Γ/ψ/₂F₁/₁F₁ depends on BigSpecial explicitly.

Owns (per §5.2 + §6 crosswalk rows dest=BigSpecial):
  - `biggamma.jl`  — gamma_big/loggamma/digamma/polygamma/pochhammer + factorial_big/fact_big
    (consolidates Teukolsky `_gamma_big`/`_pochhammer`, H1Lorenz `_fact`, Punctures `_poch`,
    WaSABI hand-rolled digamma/PolyGamma into one type-generic canonical primitive).
  - `hyp2f1_primitive.jl` — self-contained ₂F₁: Maclaurin / z→∞ (DLMF 15.8.2) / safe dispatcher /
    terminating (consolidates Teukolsky `_f21_maclaurin`/`_f21_zinf`, Punctures `_f21_term`).
    AHA: Teukolsky's `safe_2f1` z→∞ branch-SELECTION heuristic is a per-caller decision and
    stays LOCAL to Teukolsky — it is NOT pulled here.
  - `kummer.jl` — `kummer_hyp`, the HypGeo-DELEGATED ₁F₁ confluent. DISTINCT named function from
    NumericsCore's series-self-contained `kummer_1f1` (AHA: HypGeo-delegated vs series-self-contained
    are different decisions ⇒ two named fns, not one fork).

Type-generic `{T}` throughout so BigFloat precision is inherited, never hard-Float64'd.
Legacy packages whose primitives this stratum absorbs (kept green as the characterization oracle):
Teukolsky, H1Lorenz, Punctures, WaSABI.
"""
module BigSpecial

using NumericsCore
using HypergeometricFunctions
using SpecialFunctions

# Module files are added incrementally as each primitive is TDD-ported + characterization-gated.

include("pochhammer.jl")
export pochhammer

include("biggamma.jl")
export gamma_big, loggamma_big, needs_refl, fact_big, factorial_big

include("digamma.jl")
export digamma_big, PolyGamma

include("hyp2f1_primitive.jl")
export f21_maclaurin, f21_zinf, safe_2f1, f21_terminating

include("kummer.jl")
export kummer_hyp

include("legendre_p3.jl")
export legendre_P3, legendre_P3_dz

end # module BigSpecial

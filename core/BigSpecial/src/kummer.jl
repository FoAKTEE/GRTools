# Confluent ₁F₁ (Kummer M) — the HypGeo-DELEGATED variant. This is a DISTINCT named function from
# NumericsCore's self-contained series `kummer_1F1`: delegating to HypergeometricFunctions (its
# research-grade branch handling + BigFloat support via `pFq`) is a different decision from the
# vendored Maclaurin sum, so per AHA they stay two named primitives rather than one forked fork.
# It lives in BigSpecial because it needs the HypergeometricFunctions dep that the NumericsCore
# spine must not carry.

"""
    kummer_hyp(a, b, z) -> typeof(z)

Confluent hypergeometric `₁F₁(a;b;z)` (Kummer's `M`), delegated to
`HypergeometricFunctions.pFq((a,),(b,),z)` — type-generic, so `Complex{BigFloat}` precision flows
through. Distinct from `NumericsCore.kummer_1F1` (the self-contained vendored series).
"""
kummer_hyp(a, b, z) = HypergeometricFunctions.pFq((a,), (b,), z)

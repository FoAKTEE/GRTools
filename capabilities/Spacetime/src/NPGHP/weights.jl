# weights.jl — GHP boost/spin weight bookkeeping: {p, q} = (s + b, −s + b)
# (GHP 1973; TeukEvolution GHP.jl:205,232,262-263). Pure integer algebra.

"""
    ghp_weights(s, b) -> (p, q)

GHP type `{p, q} = (s + b, −s + b)` of a scalar with spin-weight `s` and boost-weight `b`
(Geroch–Held–Penrose 1973; TeukEvolution `GHP.jl:205,232,262-263`). Integer-preserving.

The falloff (peeling) power is a SEPARATE datum carried per field (R_teuk.md §2 table) and is not a
function of `(s, b)`; it enters the þ operator as `R·∂_R + falloff`, applied consumer-side.
"""
ghp_weights(s::Integer, b::Integer) = (s + b, -s + b)

"""
    ghp_prime_weights(p, q) -> (−p, −q)

GHP prime (l↔n, m↔m̄) maps type `{p, q} → {−p, −q}` — equivalently `(s, b) → (−s, −b)`.
"""
ghp_prime_weights(p::Integer, q::Integer) = (-p, -q)

"""
    ghp_bar_weights(p, q) -> (q, p)

GHP complex conjugation (bar, m↔m̄ only) maps type `{p, q} → {q, p}` — equivalently `(s,b)→(−s,b)`.
"""
ghp_bar_weights(p::Integer, q::Integer) = (q, p)

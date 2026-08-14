# Hyperboloidal-layer compactification (Paper 1 §II.D, Paper 2 §II.A) — verbatim port of
# DG-Teulolsky/src/Hyperboloidal.jl.
#
# `r_* = ρ/Ω(ρ)`,  `Ω = 1 - ((ρ-R)/(s-R))^P` for `ρ ≥ R` (else 1),  future null infinity 𝓘⁺ at
# `ρ = s`.  The layer function `H = 1 - Ω²/(Ω - ρΩ')` is `0` for `ρ < R` and `1` at `ρ = s`.
# All quantities here are finite on `[ρ_L, s]` *including* `ρ = s` (the apparent `1/Ω`
# singularities cancel) — the caller must use the finite combinations
# (e.g. `A = (1/(1+H))[[H,1],[1,H]]`), never `1/Ω` directly.

"""
    Layer(R, s, P)

Hyperboloidal layer: layer start `R`, scri location `s`, integer smoothness `P`
(Paper 1 §II.D). `Ω = 1 - ((ρ-R)/(s-R))^P` for `ρ ≥ R`, else `1`.
"""
struct Layer{T<:AbstractFloat}
    R::T
    s::T
    P::Int
end
Layer(R::Real, s::Real, P::Integer) = Layer{Float64}(Float64(R), Float64(s), Int(P))

"""
    Ωfun(L::Layer, ρ) -> Ω(ρ)

Layer compactification factor `Ω = 1 - ((ρ-R)/(s-R))^P` for `ρ ≥ R` (else `1`); `Ω(s) = 0`.
"""
@inline function Ωfun(L::Layer{T}, ρ::Real) where {T}
    ρ < L.R && return one(T)
    return one(T) - ((T(ρ) - L.R) / (L.s - L.R))^L.P
end

"""
    Ωprime(L::Layer, ρ) -> Ω'(ρ)

Derivative `Ω' = -(P/(s-R))((ρ-R)/(s-R))^(P-1)` for `ρ ≥ R` (else `0`).
"""
@inline function Ωprime(L::Layer{T}, ρ::Real) where {T}
    ρ < L.R && return zero(T)
    return -T(L.P) / (L.s - L.R) * ((T(ρ) - L.R) / (L.s - L.R))^(L.P - 1)
end

"""
    denomΩ(L::Layer, ρ) -> Ω - ρΩ'

`Ω - ρΩ'`  (the denominator in `H`; strictly positive on `[ρ_L, s]`).
"""
@inline denomΩ(L::Layer, ρ::Real) = Ωfun(L, ρ) - ρ * Ωprime(L, ρ)

"""
    Hfun(L::Layer, ρ) -> H(ρ)

Layer function `H = 1 - Ω²/(Ω - ρΩ')`  (= 0 for ρ<R, → 1 at ρ=s).
"""
@inline function Hfun(L::Layer{T}, ρ::Real) where {T}
    ρ < L.R && return zero(T)
    Ωv = Ωfun(L, ρ)
    return one(T) - Ωv^2 / (Ωv - ρ * Ωprime(L, ρ))
end

"""
    omega_factor(L::Layer, ρ) -> 1 - H

`1 - H = Ω²/(Ω - ρΩ')` (the paper's `ω`).
"""
@inline omega_factor(L::Layer{T}, ρ::Real) where {T} = one(T) - Hfun(L, ρ)

"""
    rstar_of_ρ(L::Layer, ρ) -> r_*

`r_* = ρ/Ω`  (→ ∞ at ρ=s; only use where finite).
"""
@inline rstar_of_ρ(L::Layer, ρ::Real) = ρ / Ωfun(L, ρ)

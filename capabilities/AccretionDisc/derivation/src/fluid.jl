# fluid.jl — the covariant fluid equations of arXiv:2608.17003 Eq. (3), written ONCE and
# element-type generic.
#
# The SAME four functions are evaluated twice: with `Num` scalars and `∂bg` they give the background
# system (whose vanishing is the stage-1 certificate of (6),(7),(11)); with the λ-jets of jets.jl and
# `∂jet` their O(λ) parts ARE the linearised system (15b),(15c) plus the perturbed normalisation —
# which is why no perturbed Christoffel, no δq^{μν} and no δ(u^ν∇_νu^μ) is ever written by hand.
#
# Conventions: signature (−,+,+,+), coordinates (t,r,φ,z), Γ^λ_{μν} from `christoffel_from_metric`,
# `∂(x, α)` = ∂_α x.  Continuity (3a) decouples (M:591) and is not used.

"""
    four_acceleration_up(Γ, u, ∂) -> a^μ

a^μ = u^ν∇_ν u^μ = u^ν∂_ν u^μ + Γ^μ_{νσ}u^ν u^σ.
"""
function four_acceleration_up(Γ, u, ∂)
    n = length(u)
    [sum(u[ν] * ∂(u[μ], ν) for ν in 1:n) + sum(Γ[μ, ν, σ] * u[ν] * u[σ] for ν in 1:n, σ in 1:n)
     for μ in 1:n]
end

"""
    divergence_u(Γ, u, ∂) -> ∇_α u^α

∇_α u^α = ∂_α u^α + Γ^α_{ασ}u^σ.
"""
function divergence_u(Γ, u, ∂)
    n = length(u)
    sum(∂(u[α], α) for α in 1:n) + sum(Γ[α, α, σ] * u[σ] for α in 1:n, σ in 1:n)
end

"""
    euler_up(ginv, Γ, u, pp, ee, ∂) -> ℰ^μ

The Euler equation (3c) with the index RAISED,
ℰ^μ = (ε+p) u^ν∇_ν u^μ + q^{μα}∇_α p, q^{μα} = u^μu^α + g^{μα} (2).
"""
function euler_up(ginv, Γ, u, pp, ee, ∂)
    n = length(u)
    a = four_acceleration_up(Γ, u, ∂)
    dpv = [∂(pp, α) for α in 1:n]
    [(ee + pp) * a[μ] + sum((u[μ] * u[α] + ginv[μ, α]) * dpv[α] for α in 1:n) for μ in 1:n]
end

"""
    energy_eq(Γ, u, pp, ee, ∂) -> u^μ∇_με + (ε+p)∇_αu^α

The energy-conservation equation (3b).
"""
function energy_eq(Γ, u, pp, ee, ∂)
    n = length(u)
    sum(u[μ] * ∂(ee, μ) for μ in 1:n) + (ee + pp) * divergence_u(Γ, u, ∂)
end

"""
    norm_constraint(g, u) -> g_{μν}u^μu^ν + 1

Four-velocity normalisation; its O(λ) part is `2u_αδu^α + h_{αβ}u^αu^β = 0`, step (i) of the
reduction scheme (M:604).
"""
function norm_constraint(g, u)
    n = length(u)
    sum(g[μ, ν] * u[μ] * u[ν] for μ in 1:n, ν in 1:n) + 1
end

"Lower an index: V_α = g_{αμ}V^μ."
lower(g, V) = [sum(g[α, μ] * V[μ] for μ in eachindex(V)) for α in eachindex(V)]
"Raise an index: V^α = g^{αμ}V_μ."
raise(ginv, V) = [sum(ginv[α, μ] * V[μ] for μ in eachindex(V)) for α in eachindex(V)]

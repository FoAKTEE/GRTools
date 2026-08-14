# operators.jl — the weighted NP transport operators of Li–Yunes, as an APPLICATION layer over the
# Spacetime.NPGHP coefficient kernel (PLAN §Architectural principle: NPGHP owns pointwise
# coefficients; operator application stays consumer-side).
#
# LY eq:derivative-notation (main.tex:203–212) defines four integer-weighted directional operators:
#
#   D_{[a,b,c,d]}       = D + a ε + b ε̄ + c ρ + d ρ̄            (l-leg,  D  = l^μ∂_μ)
#   Δ_{[a,b,c,d]} (bold)= Δ + a μ + b μ̄ + c γ + d γ̄            (n-leg,  Δ  = n^μ∂_μ)
#   δ_{[a,b,c,d]}       = δ + a ᾱ + b β + c π̄ + d τ            (m-leg,  δ  = m^μ∂_μ)
#   δ̄_{[a,b,c,d]}       = δ̄ + a α + b β̄ + c π + d τ̄            (m̄-leg, δ̄ = m̄^μ∂_μ)
#
# The whole transport hierarchy runs along Δ_{[…]} on the LHS (LY:319 "a first-order radial transport
# equation along the ingoing principal null direction n^μ"). This file provides:
#   • the `WeightedOp` label + `weight_combo` (the additive spin-coefficient scalar, from ANY NP
#     table — Kerr Kinnersley included, via Spacetime.NPGHP.np_spincoeff), and
#   • the static-Schwarzschild specialization (a=0): the ingoing derivative Δ reduces to −(f/2)∂_r
#     (Kinnersley n^r = −f/2, LY:736), and the coefficient combination reduces to a Laurent in r.
#
# {T}-generic; the background coefficients are SOURCED from Spacetime.NPGHP, never re-hardcoded
# (the a=0 Laurent table below is gated field-by-field against NPGHP.np_spincoeff(0,r,θ), G-E1b).

import Spacetime

"""
    WeightedOp

A weighted NP directional operator label: `leg ∈ (:D, :Δ, :δ, :δbar)` and the four integer weights
`(a,b,c,d)` of LY eq:derivative-notation. Carries no data beyond the label — application is provided
by [`weight_combo`](@ref) (the algebraic part) plus the caller's directional derivative.
"""
struct WeightedOp
    leg::Symbol
    a::Int
    b::Int
    c::Int
    d::Int
    function WeightedOp(leg::Symbol, a::Integer, b::Integer, c::Integer, d::Integer)
        leg in (:D, :Δ, :δ, :δbar) || error("WeightedOp: leg must be one of :D,:Δ,:δ,:δbar (got $leg)")
        new(leg, a, b, c, d)
    end
end

"Δ_{[a,b,c,d]} — the ingoing (n-leg) weighted operator (LY:207)."
Δw(a, b, c, d) = WeightedOp(:Δ, a, b, c, d)
"D_{[a,b,c,d]} — the outgoing (l-leg) weighted operator (LY:204)."
Dw(a, b, c, d) = WeightedOp(:D, a, b, c, d)
"δ_{[a,b,c,d]} — the m-leg weighted operator (LY:208)."
δw(a, b, c, d) = WeightedOp(:δ, a, b, c, d)
"δ̄_{[a,b,c,d]} — the m̄-leg weighted operator (LY:210)."
δbarw(a, b, c, d) = WeightedOp(:δbar, a, b, c, d)

"""
    weight_combo(op::WeightedOp, np) -> scalar

The additive spin-coefficient combination of `op` evaluated on an NP coefficient table `np` (a
NamedTuple with fields `ε, ρ, μ, γ, α, β, π, τ` — e.g. the output of
`Spacetime.NPGHP.np_spincoeff(a, r, θ)`), per the leg's quartet in LY eq:derivative-notation:

    :D     → a·ε + b·conj(ε) + c·ρ + d·conj(ρ)
    :Δ     → a·μ + b·conj(μ) + c·γ + d·conj(γ)
    :δ     → a·conj(α) + b·β + c·conj(π) + d·τ
    :δbar  → a·α + b·conj(β) + c·π + d·conj(τ)

Generic over the coefficient element type (Kerr Kinnersley `Complex{T}` or the real Schwarzschild
reduction). Returns ONLY the algebraic part; the caller supplies the directional derivative.
"""
function weight_combo(op::WeightedOp, np)
    if op.leg === :D
        op.a * np.ε + op.b * conj(np.ε) + op.c * np.ρ + op.d * conj(np.ρ)
    elseif op.leg === :Δ
        op.a * np.μ + op.b * conj(np.μ) + op.c * np.γ + op.d * conj(np.γ)
    elseif op.leg === :δ
        op.a * conj(np.α) + op.b * np.β + op.c * conj(np.π) + op.d * np.τ
    else # :δbar
        op.a * np.α + op.b * conj(np.β) + op.c * np.π + op.d * conj(np.τ)
    end
end

# ------------------------------------------------------------------------------------------------
# Static-Schwarzschild (a=0) specialization
# ------------------------------------------------------------------------------------------------
"""
    schw_f(T, M) -> Laurent{T}

The Schwarzschild metric function `f(r) = 1 − 2M/r` as a Laurent (`laurent(0=>1, -1=>-2M)`).
"""
schw_f(::Type{T}, M) where {T} = laurent(T, 0 => one(T), -1 => T(-2M))

"""
    schw_np_laurent(T, M) -> NamedTuple

The a=0 Schwarzschild Kinnersley NP spin coefficients that enter the radial (monopole) transport
chain, each as a `Laurent{T}` in r (M=area-mass):

    ρ = −1/r ,   μ = −f(r)/(2r) ,   γ = M/(2r²) ,   ψ₂ = −M/r³ ,   τ = π = ε = 0 .

These are the ANCHOR-2 closed forms (LY:743–748); they are GATED field-by-field against
`Spacetime.NPGHP.np_spincoeff(0, r, θ)` (test G-E1b) — this table is an application-side cache of the
NPGHP kernel, not an independent hardcode. `α, β` are θ-dependent (cotθ/2√2 r) and do not enter the
spherically symmetric radial chain, so are omitted here.
"""
function schw_np_laurent(::Type{T}, M) where {T}
    Mt = T(M)
    (ρ  = laurent(T, -1 => -one(T)),                        # −1/r
     μ  = laurent(T, -1 => -one(T) / 2, -2 => Mt),          # −f/2r = −1/2r + M/r²
     γ  = laurent(T, -2 => Mt / 2),                         # M/2r²
     τ  = lzero(T),
     π  = lzero(T),
     ε  = lzero(T),
     ψ₂ = laurent(T, -3 => -Mt))                            # −M/r³
end

"""
    delta_static(op::WeightedOp, X::ShellProfile, M) -> ShellProfile

Apply the static-Schwarzschild reduction of a Δ_{[a,b,c,d]} operator (`op.leg === :Δ`) to a radial
field `X`. For a STATIC field the ingoing derivative is `Δ = n^μ∂_μ = n^r ∂_r = −(f/2)∂_r`
(Kinnersley n^μ = (1/2, −f/2, 0, 0), LY:736), so

    Δ_{[a,b,c,d]} X = −(f/2)·∂_r X + (a+b)·μ·X + (c+d)·γ·X

with the real Schwarzschild μ, γ (conjugates equal). Returns a `ShellProfile` (the algebra stays in
the distributional class). This is the residual operator the shell gates apply to the reconstructed
fields to confirm they satisfy the transport equations.
"""
function delta_static(op::WeightedOp, X::ShellProfile{T}, M) where {T}
    op.leg === :Δ || error("delta_static: only the Δ (n-leg) reduction is implemented (got $(op.leg))")
    np = schw_np_laurent(T, M)
    f = schw_f(T, M)
    combo = (op.a + op.b) === 0 ? lzero(T) : lscale(op.a + op.b, np.μ)
    combo = combo + ((op.c + op.d) === 0 ? lzero(T) : lscale(op.c + op.d, np.γ))
    # −(f/2)∂_r X  +  combo·X
    pscale(-one(T) / 2, psmooth(f, pderiv(X))) + psmooth(combo, X)
end

# spincoeff_bl.jl — Kerr Kinnersley Newman–Penrose spin-coefficient table in the Boyer–Lindquist
# chart (M = 1, signature −+++, ζ = r − i a cosθ), and its GHP-primed view.
#
# Closed forms (type-D Kinnersley: κ=σ=ν=λ=ε≡0) with ζ = r − i a cosθ, Σ = ζζ̄, Δ = r²−2r+a²:
#     ρ = −1/ζ                     μ = −Δ/(2 Σ ζ)              τ = −i a sinθ /(√2 Σ)
#     π =  i a sinθ /(√2 ζ²)       β =  cotθ /(2√2 ζ̄)          α =  π − conj(β)
#     γ =  μ + (r−1)/(2 Σ)         ψ₂ = −1/ζ³
# Admitted by the first-principles oracle (test/npghp_oracle.jl, gate G-B2).
# [AXIOM] sign pins: ρ = −l_{μ;ν}m^μm̄^ν class; n-leg block plus; MTW/Wald ⇒ ψ₂=−1/ζ³ no extra sign.
# √2 = sqrt(T(2)) — a Float64 literal here would cap BigFloat at ~1e-16.

using ..Spacetime: zeta_kerr, sigma_kerr, delta_kerr

"""
    np_spincoeff(a, r, θ; coords=:BL) -> NamedTuple

Full Kerr Kinnersley Newman–Penrose spin-coefficient table at areal radius `r`, polar angle `θ`,
spin `a` (M = 1, signature −+++). Returns `(; κ, ρ, σ, τ, π, ν, μ, λ, ε, γ, β, α, ψ₂)` as
`Complex{T}` with `T` promoted from `(a, r, θ)`. The type-D vanishing set `κ=σ=ν=λ=ε` is returned
as EXACT zeros; the nonzero `ρ,τ,π,μ,γ,β,α` and the Weyl scalar `ψ₂ = −1/ζ³` are the closed forms
above. Gated by the first-principles oracle (G-B2), reduces to `ModifiedGravity.schw_spincoeff` at
a=0/equator (G-B3).

`coords ∈ (:BL, :In, :Out)`. The `:In`/`:Out` Kinnersley tetrads are pure coordinate
push-forwards of the BL legs (no boost, no phase), and NP spin coefficients are frame
scalars — the table is chart-invariant (gate G-F5 recomputes it in each chart).
"""
function np_spincoeff(a, r, θ; coords::Symbol = :BL)
    coords in (:BL, :In, :Out) ||
        throw(ArgumentError("np_spincoeff: unknown coords $(repr(coords)) (expected :BL, :In, :Out)"))
    T  = float(promote_type(typeof(a), typeof(r), typeof(θ)))
    aT = T(a); rT = T(r); θT = T(θ)
    ζ  = zeta_kerr(aT, rT, θT)                 # r − i a cosθ, Complex{T}
    ζ̄  = conj(ζ)
    Σ  = sigma_kerr(aT, rT, θT)                # real T
    Δ  = delta_kerr(aT, rT)                    # real T
    sθ = sin(θT); cθ = cos(θT)
    rt2 = sqrt(T(2))
    z0 = zero(Complex{T})

    ρ = -inv(ζ)
    μ = -Δ / (2 * Σ * ζ)
    τ = -im * aT * sθ / (rt2 * Σ)
    π = im * aT * sθ / (rt2 * ζ^2)
    β = (cθ / sθ) / (2 * rt2 * ζ̄)              # cotθ /(2√2 ζ̄)
    α = π - conj(β)
    γ = μ + (rT - one(T)) / (2 * Σ)            # μ + (r−M)/(2Σ), M=1
    ψ₂ = -inv(ζ^3)

    (; κ = z0, ρ = Complex{T}(ρ), σ = z0, τ = Complex{T}(τ), π = Complex{T}(π),
       ν = z0, μ = Complex{T}(μ), λ = z0, ε = z0, γ = Complex{T}(γ),
       β = Complex{T}(β), α = Complex{T}(α), ψ₂ = Complex{T}(ψ₂))
end

"""
    ghp_spincoeff(a, r, θ; coords=:BL) -> NamedTuple

The GHP-PRIMED spin-coefficient table — the NP coefficients of the prime tetrad (l↔n, m↔m̄) at
areal radius `r`, polar angle `θ`, spin `a` (M = 1). Returns `(; κ, ρ, σ, τ, π, ν, μ, λ, ε, γ, β,
α, ψ₂)` where each field holds the PRIMED coefficient (`.ρ` = ρ′, `.μ` = μ′, `.ε` = ε′, …), obtained
from the standard GHP prime relations applied to `np_spincoeff`:
    κ′=−ν  ρ′=−μ  σ′=−λ  τ′=−π  π′=−τ  ν′=−κ  μ′=−ρ  λ′=−σ  ε′=−γ  γ′=−ε  β′=−α  α′=−β,
and ψ₂ prime-invariant (GHP type {0,0}). The prime MAP is not trusted from the a=0 pin: G-B2 gates
this field-by-field against the LEG-SWAPPED oracle at full Kerr. Consumers using the [AXIOM-4]
convention ρ′=−μ, ε′=−γ (e.g. `ModifiedGravity.thorn_prime`) read this view directly (G-B3).
`coords ∈ (:BL, :In, :Out)` is chart-invariant (the prime tetrad is the same push-forward, F5).
"""
function ghp_spincoeff(a, r, θ; coords::Symbol = :BL)
    np = np_spincoeff(a, r, θ; coords = coords)
    (; κ = -np.ν, ρ = -np.μ, σ = -np.λ, τ = -np.π, π = -np.τ, ν = -np.κ,
       μ = -np.ρ, λ = -np.σ, ε = -np.γ, γ = -np.ε, β = -np.α, α = -np.β, ψ₂ = np.ψ₂)
end

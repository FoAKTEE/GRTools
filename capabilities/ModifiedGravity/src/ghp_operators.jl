# Native GHP þ/þ′ directional derivatives on equatorial Schwarzschild (Sourceincoordinates.nb
# SC l.368: þ′ ≡ PDNP[{2,−NP}], the n-leg derivative). Modes ∝ e^{iωt} (∂_t→iω); ε=ε†=0 and
# ε′=ε′† here, so only the weight-SUM p+q enters þ′. Radial functions travel as NumericsCore
# jet towers [g,g′,…] at r0 (∂_r = jderiv).

using NumericsCore: jadd, jscale, jmul, jderiv, jinv

"""
    schw_spincoeff(r; M=1) -> NamedTuple

Equatorial Schwarzschild GHP spin coefficients at Boyer-Lindquist radius `r` (Schwspincoeff reduction
rules, Sourceincoordinates.nb): `ρ=−1/r`, `ρ′=(r−2M)/(2r²)`, `ε′=ε′†=−M/(2r²)`, `ψ₂=−M/r³`, and the
vanishing set `κ=σ=τ=τ′=ε=ε†=0`. Type-generic in `r`/`M`.
"""
function schw_spincoeff(r::Real; M::Real=1)
    T = typeof(float(r*M))
    rT = T(r); MT = T(M)
    (ρ  = -one(T)/rT,
     ρ′ = (rT - 2MT)/(2rT^2),
     ε′ = -MT/(2rT^2),
     ε′dag = -MT/(2rT^2),
     ψ₂ = -MT/rT^3,
     κ = zero(T), σ = zero(T), τ = zero(T), τ′ = zero(T), ε = zero(T), εdag = zero(T))
end

# ── coefficient jets at r0 (reusing Operators jet algebra) ───────────────────────────────────────
"Length-(N+1) coefficient jets (r², Δ, r²/Δ, Δ/2r², ε′=−M/2r²) at r0, element type Complex{T}."
function _coeff_jets(r0::Real, M::Real, N::Integer)
    T = typeof(float(r0*M)); CT = Complex{T}
    r0T = T(r0); MT = T(M)
    rjet = zeros(CT, N+1); rjet[1] = r0T; N ≥ 1 && (rjet[2] = one(CT))
    r2   = jmul(rjet, rjet)
    Δjet = zeros(CT, N+1); Δjet[1] = r0T^2 - 2MT*r0T; N ≥ 1 && (Δjet[2] = 2r0T - 2MT); N ≥ 2 && (Δjet[3] = CT(2))
    invΔ = jinv(Δjet)
    r2overΔ = jmul(r2, invΔ)
    Δover2r2 = jmul(Δjet, jinv(jscale(CT(2), r2)))
    εpjet = jscale(CT(-MT), jinv(jscale(CT(2), r2)))    # ε′ = −M/(2r²)
    (CT=CT, rjet=rjet, r2overΔ=r2overΔ, Δover2r2=Δover2r2, εpjet=εpjet)
end

"""
    thorn(jet, p, q, r0, ω; M=1) -> jet

GHP thorn `þ = D − pε − qε̄` on the type-`{p,q}` mode whose radial Taylor tower at `r0` is `jet`
(`∝ e^{iωt}`, `∂_t→iω`). On equatorial Schwarzschild `ε=ε†=0`, so `þ = iω r²/Δ + ∂_r`; maps
`{p,q}→{p+1,q+1}`. Returns the tower of `þη` at `r0`.
"""
function thorn(jet::AbstractVector, p::Integer, q::Integer, r0::Real, ω::Real; M::Real=1)
    N = length(jet) - 1
    c = _coeff_jets(r0, M, N)
    ωC = c.CT(ω)
    jadd(jderiv(jet), jscale(im*ωC, jmul(c.r2overΔ, jet)))
end

"""
    thorn_prime(jet, p, q, r0, ω; M=1) -> jet

GHP thorn-prime `þ′ = Δ′ + pε′ + qε̄′ = PDNP[{2,−NP}]` (the n-leg derivative, Sourceincoordinates.nb
SC l.368) on the type-`{p,q}` mode whose radial Taylor tower at `r0` is `jet` (`∝ e^{iωt}`,
`∂_t→iω`). On equatorial Schwarzschild `ε′=ε̄′=−M/2r²` so only the weight-sum `p+q` enters:
`þ′ = iω/2 − (Δ/2r²)∂_r + (p+q)ε′`; maps `{p,q}→{p−1,q−1}`. Returns the tower of `þ′η` at `r0`.
"""
function thorn_prime(jet::AbstractVector, p::Integer, q::Integer, r0::Real, ω::Real; M::Real=1)
    N = length(jet) - 1
    c = _coeff_jets(r0, M, N)
    ωC = c.CT(ω); S = c.CT(p + q)
    jadd(jscale(im*ωC/2, jet),
             jadd(jscale(-one(c.CT), jmul(c.Δover2r2, jderiv(jet))),
                      jscale(S, jmul(c.εpjet, jet))))
end

export schw_spincoeff, thorn, thorn_prime

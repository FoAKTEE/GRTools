# RT-A2 velocity-reconstruction operator Q^{-1} (arXiv:2601.19123 l.379) + its inverse Q, the
# principal-symbol building block of the height-averaged reduced operator (MasterPDE). At the
# equator (θ=π/2) the perturbed fluid velocity is δu^μ = Q^μ_ν · (reconstruction source), with
#
#   Q^{-1}^ν_μ = iω̃α δ^ν_μ + Γ^ν_{βμ}U^β + ∇_μU^ν + U^ν R_μ + R^ν B_μ/(iαω̃),
#
#   U^μ = α(1,0,0,Ω) (circular equatorial 4-velocity), R^ν = ∇^νp/(ε+p) (pressure-gradient
#   accel.), B_μ = −∂_μ log c_s²/(1+c_s²). ω̃ = m(Ω_p−Ω) is the comoving frequency; Ω=Ω^k is the
#   geometric circular frequency (the fluid drift δω is O(h²), a PMLA refinement).
#
# AHA: the legacy in-package FD `kerr_christoffel` fork is DROPPED — the connection comes from the
# modern analytic `Spacetime.christoffel` (Γ^λ_{μν}, λ first, same convention). Complex{T},
# BigFloat-inherited. From DiscModes/relativistic_master.jl/build_Qinv,build_Q.

using ..Modes: RelativisticDisc, disc_OmegaK, csound, eden0
import Spacetime: kerr_metric, kerr_metric_inv, christoffel

"""
    qinv_operator(d::RelativisticDisc, r, r_p, m) -> 4×4 Complex{T}

Assemble the velocity-reconstruction operator `Q^{-1}^ν_μ` (mixed tensor `M[ν,μ]`) at the equator
for radius `r`, perturber radius `r_p`, azimuthal mode `m` (RT-A2, l.379). Reuses the analytic
`Spacetime.{kerr_metric,kerr_metric_inv,christoffel}`. Singular at corotation (ω̃→0) and `m=0`.
"""
function qinv_operator(d::RelativisticDisc{T}, r, r_p, m) where {T}
    a = d.a; θ = T(π) / 2; hh = sqrt(eps(T)) * max(abs(r), one(T))
    gi = kerr_metric_inv(a, r, θ); Γ = christoffel(a, r, θ)
    Ω = disc_OmegaK(d, r); Ωp = disc_OmegaK(d, r_p)
    αof(s) = (Ωs = disc_OmegaK(d, s); gs = kerr_metric(a, s, θ);
              sqrt(-one(T) / (gs[1, 1] + 2 * Ωs * gs[1, 4] + Ωs^2 * gs[4, 4])))
    α = αof(r); ω̃ = m * (Ωp - Ω)
    U = T[α, zero(T), zero(T), α * Ω]                          # U^μ
    # R^ν = ∇^νp/(ε+p) (equatorial radial), R_μ; p = c_s² ε_0
    pof(s) = csound(d, s)^2 * eden0(d, s); dpdr = (pof(r + hh) - pof(r - hh)) / (2hh)
    εp = eden0(d, r) + pof(r)
    Rup = zeros(T, 4); Rup[2] = gi[2, 2] * dpdr / εp           # R^r
    Rlo = zeros(T, 4); Rlo[2] = dpdr / εp                      # R_r = g_rr R^r
    # B_μ = −∂_μ log c_s²/(1+c_s²) (equatorial radial)
    cs2of(s) = csound(d, s)^2
    dlogcs2 = (log(cs2of(r + hh)) - log(cs2of(r - hh))) / (2hh)
    Blo = zeros(T, 4); Blo[2] = -dlogcs2 / (1 + cs2of(r))      # B_r
    # ∂_μ U^ν (only μ=r): ∂_r U^t = α', ∂_r U^φ = (αΩ)'
    dα = (αof(r + hh) - αof(r - hh)) / (2hh)
    dαΩ = (αof(r + hh) * disc_OmegaK(d, r + hh) - αof(r - hh) * disc_OmegaK(d, r - hh)) / (2hh)
    dU = zeros(T, 4, 4); dU[2, 1] = dα; dU[2, 4] = dαΩ         # dU[μ,ν] = ∂_μ U^ν
    M = zeros(Complex{T}, 4, 4)
    for ν in 1:4, μ in 1:4
        t1 = (ν == μ) ? im * ω̃ * α : zero(Complex{T})
        t2 = sum(Γ[ν, β, μ] * U[β] for β in 1:4)               # Γ^ν_{βμ} U^β
        covD = dU[μ, ν] + sum(Γ[ν, μ, σ] * U[σ] for σ in 1:4)  # ∇_μ U^ν
        t4 = U[ν] * Rlo[μ]
        t5 = Rup[ν] * Blo[μ] / (im * α * ω̃)
        M[ν, μ] = t1 + t2 + covD + t4 + t5
    end
    M
end

"""
    q_operator(d::RelativisticDisc, r, r_p, m) -> 4×4 Complex{T}

Velocity-reconstruction operator `Q = (Q^{-1})^{-1}` (RT-A2). `Q·Q^{-1}=I` for `m≠0` away from
corotation (singular at corotation / `m=0`).
"""
q_operator(d::RelativisticDisc, r, r_p, m) = inv(qinv_operator(d, r, r_p, m))

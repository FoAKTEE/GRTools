# Off-equatorial (z=cosθ) disc background to O(z²)~O(h²) + the RT-A2 velocity-reconstruction
# operator OFF the equator + the PMLA Gaussian height-average (design §4.10 Modes: "off-equatorial
# Q-corrections, Gaussian height-averaging"). The PMLA expansion of the Kerr disc background
# (arXiv:2601.19123 eq:backgroundVel): the poloidal (θ) structure of Ω^k, α^k, U^μ, R_μ and the
# metric is what the Gaussian height-average + IBP converts into the radial ODE's A_m/B_m, and what
# the Force pipeline's 𝒬-projected source needs off the midplane. At z=0 everything reduces to the
# verified equatorial pieces (Background.jl / Velocity.jl). M=1, BL coords.
#
# AHA: the connection comes from the analytic `Spacetime.christoffel` (Γ^λ_{μν}, λ first — the same
# convention as the legacy FD `DiscModes.kerr_christoffel`, collapsed per the Velocity.jl dedup
# decision). The legacy `F6_NOT5` env experiment knob (iter289, default-off) is NOT ported: the
# validated default (t5 pressure term KEPT) is hardwired. Type-generic {T}; BigFloat-inherited.
# From DiscModes/off_equator.jl + DiscModes/heightaverage.jl.

using ..Modes: RelativisticDisc, disc_OmegaK, csound, eden0
import Spacetime: kerr_metric, kerr_metric_inv, christoffel

# purely-geometric z² coefficients (eq:backgroundVel): Ω^k=Ω^k_eq(1+gΩ z²), α^k=α^k_eq(1−gα z²)
_gOmega(a, r) = (r^2 + 2a*sqrt(r) - 3a^2) / (2*sqrt(r)*(a + r^(3//2)))
_alphaEq(a, r) = (a + r^(3//2)) / sqrt(2a*r^(3//2) + (r - 3)*r^2)          # |α^k_eq| (U^t>0 convention)
_gAlpha(a, r) = (3a^4 - 6a^3*sqrt(r) + 4a^2*r*(r+1) - r^4) /
                (2*r^2*(2a + (r-3)*sqrt(r))*(a + r^(3//2)))

"Off-equatorial geometric angular frequency Ω^k(r,z) to O(z²) (eq:backgroundVel). z=cosθ."
disc_OmegaK_z(d::RelativisticDisc, r, z) = disc_OmegaK(d, r) * (1 + _gOmega(d.a, r)*z^2)

"Off-equatorial lapse-like normalisation α^k(r,z)=U^t to O(z²) (eq:backgroundVel; |·| with U^t>0)."
alpha_z(d::RelativisticDisc, r, z) = _alphaEq(d.a, r) * (1 - _gAlpha(d.a, r)*z^2)

"Background 4-velocity U^μ(r,z)=α{1,0,0,Ω} (geometric, BL; pressure δω,δα are O(h²) sub-leading)."
function Uvec_z(d::RelativisticDisc{T}, r, z) where {T}
    al = alpha_z(d, r, z); Om = disc_OmegaK_z(d, r, z)
    T[al, 0, 0, al*Om]
end

# Gaussian vertical pressure profile p(r,z)=p_0(r)e^{−z²/2h²} (eq:pressureprolfe1); ε(r,z)=p/c_s²
# (poloidally isothermal: c_s²=c_s²(r)). ε_{(0,1)}=ε_I(r/r_I)^{-2} equatorial.
"Gaussian vertical pressure p(r,z)=c_s²(r)·ē(r)·e^{−z²/2h₀²} (eq:pressureprolfe1). z=cosθ."
pressure_z(d::RelativisticDisc, r, z) = csound(d, r)^2 * eden0(d, r) * exp(-z^2/(2*d.h0^2))
"Gaussian vertical energy density ε(r,z)=ē(r)·e^{−z²/2h₀²} (poloidally isothermal). z=cosθ."
eden_z(d::RelativisticDisc, r, z) = eden0(d, r) * exp(-z^2/(2*d.h0^2))

"""
    qinv_operator_z(d, r, z, r_p, m) -> 4×4 Complex{T}

RT-A2 velocity-reconstruction operator `Q^{-1}` OFF the equator (the `qinv_operator` generalization,
θ=acos z): the same operator with the full poloidal (z) background — `U^μ(r,z)`, `R_μ` with `R_θ≠0`
(Gaussian pressure), `∂_θ U^ν`, `Γ(r,θ)`, `metric(r,θ)`. At z=0 it reduces to the equatorial
`qinv_operator`. Returns the mixed tensor `M[ν,μ]=(Q^{-1})^ν_μ`. From DiscModes/off_equator.jl
`build_Qinv_z` (the validated default: the O(c_s²) t5 pressure term KEPT).
"""
function qinv_operator_z(d::RelativisticDisc{T}, r, z, r_p, m) where {T}
    a = d.a; θ = acos(z); hr = sqrt(eps(T))*max(abs(r), one(T)); hθ = sqrt(eps(T))
    gi = kerr_metric_inv(a, r, θ); Γ = christoffel(a, r, θ)
    Ω  = disc_OmegaK_z(d, r, z); Ωp = disc_OmegaK(d, r_p)
    α  = alpha_z(d, r, z);  ω̃ = m*(Ωp - Ω)
    U  = Uvec_z(d, r, z)                                          # U^μ(r,z)
    # background scalar fields as functions of (r,θ) for finite-diff (z=cosθ)
    Uf(rr, θθ)  = Uvec_z(d, rr, cos(θθ))
    pf(rr, θθ)  = csound(d, rr)^2 * eden0(d, rr) * exp(-cos(θθ)^2/(2*d.h0^2))
    εpf(rr, θθ) = eden0(d, rr) * exp(-cos(θθ)^2/(2*d.h0^2)) * (1 + csound(d, rr)^2)   # ε+p
    cs2(rr)     = csound(d, rr)^2
    # R_μ = ∂_μ p/(ε+p)  (radial + poloidal; ∂_t=∂_φ=0)
    Rlo = zeros(T, 4)
    Rlo[2] = ((pf(r+hr, θ) - pf(r-hr, θ))/(2hr)) / εpf(r, θ)
    Rlo[3] = ((pf(r, θ+hθ) - pf(r, θ-hθ))/(2hθ)) / εpf(r, θ)
    Rup = gi * Rlo                                               # R^ν
    # B_μ = −∂_μ log c_s²/(1+c_s²)  (radial only; c_s²=c_s²(r))
    Blo = zeros(T, 4)
    Blo[2] = -((log(cs2(r+hr)) - log(cs2(r-hr)))/(2hr)) / (1 + cs2(r))
    # ∂_μ U^ν  (μ=r,θ)
    dU = zeros(T, 4, 4)                                          # dU[μ,ν] = ∂_μ U^ν
    dU[2, :] = (Uf(r+hr, θ) .- Uf(r-hr, θ)) ./ (2hr)
    dU[3, :] = (Uf(r, θ+hθ) .- Uf(r, θ-hθ)) ./ (2hθ)
    M = zeros(Complex{T}, 4, 4)
    for ν in 1:4, μ in 1:4
        t1   = (ν == μ) ? im*ω̃*α : zero(Complex{T})
        t2   = sum(Γ[ν, β, μ]*U[β] for β in 1:4)                 # Γ^ν_{βμ} U^β
        covD = dU[μ, ν] + sum(Γ[ν, μ, σ]*U[σ] for σ in 1:4)      # ∇_μ U^ν
        t4   = U[ν]*Rlo[μ]
        t5   = Rup[ν]*Blo[μ]/(im*α*ω̃)                            # O(c_s²) pressure term (validated default: KEPT)
        M[ν, μ] = t1 + t2 + covD + t4 + t5
    end
    M
end

"""
    q_operator_z(d, r, z, r_p, m) -> 4×4 Complex{T}

Off-equatorial velocity-reconstruction operator `Q = (Q^{-1})^{-1}` at height z=cosθ (RT-A2).
Singular at corotation / the 𝒬-resonances. From DiscModes/off_equator.jl `build_Q_z`.
"""
q_operator_z(d::RelativisticDisc, r, z, r_p, m) = inv(qinv_operator_z(d, r, z, r_p, m))

# --- PMLA Gaussian height-average (RT-A1 reduction; eq:heightprojector A10) ---------------------------------
# EXACT Gaussian moments on the O(h²) Taylor expansion (the PMLA truncates the polar DOF to the n=0
# Hermite, O(h²)). For z=cosθ and width h: ⟨z^{2k}⟩ = (2k−1)!! h^{2k}, odd moments 0, so
# ⟨f⟩ = f(0) + (h²/2) f''(0) + (h⁴/8) f''''(0) + … (Taylor × moments; exact for the PMLA truncation).

# central finite-diff derivatives at 0 (for the operator's z-Taylor coefficients; type-generic)
_d2_0(f, hz) = (f(hz) - 2f(zero(hz)) + f(-hz)) / hz^2
_d4_0(f, hz) = (f(2hz) - 4f(hz) + 6f(zero(hz)) - 4f(-hz) + f(-2hz)) / hz^4

"""
    heightavg(f, h; order=2, hz=nothing)

PMLA Gaussian height-average ⟨f⟩ to O(h²): `f(0)+(h²/2)f''(0)`; `order=4` adds `(h⁴/8)f''''(0)`.
`hz` = z-step for the central-FD z-derivatives (default `eps(T)^(1/4)`, balancing truncation vs
roundoff). The vertical projection reducing the master PDE's z-structure to the radial ODE.
From DiscModes/heightaverage.jl.
"""
function heightavg(f, h::T; order::Int = 2, hz = nothing) where {T}
    hzz = hz === nothing ? sqrt(sqrt(eps(T))) : T(hz)
    r = f(zero(T)) + (h^2/2) * _d2_0(f, hzz)
    order >= 4 && (r += (h^4/8) * _d4_0(f, hzz))
    r
end

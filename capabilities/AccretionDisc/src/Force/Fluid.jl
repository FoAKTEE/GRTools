# FLUID — reconstruct the fluid perturbation (δU^μ, δp, δε) from the master variable 𝔥, and the
# advected angular-momentum flux density (design §4.10 Force: "fluid reconstruction + advected
# ang-mom flux"). From DiscEMRIForce/fluid.jl (h1lorenz fluid_recon.jl / adv_flux.jl).
#
# Physics (arXiv:2601.19123 eq:PresEnthrecon / denEnthecon / VelocityReconstruction):
#   p^{1,1} = (ε+p) 𝔥
#   δU^μ    = 𝒬^μ_ν[ prj^ν_α f^α/(ε+p) − prj^{να}∂_α 𝔥 − 𝔥 B^ν ]
#   ε^{1,1} = (ε+p)[ 𝔥/c_s² − B_ν δU^ν/(i ω̃ α) ]
#   J̇^r_Adv = (ε+p) Re[δu^r δu_φ*]              (the advected ang-mom flux density; one side of F8)
#
# `δU_wave` (drop the direct metric force f) is the wave-propagating part used in the F8 flux
# cross-check. TYPE-GENERIC {T} throughout.

using ..Modes: RelativisticDisc, disc_OmegaK, csound, eden0, Uvec_z, alpha_z, q_operator_z
import Spacetime: kerr_metric, kerr_metric_inv

"""
    FluidVars{T}

Reconstructed fluid perturbation at a radius: full velocity `δU` (Complex{T} 4-vector, upper index),
`δU_wave` (the 𝔥-response part with the direct metric force dropped), pressure `δp`, energy `δε`.
"""
struct FluidVars{T<:Real}
    δU::Vector{Complex{T}}; δU_wave::Vector{Complex{T}}; δp::Complex{T}; δε::Complex{T}
end

"""
    reconstruct(d, modes, k, m, r_p, ω0, h_val, hp; lmax) -> FluidVars{T}

Reconstruct (δU, δU_wave, δp, δε) from the master variable 𝔥 = `h_val` and 𝔥' = `hp` at radial
index `k`. The projector prj (=g^{να}+U^νU^α, upper-upper) acts on the LOWER force f_α (the
iter287-fixed q^ν_α f^α index placement). From DiscEMRIForce/fluid.jl (h1lorenz reconstruct).
"""
function reconstruct(d::RelativisticDisc, modes::Dict{Int,H1Mode{T}}, k::Integer, m::Integer, r_p, ω0,
                     h_val, hp; lmax::Integer) where {T}
    a = T(d.a); r = modes[m].r_right[k]; θ = T(π)/2
    gi = kerr_metric_inv(a, r, θ); U = Uvec_z(d, r, zero(T)); prj = gi .+ U*U'
    Q = q_operator_z(d, r, zero(T), r_p, m); cs2 = csound(d, r)^2; εp = eden0(d, r)*(1 + cs2)
    # B^ν : B_r = −∂_r log c_s²/(1+c_s²), raised   (validated FD step 1e-6)
    dr = T(1)/T(10)^6; cs2f(s) = csound(d, s)^2
    Blo = zeros(Complex{T}, 4); Blo[2] = -((log(cs2f(r+dr)) - log(cs2f(r-dr)))/(2dr))/(1 + cs2); Bup = gi*Blo
    # ∂_α 𝔥 : ∂_t=−iω₀𝔥, ∂_r=𝔥', ∂_θ=0 (n=0), ∂_φ=im𝔥
    dh = Complex{T}[-im*ω0*h_val, hp, 0, im*m*h_val]
    # direct metric force f (sum over l), parity-aware — reuse the Source-module force assembly
    f = zeros(Complex{T}, 4)
    for l in m:lmax
        haskey(modes, l) || continue
        f .+= force_mode(d, modes[l], k, l, m, ω0, zero(T))
    end
    δU_full = Q*((prj*f)/εp .- (prj*dh) .- h_val .* Bup)
    δU_wave = Q*(            (-(prj*dh)) .- h_val .* Bup)
    δp = εp*h_val
    ω̃ = ω0 - m*disc_OmegaK(d, r); α = alpha_z(d, r, zero(T))   # spin-general U^t (a→0 ⇒ 1/√(1−3/r))
    δε = εp*(h_val/cs2 - sum(Blo[ν]*δU_full[ν] for ν in 1:4)/(im*ω̃*α))
    FluidVars{T}(δU_full, δU_wave, δp, δε)
end

"""
    Jadv(d, modes, k, m, r_p, ω0, h_val, hp; lmax, wave=false) -> T

Advected angular-momentum flux density J̇^r_Adv = (ε+p) Re[δu^r δu_φ*] from the reconstructed fluid.
`wave=true` uses the wave-propagating part `δU_wave` (the conserved-flux tracker of the F8
cross-check). From DiscEMRIForce/fluid.jl.
"""
function Jadv(d::RelativisticDisc, modes::Dict{Int,H1Mode{T}}, k::Integer, m::Integer, r_p, ω0,
              h_val, hp; lmax::Integer, wave::Bool=false) where {T}
    a = T(d.a); r = modes[m].r_right[k]; gm = kerr_metric(a, r, T(π)/2)
    R = reconstruct(d, modes, k, m, r_p, ω0, h_val, hp; lmax=lmax)
    δU = wave ? R.δU_wave : R.δU
    δur = δU[2]; δuφ = sum(gm[4,ν]*δU[ν] for ν in 1:4)
    εp = eden0(d, r)*(1 + csound(d, r)^2)
    εp*real(δur*conj(δuφ))
end

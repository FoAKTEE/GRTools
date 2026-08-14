# =============================================================================================================
# WeylScalars — the GAUGE-INVARIANT linearized Weyl scalars ψ₀⁽¹⁾[h], ψ₄⁽¹⁾[h] REGENERATED from the CCK-
# reconstructed coordinate metric h (legacy MetricReconstruction/weyl_scalars.jl). ψ₀ and ψ₄ are gauge-invariant
# for a vacuum perturbation about a Petrov-type-D background (Kerr, Kinnersley tetrad: background ψ₀=ψ₁=ψ₃=ψ₄=0),
# so they are the natural CROSS-GAUGE quantity — and the H10 physics gate that validates the CCK preset CHOICE.
#
# Construction (vacuum background + vacuum perturbation ⇒ linearized Weyl = linearized Riemann):
#   ∇_μ h_νρ                                            (Ch)
#   δΓ^a_bc  = ½ g^{ad}(∇_b h_dc + ∇_c h_db − ∇_d h_bc)  (dG)
#   δR^a_bcd = ∇_c δΓ^a_bd − ∇_d δΓ^a_bc
#   δR_abcd = g_ae δR^e_bcd + h_ae R⁽⁰⁾^e_bcd
#   ψ₀ = δR_abcd l^a m^b l^c m^d ,   ψ₄ = δR_abcd n^a m̄^b n^c m̄^d
# Everything is assembled from the EXACT analytic (h,∂h,∂∂h) jets + analytic Christoffels/∂Christoffels (Einstein.jl,
# green) — no finite differences. By construction the CCK ORG reconstruction reproduces ψ₄ = (const)·R₋₂(r)·S₋₂(θ)·ρ⁴
# and IRG reproduces ψ₀ = (const)·R₊₂(r)·S₊₂(θ); the proportionality constant is (r,θ)-independent (H10 gate).
# Background g/g⁻¹/Γ and analytic ∂g DELEGATE to Spacetime (via Einstein). Type-generic {T} (BigFloat-inherited).
# =============================================================================================================

import Spacetime: kerr_metric, kerr_metric_inv, kinnersley_tetrad

export weyl_scalar, lin_riemann_lowered

# Metric compatibility gives ∂_α g^{ab} directly from g⁻¹ and Γ. Keeping this adapter independent of
# Kerr coordinates makes the curvature kernel testable with arbitrary analytic background jets.
function _dginv_from_connection(::Type{C}, gi, Γ) where {C}
    dgi = zeros(C, 4, 4, 4)   # dgi[α,a,b]
    @inbounds for α in 1:4, a_ in 1:4, b in 1:4
        acc = zero(C)
        for s in 1:4
            acc -= Γ[a_,α,s]*gi[s,b] + Γ[b,α,s]*gi[a_,s]
        end
        dgi[α,a_,b] = acc
    end
    dgi
end

"""
    lin_riemann_lowered(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=:Up) -> δR[a,b,c,d]  (all indices DOWN)

Linearized Riemann tensor δR_abcd of the reconstructed perturbation at the BL point, fully lowered. For a vacuum
reconstruction this equals the linearized Weyl tensor (consumed by `weyl_scalar`). Uses the EXACT analytic h-jets
(`_h_jets_analytic`, Projection RTJet) + analytic Christoffels/∂Christoffels (`christoffel_BL`,
`_dchristoffel_BL_analytic`, Einstein) — no finite differences.
"""
function lin_riemann_lowered(gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ, t, φ, weyl; bc::Symbol=:Up)
    g  = kerr_metric(a, r, θ; coords=:BL)
    gi = kerr_metric_inv(a, r, θ; coords=:BL)
    Γ  = christoffel_BL(a, r, θ)
    dΓ = _dchristoffel_BL_analytic(a, r, θ)
    h, dh, ddh = _h_jets_analytic(gauge, a, ω, l, m, r, θ, t, φ, weyl; bc=bc)
    _lin_riemann_lowered_from_jets(g, gi, Γ, dΓ, h, dh, ddh)
end

"All-down linearized Riemann kernel for an analytic background and perturbation 2-jet."
function _lin_riemann_lowered_from_jets(g, gi, Γ, dΓ, h, dh, ddh)
    C = eltype(h); gC = C.(g); giC = C.(gi)
    dgi = _dginv_from_connection(C, giC, Γ)
    # ∇_μ h_νρ
    Ch = zeros(C,4,4,4)
    @inbounds for μ in 1:4, ν in 1:4, ρ in 1:4
        acc = dh[μ,ν,ρ]
        for σ in 1:4; acc -= Γ[σ,μ,ν]*h[σ,ρ] + Γ[σ,μ,ρ]*h[ν,σ]; end
        Ch[μ,ν,ρ] = acc
    end
    # ∂_α(∇_μ h_νρ)
    dCh = zeros(C,4,4,4,4)
    @inbounds for α in 1:4, μ in 1:4, ν in 1:4, ρ in 1:4
        acc = ddh[α,μ,ν,ρ]
        for σ in 1:4
            acc -= dΓ[α,σ,μ,ν]*h[σ,ρ] + Γ[σ,μ,ν]*dh[α,σ,ρ]
            acc -= dΓ[α,σ,μ,ρ]*h[ν,σ] + Γ[σ,μ,ρ]*dh[α,ν,σ]
        end
        dCh[α,μ,ν,ρ] = acc
    end
    # δΓ^a_bc and ∂_α δΓ^a_bc
    dG  = zeros(C,4,4,4)      # dG[a,b,c]
    ddG = zeros(C,4,4,4,4)    # ddG[α,a,b,c]
    @inbounds for a_ in 1:4, b in 1:4, c in 1:4
        acc = zero(C)
        for d in 1:4; acc += giC[a_,d]*(Ch[b,d,c] + Ch[c,d,b] - Ch[d,b,c]); end
        dG[a_,b,c] = acc/2
    end
    @inbounds for α in 1:4, a_ in 1:4, b in 1:4, c in 1:4
        acc = zero(C)
        for d in 1:4
            acc += dgi[α,a_,d]*(Ch[b,d,c] + Ch[c,d,b] - Ch[d,b,c])
            acc += giC[a_,d]*(dCh[α,b,d,c] + dCh[α,c,d,b] - dCh[α,d,b,c])
        end
        ddG[α,a_,b,c] = acc/2
    end
    # ∇_μ δΓ^a_bc = ∂_μ δΓ + Γ^a_μσ δΓ^σ_bc − Γ^σ_μb δΓ^a_σc − Γ^σ_μc δΓ^a_bσ
    CdG = zeros(C,4,4,4,4)   # CdG[μ,a,b,c]
    @inbounds for μ in 1:4, a_ in 1:4, b in 1:4, c in 1:4
        acc = ddG[μ,a_,b,c]
        for σ in 1:4
            acc += Γ[a_,μ,σ]*dG[σ,b,c] - Γ[σ,μ,b]*dG[a_,σ,c] - Γ[σ,μ,c]*dG[a_,b,σ]
        end
        CdG[μ,a_,b,c] = acc
    end
    # Fully lowered variation:
    #   δR_abcd = g_ae δR^e_bcd + h_ae R⁽⁰⁾^e_bcd.
    # The second term is the variation of the lowering metric.  Omitting it breaks the
    # first-pair antisymmetry R_abcd = -R_bacd on curved backgrounds.
    Rout = zeros(C,4,4,4,4)
    @inbounds for a_ in 1:4, b in 1:4, c in 1:4, d in 1:4
        acc = zero(C)
        for e in 1:4
            δRu = CdG[c,e,b,d] - CdG[d,e,b,c]
            R0u = dΓ[c,e,d,b] - dΓ[d,e,c,b]
            for s in 1:4
                R0u += Γ[e,c,s]*Γ[s,d,b] - Γ[e,d,s]*Γ[s,c,b]
            end
            acc += gC[a_,e]*δRu + h[a_,e]*R0u
        end
        Rout[a_,b,c,d] = acc
    end
    Rout
end

"Repeated-null Riemann/Weyl contraction used by the Newman–Penrose ψ₀/ψ₄ definitions."
function _repeated_null_contraction(Rt, u1, u2)
    C = eltype(Rt)
    ψ = zero(C)
    @inbounds for a_ in 1:4, b in 1:4, c in 1:4, d in 1:4
        ψ += Rt[a_,b,c,d]*u1[a_]*u2[b]*u1[c]*u2[d]
    end
    ψ
end

"""
    weyl_scalar(which, gauge, a, ω, l, m, r, θ; weyl, t=0, φ=0, bc=:Up) -> Complex

Gauge-invariant linearized Weyl scalar of the reconstructed metric. `which` ∈ (:psi0, :psi4).
  ψ₀ = δR_abcd l^a m^b l^c m^d        (peels the s=+2 Teukolsky content)
  ψ₄ = δR_abcd n^a m̄^b n^c m̄^d       (peels the s=−2 Teukolsky content)
Kinnersley tetrad up-legs (l,n,m,m̄) DELEGATED to Spacetime. For a CCK reconstruction, ψ₀[h_IRG] ∝ R₊₂·S₊₂ and
ψ₄[h_ORG] ∝ R₋₂·S₋₂·ρ⁴ — gauge-invariant, hence the Teukolsky-shaped H10 physics gate. Default `weyl` is
`weyl_amplitudes(:unit,…)`; pass `:psi0_IRG`/`:psi4_ORG` presets for the value-grounded reconstruction.
"""
function weyl_scalar(which::Symbol, gauge::Symbol, a, ω, l::Integer, m::Integer, r, θ;
                     weyl=nothing, t=0, φ=0, bc::Symbol=:Up)
    w = weyl === nothing ? weyl_amplitudes(:unit, ω, l, m, a) : weyl
    Rt = lin_riemann_lowered(gauge, a, ω, l, m, r, θ, t, φ, w; bc=bc)
    legs = kinnersley_tetrad(a, r, θ; coords=:BL)   # (l, n, m, m̄) up-vectors
    C = eltype(Rt)
    lup = C.(legs[1]); nup = C.(legs[2]); mup = C.(legs[3]); mbup = C.(legs[4])
    if which === :psi0
        u1, u2 = lup, mup                # l^a m^b l^c m^d
    elseif which === :psi4
        u1, u2 = nup, mbup               # n^a m̄^b n^c m̄^d
    else
        error("weyl_scalar: which must be :psi0 or :psi4")
    end
    _repeated_null_contraction(Rt, u1, u2)
end

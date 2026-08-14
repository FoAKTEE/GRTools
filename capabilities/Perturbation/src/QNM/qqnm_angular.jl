# QNM/qqnm_angular.jl — the general-a angular reduction: the four Paper-I:591 angular structures
#     𝒜₁ = [L₂†S]²,  𝒜₂ = S·L₁†L₂†S,  𝒜₃ = 3ia sinθ·S·L₂†S,  𝒜₄ = (a sinθ·S)²
# over complex-γ spheroidal values + θ-derivatives, projected onto the child spheroidal
# ₋₂S_{l_Q m_Q}(θ; γ_Q) at γ_Q = aω_Q = a(ω₁+ω₂) by Gauss–Legendre quadrature (Paper I:723,840).
#
# ─────────────────────────────── THE L†_n LADDER (Paper I:218, 232-237) ─────────────────────────────
# The Chandrasekhar angular ladder operator, VERBATIM (Paper I:232-237):
#     L†_n = ∂_θ − P + n cotθ ,   P = −aω sinθ + m cscθ      ⟹    L†_n = ∂_θ + γ sinθ − m cscθ + n cotθ
# with γ = aω the parent spheroidicity (COMPLEX for a QNM ω).  Acts on S = ₋₂S_{l_L m_L}(θ; γ_L),
# carrying the PARENT (m_L, γ_L=aω_L).  Verified against the a=0 Schwarzschild reduction (Paper I:293,
# 749):  L₂†(₋₂Y_{lm}) = −ð(₋₂Y) = −√((l+2)(l−1)) ₋₁Y  and  L₁†L₂†(₋₂Y) = √((l+2)(l+1)l(l−1)) ₀Y
# (test: 10-digit match at a=0).  As first-order θ-operators on S with S0,S1,S2 = S, ∂_θS, ∂²_θS:
#     L₂†S      = S1 + P₂·S0 ,                         P₂ = γ sinθ − m cscθ + 2 cotθ
#     ∂_θ(L₂†S) = S2 + P₂·S1 + P₂′·S0 ,                P₂′ = γ cosθ + m cscθ cotθ − 2 csc²θ
#     L₁†L₂†S   = ∂_θ(L₂†S) + P₁·(L₂†S) ,              P₁ = γ sinθ − m cscθ + 1 cotθ
#
# PRECISION TIER: F64 angular × {T}/BigFloat radial — the Harmonics value tier is Float64-capped
# by design, so the joint projection carries {T} precision with an F64 angular ceiling.
#
# ─────────────────────────────── CHILD PROJECTION ORTHOGONALITY ─────────────────────────────────────
# Complex-γ spheroidals are COMPLEX-symmetric (not Hermitian) ⇒ bi-orthogonal WITHOUT conjugation:
#     ∫ ₋₂S_{l m}(θ;γ) ₋₂S_{l′ m}(θ;γ) sinθ dθ ∝ δ_{ll′}   (NO c.c.; Paper I:715).
# So the projection sum uses S_Q, not conj(S_Q).  GL rule on y=−cosθ, dy=sinθ dθ (SphereCollocation).

# ═══════════════════════════════ the four angular structures ══════════════════════════════════════
"""
    qqnm_angular_structures(l_L, m_L, γ_L, a, θ) -> (𝒜₁,𝒜₂,𝒜₃,𝒜₄)

The four Paper-I:591 angular structures at colatitude `θ` for the parent spheroidal
`S = ₋₂S_{l_L m_L}(·; γ_L)` (`γ_L = aω_L`, complex).  `𝒜₁=[L₂†S]²`, `𝒜₂=S·L₁†L₂†S`,
`𝒜₃=3ia sinθ·S·L₂†S`, `𝒜₄=(a sinθ·S)²` (`𝒜₃,𝒜₄` carry explicit `a`, vanish at `a=0`).  Built from
the complex-γ spheroidal value + its ∂_θ,∂²_θ; the L†_n ladders act as first-order θ-operators
(file header).  Returns `ComplexF64`.
"""
function qqnm_angular_structures(l_L::Integer, m_L::Integer, γ_L::Number, a::Real, θ::Real)
    S = SpinWeightedSpheroidalHarmonicS(-2, l_L, m_L, complex(γ_L); method=:SphericalExpansion)
    θf = Float64(θ); γ = ComplexF64(γ_L); mL = Float64(m_L); af = Float64(a)
    S0 = swsh_derivative(S, 0, 0, θf, 0.0)
    S1 = swsh_derivative(S, 1, 0, θf, 0.0)
    S2 = swsh_derivative(S, 2, 0, θf, 0.0)
    st = sin(θf); ct = cos(θf); csc = 1 / st; cot = ct / st
    P2  = γ * st - mL * csc + 2 * cot
    P2p = γ * ct + mL * csc * cot - 2 * csc^2
    P1  = γ * st - mL * csc + 1 * cot
    L2S   = S1 + P2 * S0
    dL2S  = S2 + P2 * S1 + P2p * S0
    L1L2S = dL2S + P1 * L2S
    𝒜1 = L2S^2
    𝒜2 = S0 * L1L2S
    𝒜3 = 3im * af * st * S0 * L2S
    𝒜4 = (af * st)^2 * S0^2
    return (𝒜1, 𝒜2, 𝒜3, 𝒜4)
end

"""
    qqnm_angular_projections(l_L, m_L, l_Q, m_Q, ω_L, a; ny=96) -> (Π₁,Π₂,Π₃,Π₄)

The four child-projection coefficients `Π_ν = ∫ 𝒜_ν(θ) S_Q(θ) sinθ dθ / ∫ S_Q(θ)² sinθ dθ`, with the
child spheroidal `S_Q = ₋₂S_{l_Q m_Q}(·; γ_Q)` at `γ_Q = a·ω_Q = a·(2ω_L)` (Paper I:1123 — the
quadratic-mode basis, DISTINCT from the linear child).  Bi-orthogonal projection (NO c.c.), GL
quadrature over `ny` points.  These are the general-a analogue of the a=0 3-j pair: as `a→0`,
`Π₁ → (l_L+2)(l_L−1)·P_B` and `Π₂ → √((l_L+2)(l_L+1)l_L(l_L−1))·P_A` (the certified
`_qqnm_proj_a0_analytic`), while `Π₃,Π₄ = O(a),O(a²) → 0`.  ComplexF64 (F64 value tier).
"""
function qqnm_angular_projections(l_L::Integer, m_L::Integer, l_Q::Integer, m_Q::Integer,
                                  ω_L::Number, a::Real; ny::Integer=96)
    γ_L = a * ω_L
    ω_Q = 2 * ω_L
    γ_Q = a * ω_Q
    yv, wv = gl_nodes_weights(Float64, ny)                 # y=−cosθ, dy=sinθdθ
    SQ = SpinWeightedSpheroidalHarmonicS(-2, l_Q, m_Q, complex(γ_Q); method=:SphericalExpansion)
    num = zeros(ComplexF64, 4); den = zero(ComplexF64)
    for j in 1:ny
        θ = acos(-yv[j])
        A = qqnm_angular_structures(l_L, m_L, γ_L, a, θ)
        sq = swsh_derivative(SQ, 0, 0, θ, 0.0)
        @inbounds for ν in 1:4
            num[ν] += wv[j] * A[ν] * sq
        end
        den += wv[j] * sq^2
    end
    return (num[1] / den, num[2] / den, num[3] / den, num[4] / den)
end

# ═══════════════════════════════ general-a joint (Γ(θ)-coupled) source assembly ═══════════════════
# For a≠0 the radial coefficients a_ij^(ν)(Γ,Γ̄,r) depend on θ through Γ=r+ia cosθ, Γ̄=r−ia cosθ, so
# the child projection does NOT factorize into (radial contraction)×(angular projection) as it does at
# a=0.  The joint θ-integral is
#     W_ab(r) = −(Δ⁶/16) · Σ_ν ∫ 𝒜_ν(θ) · C_ab^(ν)(Γ(θ),Γ̄(θ),r) · S_Q(θ) sinθ dθ / ∫ S_Q(θ)² sinθ dθ,
# where C_ab^(ν) is the (α_k,β_k)-folded radial contraction of a_ij^(ν) (identical to the a=0
# `_qqnm_Cab_pp`, now with Γ=Γ(θ)).  The parent chains (α_k,β_k) are θ-INDEPENDENT (X_i=D₂†ⁱℛ radial),
# computed ONCE per r-node.  As a→0: Γ(θ)→r, 𝒜₃,𝒜₄→0, and W_ab → the certified a=0 coefficients.
"""
    qqnm_source_coeffs_gl(parent1, parent2, l_Q; a, r, T=ComplexF64, ny=96)

The general-a (`a≠0`) ++ source coefficients `(W00,W01,W10,W11)` of `Q(r)=Σ W_ab ℛ₁^{(a)}ℛ₂^{(b)}`
via the Γ(θ)-coupled joint child projection (file header).  Reduces to `qqnm_source_coeffs(...;a=0)`
as `a→0` (continuity).  Parent chains at BigFloat if `T=Complex{BigFloat}` (radial tier); the angular
structures + child spheroidal stay ComplexF64 (F64 value tier).  Δ = r²−2r+a².
"""
function qqnm_source_coeffs_gl(parent1, parent2, l_Q::Integer;
                               a::Real, r, T::Type=ComplexF64, ny::Integer=96)
    C = T <: Complex ? T : Complex{T}
    rs = r isa AbstractVector ? r : [r]
    m_Q = parent1.m + parent2.m
    identical = (parent1 === parent2) ||
                (parent1.l == parent2.l && parent1.m == parent2.m &&
                 parent1.n == parent2.n && parent1.ω == parent2.ω)
    ω_L = parent1.ω
    γ_L = a * ω_L
    γ_Q = a * (parent1.ω + parent2.ω)
    yv, wv = gl_nodes_weights(Float64, ny)
    SQ = SpinWeightedSpheroidalHarmonicS(-2, l_Q, m_Q, complex(γ_Q); method=:SphericalExpansion)
    # angular structures + child value on the GL grid (θ-only, r-independent) — precompute once
    θv = [acos(-yv[j]) for j in 1:ny]
    Amat = [qqnm_angular_structures(parent1.l, parent1.m, γ_L, a, θv[j]) for j in 1:ny]  # ν-tuples
    sqv  = ComplexF64[swsh_derivative(SQ, 0, 0, θv[j], 0.0) for j in 1:ny]
    den  = zero(ComplexF64)
    for j in 1:ny
        den += wv[j] * sqv[j]^2
    end
    W00 = zeros(C, length(rs)); W01 = zeros(C, length(rs))
    W10 = zeros(C, length(rs)); W11 = zeros(C, length(rs))
    for (k, rr) in enumerate(rs)
        rrC = C(rr)
        (α1, β1) = _qqnm_parent_chain(parent1, a, rr, C)
        (α2, β2) = identical ? (α1, β1) : _qqnm_parent_chain(parent2, a, rr, C)
        Δ = C(rr^2 - 2rr + a^2); pref = -(Δ^6) / 16
        acc00 = zero(C); acc01 = zero(C); acc10 = zero(C); acc11 = zero(C)
        for j in 1:ny
            Γ  = rrC + im * C(a) * C(cos(θv[j]))           # Γ = r + ia cosθ
            Γb = rrC - im * C(a) * C(cos(θv[j]))           # Γ̄ = r − ia cosθ
            wj = C(wv[j]) * C(sqv[j])                       # GL weight × child value (bi-orthogonal)
            for ν in 1:4
                (C00, C01, C10, C11) = _qqnm_Cab_pp(ν, α1, β1, α2, β2, Γ, Γb, rr, C)
                Aν = C(Amat[j][ν])
                acc00 += wj * Aν * C00; acc01 += wj * Aν * C01
                acc10 += wj * Aν * C10; acc11 += wj * Aν * C11
            end
        end
        denC = C(den)
        W00[k] = pref * acc00 / denC; W01[k] = pref * acc01 / denC
        W10[k] = pref * acc10 / denC; W11[k] = pref * acc11 / denC
    end
    if identical
        return (W00=W00, W01=(W01 .+ W10), W11=W11, r=collect(rs), m_Q=m_Q)
    end
    return (W00=W00, W01=W01, W10=W10, W11=W11, r=collect(rs), m_Q=m_Q)
end

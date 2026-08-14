# QNM/qqnm_source.jl — the second-order QQNM source Q(r) of the projected equation
#   (Δ𝒟_{-1}^{†(2)}𝒟_0^{(2)} − 12iω_L r − λ_{l_Q}^{(2)})R_{l_Q}^{(2)}(r) = Q(r)
# (Paper I:722; owner:1198), from the baked coefficient tables, in three stages:
#  (1) radial quotient fold: A_ν = Σ a_ij^(ν) X_i X_j (X_i = D₂†ⁱℛ) reduced via the spin-+2 field
#      equation to A_ν = C₀₀ℛ² + C₀₁ℛℛ′ + C₁₁ℛ′², X_k = α_kℛ + β_kℛ′.  CONSTRAINT: the fold
#      eigenvalue is λ₊ = λ_Teuk(s=+2) + 4, not λ_Teuk (the +4 from the D₂† ladder; the wrong
#      choice fails the Paper-I:1427 transcription).
#  (2) angular structures (Paper I:591) over the parent spheroidal at γ=aω_L; at a=0 the ladders
#      collapse to L₂†S=−√((l+2)(l−1))₋₁Y, L₁†L₂†S=√((l+2)(l+1)l(l−1))₀Y (Paper I:749).
#  (3) child projection (Paper I:723; owner:1204) with the (16/Δ⁶) bookkeeping; child basis at
#      aω_Q = a(ω₁+ω₂) (Paper I:1123). Routes: a=0 analytic 3-j; general-a Gauss–Legendre.
# Channels: ++ symmetric fold; +− ordered bilinear with the 𝔠_- mirror carry; −− ≡ 0.
# {T}/BigFloat-inherited (the 3-j route via exact BigInt factorials).

# ═══════════════════════════════ STAGE 1 — radial quotient fold ════════════════════════════════════
# App-B coefficient jets U,V,T,F,G as length-LEN r-jets at r. λ is the pinned λ₊ (=λ_Teuk+4).
# {T}-generic via NumericsCore jet ops.
function _qqnm_coefjets(::Type{C}, a, r, ω, m, λ, LEN::Int) where {C}
    K = zeros(C, LEN); K[1] = (r^2 + a^2)*ω - a*m; K[2] = 2r*ω; K[3] = C(2ω)
    Δ = zeros(C, LEN); Δ[1] = r^2 - 2r + a^2;      Δ[2] = C(2r - 2); Δ[3] = C(2)
    rm1 = zeros(C, LEN); rm1[1] = C(r - 1); rm1[2] = one(C)
    invΔ = jinv(Δ)
    U = jmul(( im .* K) .+ (4 .* rm1), invΔ)          # D₂† coefficient (F7-3a double-entry certified)
    V = jmul((-im .* K) .+ (2 .* rm1), invΔ)          # D₁  coefficient
    Tn = zeros(C, LEN); Tn[1] = 6*im*ω*r - λ; Tn[2] = C(6*im*ω)
    T = jmul(Tn, invΔ)                                 # (6iωr − λ₊)/Δ
    F = (-6) .* jmul(rm1, invΔ)                        # −6(r−1)/Δ  (R″=FR′+GR)
    G = (-1) .* jderiv(U) .- jmul(V, U) .- T
    (U=U, V=V, T=T, F=F, G=G)
end

# The (α_k,β_k) chain: X_k = P^k ℛ = α_k ℛ + β_k ℛ′.
# α_{k+1}=α_k′+Uα_k+Gβ_k, β_{k+1}=α_k+β_k′−Vβ_k; (α_0,β_0)=(1,0).  Returns the VALUES α_k(r),β_k(r)
# (jet [1] component) as vectors indexed 0..kmax.
function _qqnm_alphabeta_vals(cj, LEN::Int, kmax::Int, ::Type{C}) where {C}
    α = Vector{Vector{C}}(undef, kmax+1); β = Vector{Vector{C}}(undef, kmax+1)
    α[1] = zeros(C, LEN); α[1][1] = one(C); β[1] = zeros(C, LEN)
    for k in 1:kmax
        α[k+1] = jderiv(α[k]) .+ jmul(cj.U, α[k]) .+ jmul(cj.G, β[k])
        β[k+1] = α[k] .+ jderiv(β[k]) .- jmul(cj.V, β[k])
    end
    (C[α[k][1] for k in 1:kmax+1], C[β[k][1] for k in 1:kmax+1])
end

"""
    qqnm_lambda_plus(l, m, a, ω) -> Complex

The fold eigenvalue `λ₊ = λ_Teuk(s=+2,l,m,aω) + 4` for the spin-+2 Hertz radial (file-header
constraint).  `λ_Teuk` is `Harmonics.SpinWeightedSpheroidalEigenvalue`.  {T}/BigFloat-inherited.
"""
function qqnm_lambda_plus(l::Integer, m::Integer, a::Real, ω::Number)
    γ = complex(a * ω)
    # a=0 (γ=0): the spherical eigenvalue is exactly l(l+1)−s(s+1), s=+2 (avoids any eigensolve).
    iszero(γ) && return (l*(l+1) - 2*(2+1)) + 4
    return SpinWeightedSpheroidalEigenvalue(2, l, m, γ) + 4
end

# Per-parent radial chain (α_k,β_k values at r).  `parent = (l=,m=,n=,ω=)`.  s=+2 Hertz radial.
function _qqnm_parent_chain(parent, a, r, ::Type{C}; kmax::Int=6, LEN::Int=8) where {C}
    λ = C(qqnm_lambda_plus(parent.l, parent.m, real(a), parent.ω))
    cj = _qqnm_coefjets(C, C(a), C(r), C(parent.ω), parent.m, λ, LEN)
    _qqnm_alphabeta_vals(cj, LEN, kmax, C)
end

# Mirror-channel chain for the +− `−` factor (owner:1263 eq:app-Cminus-channel-factor):
# 𝔠_-: (ω_+,m_+,λ_+,Γ,Γ̄,i) ↦ (ω_-,m_-,λ_-,Γ̄,Γ,−i), ω_-=−conj(ω_{L̄}), m_-=−m_{L̄}.  The i↦−i flips
# the sign of the iK term in U,V (= conjugating the operator's explicit i, NOT the fields).
function _qqnm_parent_chain_mirror(parent, a, r, ::Type{C}; kmax::Int=6, LEN::Int=8) where {C}
    ω_m = -conj(C(parent.ω)); m_m = -parent.m
    λ_m = C(qqnm_lambda_plus(parent.l, m_m, real(a), ω_m))
    # i↦−i: build the coef jets with conjugated explicit i by negating the imaginary unit in U,V,T.
    K = zeros(C, LEN); K[1] = (r^2 + a^2)*ω_m - a*m_m; K[2] = 2r*ω_m; K[3] = C(2ω_m)
    Δ = zeros(C, LEN); Δ[1] = r^2 - 2r + a^2; Δ[2] = C(2r - 2); Δ[3] = C(2)
    rm1 = zeros(C, LEN); rm1[1] = C(r - 1); rm1[2] = one(C)
    invΔ = jinv(Δ); mi = -im
    U = jmul((mi .* K) .+ (4 .* rm1), invΔ)
    V = jmul((-mi .* K) .+ (2 .* rm1), invΔ)
    Tn = zeros(C, LEN); Tn[1] = 6*mi*ω_m*r - λ_m; Tn[2] = C(6*mi*ω_m)
    T = jmul(Tn, invΔ)
    G = (-1) .* jderiv(U) .- jmul(V, U) .- T
    cj = (U=U, V=V, T=T, F=(-6) .* jmul(rm1, invΔ), G=G)
    _qqnm_alphabeta_vals(cj, LEN, kmax, C)
end

# ── ++ C_ab: general-parent-pair bilinear symmetrization of the a_ij contraction ────────────────────
# A_ν = Σ_i a_ii X_i^1 X_i^2 + Σ_{i<j}(a_ij/2)(X_i^1 X_j^2 + X_j^1 X_i^2),  X_k^p = α_k^p ℛ_p + β_k^p ℛ_p′.
# ⇒ A_ν = C₀₀ ℛ₁ℛ₂ + C₀₁ ℛ₁ℛ₂′ + C₁₀ ℛ₁′ℛ₂ + C₁₁ ℛ₁′ℛ₂′.  Identical parents ⇒ printed Paper-I:591
# form (C₀₁=C₁₀ merge, upper-triangular Σ a_ij α_iα_j …).
function _qqnm_Cab_pp(ν::Integer, α1, β1, α2, β2, Γ, Γb, r, ::Type{C}) where {C}
    aij = qqnm_pp_aij(ν, C, Γ, Γb, r)
    C00 = zero(C); C01 = zero(C); C10 = zero(C); C11 = zero(C)
    for ((i,j), a) in aij
        if i == j
            C00 += a*α1[i+1]*α2[i+1]; C01 += a*α1[i+1]*β2[i+1]
            C10 += a*β1[i+1]*α2[i+1]; C11 += a*β1[i+1]*β2[i+1]
        else
            h = a/2
            C00 += h*(α1[i+1]*α2[j+1] + α1[j+1]*α2[i+1])
            C01 += h*(α1[i+1]*β2[j+1] + α1[j+1]*β2[i+1])
            C10 += h*(β1[i+1]*α2[j+1] + β1[j+1]*α2[i+1])
            C11 += h*(β1[i+1]*β2[j+1] + β1[j+1]*β2[i+1])
        end
    end
    (C00, C01, C10, C11)
end

# ── +− C_ab per mixed atom μ: ordered bilinear (owner:1560) ─────────────────────────────────────────
# B_μ = Σ_{i,j} b_ij^(μ) X_i^+ X_j^-  ⇒  C₀₀ℛ₊ℛ₋ + C₀₁ℛ₊ℛ₋′ + C₁₀ℛ₊′ℛ₋ + C₁₁ℛ₊′ℛ₋′.
function _qqnm_Cab_pm(μ::Integer, αp, βp, αm, βm, Γ, Γb, r, ::Type{C}) where {C}
    bij = qqnm_pm_bij(μ, C, Γ, Γb, r)
    C00 = zero(C); C01 = zero(C); C10 = zero(C); C11 = zero(C)
    for ((i,j), b) in bij
        C00 += b*αp[i+1]*αm[j+1]; C01 += b*αp[i+1]*βm[j+1]
        C10 += b*βp[i+1]*αm[j+1]; C11 += b*βp[i+1]*βm[j+1]
    end
    (C00, C01, C10, C11)
end

# ═══════════════════════════════ STAGE 3 — angular projection (a=0 analytic 3-j) ═══════════════════
# Wigner 3-j (Racah single-sum), exact via BigInt factorials, cast to T.  {T}/BigFloat-inherited.
function _wigner3j(::Type{T}, j1, j2, j3, m1, m2, m3) where {T}
    (m1 + m2 + m3 == 0) || return zero(T)
    (abs(m1) ≤ j1 && abs(m2) ≤ j2 && abs(m3) ≤ j3) || return zero(T)
    (j3 ≥ abs(j1 - j2) && j3 ≤ j1 + j2) || return zero(T)
    fac(n) = factorial(big(n))
    Δtri = fac(j1+j2-j3)*fac(j1-j2+j3)*fac(-j1+j2+j3)//fac(j1+j2+j3+1)
    pref = fac(j1+m1)*fac(j1-m1)*fac(j2+m2)*fac(j2-m2)*fac(j3+m3)*fac(j3-m3)
    kmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    kmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    s = big(0)//1
    for k in kmin:kmax
        den = fac(k)*fac(j1+j2-j3-k)*fac(j1-m1-k)*fac(j2+m2-k)*fac(j3-j2+m1+k)*fac(j3-j1-m2+k)
        s += (iseven(k) ? 1 : -1) * (big(1)//den)
    end
    sq = Δtri * pref            # rational; the 3-j is (−1)^{j1−j2−m3} √sq · s
    val = T(s) * sqrt(T(sq.num) / T(sq.den))
    return (iseven(j1 - j2 - m3) ? 1 : -1) * val
end

# The two Schwarzschild angular projections (Paper I:793 eq:3j_symbol_A, :809 eq:3j_symbol_B),
# = P_A, P_B in the unit-normalized child convention (∫S²_{l_Q}sinθdθ cancels).
function _qqnm_proj_a0_analytic(::Type{T}, l_L, m_L, l_Q, m_Q) where {T}
    common = T(-1)^m_Q * sqrt(T((2l_L+1)^2 * (2l_Q+1)) / (4*T(pi))) *
             _wigner3j(T, l_L, l_L, l_Q, m_L, m_L, -m_Q)
    P_A = common * _wigner3j(T, l_L, l_L, l_Q, 2, 0, -2)     # ₋₂Y₀Y   ← pairs with A₂
    P_B = common * _wigner3j(T, l_L, l_L, l_Q, 1, 1, -2)     # ₋₁Y²    ← pairs with A₁
    (P_A, P_B)
end

# GL-quadrature route (a=0): the SAME projections by Gauss–Legendre over the spin-weighted spherical
# θ-functions on the SphereCollocation grid.  dy = sinθ dθ ⇒ `inner_product` IS the θ-measure.
_sph_theta(::Type{T}, s, l, m, θ) where {T} =
    l < max(abs(m), abs(s)) ? zero(T) :        # nonexistent harmonic ≡ 0 (matches the 3-j guard)
    real(SpinWeightedSphericalHarmonicY(s, l, m, T(θ), zero(T)))
function _qqnm_proj_a0_gl(::Type{T}, l_L, m_L, l_Q, m_Q, ny::Integer) where {T}
    θv = acos.(cos_vals(T, ny))
    y1 = T[_sph_theta(T, -1, l_L, m_L, θ) for θ in θv]
    y2m2 = T[_sph_theta(T, -2, l_L, m_L, θ) for θ in θv]
    y0 = T[_sph_theta(T, 0, l_L, m_L, θ) for θ in θv]
    yQ = T[_sph_theta(T, -2, l_Q, m_Q, θ) for θ in θv]
    den = inner_product(yQ, yQ)
    P_B = inner_product(y1 .* y1, yQ) / den
    P_A = inner_product(y2m2 .* y0, yQ) / den
    (P_A, P_B)
end

# ═══════════════════════ +− CHILD PROJECTION (a=0) — the per-atom angular reduction [Q4-H-2] ════════
# The +− source is  Ŝ_{+-} = e^{i((m₁+m₂)φ−ω_Q t)} Σ_{μ=1}^{16} B_μ 𝒵_μ  (owner:1536), with the 16 mixed
# angular atoms 𝒵_μ (owner:1516, `QQNM_PM_ATOMS`) built from the + Weyl harmonic S_A=₋₂S_{l₁m₁}, the −
# reflected CCK metric slots  Ym0=S_B , Ym1=L₂†S_B , Ym2=L₁†L₂†S_B  (S_B=₋₂S_{l₂m₂}; the L†_n
# Chandrasekhar ladders — PINNED by the ++ hll↔`derive-spin2-plus-minus` metric-slot analogy
# hll=xm[0](Γ̄²/2·Ym2 − iΓ̄ s_a·Ym1)  vs  ++ hll=(Γ²/2)S2 + iΓ s_a·S1 with S1=L₂†S, S2=L₁†L₂†S) and the
# LOCAL angular derivative ℒ₀=∂_θ.  The child projection is the SAME bi-orthogonal GL θ-integral as the
# ++ round (`_qqnm_proj_a0_gl`, validated ≡ the analytic 3-j cA1,cA2 to machine ε).
#
# THE a→0 LIMIT.  Reorganize Ŝ_{+-}=Σ_{(i,j)} X_i^+ X_j^- g_{ij}(θ), g_{ij}(θ)=Σ_μ b_{ij}^{(μ)} 𝒵_μ(θ).
# At a=0 the s_a^{+k} atoms vanish and the s_a^{−k} atoms have b_{ij}∼a^k, so the joint limit is
# FINITE (=cot^kθ·radial).  CONSTRAINT: this joint per-atom projection is what makes the +− RHS
# horizon-analytic — a flat atom sum leaves a genuine (r−2)⁻² pole.  Evaluated at a=`aeps` (tiny).
# {T}/BigFloat radial; F64 angular value tier.
function _qqnm_pm_child_gij(parent1, parent2, l_Q::Integer, rs::AbstractVector, ::Type{C};
                            ny::Integer=200, aeps=nothing) where {C}
    RT = real(C)
    l1 = parent1.l; m1 = parent1.m; l2 = parent2.l; m2 = parent2.m; m_Q = m1 + m2
    a = aeps === nothing ? RT(10)^(-40) : RT(aeps)         # tiny a: the analytic s_a^{−k} limit
    yv, wv = gl_nodes_weights(Float64, ny)                 # y=−cosθ, dy=sinθ dθ
    θv = acos.(-yv)
    SA = SpinWeightedSpheroidalHarmonicS(-2, l1, m1, complex(0.0); method=:SphericalExpansion)
    SB = SpinWeightedSpheroidalHarmonicS(-2, l2, m2, complex(0.0); method=:SphericalExpansion)
    SQ = SpinWeightedSpheroidalHarmonicS(-2, l_Q, m_Q, complex(0.0); method=:SphericalExpansion)
    dA(k,θ) = real(swsh_derivative(SA, k, 0, θ, 0.0))
    dB(k,θ) = real(swsh_derivative(SB, k, 0, θ, 0.0))
    # − CCK ladder tower on S_B (L†_n = ∂_θ − m₂cscθ + n cotθ at a=0): Ym1=L₂†S_B, Ym2=L₁†L₂†S_B
    L2(θ)   = (c=1/sin(θ); ct=cos(θ)/sin(θ); dB(1,θ) - m2*c*dB(0,θ) + 2*ct*dB(0,θ))
    dL2(θ)  = (c=1/sin(θ); ct=cos(θ)/sin(θ); dB(2,θ) - m2*(c*dB(1,θ) - c*ct*dB(0,θ)) + 2*(ct*dB(1,θ) - c^2*dB(0,θ)))
    L1L2(θ) = (c=1/sin(θ); ct=cos(θ)/sin(θ); dL2(θ) - m2*c*L2(θ) + ct*L2(θ))
    yQ = Float64[real(swsh_derivative(SQ, 0, 0, θ, 0.0)) for θ in θv]
    den = inner_product(yQ, yQ)
    # the 16 atom angular words 𝒵word_μ(θ) (s_a-power stripped) + their s_a power (QQNM_PM_ATOMS)
    zword = Vector{Function}(undef, 16)
    zword[1]  = θ-> dA(1,θ)*dB(1,θ) + dA(0,θ)*dB(2,θ)      # ∂θ(S_A ∂θ S_B)
    zword[2]  = θ-> dA(1,θ)*L2(θ)   + dA(0,θ)*dL2(θ)       # ∂θ(S_A L₂†S_B)
    zword[3]  = θ-> dB(1,θ)*dA(1,θ) + dB(0,θ)*dA(2,θ)      # ∂θ(S_B ∂θ S_A)
    zword[4]  = θ-> dA(0,θ)*dL2(θ)                          # S_A ∂θ(L₂†S_B)
    zword[5]  = θ-> dA(0,θ)*dB(0,θ)                          # S_A S_B
    zword[6]  = θ-> dA(0,θ)*L1L2(θ)                         # S_A L₁†L₂†S_B
    zword[7]  = θ-> dA(1,θ)*dB(0,θ) + dA(0,θ)*dB(1,θ)       # ∂θ(S_A S_B)  [s_a^{-1}]
    zword[8]  = zword[7]                                     # [s_a^{+1}] → 0
    zword[9]  = θ-> dA(0,θ)*dB(1,θ)                          # S_A ∂θ S_B  [s_a^{-1}]
    zword[10] = zword[9]                                     # [s_a^{+1}] → 0
    zword[11] = θ-> dA(0,θ)*L2(θ)                            # S_A L₂†S_B  [s_a^{-1}]
    zword[12] = zword[11]                                    # [s_a^{+1}] → 0
    zword[13] = θ-> dB(0,θ)*dA(1,θ)                          # S_B ∂θ S_A  [s_a^{-1}]
    zword[14] = zword[13]                                    # [s_a^{+1}] → 0
    zword[15] = zword[5]                                     # S_A S_B     [s_a^{-2}]
    zword[16] = zword[5]                                     # [s_a^{+2}] → 0
    spow = [a.sapow for a in QQNM_PM_ATOMS]
    # precompute the F64 angular products 𝒵word_μ(θ)·yQ(θ)·w(θ) per atom
    zθ = [Float64[zword[μ](θv[j]) for j in 1:ny] for μ in 1:16]
    # g_{ij}(r) = Σ_μ Σ_θ w·[b_{ij}^{(μ)}(Γ(θ),Γ̄(θ),r)·(a sinθ)^{spow}]·𝒵word_μ·yQ / den
    pairs = QQNM_PM_SUPPORT
    G = Dict{Tuple{Int,Int},Vector{C}}(pr => zeros(C, length(rs)) for pr in pairs)
    for (k, rr) in enumerate(rs)
        rC = C(rr)
        for μ in 1:16
            p = spow[μ]
            for (jθ, θ) in enumerate(θv)
                cθ = cos(θ); sθ = sin(θ)
                Γ = rC + im*C(a)*C(cθ); Γb = rC - im*C(a)*C(cθ)
                bij = qqnm_pm_bij(μ, C, Γ, Γb, rC)
                isempty(bij) && continue
                sa_p = (C(a)*C(sθ))^p
                wz = C(wv[jθ]) * C(zθ[μ][jθ]) * C(yQ[jθ]) * sa_p
                for (pr, b) in bij
                    G[pr][k] += b * wz
                end
            end
        end
        for pr in pairs; G[pr][k] /= C(den); end
    end
    G
end

# ═══════════════════════════════ ASSEMBLY — the source coefficients W_ab(r) ════════════════════════
"""
    qqnm_source_coeffs(parent1, parent2, channel, l_Q; a, r, T=ComplexF64, ny=64)

The assembled radial source coefficients of `Q(r)` for the excitation channel, on the r-point(s) `r`
(scalar or vector).  Returns a NamedTuple:
  * `channel=:pp` (or `:plusplus`): `(W00, W01, W11)` with `Q(r) = W00 ℛ² + W01 ℛℛ′ + W11 ℛ′²`
    (identical parents) — or `(W00, W01, W10, W11)` for a distinct-parent pair
    (`Q = W00 ℛ₁ℛ₂ + W01 ℛ₁ℛ₂′ + W10 ℛ₁′ℛ₂ + W11 ℛ₁′ℛ₂′`).
  * `channel=:pm` (`:plusminus`): the 4 mixed coefficients (mirror carry on parent2).
  * `channel=:mm` (`:minusminus`): all-zero (structural).
`parent = (l=,m=,n=,ω=)`.  `a=0` uses the analytic-3-j projection (`method=:analytic`), else the GL
route over spheroidals.  {T}/BigFloat-inherited.  Δ = r²−2r+a².
"""
function qqnm_source_coeffs(parent1, parent2, channel::Symbol, l_Q::Integer;
                            a=0.0, r, T=ComplexF64, ny::Integer=64, method::Symbol=:auto)
    C = T <: Complex ? T : Complex{T}
    RT = real(C)
    a == 0 || method == :gl || (method === :auto) ||
        throw(ArgumentError("a≠0 requires method=:gl"))
    rs = r isa AbstractVector ? r : [r]
    m_Q = parent1.m + parent2.m
    identical = (parent1 === parent2) ||
                (parent1.l==parent2.l && parent1.m==parent2.m && parent1.n==parent2.n &&
                 parent1.ω==parent2.ω)

    if channel in (:mm, :minusminus)
        z = zeros(C, length(rs))
        return (W00=z, W01=z, W10=copy(z), W11=copy(z), r=collect(rs), zero=qqnm_mm_is_zero())
    end

    # General-a ++ (a≠0): delegate to the Q4-H Deliverable-A Γ(θ)-coupled GL joint projection
    # (qqnm_angular.jl). Reduces to the a=0 analytic-3-j branch below as a→0 (gate G-Q4H-3).
    if a != 0 && channel in (:pp, :plusplus)
        return qqnm_source_coeffs_gl(parent1, parent2, l_Q; a=a, r=rs, T=C, ny=ny)
    end

    # projection coefficients (a=0 → analytic 3-j for ν=1,2; A₃,A₄ drop).  General-a is [FUTURE-gl].
    use_analytic = (a == 0)
    if use_analytic
        (P_A, P_B) = _qqnm_proj_a0_analytic(RT, parent1.l, parent1.m, l_Q, m_Q)
        cA1 = C((parent1.l+2)*(parent1.l-1)) * C(P_B)         # A₁ ↦ [L₂†S]² ↦ ₋₁Y² ↦ 3jB
        cA2 = C(sqrt(RT((parent1.l+2)*(parent1.l+1)*parent1.l*(parent1.l-1)))) * C(P_A)  # A₂ ↦ 3jA
    else
        throw(ArgumentError("general-a projection is method=:gl (Deliverable-B [PRELIMINARY]); " *
                            "call qqnm_source_coeffs_gl for the spheroidal GL route"))
    end

    W00 = zeros(C, length(rs)); W01 = zeros(C, length(rs))
    W10 = zeros(C, length(rs)); W11 = zeros(C, length(rs))
    # +− (a=0): the per-atom child projection g_{ij}(r) (Q4-H-2), computed once for all r.  This is the
    # DECISIVE angular reduction that cancels the (r−2)⁻² pole (owner CSV ⇒ finite +−); it REPLACES the
    # Q3 placeholder (16 atoms weight 1).  Both radial chains are the D₂† chain (per the 𝔠_- horizon
    # analysis: at a=0 the mirror i↦−i lives in the b_ij explicit `im`, NOT the radial ladder; the −
    # radial ℛ_-=₊₂R_{lmω₂} shares the D₂† chain, ω₂=parent2.ω — the child sits at ω_Q=ω₁+ω₂).
    Gpm = channel in (:pm, :plusminus) ? _qqnm_pm_child_gij(parent1, parent2, l_Q, rs, C; ny=max(200, ny)) : nothing
    for (k, rr) in enumerate(rs)
        Γ = C(rr + im*a*0); Γb = C(rr)                        # a=0: Γ=Γb=r
        Δ = C(rr^2 - 2rr + a^2); pref = -(Δ^6)/16
        if channel in (:pp, :plusplus)
            (α1, β1) = _qqnm_parent_chain(parent1, a, rr, C)
            (α2, β2) = identical ? (α1, β1) : _qqnm_parent_chain(parent2, a, rr, C)
            c00=zero(C); c01=zero(C); c10=zero(C); c11=zero(C)
            for (ν, cν) in ((1, cA1), (2, cA2))               # a=0: only ν=1,2 contribute
                (C00, C01, C10, C11) = _qqnm_Cab_pp(ν, α1, β1, α2, β2, Γ, Γb, rr, C)
                c00 += cν*C00; c01 += cν*C01; c10 += cν*C10; c11 += cν*C11
            end
            W00[k]=pref*c00; W01[k]=pref*c01; W10[k]=pref*c10; W11[k]=pref*c11
        elseif channel in (:pm, :plusminus)
            # ORDERED bilinear: X_i^+ (parent1 D₂† chain) × X_j^- (parent2 D₂† chain), weighted by the
            # per-atom child projection g_{ij}(r) (`_qqnm_pm_child_gij`).  Support (i,j)∈QQNM_PM_SUPPORT.
            (αp, βp) = _qqnm_parent_chain(parent1, a, rr, C)
            (αm, βm) = identical ? (αp, βp) : _qqnm_parent_chain(parent2, a, rr, C)
            c00=zero(C); c01=zero(C); c10=zero(C); c11=zero(C)
            for (pr, gv) in Gpm
                i, j = pr; g = gv[k]
                c00 += g*αp[i+1]*αm[j+1]; c01 += g*αp[i+1]*βm[j+1]
                c10 += g*βp[i+1]*αm[j+1]; c11 += g*βp[i+1]*βm[j+1]
            end
            W00[k]=pref*c00; W01[k]=pref*c01; W10[k]=pref*c10; W11[k]=pref*c11
        else
            throw(ArgumentError("unknown channel $channel"))
        end
    end
    if channel in (:pp, :plusplus) && identical
        return (W00=W00, W01=(W01 .+ W10), W11=W11, r=collect(rs), m_Q=m_Q)
    end
    (W00=W00, W01=W01, W10=W10, W11=W11, r=collect(rs), m_Q=m_Q)
end

"""
    qqnm_source(parent1, parent2, channel, l_Q; a, r, R, dR, T=ComplexF64) -> Vector

The projected second-order source `Q(r)` on the collocation grid `r`, given the parent Hertz radial
`ℛ = ₊₂R_{lmω}` and its derivative `ℛ′` (from the Q2-H s=+2 Hertz path, `ℛ(scri)=1` pin) sampled at
`r`.  For identical parents `Q = W00 ℛ² + W01 ℛℛ′ + W11 ℛ′²`; for a distinct pair `Q = W00 ℛ₁ℛ₂ +
W01 ℛ₁ℛ₂′ + W10 ℛ₁′ℛ₂ + W11 ℛ₁′ℛ₂′` (pass `R`,`dR` as 2-tuples of vectors).  If `R` is omitted only
the coefficient NamedTuple is returned (the R-independent certified object, gate G-Q3-1).
{T}/BigFloat-inherited.
"""
function qqnm_source(parent1, parent2, channel::Symbol, l_Q::Integer;
                     a=0.0, r, R=nothing, dR=nothing, T=ComplexF64, ny::Integer=64)
    co = qqnm_source_coeffs(parent1, parent2, channel, l_Q; a=a, r=r, T=T, ny=ny)
    R === nothing && return co
    if haskey(co, :W10)                 # distinct-parent pair
        (R1, R2) = R; (dR1, dR2) = dR
        return co.W00 .* (R1 .* R2) .+ co.W01 .* (R1 .* dR2) .+
               co.W10 .* (dR1 .* R2) .+ co.W11 .* (dR1 .* dR2)
    else                                # identical parents
        return co.W00 .* (R .^ 2) .+ co.W01 .* (R .* dR) .+ co.W11 .* (dR .^ 2)
    end
end

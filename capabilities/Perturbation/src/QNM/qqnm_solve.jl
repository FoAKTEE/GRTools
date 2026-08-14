# QNM/qqnm_solve.jl — the sourced hyperboloidal solve: L(ω_Q)·ψ° = 𝒮_hyp on the collocation
# grid (s=−2) at the fixed quadratic frequency; the scri readout ψ°(R=0) = E = [rΨ₄⁽²⁾]^{(∞)}
# (2410 supp:187; Paper I:987).
#
# The source transform is DERIVED, not fitted: L_hphc = −½·F·𝒪_BL·F⁻¹ (operator conjugation;
# 𝒪_BL per Paper I:722) with F = e^{−iω_Q r}·r^{−1−4iω_Q}·(r−2)^{−2+2iω_Q} at a=0, so
# 𝒮_hyp = −½·F·Q.  With the Hertz prefactors (Paper I:335,407) the universal factor
# F·E_pref₁·E_pref₂ is exponential-free: r⁻⁷(r−2)⁻⁶ at a=0; at general a the branches cancel
# because the child superradiant factor is the sum of the parents' (σ_Q = σ₁+σ₂), leaving the
# INTEGER horizon exponent e₊ = −6 for all a — the RHS is horizon-analytic.
#
# Two load-bearing constraints:
#  • the horizon node r=r₊ is 0·∞ in the closed form but analytic-finite — it MUST be filled
#    with the δ→0 extrapolated limit (zeroing it makes the solve diverge with n);
#  • in the +− channel BOTH parent radials are +ω Hertz series (the mirror bookkeeping lives in
#    the b_ij phases, not the radial): any −ω̄/−ω/conj variant re-introduces a horizon pole.
# R_+− carries an unresolved overall normalization factor vs the reference tabs (see docs).
# Radial tier {T}/BigFloat; F64 fails the solve at n≥48 (measured cond ~ n⁴) — use qqnm_hphc_cond.

# ═══════════════════════════════ channel ω_Q map ══════════════════════════════════════════════════
"""
    qqnm_channel_omega(ω1, ω2, channel) -> ω_Q

The quadratic-daughter frequency (2410 main.tex:238-260 eq:AQ_Al).  In the daughter-fixed
(nonprecessing-tab) perspective all four channels feed the SAME QQNM `(ω₁+ω₂, m₁+m₂)` and
differ only in their source: the `+−` reflected mode enters through `−ω̄_{ℓ₂,−m₂}`, which is
numerically `ω₂` (main.tex:255).  `−−` is the structural zero.
"""
function qqnm_channel_omega(ω1::Number, ω2::Number, channel::Symbol)
    channel in (:pp, :plusplus, :pm, :plusminus, :mp, :minusplus, :mm, :minusminus) ||
        throw(ArgumentError("qqnm_channel_omega: unknown channel $channel"))
    return ω1 + ω2
end

# ═══════════════════════════════ Leaver Hertz radial (s=+2) + Σaₙ ══════════════════════════════════
"""
    qqnm_hertz_leaver(l, m, a, ω; N=200, T=Complex{Float64}) -> (aₙ, Σaₙ)

The Leaver minimal-solution coefficients `aₙ` (n=0..N, a₀=1) of the s=+2 Hertz radial
`₊₂R_{lmω} = X_L(r-r₋)⁻¹Δ⁻²Σaₙρⁿ` (eq:leaver_hertz) and their sum `Σaₙ` (the linear-Ψ₄ scri
normalization `[rΨ₄⁽¹⁾]^{(∞)}=(ω⁴/2)Σaₙ`, Paper I:542).  The minimal (recessive) solution ratio is
`ρₙ=aₙ/a_{n−1}=−γₙ/(βₙ+αₙρ_{n+1})` by backward recurrence (seed ρ_{N+1}=0) over the Leaver
recurrence (`_qnm_leaver_c`, s=+2).  {T}/BigFloat-inherited (λ at the request precision).
"""
function qqnm_hertz_leaver(l::Integer, m::Integer, a::Real, ω::Number;
                           N::Integer=200, T::Type=Complex{Float64})
    RT = real(T)
    ωT = T(ω); aT = RT(a)
    λ = T(_qnm_angular_lambda(2, l, m, aT, ωT))          # s=+2 spheroidal eigenvalue at γ=aω
    c0, c1, c2, c3, c4 = _qnm_leaver_c(2, m, aT, ωT, λ)
    αn(n) = n^2 + (c0 + 1) * n + c0
    βn(n) = -2 * n^2 + (2 + c1) * n + c3
    γn(n) = n^2 + (c2 - 3) * n + c4 - c2 + 2
    ρ = zeros(T, N + 2)                                   # ρ[n+1] = aₙ/a_{n−1}; ρ[N+2]=0 seed
    for n in N:-1:1
        ρ[n+1] = -γn(n) / (βn(n) + αn(n) * ρ[n+2])
    end
    aₙ = Vector{T}(undef, N + 1); aₙ[1] = one(T)
    for n in 1:N
        aₙ[n+1] = aₙ[n] * ρ[n+1]
    end
    return aₙ, sum(aₙ)
end

# The Leaver series P(ρ)=Σaₙρⁿ and P′(ρ)=Σ n aₙρⁿ⁻¹ at ρ=(r−2)/r (a=0).
@inline function _qqnm_leaver_series(aₙ::AbstractVector{T}, ρ::Number) where {T}
    P = zero(T); Pp = zero(T); ρk = one(T)
    @inbounds for k in 0:(length(aₙ) - 1)
        P += aₙ[k+1] * ρk
        k >= 1 && (Pp += k * aₙ[k+1] * ρk / ρ)
        ρk *= ρ
    end
    return P, Pp
end

# The Leaver-gauge log-derivative ġ_j(r)=E_pref_j′/E_pref_j (a=0), E_pref=X_L(ω)(r−r₋)⁻¹Δ⁻² (Paper I:407).
@inline _qqnm_gdot(r::Number, ω::Number) = im*ω - (2im*ω + 2)/(r - 2) + (4im*ω - 3)/r

# General-a (Kerr) ġ_j = E_pref_j′/E_pref_j, E_pref_j = X_L(ω_j,m_j)(r−r₋)⁻¹Δ⁻² (Paper I:335,407),
# X_L = e^{iω_j r}(r−r₊)^{−iσ₊}(r−r₋)^{2iω_j+iσ₊}, σ₊=(ω_j r₊ − a m_j/2)/√(1−a²) (== _qnm_leaver_c σp):
#   ġ = iω_j − (iσ₊+2)/(r−r₊) + (2iω_j+iσ₊−3)/(r−r₋).  a=0 (σ₊=2ω,r₊=2,r₋=0) ⇒ the a=0 form above.
@inline _qqnm_gdot(r::Number, ω::Number, σp::Number, rp::Number, rm::Number) =
    im*ω - (im*σp + 2)/(r - rp) + (2im*ω + im*σp - 3)/(r - rm)

# 𝒮_hyp(r) = −½ F Q at a=0, on the collocation rvals.  Scri node (r=∞) → 0 (𝒮 ~ R¹³);
# horizon node → the analytic limit (file-header constraint).  Ṗ_j = ġ_j P_j + P_j′ρ′ with
# ġ_j = E_pref_j′/E_pref_j.
function _qqnm_rhs_hyp(parent1, parent2, channel::Symbol, l_Q::Integer,
                       rvals::AbstractVector, T::Type; a::Real=0,
                       N_leaver::Integer=400, hertz_scale::Number=1)
    C = T <: Complex ? T : Complex{T}
    # a≠0 routes to the general-a transform; the a=0 body below is the regression-gated reference path.
    iszero(a) || return _qqnm_rhs_hyp_kerr(parent1, parent2, channel, l_Q, rvals, C, real(C)(a);
                                           N_leaver=N_leaver, hertz_scale=hertz_scale)
    ω1 = C(parent1.ω)
    a1, _ = qqnm_hertz_leaver(parent1.l, parent1.m, 0, ω1; N=N_leaver, T=C)
    a1 = C(hertz_scale) .* a1
    if channel in (:pm, :plusminus)                       # − factor at ω₂=parent2.ω (+ω, NOT −ω̄ mirror)
        ω2 = C(parent2.ω)
        if ω2 == ω1
            a2 = a1
        else
            a2, _ = qqnm_hertz_leaver(parent2.l, parent2.m, 0, ω2; N=N_leaver, T=C)
            a2 = C(hertz_scale) .* a2
        end
    else
        ω2 = ω1; a2 = a1
    end
    # per-r builder over a finite vector r>r₊ (r₊=2 at a=0)
    function build(rC::AbstractVector)
        co = qqnm_source_coeffs(parent1, parent2, channel, l_Q; a=0, r=rC, T=C)
        W00 = co.W00; W01 = co.W01; W11 = co.W11
        W10 = haskey(co, :W10) ? co.W10 : zero.(W01)      # ++ identical: W10 already folded into W01
        out = zeros(C, length(rC))
        @inbounds for k in eachindex(rC)
            r = C(rC[k]); ρ = (r - 2)/r; ρp = C(2)/r^2
            P1, dP1 = _qqnm_leaver_series(a1, ρ); Ṗ1 = _qqnm_gdot(r, ω1)*P1 + dP1*ρp
            P2, dP2 = _qqnm_leaver_series(a2, ρ); Ṗ2 = _qqnm_gdot(r, ω2)*P2 + dP2*ρp
            pref = -C(1)/2 * r^(-7) * (r - 2)^(-6)
            out[k] = pref*(W00[k]*P1*P2 + W01[k]*P1*Ṗ2 + W10[k]*Ṗ1*P2 + W11[k]*Ṗ1*Ṗ2)
        end
        out
    end
    RT = real(C)
    𝒮 = zeros(C, length(rvals))
    interior = Int[k for k in eachindex(rvals) if isfinite(rvals[k]) && rvals[k] > 2]  # r>r₊
    horizon  = Int[k for k in eachindex(rvals) if isfinite(rvals[k]) && rvals[k] <= 2] # r==r₊
    if !isempty(interior)
        Sv = build(RT[RT(rvals[k]) for k in interior])
        for (j, k) in enumerate(interior); 𝒮[k] = Sv[j]; end
    end
    # horizon node r=r₊: 𝒮 is 0·∞ in closed form ⇒ analytic limit by 3-pt quadratic extrapolation δ→0
    for k in horizon
        rp = RT(rvals[k])                                 # r₊ (=2 at a=0)
        δs = (RT(1)/100, RT(1)/1000, RT(1)/10000)
        Sδ = C[build(RT[rp + δ])[1] for δ in δs]
        𝒮[k] = _qqnm_extrap0(collect(δs), Sδ)             # Lagrange extrapolation to δ=0
    end
    return 𝒮                                              # scri nodes (r=∞) stay 0 (𝒮 ~ R¹³ → 0)
end

# 𝒮_hyp at general a: the conjugation F acquires an r=0 branch r^{1−4iω_Q} and superradiant
# exponents p₊,p₋ at r₊,r₋ (in code below); σ_Q = σ₁+σ₂ cancels the branches so
# U = F·E_pref₁·E_pref₂ = r^{e₀}(r−r₊)^{e₊}(r−r₋)^{e₋} is exponential-free with e₊ = −6 exactly
# (Hertz E_pref exponents per Paper I:335,407).  Same horizon fill as a=0.  ++ only; +− [FUTURE].
function _qqnm_rhs_hyp_kerr(parent1, parent2, channel::Symbol, l_Q::Integer,
                           rvals::AbstractVector, C::Type, aR::Real;
                           N_leaver::Integer=400, hertz_scale::Number=1)
    RT = real(C)
    b  = sqrt(one(RT) - aR^2)                              # √(1−a²)  (M=1)
    rp = one(RT) + b                                        # r₊ = 1 + √(1−a²)
    rm = one(RT) - b                                        # r₋ = 1 − √(1−a²)
    m1 = parent1.m; ω1 = C(parent1.ω)
    σ1 = (ω1*rp - aR*m1/2)/b                                # σ₊^{(1)} (Paper I:335; == _qnm_leaver_c σp)
    a1, _ = qqnm_hertz_leaver(parent1.l, m1, aR, ω1; N=N_leaver, T=C)
    a1 = C(hertz_scale) .* a1
    if channel in (:pm, :plusminus)
        m2 = parent2.m; ω2 = C(parent2.ω); σ2 = (ω2*rp - aR*m2/2)/b
        if ω2 == ω1 && m2 == m1
            a2 = a1
        else
            a2, _ = qqnm_hertz_leaver(parent2.l, m2, aR, ω2; N=N_leaver, T=C); a2 = C(hertz_scale) .* a2
        end
    else
        m2 = m1; ω2 = ω1; σ2 = σ1; a2 = a1
    end
    ωQ = ω1 + ω2; mQ = m1 + parent2.m; σQ = σ1 + σ2
    # F exponents (operator-conjugation CAS) at the child (ω_Q,m_Q):
    Fzero  = 1 - 4im*ωQ                                     # r=0 residue (a≠0 feature)
    Fplus  = (im*aR^2*ωQ - 2aR^2 + im*aR*mQ*b/2 + im*aR*mQ/2 - 2im*ωQ*b - 2im*ωQ + 2b + 2)/(aR^2 - b - 1)
    Fminus = (im*aR^2*ωQ - 2aR^2 - im*aR*mQ*b/2 + im*aR*mQ/2 + 2im*ωQ*b - 2im*ωQ - 2b + 2)/(aR^2 + b - 1)
    # U = F·E_pref₁·E_pref₂ = r^{e₀}(r−r₊)^{e₊}(r−r₋)^{e₋} (exp-free; the child σ_Q = σ₁+σ₂ kills the branch)
    e0 = Fzero
    ep = Fplus  + (-2 - im*σ1) + (-2 - im*σ2)
    em = Fminus + (-3 + 2im*ω1 + im*σ1) + (-3 + 2im*ω2 + im*σ2)
    function build(rC::AbstractVector)
        co = qqnm_source_coeffs(parent1, parent2, channel, l_Q; a=aR, r=rC, T=C)
        W00 = co.W00; W01 = co.W01; W11 = co.W11
        W10 = haskey(co, :W10) ? co.W10 : zero.(W01)
        out = zeros(C, length(rC))
        @inbounds for k in eachindex(rC)
            r = C(rC[k]); ρ = (r - rp)/(r - rm); ρp = C(rp - rm)/(r - rm)^2
            P1, dP1 = _qqnm_leaver_series(a1, ρ); Ṗ1 = _qqnm_gdot(r, ω1, σ1, rp, rm)*P1 + dP1*ρp
            P2, dP2 = _qqnm_leaver_series(a2, ρ); Ṗ2 = _qqnm_gdot(r, ω2, σ2, rp, rm)*P2 + dP2*ρp
            pref = -C(1)/2 * r^e0 * (r - rp)^ep * (r - rm)^em
            out[k] = pref*(W00[k]*P1*P2 + W01[k]*P1*Ṗ2 + W10[k]*Ṗ1*P2 + W11[k]*Ṗ1*Ṗ2)
        end
        out
    end
    thresh = rp*(1 + 4096*eps(RT))                          # horizon-node band (>rp roundtrip 1/Rp safe)
    𝒮 = zeros(C, length(rvals))
    interior = Int[k for k in eachindex(rvals) if isfinite(rvals[k]) && rvals[k] > thresh]
    horizon  = Int[k for k in eachindex(rvals) if isfinite(rvals[k]) && rvals[k] <= thresh]
    if !isempty(interior)
        Sv = build(RT[RT(rvals[k]) for k in interior])
        for (j, k) in enumerate(interior); 𝒮[k] = Sv[j]; end
    end
    for k in horizon                                        # analytic-limit fill (Richardson δ→0)
        rph = RT(rvals[k])
        δs = (RT(1)/100, RT(1)/1000, RT(1)/10000)
        Sδ = C[build(RT[rph + δ])[1] for δ in δs]
        𝒮[k] = _qqnm_extrap0(collect(δs), Sδ)
    end
    return 𝒮
end

# 3-pt Lagrange extrapolation of y(x) to x=0 (the analytic horizon limit; y analytic in δ=r−r₊).
@inline function _qqnm_extrap0(xs::AbstractVector, ys::AbstractVector{C}) where {C}
    n = length(xs); L0 = zero(C)
    @inbounds for i in 1:n
        li = one(real(C))
        for j in 1:n
            i == j && continue
            li *= (0 - xs[j]) / (xs[i] - xs[j])
        end
        L0 += C(li) * ys[i]
    end
    return L0
end

# ═══════════════════════════════ conditioning of L(ω_Q) ════════════════════════════════════════════
"""
    qqnm_hphc_cond(l_Q, m_Q, a, ω_Q; n=64, T=Float64) -> (cond, absdet_over_cond)

The 2-norm condition number `cond₂(L(ω_Q))` of the s=−2 hyperboloidal operator at the FIXED quadratic
frequency (the sourced-solve matrix, NOT the eigen QEP).  Also returns `|det L|^{1/N}/cond` as the
resonance-proximity sanity metric (small ⇒ near a child QNM ⇒ resonance).  `cond` via SVD at Float64
(the eigen-study's companion metric does not cover this fixed-ω solve — measured here directly).
"""
function qqnm_hphc_cond(l_Q::Integer, m_Q::Integer, a::Real, ω_Q::Number;
                        n::Integer=64, T::Type=Float64)
    λ_Q = _qnm_angular_lambda(-2, l_Q, m_Q, T(a), ω_Q)
    op = hphc_radial_operator(-2, m_Q, T(a), ω_Q, λ_Q; n=n, T=T)
    L = op.L
    svals = svdvals(Matrix{ComplexF64}(L))
    condL = svals[1] / svals[end]
    N = size(L, 1)
    logdet_mag = sum(log, svals) / N                      # (Σ log σ)/N = log |det|^{1/N}
    return (cond=condL, det_metric=exp(logdet_mag) / condL, smin=svals[end], smax=svals[1])
end

# ═══════════════════════════════ the sourced solve ════════════════════════════════════════════════
"""
    qqnm_solve(parent1, parent2, channel, l_Q; a=0, n=64, T=Complex{BigFloat}, N_leaver=400) -> NamedTuple

Solve `L(ω_Q)·ψ° = 𝒮_hyp` (s=−2) on the collocation grid and read the scri amplitude
`E = ψ°(R=0)`, via the derived transform `𝒮_hyp = −½ F(r) Q(r)` (file header) with the
horizon-node analytic-limit fill.  `parent = (l=,m=,n=,ω=)`; `channel ∈ {:pp,:pm,:mp,:mm}`.
`−−` returns `E=0` structurally.  `R_+−` carries an unresolved overall normalization factor
vs the reference tabs (see docs).  `++` supports general `a`; `+−` at `a≠0` is [FUTURE].
Returns `(; E, omega_Q, lambda_Q, psi, R, rvals, cond, Sigma_an, resonance_ok, n)`.
{T}/BigFloat-inherited (radial tier).
"""
function qqnm_solve(parent1, parent2, channel::Symbol, l_Q::Integer;
                    a::Real=0, n::Integer=64, T::Type=Complex{BigFloat}, N_leaver::Integer=400,
                    hertz_scale::Number=1)
    C = T <: Complex ? T : Complex{T}
    RT = real(C)
    m_Q = parent1.m + parent2.m
    ω_Q = qqnm_channel_omega(C(parent1.ω), C(parent2.ω), channel)
    # −− structural zero: no solve.
    if channel in (:mm, :minusminus)
        return (; E=zero(C), omega_Q=ω_Q, lambda_Q=zero(C), psi=C[], R=RT[], rvals=RT[],
                cond=RT(NaN), Sigma_an=zero(C), resonance_ok=true, n=n, zero=qqnm_mm_is_zero())
    end
    # ++ supports general a; the +− Kerr source is [FUTURE].
    (a == 0 || channel in (:pp, :plusplus)) ||
        throw(ArgumentError("qqnm_solve: Kerr (a≠0) source injection is [FUTURE] for channel $channel " *
                            "(only ++ is derived at general a — Q5-STAGE-1)"))
    λ_Q = C(_qnm_angular_lambda(-2, l_Q, m_Q, RT(a), ω_Q))
    op = hphc_radial_operator(-2, m_Q, RT(a), ω_Q, λ_Q; n=n, T=RT)
    R = op.R; rvals = [iszero(Rj) ? RT(Inf) : one(RT) / Rj for Rj in R]
    # Σaₙ (parent-1 Hertz Leaver); `hertz_scale` must flow into BOTH the RHS and Σaₙ so the
    # bilinear-over-product R stays invariant under ℛ→cℛ.
    _, Σaₙ = qqnm_hertz_leaver(parent1.l, parent1.m, a, C(parent1.ω); N=N_leaver, T=C)
    Σaₙ = C(hertz_scale) * Σaₙ
    𝒮 = _qqnm_rhs_hyp(parent1, parent2, channel, l_Q, rvals, C; a=a, N_leaver=N_leaver, hertz_scale=hertz_scale)
    ψ = op.L \ 𝒮
    E = ψ[op.iscri]
    # conditioning + resonance sanity (Float64 SVD)
    cm = qqnm_hphc_cond(l_Q, m_Q, a, ComplexF64(ω_Q); n=n, T=Float64)
    resonance_ok = cm.det_metric > 1e-8                  # not near a child QNM (|det|^{1/N}/cond)
    return (; E=E, omega_Q=ω_Q, lambda_Q=λ_Q, psi=ψ, R=R, rvals=rvals, cond=RT(cm.cond),
            Sigma_an=Σaₙ, resonance_ok=resonance_ok, n=n)
end

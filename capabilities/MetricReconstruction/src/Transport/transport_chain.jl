# transport_chain.jl — the seven generalized transport equations B1…B7 (LY:245–318) + the S_μ
# auxiliary, specialized to the static Schwarzschild shell, and the Fig.-1 DAG driver.
#
#   B1 lambda_transport   Δ_{[1,1,3,-1]} λ      = −Ψ4                              (LY:246)
#   B2 hmm_transport      Δ_{[-1,1,2,-2]} h_m̄m̄ = 2λ                               (LY:250)
#   B3 psi3_transport     Δ_{[4,0,2,0]} Ψ3       = δ_{[0,4,0,-1]}Ψ4 + R_h          (LY:257)
#   B4 pi_transport       Δ_{[0,0,1,-1]} π       = −Ψ3 − Φ21 − …                    (LY:262)
#   B5 hlm_transport      Δ_{[0,1,0,-2]} h_lm̄   = 2π − τ̄h_mm̄ − τh_m̄m̄            (LY:276)
#   B6 psi2_transport     Δ_{[3,0,0,0]} Ψ2       = −3Ψ2μ^{(1)} + δ_{[0,2,0,-2]}Ψ3 + R_g  (LY:284)
#   B7 hll_transport      Δ_{[-1,1,-2,-2]} h_ll = … + (2/μ)S_μ                     (LY:311)
#
# For the shell (LY:771) the spin-2 (λ, h_m̄m̄) and spin-1 (Ψ3, π, h_lm̄) sectors VANISH; only the
# spin-0 Ψ2 (B6) and h_ll (B7) survive. This file provides the vanishing check, the Ψ2 solve, the
# S_μ auxiliary + its δ-cancellation, and the h_ll integrating-factor solve.

# ================================================================================================
# The Δ-operator quartets of the chain (LY labels) — the equation-level structural surface.
# ================================================================================================
"The LHS Δ_{[a,b,c,d]} weight quartets of the transport chain, keyed by LY equation label."
const CHAIN_DELTA_OPS = (
    mu_transport_final = Δw(2, 0, 1, 1),     # LY:230
    trace_transport    = Δw(-1, 1, 0, 0),    # LY:235
    lambda_transport   = Δw(1, 1, 3, -1),    # LY:246  (B1)
    hmm_transport      = Δw(-1, 1, 2, -2),   # LY:250  (B2)
    psi3_transport     = Δw(4, 0, 2, 0),     # LY:257  (B3)
    pi_transport       = Δw(0, 0, 1, -1),    # LY:262  (B4)
    hlm_transport      = Δw(0, 1, 0, -2),    # LY:276  (B5)
    psi2_transport     = Δw(3, 0, 0, 0),     # LY:284  (B6)
    hll_transport      = Δw(-1, 1, -2, -2),  # LY:311  (B7)
)

# ================================================================================================
# ANCHOR-5 vanishing (LY:771): the spin-2 / spin-1 sectors collapse to zero for the shell.
# ================================================================================================
"""
    vanishing_shell(T, M, R; Ψ4, Φ21, Rh) -> NamedTuple

The spin-2 / spin-1 collapse of LY:771. Given the shell inputs `Ψ4=0` (no stationary Teukolsky
source, homogeneous spherical), `Φ21=0`, and `Rh=0` (monopole, [`Rh_shell`](@ref)), every source of
B1…B5 vanishes and spherical symmetry forces the homogeneous solution to zero:

    Ψ4=0            ⇒ (B1) λ=0        ⇒ (B2) h_m̄m̄=0
    δΨ4=0, Rh=0     ⇒ (B3) Ψ3=0
    Ψ3=Φ21=0, π_bg=τ_bg=0 ⇒ (B4) π=0 ⇒ (B5) h_lm̄=0

Returns the five fields as zero [`ShellProfile`](@ref)s together with the boolean `all_zero` (the
G-E2 vanishing certificate). The check is structural: it verifies each B-equation SOURCE is zero
given the inputs, so the spin-s homogeneous+spherical solution is zero.
"""
function vanishing_shell(::Type{T}, M, R; Ψ4, Φ21, Rh) where {T}
    z = ShellProfile{T}(lzero(T), lzero(T), zero(T), zero(T), T(R))
    # B1 source = −Ψ4 = 0 ; B3 source = δΨ4 + Rh = 0 ; B4 source = −Ψ3−Φ21−… = 0 (Ψ3=Φ21=0, τ_bg=π_bg=0)
    src_lambda_zero = _iszeroprofile(Ψ4)
    src_psi3_zero   = _iszeroprofile(Rh)               # δΨ4=0 (Ψ4=0), Rh=0
    src_pi_zero     = _iszeroprofile(Φ21)              # Ψ3=0, Φ21=0, τ_bg=π_bg=0 ⇒ src=0
    all_zero = src_lambda_zero && src_psi3_zero && src_pi_zero
    (λ = z, h_m̄m̄ = z, Ψ3 = z, π = z, h_lm̄ = z, all_zero = all_zero)
end

_iszeroprofile(x::ShellProfile) = isempty(x.inside.c) && isempty(x.outside.c) && x.dc == 0 && x.dpc == 0

# ================================================================================================
# B6 Ψ2 (LY:283–287; shell eq:psi2_eq_shell LY:805) — spin-0 mass-shift completion Weyl scalar.
# ================================================================================================
"""
    solve_psi2_shell(T, M, R, h_mm̄, μ, Rg) -> ShellProfile

Solve the B6 transport for `Ψ2^{(1)}` on the static Schwarzschild shell. From B6 with `Ψ3=0`,
`Δ_{[3,0,0,0]}=−(f/2)∂_r+3μ_bg`, one has (eq:psi2_eq_shell, LY:805)

    ∂_r(r³Ψ2^{(1)}) = −(2r³/f) R_g − (6M/f) μ^{(1)},

where the μ-term is f-free via the trace relation `−(6M/f)μ^{(1)} = (3M/2)∂_r h_mm̄^{(1)}` (μ = −f/4
∂_r h_mm̄) and the `R_g` term is scaled by the smooth `−2r³/f` (value+derivative at R only, since
`R_g` is pure δ′/δ). Integration constants are the **[AXIOM-E]** mass-shift boundary conditions
`Ψ2,reg⁺=−M/r³` and `[Ψ2,reg]=−M/R³` (LY:833). Reproduces eq:psi2_shell (LY:820–828).
"""
function solve_psi2_shell(::Type{T}, M, R, h_mm̄::ShellProfile{T}, μ::ShellProfile{T},
                          Rg::ShellProfile{T}) where {T}
    Mt, Rt = T(M), T(R)
    fR = 1 - 2Mt / Rt
    fpR = 2Mt / Rt^2
    # μ-term:  −(6M/f)μ = (3M/2)∂_r h_mm̄   (Θ part only; f cancels — genuine, no 1/f)
    src_mu = pscale(3Mt / 2, pderiv(h_mm̄))
    # R_g-term:  −(2r³/f)R_g  (R_g pure δ′/δ ⇒ use value/deriv of −2r³/f at R)
    smval = -2Rt^3 / fR                                   # (−2r³/f)|_R
    smder = -2 * (3Rt^2 * fR - Rt^3 * fpR) / fR^2          # d/dr(−2r³/f)|_R
    src_Rg = psmooth_vd(smval, smder, Rg)
    src = src_mu + src_Rg                                  # = ∂_r(r³Ψ2)
    # integrate  r³Ψ2 = ∫src, then Ψ2 = (∫src)/r³, with the mass-shift BCs.
    Y_in = lintegrate(src.inside)                          # ∫ interior Θ part
    # outside: no interior source integrated here; BC r³Ψ2⁺ = −M  (⇒ Ψ2,reg⁺ = −M/r³)
    numer_out = laurent(T, 0 => -Mt)
    # inside constant from the jump BC [Ψ2,reg] = −M/R³  ⇒  Y_in(R)+c_in = 0
    c_in = -leval(Y_in, Rt)
    numer_in = Y_in + laurent(T, 0 => c_in)
    inv_r3 = laurent(T, -3 => one(T))
    Ψ2_in = numer_in * inv_r3
    Ψ2_out = numer_out * inv_r3
    # δ of Ψ2:  the δ′ in r³Ψ2 (∫δ′=δ, coeff src.dpc) becomes (src.dpc/R³)δ  [ (1/r³)δ = (1/R³)δ ]
    Ψ2_dc = src.dpc / Rt^3
    ShellProfile{T}(Ψ2_in, Ψ2_out, Ψ2_dc, zero(T), Rt)
end

# ================================================================================================
# S_μ auxiliary (LY:289–307) + the δ-cancellation (ANCHOR-11, LY:847–853).
# ================================================================================================
"""
    smu_cancel_shell(T, μ, Ψ2, Λ) -> ShellProfile

The core of the S_μ auxiliary reduced for the shell (LY:849): `∂_r μ^{(1)} − Ψ2^{(1)} − 2Λ^{(1)}`.
Feeding it the reconstructed `μ, Ψ2` and the source `Λ`, the three distributional (δ) contributions
must CANCEL EXACTLY (`.dc ≈ 0`) — the stringent ANCHOR-11 internal-consistency check — leaving the
smooth `M(R−3M)/(f(R)Rr³)Θ(R−r) + M/r³ Θ(r−R)`.
"""
function smu_cancel_shell(::Type{T}, μ::ShellProfile{T}, Ψ2::ShellProfile{T}, Λ::ShellProfile{T}) where {T}
    pderiv(μ) - Ψ2 - pscale(T(2), Λ)
end

# ================================================================================================
# B7 h_ll (LY:309–318; shell eq:hll_eq_shell LY:836) — the spin-0 radial completion.
# ================================================================================================
"""
    HllField{T}

The reconstructed `h_ll^{(1)}` on the shell, stored as interior/exterior NUMERATOR Laurents with an
implicit `1/f²` factor: `h_ll(r) = numer(r)/f(r)²` (the integrating factor of the B7 transport is
`f²`, LY:863). Evaluate with [`hll_eval`](@ref).
"""
struct HllField{T}
    numer_in::Laurent{T}
    numer_out::Laurent{T}
    R::T
    M::T
end

"Evaluate `h_ll^{(1)}(r) = numer(r)/f(r)²` (interior branch for r<R, exterior for r>R)."
function hll_eval(h::HllField{T}, r) where {T}
    fr = 1 - 2h.M / r
    numer = r < h.R ? h.numer_in : h.numer_out
    leval(numer, r) / fr^2
end

"""
    solve_hll_shell(T, M, R, μ, cancel; Cp=0) -> HllField

Solve the B7 transport for `h_ll^{(1)}` on the shell. The equation (eq:hll_eq_shell_1, LY:836) is
`(f/2 ∂_r + 2M/r²)h_ll = S`, `S = −8μ^{(1)}/f + (2/μ_bg)(∂_r μ^{(1)} − Ψ2 − 2Λ)`. Its integrating
factor is `f²` (∫4M/(r²f) = 2 ln f), so `h_ll = f^{-2}(C + ∫2fS dr)` with `2fS = −16μ^{(1)} − 8r·
cancel` a Laurent (all `1/f` cancels; `cancel` = [`smu_cancel_shell`](@ref)). The exterior amplitude
`Cp` is the **[AXIOM-E]** outside-Schwarzschild input (paper default 0), and continuity of `h_ll`
across the shell fixes the interior constant. Reproduces eq:hll_shell (LY:863–869) + continuity
(ANCHOR-13, LY:870).
"""
function solve_hll_shell(::Type{T}, M, R, μ::ShellProfile{T}, cancel::ShellProfile{T};
                         Cp = zero(T)) where {T}
    Mt, Rt = T(M), T(R)
    rL = laurent(T, 1 => one(T))                          # r
    # 2fS = −16μ − 8r·cancel   (per branch, Laurent)
    twofS_in  = lscale(-16, μ.inside)  + lscale(-8, rL * cancel.inside)
    twofS_out = lscale(-16, μ.outside) + lscale(-8, rL * cancel.outside)
    I_in  = lintegrate(twofS_in)
    I_out = lintegrate(twofS_out)
    numer_out = laurent(T, 0 => Cp) + I_out               # C⁺ + ∫2fS⁺
    # continuity  numer_in(R)/f² = numer_out(R)/f²  ⇒  numer_in(R)=numer_out(R)
    Cm = leval(numer_out, Rt) - leval(I_in, Rt)
    numer_in = laurent(T, 0 => Cm) + I_in
    HllField{T}(numer_in, numer_out, Rt, Mt)
end

# ================================================================================================
# The Fig.-1 DAG driver + the [FUTURE] smooth-source seam.
# ================================================================================================
"""
    ShellSource(; M, R)

The static Schwarzschild thin-shell source (LY Supplemental Material). Element type `T` is chosen at
`transport_solve` time.
"""
struct ShellSource
    M::Float64
    R::Float64
end
ShellSource(; M, R) = ShellSource(Float64(M), Float64(R))

"""
    transport_solve(src::ShellSource, ::Type{T}=Float64) -> NamedTuple

Run the full Fig.-1 transport DAG (topological order — trace branch BEFORE h_ll, LY:414–422) for the
static Schwarzschild shell, returning the reconstructed fields:
`(; Φ22, Λ, Rg, h_mm̄, μ, γ, Ψ2, cancel, h_ll, sources…)`, each a [`ShellProfile`](@ref) (or
[`HllField`](@ref) for h_ll). This is the object the G-E2…G-E6 gates interrogate.
"""
function transport_solve(src::ShellSource, ::Type{T} = Float64) where {T}
    M, R = T(src.M), T(src.R)
    # sources (HOLE-6): T_ab → Φ,Λ + Bianchi R_g/R_h
    ric = shell_ricci_profiles(T, M, R)
    Rg = Rg_shell(T, M, R)
    Rh = Rh_shell(T, M, R)
    zΨ4 = ShellProfile{T}(lzero(T), lzero(T), zero(T), zero(T), R)   # Ψ4^{(1)}=0 (spherical)
    zΦ21 = ShellProfile{T}(lzero(T), lzero(T), zero(T), zero(T), R)  # Φ21^{(1)}=0
    # spin-2 / spin-1 vanishing (ANCHOR-5)
    van = vanishing_shell(T, M, R; Ψ4 = zΨ4, Φ21 = zΦ21, Rh = Rh)
    # trace branch (A1/A2): h_mm̄ ← Φ22, then μ, γ
    h_mm̄ = solve_hmm_shell(T, M, R, ric.Φ22)
    μ = mu_from_trace(T, M, h_mm̄)
    γ = gamma_from_trace(T, M, h_mm̄)
    # Ψ2 (B6) ← μ, R_g ; S_μ δ-cancellation ; h_ll (B7)
    Ψ2 = solve_psi2_shell(T, M, R, h_mm̄, μ, Rg)
    cancel = smu_cancel_shell(T, μ, Ψ2, ric.Λ)
    h_ll = solve_hll_shell(T, M, R, μ, cancel)
    (Φ00 = ric.Φ00, Φ22 = ric.Φ22, Φ11 = ric.Φ11, Λ = ric.Λ, Rg = Rg, Rh = Rh,
     vanishing = van, h_mm̄ = h_mm̄, μ = μ, γ = γ, Ψ2 = Ψ2, cancel = cancel, h_ll = h_ll)
end

# `SmoothSource` + `transport_solve(::SmoothSource, T)` (the F4 DP5 smooth-source path) are defined
# in `smooth_source.jl` (included after this file), which fills the seam declared here.

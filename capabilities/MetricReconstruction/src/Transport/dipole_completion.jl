# dipole_completion.jl — the ℓ=1 angular-momentum (Lense–Thirring) completion; derivation
# F6B_DERIVATION.md, machine referee F6B_referee.md. Setting: static, vacuum, a=0, M=1, ℓ=1 m=0
# axial; Ψ₄⁽¹⁾=0 (the s=−2 ladder vanishes at ℓ=1). Closed-form chain (residual gates in
# test/transport_dipole.jl):
#   B3  Δ_{[4,0,2,0]}Ψ₃ = 0           ⇒  Ψ₃ = C₃ f/r⁴
#   B4  Δ_{[0,0,1,-1]}π = −Ψ₃         ⇒  π  = −(2C₃/3)/r³
#   B5  Δ_{[0,1,0,-2]}h_lm̄ = +2π      ⇒  h_lm̄^{ORG} = −(8C₃/3)/(r²f)
#   B6  Δ_{[3,0,0,0]}Ψ₂ = (δ+2β−2τ)Ψ₃ ⇒  Ψ₂ = 2√2 C₃ cosθ/r⁴
# CRITICAL SIGN: the B5 source is +2π (LY:276), NOT the −2π printed in RLGP-F:1215 — only +2π lets
# ONE amplitude C₃ = −3√2 i·δJ/4 reproduce both referee targets (h_tφ AND the frame-invariant Ψ₂);
# do not "fix" it back. The BL representative is recovered by the axial gauge vector ξ = c(r)sin²θ dφ.
# Background NP/tetrad SOURCED from Spacetime, never re-hardcoded. {T}-generic.

import Spacetime

# ================================================================================================
# DipoleCompletion — the static ℓ=1 m=0 axial angular-momentum completion.
# ================================================================================================
"""
    DipoleCompletion{T}

The static (ω=0) ℓ=1 m=0 AXIAL angular-momentum (Lense–Thirring) completion of the Li–Yunes
traceful-ORG transport chain on the a=0 Schwarzschild background (M=1). Constructed from the single
physical input `δJ` (= the Kerr spin `a`, the angular-momentum shift); every chain field is a
closed-form function of the derived homogeneous amplitude `C₃`.

Fields:
  - `δJ  :: T`          — **[AXIOM-F6B-4]** the angular-momentum shift; sets the metric normalization
                          `h_tφ = −2δJ sin²θ/r` (the exact O(a) Kerr Lense–Thirring term).
  - `C₃  :: Complex{T}` — the B3 homogeneous amplitude, `C₃ = −3√2 i δJ/4`, fixed by δJ (derived in
                          the file header; gate-verified by G-F6B-1 metric + G-F6B-2 Weyl).
  - `Kc  :: T`          — **[AXIOM-F6B-5]** the residual axial-gauge constant (the rigid `r²`
                          homogeneous mode of the gauge ODE); paper default `0`.

**[AXIOM] ledger** (choices not derivable from the transport ODEs):
  - **[AXIOM-F6B-1]** B3 branch: the (single, first-order) homogeneous solution `f/r⁴` — horizon-
    regular (f→0) AND decaying (∼r⁻⁴); the caller supplies no second branch.
  - **[AXIOM-F6B-2]** `K_π = 0`: π decays at infinity (no additive constant in the B4 quadrature).
  - **[AXIOM-F6B-3]** axial parity: the trace sector (`h_mm̄`, `μ⁽¹⁾`) and `h_ll` vanish — an ℓ=1
    axial perturbation carries no polar spin-0 metric content (referee ITEM-2: h_mm̄=h_ll=0).
    (Refinement: B7 DOES have a nonzero-ladder `ð h_lm̄` source; `h_ll=0` holds because the full
    source cancels in-combination, not termwise — see the file header / G-F6B-3.)
  - **[AXIOM-F6B-5]** residual gauge constant `Kc` (default 0), and the h_lm̄ decay BC `K₅=0`.
"""
struct DipoleCompletion{T}
    δJ::T
    C₃::Complex{T}
    Kc::T
end

"""
    DipoleCompletion(δJ; Kc=0) -> DipoleCompletion

Build the completion from the angular-momentum shift `δJ`. Derives `C₃ = −3√2 i δJ/4` (the amplitude
fixing `h_tφ = −2δJ sin²θ/r`; equivalently `Ψ₂⁽¹⁾ = −3i δJ cosθ/r⁴`). Element type `T = float(δJ)`
(BigFloat inherited). `Kc` is the **[AXIOM-F6B-5]** residual axial-gauge constant.
"""
function DipoleCompletion(δJ; Kc = 0)
    T = float(typeof(δJ))
    δJt = T(δJ)
    C₃ = Complex{T}(0, -3 * sqrt(T(2)) * δJt / 4)          # = −3√2 i δJ/4 = −3i δJ/(2√2)
    DipoleCompletion{T}(δJt, C₃, T(Kc))
end

# ------------------------------------------------------------------------------------------------
# f(r) helper (a=0 Schwarzschild, M=1), {T}-generic.
# ------------------------------------------------------------------------------------------------
_f(::Type{T}, r) where {T} = one(T) - 2 / T(r)

# ================================================================================================
# The chain fields as closed forms (each cited to its transport equation; all RE-VERIFIED).
# ================================================================================================
"""
    dipole_psi3_radial(dc::DipoleCompletion) -> Laurent

The radial part of `Ψ₃⁽¹⁾` (spin −1): `C₃ f/r⁴ = C₃(r⁻⁴ − 2r⁻⁵)`, the unique horizon-regular &
decaying homogeneous solution of B3 `Δ_{[4,0,2,0]}Ψ₃ = 0` (**[AXIOM-F6B-1]**). Full field is
[`dipole_psi3`](@ref) `= radial·sinθ`.
"""
dipole_psi3_radial(dc::DipoleCompletion{T}) where {T} =
    laurent(Complex{T}, -4 => dc.C₃, -5 => -2 * dc.C₃)

"`Ψ₃⁽¹⁾(r,θ) = C₃ (f/r⁴) sinθ` (spin −1; ₋₁Y₁₀ ∝ sinθ)."
dipole_psi3(dc::DipoleCompletion{T}, r, θ) where {T} =
    leval(dipole_psi3_radial(dc), r) * sin(T(θ))

"""
    dipole_pi_radial(dc::DipoleCompletion) -> Laurent

The radial part of `π⁽¹⁾` (spin −1): `−(2C₃/3) r⁻³`, from the B4 quadrature `Δ_{[0,0,1,-1]}π = −Ψ₃`
(alg coeff 0 ⇒ `π′ = 2Ψ₃/f = 2C₃/r⁴`), with `K_π=0` (**[AXIOM-F6B-2]**).
"""
dipole_pi_radial(dc::DipoleCompletion{T}) where {T} =
    laurent(Complex{T}, -3 => -2 * dc.C₃ / 3)

"`π⁽¹⁾(r,θ) = −(2C₃/3) r⁻³ sinθ` (spin −1)."
dipole_pi(dc::DipoleCompletion{T}, r, θ) where {T} =
    leval(dipole_pi_radial(dc), r) * sin(T(θ))

"""
    dipole_hlmbar_radial(dc::DipoleCompletion, r) -> value

The radial part of the ORG `h_lm̄⁽¹⁾` (spin −1): `−(8C₃/3)/(r²f)`, from the B5 integrating-factor
solve `Δ_{[0,1,0,-2]}h_lm̄ = +2π` (alg coeff `μ−2γ = −1/(2r)`; IF `(r−2)`; `K₅=0` decay). NOT a
Laurent (carries the `1/f` horizon pole). SOURCE SIGN `+2π` (transport_chain.jl/LY:276) — see header.
"""
function dipole_hlmbar_radial(dc::DipoleCompletion{T}, r) where {T}
    -(8 * dc.C₃ / 3) / (T(r)^2 * _f(T, r))
end

"`h_lm̄⁽¹⁾(r,θ) = −(8C₃/3)/(r²f)·sinθ` (ORG representative, spin −1; = 2√2 i δJ sinθ/(r²f))."
dipole_hlmbar(dc::DipoleCompletion{T}, r, θ) where {T} =
    dipole_hlmbar_radial(dc, r) * sin(T(θ))

"""
    dipole_psi2_radial(dc::DipoleCompletion) -> Laurent

The radial part of the Weyl scalar `Ψ₂⁽¹⁾` (spin 0): `2√2 C₃ r⁻⁴`, from the B6 ð-sourced solve
`Δ_{[3,0,0,0]}Ψ₂ = (δ+2β−2τ)Ψ₃` (source `√2 C₃ f cosθ/r⁵`; alg `3μ = −3/(2r)+3/r²`; the r⁻³ mass-
shift homogeneous mode is O(a²)-even ⇒ absent at O(a)). Equals the frame-invariant referee target
`−3i δJ cosθ/r⁴`.
"""
dipole_psi2_radial(dc::DipoleCompletion{T}) where {T} =
    laurent(Complex{T}, -4 => 2 * sqrt(T(2)) * dc.C₃)

"`Ψ₂⁽¹⁾(r,θ) = 2√2 C₃ r⁻⁴ cosθ = −3i δJ cosθ/r⁴` (spin 0, Y₁₀ ∝ cosθ; frame-invariant, G-F6B-2)."
dipole_psi2(dc::DipoleCompletion{T}, r, θ) where {T} =
    leval(dipole_psi2_radial(dc), r) * cos(T(θ))

# ================================================================================================
# Zero sectors (F6B_DERIVATION step-1/5; referee ITEM-2/4) — spin-2 and polar spin-0 vanish at ℓ=1.
# ================================================================================================
"`λ⁽¹⁾ = h_m̄m̄⁽¹⁾ ≡ 0` (spin-2; ℓ=1<|s|=2, s=−2 raise ladder (ℓ−s)(ℓ+s+1)=0). Returns `0`."
dipole_lambda(dc::DipoleCompletion{T}, r, θ) where {T} = zero(Complex{T})
dipole_hmbarmbar(dc::DipoleCompletion{T}, r, θ) where {T} = zero(Complex{T})
"`h_mm̄⁽¹⁾ ≡ 0` (trace/polar spin-0; **[AXIOM-F6B-3]** axial parity). Returns `0`."
dipole_hmmbar(dc::DipoleCompletion{T}, r, θ) where {T} = zero(Complex{T})
"`h_ll⁽¹⁾ ≡ 0` (polar spin-0; **[AXIOM-F6B-3]**; referee ITEM-2, gauge-invariant under axial ξ). Returns `0`."
dipole_hll(dc::DipoleCompletion{T}, r, θ) where {T} = zero(Complex{T})

# ================================================================================================
# Traceful-ORG tetrad assembly  h_ab → h_μν   (SOURCED from Spacetime.kinnersley_tetrad/kerr_metric)
# ================================================================================================
"""
    _a0_tetrad_lowered(::Type{T}, r, θ) -> (l, n, m, mb)  (covariant legs)

The a=0 Kinnersley legs from `Spacetime.kinnersley_tetrad(0,r,θ;:BL)`, index-lowered with
`Spacetime.kerr_metric(0,r,θ;:BL)`. Grounded (not re-hardcoded); {T}-generic.
"""
function _a0_tetrad_lowered(::Type{T}, r, θ) where {T}
    lu, nu, mu, mbu = Spacetime.kinnersley_tetrad(zero(T), T(r), T(θ); coords = :BL)
    g = Spacetime.kerr_metric(zero(T), T(r), T(θ); coords = :BL)
    low(v) = Complex{T}[sum(g[μ, ν] * v[ν] for ν in 1:4) for μ in 1:4]
    (low(lu), low(nu), low(mu), low(mbu))
end

"""
    dipole_assemble_org(dc::DipoleCompletion, r, θ) -> 4×4 Complex{T}

The traceful-ORG metric perturbation `h_μν⁽¹⁾` assembled from the sole nonzero spin-1 tetrad
component `h_lm̄` (and `h_lm = conj(h_lm̄)`); the n-row `h_nm̄=h_nm=0` (ORG). Using the co-basis
`ω^{(l)}=−n, ω^{(n)}=−l, ω^{(m)}=m̄, ω^{(m̄)}=m` (null-tetrad duals, l·n=−1, m·m̄=1):

    h_μν = −h_lm̄(n_μ m_ν + m_μ n_ν) − h_lm(n_μ m̄_ν + m̄_μ n_ν).

Yields the ORG representative `h_tφ = −2δJ sin²θ/r`, `h_rφ = −2δJ sin²θ/(rf)` (= h_tφ/f, the
n·h=0 condition), all else 0. The axial gauge (below) removes `h_rφ`.
"""
function dipole_assemble_org(dc::DipoleCompletion{T}, r, θ) where {T}
    l, n, m, mb = _a0_tetrad_lowered(T, r, θ)
    hlmb = dipole_hlmbar(dc, r, θ)
    hlm = conj(hlmb)
    H = zeros(Complex{T}, 4, 4)
    for μ in 1:4, ν in 1:4
        H[μ, ν] = -hlmb * (n[μ] * m[ν] + m[μ] * n[ν]) -
                   hlm * (n[μ] * mb[ν] + mb[μ] * n[ν])
    end
    H
end

# ================================================================================================
# (ANCHOR-14) the static axial gauge vector  ξ = c(r) sin²θ dφ   (removes the ORG h_rφ, keeps h_tφ).
# ================================================================================================
"""
    dipole_gauge_c(dc::DipoleCompletion, r) -> value

The axial gauge generator `c(r)` of `ξ_φ = c(r) sin²θ` (ANCHOR-14 ξ_t log pattern), solving the
gauge ODE `c′ − 2c/r = 2δJ/(rf)` (the ORG→BL direction, which removes the residual `h_rφ`):

    c(r) = r²[ Kc + 2δJ(−¼ln r + 1/(2r) + ¼ln(r−2)) ].

(F6B_DERIVATION step-6 writes the reverse BL→ORG sign `−2δJ`; the log structure is identical.)
Only `c′ − 2c/r` (below) enters the metric; `Kc` is the **[AXIOM-F6B-5]** residual rigid-`r²` mode.
"""
function dipole_gauge_c(dc::DipoleCompletion{T}, r) where {T}
    rt = T(r)
    rt^2 * (dc.Kc + 2 * dc.δJ * (-log(rt) / 4 + 1 / (2 * rt) + log(rt - 2) / 4))
end

"""
    dipole_gauge_hrphi(dc::DipoleCompletion, r, θ) -> value

The gauge-generated `δh_rφ = (c′ − 2c/r) sin²θ = 2δJ sin²θ/(rf)`. This is the ONLY component the
static axial ξ changes (`δh_tφ = δh_θφ = δh_tt = … = 0`: `g_tφ` is gauge-rigid, F6B_DERIVATION
step-6). Adding it to the ORG assembly cancels `h_rφ^{ORG} = −2δJ sin²θ/(rf)`, recovering the BL
`d/da kerr_metric` (only `h_tφ` nonzero). Closed form (independent of `Kc`; the `r²` mode drops).
"""
function dipole_gauge_hrphi(dc::DipoleCompletion{T}, r, θ) where {T}
    rt = T(r)
    2 * dc.δJ * sin(T(θ))^2 / (rt * _f(T, rt))          # = (c′−2c/r) sin²θ,  c′−2c/r = 2δJ/(rf)
end

# ================================================================================================
# The assembled BL metric perturbation:  h_μν^{ORG} + gauge  ==  d/da kerr_metric|₀  (G-F6B-1).
# ================================================================================================
"""
    dipole_h(dc::DipoleCompletion, r, θ) -> 4×4 Complex{T}

The FULL reconstructed ℓ=1 metric perturbation in BL coordinates: the traceful-ORG assembly
[`dipole_assemble_org`](@ref) PLUS the axial-gauge Lie term (only the `rφ`/`φr` slot,
[`dipole_gauge_hrphi`](@ref)). Equals `d/da kerr_metric(δJ,r,θ;:BL)|₀` — i.e. only `h_tφ = h_φt =
−2δJ sin²θ/r` survives; every other component (including `h_rφ`) is 0. This is the DECISIVE G-F6B-1
object. Index order (t,r,θ,φ)=(1,2,3,4).
"""
function dipole_h(dc::DipoleCompletion{T}, r, θ) where {T}
    H = dipole_assemble_org(dc, r, θ)
    δhrφ = dipole_gauge_hrphi(dc, r, θ)
    H[2, 4] += δhrφ
    H[4, 2] += δhrφ
    H
end

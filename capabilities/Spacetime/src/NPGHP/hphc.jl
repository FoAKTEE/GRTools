# hphc.jl — the hyperboloidally-compactified, horizon-penetrating (HPHC) chart of TeukEvolution.jl
# (Ripley et al., arXiv:2010.00162 §V.C): radial R = cl²/r, angular (cy,sy) = (cosθ,sinθ), y = −cosθ.
# The evolved NP scalar carries its peeling factor stripped: evolved = physical × (r/cl²)^falloff =
# physical × R^(−falloff), so all coefficients are finite at scri⁺ (R=0). M = 1 (bhm ≡ 1) canon.
#
# Provides the SIX background NP_0 scalars (BackgroundNP.jl:41-57) and the ELEVEN GHP pre_* operator
# coefficient fields (GHP.jl:82-134), transcribed VERBATIM from the upstream source and independently
# cross-checked against R_teuk.md §3a/§3b (double-entry). Every form is {T}-generic with the quirk-4
# fix √2 = sqrt(T(2)) (DEVIATION-Q4: the TRUE closed form, not the upstream Float64-√2 artifact) and
# the quirk-5 fix (all of a/cl promoted into T — no `::Real` precision leak). NP_0 gated by G-B1
# (baked fixture) + G-B4 (ψ₂ chart identity); the GHP pre_* application (þ/þ′/ð/ð′ operators) stays
# consumer-side per the coefficients-vs-application principle.

using ..Spacetime: horizons, delta_kerr

# ---- chart maps R ↔ r (R = cl²/r) ----------------------------------------------------------------
"HPHC compactified radius from Boyer–Lindquist areal radius: `R = cl²/r` (Id.jl:61). R=0 ⟺ scri⁺."
hphc_R_of_r(r; cl = 1) = (T = float(promote_type(typeof(r), typeof(cl))); T(cl)^2 / T(r))

"Boyer–Lindquist areal radius from the HPHC compactified radius: `r = cl²/R`."
hphc_r_of_R(R; cl = 1) = (T = float(promote_type(typeof(R), typeof(cl))); T(cl)^2 / T(R))

"""
    hphc_horizon_R(a; cl=1) -> R₊

Outer horizon in the HPHC chart: `R₊ = cl²/r₊` with `r₊ = 1 + √(1−a²)` (M = 1). This is the
Δ-zero of the radial characteristic `A_pq` (Evolution.jl:93-94), the horizon-penetrating edge.
"""
hphc_horizon_R(a; cl = 1) = (T = float(promote_type(typeof(a), typeof(cl))); T(cl)^2 / horizons(T(a))[1])

# ---- HPHC ↔ BL-Kinnersley GHP boost-map factors (F1_DERIVATION.md) --------------------------------
"""
    hphc_boost(a, R, cy; cl=1) -> A

The GHP **BOOST** factor `A = Δ/(2Σ)` of the HPHC ↔ BL-Kinnersley tetrad map, evaluated at
Boyer–Lindquist radius `r = cl²/R` with `Δ = r²−2r+a²`, `Σ = r² + (a·cy)²` (M = 1, `cy = cosθ`).
Real; positive outside the horizon. This is the `l → A l` / `n → A⁻¹ n` rescaling
(`l_H = A·(J·l_K)`, `n_H = A⁻¹·(J·n_K)`; F1_DERIVATION.md "The law", tex:3780-3786) and the boost
weight in the NP_0 identities `ρ₀ = (A/R)·ρ_BL`, `μ₀ = μ_BL/(A·R)` (F1_DERIVATION.md
"The six NP_0 identities"). `{T}`-generic: `T` promoted from `(a, R, cy, cl)`.
"""
@inline function hphc_boost(a, R, cy; cl = 1)
    T = float(promote_type(typeof(a), typeof(R), typeof(cy), typeof(cl)))
    r = T(cl)^2 / T(R)
    Δ = delta_kerr(T(a), r)                 # r²−2r+a² (M=1)
    Σ = r^2 + (T(a) * T(cy))^2              # r²+a²cos²θ = sigma_kerr(a,r,θ), cy=cosθ
    Δ / (2 * Σ)
end

"""
    hphc_spin_phase(a, R, cy; cl=1) -> e^{iλ}

The GHP pure **SPIN-ROTATION** phase `e^{iλ} = −ζ̄/ζ` (unit modulus) of the HPHC ↔ BL-Kinnersley
map, with `ζ = r − i·a·cy` at `r = cl²/R` (`cy = cosθ`). This is the `m → e^{iλ} m` rotation
(`m_H = e^{iλ}·(J·m_K)`; F1_DERIVATION.md "The law", tex:3788-3792). Its two conjugate spin-weight
partners drive the NP_0 identities `τ₀ = −(ζ̄/ζ)·τ_BL/R²` and `π₀ = −(ζ/ζ̄)·π_BL/R²` (= conj law;
F1_DERIVATION.md "The six NP_0 identities"). On the equator (`cy = 0`) it is exactly `−1` for every
`a`. `{T}`-generic; returns `Complex{T}` with `T` promoted from `(a, R, cy, cl)`.
"""
@inline function hphc_spin_phase(a, R, cy; cl = 1)
    T = float(promote_type(typeof(a), typeof(R), typeof(cy), typeof(cl)))
    r = T(cl)^2 / T(R)
    ζ = complex(r, -T(a) * T(cy))           # ζ = r − i a cy = zeta_kerr(a, cl²/R, θ)
    -conj(ζ) / ζ
end

# ---- six background NP_0 scalars (BackgroundNP.jl:41-57 ≡ R_teuk.md §3b) --------------------------
"""
    np_spincoeff_hphc(R, cy, sy; a, cl=1) -> NamedTuple

The six background Newman–Penrose scalars of the HPHC chart — `(; mu_0, tau_0, pi_0, rho_0, eps_0,
psi2_0)` as `Complex{T}` — at compactified radius `R`, `cy = cosθ`, `sy = sinθ`, spin `a` (= bhs),
compactification `cl` (M = 1). These are the peeling-stripped (evolved) NP background of
TeukEvolution.jl (`BackgroundNP.NP_0`), reproduced as closed forms. `eps_0 ≠ 0` encodes the GHP
boost relative to strict Kinnersley. `√2 = sqrt(T(2))` (quirk-4 TRUE form). Gated by G-B1 (fixture)
and G-B4 (`psi2_0 == ψ₂_BL(cl²/R,θ)·(r/cl²)³`).
"""
function np_spincoeff_hphc(R, cy, sy; a, cl = 1)
    T = float(promote_type(typeof(R), typeof(cy), typeof(sy), typeof(a), typeof(cl)))
    RT = T(R); cyT = T(cy); syT = T(sy); bhs = T(a); clT = T(cl); bhm = one(T)
    rt2 = sqrt(T(2))
    zm = clT^2 - im * bhs * RT * cyT               # cl² − i·bhs·R·cosθ
    zp = clT^2 + im * bhs * RT * cyT               # cl² + i·bhs·R·cosθ
    d2 = clT^4 + (bhs * cyT * RT)^2                 # = real(zm·zp)

    mu_0   = 1 / (-(clT^2) + im * bhs * RT * cyT)
    tau_0  = (im * bhs * syT / rt2) / (zm^2)
    pi_0   = -(im * bhs * syT / rt2) / d2
    rho_0  = -(clT^4 - 2 * (clT^2) * bhm * RT + (bhs * RT)^2) / (2 * zm^2 * zp)
    eps_0  = (clT^2 * bhm - bhs^2 * RT - im * bhs * (clT^2 - bhm * RT) * cyT) / (2 * zm^2 * zp)
    psi2_0 = -bhm / (zm^3)
    (; mu_0 = Complex{T}(mu_0), tau_0 = Complex{T}(tau_0), pi_0 = Complex{T}(pi_0),
       rho_0 = Complex{T}(rho_0), eps_0 = Complex{T}(eps_0), psi2_0 = Complex{T}(psi2_0))
end

# ---- eleven GHP pre_* operator coefficient fields (GHP.jl:82-134 ≡ R_teuk.md §3a) -----------------
"""
    hphc_ghp_coeffs(R, cy, sy; a, cl=1) -> NamedTuple

The eleven GHP operator PRE-coefficient fields of the HPHC chart (`GHP.GHP_ops`):
`(; pre_edth_DT, pre_edth_raised, pre_edth, pre_edth_prime_DT, pre_edth_prime_lowered,
pre_edth_prime, pre_thorn_DT, pre_thorn_DR, pre_thorn, pre_thorn_prime_DT, pre_thorn_prime_DR)`,
`Complex{T}`, at `R`, `cy=cosθ`, `sy=sinθ`, spin `a`, compactification `cl` (M = 1). These are the
pointwise coefficient multipliers in the ð/ð′/þ/þ′ operator assemblies (GHP.jl:211-272); the
operator APPLICATION (SWSH raise/lower, `R·∂_R + falloff`, mode iω) is consumer-side. `√2 =
sqrt(T(2))` (quirk-4). Note the upstream `set_thorn_prime!` dead-code typo (quirk-1) and the
never-called þ′ path (quirk-2) are documented in NPGHP.jl [FUTURE]; the coefficients here are the
LIVE, cross-checked (`pre_edth = −R·tau_0`) fields.
"""
function hphc_ghp_coeffs(R, cy, sy; a, cl = 1)
    T = float(promote_type(typeof(R), typeof(cy), typeof(sy), typeof(a), typeof(cl)))
    RT = T(R); cyT = T(cy); syT = T(sy); bhs = T(a); clT = T(cl); bhm = one(T)
    rt2 = sqrt(T(2))
    zm = clT^2 - im * bhs * RT * cyT
    zp = clT^2 + im * bhs * RT * cyT
    d2 = clT^4 + (bhs * RT * cyT)^2

    pre_edth_DT             = (1 / rt2) * (1 / zm) * (-im * bhs * syT)
    pre_edth_raised         = (1 / rt2) * (1 / zm)
    pre_edth                = (im * bhs * RT * syT / rt2) / ((im * (clT^2) + bhs * RT * cyT)^2)
    pre_edth_prime_DT       = (1 / rt2) * (1 / zp) * im * bhs * syT
    pre_edth_prime_lowered  = (1 / rt2) * (1 / zp)
    pre_edth_prime          = (im * bhs * RT * syT / rt2) / (zp^2)
    pre_thorn_DT            = (1 / d2) * RT * 2 * bhm * (2 * bhm - ((bhs / clT)^2) * RT)
    pre_thorn_DR            = (1 / d2) * (-(clT^2 - 2 * bhm * RT + (bhs * RT / clT)^2) / 2)
    pre_thorn               = (1 / d2) * RT * (im * bhs)
    pre_thorn_prime_DT      = 2 + 4 * bhm * RT / (clT^2)
    pre_thorn_prime_DR      = ((1 / clT)^2) * RT
    (; pre_edth_DT = Complex{T}(pre_edth_DT), pre_edth_raised = Complex{T}(pre_edth_raised),
       pre_edth = Complex{T}(pre_edth), pre_edth_prime_DT = Complex{T}(pre_edth_prime_DT),
       pre_edth_prime_lowered = Complex{T}(pre_edth_prime_lowered),
       pre_edth_prime = Complex{T}(pre_edth_prime), pre_thorn_DT = Complex{T}(pre_thorn_DT),
       pre_thorn_DR = Complex{T}(pre_thorn_DR), pre_thorn = Complex{T}(pre_thorn),
       pre_thorn_prime_DT = Complex{T}(pre_thorn_prime_DT),
       pre_thorn_prime_DR = Complex{T}(pre_thorn_prime_DR))
end

# Puncture/derivatives.jl — analytic FIRST derivatives of the m-mode puncture (ESSF node
# `puncture-derivatives`, W6-C; arXiv:2606.04998 main.tex:481-487). The r/θ dependence of
# Ψ^P_m enters only through (Δr, Δθ, ϱ²), so the derivatives are exact chain-rule
# assemblies of the same closed forms as the puncture itself:
#
#   * ∂J_n/∂ϱ² by the exact POWER REDUCTION ∂R^n/∂ϱ² = (n/2) R^{n−2} under the integral
#     sign (main.tex:483 "∇_α ϱ^n = (n/2) ϱ^{n−2} ∇_α ϱ²") — this extends the
#     `mode_integral_Rn` range to n = −9 and, unlike the d₃P/dz-recurrence form
#     (main.tex:486, the C oracle's DrhoR route), has NO z̃-cancellation either near the
#     worldline (z̃ → ∞) or far from it (z̃ → 1⁺); the recurrence form is the
#     cross-oracle in test_puncture_derivatives.jl.
#   * ∂A_q/∂Δr, ∂A_q/∂Δθ as verbatim transcriptions of the C oracle's coefficient-group
#     derivatives (kerr-circular-comoving.c:384-388, 398-402) — the same mechanical
#     A_ijk provenance as field4d.jl/mmode_integrals.jl.
#   * t and φ̄ derivatives are the trivial ansatz factors −imΩ and +im (main.tex:487),
#     exposed on the same API so the effective-source node consumes one bundle.

"""
    mode_integral_Rn_drho2(n::Integer, m::Integer, rho2, zc2) -> real

Exact ∂/∂ϱ² of [`mode_integral_Rn`](@ref) at fixed `z_c²`, via the power reduction
`∂R^n/∂ϱ² = (n/2) R^{n−2}` under the integral sign (main.tex:483):

    ∂J_n/∂ϱ² = (n/2) · J_{n−2},   J_n = ∫_{−π}^{π} R^n e^{−imΔφ} dΔφ,

so the derivative needs `mode_integral_Rn` down to n = −9 (puncture n ∈ {−7,…,1}).
Well-conditioned everywhere the puncture is (no z̃ cancellation — the equivalent
d₃P/dz-recurrence form loses ~z̃² digits near the worldline and (z̃²−1)⁻¹ far from it;
both forms are cross-gated in the node tests). REAL, EVEN in `m`, requires `ϱ > 0`;
type-generic (BigFloat inherits full precision).
"""
mode_integral_Rn_drho2(n::Integer, m::Integer, rho2::Real, zc2::Real) =
    (mode_integral_Rn(n - 2, m, rho2, zc2) * n) / 2

# ∂A_q/∂Δr groups — VERBATIM transcription of kerr-circular-comoving.c:384-388 (dA_dr[0..4];
# same 50-coefficient provenance as `_puncture_A_groups`, pinned SHA 1cb9834f).
function _puncture_dA_dr_groups(p::Puncture, dr::Real, dtheta::Real)
    (; A106, A108, A124, A126, A142, A144, A160, A162, A180,
       A204, A206, A222, A224, A240, A242, A260,
       A304, A306, A322, A324, A340, A342, A360,
       A402, A404, A420, A422, A440,
       A502, A504, A520, A522, A540,
       A600, A602, A620, A700, A702, A720, A800, A900) = p

    dr2 = dr * dr
    dr3 = dr2 * dr
    dr4 = dr2 * dr2
    dr5 = dr3 * dr2
    dr6 = dr3 * dr3
    dr7 = dr4 * dr3
    dr8 = dr4 * dr4
    dtheta2 = dtheta * dtheta
    dtheta4 = dtheta2 * dtheta2
    dtheta6 = dtheta4 * dtheta2
    dtheta8 = dtheta4 * dtheta4

    d0 = 6 * A600 * dr5 + 7 * A700 * dr6 + 8 * A800 * dr7 + 9 * A900 * dr8 +
         4 * A420 * dr3 * dtheta2 + 5 * A520 * dr4 * dtheta2 + 6 * A620 * dr5 * dtheta2 +
         7 * A720 * dr6 * dtheta2 + 2 * A240 * dr * dtheta4 + 3 * A340 * dr2 * dtheta4 +
         4 * A440 * dr3 * dtheta4 + 5 * A540 * dr4 * dtheta4 + A160 * dtheta6 +
         2 * A260 * dr * dtheta6 + 3 * A360 * dr2 * dtheta6 + A180 * dtheta8
    d1 = 4 * A402 * dr3 + 5 * A502 * dr4 + 6 * A602 * dr5 + 7 * A702 * dr6 +
         2 * A222 * dr * dtheta2 + 3 * A322 * dr2 * dtheta2 + 4 * A422 * dr3 * dtheta2 +
         5 * A522 * dr4 * dtheta2 + A142 * dtheta4 + 2 * A242 * dr * dtheta4 +
         3 * A342 * dr2 * dtheta4 + A162 * dtheta6
    d2 = 2 * A204 * dr + 3 * A304 * dr2 + 4 * A404 * dr3 + 5 * A504 * dr4 +
         A124 * dtheta2 + 2 * A224 * dr * dtheta2 + 3 * A324 * dr2 * dtheta2 + A144 * dtheta4
    d3 = A106 + 2 * A206 * dr + 3 * A306 * dr2 + A126 * dtheta2
    d4 = A108
    return (d0, d1, d2, d3, d4)
end

# ∂A_q/∂Δθ groups — VERBATIM transcription of kerr-circular-comoving.c:398-402
# (dA_dtheta[0..4]). Every term carries an explicit Δθ factor (the A_ijk are even in Δθ),
# so the equatorial value is an EXACT zero — gated in the node tests.
function _puncture_dA_dtheta_groups(p::Puncture, dr::Real, dtheta::Real)
    (; A024, A026, A042, A044, A060, A062, A080,
       A124, A126, A142, A144, A160, A162, A180,
       A222, A224, A240, A242, A260,
       A322, A324, A340, A342, A360,
       A420, A422, A440,
       A520, A522, A540, A620, A720) = p

    dr2 = dr * dr
    dr4 = dr2 * dr2
    dtheta2 = dtheta * dtheta
    dtheta3 = dtheta2 * dtheta
    dtheta4 = dtheta2 * dtheta2
    dtheta5 = dtheta3 * dtheta2
    dtheta7 = dtheta4 * dtheta3

    d0 = 2 * (A420 + dr * (A520 + dr * (A620 + A720 * dr))) * dr4 * dtheta +
         4 * ((A240 + A340 * dr) * dr2 + (A440 + A540 * dr) * dr4 +
              (A060 + dr * (A160 + dr * (A260 + A360 * dr))) * dtheta2) * dtheta3 +
         2 * (A060 + dr * (A160 + dr * (A260 + A360 * dr))) * dtheta5 +
         8 * (A080 + A180 * dr) * dtheta7
    d1 = 2 * dtheta * ((A222 + A322 * dr) * dr2 + (A422 + A522 * dr) * dr4 +
                       2 * dtheta2 * (A042 + A142 * dr + (A242 + A342 * dr) * dr2 +
                                      (A062 + A162 * dr) * dtheta2) +
                       (A062 + A162 * dr) * dtheta4)
    d2 = 2 * dtheta * (A024 + A124 * dr + (A224 + A324 * dr) * dr2 +
                       2 * (A044 + A144 * dr) * dtheta2)
    d3 = 2 * (A026 + A126 * dr) * dtheta
    d4 = zero(d3)
    return (d0, d1, d2, d3, d4)
end

# binomial elimination c_j of Q² = (R² − ϱ²)/z_c² for ONE A-vector — the same combination
# puncture_phi_m inlines for the value (kerr-circular-comoving.c:171-176); local helper here
# because this file applies it to three A-vectors (A, ∂A/∂Δr, ∂A/∂Δθ; the landed value path
# is untouched per node scope).
function _puncture_cs(al, be, A0, A1, A2, A3, A4)
    c0 = be^4 * A0 - al * be^3 * A1 + al^2 * be^2 * A2 - al^3 * be * A3 + al^4 * A4
    c1 = be^3 * A1 - 2 * al * be^2 * A2 + 3 * al^2 * be * A3 - 4 * al^3 * A4
    c2 = be^2 * A2 - 3 * al * be * A3 + 6 * al^2 * A4
    c3 = be * A3 - 4 * al * A4
    c4 = A4
    return (c0, c1, c2, c3, c4)
end

"""
    puncture_dphi_m(p::Puncture, m::Integer, dr, dtheta) -> (phi, dr, dtheta)

The BARE azimuthal-integral bundle: `phi` = I_m(Δr, Δθ) exactly as
[`puncture_phi_m`](@ref), plus its analytic first derivatives `dr` = ∂I_m/∂Δr and
`dtheta` = ∂I_m/∂Δθ (all REAL and EVEN in `m`; the C oracle's `effsource_calc_m`
dPhiS_dr/dPhiS_dth up to the e^{−imφ_p} rotation). Chain rule with x ∈ {Δr, Δθ}
(kerr-circular-comoving.c:556-559 structure, with ∂/∂ϱ² in place of its ∂/∂ϱ):

    ∂I_m/∂x = z_c^{−8} Σ_j [ c_j(∂A/∂x) J_{n_j}  +  (∂c_j/∂ϱ² · J_{n_j}
              + c_j · ∂J_{n_j}/∂ϱ²) ∂ϱ²/∂x ],      n_j = 2j − 7,

with `∂ϱ²/∂Δr = 2 (r_p²/Δ) Δr`, `∂ϱ²/∂Δθ = 2 r_p² Δθ` (main.tex:452) and
`∂J_n/∂ϱ²` = [`mode_integral_Rn_drho2`](@ref) (n down to −9). On the equator
(`Δθ = 0`) `dtheta` is an EXACT zero (evenness of I_m). Requires `ϱ > 0`;
type-generic (BigFloat inherits full precision).
"""
function puncture_dphi_m(p::Puncture, m::Integer, dr::Real, dtheta::Real)
    o = p.orbit
    al = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dtheta^2    # ϱ² (main.tex:452)
    be = o.zc2                                            # z_c² (main.tex:453)

    A = _puncture_A_groups(p, dr, dtheta)
    c = _puncture_cs(al, be, A...)
    cr = _puncture_cs(al, be, _puncture_dA_dr_groups(p, dr, dtheta)...)
    ct = _puncture_cs(al, be, _puncture_dA_dtheta_groups(p, dr, dtheta)...)

    Jm7 = mode_integral_Rn(-7, m, al, be)
    Jm5 = mode_integral_Rn(-5, m, al, be)
    Jm3 = mode_integral_Rn(-3, m, al, be)
    Jm1 = mode_integral_Rn(-1, m, al, be)
    Jp1 = mode_integral_Rn(1, m, al, be)

    phi = (c[1] * Jm7 + c[2] * Jm5 + c[3] * Jm3 + c[4] * Jm1 + c[5] * Jp1) / be^4

    # explicit A-coefficient parts (kerr-circular-comoving.c:470-489)
    dI_dA_r = (cr[1] * Jm7 + cr[2] * Jm5 + cr[3] * Jm3 + cr[4] * Jm1 + cr[5] * Jp1) / be^4
    dI_dA_t = (ct[1] * Jm7 + ct[2] * Jm5 + ct[3] * Jm3 + ct[4] * Jm1 + ct[5] * Jp1) / be^4

    # ∂c_j/∂ϱ² at fixed A_q (kerr-circular-comoving.c:499-503 dRePhiS_dalpha; ∂c₄/∂ϱ² = 0)
    dc0 = -be^3 * A[2] + 2 * al * be^2 * A[3] - 3 * al^2 * be * A[4] + 4 * al^3 * A[5]
    dc1 = -2 * be^2 * A[3] + 6 * al * be * A[4] - 12 * al^2 * A[5]
    dc2 = -3 * be * A[4] + 12 * al * A[5]
    dc3 = -4 * A[5]

    # full ∂I_m/∂ϱ²: coefficient part + power-reduced integral part (c.f. c:505-509 drho)
    dI_du = (dc0 * Jm7 + dc1 * Jm5 + dc2 * Jm3 + dc3 * Jm1 +
             c[1] * mode_integral_Rn_drho2(-7, m, al, be) +
             c[2] * mode_integral_Rn_drho2(-5, m, al, be) +
             c[3] * mode_integral_Rn_drho2(-3, m, al, be) +
             c[4] * mode_integral_Rn_drho2(-1, m, al, be) +
             c[5] * mode_integral_Rn_drho2(1, m, al, be)) / be^4

    du_dr = 2 * o.rho2_dr2 * dr
    du_dt = 2 * o.rho2_dtheta2 * dtheta

    return (phi = phi, dr = dI_dA_r + dI_du * du_dr, dtheta = dI_dA_t + dI_du * du_dt)
end

"""
    puncture_dpsi_m(p::Puncture, m::Integer, r, z; varphi_p = 0)
        -> (psi, dt, dr, dtheta, dphi)

THE m-mode puncture derivative bundle the effective source consumes: `psi` is EXACTLY the
landed [`puncture_psi_m`](@ref) (same call, bitwise — single knowledge home of the
assembly), and the four first derivatives follow the LOCKED mode-ansatz convention
(main.tex:481-487 applied to `Ψ^P_m = r e^{−imφ̄_p}/(2π sin^{|m|}θ) · I_m`):

    dt     = −imΩ · psi,                dphi = +im · psi        (trivial ansatz factors),
    dr     = e^{−imφ̄_p}/(2π sin^{|m|}θ) · (I_m + r ∂I_m/∂Δr),
    dtheta = r e^{−imφ̄_p}/(2π) · (∂I_m/∂Δθ / sin^{|m|}θ + I_m sinθ · d(sin^{|m|}θ)/dz / sin^{2|m|}θ),

where the `dtheta` ansatz-factor term uses `d(sin^{|m|}θ)/dθ = −sinθ · sinm_theta_dz`
(single home of the ansatz-factor derivative; exact zero at m = 0). `dt`/`dphi` are
derivatives of the FULL ansatz term ∂_t, ∂_φ̄ brought to mode level — time-independence
of the co-rotating puncture (φ̄_p(t) = φ̄_p + Ωt) is inherited from `puncture_psi_m`.
Conjugate symmetry: every slot of the `−m` bundle is the conjugate of the `+m` slot.
Valid at `|z| < 1` for m ≠ 0 and off the worldline circle (ϱ > 0); type-generic
(BigFloat inputs give `Complex{BigFloat}` slots).
"""
function puncture_dpsi_m(p::Puncture, m::Integer, r::Real, z::Real; varphi_p::Real = 0)
    o = p.orbit
    b = puncture_dphi_m(p, m, r - o.r0, -asin(z))         # Δθ = θ − π/2 = −asin z (locked)
    psi = puncture_psi_m(p, m, r, z; varphi_p = varphi_p)
    pref = cis(-m * varphi_p) / (2 * oftype(b.phi, pi))
    s = sinm_theta(m, z)
    dpsi_dr = pref * (b.phi + r * b.dr) / s
    dpsi_dtheta = r * pref * (b.dtheta / s +
                              b.phi * sintheta_of_z(z) * sinm_theta_dz(m, z) / s^2)
    return (psi = psi, dt = -im * m * o.Omega * psi, dr = dpsi_dr,
            dtheta = dpsi_dtheta, dphi = im * m * psi)
end

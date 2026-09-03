# Puncture/second_derivatives.jl — analytic SECOND derivatives of the m-mode puncture
# (ESSF node `effective-source-direct`, W7-C; additive extension of the W6-C first-derivative
# layer; arXiv:2606.04998 main.tex:481-490). The effective source is assembled from
# ∂²Ψ^P_m/∂r² and ∂²Ψ^P_m/∂z² (the first-order flux form is DIAGONAL in (r, z), so no mixed
# second derivative is ever needed — the same five slots the C oracle records as NaN,
# kerr-circular-comoving.c:570-574).
#
# The bare azimuthal integral is I_m(Δr, Δθ) = z_c^{−8} Σ_{j=0}^{4} c_j J_{n_j}, n_j = 2j−7,
# with the binomial-elimination coefficients (mmode_integrals.jl `_puncture_cs`)
#
#     c_j(B, u) = Σ_{q≥j} binom(q, j) (−u)^{q−j} z_c^{2(4−q)} B_q ,   u ≡ ϱ² ,
#
# so I_m = Ĝ(A(Δr, Δθ), u(Δr, Δθ)) with Ĝ LINEAR in the A-vector. Two exact identities carry
# the whole second-derivative chain, and neither adds new coefficient algebra:
#
#   (i)  ∂c_j/∂u = −(j+1) c_{j+1}   (c_5 ≡ 0), hence ∂²c_j/∂u² = (j+1)(j+2) c_{j+2}.
#        Proof: binom(q,j)(q−j) = (j+1) binom(q,j+1) termwise. The SEALED `_puncture_cs` is
#        therefore the only home of the c_j algebra — the W6-C `puncture_dphi_m` ∂c_j/∂ϱ²
#        block is this identity written out.
#   (ii) ∂J_n/∂u = (n/2) J_{n−2} and ∂²J_n/∂u² = (n/2)((n−2)/2) J_{n−4}, the exact power
#        reduction ∂R^n/∂ϱ² = (n/2)R^{n−2} under the integral sign (main.tex:483) already
#        sealed at first order as `mode_integral_Rn_drho2`. The second order extends
#        `mode_integral_Rn` to n = −11 (W4-D reached −7, W6-C −9).
#
# With Ĝ_u, Ĝ_uu the u-derivatives of Ĝ at FIXED A-vector, x ∈ {Δr, Δθ}:
#
#     ∂²I_m/∂x² = Ĝ(∂²A/∂x², u) + 2 Ĝ_u(∂A/∂x, u) u_x + Ĝ_uu(A, u) u_x² + Ĝ_u(A, u) u_xx ,
#
# with u_x = 2·(ϱ²-prefactor)·x and u_xx = 2·(ϱ²-prefactor) exactly (main.tex:452, the
# prefactors cached on `CircularOrbit`). The only new transcription is the A-group second
# derivatives, taken VERBATIM from the C oracle (pinned SHA 1cb9834f) exactly as
# field4d.jl / mmode_integrals.jl / derivatives.jl took the value and first-derivative groups.

"""
    mode_integral_Rn_d2rho2(n::Integer, m::Integer, rho2, zc2) -> real

Exact ∂²/∂(ϱ²)² of [`mode_integral_Rn`](@ref) at fixed `z_c²`, by applying the power
reduction `∂R^n/∂ϱ² = (n/2) R^{n−2}` (main.tex:483) twice under the integral sign:

    ∂²J_n/∂(ϱ²)² = (n/2)·((n−2)/2) · J_{n−4},   J_n = ∫_{−π}^{π} R^n e^{−imΔφ} dΔφ,

so the effective source needs `mode_integral_Rn` down to `n = −11` (the puncture itself uses
n ∈ {−7,…,1}, its first derivatives −9). Companion of the sealed
[`mode_integral_Rn_drho2`](@ref) and equally free of z̃ cancellation. REAL, EVEN in `m`,
requires `ϱ > 0`; type-generic (BigFloat inherits full precision).
"""
mode_integral_Rn_d2rho2(n::Integer, m::Integer, rho2::Real, zc2::Real) =
    (mode_integral_Rn(n - 4, m, rho2, zc2) * (n * (n - 2))) / 4

# ∂²A_q/∂Δr² groups — VERBATIM transcription of kerr-circular-comoving.c:391-395
# (d2A_dr2[0..4]; same 50-coefficient provenance as `_puncture_A_groups`, SHA 1cb9834f).
function _puncture_d2A_dr2_groups(p::Puncture, dr::Real, dtheta::Real)
    (; A204, A206, A222, A224, A240, A242, A260,
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
    dtheta2 = dtheta * dtheta
    dtheta4 = dtheta2 * dtheta2
    dtheta6 = dtheta4 * dtheta2

    d0 = 2 * (15 * A600 * dr4 + 21 * A700 * dr5 + 28 * A800 * dr6 + 36 * A900 * dr7 +
              6 * A420 * dr2 * dtheta2 + 10 * A520 * dr3 * dtheta2 +
              15 * A620 * dr4 * dtheta2 + 21 * A720 * dr5 * dtheta2 +
              A240 * dtheta4 + 3 * A340 * dr * dtheta4 + 6 * A440 * dr2 * dtheta4 +
              10 * A540 * dr3 * dtheta4 + A260 * dtheta6 + 3 * A360 * dr * dtheta6)
    d1 = 2 * (6 * A402 * dr2 + 10 * A502 * dr3 + 15 * A602 * dr4 + 21 * A702 * dr5 +
              A222 * dtheta2 + 3 * A322 * dr * dtheta2 + 6 * A422 * dr2 * dtheta2 +
              10 * A522 * dr3 * dtheta2 + A242 * dtheta4 + 3 * A342 * dr * dtheta4)
    d2 = 2 * (A204 + 3 * A304 * dr + 6 * A404 * dr2 + 10 * A504 * dr3 +
              A224 * dtheta2 + 3 * A324 * dr * dtheta2)
    d3 = 2 * (A206 + 3 * A306 * dr)
    d4 = zero(d3)
    return (d0, d1, d2, d3, d4)
end

# ∂²A_q/∂Δθ² groups — VERBATIM transcription of kerr-circular-comoving.c:405-409
# (d2A_dtheta2[0..4]). Unlike the FIRST θ-derivative these do NOT all carry an explicit Δθ
# factor, so the equatorial value is generally NONZERO — it carries the leading angular
# curvature of the effective source on the equator.
function _puncture_d2A_dtheta2_groups(p::Puncture, dr::Real, dtheta::Real)
    (; A024, A026, A042, A044, A060, A062, A080,
       A124, A126, A142, A144, A160, A162, A180,
       A222, A224, A240, A242, A260,
       A322, A324, A340, A342, A360,
       A420, A422, A440, A520, A522, A540, A620, A720) = p

    dr2 = dr * dr
    dr4 = dr2 * dr2
    dtheta2 = dtheta * dtheta
    dtheta4 = dtheta2 * dtheta2
    dtheta6 = dtheta4 * dtheta2

    d0 = 2 * ((A420 + dr * (A520 + dr * (A620 + A720 * dr))) * dr4 +
              6 * dtheta2 * ((A240 + A340 * dr) * dr2 + (A440 + A540 * dr) * dr4 +
                             (A060 + dr * (A160 + dr * (A260 + A360 * dr))) * dtheta2) +
              9 * (A060 + dr * (A160 + dr * (A260 + A360 * dr))) * dtheta4 +
              28 * (A080 + A180 * dr) * dtheta6)
    d1 = 2 * ((A222 + A322 * dr) * dr2 + (A422 + A522 * dr) * dr4 +
              6 * dtheta2 * (A042 + A142 * dr + (A242 + A342 * dr) * dr2 +
                             (A062 + A162 * dr) * dtheta2) +
              9 * (A062 + A162 * dr) * dtheta4)
    d2 = 2 * (A024 + A124 * dr + (A224 + A324 * dr) * dr2 +
              6 * (A044 + A144 * dr) * dtheta2)
    d3 = 2 * (A026 + A126 * dr)
    d4 = zero(d3)
    return (d0, d1, d2, d3, d4)
end

# Ĝ(B, u) and its two u-derivatives at FIXED A-vector B, from ONE evaluation of the J ladder.
# `c` = `_puncture_cs(u, z_c², B...)`; `J[k]` = J_{n_k}, `Jd[k]` = ∂J_{n_k}/∂u,
# `Jdd[k]` = ∂²J_{n_k}/∂u² with n_k = 2(k−1)−7. The ∂c_j/∂u = −(j+1)c_{j+1} ladder (file
# header identity (i)) is what makes the u-derivatives pure re-indexing of the SAME c vector.
_puncture_ghat(c, J, be) =
    (c[1] * J[1] + c[2] * J[2] + c[3] * J[3] + c[4] * J[4] + c[5] * J[5]) / be^4

_puncture_ghat_u(c, J, Jd, be) =
    (-c[2] * J[1] - 2 * c[3] * J[2] - 3 * c[4] * J[3] - 4 * c[5] * J[4] +
     c[1] * Jd[1] + c[2] * Jd[2] + c[3] * Jd[3] + c[4] * Jd[4] + c[5] * Jd[5]) / be^4

_puncture_ghat_uu(c, J, Jd, Jdd, be) =
    (2 * c[3] * J[1] + 6 * c[4] * J[2] + 12 * c[5] * J[3] -
     2 * (c[2] * Jd[1] + 2 * c[3] * Jd[2] + 3 * c[4] * Jd[3] + 4 * c[5] * Jd[4]) +
     c[1] * Jdd[1] + c[2] * Jdd[2] + c[3] * Jdd[3] + c[4] * Jdd[4] + c[5] * Jdd[5]) / be^4

"""
    puncture_d2phi_m(p::Puncture, m::Integer, dr, dtheta)
        -> (phi, dr, dtheta, drr, dthetatheta)

The FUSED bare azimuthal-integral bundle of the m-mode puncture through SECOND order — the
effective source's single puncture call (main.tex:459-490). All five slots are REAL and EVEN
in `m`:

  * `phi`, `dr`, `dtheta` — bitwise/rounding-identical to the sealed [`puncture_phi_m`](@ref)
    and [`puncture_dphi_m`](@ref) (the value slot is the SAME expression over the same `c_j`
    and `J_n`; the two first derivatives regroup the ∂/∂ϱ² sum through the c_j ladder and so
    agree to rounding — gated at 1e-13 in the node tests);
  * `drr` = ∂²I_m/∂Δr², `dthetatheta` = ∂²I_m/∂Δθ² — the new content, assembled as

        ∂²I_m/∂x² = Ĝ(∂²A/∂x²) + 2 Ĝ_u(∂A/∂x) u_x + Ĝ_uu(A) u_x² + Ĝ_u(A) u_xx

    (file header; `u = ϱ²`, `u_x`/`u_xx` exact from main.tex:452) over ONE evaluation of the
    seven-member J ladder `J_{−11..1}`, so the whole bundle costs the same seven
    [`mode_integral_Rn`](@ref) calls the value alone would.

NO MIXED second derivative is produced: the first-order flux form is diagonal in `(r, z)`, so
`∂²/∂Δr∂Δθ` never enters the effective source — the same five slots the C oracle leaves as
NaN. On the equator (`dtheta == 0`) `dtheta` is an EXACT zero (θ-evenness of `I_m`) while
`dthetatheta` is not. Requires `ϱ > 0`; type-generic (BigFloat inherits full precision).
"""
function puncture_d2phi_m(p::Puncture, m::Integer, dr::Real, dtheta::Real)
    o = p.orbit
    al = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dtheta^2    # ϱ² (main.tex:452)
    be = o.zc2                                            # z_c² (main.tex:453)

    c   = _puncture_cs(al, be, _puncture_A_groups(p, dr, dtheta)...)
    cr  = _puncture_cs(al, be, _puncture_dA_dr_groups(p, dr, dtheta)...)
    ct  = _puncture_cs(al, be, _puncture_dA_dtheta_groups(p, dr, dtheta)...)
    crr = _puncture_cs(al, be, _puncture_d2A_dr2_groups(p, dr, dtheta)...)
    ctt = _puncture_cs(al, be, _puncture_d2A_dtheta2_groups(p, dr, dtheta)...)

    # the J ladder: the value needs n ∈ {−7,−5,−3,−1,1}; ∂/∂ϱ² pulls n−2, ∂²/∂(ϱ²)² pulls n−4
    Jm11 = mode_integral_Rn(-11, m, al, be)
    Jm9  = mode_integral_Rn(-9, m, al, be)
    Jm7  = mode_integral_Rn(-7, m, al, be)
    Jm5  = mode_integral_Rn(-5, m, al, be)
    Jm3  = mode_integral_Rn(-3, m, al, be)
    Jm1  = mode_integral_Rn(-1, m, al, be)
    Jp1  = mode_integral_Rn(1, m, al, be)
    J   = (Jm7, Jm5, Jm3, Jm1, Jp1)
    Jd  = ((-7 * Jm9) / 2, (-5 * Jm7) / 2, (-3 * Jm5) / 2, -Jm3 / 2, Jm1 / 2)
    Jdd = ((63 * Jm11) / 4, (35 * Jm9) / 4, (15 * Jm7) / 4, (3 * Jm5) / 4, -Jm3 / 4)

    u_r  = 2 * o.rho2_dr2 * dr
    u_rr = 2 * o.rho2_dr2
    u_t  = 2 * o.rho2_dtheta2 * dtheta
    u_tt = 2 * o.rho2_dtheta2

    gu  = _puncture_ghat_u(c, J, Jd, be)
    guu = _puncture_ghat_uu(c, J, Jd, Jdd, be)

    return (phi = _puncture_ghat(c, J, be),
            dr = _puncture_ghat(cr, J, be) + gu * u_r,
            dtheta = _puncture_ghat(ct, J, be) + gu * u_t,
            drr = _puncture_ghat(crr, J, be) +
                  2 * _puncture_ghat_u(cr, J, Jd, be) * u_r + guu * u_r^2 + gu * u_rr,
            dthetatheta = _puncture_ghat(ctt, J, be) +
                          2 * _puncture_ghat_u(ct, J, Jd, be) * u_t + guu * u_t^2 +
                          gu * u_tt)
end

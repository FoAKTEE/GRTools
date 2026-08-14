# Puncture/mmode_integrals.jl — the m-mode puncture via the CRP ₃P route (ESSF node
# `mmode-puncture-integrals`, W4-D; arXiv:2606.04998 Sec. II.C "Mode decomposition",
# main.tex:459-486; Bourg et al. arXiv:2404.10082 Eq. (115) = eq. I00nm there).
#
# Route: each A_ijk Δr^i Δθ^j Q^k / R⁷ term of the 4D puncture (eq:puncture) has its
# azimuthal integral reduced by eliminating Q via Q² = (R² − ϱ²)/z_c² (eq:s2), leaving
# integrals ∫ R^n e^{−imΔφ} dΔφ with n ODD ∈ {−7,−5,−3,−1,1}, each evaluated by the exact
# closed form of Bourg Eq. (115) over `BigSpecial.legendre_P3` (half-integer degree; the
# CRP's own z ≥ 0.995 DLMF 15.8.10/15.8.12 near-worldline switch lives inside it) and
# `BigSpecial.pochhammer`. This replaces the Thornburg (2017) complete-elliptic-integral ×
# m-dependent-polynomial route, whose polynomials lose digits far from the worldline at
# high m (main.tex:465-466) — adjudicated in test_mmode_puncture.jl: at (Δr = 10, m = 20)
# the E/K-polynomial C oracle is off by ≥ 7e1 relative while this route tracks a
# BigFloat(256) quadrature of the 4D puncture to ≤ 5.7e-11.

"""
    mode_integral_Rn(n::Integer, m::Integer, rho2, zc2) -> real

The azimuthal mode integral of an odd power of the puncture distance function
(arXiv:2606.04998 main.tex:461-473; Bourg arXiv:2404.10082 Eq. (115)):

    ∫_{−π}^{π} R^n e^{−i m Δφ} dΔφ
      = 2π (−1)^m / (1 + n/2)_m · [ϱ √(ϱ² + z_c²)]^{n/2} · ₃P^m_{n/2}(z̃),

with `R² = ϱ² + z_c² Q²`, `Q = sin(Δφ/2)` (eq:s2), `rho2 = ϱ² > 0`, `zc2 = z_c² > 0`, and
`z̃ = (ϱ² + z_c²/2)/(ϱ√(ϱ² + z_c²)) > 1` (AM–GM, margin z_c⁴/4). The result is REAL (the
integrand's imaginary part is odd in Δφ) and EVEN in `m` — the `m < 0` route through the
reciprocal-falling Pochhammer and the negative-order `legendre_P3` connection reproduces
the `+m` value exactly (Γ-function identity, gated at 1e-13). `n` is any odd integer
(`1 + n/2` half-integer ⇒ the Pochhammer never vanishes); the puncture needs
n ∈ {−7,…,1}, its r/θ derivatives push to n = −9…3 and the identity is tested on
[−7, 9]. Negative half-integer degrees route through ₃P_ν = ₃P_{−ν−1} inside
`legendre_P3`. Requires `ϱ > 0`: on the worldline circle the closed form (and for
n ≤ −1 the integral itself) is singular. Type-generic: BigFloat arguments inherit full
precision through the BigSpecial primitives.
"""
function mode_integral_Rn(n::Integer, m::Integer, rho2::Real, zc2::Real)
    T = float(promote_type(typeof(rho2), typeof(zc2)))
    rho2T, zc2T = T(rho2), T(zc2)
    rR0 = sqrt(rho2T) * sqrt(rho2T + zc2T)          # ϱ√(ϱ² + z_c²)
    zt = (rho2T + zc2T / 2) / rR0
    sgn = iseven(m) ? one(T) : -one(T)
    half_n = T(n) / 2
    return 2 * T(pi) * sgn * rR0^half_n * legendre_P3(half_n, T(m), zt) /
           pochhammer(1 + half_n, m)
end

# Q^{2q}-grouped A_ijk polynomial: A(Δr, Δθ, Q) = Σ_{q=0}^{4} A_q(Δr, Δθ) Q^{2q} — the
# eq:puncture numerator regrouped by azimuthal power (VERBATIM regrouping of the C oracle,
# kerr-circular-comoving.c:158-162; same 50 coefficients as `puncture_phi`, consistency
# gated in test_mmode_puncture.jl at 1e-13). Second home of the 50-coefficient polynomial
# knowledge (with puncture_phi's flat form) — kept as the C-auditable transcription rather
# than abstracted, per the field4d.jl precedent.
function _puncture_A_groups(p::Puncture, dr::Real, dtheta::Real)
    (; A006, A008, A024, A026, A042, A044, A060, A062, A080,
       A106, A108, A124, A126, A142, A144, A160, A162, A180,
       A204, A206, A222, A224, A240, A242, A260,
       A304, A306, A322, A324, A340, A342, A360,
       A402, A404, A420, A422, A440,
       A502, A504, A520, A522, A540,
       A600, A602, A620, A700, A702, A720, A800, A900) = p

    dr2 = dr * dr
    dr3 = dr2 * dr
    dr4 = dr2 * dr2
    dr6 = dr3 * dr3
    dr8 = dr4 * dr4
    dtheta2 = dtheta * dtheta
    dtheta4 = dtheta2 * dtheta2
    dtheta8 = dtheta4 * dtheta4

    A0 = dr6 * (A600 + A700 * dr) + dr8 * (A800 + A900 * dr) +
         dr4 * (A420 + A520 * dr + dr2 * (A620 + A720 * dr)) * dtheta2 +
         (A080 + A180 * dr) * dtheta8 +
         dtheta4 * (dr2 * (A240 + A340 * dr) + dr4 * (A440 + A540 * dr) +
                    (A060 + A160 * dr + dr2 * (A260 + A360 * dr)) * dtheta2)
    A1 = dr4 * (A402 + A502 * dr + dr2 * (A602 + A702 * dr)) +
         (dr2 * (A222 + A322 * dr) + dr4 * (A422 + A522 * dr)) * dtheta2 +
         dtheta4 * (A042 + A142 * dr + dr2 * (A242 + A342 * dr) +
                    (A062 + A162 * dr) * dtheta2)
    A2 = dr2 * (A204 + A304 * dr) + dr4 * (A404 + A504 * dr) +
         (A024 + A124 * dr + dr2 * (A224 + A324 * dr)) * dtheta2 +
         (A044 + A144 * dr) * dtheta4
    A3 = A006 + A106 * dr + dr2 * (A206 + A306 * dr) + (A026 + A126 * dr) * dtheta2
    A4 = A008 + A108 * dr
    return (A0, A1, A2, A3, A4)
end

"""
    puncture_phi_m(p::Puncture, m::Integer, dr, dtheta) -> real

The BARE azimuthal integral of the 4D Detweiler–Whiting puncture — NOT yet the solver's
mode function (that is [`puncture_psi_m`](@ref)):

    I_m(Δr, Δθ) = ∫_{−π}^{π} Φ^P(Δr, Δθ, Δφ) e^{−i m Δφ} dΔφ,

REAL (Φ^P is even in Δφ) and EVEN in `m`; carries NO 1/(2π), NO e^{−imφ̄_p} particle
phase, and NO sin^{|m|}θ/r ansatz factors. This is exactly the C oracle's
`effsource_PhiS_m` output up to its e^{−imφ_p} rotation (both variants; fixture headers
in test/fixtures/effsource/).

Evaluation (main.tex:459-473): with the numerator regrouped as Σ_q A_q Q^{2q}
([`_puncture_A_groups`](@ref)), `Q² = (R² − ϱ²)/z_c²` (eq:s2) turns each term into a
binomial combination of [`mode_integral_Rn`](@ref) values with n odd ∈ {−7,…,1}:

    I_m = z_c^{−8} Σ_{j=0}^{4} c_j ∫R^{2j−7}e^{−imΔφ}dΔφ,
    c_j = Σ_{q≥j} binom(q,j) (−ϱ²)^{q−j} z_c^{2(4−q)} A_q,

with ϱ², z_c² consumed from the orbit's solid puncture scales (`rho2_dr2`,
`rho2_dtheta2`, `zc2`). Conditioning (measured, node evidence): the alternating c_j
combination is benign near the worldline (ϱ² ≪ z_c²; loss ~1e-16) and loses at most
~6e-11 at the grid extreme ϱ²/z_c² ~ 1e2 with m = 20 — far better than either C fixture
route in its bad window. Requires ϱ > 0 (off the worldline circle). Type-generic; the
worldline-crossing sign structure (odd Δr powers) is inherited from the A_ijk table.
"""
function puncture_phi_m(p::Puncture, m::Integer, dr::Real, dtheta::Real)
    A0, A1, A2, A3, A4 = _puncture_A_groups(p, dr, dtheta)
    o = p.orbit
    al = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dtheta^2    # ϱ² (main.tex:452)
    be = o.zc2                                            # z_c² (main.tex:453)

    # binomial elimination of Q² = (R² − ϱ²)/z_c² (kerr-circular-comoving.c:171-176 order)
    c0 = be^4 * A0 - al * be^3 * A1 + al^2 * be^2 * A2 - al^3 * be * A3 + al^4 * A4
    c1 = be^3 * A1 - 2 * al * be^2 * A2 + 3 * al^2 * be * A3 - 4 * al^3 * A4
    c2 = be^2 * A2 - 3 * al * be * A3 + 6 * al^2 * A4
    c3 = be * A3 - 4 * al * A4
    c4 = A4

    return (c0 * mode_integral_Rn(-7, m, al, be) +
            c1 * mode_integral_Rn(-5, m, al, be) +
            c2 * mode_integral_Rn(-3, m, al, be) +
            c3 * mode_integral_Rn(-1, m, al, be) +
            c4 * mode_integral_Rn(1, m, al, be)) / be^4
end

"""
    puncture_psi_m(p::Puncture, m::Integer, r, z; varphi_p = 0) -> Complex

THE m-mode puncture the elliptic solver consumes: `Ψ^P_m(r, z)` per the mode convention
LOCKED in MModeOperator/ansatz.jl (eq:modeIntegral, main.tex:233-236, applied to Φ^P) —
the sin^{|m|}θ/r ansatz factoring IS included (stripped through `mode_amplitude`), so the
output is the mode FUNCTION, not the bare azimuthal integral [`puncture_phi_m`](@ref):

    Ψ^P_m(r, z) = r / (2π sin^{|m|}θ) · e^{−i m φ̄_p} · I_m(r − r_p, θ − π/2),

with `z = cosθ` (θ − π/2 = −asin z, the locked branch) and `varphi_p` = the particle's
azimuth in the horizon-regular chart φ̄ at t = 0 (the baked fixtures use φ_p = 0.31).
Time-independence: the co-rotating puncture has φ̄_p(t) = φ̄_p + Ωt, whose e^{−imΩt}
cancels eq:modeIntegral's e^{+imΩt} — gated against the ts ≠ 0 assembly in the tests.
Slicing: the puncture lives in the worldtube/source region where H = 0 and h = 0
(main.tex:280), so Ψ^P_m = Ψ̃^P_m there and no e^{−imΩh} factor applies.

Conjugate symmetry (real Φ^P): `Ψ^P_{−m} = conj(Ψ^P_m)` — automatic from the m-evenness
of I_m plus the |m| ansatz factors. Valid at |z| < 1 for m ≠ 0 (poles) and off the
worldline circle (ϱ > 0); diverges logarithmically towards it (main.tex:424).
Type-generic: BigFloat inputs yield `Complex{BigFloat}`.
"""
function puncture_psi_m(p::Puncture, m::Integer, r::Real, z::Real; varphi_p::Real = 0)
    o = p.orbit
    Im_ = puncture_phi_m(p, m, r - o.r0, -asin(z))
    cm = Im_ * cis(-m * varphi_p) / (2 * oftype(Im_, pi))
    return mode_amplitude(m, cm; Omega = o.Omega, ts = zero(Im_), r = r, z = z)
end

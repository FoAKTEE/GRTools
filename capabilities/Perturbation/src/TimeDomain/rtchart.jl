# Rácz–Tóth (RT) horizon-penetrating hyperboloidal chart for the GENERIC-SPIN time-domain
# Teukolsky equation — Harms, Bernuzzi & Brügmann [arXiv:1301.1591] (DG-P4 chart decision).
#
# WHY THIS CHART: in the s=0 hyperboloidal-layer chart the s≠0 frictions decay only as 2s/r and
# the E-inversion produces coefficients ~ r, DIVERGENT at scri — the layer chart cannot carry
# s≠0. The RT chart with Ψ = Δ⁻ˢ r⁻¹ ψ [1301.1591 Eq. (12)] keeps every coefficient a rational
# polynomial regular on R ∈ [R₊, 1], and both boundaries are pure outflow (C_RR ∝ Δ·(R²−1)⁴
# vanishes at both ends; §2.3 Eqs. (16)–(19)). The s=0 layer-chart path stays bit-untouched.
#
# CHART [arXiv:1301.1591, Eqs. (3)–(6); classic RT parameters κ = S = 1 of Eqs. (8)–(9),
# which is what the appendix coefficients assume]:
#   ingoing-Kerr:  τ = t − r + ∫dr (r²+a²)/Δ,   φ̃ = φ + a∫dr/Δ        [Eqs. (3),(4)]
#   RT:            τ = T − 4M log(1−R²) + (1+R²)/(1−R²),  r = 2R/(1−R²)  [Eqs. (5),(6)]
# scri⁺ at R = 1, horizon at R₊ [Eq. (7)]. Master-variable rescaling and m-mode:
#   Ψ = Δ⁻ˢ r⁻¹ ψ,   ψ(T,R,θ,φ̃) = Σₘ ψₘ(T,R,θ) e^{imφ̃}                 [Eqs. (12),(15)]
# (Same asymptotic weight r^{−2s−1} as Barack–Giudice's (rΔˢ)⁻¹ [arXiv:1702.04204, Eq. (27)];
# at a=0 the two mode variables COINCIDE, ψ_RT,ℓm = φ_BG,sℓm — used by the Bardeen–Press
# gate D.)
#
# COEFFICIENT TRANSCRIPTION: every C_X below is the appendix entry of [arXiv:1301.1591],
# with the θ-dependence split off analytically:
#   C_T  = rt_CT_flat(R)  + 2ias·cosθ/R                       (cosθ ⇒ ℓ±1 coupling)
#   C_TT = rt_CTT_flat(R) + a²·cos2θ/(2R)                     (cos2θ = 1−2sin²θ ⇒ ℓ,ℓ±1,ℓ±2)
#   C_φ  = a(1−R²)/R² − 2is·cotθcscθ/R                        (Im part ⇒ angular operator)
#   C_0  = rt_C0_flat(R) + (s²cot²θ − s)/R                    (θ part ⇒ angular operator)
#   C_θ = −cotθ/R, C_θθ = −1/R, C_φφ = −csc²θ/R               (⇒ angular operator)
# The transcription (including the splits) is gated against an independent sympy re-derivation
# (test/fixtures/rt_coeffs_reference.jl). NOTE: [1301.1591 Eq. (10)] has a typo — "(r²+a²)Δ⁻¹"
# should read "(r²+a²)²Δ⁻¹".
#
# 1+1 ℓ-COUPLED PROJECTION (not in print — Barack–Giudice re-expansions [1702.04204,
# Eqs. (32)–(35)] applied to the 2+1 m-mode equation [1301.1591, Eq. (15)]):
# insert ψₘ = Σ_ℓ ψ_{ℓm}(T,R) ₛYₗₘ(θ) (φ̃-stripped) and project on ₛȲ_{ℓ′m}. The angular
# collection −(1/R)[∂θθ + cotθ∂θ − (m²+2ms·cosθ)/sin²θ − s²cot²θ + s] acting on ₛYₗₘ gives
# +(ℓ−s)(ℓ+s+1)/R by the spin-weighted-harmonic equation [arXiv:1702.04204, Eq. (25)];
# the cosθ / sin²θ factors give the 𝒦 / 𝒮 matrices of spinweight.jl. Result, per (s,m):
#
#   𝒜_TT ∂_TT ψ⃗ + C_TR ∂_TR ψ⃗ + C_RR ∂_RR ψ⃗ + 𝒜_T ∂_T ψ⃗ + C̃_R ∂_R ψ⃗ + diag(C̄₀^{(ℓ)}) ψ⃗ = 0
#
#   𝒜_TT = (rt_CTT_flat + a²/(2R))·I − (a²/R)·𝒮          (real, symmetric, PENTADIAGONAL)
#   𝒜_T  = (rt_CT_flat + im·C_Tφ)·I + (2ias/R)·𝒦         (complex, tridiagonal)
#   C̃_R  = C_R + im·C_Rφ,  C̃_TR = C_TR,  C̃_RR = C_RR     (scalars; principal part REAL)
#   C̄₀^{(ℓ)} = rt_C0_flat + im·a(1−R²)/R² + (ℓ−s)(ℓ+s+1)/R
#
# s=0 reduction: the ℓ±1 bands vanish, 𝒮 reduces to the sin²θ-coupling, and the equation is the
# RT scalar equation [1301.1591 §2.2] — gated against the layer-chart KerrModes path.

# rt_r_of_R and rt_Rplus are pure chart maps homed in Spacetime.RTChart, re-imported by the
# TimeDomain module header; call sites below resolve to the same Spacetime objects.

"""
    rt_C0_flat(R; M=1, a=0, s=0)

θ-independent part of the appendix coefficient C₀ of [arXiv:1301.1591, Eq. (14)]:
`(−a²(R²−1)²/2 − M(R²−1)R(s+1))/R³`. The remaining `R²s(s·cot²θ−1)/R³ = (s²cot²θ−s)/R`
belongs to the angular operator (consumed by the ₛYₗₘ eigenvalue, header).
"""
rt_C0_flat(R; M=1.0, a=0.0, s::Integer=0) =
    (-(a^2 / 2) * (R^2 - 1)^2 - M * (R^2 - 1) * R * (s + 1)) / R^3

"""
    rt_CT_flat(R; M=1, a=0, s=0)

θ-independent part of the appendix coefficient C_T of [arXiv:1301.1591, Eq. (14)] (the
`2ia(R+1)(R²+1)³s·cosθ` numerator term — the ℓ±1 coupling `2ias·cosθ/R` — is split off).
"""
function rt_CT_flat(R; M=1.0, a=0.0, s::Integer=0)
    pref = 1 / (R * (R + 1) * (R^2 + 1)^3)
    t1 = -a^2 * (R - 1) * (R + 1)^2 * (2M * (3R^6 + 5R^4 - 7R^2 - 1) - R^4 - 6R^2 - 1)
    t2 = 4M^2 * (R - 1) * R * (R + 1)^2 * (R^4 * (s + 2) + 2R^2 * (s + 3) + s)
    t3 = M * (R + 1) * (R^6 * (7s + 3) - 2R^5 * s + 13R^4 * (s + 1) - 4R^3 * (s + 2) +
              R^2 * (5s - 7) - 2R * s - s - 1)
    t4 = 2 * (R - 1) * R * (R^4 * s + R^3 + 2R^2 * (s + 1) + R + s)
    return pref * (t1 - 2 * (t2 + t3 + t4))
end

"""
    rt_CR(R; M=1, a=0, s=0)

Appendix coefficient C_R of [arXiv:1301.1591, Eq. (14)] (θ-independent as printed).
"""
function rt_CR(R; M=1.0, a=0.0, s::Integer=0)
    pref = -(R^2 - 1) / (2R^2 * (R^2 + 1)^3)
    t1 = a^2 * (2R^4 + 5R^2 + 1) * (R^2 - 1)^2
    t2 = M * (R^2 - 1) * (R^4 * (s + 3) + 2R^2 * (s + 4) + s + 1) +
         2R * (R^4 * (s + 1) + R^2 * (2s + 3) + s)
    return pref * (t1 + 2R * t2)
end

"""
    rt_CTT_flat(R; M=1, a=0)

θ-independent part of the appendix coefficient C_TT of [arXiv:1301.1591, Eq. (14)] (the
`a²(R+1)(R²+1)²cos2θ` numerator term — i.e. `a²cos2θ/(2R) = a²/(2R) − (a²/R)sin²θ`, the
ℓ,ℓ±1,ℓ±2 mass coupling — is split off). s-independent.
"""
function rt_CTT_flat(R; M=1.0, a=0.0)
    pref = 1 / (2R * (R + 1) * (R^2 + 1)^2)
    t1 = -a^2 * (R + 1) * (32M^2 * R^6 - (64M^2 + 32M + 1) * R^4 + (32M^2 + 32M + 6) * R^2 - 1)
    t2 = 16M^3 * (R - 1) * R^2 * (R + 1)^2 + 8M^2 * R * (R^3 - R^2 - 3R - 1) +
         M * (R^3 - 11R^2 - 5R - 1) - R * (R + 1)
    return pref * (t1 - 8R * t2)
end

"""
    rt_CRR(R; M=1, a=0)

Appendix coefficient C_RR of [arXiv:1301.1591, Eq. (14)]:
`−(R²−1)²(a²(R²−1)² + 4R(M(R²−1)+R))/(4R(R²+1)²)`. The second factor equals Δ(r(R))·(1−R²)²,
so C_RR vanishes at BOTH R₊ (Δ=0) and scri (R=1) — the outflow-boundary property
[arXiv:1301.1591, §2.3]. s-independent.
"""
rt_CRR(R; M=1.0, a=0.0) =
    -(R^2 - 1)^2 * (a^2 * (R^2 - 1)^2 + 4R * (M * (R^2 - 1) + R)) / (4R * (R^2 + 1)^2)

"""
    rt_CTR(R; M=1, a=0)

Appendix coefficient C_TR of [arXiv:1301.1591, Eq. (14)]. s-independent.
"""
rt_CTR(R; M=1.0, a=0.0) =
    -(2 / (R^2 + 1)^2) * (a^2 * (R^2 - 1)^2 * (2M * (R^2 - 1) - 1) +
                          2 * (4M^2 * R * (R^2 - 1)^2 + M * (3R^2 + 4R + 1) * (R - 1)^2 - 2R^2))

"""
    rt_CTphi(R; M=1, a=0)

Appendix coefficient C_Tφ of [arXiv:1301.1591, Eq. (14)]: `4a(1−2M(R²−1))/(R²+1)`
(enters the m-mode as `C̃_T = C_T + im·C_Tφ`, appendix tilde table).
"""
rt_CTphi(R; M=1.0, a=0.0) = 4a * (1 - 2M * (R^2 - 1)) / (R^2 + 1)

"""
    rt_CRphi(R; M=1, a=0)

Appendix coefficient C_Rφ of [arXiv:1301.1591, Eq. (14)]: `−a(R²−1)²/(R³+R)`
(enters the m-mode as `C̃_R = C_R + im·C_Rφ`, appendix tilde table).
"""
rt_CRphi(R; M=1.0, a=0.0) = -a * (R^2 - 1)^2 / (R^3 + R)

"""
    rt_Cphi_re(R; a=0)

Real (θ-independent) part of the appendix coefficient C_φ of [arXiv:1301.1591, Eq. (14)]:
`a(1−R²)/R²` (enters the m-mode as `C̃₀ ⊃ im·a(1−R²)/R²`, appendix tilde table; the imaginary
part `−2s·cotθcscθ/R` is consumed by the angular operator, header).
"""
rt_Cphi_re(R; a=0.0) = a * (1 - R^2) / R^2

# ---- projected 1+1 blocks (DERIVED-RT, header) ---------------------------------------------------

"""
    rt_ATT(R, s, ℓs, m; M=1, a=0) -> Matrix (real, symmetric)

Projected ∂_TT block `𝒜_TT = (rt_CTT_flat + a²/(2R))·I − (a²/R)·𝒮` over the mode list `ℓs`
(header DERIVED-RT; 𝒮 = sin²θ matrix [arXiv:1702.04204, Eqs. (33),(35)]). PENTADIAGONAL for
s·m ≠ 0; parity-split tridiagonal-in-ℓ±2 otherwise. Symmetric (C₊^ℓ = C₋^{ℓ+1},
C₊₊^ℓ = C₋₋^{ℓ+2} [arXiv:1702.04204, Eq. (34)]; dossier §3.3).
"""
function rt_ATT(R, s::Integer, ℓs::Vector{Int}, m::Integer; M=1.0, a=0.0)
    n = length(ℓs)
    A = -(a^2 / R) .* _spin_sin2_matrix(s, ℓs, m)
    d = rt_CTT_flat(R; M=M, a=a) + a^2 / (2R)
    for i in 1:n
        A[i, i] += d
    end
    return A
end

"""
    rt_AT(R, s, ℓs, m; M=1, a=0) -> Matrix (complex)

Projected ∂_T block `𝒜_T = (rt_CT_flat + im·C_Tφ)·I + (2ias/R)·𝒦` over the mode list `ℓs`
(header DERIVED-RT; 𝒦 = cosθ matrix [arXiv:1702.04204, Eqs. (32),(34)]). The imaginary
tridiagonal `(2ias/R)𝒦` band is the ℓ±1 friction coupling β± (dossier §3.1 item 5) — nonzero
whenever a·s ≠ 0 (for EVERY m), which is what breaks the even/odd-ℓ parity split.
"""
function rt_AT(R, s::Integer, ℓs::Vector{Int}, m::Integer; M=1.0, a=0.0)
    n = length(ℓs)
    A = Complex{typeof(float(R))}.((2im * a * s / R) .* _spin_cos_matrix(s, ℓs, m))
    d = rt_CT_flat(R; M=M, a=a, s=s) + im * m * rt_CTphi(R; M=M, a=a)
    for i in 1:n
        A[i, i] += d
    end
    return A
end

"""
    rt_C0diag(R, s, l, m; M=1, a=0) -> Complex

Projected zeroth-order diagonal `C̄₀^{(ℓ)} = rt_C0_flat + im·a(1−R²)/R² + (ℓ−s)(ℓ+s+1)/R`
(header DERIVED-RT; the eigenvalue term from [arXiv:1702.04204, Eq. (25)] — reduces to
ℓ(ℓ+1) at s=0).
"""
rt_C0diag(R, s::Integer, l::Integer, m::Integer; M=1.0, a=0.0) =
    rt_C0_flat(R; M=M, a=a, s=s) + im * m * rt_Cphi_re(R; a=a) + (l - s) * (l + s + 1) / R

"""
    rt_CtildeR(R, m; M=1, a=0, s=0) -> Complex

m-mode radial-advection scalar `C̃_R = C_R + im·C_Rφ` [arXiv:1301.1591, appendix tilde table].
"""
rt_CtildeR(R, m::Integer; M=1.0, a=0.0, s::Integer=0) =
    rt_CR(R; M=M, a=a, s=s) + im * m * rt_CRphi(R; M=M, a=a)

# ---- chart map to Boyer–Lindquist (a=0; Bardeen–Press gate D) ------------------------------------

# rt_F_of_R (F = t − T, the a=0 RT↔BL height composition, Eqs. (3),(5),(6)) is a PURE chart map —
# RE-HOMED to Spacetime.RTChart (F5 Part B); re-imported by the TimeDomain module header (same
# object, ===). Used by gate D (Bardeen–Press push-forward) and the sourced-sector jump amplitudes;
# both call sites now resolve to Spacetime.RTChart.rt_F_of_R. Def absent here.

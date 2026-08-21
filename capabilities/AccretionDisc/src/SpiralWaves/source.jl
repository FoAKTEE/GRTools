# SpiralWaves/source.jl — the relativistic gravitational potential of the small object, arXiv:2608.17003
# Sec. II.D (main.tex:653-711) + appendix B (appendix.tex:27-53): the leading-order singular field
#
#     h_μν(x) = (2m•/s0) [g_μν(w) + 2 v_μ v_ν]                                             (52)
#
# of the circular equatorial EMRI worldline w^μ(t) = (t, a•, Ω•t, 0) (53) with four-velocity
# v^μ = dt/dτ (1, 0, Ω•, 0) (54)-(55), the coordinate-distance approximation s0² (58) in its printed
# equatorial closed form (B1) (a function of r and Ψ = φ − Ω•t only, (B2)), the Plummer softening
# ε_soft = a• h0 b (59) (`softening_length` of background.jl), and the cosine Fourier projection (60)
#
#     ĥ_m(r) = (1/2π) ∫₀^{2π} h(r, Ψ) cos(mΨ) dΨ
#
# by the PERIODIC TRAPEZOID RULE (spectrally accurate for the smooth softened integrand; the node
# count is derived from the exact nearest complex zero of s0², see `strip_halfwidth`). Because g_μν(w)
# and v_μ are evaluated ON THE WORLDLINE, every nonzero component is a CONSTANT tensor K_μν =
# 2(g_μν(w) + 2v_μv_ν) times the one scalar kernel 1/s0(r, Ψ): ĥ_m,μν(r) = m• K_μν I_m(r) with
# I_m = (1/2π)∫cos(mΨ)/s0 dΨ, and ∂_r ĥ_m,μν = m• K_μν I_m'(r) with I_m' = −(1/4π)∫cos(mΨ) ∂_r s0²/s0³ dΨ
# (differentiation under the integral; s0² is a quadratic polynomial in r at fixed Ψ). Nonzero
# components: h_tt, h_tφ, h_φφ, h_rr, h_zz (h_tr = h_rφ = h_φz = 0 — the paper's `circmetriczero`).
#
# CONVENTIONS (EGD::deriv-singular-field CAS certificate, 26/26 PASS): covariant modes ĥ_μν are the
# PRIMARY output — exactly what the master system consumes
# (`master_inputs(...; hmodes)`, property names htt, htt_r, htphi, htphi_r, hphiphi, hphiphi_r, hrr, hzz);
# (60) is printed for h^{μν}, and raising with the Ψ-independent background metric g^{μν}(r) at the FIELD
# POINT commutes with the projection, so `raise_modes` gives ĥ^{μν}_m from ĥ_m,μν exactly. Newtonian
# limit (M/a• → 0): h_tt → +2m•/s0, i.e. the linearised g_tt = −(1 + 2Φ) gives Φ = −h_tt/2 = −m•/s0 < 0,
# the attractive Plummer-softened potential with the sign used in (50f,g) — no sign flip anywhere.
# Per-unit-mass (`mEMRI = 1` default): every output is ∝ m• = q M. Type-generic `{T}`; BigFloat
# precision is inherited (the node count grows with −log eps(T)). Kerr geometry is reached ONLY via
# `Spacetime.kerr_equatorial_potentials(_dr)` / `kerr_circular_dtdtau`. Additive (EGD::sw-singular-source).

using Spacetime: kerr_circular_dtdtau

"""
    worldline_velocity(bg::DiscBackground) -> (ut, uphi, vt, vphi)

Four-velocity of the circular equatorial EMRI worldline (53)-(55), M = 1: contravariant
`ut` = dt/dτ = (a⋆^{3/2}+χ)/(a⋆^{3/4}√(√a⋆(a⋆−3)+2χ)) (`Spacetime.kerr_circular_dtdtau(χ, a•)`, the
named accessor of the existing closed form), `uphi` = Ω• dt/dτ (Ω• from (13), `Omega_emri`), and the
COVARIANT components v_t = g_tt u^t + g_tφ u^φ (= −E), v_φ = g_tφ u^t + g_φφ u^φ (= +L) lowered with
the Page–Thorne metric at r = a• (g_tt = −e^{2ν} + e^{2ψ}ω̃², g_tφ = −e^{2ψ}ω̃, g_φφ = e^{2ψ}). χ < 0
is the retrograde orbit (motion in +φ, spin mirrored), as everywhere in SpiralWaves.
"""
function worldline_velocity(bg::DiscBackground{T}) where {T}
    a = bg.aEMRI; χ = bg.χ
    pot = kerr_equatorial_potentials(χ, a)
    e2ν = exp(2 * pot.ν); e2ψ = exp(2 * pot.ψ); ω̃ = pot.ω̃
    gtt = -e2ν + e2ψ * ω̃^2
    gtφ = -e2ψ * ω̃
    gφφ = e2ψ
    ut = convert(T, kerr_circular_dtdtau(χ, a))
    uφ = Omega_emri(bg) * ut
    (ut = ut, uphi = uφ, vt = gtt * ut + gtφ * uφ, vphi = gtφ * ut + gφφ * uφ)
end

"""
    SingularFieldCoeffs{T}

Worldline data of the singular field (52) and of s0² (B1), all independent of the field point
(computed once per background by [`singular_field_coeffs`](@ref)):

  - `A0, A2, B1, C2`: s0²(r, Ψ) = A0 + A2 r² + B1 r cosΨ + C2 r² cos2Ψ — the printed (B1) regrouped
    by its r-dependence (A0 = ε² + a⁴/Δ_a, B1 = −2a³/Δ_a, A2 = (dt/dτ)²Q²/(2a⁴) + ½[(a+2)χ²/a³ + a²/Δ_a + 1],
    C2 = ½[−(a+2)χ²/a³ + a²/Δ_a − 1] − (dt/dτ)²Q²/(2a⁴), Δ_a = (a−2)a + χ², Q = Ω•(a³+(a+2)χ²) − 2χ);
  - `Pr, Pt`: the SAME (B1) in its cancellation-free radial/tangential form (exact trigonometric
    rearrangement, the structure of (58): P_rr = g_rr(a) radial weight, P_φφ/a² tangential weight, no cross
    term — CAS certificate of EGD::deriv-singular-field),

        s0² = ε² + P_r (r cosΨ − a)² + P_t r² sin²Ψ,    P_r = a²/Δ_a = −B1/(2a),   P_t = A2 − C2,

    which is what every kernel evaluates (with r cosΨ − a = (r − a) − 2r sin²(Ψ/2)): the regrouped
    A0 + … form is a difference of O(a⁴/Δ_a) terms and loses ~(a/ε_soft)² eps relative accuracy near the
    worldline (1e-13 for a• = 20, b h0 = 0.015), the radial/tangential form is a sum of non-negative terms;
  - `Ktt, Ktphi, Kphiphi, Krr, Kzz`: K_μν = 2(g_μν(w) + 2v_μv_ν), the constant tensor of (52) per unit
    mass (h_μν = m• K_μν/s0);
  - `ut, uphi, vt, vphi`: the worldline four-velocity (54)-(55);
  - `εsoft`: the Plummer softening length (59); `a`, `χ` for reference.
"""
struct SingularFieldCoeffs{T<:Real}
    a::T; χ::T; εsoft::T
    A0::T; A2::T; B1::T; C2::T
    Pr::T; Pt::T
    Ktt::T; Ktphi::T; Kphiphi::T; Krr::T; Kzz::T
    ut::T; uphi::T; vt::T; vphi::T
end

"""
    singular_field_coeffs(bg::DiscBackground) -> SingularFieldCoeffs

Evaluate the worldline data of [`SingularFieldCoeffs`](@ref): the (B1) coefficients TRANSCRIBED
from appendix.tex:31-44 (M = 1, a⋆ = a•, dt/dτ from (55), Ω• from (13)), and K_μν = 2(g_μν(a•) + 2v_μv_ν)
with the Page–Thorne metric at r = a• and the covariant four-velocity of [`worldline_velocity`](@ref).
"""
function singular_field_coeffs(bg::DiscBackground{T}) where {T}
    a = bg.aEMRI; χ = bg.χ
    ε = softening_length(bg)
    w = worldline_velocity(bg)
    τ = w.ut
    ΩE = Omega_emri(bg)
    Δa = (a - 2) * a + χ^2
    Q = ΩE * (a^3 + (a + 2) * χ^2) - 2 * χ
    τQ = τ^2 * Q^2 / (2 * a^4)
    A0 = ε^2 + a^4 / Δa
    B1 = -2 * a^3 / Δa
    A2 = τQ + ((a + 2) * χ^2 / a^3 + a^2 / Δa + 1) / 2
    C2 = (-(a + 2) * χ^2 / a^3 + a^2 / Δa - 1) / 2 - τQ
    pot = kerr_equatorial_potentials(χ, a)
    e2ν = exp(2 * pot.ν); e2ψ = exp(2 * pot.ψ); ω̃ = pot.ω̃
    gtt = -e2ν + e2ψ * ω̃^2
    gtφ = -e2ψ * ω̃
    gφφ = e2ψ
    grr = exp(2 * pot.μ)
    SingularFieldCoeffs{T}(a, χ, ε, A0, A2, B1, C2, a^2 / Δa, A2 - C2,
                           2 * (gtt + 2 * w.vt^2), 2 * (gtφ + 2 * w.vt * w.vphi), 2 * (gφφ + 2 * w.vphi^2),
                           2 * grr, T(2), w.ut, w.uphi, w.vt, w.vphi)
end

# the allocation-free kernels: s0² and ∂_r s0² in the cancellation-free radial/tangential form of (B1) at
# (r, sh² = sin²(Ψ/2)); cosΨ = 1 − 2sh², sin²Ψ = 4sh²(1 − sh²), r cosΨ − a = (r − a) − 2r sh²
@inline function _s0sq(K::SingularFieldCoeffs, r, sh2)
    d = (r - K.a) - 2 * r * sh2
    K.εsoft^2 + K.Pr * d * d + K.Pt * r * r * (4 * sh2 * (1 - sh2))
end
@inline function _s0sq_dr(K::SingularFieldCoeffs, r, sh2)
    d = (r - K.a) - 2 * r * sh2
    2 * K.Pr * d * (1 - 2 * sh2) + 2 * K.Pt * r * (4 * sh2 * (1 - sh2))
end
# the printed regrouping A0 + A2 r² + B1 r cosΨ + C2 r² cos2Ψ (for the form = :printed cross-check only)
@inline _s0sq_printed(K::SingularFieldCoeffs, r, c1, c2) = K.A0 + r * (K.B1 * c1 + r * (K.A2 + K.C2 * c2))
@inline _s0sq_printed_dr(K::SingularFieldCoeffs, r, c1, c2) = K.B1 * c1 + 2 * r * (K.A2 + K.C2 * c2)

"""
    s0_squared(bg::DiscBackground, r, Ψ; form = :stable) -> s0²

The equatorial coordinate-distance approximation s0²(r, Ψ) of (58) in its printed closed form (B1)
(appendix.tex:31-44), Ψ = φ − Ω•t (B2), including the softening +ε_soft² (59); M = 1, a⋆ = a•. Even and
2π-periodic in Ψ; zero on the unsoftened worldline (r = a•, Ψ = 0). `form = :stable` (default) evaluates
the cancellation-free radial/tangential rearrangement ε² + P_r(r cosΨ − a)² + P_t r² sin²Ψ used by every
kernel; `form = :printed` evaluates the literal regrouped (B1) A0 + A2r² + B1 r cosΨ + C2 r² cos2Ψ (identical
algebra; loses ~(a•/ε_soft)² eps near the worldline — see [`SingularFieldCoeffs`](@ref)). Type-generic.
"""
function s0_squared(bg::DiscBackground{T}, r, Ψ; form::Symbol = :stable) where {T}
    S = promote_type(T, typeof(r), typeof(Ψ))
    K = singular_field_coeffs(bg)
    rr = convert(S, r); ψ = convert(S, Ψ)
    if form === :stable
        return _s0sq(K, rr, sin(ψ / 2)^2)
    elseif form === :printed
        return _s0sq_printed(K, rr, cos(ψ), cos(2 * ψ))
    end
    throw(ArgumentError("s0_squared: form must be :stable or :printed, got $(repr(form))"))
end

"Analytic ∂_r of [`s0_squared`](@ref) (same `form` keyword): 2P_r(r cosΨ − a)cosΨ + 2P_t r sin²Ψ ≡ B1 cosΨ + 2r(A2 + C2 cos2Ψ)."
function s0_squared_dr(bg::DiscBackground{T}, r, Ψ; form::Symbol = :stable) where {T}
    S = promote_type(T, typeof(r), typeof(Ψ))
    K = singular_field_coeffs(bg)
    rr = convert(S, r); ψ = convert(S, Ψ)
    if form === :stable
        return _s0sq_dr(K, rr, sin(ψ / 2)^2)
    elseif form === :printed
        return _s0sq_printed_dr(K, rr, cos(ψ), cos(2 * ψ))
    end
    throw(ArgumentError("s0_squared_dr: form must be :stable or :printed, got $(repr(form))"))
end

"""
    singular_field(bg::DiscBackground, r, Ψ; mEMRI = 1) -> (htt, htphi, hphiphi, hrr, hzz)

The pointwise singular field (52) in the equatorial plane, covariant components at (r, Ψ = φ − Ω•t):
h_μν = m• K_μν / s0(r, Ψ) with K_μν of [`singular_field_coeffs`](@ref) (the vanishing h_tr, h_rφ, h_φz,
h_tz, h_rz are not returned). `mEMRI` = m• (default 1: per unit mass).
"""
function singular_field(bg::DiscBackground{T}, r, Ψ; mEMRI = 1) where {T}
    S = promote_type(T, typeof(r), typeof(Ψ), typeof(mEMRI))
    K = singular_field_coeffs(bg)
    rr = convert(S, r); ψ = convert(S, Ψ)
    f = convert(S, mEMRI) / sqrt(_s0sq(K, rr, sin(ψ / 2)^2))
    (htt = K.Ktt * f, htphi = K.Ktphi * f, hphiphi = K.Kphiphi * f, hrr = K.Krr * f, hzz = K.Kzz * f)
end

# ---------------------------------------------------------------------------------------------
# Fourier projection (60): periodic trapezoid rule, node count from the analytic strip of 1/s0
# ---------------------------------------------------------------------------------------------

"""
    strip_halfwidth(bg::DiscBackground, r) -> δ

Half-width δ of the strip |Im Ψ| < δ in which the integrand 1/s0(r, Ψ) of (60) is analytic: s0² =
A + B w + C(2w² − 1) with w = cosΨ, so its complex zeros are Ψ = ±arccos(w_±) + 2πk with w_± the roots of
2C w² + B w + (A − C) = 0, and δ = min |Im arccos(w_±)|. The periodic trapezoid rule then converges like
e^{−Nδ} (aliasing of the m-th cosine coefficient ≈ e^{−(N−2m)δ} relative). At r = a• and for small
softening δ ≈ ε_soft/(a•√P_tan) ≈ h0 b — the softening b sets the resolution: b → 0 needs N → ∞ (the
unsoftened kernel has a real singularity). Away from a• the strip widens (δ ~ |ln(r/a•)|). Evaluated
directly from the (B1) coefficients (`singular_field_coeffs`); type-generic.
"""
function strip_halfwidth(bg::DiscBackground{T}, r) where {T}
    K = singular_field_coeffs(bg)
    _strip_halfwidth(K, convert(promote_type(T, typeof(r)), r))
end

function _strip_halfwidth(K::SingularFieldCoeffs, r::S) where {S<:Real}
    A = K.A0 + K.A2 * r^2
    B = K.B1 * r
    C = K.C2 * r^2
    if iszero(C)
        w = Complex(-A / B, zero(S))
        return abs(imag(acos(w)))
    end
    disc = sqrt(Complex(B^2 - 8 * C * (A - C), zero(S)))
    w1 = (-B + disc) / (4 * C)
    w2 = (-B - disc) / (4 * C)
    min(abs(imag(acos(w1))), abs(imag(acos(w2))))
end

"""
    trapezoid_nodes(bg::DiscBackground, m, r; T = float type) -> N

Default (even) node count of the periodic trapezoid projection (60) at radius `r` for azimuthal
number `m`: N = max(256, 8m + 64, 2m + ⌈(−ln eps(T) + ln(1/δ) + 8)/δ⌉) rounded up to even, with δ =
[`strip_halfwidth`](@ref)`(bg, r)`. The third term bounds the aliasing error of the m-th coefficient
(≈ e^{−(N−2m)δ}, relative to the coefficient itself, times the strip growth ∝ 1/δ) by ~ e^{−8} eps(T);
the first two are the floor of the task specification (and dominate away from r ≈ a•). Grows like
1/(h0 b) near the worldline and like −ln eps(T) in BigFloat. Pass `N` explicitly to
[`singular_field_mode`](@ref) to override.
"""
function trapezoid_nodes(bg::DiscBackground{T}, m::Integer, r) where {T}
    S = promote_type(T, typeof(r))
    _trapezoid_nodes(singular_field_coeffs(bg), m, convert(S, r))
end

function _trapezoid_nodes(K::SingularFieldCoeffs, m::Integer, r::S) where {S<:Real}
    δ = _strip_halfwidth(K, r)
    δ > 0 || throw(DomainError(r, "trapezoid_nodes: 1/s0 is singular on the real Ψ axis (b = 0 on the worldline?)"))
    need = 2 * m + ceil(Int, (-log(eps(S)) + log(1 / δ) + 8) / δ)
    N = max(256, 8 * m + 64, need)
    N + (N & 1)                                                    # even
end

# Neumaier (Kahan–Babuška) compensated accumulation step: returns the updated (sum, compensation).
# Keeps the roundoff floor of the N-term trapezoid sums at O(eps) instead of O(√N eps) (N ~ 10³).
@inline function _kbn_add(s, c, x)
    t = s + x
    c += abs(s) ≥ abs(x) ? (s - t) + x : (x - t) + s
    (t, c)
end

# the cosine-projected kernel integrals I_m = (1/2π)∫cos(mΨ)/s0 dΨ and I_m' = −(1/4π)∫cos(mΨ) ∂_r s0²/s0³ dΨ
# by the N-node periodic trapezoid rule on [0, 2π), using the evenness of s0² in Ψ (nodes k and N−k coincide):
# (1/N)[f(0) + f(π) + 2 Σ_{k=1}^{N/2−1} f(2πk/N)]. Allocation-free; every node is evaluated in the working
# type S (cos of the node angle in S for BigFloat fidelity); compensated summation.
function _kernel_integrals(K::SingularFieldCoeffs, m::Integer, r::S, N::Integer) where {S<:Real}
    iseven(N) && N ≥ 2 || throw(ArgumentError("trapezoid node count must be even and ≥ 2, got N = $N"))
    twoπ = 2 * S(π)
    # endpoints Ψ = 0 (sin²(Ψ/2) = 0) and Ψ = π (sin²(Ψ/2) = 1)
    s0 = sqrt(_s0sq(K, r, zero(S)))
    acc = 1 / s0; cacc = zero(S)
    accd = _s0sq_dr(K, r, zero(S)) / s0^3; caccd = zero(S)
    s0 = sqrt(_s0sq(K, r, one(S)))
    cm = isodd(m) ? -one(S) : one(S)
    acc, cacc = _kbn_add(acc, cacc, cm / s0)
    accd, caccd = _kbn_add(accd, caccd, cm * _s0sq_dr(K, r, one(S)) / s0^3)
    half = N ÷ 2
    @inbounds for k in 1:half-1
        sh = sin(S(π) * k / N)                                       # sin(Ψ/2), Ψ = 2πk/N
        sh2 = sh * sh
        cmk = cos(twoπ * ((m * k) % N) / N)                          # exact reduction of mΨ mod 2π
        s0 = sqrt(_s0sq(K, r, sh2))
        inv3 = 1 / (s0 * s0 * s0)
        acc, cacc = _kbn_add(acc, cacc, 2 * cmk / s0)
        accd, caccd = _kbn_add(accd, caccd, 2 * cmk * _s0sq_dr(K, r, sh2) * inv3)
    end
    ((acc + cacc) / N, -(accd + caccd) / (2 * N))
end

"""
    singular_field_mode(bg::DiscBackground, m, r; N = trapezoid_nodes(bg, m, r), mEMRI = 1) -> NamedTuple

The cosine Fourier modes (60) of every nonzero covariant component of the singular field (52) at radius
`r` and their analytic r-derivatives, as the NamedTuple

    (htt, htt_r, htphi, htphi_r, hphiphi, hphiphi_r, hrr, hrr_r, hzz, hzz_r, I, I_r, N)

with ĥ_m,μν = m• K_μν I_m(r), ∂_r ĥ_m,μν = m• K_μν I_m'(r), I_m = (1/2π)∫₀^{2π} cos(mΨ)/s0(r,Ψ) dΨ the
per-unit-mass kernel mode (so that the Newtonian-limit potential mode is Φ̂_m = −ĥ_tt,m/2 → −m• I_m, the
softened Plummer potential), I_m' by differentiation of s0^{-1} under the integral, and `N` the node
count of the periodic trapezoid rule used (default [`trapezoid_nodes`](@ref); pass `N` to override).
The property names `htt … hzz` are exactly those consumed by `master_inputs(...; hmodes)` (COVARIANT,
the paper's convention, ERRATUM E-M2); use [`raise_modes`](@ref) for ĥ^{μν}. `mEMRI` = m• (default 1 → results ∝ q).
Allocation-free, type-generic (BigFloat inherits); the result is a deterministic function of
(bg, m, r, N) — identical bits from [`source_modes`](@ref) on any thread.
"""
function singular_field_mode(bg::DiscBackground{T}, m::Integer, r; N::Union{Nothing,Integer} = nothing, mEMRI = 1) where {T}
    S = promote_type(T, typeof(r), typeof(mEMRI))
    K = singular_field_coeffs(bg)
    _singular_field_mode(K, m, convert(S, r), N, convert(S, mEMRI))
end

function _singular_field_mode(K::SingularFieldCoeffs, m::Integer, r::S, N::Union{Nothing,Integer}, mE::S) where {S<:Real}
    m ≥ 0 || throw(DomainError(m, "singular_field_mode: need m ≥ 0 (the cosine modes are even in m)"))
    NN = N === nothing ? _trapezoid_nodes(K, m, r) : Int(N)
    I, Ir = _kernel_integrals(K, m, r, NN)
    f = mE * I; fr = mE * Ir
    (htt = K.Ktt * f, htt_r = K.Ktt * fr,
     htphi = K.Ktphi * f, htphi_r = K.Ktphi * fr,
     hphiphi = K.Kphiphi * f, hphiphi_r = K.Kphiphi * fr,
     hrr = K.Krr * f, hrr_r = K.Krr * fr,
     hzz = K.Kzz * f, hzz_r = K.Kzz * fr,
     I = I, I_r = Ir, N = NN)
end

"""
    source_modes(bg::DiscBackground, m, rgrid; threads = false, N = nothing, mEMRI = 1) -> Vector{NamedTuple}

[`singular_field_mode`](@ref) at every radius of `rgrid` (the default N is chosen per radius), returned
as a `Vector` of the same NamedTuples (element `i` feeds `master_inputs(bg, m, rgrid[i]; hmodes = out[i])`).
`threads = true` distributes the radii over `Threads.@threads` (default scheduler, safe to nest inside an
outer threaded m-loop); every point is an independent deterministic function of (bg, m, r, N), so the
threaded result is BITWISE identical to the serial one. The worldline coefficients are computed once; the
per-point work is allocation-free. See [`source_modes!`](@ref) for the in-place variant.
"""
function source_modes(bg::DiscBackground{T}, m::Integer, rgrid::AbstractVector; threads::Bool = false,
                      N::Union{Nothing,Integer} = nothing, mEMRI = 1) where {T}
    S = promote_type(T, eltype(rgrid), typeof(mEMRI))
    K = singular_field_coeffs(bg)
    # element type = the concrete NamedTuple of one (cheap, N = 2) evaluation
    proto = _singular_field_mode(K, m, convert(S, first(rgrid)), 2, convert(S, mEMRI))
    out = Vector{typeof(proto)}(undef, length(rgrid))
    source_modes!(out, bg, m, rgrid; threads = threads, N = N, mEMRI = mEMRI)
end

"""
    source_modes!(out, bg::DiscBackground, m, rgrid; threads = false, N = nothing, mEMRI = 1) -> out

In-place [`source_modes`](@ref): fills `out[i] = singular_field_mode(bg, m, rgrid[i]; N, mEMRI)` (no
allocation beyond the one-off worldline coefficients).
"""
function source_modes!(out::AbstractVector, bg::DiscBackground{T}, m::Integer, rgrid::AbstractVector;
                       threads::Bool = false, N::Union{Nothing,Integer} = nothing, mEMRI = 1) where {T}
    length(out) == length(rgrid) || throw(DimensionMismatch("source_modes!: length(out) = $(length(out)) ≠ length(rgrid) = $(length(rgrid))"))
    S = promote_type(T, eltype(rgrid), typeof(mEMRI))
    K = singular_field_coeffs(bg)
    mE = convert(S, mEMRI)
    if threads
        Threads.@threads for i in eachindex(rgrid)
            @inbounds out[i] = _singular_field_mode(K, m, convert(S, rgrid[i]), N, mE)
        end
    else
        @inbounds for i in eachindex(rgrid)
            out[i] = _singular_field_mode(K, m, convert(S, rgrid[i]), N, mE)
        end
    end
    out
end

# ---------------------------------------------------------------------------------------------
# contravariant modes ĥ^{μν} = g^{μα}(r) g^{νβ}(r) ĥ_αβ (the index position of the printed (60))
# ---------------------------------------------------------------------------------------------

"""
    raise_modes(bg::DiscBackground, r, hmodes) -> NamedTuple

Contravariant modes ĥ^{μν}_m(r) = g^{μα}(r) g^{νβ}(r) ĥ_m,αβ(r) (the index position in which (60) is
printed) and their r-derivatives from the covariant `hmodes` (any object with properties
`htt, htt_r, htphi, htphi_r, hphiphi, hphiphi_r, hrr, hrr_r, hzz, hzz_r`, e.g. [`singular_field_mode`](@ref)),
raised with the Page–Thorne background at the FIELD POINT r (Ψ-independent, so raising commutes with the
cosine projection): g^{tt} = −e^{−2ν}, g^{tφ} = −ω̃e^{−2ν}, g^{φφ} = e^{−2ψ} − ω̃²e^{−2ν}, g^{rr} = e^{−2μ},
g^{zz} = 1, and d/dr of each from `Spacetime.kerr_equatorial_potentials_dr` (product rule). Returns
`(htt, htt_r, htphi, htphi_r, hphiphi, hphiphi_r, hrr, hrr_r, hzz, hzz_r)` — the SAME property names, now
meaning ĥ^{tt}, ĥ^{tφ}, … (do NOT feed these to `master_inputs`, which takes the covariant modes).
"""
function raise_modes(bg::DiscBackground{T}, r, hmodes) where {T}
    S = promote_type(T, typeof(r), typeof(hmodes.htt))
    rr = convert(S, r); χ = convert(S, bg.χ)
    pot = kerr_equatorial_potentials(χ, rr)
    d = kerr_equatorial_potentials_dr(χ, rr)
    em2ν = exp(-2 * pot.ν); em2ψ = exp(-2 * pot.ψ); em2μ = exp(-2 * pot.μ); ω̃ = pot.ω̃
    gtt = -em2ν
    gtφ = -ω̃ * em2ν
    gφφ = em2ψ - ω̃^2 * em2ν
    grr = em2μ
    dgtt = 2 * d.dν * em2ν
    dgtφ = -d.dω̃ * em2ν + 2 * d.dν * ω̃ * em2ν
    dgφφ = -2 * d.dψ * em2ψ - 2 * ω̃ * d.dω̃ * em2ν + 2 * d.dν * ω̃^2 * em2ν
    dgrr = -2 * d.dμ * em2μ
    htt = convert(S, hmodes.htt); htt_r = convert(S, hmodes.htt_r)
    htφ = convert(S, hmodes.htphi); htφ_r = convert(S, hmodes.htphi_r)
    hφφ = convert(S, hmodes.hphiphi); hφφ_r = convert(S, hmodes.hphiphi_r)
    hrr = convert(S, hmodes.hrr); hrr_r = convert(S, hmodes.hrr_r)
    hzz = convert(S, hmodes.hzz); hzz_r = convert(S, hmodes.hzz_r)
    # H^{ab} = G^{ac} H_cd G^{db} on the (t, φ) block
    Htt = gtt * gtt * htt + 2 * gtt * gtφ * htφ + gtφ * gtφ * hφφ
    Htφ = gtt * gtφ * htt + (gtt * gφφ + gtφ * gtφ) * htφ + gtφ * gφφ * hφφ
    Hφφ = gtφ * gtφ * htt + 2 * gtφ * gφφ * htφ + gφφ * gφφ * hφφ
    Htt_r = 2 * gtt * dgtt * htt + gtt * gtt * htt_r +
            2 * (dgtt * gtφ + gtt * dgtφ) * htφ + 2 * gtt * gtφ * htφ_r +
            2 * gtφ * dgtφ * hφφ + gtφ * gtφ * hφφ_r
    Htφ_r = (dgtt * gtφ + gtt * dgtφ) * htt + gtt * gtφ * htt_r +
            (dgtt * gφφ + gtt * dgφφ + 2 * gtφ * dgtφ) * htφ + (gtt * gφφ + gtφ * gtφ) * htφ_r +
            (dgtφ * gφφ + gtφ * dgφφ) * hφφ + gtφ * gφφ * hφφ_r
    Hφφ_r = 2 * gtφ * dgtφ * htt + gtφ * gtφ * htt_r +
            2 * (dgtφ * gφφ + gtφ * dgφφ) * htφ + 2 * gtφ * gφφ * htφ_r +
            2 * gφφ * dgφφ * hφφ + gφφ * gφφ * hφφ_r
    (htt = Htt, htt_r = Htt_r, htphi = Htφ, htphi_r = Htφ_r, hphiphi = Hφφ, hphiphi_r = Hφφ_r,
     hrr = grr * grr * hrr, hrr_r = 2 * grr * dgrr * hrr + grr * grr * hrr_r,
     hzz = hzz, hzz_r = hzz_r)
end

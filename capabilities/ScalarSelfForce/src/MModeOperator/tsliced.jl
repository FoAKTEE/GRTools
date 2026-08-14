# MModeOperator/tsliced.jl — the m-mode elliptic operator COEFFICIENT layer in first-order
# flux form (ESSF node `elliptic-operator-tsliced`, W4-C; arXiv:2606.04998 Sec. II.B.3-4,
# eq:flux main.tex:366-368, eq:first-order main.tex:370-373, β_m/γ^i_m main.tex:375-383,
# z² equatorial variant main.tex:389-394). The PDE in 2D coordinates q^i = (r, z) is
#
#   △_m Ψ_m = −∂_i 𝓕^i_m + β_m Ψ_m + γ^i_m ∂_i Ψ_m ,
#   𝓕^i_m = (A^rr ∂_r Ψ_m, A^zz ∂_z Ψ_m),  A^rr = Δ/(r²+a²),  A^zz = (1−z²)/(r²+a²),
#
# in horizon-penetrating Boyer–Lindquist r (∂_{r_*} = Δ/(r²+a²) ∂_r, main.tex:354-358) with a
# CALLER-supplied boost (H, H′ ≡ dH/dr_*) evaluated at the caller's point, so the SAME
# coefficients serve every slicing: t-slicing is the (H, H′) = (0, 0) default (this node);
# the vtu piecewise H ∈ {−1, 0, +1} assembly and its junction conditions are a LATER node
# (vtu-slicing-jumps). The flux coefficients A^ii carry NO H-dependence — all slicing
# dependence lives in β_m and γ^r_m, pointwise algebraic in (H, H′) via the replacement rules
# ∂_{r_*} → ∂_{r_*} + imΩH, ∂²_{r_*} → ∂²_{r_*} + 2imΩH∂_{r_*} + imΩH′ − m²Ω²H²
# (main.tex:306-308).
#
# GENERAL-m (|m|) READING — the single deviation from the printed forms: the paper's
# m(m+1) (in β_m) and 2mz (in γ^z_m) originate from the angular reduction of the ansatz
# factor sin^{|m|}θ, which produces (1/sinθ)∂_θ(sinθ∂_θ) − m²/sin²θ acting on sin^{|m|}θΨ =
# sin^{|m|}θ·[∂_z(1−z²)∂_z − 2|m|z∂_z − |m|(|m|+1)]Ψ (the m²−|m|² singular term cancels);
# the printed forms are the m ≥ 0 branch. All other m-dependence is purely imaginary-odd
# (im·a, imΩH, imΩH′) or even (m²). Hence coefficient_{−m} = conj(coefficient_m) exactly —
# consistent with the ansatz layer's conjugate symmetry Ψ_{−m} = conj(Ψ_m) (ansatz.jl).
# Production uses m ≥ 0 (module convention), where this is identical to the printed text.
#
# Pure pointwise functions of (r, z) given (orbit::CircularOrbit, m, H, H′); type-generic
# (BigFloat precision inherited through promotion, never down-cast); complex-valued exactly
# where m makes them complex (β_m, γ^r_m); Δ is evaluated in the factored form
# (r − r₊)(r − r₋) off the orbit's cached horizons — exact zero AT r₊ and no cancellation
# error near it (the horizon-regularity property main.tex:386 is gate (b) of the tests).

"""
    flux_coeff_r(orbit, r) -> A^rr = Δ/(r² + a²)

Radial flux coefficient of the first-order form (eq:flux, main.tex:366-368):
`𝓕^r_m = A^rr ∂_r Ψ_m`. Δ is computed as `(r − r₊)(r − r₋)` from the orbit's cached
horizons — exactly zero AT the horizon, which is why no explicit horizon boundary
condition is ever needed (main.tex:385 + 399-400). Positive for r > r₊ and shares its
sign with [`flux_coeff_z`](@ref) (ellipticity, main.tex:263). Slicing-independent (no H).
"""
flux_coeff_r(orbit::CircularOrbit, r) =
    (r - orbit.rp) * (r - orbit.rm) / (r^2 + orbit.a^2)

"""
    flux_coeff_z(orbit, r, z) -> A^zz = (1 − z²)/(r² + a²)

Angular flux coefficient of the first-order form (eq:flux): `𝓕^z_m = A^zz ∂_z Ψ_m`.
Computed as `(1 − z)(1 + z)` (pole-accurate product form, exactly zero at z = ±1): the
sin^{|m|}θ factoring of the ansatz makes the polar boundaries pure regularity conditions
with no explicit boundary data (main.tex:385 + 400).
"""
flux_coeff_z(orbit::CircularOrbit, r, z) =
    (1 - z) * (1 + z) / (r^2 + orbit.a^2)

"""
    flux_coeff_z2(orbit, r, z) -> A^{z²} = 4z²·A^zz

Angular flux coefficient of the equatorial-symmetry variant (main.tex:389-391): with the
polar coordinate w = z² on the half-domain w ∈ [0, 1], `𝓕^{z²}_m = A^{z²} ∂_w Ψ_m` and
`𝓕^{z²}_m = 4z² 𝓕^z_m` at coefficient level (the printed relation, with the derivative
slot reinterpreted as ∂_w = ∂_z/(2z)). The argument is still z — callers on a w-grid pass
`z = √w`. Exactly zero at BOTH ends: the pole w = 1 (regularity as above) and the equator
w = 0, so the even sector needs no explicit boundary data either.
"""
flux_coeff_z2(orbit::CircularOrbit, r, z) = 4 * z^2 * flux_coeff_z(orbit, r, z)

"""
    beta_m(orbit, m, r, z, H = zero(r), Hprime = zero(r)) -> Complex

Zeroth-order coefficient β_m of the first-order flux form (eq:first-order,
main.tex:370-379):

    β_m = [ |m|(|m|+1) + (2M/r)(1 − a²/(Mr)) + (2iam/r)(1 + aΩH) ] / (r²+a²)
        + (r²+a²)/Δ · [ m²Ω²(H² − 1 + a²Δ(1−z²)/(r²+a²)²)
                        + (2am²Ω/(r²+a²))(2Mr/(r²+a²) + H) − imΩH′ ]

with `(H, Hprime) = (H(r_*), dH/dr_*)` supplied by the caller at this point (t-slicing =
the (0, 0) defaults; vtu slabs pass their constant H with H′ = 0). The |m| in the first
term is the general-m reading of the printed m(m+1) (see the file header): β_{−m} =
conj(β_m) exactly. Horizon regularity (main.tex:386): at H = −1 the 1/Δ pieces cancel —
H² − 1 = 0 exactly, and (2Mr/(r²+a²) + H)/Δ → −1/(r²+a²) via 2Mr₊ − r₊² − a² = −Δ(r₊) —
so β_m is finite at r → r₊ (with H′ = 0); the residual near-horizon roundoff of the
as-printed evaluation is O(eps/Δ), gated at 1e-3 in test (b). With H = 0 it genuinely
diverges ~1/Δ (the t-sliced form is not horizon-regular; the v domain is).
"""
function beta_m(orbit::CircularOrbit, m::Integer, r, z, H = zero(r), Hprime = zero(r))
    a, M, Om = orbit.a, orbit.M, orbit.Omega
    D = (r - orbit.rp) * (r - orbit.rm)
    w = r^2 + a^2
    am = abs(m)
    (am * (am + 1) + (2M / r) * (1 - a^2 / (M * r)) + (2im * a * m / r) * (1 + a * Om * H)) / w +
        (w / D) * (m^2 * Om^2 * (H^2 - 1 + a^2 * D * (1 - z) * (1 + z) / w^2) +
                   (2a * m^2 * Om / w) * (2M * r / w + H) - im * m * Om * Hprime)
end

"""
    gamma_r_m(orbit, m, r, H = zero(r)) -> Complex

Radial first-order coefficient γ^r_m (main.tex:380):

    γ^r_m = 2[ a²Δ/(r(r²+a²)²) − ima/(r²+a²) − imΩH ]

z-independent; the −2imΩH piece is the boost contribution of the replacement rules
(t-slicing default H = 0). Purely real for m = 0; γ^r_{−m} = conj(γ^r_m) exactly.
"""
function gamma_r_m(orbit::CircularOrbit, m::Integer, r, H = zero(r))
    a, Om = orbit.a, orbit.Omega
    D = (r - orbit.rp) * (r - orbit.rm)
    w = r^2 + a^2
    2 * (a^2 * D / (r * w^2) - im * m * a / w - im * m * Om * H)
end

"""
    gamma_z_m(orbit, m, r, z) -> Real

Angular first-order coefficient γ^z_m = 2|m|z/(r² + a²) (main.tex:381-383, |m| reading of
the printed 2mz — see the file header). Real, even in m, slicing-independent.
"""
gamma_z_m(orbit::CircularOrbit, m::Integer, r, z) = 2 * abs(m) * z / (r^2 + orbit.a^2)

"""
    gamma_z2_m(orbit, m, r, z) -> Real

Angular first-order coefficient of the equatorial-symmetry variant (main.tex:392):
`γ^{z²}_m = 2z γ^z_m + 2(1−z²)/(r²+a²)`, multiplying ∂_w Ψ_m in the w = z² coordinate
(argument still z; callers on a w-grid pass z = √w). Real, even in m, slicing-independent.
"""
gamma_z2_m(orbit::CircularOrbit, m::Integer, r, z) =
    2 * z * gamma_z_m(orbit, m, r, z) + 2 * (1 - z) * (1 + z) / (r^2 + orbit.a^2)

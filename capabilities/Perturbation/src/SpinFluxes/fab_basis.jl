# Shared fab^(i) basis + angular layer + Carter-K helpers of the FromCorrection family — faithful
# 1:1 of TeukolskySpinFluxes/src/from_correction.jl (helper half) + the _dfabidS/_dfabidω analytic
# ω-derivative bases (spin_mode_spherical_correction*.jl). The fab^(i) basis is SOURCE-EXTRACTED+
# TRANSLATED from TeukolskySpinFluxes.wl L48-85 and characterization-verified vs the BHPToolkit
# oracle (REF_FABI_BASIS). Consumed by from_correction/spin_mode_num/spin_mode_correction/spherical*.

# --- CarterConstantK flux helpers (.wl L1606-1653). dJzd* = ∂J_z/∂{K,En,Lz} along the polar motion
# (geodesic + s·spin), KCorrection = the linear-in-spin correction Kc1 to the Carter-Killing K. ---

# z1,z2 (polar turning points), kzsq=z2²/z1² shared by all three dJzd* (.wl L1607-1609).
@inline function _dJz_zroots(a, Eso, Lso, Kso)
    disc = sqrt(a^4 - 2 * a^2 * Kso + Kso^2 - 4 * a^3 * Eso * Lso + 4 * a * Eso * Kso * Lso + 4 * a^2 * Lso^2)
    base = a^2 - 2 * a^2 * Eso^2 + Kso + 2 * a * Eso * Lso
    den  = 2 * (a^2 - a^2 * Eso^2)
    z2 = sqrt((base - disc) / den)
    z1 = sqrt((base + disc) / den)
    return z1, z2, z2^2 / z1^2
end

"""`_dJzdK(a, s, En, Lz, Kc)` — ∂J_z/∂K = dJz0dKso + s·dJz1dKso (.wl `dJzdK`, L1606-1613)."""
function _dJzdK(a, s, Eso, Lso, Kso)
    z1, z2, kzsq = _dJz_zroots(a, Eso, Lso, Kso)
    EK = real(ellipK(kzsq)); EE = real(ellipE(kzsq))
    dJz0dKso = EK / (a * sqrt(1 - Eso^2) * pi * z1)
    dJz1dKso = -(1 / (a * sqrt(1 - Eso^2) * sqrt(Kso) * pi * z1)) * (Eso * EK +
        ((a^2 * Eso - Eso * Kso - a * Lso) * (a^2 * z1^2 * EE + (Kso - a^2 * z1^2) * EK)) /
            ((Kso - a^2 * z1^2) * (Kso - a^2 * z2^2)) -
        (Kso * (-(((1 - z1^2) * (Lso - a * Eso * (1 - z1^2)) * EE) / ((1 - kzsq) * z1^2 * (Kso - a^2 * z1^2))) +
                ((1 - z2^2) * (Lso - a * Eso * (1 - z2^2)) * (-(EE / (1 - kzsq)) + EK)) / (z2^2 * (Kso - a^2 * z2^2)))) /
            (a * (1 - Eso^2) * (-z1^2 + z2^2)))
    return dJz0dKso + s * dJz1dKso
end

"""`_dJzdEn(a, s, En, Lz, Kc)` — ∂J_z/∂En = dJz0dEso + s·dJz1dEso (.wl `dJzdEn`, L1615-1622)."""
function _dJzdEn(a, s, Eso, Lso, Kso)
    z1, z2, kzsq = _dJz_zroots(a, Eso, Lso, Kso)
    EK = real(ellipK(kzsq)); EE = real(ellipE(kzsq))
    dJz0dEso = 1 / (sqrt(1 - Eso^2) * pi * z1) * 2 *
               (-a * Eso * z1^2 * EE + (Lso - a * Eso * (1 - z1^2)) * EK)
    dJz1dEso = (2 * sqrt(Kso) * ((z1^2 + z2^2 * (1 + z1^4 + z1^2 * (-4 + z2^2))) * EE +
                (-z1 + z2) * (z1 + z2) * (1 + (-2 + z1^2) * z2^2) * EK)) /
               (a * (1 - Eso^2)^(3 / 2) * pi * z1 * z2^2 * (-z1^2 + z2^2)^2)
    return dJz0dEso + s * dJz1dEso
end

"""`_dJzdLz(a, s, En, Lz, Kc)` — ∂J_z/∂Lz = dJz0dLso + s·dJz1dLso (.wl `dJzdLz`, L1624-1631)."""
function _dJzdLz(a, s, Eso, Lso, Kso)
    z1, z2, kzsq = _dJz_zroots(a, Eso, Lso, Kso)
    EK = real(ellipK(kzsq)); EE = real(ellipE(kzsq)); EPi = real(ellipPi(z2^2, kzsq))
    dJz0dLso = -((2 * (-a * Eso * EK + Lso * EPi)) / (a * sqrt(1 - Eso^2) * pi * z1))
    dJz1dLso = (2 * sqrt(Kso) * ((-z1^2 + (-1 + 2 * z1^2) * z2^2) * EE +
                (-z1 + z2) * (z1 + z2) * (-1 + z2^2) * EK)) /
               (a^2 * (1 - Eso^2)^(3 / 2) * pi * z1 * z2^2 * (-z1^2 + z2^2)^2)
    return dJz0dLso + s * dJz1dLso
end

"""
    _KCorrection(a, p, e, x, En, Lz, K) -> Real

Linear-in-spin correction Kc1 to the Carter-Killing constant K (.wl `KCorrection`, L1633-1653).
`(En,Lz,K)` are the GEODESIC constants; solves the 3×3 system M·{ES,LzS,KS} = −v built from the
radial/polar turning-point integrals, returns `KS + 2((K+a²)En − a·Lz)/√K`.
"""
function _KCorrection(a, p, e, x, Eng, Lzg, Kg)
    # radial roots of ((r²+a²)En − a·Lz)² − (r²−2r+a²)(K+r²): identical polynomial to radial_roots
    # (since K = Q + (Lz−aE)²), DESCENDING order r1≥r2≥r3≥r4 (.wl Reverse[Sort]). r1,r2 = apo/peri.
    Q = Kg - (Lzg - a * Eng)^2
    r1, r2, r3, r4 = radial_roots(a, p, e, x, Eng, Q)
    z1 = sqrt(1 - x^2)
    z2 = sqrt(a^2 * (1 - Eng^2) + Lzg^2 / x^2)
    kz = a^2 * (1 - Eng^2) * z1^2 / z2^2
    kr = real((r1 - r2) * (r3 - r4) / ((r1 - r3) * (r2 - r4)))

    IntRξtφ(r) = -2 * (a^2 * Eng - a * Lzg + Eng * r^2) / sqrt(Kg) *
        ((r * ((r - 1) * (Kg + r^2) + r * (r^2 - 2r + a^2)) - 2 * r^2 * Eng * ((r^2 + a^2) * Eng - a * Lzg)) *
            real(ellipPi(-a^2 * z1^2 / r^2, kz)) / real(ellipK(kz)) / r^2 +
         (Eng * ((r^2 + a^2) * Eng - a * Lzg) - r^2 + r))
    dRdEn(r) = 2 * (a^2 + r^2) * (-a * Lzg + Eng * (a^2 + r^2))
    dRdLz(r) = -2 * a * (-a * Lzg + Eng * (a^2 + r^2))
    dRdK(r)  = -a^2 - (-2 + r) * r
    dZdEn(z) = -2 * a * (-1 + z^2) * (Lzg + a * Eng * (-1 + z^2))
    dZdLz(z) = -2 * (Lzg + a * Eng * (-1 + z^2))
    dZdK(z)  = 1 - z^2

    az12IntΣ = (a * z1^2) / (r3^2 + a^2 * z1^2) +
        2 * real((im * (-r2 + r3) * z1) / (2 * (r2 + im * a * z1) * (r3 + im * a * z1)) *
                 ellipPi(((r1 - r2) * (r3 + im * a * z1)) / ((r1 - r3) * (r2 + im * a * z1)), kr)) /
        real(ellipK(kr))
    IntZξtφ = -2 * (Lzg + a * Eng * (-1 + z1^2)) / sqrt(Kg) *
        ((Eng * ((1 - z1^2) * a * Eng - Lzg) + a * z1^2) +
         ((Kg + a^2 * (1 - 2 * z1^2)) - 2 * a * Eng * ((1 - z1^2) * a * Eng - Lzg)) * az12IntΣ)

    # M·{ES,LzS,KS} = −v (.wl LinearSolve). M is real 3×3; KS via Cramer's rule (only KS is needed).
    a11 = real(dRdEn(r1)); a12 = real(dRdLz(r1)); a13 = real(dRdK(r1))
    a21 = real(dRdEn(r2)); a22 = real(dRdLz(r2)); a23 = real(dRdK(r2))
    a31 = real(dZdEn(z1)); a32 = real(dZdLz(z1)); a33 = real(dZdK(z1))
    b1 = -real(IntRξtφ(r1)); b2 = -real(IntRξtφ(r2)); b3 = -real(IntZξtφ)
    detM = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31)
    KS = (a11 * (a22 * b3 - b2 * a32) - a12 * (a21 * b3 - b2 * a31) + b1 * (a21 * a32 - a22 * a31)) / detM
    return KS + 2 * ((Kg + a^2) * Eng - a * Lzg) / sqrt(Kg)
end

# h1/h2/h3: projections of s^{μν} to ξ_{μ;ν} (Eq. 21 in 2303.16798); identical to the verified
# _oc_h* (spinning_orbit_correction_local.jl) — kept as the legacy `_fc_h*` names the mode files use.
_fc_h1(r, z, K, a) = (r * (-3 * a^2 * r^2 * z^2 + a^4 * z^4 + K * (r^2 - 3 * a^2 * z^2))) /
                     (sqrt(K) * (r^2 + a^2 * z^2)^3)
_fc_h2(r, z, K, En, Lz, a) = (-(En * Lz * r^6) + a^4 * En * Lz * r^2 * z^4 + a^2 * En * Lz * r^4 * (-2 + z^2) -
    a^6 * En * Lz * z^4 * (-2 + z^2) + a^7 * En^2 * z^4 * (-1 + z^2) +
    a * r^3 * (Lz^2 * r + K * (-1 + z^2) + K * r * (-1 + 2 * z^2) + r^3 * (z^2 - En^2 * (-1 + z^2))) +
    a^3 * (-(En^2 * r^4 * (-1 + z^2)) + r * z^2 * (2 * r^3 + 2 * K * r * z^2 - 3 * K * (-1 + z^2) - 3 * r^2 * (-1 + z^2))) +
    a^5 * z^4 * (K - Lz^2 + r * (-1 + z^2 + r * (2 - z^2 + En^2 * (-1 + z^2))))) /
    (sqrt(K) * (r^2 + a^2 * z^2)^3)
_fc_h3(r, z, K, a) = (-2 * a * r * z) / (sqrt(K) * (r^2 + a^2 * z^2)^2)

# --- fab^(i) basis (geodesic source structure); SOURCE-EXTRACTED+TRANSLATED from .wl L48-85. ---
# Returns (fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2).
function _fabi(zeta, zetabar, sinth, Del, dDel, KK, dKK, SS, L2SS, L1L2SS, aa)
    return (
        (-2 * (L1L2SS - ((2 * im) * aa * L2SS * sinth) / zeta) * zeta^2) / Del^2,
        (2 * sqrt(2) * zeta^2 * (-((aa * KK * sinth * SS * (-zeta^(-1) + zetabar^(-1))) / Del) + L2SS * ((im * KK) / Del + zeta^(-1) + zetabar^(-1)))) / (Del * zetabar),
        (2 * sqrt(2) * zeta^2 * (L2SS + im * aa * sinth * SS * (-zeta^(-1) + zetabar^(-1)))) / (Del * zetabar),
        (SS * (KK^2 + im * (Del * dKK - dDel * KK) - ((2 * im) * Del * KK) / zeta) * zeta^2) / (Del^2 * zetabar^2),
        (-2 * SS * ((im * KK) / Del + zeta^(-1)) * zeta^2) / zetabar^2,
        -((SS * zeta^2) / zetabar^2),
    )
end

# ∂/∂r of the fab^(i) basis (.wl dfabidr, L71-77).
function _dfabidr(zeta, zetabar, sinth, Del, dDel, KK, dKK, SS, L2SS, L1L2SS, aa, ww)
    return (
        (Del * ((4 * im) * aa * L2SS * sinth - 4 * L1L2SS * zeta) + 4 * dDel * zeta * ((-2 * im) * aa * L2SS * sinth + L1L2SS * zeta)) / Del^3,
        (2 * sqrt(2) * (aa * sinth * SS * (zeta - zetabar) * (2 * dDel * KK * zeta * zetabar - Del * (dKK * zeta * zetabar + KK * (-2 * zeta + zetabar))) + L2SS * ((-im) * KK * zeta * zetabar * (Del * (zeta - 2 * zetabar) + 2 * dDel * zeta * zetabar) + Del * (-(Del * (zeta - zetabar) * (2 * zeta + zetabar)) + zeta * zetabar * (im * dKK * zeta * zetabar - dDel * (zeta + zetabar)))))) / (Del^3 * zetabar^3),
        (2 * sqrt(2) * (-(L2SS * zeta * zetabar * (Del * (zeta - 2 * zetabar) + dDel * zeta * zetabar)) - im * aa * sinth * SS * (zeta - zetabar) * (Del * (2 * zeta - zetabar) + dDel * zeta * zetabar))) / (Del^2 * zetabar^3),
        ((-2 * im) * SS * ((-im) * KK^2 * zeta * (Del * (zeta - zetabar) + dDel * zeta * zetabar) + KK * (-(dDel^2 * zeta^2 * zetabar) + Del^2 * (-2 * zeta + zetabar) + Del * zeta^2 * (-dDel + zetabar + im * dKK * zetabar)) + Del * zeta^2 * (dDel * dKK * zetabar + Del * (dKK - ww * zetabar)))) / (Del^3 * zetabar^3),
        (2 * SS * (im * KK * zeta * (2 * Del * (zeta - zetabar) + dDel * zeta * zetabar) + Del * (Del * (2 * zeta - zetabar) - im * dKK * zeta^2 * zetabar))) / (Del^2 * zetabar^3),
        (2 * SS * zeta * (zeta - zetabar)) / zetabar^3,
    )
end

# ∂/∂θ of the fab^(i) basis (.wl dfabidθ, L79-85).
function _dfabidθ(zeta, zetabar, sinth, Del, dDel, KK, dKK, SS, L2SS, L1L2SS, dSdth, dL2Sdth, dL1L2Sdth, aa)
    return (
        (-4 * aa^2 * L2SS * sinth^2 + (4 * im) * aa * (dL2Sdth - L1L2SS) * sinth * zeta + 2 * zeta * (-((dL1L2Sdth + L2SS) * zeta) + L2SS * zetabar)) / Del^2,
        (sqrt(2) * ((-2 * im) * aa^2 * KK * sinth^2 * SS * (2 * zeta - zetabar) * (zeta + zetabar) + zeta * zetabar * ((-im) * KK * SS * (zeta - zetabar)^2 + (2 * im) * dL2Sdth * KK * zeta * zetabar + 2 * Del * dL2Sdth * (zeta + zetabar)) + 2 * aa * sinth * (dSdth * KK * zeta * zetabar * (-zeta + zetabar) + im * Del * L2SS * (zeta + zetabar) * (2 * zeta + zetabar) - KK * L2SS * zeta * zetabar * (zeta + 2 * zetabar)))) / (Del^2 * zetabar^3),
        (sqrt(2) * (-2 * aa^2 * sinth^2 * SS * (2 * zeta - zetabar) * (zeta + zetabar) - zeta * zetabar * (SS * (zeta - zetabar)^2 - 2 * dL2Sdth * zeta * zetabar) + (2 * im) * aa * sinth * zeta * zetabar * (dSdth * (zeta - zetabar) + L2SS * (zeta + 2 * zetabar)))) / (Del * zetabar^3),
        (Del * dKK * zeta * (im * dSdth * zeta * zetabar - 2 * aa * sinth * SS * (zeta + zetabar)) + KK^2 * zeta * (dSdth * zeta * zetabar + (2 * im) * aa * sinth * SS * (zeta + zetabar)) + KK * (dDel * zeta * ((-im) * dSdth * zeta * zetabar + 2 * aa * sinth * SS * (zeta + zetabar)) + 2 * Del * ((-im) * dSdth * zeta * zetabar + aa * sinth * SS * (2 * zeta + zetabar)))) / (Del^2 * zetabar^3),
        (-2 * (KK * zeta * (im * dSdth * zeta * zetabar - 2 * aa * sinth * SS * (zeta + zetabar)) + Del * (dSdth * zeta * zetabar + im * aa * sinth * SS * (2 * zeta + zetabar)))) / (Del * zetabar^3),
        -((zeta * (dSdth * zeta * zetabar + (2 * im) * aa * sinth * SS * (zeta + zetabar))) / zetabar^3),
    )
end

# ∂_S (ω-linear, S→S1) of the fab^(i) basis (.wl dfabidS, L55-61); analytic fab-basis ω-derivative.
# Returns (fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21).
function _dfabidS(ζ, ζbar, sinθ, Δ, dΔ, K, dK, K1, dK1, S, L2S, L1L2S, S1, L2S1, L1L2S1, a)
    return (
        -2 * ζ^2 / Δ^2 * (L1L2S1 - 2im * a / ζ * sinθ * L2S1),
        4 / sqrt(2) * ζ^2 / ζbar / Δ * ((im * K / Δ + 1 / ζ + 1 / ζbar) * L2S1 + (im * K1 / Δ) * L2S -
            a * sinθ * K / Δ * (1 / ζbar - 1 / ζ) * S1 - a * sinθ * K1 / Δ * (1 / ζbar - 1 / ζ) * S),
        4 / sqrt(2) * ζ^2 / ζbar / Δ * (L2S1 + im * a * sinθ * (1 / ζbar - 1 / ζ) * S1),
        ζ^2 / ζbar^2 * ((im * (dK * Δ - K * dΔ) - 2 / ζ * im * K * Δ + K^2) / Δ^2 * S1 +
            (im * (dK1 * Δ - K1 * dΔ) - 2 / ζ * im * K1 * Δ + 2 * K * K1) / Δ^2 * S),
        -2 * ζ^2 / ζbar^2 * ((1 / ζ + im * K / Δ) * S1 + (im * K1 / Δ) * S),
        -ζ^2 / ζbar^2 * S1,
    )
end

# ω-derivative of the fab^(i) basis (.wl dfabidω, L62-69); the analytic counterpart of `_dfabidS`.
# Returns (dfnn0dω,dfnmb0dω,dfnmb1dω,dfmbmb0dω,dfmbmb1dω,dfmbmb2dω).
function _dfabidω(ζ, ζbar, sinθ, Δ, dΔ, K, dKdr, dKdω, d2Kdrdω, S, L2S, L1L2S, dSdω, dL2Sdω, dL1L2Sdω, a)
    return (
        -2 * ζ^2 / Δ^2 * (dL1L2Sdω - 2im * a / ζ * sinθ * dL2Sdω),
        4 / sqrt(2) * ζ^2 / ζbar / Δ * ((im * K / Δ + 1 / ζ + 1 / ζbar) * dL2Sdω + (im * dKdω / Δ) * L2S -
            a * sinθ * K / Δ * (1 / ζbar - 1 / ζ) * dSdω - a * sinθ * dKdω / Δ * (1 / ζbar - 1 / ζ) * S),
        4 / sqrt(2) * ζ^2 / ζbar / Δ * (dL2Sdω + im * a * sinθ * (1 / ζbar - 1 / ζ) * dSdω),
        ζ^2 / ζbar^2 * ((im * (dKdr * Δ - K * dΔ) - 2 / ζ * im * K * Δ + K^2) / Δ^2 * dSdω +
            (im * (d2Kdrdω * Δ - dKdω * dΔ) + 2 * (-1 / ζ * im * Δ + K) * dKdω) / Δ^2 * S),
        -2 * ζ^2 / ζbar^2 * ((1 / ζ + im * K / Δ) * dSdω + (im * dKdω / Δ) * S),
        -ζ^2 / ζbar^2 * dSdω,
    )
end

"""
    _angular_fromcorr(SWSH, λ, a, m, ω, zp) -> NamedTuple

Polar functions and angular operators along the geodesic, in the convention used by
TeukolskySpinModeFromCorrection (.wl L1733-1742) — the OLDER `L2S`/`L1L2S` convention (differs from
`_angular_generic`): `(S, dSdθ, d2Sdθ2, d3Sdθ3, L2S, dL2Sdθ, L1L2S, dL1L2Sdθ)`.
`SWSH` is a callable `SpinWeightedSpheroidalHarmonicS(-2,l,m,aω)` object.
"""
function _angular_fromcorr(SWSH, λ, a::Real, m::Integer, ω::Real, zp::Real)
    θp = acos(zp); sinθp = sqrt(1 - zp^2)
    S = SWSH(θp, 0.0)
    dSdθ = swsh_derivative(SWSH, 1, 0, θp, 0.0)
    d2Sdθ2 = -(-(a * ω)^2 * (1 - zp^2) - (m - 2 * zp)^2 / (1 - zp^2) + 4 * a * ω * zp - 2 + 2 * m * a * ω + λ) * S -
             zp / sinθp * dSdθ
    d3Sdθ3 = -(1 / sinθp^3) * 2 * (-2 + m * zp + a * zp * (-1 + zp^2) * ω) * (m + a * ω - zp * (2 + a * zp * ω)) * S -
             (-(a * ω)^2 * (1 - zp^2) - (m - 2 * zp)^2 / (1 - zp^2) + 4 * a * ω * zp - 2 + 2 * m * a * ω + λ - 1 / (1 - zp^2)) * dSdθ -
             zp / sinθp * d2Sdθ2
    L2S = dSdθ - (S * (m - 2 * zp + a * (-1 + zp^2) * ω)) / sinθp
    dL2Sdθ = d2Sdθ2 + 1 / (-1 + zp^2) * (dSdθ * sinθp * (m - 2 * zp + a * (-1 + zp^2) * ω) +
              S * (2 - m * zp + a * zp * (-1 + zp^2) * ω))
    L1L2S = d2Sdθ2 + (dSdθ * (-2 * m + 3 * zp - 2 * a * (-1 + zp^2) * ω)) / sinθp +
            S * (-2 - (m * (m - 2 * zp)) / (-1 + zp^2) - 2 * a * (m - 2 * zp) * ω - a^2 * (-1 + zp^2) * ω^2)
    dL1L2Sdθ = d3Sdθ3 +
        1 / (-1 + zp^2) * dSdθ * (5 - m^2 - 2 * zp^2 - 2 * a * (m - 3 * zp) * (-1 + zp^2) * ω - a^2 * (-1 + zp^2)^2 * ω^2) +
        1 / sinθp^3 * (d2Sdθ2 * (-1 + zp^2) * (2 * m - 3 * zp + 2 * a * (-1 + zp^2) * ω) +
              2 * S * (-m^2 * zp + m * (1 + zp^2) + a * (-1 + zp^2)^2 * ω * (-2 + a * zp * ω)))
    return (S = S, dSdθ = dSdθ, d2Sdθ2 = d2Sdθ2, d3Sdθ3 = d3Sdθ3,
            L2S = L2S, dL2Sdθ = dL2Sdθ, L1L2S = L1L2S, dL1L2Sdθ = dL1L2Sdθ)
end

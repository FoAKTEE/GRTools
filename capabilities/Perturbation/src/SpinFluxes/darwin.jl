# Linear-in-spin GW flux, ECCENTRIC EQUATORIAL s=-2, DARWIN χr parametrization — faithful port of
# TeukolskySpinFluxes.wl TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin (L1004-1199). The
# Darwin-reparametrized sibling of equatorial_spin_flux (Mino-phase): same FabEq/G_a/∂_ωF building
# blocks (REUSED), but the χr grid {π/s,3π/s,…,(s-1)π/s} with the Jacobian dλ/dχr and the prefactor
# Ωr=Υr/Υt replacing the Mino 1/Υt. χr→qr map + deltas from the LOCAL `_darwin_trajectory_eq`
# (1:1 of KerrGeodesics.KerrGeoTrajectoryEquatorialDarwin / KerrGeoOrbit.m L149-186 — Spacetime
# ships no Darwin-parametrized trajectory, so it lives here; see the FIREWALL NOTE in
# spinning_secondary_local.jl). R3.

"""
    _darwin_trajectory_eq(a, p, e) -> Function

Equatorial geodesic trajectory in the **Darwin** χ-parametrization (M=1). Returns
`χ -> (; r, t, φ, Δtr, Δφr, λ, qr)` with `r(χ)=p/(1+e cosχ)`; `t`,`φ` the full BL coordinates and
`Δtr=t−Υt·λ(χ)`, `Δφr=φ−Υφ·λ(χ)` the oscillating deviations (KGO convention); `λ(χ)` the Mino time
and `qr=Υr·λ` the radial Mino phase. Faithful 1:1 of `KerrGeoOrbit.m` `KerrGeoOrbitEquatorialDarwin`
(L149-186): the Mino trajectory composed with the reparametrisation `λ(χ)` (closed-form `EllipticF`).
"""
function _darwin_trajectory_eq(a, p, e)
    En, L, Q = constants_of_motion(a, p, e, 1)
    r1c, r2c, r3c, r4c = radial_roots(a, p, e, 1, En, Q)
    r1 = real(r1c); r2 = real(r2c); r3 = real(r3c); r4 = real(r4c)
    fr = kerr_frequencies(a, p, e, 1; Time="Mino"); Υr = fr.Ur; Υt = fr.Ut; Υφ = fr.Uph
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    Λr = 2pi / Υr
    pref = 2 / (sqrt(1 - En^2) * sqrt((r1 - r3) * (r2 - r4)))
    yr(r) = e == 0 ? 1.0 : sqrt((r1 - r3) / (r1 - r2) * (r - r2) / (r - r3))
    λ0r(r) = pref * real(ellipFinc(asin(clamp(yr(r), -1.0, 1.0)), kr))
    Λr1 = λ0r(r2)
    rad = geodesic_trajectory_radial(a, p, e, 1)          # KGO-convention Δtr/Δφr (deviation)
    return function (χ)
        rχ = p / (1 + e * cos(χ))
        mχ = mod(χ, 2pi)
        λ = Λr * floor(χ / (2pi)) + (mχ <= pi ? λ0r(rχ) - Λr1 : Λr - λ0r(rχ))
        qr = Υr * λ                                      # Mino phase at this χ
        d = rad(qr)
        return (r = rχ, t = Υt * λ + d.Δtr, φ = Υφ * λ + d.Δφr, Δtr = d.Δtr, Δφr = d.Δφr,
                λ = λ, qr = qr)
    end
end

"""
    _spin_mode_eq_correction_darwin(l, m, n, a, p, e; Nmax=80, h=3e-4) -> NamedTuple

Linear-in-spin GW flux correction (dE¹,dLz¹ at ℐ/ℋ) for an ECCENTRIC EQUATORIAL s=-2 point particle
in the DARWIN χr parametrization, per `TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin`
(.wl L1004-1199). The χr sum {π/s,…,(s-1)π/s} with the Jacobian dλ/dχr + the two velocity branches;
∂_ωF/∂_ωW via central FD. Returns `(; dE1I,dE1H,dL1I,dL1H, CPlus0,CMinus0,CPlus1,CMinus1, omega,
stepsr)`. R3.
"""
function _spin_mode_eq_correction_darwin(l::Integer, m::Integer, n::Integer, a::Real, p::Real, e::Real;
        Nmax::Integer=80, h::Real=3e-4)
    En0 = geodesic_energy(a, p, e, 1); Lz0 = geodesic_angular_momentum(a, p, e, 1)
    Q0 = carter_constant(a, p, e, 1); K0 = Q0 + (Lz0 - a * En0)^2; sK = sqrt(K0)
    fr = kerr_frequencies(a, p, e, 1; Time="Mino")          # Υr, Υθ, Υφ, Υt (geodesic)
    Uτ = proper_frequency_factor(a, p, e, 1)
    Υr = fr.Ur; Υt = fr.Ut; Υφ = fr.Uph
    Ωr = Υr / Υt; Ωφ = Υφ / Υt
    ω = m * Ωφ + n * Ωr
    ω1 = 3 * Uτ * ω / (2 * sK * Υt)
    Υt1 = -3 * Uτ / (2 * sK)
    En, _, Q = constants_of_motion(a, p, e, 1)          # 3rd/4th radial roots (.wl L1026)
    roots = radial_roots(a, p, e, 1, En, Q); r3 = real(roots[3]); r4 = real(roots[4])
    darwin = _darwin_trajectory_eq(a, p, e)       # χ→(t,φ,λ,qr) map (closed-form EllipticF)
    eqtraj = geodesic_trajectory_equatorial(a, p, e)             # Mino Δτr(qr) for the proper-time deviation
    # Darwin-angle secular-subtracted deviations (.wl orbit["TrajectoryDeltas"] for the Darwin orbit):
    # the radial harmonic phase is n·χ, so the secular drift is removed PER χ — Δtr_D=t−(Υt/Υr)χ etc.
    # This makes ω·Δtr_D − m·Δφr_D + n·χ reconstruct the exact physical phase ω·t − m·φ (verified).
    function deltas_D(χ)
        dv = darwin(χ); λ = dv.λ
        τ = Uτ * λ + eqtraj(dv.qr).Δτr                       # proper time τ(χ)=Υτ·λ + Δτr_mino(qr)
        return (Δtr = dv.t - (Υt / Υr) * χ, Δφr = dv.φ - (Υφ / Υr) * χ, Δτr = τ - (Uτ / Υr) * χ)
    end
    # stepsr power-of-2 heuristic (.wl L1071-1072: 2^(5+⌈log2|ω·Δtr'−m·Δφr'+n|⌉) at χ=π,0, ≥32)
    dphase(χ) = (t1 = deltas_D(χ + 1e-6); t0 = deltas_D(χ - 1e-6);
                 ω * (t1.Δtr - t0.Δtr) / 2e-6 - m * (t1.Δφr - t0.Δφr) / 2e-6 + n)
    stepterm(χ) = (d = abs(dphase(χ)); d < 1e-12 ? 0.0 : 2.0^(5 + ceil(log2(d))))
    stepsr = Int(max(stepterm(pi), stepterm(0.0), 32.0))
    b0 = _eqflux_bundle(l, m, a, ω; Nmax=Nmax)
    bp = _eqflux_bundle(l, m, a, ω + h; Nmax=Nmax)
    bm = _eqflux_bundle(l, m, a, ω - h; Nmax=Nmax)
    vmb = -im * (a * En0 - Lz0)
    contractF(vn, F) = vn^2 * F.Fnn + vn * vmb * F.Fnmb + vmb^2 * F.Fmbmb
    sumP0 = zero(ComplexF64); sumP1 = zero(ComplexF64); sumM0 = zero(ComplexF64); sumM1 = zero(ComplexF64)
    W = zero(ComplexF64); dWdω = zero(ComplexF64)
    npts = stepsr ÷ 2
    for i in 1:npts
        χr = (2i - 1) * pi / stepsr                           # Darwin grid {π/s,3π/s,…,(s-1)π/s}
        cosχ = cos(χr); sinχ = sqrt(1 - cosχ^2)
        rp = p / (1 + e * cosχ)
        dλdχr = sqrt((1 - e^2) / ((1 - En0^2) * (-p + r3 + e * r3 * cosχ) * (-p + r4 + e * r4 * cosχ)))
        Urp = p * e * sinχ / (1 + e * cosχ)^2 / dλdχr         # geodesic dr/dλ (.wl L1083)
        d = deltas_D(χr)                                       # Darwin-angle secular-subtracted deltas
        vnp = -((rp^2 + a^2) * En0 - a * Lz0 + Urp) / (rp^2 - 2rp + a^2)
        vnm = -((rp^2 + a^2) * En0 - a * Lz0 - Urp) / (rp^2 - 2rp + a^2)
        expr  = exp(im * (ω * d.Δtr - m * d.Δφr + n * χr))    # n·χr (Darwin angle), with Darwin-subtracted Δtr/Δφr
        expr1 = im * (ω1 * d.Δtr + ω * (-3 / (2 * sK)) * d.Δτr)
        FinP = _eqflux_fab(b0.R.In, rp, a, m, ω, b0.λ, b0.ang)
        FupP = _eqflux_fab(b0.R.Up, rp, a, m, ω, b0.λ, b0.ang)
        dFab(Rfp, λp, angp, Rfm, λm, angm) = let
            Fp = _eqflux_fab(Rfp, rp, a, m, ω + h, λp, angp); Fm = _eqflux_fab(Rfm, rp, a, m, ω - h, λm, angm)
            (Fnn = (Fp.Fnn - Fm.Fnn) / (2h), Fnmb = (Fp.Fnmb - Fm.Fnmb) / (2h), Fmbmb = (Fp.Fmbmb - Fm.Fmbmb) / (2h))
        end
        dFin = dFab(bp.R.In, bp.λ, bp.ang, bm.R.In, bm.λ, bm.ang)
        dFup = dFab(bp.R.Up, bp.λ, bp.ang, bm.R.Up, bm.λ, bm.ang)
        # geodesic (Jacobian-weighted; +χr branch ×expr, conjugate branch /expr)
        sumP0 += dλdχr * (contractF(vnp, FinP) * expr + contractF(vnm, FinP) / expr)
        sumM0 += dλdχr * (contractF(vnp, FupP) * expr + contractF(vnm, FupP) / expr)
        # linear-in-spin: G_a source + phase shift (±expr1) + ∂_ωF·ω1, per velocity branch
        GPp = (vnp * FinP.Gn + vmb * FinP.Gmb) / sK; GPm = (vnm * FinP.Gn + vmb * FinP.Gmb) / sK
        GUp = (vnp * FupP.Gn + vmb * FupP.Gmb) / sK; GUm = (vnm * FupP.Gn + vmb * FupP.Gmb) / sK
        sumP1 += dλdχr * ((GPp + contractF(vnp, FinP) * expr1 + contractF(vnp, dFin) * ω1) * expr +
                          (GPm + contractF(vnm, FinP) * (-expr1) + contractF(vnm, dFin) * ω1) / expr)
        sumM1 += dλdχr * ((GUp + contractF(vnp, FupP) * expr1 + contractF(vnp, dFup) * ω1) * expr +
                          (GUm + contractF(vnm, FupP) * (-expr1) + contractF(vnm, dFup) * ω1) / expr)
        if i == 1                                              # W at the FIRST χr point (.wl L1161)
            W = _eqflux_W(b0.R, rp, a)
            dWdω = (_eqflux_W(bp.R, rp, a) - _eqflux_W(bm.R, rp, a)) / (2h)
        end
    end
    # Darwin normalisation: ·Ωr / (W·stepsr) — the dλdχr Jacobian + Ωr replaces the Mino 1/Υt (.wl L1163)
    CPlus0  = -2pi * sumP0 * Ωr / (W * stepsr); CMinus0 = -2pi * sumM0 * Ωr / (W * stepsr)
    W1 = dWdω * ω1
    pref = En0 / (2 * sK) - Υt1 / Υt - W1 / W
    CPlus1  = -2pi * sumP1 * Ωr / (W * stepsr) + pref * CPlus0
    CMinus1 = -2pi * sumM1 * Ωr / (W * stepsr) + pref * CMinus0
    α = _alpha_circ(a, m, ω, b0.λ)
    dαdω = (_alpha_circ(a, m, ω + h, bp.λ) - _alpha_circ(a, m, ω - h, bm.λ)) / (2h)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    return (dE1I=dE1I, dE1H=dE1H, dL1I=dL1I, dL1H=dL1H, CPlus0=CPlus0, CMinus0=CMinus0,
            CPlus1=CPlus1, CMinus1=CMinus1, omega=ω, stepsr=stepsr)
end

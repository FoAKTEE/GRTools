# Spherical-orbit linear-in-spin flux CORRECTION (geodesic mode + dE¹/dLz¹ + CarterConstantK), ANALYTIC
# ω-derivative path — TeukolskySpinFluxes.wl TeukolskySpinModeSphericalCorrectionAnalytical (L4472-4672).
# The θ-only restriction of inclined_correction.jl: a SINGLE wz (polar Mino) sum at fixed r≡p (Ur=0,
# radials + 𝒟-operators computed ONCE outside the loop) over 2 polar quadrants. Reuses the inclined
# _Fab / _incl_radialD / _angular_generic / _incl_bundle / _alpha_circ / _dδ helpers verbatim; ω1=
# 3Υτω/(2√K Υt), Υt1=−3Υτ/(2√K), polar phase expz1; ∂_ωF/∂_ωW/∂_ωα via central FD of the ω-bundle. R3.

"""
    _spin_mode_spherical_correction_analytical(l, m, k, a, p, x; Nmax=80, h=3e-4) -> NamedTuple

Spherical (e=0) linear-in-spin s=−2 flux CORRECTION — TeukolskySpinModeSphericalCorrectionAnalytical
(.wl L4472-4672). Geodesic orbit built internally; the analytic angparNew/HS1spin ω-derivatives are
replaced by central FD of the ω-bundle (`_incl_bundle` at ω±h), matching inclined_correction.jl.
Returns `(; l,m,k,n=0, omega, CPlus0,CMinus0,CPlus1,CMinus1, α, Spi2, EnergyI,EnergyH,AngMomI,AngMomH,
dE1I,dE1H,dL1I,dL1H, CarterKI,CarterKH, dCarterKI,dCarterKH, stepsz)`. R3.
"""
function _spin_mode_spherical_correction_analytical(l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, x::Real; Nmax::Integer=80, h::Real=3e-4)
    (l < 2 || abs(m) > l) && error("TeukolskySpinModeSphericalCorrectionAnalytical: l<2 or |m|>l")
    e = 0.0
    En0, Lz0, Q0 = constants_of_motion(a, p, e, x)
    K0 = Q0 + (Lz0 - a * En0)^2; sK = sqrt(K0)
    fr = kerr_frequencies(a, p, e, x; Time="Mino")
    pt = geodesic_proper_time(a, p, e, x)
    Υr = fr.Ur; Υz = fr.Uth; Υφ = fr.Uph; Υt = fr.Ut; Υτ = pt.ProperTimeFrequency
    Ωz = Υz / Υt; Ωφ = Υφ / Υt
    ω = m * Ωφ + k * Ωz
    ω1 = 3 * Υτ * ω / (2 * sK * Υt); Υt1 = -3 * Υτ / (2 * sK)
    # Carter-K theta Mino-frequency pieces (.wl L4498-4500): geodesic Q0 = K0 − (Lz0 − aEn0)²
    zp_r, zm_r = polar_roots(a, p, e, x)
    Υtz = Spacetime._mino_freq_t_theta(a, En0, Q0, zp_r, zm_r, x) + a * Lz0
    Υφz = Spacetime._mino_freq_phi_theta(a, En0, Lz0, zp_r, zm_r) - a * En0
    Υτz = pt.Υτz
    pol = geodesic_trajectory_polar(a, p, e, x)
    z = qz -> cos(pol(qz).θ)
    Δtθ = qz -> pol(qz).Δtz; Δφθ = qz -> pol(qz).Δφz
    # θ-only step heuristic (.wl L4527-4528)
    stepsz = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/4) - m * _dδ(Δφθ, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npz = stepsz ÷ 2

    # radial bundle at ω and ω±h (MST-analytic radial; ω-derivative via central FD)
    b0 = _incl_bundle(l, m, a, ω; Nmax=Nmax)
    bp = _incl_bundle(l, m, a, ω + h; Nmax=Nmax); bm = _incl_bundle(l, m, a, ω - h; Nmax=Nmax)

    # --- radials + 𝒟-operators computed ONCE at rp=p, OUTSIDE the wz loop (.wl L4530-4558; Ur=0) ---
    rp = float(p); Δ = rp^2 - 2rp + a^2
    Din0 = _incl_radialD(b0.R.In, rp, a, m, ω, b0.λ); Dup0 = _incl_radialD(b0.R.Up, rp, a, m, ω, b0.λ)
    Dinp = _incl_radialD(bp.R.In, rp, a, m, ω + h, bp.λ); Dinm = _incl_radialD(bm.R.In, rp, a, m, ω - h, bm.λ)
    Dupp = _incl_radialD(bp.R.Up, rp, a, m, ω + h, bp.λ); Dupm = _incl_radialD(bm.R.Up, rp, a, m, ω - h, bm.λ)
    # invariant Wronskian + its ω-linear part (.wl L4611-4612), FD in ω
    W  = (Din0.R * Dup0.dR - Din0.dR * Dup0.R) / Δ
    Wp = (Dinp.R * Dupp.dR - Dinp.dR * Dupp.R) / Δ
    Wm = (Dinm.R * Dupm.dR - Dinm.dR * Dupm.R) / Δ
    W1 = (Wp - Wm) / (2h) * ω1

    # four-velocity radial component (scalar; Ur=0, .wl L4586)
    vn = -((rp^2 + a^2) * En0 - a * Lz0) / Δ

    sumP0 = sumP1 = sumM0 = sumM1 = zero(ComplexF64)
    integ0(vmb, F) = vn^2 * F.Fnn + vn * vmb * F.Fnmb + vmb^2 * F.Fmbmb
    Gop(vmb, F) = (vn * F.Gn + vmb * F.Gmb) / sK
    dF0(vmb, dF) = vn^2 * dF.Fnn + vn * vmb * dF.Fnmb + vmb^2 * dF.Fmbmb
    for iz in 1:npz
        wz = (iz - 1//2) * 2pi / stepsz
        zp = z(wz); sinθp = sqrt(1 - zp^2)
        sqrtZ = sqrt(max(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (K0 - a^2 * zp^2), 0.0))
        expz = exp(im * (ω * Δtθ(wz) - m * Δφθ(wz) + 2pi * k * (iz - 1//2) / stepsz))
        # polar linear phase (.wl L4566): basez carries the {1,-1} quadrant sign
        basez = 3 * im * ω / (2 * sK) * (Υτ / Υt * Δtθ(wz) - pt.Δτz(wz))
        ζ = rp - im * a * zp; ζbar = rp + im * a * zp
        # angular layer at ω and ω±h
        a0 = _angular_generic(l, m, a, ω, zp; Nmax=Nmax)
        ap = _angular_generic(l, m, a, ω + h, zp; Nmax=Nmax)
        am = _angular_generic(l, m, a, ω - h, zp; Nmax=Nmax)
        FIn = _Fab(ζ, ζbar, a, sinθp, Din0.R, Din0.DR, Din0.DDR, a0.S, a0.L2S, a0.L1L2S)
        FUp = _Fab(ζ, ζbar, a, sinθp, Dup0.R, Dup0.DR, Dup0.DDR, a0.S, a0.L2S, a0.L1L2S)
        # ∂_ω F components via central FD (radial + angular both at ω±h)
        FInp = _Fab(ζ, ζbar, a, sinθp, Dinp.R, Dinp.DR, Dinp.DDR, ap.S, ap.L2S, ap.L1L2S)
        FInm = _Fab(ζ, ζbar, a, sinθp, Dinm.R, Dinm.DR, Dinm.DDR, am.S, am.L2S, am.L1L2S)
        FUpp = _Fab(ζ, ζbar, a, sinθp, Dupp.R, Dupp.DR, Dupp.DDR, ap.S, ap.L2S, ap.L1L2S)
        FUpm = _Fab(ζ, ζbar, a, sinθp, Dupm.R, Dupm.DR, Dupm.DDR, am.S, am.L2S, am.L1L2S)
        dFIn = (Fnn=(FInp.Fnn-FInm.Fnn)/(2h), Fnmb=(FInp.Fnmb-FInm.Fnmb)/(2h), Fmbmb=(FInp.Fmbmb-FInm.Fmbmb)/(2h))
        dFUp = (Fnn=(FUpp.Fnn-FUpm.Fnn)/(2h), Fnmb=(FUpp.Fnmb-FUpm.Fnmb)/(2h), Fmbmb=(FUpp.Fmbmb-FUpm.Fmbmb)/(2h))
        # 2 polar quadrants (.wl: Uzp={1,-1}·(−√Z), expz^{1,-1}, expz1={1,-1}·basez)
        vmbbase = -im * (a * sinθp^2 * En0 - Lz0)
        for (vmb, ef, p1) in (((vmbbase - sqrtZ) / sinθp^2, expz,   basez),
                              ((vmbbase + sqrtZ) / sinθp^2, 1/expz, -basez))
            sumP0 += integ0(vmb, FIn) * ef
            sumM0 += integ0(vmb, FUp) * ef
            sumP1 += (Gop(vmb, FIn) + integ0(vmb, FIn) * p1 + dF0(vmb, dFIn) * ω1) * ef
            sumM1 += (Gop(vmb, FUp) + integ0(vmb, FUp) * p1 + dF0(vmb, dFUp) * ω1) * ef
        end
    end

    den = Υt * W * stepsz
    CPlus0  = -2pi * sumP0 / den; CMinus0 = -2pi * sumM0 / den
    pref = En0 / (2 * sK) - Υt1 / Υt - W1 / W
    CPlus1  = -2pi * sumP1 / den + pref * CPlus0
    CMinus1 = -2pi * sumM1 / den + pref * CMinus0

    α  = _alpha_circ(a, m, ω, b0.λ)
    dα = (_alpha_circ(a, m, ω + h, bp.λ) - _alpha_circ(a, m, ω - h, bm.λ)) / (2h)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    # geodesic Fluxes (.wl L4639-4651)
    EnergyI = abs2(CPlus0) / (4pi * ω^2);  EnergyH = α * abs2(CMinus0) / (4pi * ω^2)
    AngMomI = abs2(CPlus0) * m / (4pi * ω^3); AngMomH = α * abs2(CMinus0) * m / (4pi * ω^3)
    cK = k * Υz - ω * Υtz + m * Υφz
    CarterKI = cK * abs2(CPlus0) / (2pi * ω^3)
    CarterKH = cK * α * abs2(CMinus0) / (2pi * ω^3)
    # FluxesCorrection dE¹/dLz¹ (.wl L4656-4662)
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (dα * ω1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (dα * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    # CarterConstantK flux correction (.wl L4666-4667)
    prefK = (En0 * cK - a * (m - a * ω) - 3 * ω * (-Υτz + Υτ / Υt * Υtz)) / (2 * sK)
    dCarterKI = prefK * abs2(CPlus0) / (2pi * ω^3) + cK * (reP - 3 * abs2(CPlus0) * ω1 / ω) / (2pi * ω^3)
    dCarterKH = prefK * α * abs2(CMinus0) / (2pi * ω^3) +
                cK * α * (dα * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) / (2pi * ω^3)

    Sc = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω; numterms = Nmax > 60 ? :auto : Nmax)
    return (l=l, m=m, k=k, n=0, omega=ω, CPlus0=CPlus0, CMinus0=CMinus0, CPlus1=CPlus1, CMinus1=CMinus1,
            α=α, Spi2=Sc(pi/2, 0.0),
            EnergyI=EnergyI, EnergyH=EnergyH, AngMomI=AngMomI, AngMomH=AngMomH,
            CarterKI=CarterKI, CarterKH=CarterKH,
            dE1I=dE1I, dE1H=dE1H, dL1I=dL1I, dL1H=dL1H, dCarterKI=dCarterKI, dCarterKH=dCarterKH,
            stepsz=stepsz)
end

# Generic (inclined) linear-in-spin flux source layer + full mode + flux correction — faithful 1:1
# of TeukolskySpinFluxes/src/{inclined_flux,inclined_correction}.jl (TeukolskySpinFluxes.wl
# TeukolskySpinModeAnalytical L1879-2015, TeukolskySpinModeCorrectionAnalytical L2864-3090),
# rebuilt on the Perturbation Radial/AmplitudesFluxes machinery + Spacetime geodesics. The
# z-dependent angular operators L2S/dL2Sdθ/L1L2S, the generic Fab (Fnn,Fnmb,Fmbmb geodesic +
# Gn,Gmb spin-dipole), the wr×wz Mino double-sum, and the DH P5 correction path. M = 1.

"""
    _angular_generic(l, m, a, omega, zp; Nmax=80) -> NamedTuple

s=-2 spin-weighted angular operators at polar position `zp=cosθ` along the geodesic, per
TeukolskySpinFluxes.wl L1953-1958: `(S, dSdθ, d2Sdθ2, L2S, dL2Sdθ, L1L2S)`. `d2Sdθ2` is the
Teukolsky angular-equation form (eigenvalue λ); `Kθ=aω(1−z²)−m`.
"""
function _angular_generic(l::Integer, m::Integer, a::Real, omega::Real, zp::Real; Nmax::Integer=80)
    s = -2
    θp = acos(zp); sinθp = sqrt(1 - zp^2)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms = Nmax > 60 ? :auto : Nmax)
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
    S = Sc(θp, 0.0)
    dSdθ = swsh_derivative(Sc, 1, 0, θp, 0.0)
    Kθ = a * omega * (1 - zp^2) - m
    d2Sdθ2 = -(-(a * omega)^2 * (1 - zp^2) - (m - 2 * zp)^2 / (1 - zp^2) + 4 * a * omega * zp
                - 2 + 2 * m * a * omega + λ) * S - zp / sinθp * dSdθ
    L2S = dSdθ + (Kθ + 2 * zp) * S / sinθp
    dL2Sdθ = d2Sdθ2 + (Kθ + 2 * zp) * dSdθ / sinθp +
             (2 * (a * omega * zp - 1) - (Kθ + 2 * zp) / sinθp^2 * zp) * S
    L1L2S = dL2Sdθ + (Kθ + zp) * L2S / sinθp
    return (S = S, dSdθ = dSdθ, d2Sdθ2 = d2Sdθ2, L2S = L2S, dL2Sdθ = dL2Sdθ, L1L2S = L1L2S)
end

"""
    _Fab(ζ, ζbar, a, sinθp, RIn, DRIn, DDRIn, S, L2S, L1L2S) -> NamedTuple

Generic source projection (TeukolskySpinFluxes.wl `Fab`, L87-91): the three geodesic components
`(Fnn,Fnmb,Fmbmb)` + the two spin-dipole components `(Gn,Gmb)`, given the GHP radial 𝒟-operator
values `DRIn=𝒟R`, `DDRIn=𝒟²R` and angular `S,L2S,L1L2S`. Reduces to `_FabEq` at z=0 (ζ=ζbar=r).
"""
function _Fab(ζ, ζbar, a, sinθp, RIn, DRIn, DDRIn, S, L2S, L1L2S)
    Fnn   = -ζ / (2 * ζbar) * (L1L2S - 2 * im * a * sinθp * L2S / ζ) * RIn
    Fnmb  = (sinθp * ζ) / ζbar * L2S * DRIn - (sinθp * (ζ + ζbar)) / ζbar^2 * L2S * RIn +
            (im * a * sinθp^2 * (ζ - ζbar)) / ζbar^2 * S * DRIn
    Fmbmb = -ζ / (2 * ζbar) * S * sinθp^2 * (DDRIn - 2 / ζ * DRIn)
    Gn    = 1//2 * im * DRIn * (a * L2S * sinθp + im * L1L2S * ζ)
    Gmb   = -1//2 * L2S * sinθp * (DRIn - DDRIn * ζ)
    return (Fnn = Fnn, Fnmb = Fnmb, Fmbmb = Fmbmb, Gn = Gn, Gmb = Gmb)
end

"GHP 𝒟R=dR−iK/Δ·R and 𝒟²R at radius rp for a radial branch (In/Up), frequency ω, eigenvalue λ."
function _incl_radialD(Rbranch, rp, a, m, ω, λ)
    Δ = rp^2 - 2rp + a^2; K = (rp^2 + a^2) * ω - a * m; dKdr = 2 * rp * ω; dΔdr = 2 * rp - 2
    Rv = Rbranch(rp); dRv = teukolsky_radial_deriv(Rbranch, rp)
    d2Rv = _d2R_ode(Rv, dRv, -2, m, a, rp, ω, λ)
    DR = dRv - im * K / Δ * Rv
    dDRdr = d2Rv - im * K / Δ * dRv - im * (dKdr * Δ - K * dΔdr) / Δ^2 * Rv
    DDR = dDRdr - im * K / Δ * DR
    return (R = Rv, DR = DR, DDR = DDR, dR = dRv)
end

"central-FD derivative of an orbit TrajectoryDeltas closure (for the stepsr/stepsz heuristic)."
_dδ(f, q; h=1e-6) = (f(q + h) - f(q - h)) / (2h)

"""
    _inclined_spin_mode(orbit, l, m, n, k; Nmax=80) -> NamedTuple

Full inclined linear-in-spin s=−2 mode — TeukolskySpinModeAnalytical (L1879-2013): the wr×wz Mino
double-sum over 4 velocity quadrants of (vᵃvᵇF_ab + σ/√K̃·vᵃG_a)·exp phases, giving CPlus/CMinus and
the ℐ/ℋ Energy/AngularMomentum/Carter fluxes. `orbit` = inclined `KerrSpinOrbit`.
"""
function _inclined_spin_mode(orbit, l::Integer, m::Integer, n::Integer, k::Integer; Nmax::Integer=80)
    s = -2
    a = orbit.a; spar = orbit.spar
    EnT = orbit.EnTilde; LzT = orbit.JzTilde; KT = orbit.KTilde; sKT = sqrt(KT)
    Ωr, Ωz, Ωφ = orbit.BLFrequencies
    Υt = orbit.Frequencies.Ut; Υz = orbit.Frequencies.Uz; Υtz = orbit.Upsilontz; Υφz = orbit.Upsilonphiz
    rfun = orbit.trajectory[1]; zfun = orbit.trajectory[2]
    Δtr = orbit.TrajectoryDeltas.Δtr; Δtz = orbit.TrajectoryDeltas.Δtz
    Δφr = orbit.TrajectoryDeltas.Δφr; Δφz = orbit.TrajectoryDeltas.Δφz
    ω = m * Ωφ + n * Ωr + k * Ωz
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, a * ω)
    ν = RenormalizedAngularMomentum(s, l, m, a, ω, λ)
    R = TeukolskyRadial(s, l, m, a, ω; lambda=λ, nu=ν, Nmax=Nmax)
    α = _alpha_circ(a, m, ω, λ)
    # step counts (32·⌈|ωΔt'−mΔφ'+idx|⌉, ≥32) — spectral midpoint sums converge well above threshold
    stepsr = max(32 * ceil(Int, abs(ω * _dδ(Δtr, pi)  - m * _dδ(Δφr, pi)  + n)),
                 32 * ceil(Int, abs(ω * _dδ(Δtr, 0.0) - m * _dδ(Δφr, 0.0) + n)), 32)
    stepsz = max(32 * ceil(Int, abs(ω * _dδ(Δtz, pi/4) - m * _dδ(Δφz, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtz, 0.0)  - m * _dδ(Δφz, 0.0)  + k)), 32)
    npr = stepsr ÷ 2; npz = stepsz ÷ 2
    # cache the z-only quantities over the wz grid (built on the first ir pass)
    zc = Vector{NamedTuple}(undef, npz)
    sumP = zero(ComplexF64); sumM = zero(ComplexF64); W = zero(ComplexF64)
    for ir in 1:npr
        wr = (ir - 1//2) * 2pi / stepsr
        rp = rfun(wr)
        Δ = rp^2 - 2rp + a^2
        sqrtR = sqrt(max(((rp^2 + a^2) * EnT - a * LzT)^2 - Δ * (rp^2 + KT), 0.0))
        expr = exp(im * (ω * Δtr(wr) - m * Δφr(wr) + 2pi * n * (ir - 1//2) / stepsr))
        Din = _incl_radialD(R.In, rp, a, m, ω, λ); Dup = _incl_radialD(R.Up, rp, a, m, ω, λ)
        ir == 1 && (W = (Din.R * Dup.dR - Din.dR * Dup.R) / Δ)
        vn_p = -((rp^2 + a^2) * EnT - a * LzT + sqrtR) / Δ   # Ur = +√R
        vn_m = -((rp^2 + a^2) * EnT - a * LzT - sqrtR) / Δ   # Ur = −√R
        for iz in 1:npz
            if ir == 1
                wz = (iz - 1//2) * 2pi / stepsz
                zp = zfun(wz); sinθp = sqrt(1 - zp^2)
                sqrtZ = sqrt(max(-((1 - zp^2) * a * EnT - LzT)^2 + (1 - zp^2) * (KT - a^2 * zp^2), 0.0))
                expz = exp(im * (ω * Δtz(wz) - m * Δφz(wz) + 2pi * k * (iz - 1//2) / stepsz))
                ang = _angular_generic(l, m, a, ω, zp; Nmax=Nmax)
                zc[iz] = (zp=zp, sinθp=sinθp, sqrtZ=sqrtZ, expz=expz, S=ang.S, L2S=ang.L2S, L1L2S=ang.L1L2S)
            end
            zq = zc[iz]; zp = zq.zp; sinθp = zq.sinθp; expz = zq.expz
            ζ = rp - im * a * zp; ζbar = rp + im * a * zp
            vmb_m = (-im * (a * sinθp^2 * EnT - LzT) - zq.sqrtZ) / sinθp^2   # Uz = −√Z
            vmb_p = (-im * (a * sinθp^2 * EnT - LzT) + zq.sqrtZ) / sinθp^2   # Uz = +√Z
            FIn = _Fab(ζ, ζbar, a, sinθp, Din.R, Din.DR, Din.DDR, zq.S, zq.L2S, zq.L1L2S)
            FUp = _Fab(ζ, ζbar, a, sinθp, Dup.R, Dup.DR, Dup.DDR, zq.S, zq.L2S, zq.L1L2S)
            integ(vn, vmb, F) = vn^2 * F.Fnn + vn * vmb * F.Fnmb + vmb^2 * F.Fmbmb +
                                spar / sKT * (vn * F.Gn + vmb * F.Gmb)
            # 4 quadrants: (Ur,Uz,expr^±,expz^±) = (+,−,+,+),(−,−,−,+),(+,+,+,−),(−,+,−,−)
            sumP += integ(vn_p, vmb_m, FIn) * expr * expz + integ(vn_m, vmb_m, FIn) / expr * expz +
                    integ(vn_p, vmb_p, FIn) * expr / expz + integ(vn_m, vmb_p, FIn) / expr / expz
            sumM += integ(vn_p, vmb_m, FUp) * expr * expz + integ(vn_m, vmb_m, FUp) / expr * expz +
                    integ(vn_p, vmb_p, FUp) * expr / expz + integ(vn_m, vmb_p, FUp) / expr / expz
        end
    end
    pre = (1 + spar * EnT / (2 * sKT)) * 2pi / (Υt * W * stepsr * stepsz)
    CPlus = -pre * sumP; CMinus = -pre * sumM
    # Carter-K flux weight (.wl TeukolskySpinModeAnalytical CarterConstantK, L2006)
    cK = (1 + spar * EnT / (2 * sKT)) * (k * Υz - ω * Υtz + m * Υφz) + a * spar / (2 * sKT) * (a * ω - m)
    return (l=l, m=m, n=n, k=k, omega=ω, α=α, CPlus=CPlus, CMinus=CMinus,
            EnergyI = abs2(CPlus) / (4pi * ω^2), EnergyH = α * abs2(CMinus) / (4pi * ω^2),
            AngMomI = abs2(CPlus) * m / (4pi * ω^3), AngMomH = α * abs2(CMinus) * m / (4pi * ω^3),
            CarterKI = cK * abs2(CPlus) / (2pi * ω^3), CarterKH = cK * α * abs2(CMinus) / (2pi * ω^3),
            stepsr=stepsr, stepsz=stepsz)
end

"radial bundle at ω': λ, ν, R (TeukolskyRadial). Angular is z-dependent → built per wz."
function _incl_bundle(l, m, a, ω; Nmax)
    s = -2
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, a * ω)
    ν = RenormalizedAngularMomentum(s, l, m, a, ω, λ)
    R = TeukolskyRadial(s, l, m, a, ω; lambda=λ, nu=ν, Nmax=Nmax)
    return (λ=λ, R=R)
end

"""
    inclined_spin_flux(l, m, n, k, a, p, e, x; Nmax=80, h=3e-4) -> NamedTuple

Linear-in-spin GW flux correction (dE¹,dLz¹ at ℐ/ℋ) for a GENERIC (inclined) s=−2 point particle,
polar/radial harmonics (n,k), per `TeukolskySpinModeCorrectionAnalytical` (.wl L2864-3090). wr×wz
midpoint double-sum over 4 velocity quadrants; ∂_ωF/∂_ωW/∂_ωα via central FD. Returns
`(; dE1I,dE1H,dL1I,dL1H, CPlus0,CMinus0,CPlus1,CMinus1, omega, stepsr,stepsz)`. The DH P5 path.
"""
function inclined_spin_flux(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real; Nmax::Integer=80, h::Real=3e-4)
    En0, Lz0, Q0 = constants_of_motion(a, p, e, x)
    K0 = Q0 + (Lz0 - a * En0)^2; sK = sqrt(K0)
    fr = kerr_frequencies(a, p, e, x; Time="Mino")
    pt = geodesic_proper_time(a, p, e, x)
    Υt = fr.Ut; Υτ = pt.ProperTimeFrequency
    Ωr = fr.Ur / Υt; Ωz = fr.Uth / Υt; Ωφ = fr.Uph / Υt
    ω = m * Ωφ + n * Ωr + k * Ωz
    ω1 = 3 * Υτ * ω / (2 * sK * Υt); Υt1 = -3 * Υτ / (2 * sK)
    rad = geodesic_trajectory_radial(a, p, e, x); pol = geodesic_trajectory_polar(a, p, e, x)
    Δtr = qr -> rad(qr).Δtr; Δφr = qr -> rad(qr).Δφr
    Δtz = qz -> pol(qz).Δtz; Δφz = qz -> pol(qz).Δφz
    stepsr = max(32 * ceil(Int, abs(ω * _dδ(Δtr, pi)  - m * _dδ(Δφr, pi)  + n)),
                 32 * ceil(Int, abs(ω * _dδ(Δtr, 0.0) - m * _dδ(Δφr, 0.0) + n)), 32)
    stepsz = max(32 * ceil(Int, abs(ω * _dδ(Δtz, pi/4) - m * _dδ(Δφz, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtz, 0.0)  - m * _dδ(Δφz, 0.0)  + k)), 32)
    npr = stepsr ÷ 2; npz = stepsz ÷ 2
    b0 = _incl_bundle(l, m, a, ω; Nmax=Nmax)
    bp = _incl_bundle(l, m, a, ω + h; Nmax=Nmax); bm = _incl_bundle(l, m, a, ω - h; Nmax=Nmax)
    # z-grid cache (zp, √Z, expz, base_z phase, angular at ω/ω±h) — built on the first ir
    zc = Vector{NamedTuple}(undef, npz)
    sumP0 = sumP1 = sumM0 = sumM1 = zero(ComplexF64); W = W1 = zero(ComplexF64)
    integ0(vn, vmb, F) = vn^2 * F.Fnn + vn * vmb * F.Fnmb + vmb^2 * F.Fmbmb
    for ir in 1:npr
        wr = (ir - 1//2) * 2pi / stepsr
        g = rad(wr); rp = g.r; Δ = rp^2 - 2rp + a^2
        sqrtR = sqrt(max(((rp^2 + a^2) * En0 - a * Lz0)^2 - Δ * (rp^2 + K0), 0.0))
        expr = exp(im * (ω * g.Δtr - m * g.Δφr + 2pi * n * (ir - 1//2) / stepsr))
        baser = 3 * im * ω / (2 * sK) * (Υτ / Υt * g.Δtr - pt.Δτr(wr))   # radial linear phase (×Ur sign)
        Din0 = _incl_radialD(b0.R.In, rp, a, m, ω, b0.λ); Dup0 = _incl_radialD(b0.R.Up, rp, a, m, ω, b0.λ)
        Dinp = _incl_radialD(bp.R.In, rp, a, m, ω + h, bp.λ); Dinm = _incl_radialD(bm.R.In, rp, a, m, ω - h, bm.λ)
        Dupp = _incl_radialD(bp.R.Up, rp, a, m, ω + h, bp.λ); Dupm = _incl_radialD(bm.R.Up, rp, a, m, ω - h, bm.λ)
        if ir == 1
            W  = (Din0.R * Dup0.dR - Din0.dR * Dup0.R) / Δ
            Wp = (Dinp.R * Dupp.dR - Dinp.dR * Dupp.R) / Δ
            Wm = (Dinm.R * Dupm.dR - Dinm.dR * Dupm.R) / Δ
            W1 = (Wp - Wm) / (2h) * ω1
        end
        vn_p = -((rp^2 + a^2) * En0 - a * Lz0 + sqrtR) / Δ   # Ur=+√R
        vn_m = -((rp^2 + a^2) * En0 - a * Lz0 - sqrtR) / Δ   # Ur=−√R
        for iz in 1:npz
            if ir == 1
                wz = (iz - 1//2) * 2pi / stepsz
                h2 = pol(wz); zp = h2.z; sinθp = sqrt(1 - zp^2)
                sqrtZ = sqrt(max(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (K0 - a^2 * zp^2), 0.0))
                expz = exp(im * (ω * h2.Δtz - m * h2.Δφz + 2pi * k * (iz - 1//2) / stepsz))
                basez = 3 * im * ω / (2 * sK) * (Υτ / Υt * h2.Δtz - pt.Δτz(wz))   # polar linear phase (×Uz sign)
                a0 = _angular_generic(l, m, a, ω, zp; Nmax=Nmax)
                ap = _angular_generic(l, m, a, ω + h, zp; Nmax=Nmax)
                am = _angular_generic(l, m, a, ω - h, zp; Nmax=Nmax)
                zc[iz] = (zp=zp, sinθp=sinθp, sqrtZ=sqrtZ, expz=expz, basez=basez, a0=a0, ap=ap, am=am)
            end
            zq = zc[iz]; zp = zq.zp; sinθp = zq.sinθp; expz = zq.expz; basez = zq.basez
            ζ = rp - im * a * zp; ζbar = rp + im * a * zp
            vmb_m = (-im * (a * sinθp^2 * En0 - Lz0) - zq.sqrtZ) / sinθp^2   # Uz=−√Z
            vmb_p = (-im * (a * sinθp^2 * En0 - Lz0) + zq.sqrtZ) / sinθp^2   # Uz=+√Z
            FIn = _Fab(ζ, ζbar, a, sinθp, Din0.R, Din0.DR, Din0.DDR, zq.a0.S, zq.a0.L2S, zq.a0.L1L2S)
            FUp = _Fab(ζ, ζbar, a, sinθp, Dup0.R, Dup0.DR, Dup0.DDR, zq.a0.S, zq.a0.L2S, zq.a0.L1L2S)
            # ∂_ω F components via central FD (radial + angular both at ω±h)
            FInp = _Fab(ζ, ζbar, a, sinθp, Dinp.R, Dinp.DR, Dinp.DDR, zq.ap.S, zq.ap.L2S, zq.ap.L1L2S)
            FInm = _Fab(ζ, ζbar, a, sinθp, Dinm.R, Dinm.DR, Dinm.DDR, zq.am.S, zq.am.L2S, zq.am.L1L2S)
            FUpp = _Fab(ζ, ζbar, a, sinθp, Dupp.R, Dupp.DR, Dupp.DDR, zq.ap.S, zq.ap.L2S, zq.ap.L1L2S)
            FUpm = _Fab(ζ, ζbar, a, sinθp, Dupm.R, Dupm.DR, Dupm.DDR, zq.am.S, zq.am.L2S, zq.am.L1L2S)
            dFIn = (Fnn=(FInp.Fnn-FInm.Fnn)/(2h), Fnmb=(FInp.Fnmb-FInm.Fnmb)/(2h), Fmbmb=(FInp.Fmbmb-FInm.Fmbmb)/(2h))
            dFUp = (Fnn=(FUpp.Fnn-FUpm.Fnn)/(2h), Fnmb=(FUpp.Fnmb-FUpm.Fnmb)/(2h), Fmbmb=(FUpp.Fmbmb-FUpm.Fmbmb)/(2h))
            G(vn, vmb, F) = (vn * F.Gn + vmb * F.Gmb) / sK
            dF0(vn, vmb, dF) = vn^2 * dF.Fnn + vn * vmb * dF.Fnmb + vmb^2 * dF.Fmbmb
            # 4 quadrants: (Ur,Uz, expr^±, expz^±, expr1 sign{1,-1,1,-1}, expz1 sign{1,1,-1,-1})
            for (vn, vmb, ef, p1) in (
                    (vn_p, vmb_m, expr * expz,        baser + basez),
                    (vn_m, vmb_m, expz / expr,       -baser + basez),
                    (vn_p, vmb_p, expr / expz,        baser - basez),
                    (vn_m, vmb_p, 1 / (expr * expz), -baser - basez))
                sumP0 += integ0(vn, vmb, FIn) * ef
                sumM0 += integ0(vn, vmb, FUp) * ef
                sumP1 += (G(vn, vmb, FIn) + integ0(vn, vmb, FIn) * p1 + dF0(vn, vmb, dFIn) * ω1) * ef
                sumM1 += (G(vn, vmb, FUp) + integ0(vn, vmb, FUp) * p1 + dF0(vn, vmb, dFUp) * ω1) * ef
            end
        end
    end
    den = Υt * W * stepsr * stepsz
    CPlus0  = -2pi * sumP0 / den; CMinus0 = -2pi * sumM0 / den
    pref = En0 / (2 * sK) - Υt1 / Υt - W1 / W
    CPlus1  = -2pi * sumP1 / den + pref * CPlus0
    CMinus1 = -2pi * sumM1 / den + pref * CMinus0
    α  = _alpha_circ(a, m, ω, b0.λ)
    dα = (_alpha_circ(a, m, ω + h, bp.λ) - _alpha_circ(a, m, ω - h, bm.λ)) / (2h)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (dα * ω1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (dα * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    return (dE1I=dE1I, dE1H=dE1H, dL1I=dL1I, dL1H=dL1H, CPlus0=CPlus0, CMinus0=CMinus0,
            CPlus1=CPlus1, CMinus1=CMinus1, omega=ω, stepsr=stepsr, stepsz=stepsz)
end

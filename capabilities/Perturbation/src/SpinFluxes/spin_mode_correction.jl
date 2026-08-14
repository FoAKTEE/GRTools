# Numeric-orbit-correction linear-in-spin flux CORRECTION (geodesic mode + dE¹/dLz¹) —
# TeukolskySpinFluxes.wl TeukolskySpinModeCorrection (L2609-2860). Same wr×wθ Mino double-sum as
# TeukolskySpinModeFromCorrection, but here the s⁰ geodesic monopole (sumPlus0/sumMinus0) and its
# linear-in-spin correction (sumPlus1/sumMinus1) are kept as SEPARATE accumulators. sumPlus1 = two
# pieces: (1) the explicit-s¹ source correction (Ann0S = (Σ1/Σ)·vn²+2·vn1·vn+Adnn, AnntφS/AnnrS/…
# spin-dipole, all at the GEODESIC point), and (2) the ω-derivative ω1·∂_ω[geodesic-monopole-
# integrand] (the .wl RInrp1/fnn01 terms) via central FD over ω (fabi+radials at ω±h). W1, Γ1, α1
# enter the final assembly. The fab^(i) basis, angular layer, h1/h2/h3, _alpha_circ and _dδ are
# SHARED with from_correction.jl (same module). FluxesCorrection = dE¹/dLz¹ at ℐ/ℋ. R3.
#
# vs from_correction.jl: that path linearises via `(Σ+s·Σ1)` × `(vn+s·vn1)²` and one accumulator;
# here the .wl writes the s⁰ and s¹ pieces explicitly (Ann0S has the explicit monopole-correction
# (Σ1/Σ)vn²+2vn1·vn). The geodesic monopole (sumPlus0/CPlus0/Fluxes) equals FromCorrection's s=0.

"""
    _spin_mode_correction(oc, l, m, n, k; Nmax=80, h=3e-4) -> NamedTuple

Core port of TeukolskySpinModeCorrection (.wl L2609-2860). `oc` is the SpinningOrbit
`KerrSpinOrbitCorrection` NamedTuple. The .wl `{angparNew,RCorrection}` ω-derivative callbacks are
replaced by internal central FD over ω (`h`). Returns `(; l,m,n,k,omega, CPlus0,CMinus0,
CPlus1,CMinus1, α, Spi2, EnergyI,EnergyH,AngMomI,AngMomH, dE1I,dE1H,dL1I,dL1H, stepsr,stepsθ)`. R3.
"""
function _spin_mode_correction(oc, l::Integer, m::Integer, n::Integer, k::Integer;
        Nmax::Integer=80, h::Real=3e-4)
    (l < 2 || abs(m) > l) && error("TeukolskySpinModeCorrection: l<2 or |m|>l")
    a = oc.a; p = oc.p; e = oc.e
    En0 = oc.Ehat; Lz0 = oc.Lzhat; Kc0 = oc.Khat        # geodesic constants of motion
    En1 = oc.ES; Lz1 = oc.LS                             # linear corrections to constants of motion
    correction = oc.OrbitCorrection                      # (wr,wz) -> (ΔtS,rS,zS,ΔφS,UrS,UzS)
    Ωr, Ωθ, Ωφ = oc.BLFrequenciesGeo                     # geodesic BL frequencies
    Ωr1, Ωθ1, Ωφ1 = oc.BLFrequenciesCorrection           # linear corrections to the BL frequencies
    rfun = oc.TrajectoryGeo[1]; θfun = oc.TrajectoryGeo[2]
    rgeo(wr) = rfun(wr).r
    zgeo(wz) = cos(θfun(wz).θ)
    Γ  = oc.MinoFrequenciesGeo[1]
    Γ1 = oc.MinoFrequenciesCorrection[1]
    ω  = m * Ωφ + n * Ωr + k * Ωθ
    isfinite(ω) || error("TeukolskySpinModeCorrection: non-real ω")
    ω1 = m * Ωφ1 + n * Ωr1 + k * Ωθ1                     # linear correction to the mode frequency

    # geodesic + ω±h radial/angular bundles (FD over ω captures ∂_ωR / ∂_ωfab / ∂_ωddR / ∂_ωα).
    SWSH = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω)
    λ = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ω)
    ν = RenormalizedAngularMomentum(-2, l, m, a, ω, λ)
    R = TeukolskyRadial(-2, l, m, a, ω; lambda = λ, nu = ν, Nmax = Nmax)
    SWSHp = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * (ω + h)); λp = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * (ω + h))
    SWSHm = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * (ω - h)); λm = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * (ω - h))
    Rp = TeukolskyRadial(-2, l, m, a, ω + h; lambda = λp, nu = RenormalizedAngularMomentum(-2, l, m, a, ω + h, λp), Nmax = Nmax)
    Rm = TeukolskyRadial(-2, l, m, a, ω - h; lambda = λm, nu = RenormalizedAngularMomentum(-2, l, m, a, ω - h, λm), Nmax = Nmax)
    α  = _alpha_circ(a, m, ω, λ)
    α1 = (_alpha_circ(a, m, ω + h, λp) - _alpha_circ(a, m, ω - h, λm)) / (2h) * ω1   # α1 = ω1·dα/dω

    # step counts (.wl L2663-2666): 32·⌈|ωΔt'−mΔφ'+idx|⌉ at the turning points, ≥32.
    Δtr = oc.TrajectoryDeltas.Δtr; Δφr = oc.TrajectoryDeltas.Δφr
    Δtθ = oc.TrajectoryDeltas.Δtz; Δφθ = oc.TrajectoryDeltas.Δφz
    stepsr = max(32 * ceil(Int, abs(ω * _dδ(Δtr, pi)  - m * _dδ(Δφr, pi)  + n)),
                 32 * ceil(Int, abs(ω * _dδ(Δtr, 0.0) - m * _dδ(Δφr, 0.0) + n)), 32)
    stepsθ = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/4) - m * _dδ(Δφθ, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npr = stepsr ÷ 2; npθ = stepsθ ÷ 2

    # trajectory corrections at all (ir,iθ) midpoints (.wl `correctionp`); -iθ indexes from the end.
    correctionp = [[correction((ir - 1//2) * 2pi / stepsr, (iθ - 1//2) * 2pi / stepsθ)
                    for iθ in 1:npθ] for ir in 1:npr]

    θc = Vector{NamedTuple}(undef, npθ)                  # θ-only quantities, cached on first ir pass
    sumP0 = sumM0 = sumP1 = sumM1 = zero(ComplexF64)
    # retained final-loop radial values + Δ for the invariant Wronskian and its ω-FD (W1).
    RInrp = dRInrp = RUprp = dRUprp = zero(ComplexF64)
    RInpr = dRInpr = RUppr = dRUppr = zero(ComplexF64)   # In/Up R,dR at ω+h
    RInmr = dRInmr = RUpmr = dRUpmr = zero(ComplexF64)   # In/Up R,dR at ω−h
    rp_last = 0.0

    # geodesic-monopole-integrand piece-0 (ω-dependent factor only): the bracket multiplying exprq·expθq
    _piece0(Amnn, Amnmb, Ammbmb, f, RI, dRI, ddRI) =
        (Amnn * f[1] + Amnmb * f[2] + Ammbmb * f[4]) * RI -
        (Amnmb * f[3] + Ammbmb * f[5]) * dRI + Ammbmb * f[6] * ddRI

    for ir in 1:npr
        wr = (ir - 1//2) * 2pi / stepsr
        rp = rgeo(wr); rp_last = rp
        sqrtR = sqrt(((rp^2 + a^2) * En0 - a * Lz0)^2 - (rp^2 - 2rp + a^2) * (rp^2 + Kc0))
        Ur = (1, -1, 1, -1) .* sqrtR                                  # geodesic radial velocity
        expr = exp(im * (ω * Δtr(wr) - m * Δφr(wr) + 2pi * n * (ir - 1//2) / stepsr))
        exprq = (expr, 1 / expr, expr, 1 / expr)                       # expr^{1,-1,1,-1}
        Δ  = rp^2 - 2rp + a^2
        K  = (rp^2 + a^2) * ω - a * m
        dΔ = 2 * (rp - 1)
        dK = 2rp * ω
        V  = -(K^2 + 4im * (rp - 1) * K) / Δ + 8im * ω * rp + λ
        dV = -((2K * dK + 4im * K + 4im * (rp - 1) * dK) * Δ - (K^2 + 4im * (rp - 1) * K) * dΔ) / Δ^2 + 8im * ω
        # geodesic radials + 2nd/3rd derivatives from the Teukolsky equation (.wl L2684-2697)
        RInrp = R.In(rp); dRInrp = teukolsky_radial_deriv(R.In, rp)
        ddRInrp  = (V * RInrp + dΔ * dRInrp) / Δ
        dddRInrp = 1 / Δ * (dV * RInrp + (V + 2) * dRInrp)
        RUprp = R.Up(rp); dRUprp = teukolsky_radial_deriv(R.Up, rp)
        ddRUprp  = (V * RUprp + dΔ * dRUprp) / Δ
        dddRUprp = 1 / Δ * (dV * RUprp + (V + 2) * dRUprp)
        # ω±h radials (for the ω-FD of piece-0 and of W) — recompute V, ddR at the shifted ω.
        Kpω = (rp^2 + a^2) * (ω + h) - a * m; dKpω = 2rp * (ω + h)
        Vp  = -(Kpω^2 + 4im * (rp - 1) * Kpω) / Δ + 8im * (ω + h) * rp + λp
        Kmω = (rp^2 + a^2) * (ω - h) - a * m; dKmω = 2rp * (ω - h)
        Vm  = -(Kmω^2 + 4im * (rp - 1) * Kmω) / Δ + 8im * (ω - h) * rp + λm
        RInpr = Rp.In(rp); dRInpr = teukolsky_radial_deriv(Rp.In, rp); ddRInpr = (Vp * RInpr + dΔ * dRInpr) / Δ
        RInmr = Rm.In(rp); dRInmr = teukolsky_radial_deriv(Rm.In, rp); ddRInmr = (Vm * RInmr + dΔ * dRInmr) / Δ
        RUppr = Rp.Up(rp); dRUppr = teukolsky_radial_deriv(Rp.Up, rp); ddRUppr = (Vp * RUppr + dΔ * dRUppr) / Δ
        RUpmr = Rm.Up(rp); dRUpmr = teukolsky_radial_deriv(Rm.Up, rp); ddRUpmr = (Vm * RUpmr + dΔ * dRUpmr) / Δ

        for iθ in 1:npθ
            if ir == 1
                wθ = (iθ - 1//2) * 2pi / stepsθ
                zp = zgeo(wθ); sinθp = sqrt(1 - zp^2)
                sqrtZ = sqrt(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (Kc0 - a^2 * zp^2))
                Uz = (1, 1, -1, -1) .* (-sqrtZ)
                expθ = exp(im * (ω * Δtθ(wθ) - m * Δφθ(wθ) + 2pi * k * (iθ - 1//2) / stepsθ))
                expθq = (expθ, expθ, 1 / expθ, 1 / expθ)               # expθ^{1,1,-1,-1}
                ang  = _angular_fromcorr(SWSH,  λ,  a, m, ω,     zp)   # geodesic angular layer
                angp = _angular_fromcorr(SWSHp, λp, a, m, ω + h, zp)   # ω+h (for the piece-0 ω-FD)
                angm = _angular_fromcorr(SWSHm, λm, a, m, ω - h, zp)   # ω−h
                θc[iθ] = (zp = zp, sinθp = sinθp, Uz = Uz, expθq = expθq, Δtθwθ = Δtθ(wθ),
                          ang = ang, angp = angp, angm = angm)
            end
            t = θc[iθ]
            zp = t.zp; sinθp = t.sinθp; Uz = t.Uz; expθq = t.expθq; Δtθwθ = t.Δtθwθ
            ang = t.ang; S = ang.S; L2S = ang.L2S; L1L2S = ang.L1L2S
            dSdθ = ang.dSdθ; dL2Sdθ = ang.dL2Sdθ; dL1L2Sdθ = ang.dL1L2Sdθ

            ζ = rp - im * a * zp; ζbar = rp + im * a * zp; Σ = rp^2 + a^2 * zp^2
            (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2) =
                _fabi(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a)
            (dfnn0dr, dfnmb0dr, dfnmb1dr, dfmbmb0dr, dfmbmb1dr, dfmbmb2dr) =
                _dfabidr(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a, ω)
            (dfnn0dθ, dfnmb0dθ, dfnmb1dθ, dfmbmb0dθ, dfmbmb1dθ, dfmbmb2dθ) =
                _dfabidθ(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, dSdθ, dL2Sdθ, dL1L2Sdθ, a)

            # geodesic four-velocity in Kinnersley tetrad (4-vectors over quadrants)
            vn  = ((-1) .* ((rp^2 + a^2) * En0 - a * Lz0 .+ Ur)) ./ (2 * Σ)
            vl  = ((-1) .* ((rp^2 + a^2) * En0 - a * Lz0 .- Ur)) ./ Δ
            vm  = (im * (a * sinθp^2 * En0 - Lz0) .+ Uz) ./ (-sqrt(2) * sinθp * ζbar)
            vmb = conj.(vm)
            Amnn   = vn .^ 2
            Amnmb  = vn .* vmb
            Ammbmb = vmb .^ 2

            # spin tensor in Kinnersley tetrad (GEODESIC point; scalars in rp,zp — no explicit s here)
            Sln  = -((rp * (Kc0 - a^2 * zp^2)) / (sqrt(Kc0) * Σ))
            Snm  = (ζ / sqrt(Kc0)) .* vm .* vn
            Snmb = conj.(Snm)
            Slmb = (-(ζ / sqrt(Kc0))) .* vl .* vmb
            Smmb = (im * a * zp * (Kc0 + rp^2)) / (sqrt(Kc0) * Σ)

            # Newman-Penrose spin coefficients (opposite sign convention vs Teukolsky paper)
            ρ  = 1 / ζ
            βc = -(zp / (2 * ζbar * sqrt(2) * sinθp))
            πc = -((im * a * sinθp) / (ζ^2 * sqrt(2)))
            τc = (im * a * sinθp) / (sqrt(2) * Σ)
            μc = Δ / (2 * ζ^2 * ζbar)
            γc = (a^2 - rp + im * a * (-1 + rp) * zp) / (2 * ζ^2 * ζbar)
            αc = -((-rp * zp - im * a * (-2 + zp^2)) / (2 * ζ^2 * sqrt(2) * sinθp))

            Scdγndc  = -Sln * 2 * real(γc) - Smmb * (-μc + conj(μc))
            Scdγmbdc = -Sln * (πc + conj(τc)) .- Snmb .* conj(ρ) .- Slmb .* (-conj(γc) + γc - μc) .- Smmb * (-αc + conj(βc))
            Adnn   = (Scdγndc .* vn .- Sln * 2 * real(γc) .* vn .- 2 .* real.(Snmb .* ((conj(αc) + βc) .* vn .- μc .* vm)))
            Admbmb = (Scdγmbdc .* vmb .- Snmb .* ((-πc) .* vl) .- Slmb .* (conj(τc) .* vn .- (conj(γc) - γc) .* vmb) .+ Smmb .* ((αc - conj(βc)) .* vmb))
            Adnmb  = (Scdγndc .* vmb .+ Scdγmbdc .* vn .- Sln .* (conj(τc) .* vn .- (conj(γc) - γc) .* vmb) .-
                      Snmb .* (conj(ρ) .* vn .- μc .* vl .- (conj(αc) - βc) .* vmb) .-
                      Snm .* ((αc - conj(βc)) .* vmb) .- Snmb .* ((-conj(πc)) .* vmb .- πc .* vm) .-
                      Slmb .* (2 * real(γc) .* vn) .+ Smmb .* ((αc + conj(βc)) .* vn .- conj(μc) .* vmb)) ./ 2
            # functions from Eq. (B3) in 2303.16798 (geodesic K,Δ,Σ,ζ)
            Stφn  = -im * K / (2 * Σ) * Sln .+ (a * ω * sinθp - m / sinθp) / (sqrt(2) * Σ) .* (ζ .* Snmb .- ζbar .* Snm)
            Stφmb = -im * K .* (1 / Δ .* Snmb .+ 1 / (2 * Σ) .* Slmb) .+ (a * ω * sinθp - m / sinθp) / (sqrt(2) * ζ) .* Smmb
            Srn   = Δ / (2 * Σ) * Sln
            Srmb  = ((-1) .* Snmb) .+ Δ / (2 * Σ) .* Slmb
            Sθn   = ((-1) .* (Snmb .* ζ .+ Snm .* ζbar)) ./ (sqrt(2) * Σ)
            Sθmb  = Smmb / (sqrt(2) * ζ)

            # numeric trajectory correction at the 4 quadrants (.wl L2762-2765). cp/cm = +iθ / -iθ.
            cp = correctionp[ir][iθ]
            cm = correctionp[ir][npθ - iθ + 1]
            r1p  = (cp.rS,  cm.rS,  cm.rS,  cp.rS)
            z1p  = (cp.zS,  -cm.zS, -cm.zS, cp.zS)
            Ur1p = (cp.UrS, -cm.UrS, cm.UrS, -cp.UrS)
            Uz1p = (cp.UzS, cm.UzS, -cm.UzS, -cp.UzS)
            Σ1   = 2 .* (rp .* r1p .+ a^2 * zp .* z1p)
            # linear parts of the Δt/Δφ phase (.wl L2767-2770): O(s) trajectory shift + ω1·(geodesic Δtr+Δtθ)
            Δtrwr = Δtr(wr)
            exp1 = (im * ((ω * cp.ΔtS - m * cp.ΔφS) + ω1 * Δtrwr + ω1 * Δtθwθ),
                    im * (-(ω * cm.ΔtS - m * cm.ΔφS) - ω1 * Δtrwr + ω1 * Δtθwθ),
                    im * ((ω * cm.ΔtS - m * cm.ΔφS) + ω1 * Δtrwr - ω1 * Δtθwθ),
                    im * (-(ω * cp.ΔtS - m * cp.ΔφS) - ω1 * Δtrwr - ω1 * Δtθwθ))

            h1v = _fc_h1(rp, zp, Kc0, a); h2v = _fc_h2(rp, zp, Kc0, En0, Lz0, a); h3v = _fc_h3(rp, zp, Kc0, a)
            # linear parts of the four-velocity in Kinnersley tetrad (.wl L2771-2774)
            vn1 = ((-1) .* (((2 * rp * En0) .* r1p) .- ((rp^2 + a^2) * En0 - a * Lz0 .+ Ur) ./ Σ .* Σ1 .+
                     ((rp^2 + a^2) * (En1 - h1v) .- a .* (Lz1 .+ h2v .+ h3v .* Ur .* Uz) .+ Ur1p))) ./ (2 * Σ)
            vmb1 = (((-im * (-2 * a * zp * En0)) .* z1p) .-
                     (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* z1p .+ (r1p .- im * a .* z1p) ./ ζ) .+
                     (-im .* (a * sinθp^2 * (En1 - h1v) .- (Lz1 .+ h2v .+ h3v .* Ur .* Uz)) .+ Uz1p)) ./ (-sqrt(2) * sinθp * ζ)

            # explicit-s¹ source amplitude pieces (.wl L2775-2786): monopole-correction + spin-dipole
            Ann0S    = ((Σ1 ./ Σ) .* vn .+ 2 .* vn1) .* vn .+ Adnn
            AnntφS   = (Stφn .+ exp1 .* vn) .* vn
            AnnrS    = (Srn .+ r1p .* vn) .* vn
            AnnθS    = (Sθn .- z1p ./ sinθp .* vn) .* vn
            Anmb0S   = ((Σ1 ./ Σ) .* vn .* vmb .+ vn1 .* vmb .+ vn .* vmb1 .+ Adnmb)
            AnmbtφS  = ((Stφn .* vmb .+ Stφmb .* vn) ./ 2 .+ exp1 .* vn .* vmb)
            AnmbrS   = ((Srn .* vmb .+ Srmb .* vn) ./ 2 .+ r1p .* vn .* vmb)
            AnmbθS   = ((Sθn .* vmb .+ Sθmb .* vn) ./ 2 .- z1p ./ sinθp .* vn .* vmb)
            Ambmb0S  = ((Σ1 ./ Σ) .* vmb .+ 2 .* vmb1) .* vmb .+ Admbmb
            AmbmbtφS = (Stφmb .+ exp1 .* vmb) .* vmb
            AmbmbrS  = (Srmb .+ r1p .* vmb) .* vmb
            AmbmbθS  = (Sθmb .- z1p ./ sinθp .* vmb) .* vmb

            # piece-2: ω1·∂_ω[geodesic-monopole-integrand]. Amnn/Amnmb/Ammbmb,Σ ω-independent ⇒ FD
            # the piece-0 bracket (fabi+radials+ddR at ω±h), × Σ, × ω1. (.wl RInrp1/fnn01 terms.)
            fp = _fabi(ζ, ζbar, sinθp, Δ, dΔ, Kpω, dKpω, t.angp.S, t.angp.L2S, t.angp.L1L2S, a)
            fm = _fabi(ζ, ζbar, sinθp, Δ, dΔ, Kmω, dKmω, t.angm.S, t.angm.L2S, t.angm.L1L2S, a)

            # quadrant sum
            for q in 1:4
                p0In = _piece0(Amnn[q], Amnmb[q], Ammbmb[q], (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2), RInrp, dRInrp, ddRInrp)
                p0Up = _piece0(Amnn[q], Amnmb[q], Ammbmb[q], (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2), RUprp, dRUprp, ddRUprp)
                phase = exprq[q] * expθq[q]
                sumP0 += Σ * p0In * phase
                sumM0 += Σ * p0Up * phase

                # part 1 (.wl L2790-2795 In, L2799-2804 Up): explicit-s¹ source × radials
                s1In = ((Ann0S[q] + AnntφS[q]) * RInrp * fnn0 + AnnrS[q] * (dRInrp * fnn0 + RInrp * dfnn0dr) + AnnθS[q] * RInrp * dfnn0dθ) +
                       ((Anmb0S[q] + AnmbtφS[q]) * RInrp * fnmb0 + AnmbrS[q] * (RInrp * dfnmb0dr + dRInrp * fnmb0) + AnmbθS[q] * RInrp * dfnmb0dθ) -
                       ((Anmb0S[q] + AnmbtφS[q]) * dRInrp * fnmb1 + AnmbrS[q] * (dRInrp * dfnmb1dr + ddRInrp * fnmb1) + AnmbθS[q] * dRInrp * dfnmb1dθ) +
                       ((Ambmb0S[q] + AmbmbtφS[q]) * RInrp * fmbmb0 + AmbmbrS[q] * (dRInrp * fmbmb0 + RInrp * dfmbmb0dr) + AmbmbθS[q] * RInrp * dfmbmb0dθ) -
                       ((Ambmb0S[q] + AmbmbtφS[q]) * dRInrp * fmbmb1 + AmbmbrS[q] * (ddRInrp * fmbmb1 + dRInrp * dfmbmb1dr) + AmbmbθS[q] * dRInrp * dfmbmb1dθ) +
                       ((Ambmb0S[q] + AmbmbtφS[q]) * ddRInrp * fmbmb2 + AmbmbrS[q] * (dddRInrp * fmbmb2 + ddRInrp * dfmbmb2dr) + AmbmbθS[q] * ddRInrp * dfmbmb2dθ)
                s1Up = ((Ann0S[q] + AnntφS[q]) * RUprp * fnn0 + AnnrS[q] * (dRUprp * fnn0 + RUprp * dfnn0dr) + AnnθS[q] * RUprp * dfnn0dθ) +
                       ((Anmb0S[q] + AnmbtφS[q]) * RUprp * fnmb0 + AnmbrS[q] * (RUprp * dfnmb0dr + dRUprp * fnmb0) + AnmbθS[q] * RUprp * dfnmb0dθ) -
                       ((Anmb0S[q] + AnmbtφS[q]) * dRUprp * fnmb1 + AnmbrS[q] * (dRUprp * dfnmb1dr + ddRUprp * fnmb1) + AnmbθS[q] * dRUprp * dfnmb1dθ) +
                       ((Ambmb0S[q] + AmbmbtφS[q]) * RUprp * fmbmb0 + AmbmbrS[q] * (dRUprp * fmbmb0 + RUprp * dfmbmb0dr) + AmbmbθS[q] * RUprp * dfmbmb0dθ) -
                       ((Ambmb0S[q] + AmbmbtφS[q]) * dRUprp * fmbmb1 + AmbmbrS[q] * (ddRUprp * fmbmb1 + dRUprp * dfmbmb1dr) + AmbmbθS[q] * dRUprp * dfmbmb1dθ) +
                       ((Ambmb0S[q] + AmbmbtφS[q]) * ddRUprp * fmbmb2 + AmbmbrS[q] * (dddRUprp * fmbmb2 + ddRUprp * dfmbmb2dr) + AmbmbθS[q] * ddRUprp * dfmbmb2dθ)

                # part 2 (.wl L2796-2797 In, L2805-2806 Up): ω1·∂_ω[piece-0] via central FD
                p0Inp = _piece0(Amnn[q], Amnmb[q], Ammbmb[q], fp, RInpr, dRInpr, ddRInpr)
                p0Inm = _piece0(Amnn[q], Amnmb[q], Ammbmb[q], fm, RInmr, dRInmr, ddRInmr)
                p0Upp = _piece0(Amnn[q], Amnmb[q], Ammbmb[q], fp, RUppr, dRUppr, ddRUppr)
                p0Upm = _piece0(Amnn[q], Amnmb[q], Ammbmb[q], fm, RUpmr, dRUpmr, ddRUpmr)
                p2In = (p0Inp - p0Inm) / (2h) * ω1
                p2Up = (p0Upp - p0Upm) / (2h) * ω1

                sumP1 += Σ * s1In * phase + Σ * p2In * phase
                sumM1 += Σ * s1Up * phase + Σ * p2Up * phase
            end
        end
    end

    Δlast = rp_last^2 - 2 * rp_last + a^2
    W  = (RInrp * dRUprp - dRInrp * RUprp) / Δlast               # invariant Wronskian (.wl L2809)
    Wp = (RInpr * dRUppr - dRInpr * RUppr) / Δlast
    Wm = (RInmr * dRUpmr - dRInmr * RUpmr) / Δlast
    W1 = (Wp - Wm) / (2h) * ω1                                   # W1 = ω1·dW/dω (.wl L2810)
    den = Γ * W * stepsr * stepsθ
    CPlus0  = -2pi * sumP0 / den
    CMinus0 = -2pi * sumM0 / den
    CPlus1  = -2pi * (sumP1 - Γ1 / Γ * sumP0 - W1 / W * sumP0) / den   # (.wl L2813)
    CMinus1 = -2pi * (sumM1 - Γ1 / Γ * sumM0 - W1 / W * sumM0) / den

    # Fluxes (.wl L2833-2843) + FluxesCorrection dE¹/dLz¹ (.wl L2845-2855)
    EnergyI = abs2(CPlus0) / (4pi * ω^2);  EnergyH = α * abs2(CMinus0) / (4pi * ω^2)
    AngMomI = abs2(CPlus0) * m / (4pi * ω^3); AngMomH = α * abs2(CMinus0) * m / (4pi * ω^3)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (α1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (α1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    return (l = l, m = m, n = n, k = k, omega = ω,
            CPlus0 = CPlus0, CMinus0 = CMinus0, CPlus1 = CPlus1, CMinus1 = CMinus1,
            α = α, Spi2 = SWSH(pi / 2, 0.0),
            EnergyI = EnergyI, EnergyH = EnergyH, AngMomI = AngMomI, AngMomH = AngMomH,
            dE1I = dE1I, dE1H = dE1H, dL1I = dL1I, dL1H = dL1H,
            stepsr = stepsr, stepsθ = stepsθ)
end

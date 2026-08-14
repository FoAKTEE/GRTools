# Numeric-orbit-correction linear-in-spin Teukolsky mode — faithful 1:1 of TeukolskySpinFluxes/src/
# from_correction.jl (TeukolskySpinFluxes.wl TeukolskySpinModeFromCorrection L1656-1872), the
# "DH-trajectory" flux path. The shared fab^(i)/angular/_fc_h*/dJz*/_KCorrection basis lives in
# fab_basis.jl; the local KerrSpinOrbitCorrection kinematics in spinning_orbit_correction_local.jl.
"""
    _from_correction_mode(oc, l, m, n, k, s) -> NamedTuple

Core port of TeukolskySpinModeFromCorrection (.wl L1656-1872). `oc` is the SpinningOrbit
`KerrSpinOrbitCorrection` NamedTuple. Returns `(; l,m,n,k,omega, CPlus, CMinus, α, Spi2,
EnergyI, EnergyH, AngMomI, AngMomH, CarterKI, CarterKH, stepsr, stepsθ)`. R3.
"""
function _from_correction_mode(oc, l::Integer, m::Integer, n::Integer, k::Integer, s::Real)
    (l < 2 || abs(m) > l) && error("TeukolskySpinModeFromCorrection: l<2 or |m|>l")
    a = oc.a; p = oc.p; e = oc.e
    En0 = oc.Ehat; Lz0 = oc.Lzhat; Kc0 = oc.Khat       # geodesic constants of motion
    En1 = oc.ES; Lz1 = oc.LS                            # linear corrections to constants of motion
    Kc1 = _KCorrection(a, p, e, cos(oc.ScriptI), En0, Lz0, Kc0)   # linear-in-spin K correction (.wl L1679)
    correction = oc.OrbitCorrection                      # (wr,wz) -> (ΔtS,rS,zS,ΔφS,UrS,UzS)
    # BL frequencies: geo + s·correction
    Ωr, Ωθ, Ωφ = oc.BLFrequenciesGeo .+ s .* oc.BLFrequenciesCorrection
    rfun = oc.TrajectoryGeo[1]; θfun = oc.TrajectoryGeo[2]
    rgeo(wr) = rfun(wr).r
    zgeo(wz) = cos(θfun(wz).θ)
    Γ  = oc.MinoFrequenciesGeo[1]
    Γ1 = oc.MinoFrequenciesCorrection[1]
    ω = m * Ωφ + n * Ωr + k * Ωθ
    isfinite(ω) || error("TeukolskySpinModeFromCorrection: non-real ω")

    SWSH = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω)
    λ = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ω)
    ν = RenormalizedAngularMomentum(-2, l, m, a, ω, λ)
    R = TeukolskyRadial(-2, l, m, a, ω; lambda = λ, nu = ν)
    α = _alpha_circ(a, m, ω, λ)

    # step counts (.wl L1700-1703): 32·⌈|ωΔt'−mΔφ'+idx|⌉ at the turning points, ≥32.
    Δtr = oc.TrajectoryDeltas.Δtr; Δφr = oc.TrajectoryDeltas.Δφr
    Δtθ = oc.TrajectoryDeltas.Δtz; Δφθ = oc.TrajectoryDeltas.Δφz
    stepsr = max(32 * ceil(Int, abs(ω * _dδ(Δtr, pi)  - m * _dδ(Δφr, pi)  + n)),
                 32 * ceil(Int, abs(ω * _dδ(Δtr, 0.0) - m * _dδ(Δφr, 0.0) + n)), 32)
    stepsθ = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/4) - m * _dδ(Δφθ, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npr = stepsr ÷ 2; npθ = stepsθ ÷ 2

    # precompute the trajectory corrections at all (ir,iθ) midpoints (.wl `correctionp`).
    # correctionp[ir][iθ] = (ΔtS,rS,zS,ΔφS,UrS,UzS) at (wr_ir, wθ_iθ). The .wl indexes the columns
    # with ±iθ, where -iθ counts from the END of the (npθ-long) row → column (npθ-iθ+1).
    correctionp = [[correction((ir - 1//2) * 2pi / stepsr, (iθ - 1//2) * 2pi / stepsθ)
                    for iθ in 1:npθ] for ir in 1:npr]

    # cache the θ-only quantities over the wθ grid (built on the first ir pass)
    θc = Vector{NamedTuple}(undef, npθ)
    sumP = zero(ComplexF64); sumM = zero(ComplexF64)
    RInrp = zero(ComplexF64); dRInrp = zero(ComplexF64); RUprp = zero(ComplexF64); dRUprp = zero(ComplexF64)
    rp_last = 0.0  # last rp / Δ retained for the invariant Wronskian (.wl uses the final loop value)

    for ir in 1:npr
        wr = (ir - 1//2) * 2pi / stepsr
        rp = rgeo(wr); rp_last = rp
        sqrtR = sqrt(((rp^2 + a^2) * En0 - a * Lz0)^2 - (rp^2 - 2rp + a^2) * (rp^2 + Kc0))
        Ur = (1, -1, 1, -1) .* sqrtR                                  # 4 quadrants
        expr = exp(im * (ω * Δtr(wr) - m * Δφr(wr) + 2pi * n * (ir - 1//2) / stepsr))
        exprq = (expr, 1 / expr, expr, 1 / expr)                       # expr^{1,-1,1,-1}
        Δ  = rp^2 - 2rp + a^2
        K  = (rp^2 + a^2) * ω - a * m
        dΔ = 2 * (rp - 1)
        dK = 2rp * ω
        V  = -(K^2 + 4im * (rp - 1) * K) / Δ + 8im * ω * rp + λ
        dV = -((2K * dK + 4im * K + 4im * (rp - 1) * dK) * Δ - (K^2 + 4im * (rp - 1) * K) * dΔ) / Δ^2 + 8im * ω
        RInrp = R.In(rp); dRInrp = teukolsky_radial_deriv(R.In, rp)
        ddRInrp  = (V * RInrp + dΔ * dRInrp) / Δ
        dddRInrp = 1 / Δ * (dV * RInrp + (V + 2) * dRInrp)
        RUprp = R.Up(rp); dRUprp = teukolsky_radial_deriv(R.Up, rp)
        ddRUprp  = (V * RUprp + dΔ * dRUprp) / Δ
        dddRUprp = 1 / Δ * (dV * RUprp + (V + 2) * dRUprp)

        for iθ in 1:npθ
            if ir == 1
                wθ = (iθ - 1//2) * 2pi / stepsθ
                zp = zgeo(wθ); sinθp = sqrt(1 - zp^2)
                sqrtZ = sqrt(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (Kc0 - a^2 * zp^2))
                Uz = (1, 1, -1, -1) .* (-sqrtZ)                        # {1,1,-1,-1}*(-1)*√Z
                expθ = exp(im * (ω * Δtθ(wθ) - m * Δφθ(wθ) + 2pi * k * (iθ - 1//2) / stepsθ))
                expθq = (expθ, expθ, 1 / expθ, 1 / expθ)               # expθ^{1,1,-1,-1}
                ang = _angular_fromcorr(SWSH, λ, a, m, ω, zp)
                θc[iθ] = (zp = zp, sinθp = sinθp, Uz = Uz, expθq = expθq,
                          S = ang.S, dSdθ = ang.dSdθ, d2Sdθ2 = ang.d2Sdθ2, d3Sdθ3 = ang.d3Sdθ3,
                          L2S = ang.L2S, dL2Sdθ = ang.dL2Sdθ, L1L2S = ang.L1L2S, dL1L2Sdθ = ang.dL1L2Sdθ)
            end
            t = θc[iθ]
            zp = t.zp; sinθp = t.sinθp; Uz = t.Uz; expθq = t.expθq
            S = t.S; L2S = t.L2S; L1L2S = t.L1L2S; dSdθ = t.dSdθ; dL2Sdθ = t.dL2Sdθ; dL1L2Sdθ = t.dL1L2Sdθ

            ζ = rp - im * a * zp; ζbar = rp + im * a * zp; Σ = rp^2 + a^2 * zp^2
            (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2) =
                _fabi(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a)
            (dfnn0dr, dfnmb0dr, dfnmb1dr, dfmbmb0dr, dfmbmb1dr, dfmbmb2dr) =
                _dfabidr(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a, ω)
            (dfnn0dθ, dfnmb0dθ, dfnmb1dθ, dfmbmb0dθ, dfmbmb1dθ, dfmbmb2dθ) =
                _dfabidθ(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, dSdθ, dL2Sdθ, dL1L2Sdθ, a)

            # four-velocity in Kinnersley tetrad (4-vectors over quadrants)
            vn  = ((-1) .* ((rp^2 + a^2) * En0 - a * Lz0 .+ Ur)) ./ (2 * Σ)
            vl  = ((-1) .* ((rp^2 + a^2) * En0 - a * Lz0 .- Ur)) ./ Δ
            vm  = (im * (a * sinθp^2 * En0 - Lz0) .+ Uz) ./ (-sqrt(2) * sinθp * ζbar)
            vmb = conj.(vm)

            # spin tensor in Kinnersley tetrad
            Sln  = -((rp * (Kc0 - a^2 * zp^2)) / (sqrt(Kc0) * Σ))
            Snm  = (ζ / sqrt(Kc0)) .* vm .* vn
            Snmb = conj.(Snm)
            Slmb = (-(ζ / sqrt(Kc0))) .* vl .* vmb
            Smmb = (im * a * zp * (Kc0 + rp^2)) / (sqrt(Kc0) * Σ)

            # Newman-Penrose spin coefficients (opposite sign convention vs Teukolsky paper)
            ρ = 1 / ζ
            βc = -(zp / (2 * ζbar * sqrt(2) * sinθp))
            πc = -((im * a * sinθp) / (ζ^2 * sqrt(2)))
            τc = (im * a * sinθp) / (sqrt(2) * Σ)
            μc = Δ / (2 * ζ^2 * ζbar)
            γc = (a^2 - rp + im * a * (-1 + rp) * zp) / (2 * ζ^2 * ζbar)
            αc = -((-rp * zp - im * a * (-2 + zp^2)) / (2 * ζ^2 * sqrt(2) * sinθp))

            Scdγndc  = -Sln * 2 * real(γc) - Smmb * (-μc + conj(μc))
            Scdγmbdc = -Sln * (πc + conj(τc)) .- Snmb .* conj(ρ) .- Slmb .* (-conj(γc) + γc - μc) .- Smmb * (-αc + conj(βc))
            # dipole terms
            Adnn  = (Scdγndc .* vn .- Sln * 2 * real(γc) .* vn .- 2 .* real.(Snmb .* ((conj(αc) + βc) .* vn .- μc .* vm)))
            Admbmb = (Scdγmbdc .* vmb .- Snmb .* ((-πc) .* vl) .- Slmb .* (conj(τc) .* vn .- (conj(γc) - γc) .* vmb) .+ Smmb .* ((αc - conj(βc)) .* vmb))
            Adnmb = (Scdγndc .* vmb .+ Scdγmbdc .* vn .- Sln .* (conj(τc) .* vn .- (conj(γc) - γc) .* vmb) .-
                     Snmb .* (conj(ρ) .* vn .- μc .* vl .- (conj(αc) - βc) .* vmb) .-
                     Snm .* ((αc - conj(βc)) .* vmb) .- Snmb .* ((-conj(πc)) .* vmb .- πc .* vm) .-
                     Slmb .* (2 * real(γc) .* vn) .+ Smmb .* ((αc + conj(βc)) .* vn .- conj(μc) .* vmb)) ./ 2
            # functions from Eq. (B3) in 2303.16798
            Stφn  = -im * K / (2 * Σ) * Sln .+ (a * ω * sinθp - m / sinθp) / (sqrt(2) * Σ) .* (ζ .* Snmb .- ζbar .* Snm)
            Stφmb = -im * K .* (1 / Δ .* Snmb .+ 1 / (2 * Σ) .* Slmb) .+ (a * ω * sinθp - m / sinθp) / (sqrt(2) * ζ) .* Smmb
            Srn   = Δ / (2 * Σ) * Sln
            Srmb  = ((-1) .* Snmb) .+ Δ / (2 * Σ) .* Slmb
            Sθn   = ((-1) .* (Snmb .* ζ .+ Snm .* ζbar)) ./ (sqrt(2) * Σ)
            Sθmb  = Smmb / (sqrt(2) * ζ)

            # trajectory correction at the 4 quadrants (.wl L1780-1786). cp/cm = +iθ / -iθ columns.
            cp = correctionp[ir][iθ]
            cm = correctionp[ir][npθ - iθ + 1]
            r1p  = (cp.rS,  cm.rS,  cm.rS,  cp.rS)
            z1p  = (cp.zS,  -cm.zS, -cm.zS, cp.zS)
            Ur1p = (cp.UrS, -cm.UrS, cm.UrS, -cp.UrS)
            Uz1p = (cp.UzS, cm.UzS, -cm.UzS, -cp.UzS)
            Σ1   = 2 .* (rp .* r1p .+ a^2 * zp .* z1p)
            exp1 = (im * (ω * cp.ΔtS - m * cp.ΔφS), -im * (ω * cm.ΔtS - m * cm.ΔφS),
                    im * (ω * cm.ΔtS - m * cm.ΔφS), -im * (ω * cp.ΔtS - m * cp.ΔφS))

            h1v = _fc_h1(rp, zp, Kc0, a); h2v = _fc_h2(rp, zp, Kc0, En0, Lz0, a); h3v = _fc_h3(rp, zp, Kc0, a)
            # linear parts of the four-velocity in Kinnersley tetrad
            vn1 = ((-1) .* (((2 * rp * En0) .* r1p) .- ((rp^2 + a^2) * En0 - a * Lz0 .+ Ur) ./ Σ .* Σ1 .+
                     ((rp^2 + a^2) * (En1 - h1v) .- a .* (Lz1 .+ h2v .+ h3v .* Ur .* Uz) .+ Ur1p))) ./ (2 * Σ)
            vmb1 = (((-im * (-2 * a * zp * En0)) .* z1p) .-
                     (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* z1p .+ (r1p .- im * a .* z1p) ./ ζ) .+
                     (-im .* (a * sinθp^2 * (En1 - h1v) .- (Lz1 .+ h2v .+ h3v .* Ur .* Uz)) .+ Uz1p)) ./ (-sqrt(2) * sinθp * ζ)

            # monopole terms
            Amnn   = (vn .+ s .* vn1) .^ 2
            Amnmb  = (vn .+ s .* vn1) .* (vmb .+ s .* vmb1)
            Ammbmb = (vmb .+ s .* vmb1) .^ 2
            # spin-dipole amplitude pieces
            Ann0S    = Adnn
            AnntφS   = (Stφn .+ exp1 .* vn) .* vn
            AnnrS    = (Srn .+ r1p .* vn) .* vn
            AnnθS    = (Sθn .- z1p ./ sinθp .* vn) .* vn
            Anmb0S   = Adnmb
            AnmbtφS  = ((Stφn .* vmb .+ Stφmb .* vn) ./ 2 .+ exp1 .* vn .* vmb)
            AnmbrS   = ((Srn .* vmb .+ Srmb .* vn) ./ 2 .+ r1p .* vn .* vmb)
            AnmbθS   = ((Sθn .* vmb .+ Sθmb .* vn) ./ 2 .- z1p ./ sinθp .* vn .* vmb)
            Ambmb0S  = Admbmb
            AmbmbtφS = (Stφmb .+ exp1 .* vmb) .* vmb
            AmbmbrS  = (Srmb .+ r1p .* vmb) .* vmb
            AmbmbθS  = (Sθmb .- z1p ./ sinθp .* vmb) .* vmb

            # F_ab radial assembly (In)
            FnnIn = RInrp * fnn0
            FnmbIn = RInrp * fnmb0 - dRInrp * fnmb1
            FmbmbIn = RInrp * fmbmb0 - dRInrp * fmbmb1 + ddRInrp * fmbmb2
            dFnnIndr = dRInrp * fnn0 + RInrp * dfnn0dr
            dFnmbIndr = dRInrp * fnmb0 + RInrp * dfnmb0dr - ddRInrp * fnmb1 - dRInrp * dfnmb1dr
            dFmbmbIndr = dRInrp * fmbmb0 + RInrp * dfmbmb0dr - ddRInrp * fmbmb1 - dRInrp * dfmbmb1dr + dddRInrp * fmbmb2 + ddRInrp * dfmbmb2dr
            dFnnIndθ = RInrp * dfnn0dθ
            dFnmbIndθ = RInrp * dfnmb0dθ - dRInrp * dfnmb1dθ
            dFmbmbIndθ = RInrp * dfmbmb0dθ - dRInrp * dfmbmb1dθ + ddRInrp * dfmbmb2dθ
            integP = (Σ .+ s .* Σ1) .* (s .* (
                ((Ann0S .+ AnntφS) .* FnnIn .+ AnnrS .* dFnnIndr .+ AnnθS .* dFnnIndθ) .+
                ((Anmb0S .+ AnmbtφS) .* FnmbIn .+ AnmbrS .* dFnmbIndr .+ AnmbθS .* dFnmbIndθ) .+
                ((Ambmb0S .+ AmbmbtφS) .* FmbmbIn .+ AmbmbrS .* dFmbmbIndr .+ AmbmbθS .* dFmbmbIndθ)) .+
                (Amnn .* FnnIn .+ Amnmb .* FnmbIn .+ Ammbmb .* FmbmbIn)) .* exprq .* expθq
            sumP += sum(integP)

            # F_ab radial assembly (Up)
            FnnUp = RUprp * fnn0
            FnmbUp = RUprp * fnmb0 - dRUprp * fnmb1
            FmbmbUp = RUprp * fmbmb0 - dRUprp * fmbmb1 + ddRUprp * fmbmb2
            dFnnUpdr = dRUprp * fnn0 + RUprp * dfnn0dr
            dFnmbUpdr = dRUprp * fnmb0 + RUprp * dfnmb0dr - ddRUprp * fnmb1 - dRUprp * dfnmb1dr
            dFmbmbUpdr = dRUprp * fmbmb0 + RUprp * dfmbmb0dr - ddRUprp * fmbmb1 - dRUprp * dfmbmb1dr + dddRUprp * fmbmb2 + ddRUprp * dfmbmb2dr
            dFnnUpdθ = RUprp * dfnn0dθ
            dFnmbUpdθ = RUprp * dfnmb0dθ - dRUprp * dfnmb1dθ
            dFmbmbUpdθ = RUprp * dfmbmb0dθ - dRUprp * dfmbmb1dθ + ddRUprp * dfmbmb2dθ
            integM = (Σ .+ s .* Σ1) .* (s .* (
                ((Ann0S .+ AnntφS) .* FnnUp .+ AnnrS .* dFnnUpdr .+ AnnθS .* dFnnUpdθ) .+
                ((Anmb0S .+ AnmbtφS) .* FnmbUp .+ AnmbrS .* dFnmbUpdr .+ AnmbθS .* dFnmbUpdθ) .+
                ((Ambmb0S .+ AmbmbtφS) .* FmbmbUp .+ AmbmbrS .* dFmbmbUpdr .+ AmbmbθS .* dFmbmbUpdθ)) .+
                (Amnn .* FnnUp .+ Amnmb .* FnmbUp .+ Ammbmb .* FmbmbUp)) .* exprq .* expθq
            sumM += sum(integM)
        end
    end

    # invariant Wronskian (.wl L1836, evaluated at the last loop rp)
    W = (RInrp * dRUprp - dRInrp * RUprp) / (rp_last^2 - 2 * rp_last + a^2)
    pre = -2pi * (1 - s * Γ1 / Γ) / (Γ * W * stepsr * stepsθ)
    CPlus = pre * sumP; CMinus = pre * sumM
    # CarterConstantK flux (.wl L1863-1866): cK from the dJzd* helpers at the SHIFTED (geo+s·spin)
    # constants; Kc1 is the spin correction to K (KCorrection, taken at the GEODESIC constants).
    Ens = En0 + s * En1; Lzs = Lz0 + s * Lz1; Kcs = Kc0 + s * Kc1
    cK = (k - ω * _dJzdEn(a, s, Ens, Lzs, Kcs) - m * _dJzdLz(a, s, Ens, Lzs, Kcs)) /
         _dJzdK(a, s, Ens, Lzs, Kcs)
    return (l = l, m = m, n = n, k = k, omega = ω, α = α, CPlus = CPlus, CMinus = CMinus,
            Spi2 = SWSH(pi / 2, 0.0),
            EnergyI = abs2(CPlus) / (4pi * ω^2), EnergyH = α * abs2(CMinus) / (4pi * ω^2),
            AngMomI = abs2(CPlus) * m / (4pi * ω^3), AngMomH = α * abs2(CMinus) * m / (4pi * ω^3),
            CarterKI = cK * abs2(CPlus) / (4pi * ω^3), CarterKH = cK * α * abs2(CMinus) / (4pi * ω^3),
            stepsr = stepsr, stepsθ = stepsθ)
end

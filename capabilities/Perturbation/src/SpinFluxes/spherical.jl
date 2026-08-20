# Spherical-orbit linear-in-spin Teukolsky mode — TeukolskySpinFluxes.wl TeukolskySpinModeSpherical
# (L3331-3513). The θ-only analog of TeukolskySpinModeFromCorrection: r≡p constant, Ur=0, a SINGLE wz
# (polar Mino) sum over 2 quadrants Uz={1,-1}·(-√Z). ω = m·Ωφ + k·Ωθ (no n·Ωr). The radials Δ/K/V and
# RInrp/dRInrp/ddRInrp/dddRInrp (+Up) are computed ONCE at rp=p OUTSIDE the loop. Reuses the shared
# from_correction.jl helpers (_oc_h1/h2, _fabi/_dfabidr/_dfabidθ, _angular_fromcorr) and the module
# _alpha_circ/_dδ. CPlus/CMinus → ℐ/ℋ Energy/AngularMomentum fluxes (geodesic monopole + s·dipole). R3.

"""
    _spin_mode_spherical(oc, l, m, k, s; Nmax=80) -> NamedTuple

Core port of TeukolskySpinModeSpherical (.wl L3331-3513). `oc` is the SpinningOrbit
`KerrSpinOrbitCorrectionSpherical` NamedTuple (spherical orbit, e=0 required). Returns
`(; l,m,k,n=0,omega, CPlus,CMinus, α, Spi2, EnergyI,EnergyH,AngMomI,AngMomH, stepsθ)`. Radial via
MST-analytic TeukolskyRadial (lambda=λ,nu=ν) vs the .wl's default NumInt ⇒ ~1e-5 ℐ-flux agreement. R3.
"""
function _spin_mode_spherical(oc, l::Integer, m::Integer, k::Integer, s::Real; Nmax::Integer=80)
    (l < 2 || abs(m) > l) && error("TeukolskySpinModeSpherical: l<2 or |m|>l")
    a = oc.a; p = oc.p; e = oc.e
    iszero(e) || error("TeukolskySpinModeSpherical: requires a spherical orbit e=0 (got e=$e)")
    En0 = oc.Ehat; Lz0 = oc.Lzhat; Kc0 = oc.Khat       # geodesic constants of motion
    En1 = oc.ES; Lz1 = oc.LS                            # linear corrections to constants of motion
    correction = oc.OrbitCorrection                      # wz -> (ΔtS,rS,zS,ΔφS,UrS,UzS)
    # BL frequencies (θ,φ only): geo + s·correction
    Ωθ, Ωφ = oc.BLFrequenciesGeo .+ s .* oc.BLFrequenciesCorrection
    zgeo(wz) = cos(oc.TrajectoryGeo(wz).θ)
    Γ  = oc.MinoFrequenciesGeo[1]
    Γ1 = oc.MinoFrequenciesCorrection[1]
    ω = m * Ωφ + k * Ωθ
    isfinite(ω) || error("TeukolskySpinModeSpherical: non-real ω")

    SWSH = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω)
    λ = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ω)
    ν = RenormalizedAngularMomentum(-2, l, m, a, ω, λ)
    R = TeukolskyRadial(-2, l, m, a, ω; lambda = λ, nu = ν)
    α = _alpha_circ(a, m, ω, λ)

    # step count (.wl L3370-3371): θ-only — 32·⌈|ωΔtθ'−mΔφθ'+k|⌉ at the turning points, ≥32.
    Δtθ = oc.TrajectoryDeltas.Δtz; Δφθ = oc.TrajectoryDeltas.Δφz
    stepsθ = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/4) - m * _dδ(Δφθ, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npθ = stepsθ ÷ 2

    # trajectory corrections at the wθ midpoints (.wl `correctionp`, length stepsθ/2).
    correctionp = [correction((iθ - 1//2) * 2pi / stepsθ) for iθ in 1:npθ]

    # --- radials computed ONCE at rp=p, OUTSIDE the wz loop (.wl L3374-3389; Ur=0) ---
    rp = float(p)
    Δ  = rp^2 - 2rp + a^2
    K  = (rp^2 + a^2) * ω - a * m
    dΔ = 2 * (rp - 1)
    dK = 2rp * ω
    V  = -(K^2 + 4im * (rp - 1) * K) / Δ + 8im * ω * rp + λ
    dV = -((2K * dK + 4im * K + 4im * (rp - 1) * dK) * Δ - (K^2 + 4im * (rp - 1) * K) * dΔ) / Δ^2 + 8im * ω
    RInrp    = R.In(rp); dRInrp = teukolsky_radial_deriv(R.In, rp)
    ddRInrp  = (V * RInrp + dΔ * dRInrp) / Δ
    dddRInrp = 1 / Δ * (dV * RInrp + (V + 2) * dRInrp)
    RUprp    = R.Up(rp); dRUprp = teukolsky_radial_deriv(R.Up, rp)
    ddRUprp  = (V * RUprp + dΔ * dRUprp) / Δ
    dddRUprp = 1 / Δ * (dV * RUprp + (V + 2) * dRUprp)

    sumP = zero(ComplexF64); sumM = zero(ComplexF64)
    for iθ in 1:npθ
        wθ = (iθ - 1//2) * 2pi / stepsθ
        zp = zgeo(wθ); sinθp = sqrt(1 - zp^2)
        Uz = (1, -1) .* (-sqrt(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (Kc0 - a^2 * zp^2)))
        expθ = exp(im * (ω * Δtθ(wθ) - m * Δφθ(wθ) + 2pi * k * (iθ - 1//2) / stepsθ))
        expθq = (expθ, 1 / expθ)                                # expθ^{1,-1} (polar quadrant sign)

        ang = _angular_fromcorr(SWSH, λ, a, m, ω, zp)
        S = ang.S; L2S = ang.L2S; L1L2S = ang.L1L2S
        dSdθ = ang.dSdθ; dL2Sdθ = ang.dL2Sdθ; dL1L2Sdθ = ang.dL1L2Sdθ

        ζ = rp - im * a * zp; ζbar = rp + im * a * zp; Σ = rp^2 + a^2 * zp^2
        (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2) =
            _fabi(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a)
        (dfnn0dr, dfnmb0dr, dfnmb1dr, dfmbmb0dr, dfmbmb1dr, dfmbmb2dr) =
            _dfabidr(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a, ω)
        (dfnn0dθ, dfnmb0dθ, dfnmb1dθ, dfmbmb0dθ, dfmbmb1dθ, dfmbmb2dθ) =
            _dfabidθ(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, dSdθ, dL2Sdθ, dL1L2Sdθ, a)

        # four-velocity in Kinnersley tetrad (Ur=0; vn/vl scalars, vm 2-vector via Uz)
        vn  = -((rp^2 + a^2) * En0 - a * Lz0) / (2 * Σ)
        vl  = -((rp^2 + a^2) * En0 - a * Lz0) / Δ
        vm  = (im * (a * sinθp^2 * En0 - Lz0) .+ Uz) ./ (-sqrt(2) * sinθp * ζbar)
        vmb = conj.(vm)

        # spin tensor in Kinnersley tetrad
        Sln  = -((rp * (Kc0 - a^2 * zp^2)) / (sqrt(Kc0) * Σ))
        Snm  = (ζ / sqrt(Kc0)) .* vm .* vn
        Snmb = conj.(Snm)
        Slmb = (-(ζ / sqrt(Kc0))) .* vl .* vmb
        Smmb = (im * a * zp * (Kc0 + rp^2)) / (sqrt(Kc0) * Σ)

        # geodesic monopole (no spin)
        Amnn   = vn .^ 2
        Amnmb  = vn .* vmb
        Ammbmb = vmb .^ 2

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

        # trajectory correction at the 2 quadrants (.wl L3443-3449). cp/cm = +iθ / -iθ columns.
        cp = correctionp[iθ]
        cm = correctionp[npθ - iθ + 1]
        rp1  = (cp.rS,  cm.rS)
        zp1  = (cp.zS,  -cm.zS)
        Urp1 = (cp.UrS, cm.UrS)
        Uzp1 = (cp.UzS, -cm.UzS)
        Σ1   = 2 .* (rp .* rp1 .+ a^2 * zp .* zp1)
        exp1 = (im * (ω * cp.ΔtS - m * cp.ΔφS), im * (ω * cm.ΔtS - m * cm.ΔφS))

        h1v = _oc_h1(rp, zp, Kc0, a); h2v = _oc_h2(rp, zp, Kc0, En0, Lz0, a)
        # linear parts of the four-velocity in Kinnersley tetrad (Ur=0, no h3 term)
        vn1 = ((-1) .* (((2 * rp * En0) .* rp1) .- ((rp^2 + a^2) * En0 - a * Lz0) ./ Σ .* Σ1 .+
                 ((rp^2 + a^2) * (En1 - h1v) - a * (Lz1 + h2v) .+ Urp1))) ./ (2 * Σ)
        vmb1 = (((-im * (-2 * a * zp * En0)) .* zp1) .-
                 (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* zp1 .+ (rp1 .- im * a .* zp1) ./ ζ) .+
                 (-im * (a * sinθp^2 * (En1 - h1v) - (Lz1 + h2v)) .+ Uzp1)) ./ (-sqrt(2) * sinθp * ζ)

        # spin-dipole amplitude pieces (.wl L3454-3465; the Σ1/Σ monopole pieces fold into the *0S terms)
        Ann0S    = (Σ1 ./ Σ .* vn .+ 2 .* vn1) .* vn .+ Adnn
        AnntφS   = (Stφn .+ exp1 .* vn) .* vn
        AnnrS    = (Srn .+ rp1 .* vn) .* vn
        AnnθS    = (Sθn .- zp1 ./ sinθp .* vn) .* vn
        Anmb0S   = (Σ1 ./ Σ .* vn .* vmb .+ vn1 .* vmb .+ vn .* vmb1 .+ Adnmb)
        AnmbtφS  = ((Stφn .* vmb .+ Stφmb .* vn) ./ 2 .+ exp1 .* vn .* vmb)
        AnmbrS   = ((Srn .* vmb .+ Srmb .* vn) ./ 2 .+ rp1 .* vn .* vmb)
        AnmbθS   = ((Sθn .* vmb .+ Sθmb .* vn) ./ 2 .- zp1 ./ sinθp .* vn .* vmb)
        Ambmb0S  = (Σ1 ./ Σ .* vmb .+ 2 .* vmb1) .* vmb .+ Admbmb
        AmbmbtφS = (Stφmb .+ exp1 .* vmb) .* vmb
        AmbmbrS  = (Srmb .+ rp1 .* vmb) .* vmb
        AmbmbθS  = (Sθmb .- zp1 ./ sinθp .* vmb) .* vmb

        # the sum (.wl L3466-3473 In / L3474-3481 Up): s·Σ·(dipole) + Σ·(geodesic monopole), ×expθ^{1,-1}
        dipP = (
            ((Ann0S .+ AnntφS) .* RInrp .* fnn0 .+ AnnrS .* (dRInrp .* fnn0 .+ RInrp .* dfnn0dr) .+ AnnθS .* RInrp .* dfnn0dθ) .+
            ((Anmb0S .+ AnmbtφS) .* RInrp .* fnmb0 .+ AnmbrS .* (RInrp .* dfnmb0dr .+ dRInrp .* fnmb0) .+ AnmbθS .* RInrp .* dfnmb0dθ) .-
            ((Anmb0S .+ AnmbtφS) .* dRInrp .* fnmb1 .+ AnmbrS .* (dRInrp .* dfnmb1dr .+ ddRInrp .* fnmb1) .+ AnmbθS .* dRInrp .* dfnmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* RInrp .* fmbmb0 .+ AmbmbrS .* (dRInrp .* fmbmb0 .+ RInrp .* dfmbmb0dr) .+ AmbmbθS .* RInrp .* dfmbmb0dθ) .-
            ((Ambmb0S .+ AmbmbtφS) .* dRInrp .* fmbmb1 .+ AmbmbrS .* (ddRInrp .* fmbmb1 .+ dRInrp .* dfmbmb1dr) .+ AmbmbθS .* dRInrp .* dfmbmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* ddRInrp .* fmbmb2 .+ AmbmbrS .* (dddRInrp .* fmbmb2 .+ ddRInrp .* dfmbmb2dr) .+ AmbmbθS .* ddRInrp .* dfmbmb2dθ)) .* expθq
        monoP = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RInrp .-
                 (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRInrp .+ (Ammbmb .* fmbmb2) .* ddRInrp) .* expθq
        sumP += sum(s * Σ .* dipP .+ Σ .* monoP)

        dipM = (
            ((Ann0S .+ AnntφS) .* RUprp .* fnn0 .+ AnnrS .* (dRUprp .* fnn0 .+ RUprp .* dfnn0dr) .+ AnnθS .* RUprp .* dfnn0dθ) .+
            ((Anmb0S .+ AnmbtφS) .* RUprp .* fnmb0 .+ AnmbrS .* (RUprp .* dfnmb0dr .+ dRUprp .* fnmb0) .+ AnmbθS .* RUprp .* dfnmb0dθ) .-
            ((Anmb0S .+ AnmbtφS) .* dRUprp .* fnmb1 .+ AnmbrS .* (dRUprp .* dfnmb1dr .+ ddRUprp .* fnmb1) .+ AnmbθS .* dRUprp .* dfnmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* RUprp .* fmbmb0 .+ AmbmbrS .* (dRUprp .* fmbmb0 .+ RUprp .* dfmbmb0dr) .+ AmbmbθS .* RUprp .* dfmbmb0dθ) .-
            ((Ambmb0S .+ AmbmbtφS) .* dRUprp .* fmbmb1 .+ AmbmbrS .* (ddRUprp .* fmbmb1 .+ dRUprp .* dfmbmb1dr) .+ AmbmbθS .* dRUprp .* dfmbmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* ddRUprp .* fmbmb2 .+ AmbmbrS .* (dddRUprp .* fmbmb2 .+ ddRUprp .* dfmbmb2dr) .+ AmbmbθS .* ddRUprp .* dfmbmb2dθ)) .* expθq
        monoM = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RUprp .-
                 (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRUprp .+ (Ammbmb .* fmbmb2) .* ddRUprp) .* expθq
        sumM += sum(s * Σ .* dipM .+ Σ .* monoM)
    end

    # invariant Wronskian at rp=p (.wl L3483)
    W = (RInrp * dRUprp - dRInrp * RUprp) / (rp^2 - 2 * rp + a^2)
    pre = -2pi * (1 - s * Γ1 / Γ) / (Γ * W * stepsθ)
    CPlus = pre * sumP; CMinus = pre * sumM
    return (l = l, m = m, k = k, n = 0, omega = ω, α = α, CPlus = CPlus, CMinus = CMinus,
            Spi2 = SWSH(pi / 2, 0.0),
            EnergyI = abs2(CPlus) / (4pi * ω^2), EnergyH = α * abs2(CMinus) / (4pi * ω^2),
            AngMomI = abs2(CPlus) * m / (4pi * ω^3), AngMomH = α * abs2(CMinus) * m / (4pi * ω^3),
            stepsθ = stepsθ)
end

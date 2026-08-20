# Spherical-orbit linear-in-spin flux CORRECTION (geodesic mode + dE¹/dLz¹), NUMERICAL ω-FD —
# TeukolskySpinFluxes.wl TeukolskySpinModeSphericalCorrectionNum (L3516-3754). The θ-only analog of
# TeukolskySpinModeCorrection: a SINGLE wz (polar Mino) sum at fixed r≡p (Ur=0) over 2 quadrants. The
# s⁰ geodesic monopole (sumPlus0) and its linear-in-spin correction (sumPlus1) are SEPARATE
# accumulators. sumPlus1 = (1) the explicit-s¹ spin-dipole source at the GEODESIC point + (2) the
# ω-derivative ω1·∂_ω[geodesic monopole], computed by EXPLICIT central FD over ω of the radial R/SWSH
# (step δω) PLUS the analytic fab-basis ω-derivative `_dfabidS`. Radials Δ/K/V + RInrp1/dRInrp1/… and
# the angular linear parts S1/L2S1/L1L2S1 enter. CPlus0/1, CMinus0/1 → Fluxes + FluxesCorrection. R3.
# (_dfabidS lives in fab_basis.jl)
"""
    _spin_mode_spherical_correction_num(oc, l, m, k, δω; Nmax=80) -> NamedTuple

Core port of TeukolskySpinModeSphericalCorrectionNum (.wl L3516-3754). `oc` is the SpinningOrbit
`KerrSpinOrbitCorrectionSpherical` NamedTuple (spherical orbit, e=0). The ω-derivatives of R/SWSH use
EXPLICIT central FD with the user step `δω` (R±δω, SWSH±δω) plus the analytic `_dfabidS`. Returns
`(; l,m,k,n=0,omega,omega1, CPlus0,CMinus0,CPlus1,CMinus1, α,α1, Spi2, EnergyI,EnergyH,AngMomI,AngMomH,
dE1I,dE1H,dL1I,dL1H, λ,λ1, stepsθ)`. R3.
"""
function _spin_mode_spherical_correction_num(oc, l::Integer, m::Integer, k::Integer, δω::Real;
        Nmax::Integer=80)
    (l < 2 || abs(m) > l) && error("TeukolskySpinModeSphericalCorrectionNum: l<2 or |m|>l")
    a = oc.a; p = oc.p; e = oc.e
    iszero(e) || error("TeukolskySpinModeSphericalCorrectionNum: requires a spherical orbit e=0 (got e=$e)")
    En0 = oc.Ehat; Lz0 = oc.Lzhat; Kc0 = oc.Khat       # geodesic constants of motion
    En1 = oc.ES; Lz1 = oc.LS                            # linear corrections to constants of motion
    correction = oc.OrbitCorrection                      # wz -> (ΔtS,rS,zS,ΔφS,UrS,UzS)
    Ωθ, Ωφ = oc.BLFrequenciesGeo                         # geodesic BL frequencies (θ,φ)
    Ωθ1, Ωφ1 = oc.BLFrequenciesCorrection                # linear corrections to the BL frequencies
    zgeo(wz) = cos(oc.TrajectoryGeo(wz).θ)
    Γ  = oc.MinoFrequenciesGeo[1]
    Γ1 = oc.MinoFrequenciesCorrection[1]
    ω = m * Ωφ + k * Ωθ
    isfinite(ω) || error("TeukolskySpinModeSphericalCorrectionNum: non-real ω")
    ω1 = m * Ωφ1 + k * Ωθ1                               # linear correction to the mode frequency

    # SWSH + radial at ω and at ω±δω (the .wl numerical ω-derivative inputs; .wl L3547-3552)
    SWSH      = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω)
    SWSHplus  = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * (ω + δω))
    SWSHminus = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * (ω - δω))
    λ  = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ω)
    λp = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * (ω + δω))
    λm = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * (ω - δω))
    λ1 = (λp - λm) / (2 * δω) * ω1                        # linear correction to the eigenvalue (.wl L3554)
    ν  = RenormalizedAngularMomentum(-2, l, m, a, ω, λ)
    R      = TeukolskyRadial(-2, l, m, a, ω;      lambda = λ,  nu = ν, Nmax = Nmax)
    Rplus  = TeukolskyRadial(-2, l, m, a, ω + δω; lambda = λp,
                             nu = RenormalizedAngularMomentum(-2, l, m, a, ω + δω, λp), Nmax = Nmax)
    Rminus = TeukolskyRadial(-2, l, m, a, ω - δω; lambda = λm,
                             nu = RenormalizedAngularMomentum(-2, l, m, a, ω - δω, λm), Nmax = Nmax)

    # TS constant 𝒞2 + horizon constant α and their linear parts (.wl L3555-3563)
    C2 = ((λ + 2)^2 + 4a * ω * (m - a * ω)) * (λ^2 + 36a * ω * (m - a * ω)) -
         (2λ + 3) * (48a * ω * (m - 2a * ω)) + 144 * ω^2 * (1 - a^2)
    C21 = 4 * λ^3 * λ1 + 4 * λ^2 * (3 * λ1 + 10 * a * (m - 2 * a * ω) * ω1) +
          8 * λ * (λ1 * (1 + 10 * a * m * ω - 10 * a^2 * ω^2) + 6 * a * (m + 2 * a * ω) * ω1) +
          48 * ω * (a * m * λ1 + 6 * ω1 - 18 * a^3 * m * ω * ω1 + 12 * a^4 * ω^2 * ω1 +
                    a^2 * (λ1 * ω + 6 * m^2 * ω1))
    rplus = 1 + sqrt(1 - a^2)
    P = ω - m * a / (2 * rplus)
    ϵ = sqrt(1 - a^2) / (4 * rplus)
    α = 256 * (2 * rplus)^5 * P * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * ω^3 / C2
    α1 = -256 * (2 * rplus)^5 * P * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * ω^3 / C2^2 * C21 +
         256 * (2 * rplus)^5 * (ω1 * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * ω^3 +
            P * (2 * P * ω1) * (P^2 + 16 * ϵ^2) * ω^3 +
            P * (P^2 + 4 * ϵ^2) * (2 * P * ω1) * ω^3 +
            P * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * 3 * ω^2 * ω1) / C2

    # step count (.wl L3567-3568): θ-only, ≥32.
    Δtθ = oc.TrajectoryDeltas.Δtz; Δφθ = oc.TrajectoryDeltas.Δφz
    stepsθ = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/4) - m * _dδ(Δφθ, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npθ = stepsθ ÷ 2

    # trajectory corrections at the wθ midpoints (.wl `correctionp`, length stepsθ/2).
    correctionp = [correction((iθ - 1//2) * 2pi / stepsθ) for iθ in 1:npθ]

    # --- radials computed ONCE at rp=p, OUTSIDE the wz loop (.wl L3571-3595; Ur=0) ---
    rp = float(p)
    Δ  = rp^2 - 2rp + a^2
    K  = (rp^2 + a^2) * ω - a * m
    K1 = (rp^2 + a^2) * ω1
    dΔ = 2 * (rp - 1)
    dK = 2rp * ω
    dK1 = 2rp * ω1
    V  = -(K^2 + 4im * (rp - 1) * K) / Δ + 8im * ω * rp + λ
    V1 = -(2 * K * K1 + 4im * (rp - 1) * K1) / Δ + 8im * ω1 * rp + λ1
    dV = -((2K * dK + 4im * K + 4im * (rp - 1) * dK) * Δ - (K^2 + 4im * (rp - 1) * K) * dΔ) / Δ^2 + 8im * ω
    RInrp    = R.In(rp); dRInrp = teukolsky_radial_deriv(R.In, rp)
    ddRInrp  = (V * RInrp + dΔ * dRInrp) / Δ
    RInrp1   = (Rplus.In(rp) - Rminus.In(rp)) / (2 * δω) * ω1
    dRInrp1  = (teukolsky_radial_deriv(Rplus.In, rp) - teukolsky_radial_deriv(Rminus.In, rp)) / (2 * δω) * ω1
    ddRInrp1 = (V * RInrp1 + V1 * RInrp + dΔ * dRInrp1) / Δ
    dddRInrp = 1 / Δ * (dV * RInrp + (V + 2) * dRInrp)
    RUprp    = R.Up(rp); dRUprp = teukolsky_radial_deriv(R.Up, rp)
    ddRUprp  = (V * RUprp + dΔ * dRUprp) / Δ
    RUprp1   = (Rplus.Up(rp) - Rminus.Up(rp)) / (2 * δω) * ω1
    dRUprp1  = (teukolsky_radial_deriv(Rplus.Up, rp) - teukolsky_radial_deriv(Rminus.Up, rp)) / (2 * δω) * ω1
    ddRUprp1 = (V * RUprp1 + V1 * RUprp + dΔ * dRUprp1) / Δ
    dddRUprp = 1 / Δ * (dV * RUprp + (V + 2) * dRUprp)

    sumP0 = sumM0 = sumP1 = sumM1 = zero(ComplexF64)
    for iθ in 1:npθ
        wθ = (iθ - 1//2) * 2pi / stepsθ
        zp = zgeo(wθ); sinθp = sqrt(1 - zp^2)
        Uz = (1, -1) .* (-sqrt(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (Kc0 - a^2 * zp^2)))
        expθ = exp(im * (ω * Δtθ(wθ) - m * Δφθ(wθ) + 2pi * k * (iθ - 1//2) / stepsθ))
        expθq = (expθ, 1 / expθ)                               # expθ^{1,-1} (polar quadrant sign)
        Δtθwθ = Δtθ(wθ)
        θp = acos(zp)

        # geodesic angular layer (.wl L3602,L3604,L3606,L3608,L3610,L3612-3613,L3616-3617)
        ang = _angular_fromcorr(SWSH, λ, a, m, ω, zp)
        S = ang.S; dSdθ = ang.dSdθ; d2Sdθ2 = ang.d2Sdθ2
        L2S = ang.L2S; L1L2S = ang.L1L2S; dL2Sdθ = ang.dL2Sdθ; dL1L2Sdθ = ang.dL1L2Sdθ

        # ω-linear angular parts (.wl L3603,L3605,L3607,L3611,L3614-3615)
        S1   = (SWSHplus(θp, 0.0) - SWSHminus(θp, 0.0)) / (2 * δω) * ω1
        dS1dθ = (swsh_derivative(SWSHplus, 1, 0, θp, 0.0) - swsh_derivative(SWSHminus, 1, 0, θp, 0.0)) / (2 * δω) * ω1
        d2S1dθ2 = -(-(a * ω)^2 * (1 - zp^2) - (m - 2 * zp)^2 / (1 - zp^2) + 4a * ω * zp - 2 + 2 * m * a * ω + λ) * S1 -
                  (-a^2 * 2 * ω * ω1 * (1 - zp^2) + 4 * a * zp * ω1 + 2 * m * a * ω1 + λ1) * S -
                  zp / sinθp * dS1dθ
        L2S1 = dS1dθ - (S1 * (m - 2 * zp + a * (-1 + zp^2) * ω) + S * (a * (-1 + zp^2) * ω1)) / sinθp
        L1L2S1 = d2S1dθ2 +
                 (dS1dθ * (-2 * m + 3 * zp - 2 * a * (-1 + zp^2) * ω) + dSdθ * (-2 * a * (-1 + zp^2) * ω1)) / sinθp +
                 S1 * (-2 - (m * (m - 2 * zp)) / (-1 + zp^2) - 2 * a * (m - 2 * zp) * ω - a^2 * (-1 + zp^2) * ω^2) +
                 S * (-2 * a * (m - 2 * zp) * ω1 - a^2 * (-1 + zp^2) * 2 * ω * ω1)

        ζ = rp - im * a * zp; ζbar = rp + im * a * zp; Σ = rp^2 + a^2 * zp^2
        (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2) =
            _fabi(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a)
        (dfnn0dr, dfnmb0dr, dfnmb1dr, dfmbmb0dr, dfmbmb1dr, dfmbmb2dr) =
            _dfabidr(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, a, ω)
        (dfnn0dθ, dfnmb0dθ, dfnmb1dθ, dfmbmb0dθ, dfmbmb1dθ, dfmbmb2dθ) =
            _dfabidθ(ζ, ζbar, sinθp, Δ, dΔ, K, dK, S, L2S, L1L2S, dSdθ, dL2Sdθ, dL1L2Sdθ, a)
        (fnn01, fnmb01, fnmb11, fmbmb01, fmbmb11, fmbmb21) =
            _dfabidS(ζ, ζbar, sinθp, Δ, dΔ, K, dK, K1, dK1, S, L2S, L1L2S, S1, L2S1, L1L2S1, a)

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

        # trajectory correction at the 2 quadrants (.wl L3656-3662). cp/cm = +iθ / -iθ columns.
        cp = correctionp[iθ]
        cm = correctionp[npθ - iθ + 1]
        rp1  = (cp.rS,  cm.rS)
        zp1  = (cp.zS,  -cm.zS)
        Urp1 = (cp.UrS, cm.UrS)
        Uzp1 = (cp.UzS, -cm.UzS)
        Σ1   = 2 .* (rp .* rp1 .+ a^2 * zp .* zp1)
        # phase linear part (.wl L3661-3662): O(s) trajectory shift + ω1·(geodesic Δtθ)
        exp1 = (im * ((ω * cp.ΔtS - m * cp.ΔφS) + ω1 * Δtθwθ),
                im * ((ω * cm.ΔtS - m * cm.ΔφS) - ω1 * Δtθwθ))

        h1v = _oc_h1(rp, zp, Kc0, a); h2v = _oc_h2(rp, zp, Kc0, En0, Lz0, a)
        # linear parts of the four-velocity in Kinnersley tetrad (Ur=0, no h3 term) (.wl L3663-3666)
        vn1 = ((-1) .* (((2 * rp * En0) .* rp1) .- ((rp^2 + a^2) * En0 - a * Lz0) ./ Σ .* Σ1 .+
                 ((rp^2 + a^2) * (En1 - h1v) - a * (Lz1 + h2v) .+ Urp1))) ./ (2 * Σ)
        vmb1 = (((-im * (-2 * a * zp * En0)) .* zp1) .-
                 (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* zp1 .+ (rp1 .- im * a .* zp1) ./ ζ) .+
                 (-im * (a * sinθp^2 * (En1 - h1v) - (Lz1 + h2v)) .+ Uzp1)) ./ (-sqrt(2) * sinθp * ζ)

        # explicit-s¹ spin-dipole source pieces (.wl L3667-3678)
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

        # geodesic monopole sum (sumPlus0/sumMinus0; .wl L3679-3680), ×expθ^{1,-1}
        monoP = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RInrp .-
                 (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRInrp .+ Ammbmb .* fmbmb2 .* ddRInrp) .* expθq
        monoM = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RUprp .-
                 (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRUprp .+ Ammbmb .* fmbmb2 .* ddRUprp) .* expθq
        sumP0 += sum(Σ .* monoP)
        sumM0 += sum(Σ .* monoM)

        # correction sum sumPlus1 (.wl L3681-3689): (1) spin-dipole source × radials + (2) ω-deriv
        dipP = (
            ((Ann0S .+ AnntφS) .* RInrp .* fnn0 .+ AnnrS .* (dRInrp .* fnn0 .+ RInrp .* dfnn0dr) .+ AnnθS .* RInrp .* dfnn0dθ) .+
            ((Anmb0S .+ AnmbtφS) .* RInrp .* fnmb0 .+ AnmbrS .* (RInrp .* dfnmb0dr .+ dRInrp .* fnmb0) .+ AnmbθS .* RInrp .* dfnmb0dθ) .-
            ((Anmb0S .+ AnmbtφS) .* dRInrp .* fnmb1 .+ AnmbrS .* (dRInrp .* dfnmb1dr .+ ddRInrp .* fnmb1) .+ AnmbθS .* dRInrp .* dfnmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* RInrp .* fmbmb0 .+ AmbmbrS .* (dRInrp .* fmbmb0 .+ RInrp .* dfmbmb0dr) .+ AmbmbθS .* RInrp .* dfmbmb0dθ) .-
            ((Ambmb0S .+ AmbmbtφS) .* dRInrp .* fmbmb1 .+ AmbmbrS .* (ddRInrp .* fmbmb1 .+ dRInrp .* dfmbmb1dr) .+ AmbmbθS .* dRInrp .* dfmbmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* ddRInrp .* fmbmb2 .+ AmbmbrS .* (dddRInrp .* fmbmb2 .+ ddRInrp .* dfmbmb2dr) .+ AmbmbθS .* ddRInrp .* dfmbmb2dθ)) .* expθq
        ωdP = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RInrp1 .-
               (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRInrp1 .+ Ammbmb .* fmbmb2 .* ddRInrp1 .+
               (Amnn .* fnn01 .+ Amnmb .* fnmb01 .+ Ammbmb .* fmbmb01) .* RInrp .-
               (Amnmb .* fnmb11 .+ Ammbmb .* fmbmb11) .* dRInrp .+ Ammbmb .* fmbmb21 .* ddRInrp) .* expθq
        sumP1 += sum(Σ .* dipP .+ Σ .* ωdP)

        dipM = (
            ((Ann0S .+ AnntφS) .* RUprp .* fnn0 .+ AnnrS .* (dRUprp .* fnn0 .+ RUprp .* dfnn0dr) .+ AnnθS .* RUprp .* dfnn0dθ) .+
            ((Anmb0S .+ AnmbtφS) .* RUprp .* fnmb0 .+ AnmbrS .* (RUprp .* dfnmb0dr .+ dRUprp .* fnmb0) .+ AnmbθS .* RUprp .* dfnmb0dθ) .-
            ((Anmb0S .+ AnmbtφS) .* dRUprp .* fnmb1 .+ AnmbrS .* (dRUprp .* dfnmb1dr .+ ddRUprp .* fnmb1) .+ AnmbθS .* dRUprp .* dfnmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* RUprp .* fmbmb0 .+ AmbmbrS .* (dRUprp .* fmbmb0 .+ RUprp .* dfmbmb0dr) .+ AmbmbθS .* RUprp .* dfmbmb0dθ) .-
            ((Ambmb0S .+ AmbmbtφS) .* dRUprp .* fmbmb1 .+ AmbmbrS .* (ddRUprp .* fmbmb1 .+ dRUprp .* dfmbmb1dr) .+ AmbmbθS .* dRUprp .* dfmbmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* ddRUprp .* fmbmb2 .+ AmbmbrS .* (dddRUprp .* fmbmb2 .+ ddRUprp .* dfmbmb2dr) .+ AmbmbθS .* ddRUprp .* dfmbmb2dθ)) .* expθq
        ωdM = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RUprp1 .-
               (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRUprp1 .+ Ammbmb .* fmbmb2 .* ddRUprp1 .+
               (Amnn .* fnn01 .+ Amnmb .* fnmb01 .+ Ammbmb .* fmbmb01) .* RUprp .-
               (Amnmb .* fnmb11 .+ Ammbmb .* fmbmb11) .* dRUprp .+ Ammbmb .* fmbmb21 .* ddRUprp) .* expθq
        sumM1 += sum(Σ .* dipM .+ Σ .* ωdM)
    end

    # invariant Wronskian + its linear part (.wl L3700-3701)
    W  = (RInrp * dRUprp - dRInrp * RUprp) / Δ
    W1 = (RInrp1 * dRUprp + RInrp * dRUprp1 - dRInrp1 * RUprp - dRInrp * RUprp1) / Δ
    den = Γ * W * stepsθ
    CPlus0  = -2pi * sumP0 / den
    CMinus0 = -2pi * sumM0 / den
    CPlus1  = -2pi * (sumP1 - Γ1 / Γ * sumP0 - W1 / W * sumP0) / den       # (.wl L3704)
    CMinus1 = -2pi * (sumM1 - Γ1 / Γ * sumM0 - W1 / W * sumM0) / den

    # Fluxes (.wl L3727-3735) + FluxesCorrection dE¹/dLz¹ (.wl L3741-3747)
    EnergyI = abs2(CPlus0) / (4pi * ω^2);  EnergyH = α * abs2(CMinus0) / (4pi * ω^2)
    AngMomI = abs2(CPlus0) * m / (4pi * ω^3); AngMomH = α * abs2(CMinus0) * m / (4pi * ω^3)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (α1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (α1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    return (l = l, m = m, k = k, n = 0, omega = ω, omega1 = ω1,
            CPlus0 = CPlus0, CMinus0 = CMinus0, CPlus1 = CPlus1, CMinus1 = CMinus1,
            α = α, α1 = α1, Spi2 = SWSH(pi / 2, 0.0),
            EnergyI = EnergyI, EnergyH = EnergyH, AngMomI = AngMomI, AngMomH = AngMomH,
            dE1I = dE1I, dE1H = dE1H, dL1I = dL1I, dL1H = dL1H,
            λ = λ, λ1 = λ1, stepsθ = stepsθ)
end

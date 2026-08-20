# Spherical-orbit linear-in-spin flux CORRECTION + r/x DERIVATIVES, ANALYTIC ω-derivatives —
# TeukolskySpinFluxes.wl TeukolskySpinModeSphericalCorrection (L3759-4110). The largest flux function:
# a SINGLE wz (polar Mino) sum at r≡p (Ur=0) accumulating EIGHT sums (geodesic monopole sumPlus0/Minus0,
# spin correction sumPlus1/Minus1, and r/x derivative sums dsumPlus/Minus dr/dx). Geodesic monopole +
# spin-dipole source are IDENTICAL to the Num version; the ω-derivative of the monopole is ANALYTIC here
# (dfabidω·ω1 + dRIndω·ω1 from TeukolskySolverHS1spin), not δω-FD. The r/x derivative sums are a full
# chain-rule projection with ∂r/∂x of every factor. → CPlus0/1, dCPlusdr/dx → Fluxes + FluxesCorrection
# (== Num version) + FluxesDerivatives (the new inspiral-evolution content). R3.
# (_dfabidω lives in fab_basis.jl)

"""
    _spin_mode_spherical_correction(oc, odv, l, m, k; Nmax=80, hh=1e-6) -> NamedTuple

Core port of TeukolskySpinModeSphericalCorrection (.wl L3759-4110). `oc` is the SpinningOrbit
`KerrSpinOrbitCorrectionSpherical` NamedTuple (spherical orbit, e=0); `odv` is the matching
`_orbit_derivatives_spherical` object (r/x FD derivatives of the orbit correction). The ω-derivatives
of R/SWSH are ANALYTIC: `dRIndω` etc. via central FD of `TeukolskyRadial(ω±hh)` (the
`TeukolskySolverHS1spin` callback), `dSdω`/`d2Sdθdω` via FD of `SWSH(ω±hh)`, and the analytic
`_dfabidω`. Returns geodesic + FluxesCorrection (dE¹/dLz¹) + FluxesDerivatives (∂r/∂x of the flux). R3.
"""
function _spin_mode_spherical_correction(oc, odv, l::Integer, m::Integer, k::Integer;
        Nmax::Integer=80, hh::Real=1e-6)
    (l < 2 || abs(m) > l) && error("TeukolskySpinModeSphericalCorrection: l<2 or |m|>l")
    a = oc.a; p = oc.p; e = oc.e
    iszero(e) || error("TeukolskySpinModeSphericalCorrection: requires a spherical orbit e=0 (got e=$e)")
    En0 = oc.Ehat; Lz0 = oc.Lzhat; Kc0 = oc.Khat        # geodesic constants of motion
    En1 = oc.ES; Lz1 = oc.LS                             # linear corrections to constants of motion
    dEndr = odv.dEndr; dLzdr = odv.dLzdr                 # r/x derivatives of constants of motion
    dEndx = odv.dEndx; dLzdx = odv.dLzdx
    correction = oc.OrbitCorrection                      # wz -> (ΔtS,rS,zS,ΔφS,UrS,UzS)
    derivatives = odv.OrbitDerivatives                   # wz -> (dzdr,dzdx,dUzdr,dUzdx,dΔtdr,dΔφdr,dΔtdx,dΔφdx)
    Ωθ, Ωφ = oc.BLFrequenciesGeo                         # geodesic BL frequencies (θ,φ)
    Ωθ1, Ωφ1 = oc.BLFrequenciesCorrection                # linear corrections to the BL frequencies
    dΩθdr, dΩφdr = odv.dBLFrequenciesdr                  # r/x derivatives of the BL frequencies
    dΩθdx, dΩφdx = odv.dBLFrequenciesdx
    zgeo(wz) = cos(oc.TrajectoryGeo(wz).θ)
    Γ  = oc.MinoFrequenciesGeo[1]
    Γ1 = oc.MinoFrequenciesCorrection[1]
    dΓdr = odv.dMinoFrequenciesdr[1]
    dΓdx = odv.dMinoFrequenciesdx[1]
    ω = m * Ωφ + k * Ωθ                                  # geodesic mode frequency + corrections/derivatives
    isfinite(ω) || error("TeukolskySpinModeSphericalCorrection: non-real ω")
    ω1   = m * Ωφ1  + k * Ωθ1
    dωdr = m * dΩφdr + k * dΩθdr
    dωdx = m * dΩφdx + k * dΩθdx

    # --- angular layer (angparNew): λ, dλdω (FD), SWSH(ω) and SWSH(ω±hh) for the analytic ω-deriv ---
    SWSH      = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω)
    SWSHplus  = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * (ω + hh))
    SWSHminus = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * (ω - hh))
    λ  = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ω)
    λp = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * (ω + hh))
    λm = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * (ω - hh))
    dλdω = (λp - λm) / (2 * hh)                           # angparNew dλdω

    # --- radials (TeukolskySolverHS1spin): R at ω and at ω±hh for the analytic dRdω/d2Rdrdω ---
    ν  = RenormalizedAngularMomentum(-2, l, m, a, ω, λ)
    R      = TeukolskyRadial(-2, l, m, a, ω;      lambda = λ,  nu = ν, Nmax = Nmax)
    Rplus  = TeukolskyRadial(-2, l, m, a, ω + hh; lambda = λp,
                             nu = RenormalizedAngularMomentum(-2, l, m, a, ω + hh, λp), Nmax = Nmax)
    Rminus = TeukolskyRadial(-2, l, m, a, ω - hh; lambda = λm,
                             nu = RenormalizedAngularMomentum(-2, l, m, a, ω - hh, λm), Nmax = Nmax)

    # TS constant 𝒞2 + horizon constant α and their ω-derivatives (.wl L3826-3834)
    C2 = ((λ + 2)^2 + 4a * ω * (m - a * ω)) * (λ^2 + 36a * ω * (m - a * ω)) -
         (2λ + 3) * (48a * ω * (m - 2a * ω)) + 144 * ω^2 * (1 - a^2)
    dC2dω = 4 * λ^3 * dλdω + 4 * λ^2 * (3 * dλdω + 10 * a * (m - 2 * a * ω)) +
            8 * λ * (dλdω * (1 + 10 * a * m * ω - 10 * a^2 * ω^2) + 6 * a * (m + 2 * a * ω)) +
            48 * ω * (a * m * dλdω + 6 - 18 * a^3 * m * ω + 12 * a^4 * ω^2 + a^2 * (dλdω * ω + 6 * m^2))
    rplus = 1 + sqrt(1 - a^2)
    P = ω - m * a / (2 * rplus)
    ϵ = sqrt(1 - a^2) / (4 * rplus)
    α = 256 * (2 * rplus)^5 * P * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * ω^3 / C2
    dαdω = -256 * (2 * rplus)^5 * P * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * ω^3 / C2^2 * dC2dω +
           256 * (2 * rplus)^5 * ((P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * ω^3 +
              P * (2 * P) * (P^2 + 16 * ϵ^2) * ω^3 + P * (P^2 + 4 * ϵ^2) * (2 * P) * ω^3 +
              P * (P^2 + 4 * ϵ^2) * (P^2 + 16 * ϵ^2) * 3 * ω^2) / C2

    # step count (.wl L3836-3837): θ-only, ≥32. (pi/2 and 0 endpoints, as in the analytic .wl)
    Δtθ = oc.TrajectoryDeltas.Δtz; Δφθ = oc.TrajectoryDeltas.Δφz
    stepsθ = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/2) - m * _dδ(Δφθ, pi/2) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npθ = stepsθ ÷ 2

    # trajectory corrections + derivatives at the wθ midpoints (.wl `correctionp`/`derivativesp`).
    correctionp = [correction((iθ - 1//2) * 2pi / stepsθ) for iθ in 1:npθ]
    derivativesp = [derivatives((iθ - 1//2) * 2pi / stepsθ) for iθ in 1:npθ]

    # --- radials computed ONCE at rp=p, OUTSIDE the wz loop (.wl L3841-3857; Ur=0) ---
    rp = float(p)
    Δ  = rp^2 - 2rp + a^2
    K  = (rp^2 + a^2) * ω - a * m
    dKdω = (rp^2 + a^2)
    dΔ = 2 * (rp - 1)
    dKdr = 2rp * ω
    d2Kdrdω = 2rp
    V  = -(K^2 + 4im * (rp - 1) * K) / Δ + 8im * ω * rp + λ
    dVdω = -(2 * K * dKdω + 4im * (rp - 1) * dKdω) / Δ + 8im * rp + dλdω
    dVdr = -((2K * dKdr + 4im * K + 4im * (rp - 1) * dKdr) * Δ - (K^2 + 4im * (rp - 1) * K) * dΔ) / Δ^2 + 8im * ω

    RIn = R.In(rp); dRIndr = teukolsky_radial_deriv(R.In, rp)
    dRIndω = (Rplus.In(rp) - Rminus.In(rp)) / (2 * hh)
    d2RIndrdω = (teukolsky_radial_deriv(Rplus.In, rp) - teukolsky_radial_deriv(Rminus.In, rp)) / (2 * hh)
    d2RIndr2 = (V * RIn + dΔ * dRIndr) / Δ
    d3RIndr3 = 1 / Δ * (dVdr * RIn + (V + 2) * dRIndr)
    d3RIndr2dω = (V * dRIndω + dVdω * RIn + dΔ * d2RIndrdω) / Δ
    RUp = R.Up(rp); dRUpdr = teukolsky_radial_deriv(R.Up, rp)
    dRUpdω = (Rplus.Up(rp) - Rminus.Up(rp)) / (2 * hh)
    d2RUpdrdω = (teukolsky_radial_deriv(Rplus.Up, rp) - teukolsky_radial_deriv(Rminus.Up, rp)) / (2 * hh)
    d2RUpdr2 = (V * RUp + dΔ * dRUpdr) / Δ
    d3RUpdr3 = 1 / Δ * (dVdr * RUp + (V + 2) * dRUpdr)
    d3RUpdr2dω = (V * dRUpdω + dVdω * RUp + dΔ * d2RUpdrdω) / Δ

    sumP0 = sumM0 = sumP1 = sumM1 = zero(ComplexF64)
    dsumPdr = dsumMdr = dsumPdx = dsumMdx = zero(ComplexF64)
    for iθ in 1:npθ
        wθ = (iθ - 1//2) * 2pi / stepsθ
        zp = zgeo(wθ); sinθp = sqrt(1 - zp^2)
        Uz = (1, -1) .* (-sqrt(-((1 - zp^2) * a * En0 - Lz0)^2 + (1 - zp^2) * (Kc0 - a^2 * zp^2)))
        expθ = exp(im * (ω * Δtθ(wθ) - m * Δφθ(wθ) + 2pi * k * (iθ - 1//2) / stepsθ))
        expθq = (expθ, 1 / expθ)                              # expθ^{1,-1} (polar quadrant sign)
        Δtθwθ = Δtθ(wθ)
        θp = acos(zp)

        # geodesic angular layer (.wl L3865-3876)
        ang = _angular_fromcorr(SWSH, λ, a, m, ω, zp)
        S = ang.S; dSdθ = ang.dSdθ; d2Sdθ2 = ang.d2Sdθ2
        L2S = ang.L2S; L1L2S = ang.L1L2S; dL2Sdθ = ang.dL2Sdθ; dL1L2Sdθ = ang.dL1L2Sdθ

        # analytic ω-derivatives of the angular functions (.wl L3866-3878).
        # dSdω, d2Sdθdω from angparNew's dSWSHdω (FD of SWSH); higher ω-derivs analytic from the TE.
        dSdω   = (SWSHplus(θp, 0.0) - SWSHminus(θp, 0.0)) / (2 * hh)
        d2Sdθdω = (swsh_derivative(SWSHplus, 1, 0, θp, 0.0) - swsh_derivative(SWSHminus, 1, 0, θp, 0.0)) / (2 * hh)
        d3Sdθ2dω = -(zp / sinθp) * d2Sdθdω -
                   (-a^2 * ω^2 * (1 - zp^2) - (m - 2 * zp)^2 / (1 - zp^2) + 4 * a * ω * zp - 2 + 2 * m * a * ω + λ) * dSdω -
                   (-2 * a^2 * ω * (1 - zp^2) + 4 * a * zp + 2 * m * a + dλdω) * S
        dL2Sdω = (d2Sdθdω - (dSdω * (m - 2 * zp + a * (-1 + zp^2) * ω) + S * (a * (-1 + zp^2))) / sinθp)
        dL1L2Sdω = (d3Sdθ2dω + (d2Sdθdω * (-2 * m + 3 * zp - 2 * a * (-1 + zp^2) * ω) + dSdθ * (-2 * a * (-1 + zp^2))) / sinθp +
                   dSdω * (-2 - (m * (m - 2 * zp)) / (-1 + zp^2) - 2 * a * (m - 2 * zp) * ω - a^2 * (-1 + zp^2) * ω^2) +
                   S * (-2 * a * (m - 2 * zp) - a^2 * (-1 + zp^2) * 2 * ω))

        ζ = rp - im * a * zp; ζbar = rp + im * a * zp; Σ = rp^2 + a^2 * zp^2
        (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2) =
            _fabi(ζ, ζbar, sinθp, Δ, dΔ, K, dKdr, S, L2S, L1L2S, a)
        (dfnn0dr, dfnmb0dr, dfnmb1dr, dfmbmb0dr, dfmbmb1dr, dfmbmb2dr) =
            _dfabidr(ζ, ζbar, sinθp, Δ, dΔ, K, dKdr, S, L2S, L1L2S, a, ω)
        (dfnn0dθ, dfnmb0dθ, dfnmb1dθ, dfmbmb0dθ, dfmbmb1dθ, dfmbmb2dθ) =
            _dfabidθ(ζ, ζbar, sinθp, Δ, dΔ, K, dKdr, S, L2S, L1L2S, dSdθ, dL2Sdθ, dL1L2Sdθ, a)
        (dfnn0dω, dfnmb0dω, dfnmb1dω, dfmbmb0dω, dfmbmb1dω, dfmbmb2dω) =
            _dfabidω(ζ, ζbar, sinθp, Δ, dΔ, K, dKdr, dKdω, d2Kdrdω, S, L2S, L1L2S, dSdω, dL2Sdω, dL1L2Sdω, a)

        # four-velocity in Kinnersley tetrad (Ur=0; vn/vl scalars, vm 2-vector via Uz) (.wl L3886-3889)
        vn  = -((rp^2 + a^2) * En0 - a * Lz0) / (2 * Σ)
        vl  = -((rp^2 + a^2) * En0 - a * Lz0) / Δ
        vm  = (im * (a * sinθp^2 * En0 - Lz0) .+ Uz) ./ (-sqrt(2) * sinθp * ζbar)
        vmb = conj.(vm)

        # spin tensor in Kinnersley tetrad (.wl L3891-3895)
        Sln  = -((rp * (Kc0 - a^2 * zp^2)) / (sqrt(Kc0) * Σ))
        Snm  = (ζ / sqrt(Kc0)) .* vm .* vn
        Snmb = conj.(Snm)
        Slmb = (-(ζ / sqrt(Kc0))) .* vl .* vmb
        Smmb = (im * a * zp * (Kc0 + rp^2)) / (sqrt(Kc0) * Σ)

        # geodesic monopole (no spin)
        Amnn   = vn .^ 2
        Amnmb  = vn .* vmb
        Ammbmb = vmb .^ 2

        # Newman-Penrose spin coefficients (.wl L3899-3905)
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
        # functions from Eq. (B3) in 2303.16798 (.wl L3912-3917)
        Stφn  = -im * K / (2 * Σ) * Sln .+ (a * ω * sinθp - m / sinθp) / (sqrt(2) * Σ) .* (ζ .* Snmb .- ζbar .* Snm)
        Stφmb = -im * K .* (1 / Δ .* Snmb .+ 1 / (2 * Σ) .* Slmb) .+ (a * ω * sinθp - m / sinθp) / (sqrt(2) * ζ) .* Smmb
        Srn   = Δ / (2 * Σ) * Sln
        Srmb  = ((-1) .* Snmb) .+ Δ / (2 * Σ) .* Slmb
        Sθn   = ((-1) .* (Snmb .* ζ .+ Snm .* ζbar)) ./ (sqrt(2) * Σ)
        Sθmb  = Smmb / (sqrt(2) * ζ)

        # trajectory correction at the 2 quadrants (.wl L3919-3928). cp/cm = +iθ / -iθ columns.
        cp = correctionp[iθ]; cm = correctionp[npθ - iθ + 1]
        dp = derivativesp[iθ]; dm = derivativesp[npθ - iθ + 1]
        rp1  = (cp.rS,  cm.rS)
        zp1  = (cp.zS,  -cm.zS)
        Urp1 = (cp.UrS, cm.UrS)
        Uzp1 = (cp.UzS, -cm.UzS)
        dzpdr = (dp.dzdr, -dm.dzdr)
        dzpdx = (dp.dzdx, -dm.dzdx)
        dθpdr = (-1 / sinθp) .* dzpdr
        dθpdx = (-1 / sinθp) .* dzpdx
        dUzpdr = (dp.dUzdr, -dm.dUzdr)
        dUzpdx = (dp.dUzdx, -dm.dUzdx)
        Σ1   = 2 .* (rp .* rp1 .+ a^2 * zp .* zp1)
        dΣdr = 2 .* (rp .+ a^2 * zp .* dzpdr)
        dΣdx = 2 * a^2 * zp .* dzpdx
        # phase linear part + r/x derivatives (.wl L3932-3937)
        exp1 = (im * ((ω * cp.ΔtS - m * cp.ΔφS) + ω1 * Δtθwθ),
                im * ((ω * cm.ΔtS - m * cm.ΔφS) - ω1 * Δtθwθ))
        dexpdr = (im * ((ω * dp.dΔtdr - m * dp.dΔφdr) + dωdr * Δtθwθ),
                  im * ((ω * dm.dΔtdr - m * dm.dΔφdr) - dωdr * Δtθwθ))
        dexpdx = (im * ((ω * dp.dΔtdx - m * dp.dΔφdx) + dωdx * Δtθwθ),
                  im * ((ω * dm.dΔtdx - m * dm.dΔφdx) - dωdx * Δtθwθ))

        h1v = _oc_h1(rp, zp, Kc0, a); h2v = _oc_h2(rp, zp, Kc0, En0, Lz0, a)
        # linear parts + r/x derivatives of the four-velocity in Kinnersley tetrad (Ur=0) (.wl L3938-3949)
        vn1 = ((-1) .* (((2 * rp * En0) .* rp1) .- ((rp^2 + a^2) * En0 - a * Lz0) ./ Σ .* Σ1 .+
                 ((rp^2 + a^2) * (En1 - h1v) - a * (Lz1 + h2v) .+ Urp1))) ./ (2 * Σ)
        dvndr = ((-1) .* ((2 * rp * En0) .- ((rp^2 + a^2) * En0 - a * Lz0) ./ Σ .* dΣdr .+
                 ((rp^2 + a^2) * dEndr - a * dLzdr))) ./ (2 * Σ)
        dvndx = ((-1) .* ((-1) .* ((rp^2 + a^2) * En0 - a * Lz0) ./ Σ .* dΣdx .+
                 ((rp^2 + a^2) * dEndx - a * dLzdx))) ./ (2 * Σ)
        vmb1 = (((-im * (-2 * a * zp * En0)) .* zp1) .-
                 (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* zp1 .+ (rp1 .- im * a .* zp1) ./ ζ) .+
                 (-im * (a * sinθp^2 * (En1 - h1v) - (Lz1 + h2v)) .+ Uzp1)) ./ (-sqrt(2) * sinθp * ζ)
        dvmbdr = (((-im * (-2 * a * zp * En0)) .* dzpdr) .-
                 (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* dzpdr .+ (1 .- im * a .* dzpdr) ./ ζ) .+
                 (-im * (a * sinθp^2 * dEndr - dLzdr) .+ dUzpdr)) ./ (-sqrt(2) * sinθp * ζ)
        dvmbdx = (((-im * (-2 * a * zp * En0)) .* dzpdx) .-
                 (-im * (a * sinθp^2 * En0 - Lz0) .+ Uz) .* (((-zp) / sinθp^2) .* dzpdx .+ ((-im * a) .* dzpdx) ./ ζ) .+
                 (-im * (a * sinθp^2 * dEndx - dLzdx) .+ dUzpdx)) ./ (-sqrt(2) * sinθp * ζ)

        # explicit-s¹ spin-dipole source pieces (.wl L3950-3961)
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

        # --- 1. geodesic monopole sums (.wl L3962-3963) ---
        monoP = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RIn .-
                 (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRIndr .+ Ammbmb .* fmbmb2 .* d2RIndr2) .* expθq
        monoM = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* RUp .-
                 (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* dRUpdr .+ Ammbmb .* fmbmb2 .* d2RUpdr2) .* expθq
        sumP0 += sum(Σ .* monoP)
        sumM0 += sum(Σ .* monoM)

        # --- 2. spin correction sums sumPlus1/Minus1 (.wl L3964-3981): dipole source + analytic ω-deriv ---
        dipP = (
            ((Ann0S .+ AnntφS) .* RIn .* fnn0 .+ AnnrS .* (dRIndr .* fnn0 .+ RIn .* dfnn0dr) .+ AnnθS .* RIn .* dfnn0dθ) .+
            ((Anmb0S .+ AnmbtφS) .* RIn .* fnmb0 .+ AnmbrS .* (RIn .* dfnmb0dr .+ dRIndr .* fnmb0) .+ AnmbθS .* RIn .* dfnmb0dθ) .-
            ((Anmb0S .+ AnmbtφS) .* dRIndr .* fnmb1 .+ AnmbrS .* (dRIndr .* dfnmb1dr .+ d2RIndr2 .* fnmb1) .+ AnmbθS .* dRIndr .* dfnmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* RIn .* fmbmb0 .+ AmbmbrS .* (dRIndr .* fmbmb0 .+ RIn .* dfmbmb0dr) .+ AmbmbθS .* RIn .* dfmbmb0dθ) .-
            ((Ambmb0S .+ AmbmbtφS) .* dRIndr .* fmbmb1 .+ AmbmbrS .* (d2RIndr2 .* fmbmb1 .+ dRIndr .* dfmbmb1dr) .+ AmbmbθS .* dRIndr .* dfmbmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* d2RIndr2 .* fmbmb2 .+ AmbmbrS .* (d3RIndr3 .* fmbmb2 .+ d2RIndr2 .* dfmbmb2dr) .+ AmbmbθS .* d2RIndr2 .* dfmbmb2dθ)) .* expθq
        ωdP = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* dRIndω .-
               (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* d2RIndrdω .+ Ammbmb .* fmbmb2 .* d3RIndr2dω .+
               (Amnn .* dfnn0dω .+ Amnmb .* dfnmb0dω .+ Ammbmb .* dfmbmb0dω) .* RIn .-
               (Amnmb .* dfnmb1dω .+ Ammbmb .* dfmbmb1dω) .* dRIndr .+ Ammbmb .* dfmbmb2dω .* d2RIndr2) .* ω1 .* expθq
        sumP1 += sum(Σ .* dipP .+ Σ .* ωdP)

        dipM = (
            ((Ann0S .+ AnntφS) .* RUp .* fnn0 .+ AnnrS .* (dRUpdr .* fnn0 .+ RUp .* dfnn0dr) .+ AnnθS .* RUp .* dfnn0dθ) .+
            ((Anmb0S .+ AnmbtφS) .* RUp .* fnmb0 .+ AnmbrS .* (RUp .* dfnmb0dr .+ dRUpdr .* fnmb0) .+ AnmbθS .* RUp .* dfnmb0dθ) .-
            ((Anmb0S .+ AnmbtφS) .* dRUpdr .* fnmb1 .+ AnmbrS .* (dRUpdr .* dfnmb1dr .+ d2RUpdr2 .* fnmb1) .+ AnmbθS .* dRUpdr .* dfnmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* RUp .* fmbmb0 .+ AmbmbrS .* (dRUpdr .* fmbmb0 .+ RUp .* dfmbmb0dr) .+ AmbmbθS .* RUp .* dfmbmb0dθ) .-
            ((Ambmb0S .+ AmbmbtφS) .* dRUpdr .* fmbmb1 .+ AmbmbrS .* (d2RUpdr2 .* fmbmb1 .+ dRUpdr .* dfmbmb1dr) .+ AmbmbθS .* dRUpdr .* dfmbmb1dθ) .+
            ((Ambmb0S .+ AmbmbtφS) .* d2RUpdr2 .* fmbmb2 .+ AmbmbrS .* (d3RUpdr3 .* fmbmb2 .+ d2RUpdr2 .* dfmbmb2dr) .+ AmbmbθS .* d2RUpdr2 .* dfmbmb2dθ)) .* expθq
        ωdM = ((Amnn .* fnn0 .+ Amnmb .* fnmb0 .+ Ammbmb .* fmbmb0) .* dRUpdω .-
               (Amnmb .* fnmb1 .+ Ammbmb .* fmbmb1) .* d2RUpdrdω .+ Ammbmb .* fmbmb2 .* d3RUpdr2dω .+
               (Amnn .* dfnn0dω .+ Amnmb .* dfnmb0dω .+ Ammbmb .* dfmbmb0dω) .* RUp .-
               (Amnmb .* dfnmb1dω .+ Ammbmb .* dfmbmb1dω) .* dRUpdr .+ Ammbmb .* dfmbmb2dω .* d2RUpdr2) .* ω1 .* expθq
        sumM1 += sum(Σ .* dipM .+ Σ .* ωdM)

        # --- 3. r-derivative sums (.wl L3982-3990 In; L3991-3999 Up) ---
        dsumPdr += sum(Σ .* (
            (dΣdr ./ Σ .+ dexpdr) .* (Amnn .* fnn0 .* RIn .+ Amnmb .* (fnmb0 .* RIn .- fnmb1 .* dRIndr) .+ Ammbmb .* (fmbmb0 .* RIn .- fmbmb1 .* dRIndr .+ fmbmb2 .* d2RIndr2)) .+
            (2 .* vn .* dvndr .* fnn0 .* RIn .+ (vn .* dvmbdr .+ dvndr .* vmb) .* (fnmb0 .* RIn .- fnmb1 .* dRIndr) .+ 2 .* vmb .* dvmbdr .* (fmbmb0 .* RIn .- fmbmb1 .* dRIndr .+ fmbmb2 .* d2RIndr2)) .+
            (Amnn .* ((dfnn0dr .+ dfnn0dθ .* dθpdr .+ dfnn0dω .* dωdr) .* RIn .+ fnn0 .* (dRIndr .+ dRIndω .* dωdr)) .+
             Amnmb .* ((dfnmb0dr .+ dfnmb0dθ .* dθpdr .+ dfnmb0dω .* dωdr) .* RIn .+ fnmb0 .* (dRIndr .+ dRIndω .* dωdr)
                  .- (dfnmb1dr .+ dfnmb1dθ .* dθpdr .+ dfnmb1dω .* dωdr) .* dRIndr .- fnmb1 .* (d2RIndr2 .+ d2RIndrdω .* dωdr)) .+
             Ammbmb .* ((dfmbmb0dr .+ dfmbmb0dθ .* dθpdr .+ dfmbmb0dω .* dωdr) .* RIn .+ fmbmb0 .* (dRIndr .+ dRIndω .* dωdr)
                  .- (dfmbmb1dr .+ dfmbmb1dθ .* dθpdr .+ dfmbmb1dω .* dωdr) .* dRIndr .- fmbmb1 .* (d2RIndr2 .+ d2RIndrdω .* dωdr)
                  .+ (dfmbmb2dr .+ dfmbmb2dθ .* dθpdr .+ dfmbmb2dω .* dωdr) .* d2RIndr2 .+ fmbmb2 .* (d3RIndr3 .+ d3RIndr2dω .* dωdr)))
            ) .* expθq)
        dsumMdr += sum(Σ .* (
            (dΣdr ./ Σ .+ dexpdr) .* (Amnn .* fnn0 .* RUp .+ Amnmb .* (fnmb0 .* RUp .- fnmb1 .* dRUpdr) .+ Ammbmb .* (fmbmb0 .* RUp .- fmbmb1 .* dRUpdr .+ fmbmb2 .* d2RUpdr2)) .+
            (2 .* vn .* dvndr .* fnn0 .* RUp .+ (vn .* dvmbdr .+ dvndr .* vmb) .* (fnmb0 .* RUp .- fnmb1 .* dRUpdr) .+ 2 .* vmb .* dvmbdr .* (fmbmb0 .* RUp .- fmbmb1 .* dRUpdr .+ fmbmb2 .* d2RUpdr2)) .+
            (Amnn .* ((dfnn0dr .+ dfnn0dθ .* dθpdr .+ dfnn0dω .* dωdr) .* RUp .+ fnn0 .* (dRUpdr .+ dRUpdω .* dωdr)) .+
             Amnmb .* ((dfnmb0dr .+ dfnmb0dθ .* dθpdr .+ dfnmb0dω .* dωdr) .* RUp .+ fnmb0 .* (dRUpdr .+ dRUpdω .* dωdr)
                  .- (dfnmb1dr .+ dfnmb1dθ .* dθpdr .+ dfnmb1dω .* dωdr) .* dRUpdr .- fnmb1 .* (d2RUpdr2 .+ d2RUpdrdω .* dωdr)) .+
             Ammbmb .* ((dfmbmb0dr .+ dfmbmb0dθ .* dθpdr .+ dfmbmb0dω .* dωdr) .* RUp .+ fmbmb0 .* (dRUpdr .+ dRUpdω .* dωdr)
                  .- (dfmbmb1dr .+ dfmbmb1dθ .* dθpdr .+ dfmbmb1dω .* dωdr) .* dRUpdr .- fmbmb1 .* (d2RUpdr2 .+ d2RUpdrdω .* dωdr)
                  .+ (dfmbmb2dr .+ dfmbmb2dθ .* dθpdr .+ dfmbmb2dω .* dωdr) .* d2RUpdr2 .+ fmbmb2 .* (d3RUpdr3 .+ d3RUpdr2dω .* dωdr)))
            ) .* expθq)

        # --- 4. x-derivative sums (.wl L4000-4008 In; L4009-4017 Up). NB no dfXdr term (dfdx=dfdθ·dθdx+dfdω·dωdx). ---
        dsumPdx += sum(Σ .* (
            (dΣdx ./ Σ .+ dexpdx) .* (Amnn .* fnn0 .* RIn .+ Amnmb .* (fnmb0 .* RIn .- fnmb1 .* dRIndr) .+ Ammbmb .* (fmbmb0 .* RIn .- fmbmb1 .* dRIndr .+ fmbmb2 .* d2RIndr2)) .+
            (2 .* vn .* dvndx .* fnn0 .* RIn .+ (vn .* dvmbdx .+ dvndx .* vmb) .* (fnmb0 .* RIn .- fnmb1 .* dRIndr) .+ 2 .* vmb .* dvmbdx .* (fmbmb0 .* RIn .- fmbmb1 .* dRIndr .+ fmbmb2 .* d2RIndr2)) .+
            (Amnn .* ((dfnn0dθ .* dθpdx .+ dfnn0dω .* dωdx) .* RIn .+ fnn0 .* (dRIndω .* dωdx)) .+
             Amnmb .* ((dfnmb0dθ .* dθpdx .+ dfnmb0dω .* dωdx) .* RIn .+ fnmb0 .* (dRIndω .* dωdx)
                  .- (dfnmb1dθ .* dθpdx .+ dfnmb1dω .* dωdx) .* dRIndr .- fnmb1 .* (d2RIndrdω .* dωdx)) .+
             Ammbmb .* ((dfmbmb0dθ .* dθpdx .+ dfmbmb0dω .* dωdx) .* RIn .+ fmbmb0 .* (dRIndω .* dωdx)
                  .- (dfmbmb1dθ .* dθpdx .+ dfmbmb1dω .* dωdx) .* dRIndr .- fmbmb1 .* (d2RIndrdω .* dωdx)
                  .+ (dfmbmb2dθ .* dθpdx .+ dfmbmb2dω .* dωdx) .* d2RIndr2 .+ fmbmb2 .* (d3RIndr2dω .* dωdx)))
            ) .* expθq)
        dsumMdx += sum(Σ .* (
            (dΣdx ./ Σ .+ dexpdx) .* (Amnn .* fnn0 .* RUp .+ Amnmb .* (fnmb0 .* RUp .- fnmb1 .* dRUpdr) .+ Ammbmb .* (fmbmb0 .* RUp .- fmbmb1 .* dRUpdr .+ fmbmb2 .* d2RUpdr2)) .+
            (2 .* vn .* dvndx .* fnn0 .* RUp .+ (vn .* dvmbdx .+ dvndx .* vmb) .* (fnmb0 .* RUp .- fnmb1 .* dRUpdr) .+ 2 .* vmb .* dvmbdx .* (fmbmb0 .* RUp .- fmbmb1 .* dRUpdr .+ fmbmb2 .* d2RUpdr2)) .+
            (Amnn .* ((dfnn0dθ .* dθpdx .+ dfnn0dω .* dωdx) .* RUp .+ fnn0 .* (dRUpdω .* dωdx)) .+
             Amnmb .* ((dfnmb0dθ .* dθpdx .+ dfnmb0dω .* dωdx) .* RUp .+ fnmb0 .* (dRUpdω .* dωdx)
                  .- (dfnmb1dθ .* dθpdx .+ dfnmb1dω .* dωdx) .* dRUpdr .- fnmb1 .* (d2RUpdrdω .* dωdx)) .+
             Ammbmb .* ((dfmbmb0dθ .* dθpdx .+ dfmbmb0dω .* dωdx) .* RUp .+ fmbmb0 .* (dRUpdω .* dωdx)
                  .- (dfmbmb1dθ .* dθpdx .+ dfmbmb1dω .* dωdx) .* dRUpdr .- fmbmb1 .* (d2RUpdrdω .* dωdx)
                  .+ (dfmbmb2dθ .* dθpdx .+ dfmbmb2dω .* dωdx) .* d2RUpdr2 .+ fmbmb2 .* (d3RUpdr2dω .* dωdx)))
            ) .* expθq)
    end

    # invariant Wronskian + its ω-derivative (.wl L4021-4022)
    W  = (RIn * dRUpdr - dRIndr * RUp) / Δ
    dWdω = (dRIndω * dRUpdr + RIn * d2RUpdrdω - d2RIndrdω * RUp - dRIndr * dRUpdω) / Δ
    den = Γ * W * stepsθ
    CPlus0  = -2pi * sumP0 / den
    CMinus0 = -2pi * sumM0 / den
    # linear-in-spin amplitude corrections (.wl L4025-4026)
    CPlus1  = -2pi * (sumP1  - Γ1 / Γ * sumP0  - dWdω * ω1 / W * sumP0)  / den
    CMinus1 = -2pi * (sumM1  - Γ1 / Γ * sumM0  - dWdω * ω1 / W * sumM0)  / den
    # r/x amplitude derivatives (.wl L4027-4030)
    dCPlusdr  = -2pi * (dsumPdr - dΓdr / Γ * sumP0 - dWdω * dωdr / W * sumP0) / den
    dCMinusdr = -2pi * (dsumMdr - dΓdr / Γ * sumM0 - dWdω * dωdr / W * sumM0) / den
    dCPlusdx  = -2pi * (dsumPdx - dΓdx / Γ * sumP0 - dWdω * dωdx / W * sumP0) / den
    dCMinusdx = -2pi * (dsumMdx - dΓdx / Γ * sumM0 - dWdω * dωdx / W * sumM0) / den

    # Fluxes (.wl L4064-4072)
    EnergyI = abs2(CPlus0) / (4pi * ω^2);  EnergyH = α * abs2(CMinus0) / (4pi * ω^2)
    AngMomI = abs2(CPlus0) * m / (4pi * ω^3); AngMomH = α * abs2(CMinus0) * m / (4pi * ω^3)
    # FluxesCorrection dE¹/dLz¹ (.wl L4076-4081)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (2 * real(CMinus1 * conj(CMinus0)) + (dαdω / α * abs2(CMinus0) - 2 * abs2(CMinus0) / ω) * ω1) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    # FluxesDerivatives r (.wl L4087-4092)
    rePr = 2 * real(dCPlusdr * conj(CPlus0)); reMr = 2 * real(dCMinusdr * conj(CMinus0))
    dEdrI = (rePr - 2 * abs2(CPlus0) * dωdr / ω) / (4pi * ω^2)
    dEdrH = α * (dαdω * dωdr / α * abs2(CMinus0) + reMr - 2 * abs2(CMinus0) * dωdr / ω) / (4pi * ω^2)
    dLdrI = (rePr - 3 * abs2(CPlus0) * dωdr / ω) * m / (4pi * ω^3)
    dLdrH = α * (dαdω * dωdr / α * abs2(CMinus0) + reMr - 3 * abs2(CMinus0) * dωdr / ω) * m / (4pi * ω^3)
    # FluxesDerivatives x (.wl L4097-4102)
    rePx = 2 * real(dCPlusdx * conj(CPlus0)); reMx = 2 * real(dCMinusdx * conj(CMinus0))
    dEdxI = (rePx - 2 * abs2(CPlus0) * dωdx / ω) / (4pi * ω^2)
    dEdxH = α * (dαdω * dωdx / α * abs2(CMinus0) + reMx - 2 * abs2(CMinus0) * dωdx / ω) / (4pi * ω^2)
    dLdxI = (rePx - 3 * abs2(CPlus0) * dωdx / ω) * m / (4pi * ω^3)
    dLdxH = α * (dαdω * dωdx / α * abs2(CMinus0) + reMx - 3 * abs2(CMinus0) * dωdx / ω) * m / (4pi * ω^3)

    return (l = l, m = m, k = k, n = 0, omega = ω, omega1 = ω1, dωdr = dωdr, dωdx = dωdx,
            CPlus0 = CPlus0, CMinus0 = CMinus0, CPlus1 = CPlus1, CMinus1 = CMinus1,
            dCPlusdr = dCPlusdr, dCMinusdr = dCMinusdr, dCPlusdx = dCPlusdx, dCMinusdx = dCMinusdx,
            α = α, dαdω = dαdω, Spi2 = SWSH(pi / 2, 0.0),
            EnergyI = EnergyI, EnergyH = EnergyH, AngMomI = AngMomI, AngMomH = AngMomH,
            dE1I = dE1I, dE1H = dE1H, dL1I = dL1I, dL1H = dL1H,
            dEdrI = dEdrI, dEdrH = dEdrH, dLdrI = dLdrI, dLdrH = dLdrH,
            dEdxI = dEdxI, dEdxH = dEdxH, dLdxI = dLdxI, dLdxH = dLdxH,
            λ = λ, dλdω = dλdω, stepsθ = stepsθ)
end

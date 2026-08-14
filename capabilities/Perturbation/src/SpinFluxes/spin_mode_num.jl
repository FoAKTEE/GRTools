# Numeric-orbit linear-in-spin Teukolsky mode — TeukolskySpinFluxes.wl TeukolskySpinMode (L1210-1404).
# This is the FULL finite-s sibling of TeukolskySpinModeFromCorrection: given a materialized numeric
# spinning orbit (SpinningOrbit`KerrSpinOrbitNum[orbitCorrection,s], which bakes the spin s into the
# trajectory + BL/Mino frequencies) it runs the SAME wr×wθ Mino double-sum, but evaluates the source
# at the PERTURBED point rp=r0+s·r1, zp=z0+s·z1 with the radial function Taylor-expanded about r0
# (RInrp+s·r1·dRInrp …) and the FULL source assembly A0/A1/A2/B1/B2/B3 + C0nn/Crnn/Cθnn (NOT the
# linearised-in-s monopole+dipole split of FromCorrection). The fab^(i) basis (fabi/dfabidr/dfabidθ),
# the angular layer (_angular_fromcorr), h1/h2/h3 (_fc_h*), _alpha_circ and _dδ are SHARED with
# from_correction.jl (same module). CPlus/CMinus → ℐ/ℋ Energy/AngularMomentum fluxes. R3.
#
# Difference from FromCorrection in the normalisation: TeukolskySpinMode uses Γ = orbit["MinoFrequencies"][1]
# (already the spin-corrected Γ̂+s·Γ̂1) directly, so there is NO (1−sΓ1/Γ) linearisation factor.

"""
    _spin_mode_num(orbit, l, m, n, k) -> NamedTuple

Core port of TeukolskySpinMode (.wl L1210-1404). `orbit` is the SpinningOrbit `KerrSpinOrbitNum`
NamedTuple (fields `a,p,e,s,Ehat,Lzhat,Khat,ES,LS,BLFrequencies,MinoFrequencies,TrajectoryGeo,
TrajectoryDeltas,OrbitCorrection`). Returns `(; l,m,n,k,omega, CPlus, CMinus, α, Spi2, EnergyI,
EnergyH, AngMomI, AngMomH, stepsr, stepsθ)`. R3.
"""
function _spin_mode_num(orbit, l::Integer, m::Integer, n::Integer, k::Integer)
    (l < 2 || abs(m) > l) && error("TeukolskySpinMode: l<2 or |m|>l")
    a = orbit.a; s = orbit.s
    En0 = orbit.Ehat; Lz0 = orbit.Lzhat; Kc0 = orbit.Khat   # geodesic constants of motion
    En1 = orbit.ES; Lz1 = orbit.LS                          # linear corrections to constants of motion
    correction = orbit.OrbitCorrection                       # (wr,wz) -> (ΔtS,rS,zS,ΔφS,UrS,UzS)
    Ωr, Ωθ, Ωφ = orbit.BLFrequencies                         # spin-corrected BL frequencies
    rfun = orbit.TrajectoryGeo[1]; θfun = orbit.TrajectoryGeo[2]
    rgeo(wr) = rfun(wr).r
    zgeo(wz) = cos(θfun(wz).θ)
    Γ = orbit.MinoFrequencies[1]                             # spin-corrected Γ (used directly)
    ω = m * Ωφ + n * Ωr + k * Ωθ
    isfinite(ω) || error("TeukolskySpinMode: non-real ω")

    SWSH = SpinWeightedSpheroidalHarmonicS(-2, l, m, a * ω)
    λ = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ω)
    ν = RenormalizedAngularMomentum(-2, l, m, a, ω, λ)
    R = TeukolskyRadial(-2, l, m, a, ω; lambda = λ, nu = ν)
    α = _alpha_circ(a, m, ω, λ)

    # step counts (.wl L1251-1254): 32·⌈|ωΔt'−mΔφ'+idx|⌉ at the turning points, ≥32.
    Δtr = orbit.TrajectoryDeltas.Δtr; Δφr = orbit.TrajectoryDeltas.Δφr
    Δtθ = orbit.TrajectoryDeltas.Δtz; Δφθ = orbit.TrajectoryDeltas.Δφz
    stepsr = max(32 * ceil(Int, abs(ω * _dδ(Δtr, pi)  - m * _dδ(Δφr, pi)  + n)),
                 32 * ceil(Int, abs(ω * _dδ(Δtr, 0.0) - m * _dδ(Δφr, 0.0) + n)), 32)
    stepsθ = max(32 * ceil(Int, abs(ω * _dδ(Δtθ, pi/4) - m * _dδ(Δφθ, pi/4) + k)),
                 32 * ceil(Int, abs(ω * _dδ(Δtθ, 0.0)  - m * _dδ(Δφθ, 0.0)  + k)), 32)
    npr = stepsr ÷ 2; npθ = stepsθ ÷ 2

    # precompute trajectory corrections at all (ir,iθ) midpoints (.wl `correctionp`).
    correctionp = [[correction((ir - 1//2) * 2pi / stepsr, (iθ - 1//2) * 2pi / stepsθ)
                    for iθ in 1:npθ] for ir in 1:npr]

    θc = Vector{NamedTuple}(undef, npθ)                      # θ-only quantities, cached on first ir pass
    sumP = zero(ComplexF64); sumM = zero(ComplexF64)
    RInrp = zero(ComplexF64); dRInrp = zero(ComplexF64); RUprp = zero(ComplexF64); dRUprp = zero(ComplexF64)
    ddRInrp = zero(ComplexF64); dddRInrp = zero(ComplexF64); ddRUprp = zero(ComplexF64); dddRUprp = zero(ComplexF64)
    r0p_last = 0.0                                           # last geodesic r0p for the invariant Wronskian

    for ir in 1:npr
        wr = (ir - 1//2) * 2pi / stepsr
        r0p = rgeo(wr); r0p_last = r0p
        sqrtR = sqrt(((r0p^2 + a^2) * En0 - a * Lz0)^2 - (r0p^2 - 2r0p + a^2) * (r0p^2 + Kc0))
        Ur0 = (1, -1, 1, -1) .* sqrtR                        # geodesic radial velocity per quadrant
        expr = exp(im * (ω * Δtr(wr) - m * Δφr(wr) + 2pi * n * (ir - 1//2) / stepsr))
        exprq = (expr, 1 / expr, expr, 1 / expr)             # expr^{1,-1,1,-1}
        Δ0 = r0p^2 - 2r0p + a^2
        K0 = (r0p^2 + a^2) * ω - a * m
        dΔ0 = 2 * (r0p - 1)
        dK0 = 2r0p * ω
        V  = -(K0^2 + 4im * (r0p - 1) * K0) / Δ0 + 8im * ω * r0p + λ
        dV = -((2K0 * dK0 + 4im * K0 + 4im * (r0p - 1) * dK0) * Δ0 - (K0^2 + 4im * (r0p - 1) * K0) * dΔ0) / Δ0^2 + 8im * ω
        RInrp = R.In(r0p); dRInrp = teukolsky_radial_deriv(R.In, r0p)
        ddRInrp  = (V * RInrp + dΔ0 * dRInrp) / Δ0
        dddRInrp = 1 / Δ0 * (dV * RInrp + (V + 2) * dRInrp)
        RUprp = R.Up(r0p); dRUprp = teukolsky_radial_deriv(R.Up, r0p)
        ddRUprp  = (V * RUprp + dΔ0 * dRUprp) / Δ0
        dddRUprp = 1 / Δ0 * (dV * RUprp + (V + 2) * dRUprp)

        for iθ in 1:npθ
            if ir == 1
                wθ = (iθ - 1//2) * 2pi / stepsθ
                z0p = zgeo(wθ); sinθ0p = sqrt(1 - z0p^2)
                sqrtZ = sqrt(-((1 - z0p^2) * a * En0 - Lz0)^2 + (1 - z0p^2) * (Kc0 - a^2 * z0p^2))
                Uz0 = (1, 1, -1, -1) .* (-sqrtZ)             # geodesic polar velocity per quadrant
                expθ = exp(im * (ω * Δtθ(wθ) - m * Δφθ(wθ) + 2pi * k * (iθ - 1//2) / stepsθ))
                expθq = (expθ, expθ, 1 / expθ, 1 / expθ)     # expθ^{1,1,-1,-1}
                ang = _angular_fromcorr(SWSH, λ, a, m, ω, z0p)
                θc[iθ] = (z0p = z0p, sinθ0p = sinθ0p, Uz0 = Uz0, expθq = expθq,
                          S = ang.S, dSdθ = ang.dSdθ, L2S = ang.L2S, dL2Sdθ = ang.dL2Sdθ,
                          L1L2S = ang.L1L2S, dL1L2Sdθ = ang.dL1L2Sdθ)
            end
            t = θc[iθ]
            z0p = t.z0p; sinθ0p = t.sinθ0p; Uz0 = t.Uz0; expθq = t.expθq
            S = t.S; dSdθ = t.dSdθ; L2S = t.L2S; dL2Sdθ = t.dL2Sdθ; L1L2S = t.L1L2S; dL1L2Sdθ = t.dL1L2Sdθ

            # fab^(i) basis at the GEODESIC point (the trajectory perturbation enters via R-Taylor + C-coeffs)
            ζ0 = r0p - im * a * z0p; ζ0bar = r0p + im * a * z0p
            (fnn0, fnmb0, fnmb1, fmbmb0, fmbmb1, fmbmb2) =
                _fabi(ζ0, ζ0bar, sinθ0p, Δ0, dΔ0, K0, dK0, S, L2S, L1L2S, a)
            (dfnn0dr, dfnmb0dr, dfnmb1dr, dfmbmb0dr, dfmbmb1dr, dfmbmb2dr) =
                _dfabidr(ζ0, ζ0bar, sinθ0p, Δ0, dΔ0, K0, dK0, S, L2S, L1L2S, a, ω)
            (dfnn0dθ, dfnmb0dθ, dfnmb1dθ, dfmbmb0dθ, dfmbmb1dθ, dfmbmb2dθ) =
                _dfabidθ(ζ0, ζ0bar, sinθ0p, Δ0, dΔ0, K0, dK0, S, L2S, L1L2S, dSdθ, dL2Sdθ, dL1L2Sdθ, a)

            # trajectory correction at the 4 quadrants (.wl L1302-1305). cp/cm = +iθ / -iθ columns.
            cp = correctionp[ir][iθ]
            cm = correctionp[ir][npθ - iθ + 1]
            r1p  = (cp.rS,  cm.rS,  cm.rS,  cp.rS)
            z1p  = (cp.zS,  -cm.zS, -cm.zS, cp.zS)
            Ur1p = (cp.UrS, -cm.UrS, cm.UrS, -cp.UrS)
            Uz1p = (cp.UzS, cm.UzS, -cm.UzS, -cp.UzS)

            # FULL perturbed trajectory + four-velocity (4-vectors over quadrants)
            rp = r0p .+ s .* r1p
            zp = z0p .+ s .* z1p
            sinθp = sqrt.(1 .- zp .^ 2)
            Urp = Ur0 .+ s .* Ur1p
            Uzp = Uz0 .+ s .* Uz1p
            h1v = _fc_h1(r0p, z0p, Kc0, a); h2v = _fc_h2(r0p, z0p, Kc0, En0, Lz0, a); h3v = _fc_h3(r0p, z0p, Kc0, a)
            udtp = -(En0 + s * En1) + s * h1v                # scalar (covariant t-component)
            udφp = (Lz0 + s * Lz1) .+ s .* (h2v .+ h3v .* Ur0 .* Uz0)  # 4-vector (covariant φ-component)

            ζ = rp .- im * a .* zp; ζbar = rp .+ im * a .* zp
            Kp = (rp .^ 2 .+ a^2) .* ω .- a * m
            Δ  = rp .^ 2 .- 2 .* rp .+ a^2
            Σ  = rp .^ 2 .+ a^2 .* zp .^ 2
            exp1 = (exp(im * s * (ω * cp.ΔtS - m * cp.ΔφS)), exp(-im * s * (ω * cm.ΔtS - m * cm.ΔφS)),
                    exp(im * s * (ω * cm.ΔtS - m * cm.ΔφS)), exp(-im * s * (ω * cp.ΔtS - m * cp.ΔφS)))

            # four-velocity in Kinnersley tetrad (full perturbed values)
            vn  = ((rp .^ 2 .+ a^2) .* udtp .+ a .* udφp .- Urp) ./ (2 .* Σ)
            vl  = ((rp .^ 2 .+ a^2) .* udtp .+ a .* udφp .+ Urp) ./ Δ
            vm  = (-im .* (a .* sinθp .^ 2 .* udtp .+ udφp) .+ Uzp) ./ (-sqrt(2) .* sinθp .* ζbar)
            vmb = conj.(vm)

            # spin tensor in Kinnersley tetrad (carries explicit s; full perturbed Σ,ζ).
            # NB: ζ/ζbar/Σ/Δ/Kp and the NP coeffs are 4-vectors here (perturbed rp,zp per quadrant),
            # so unary minus must broadcast: `(.-x)` (Julia has no unary `-` on tuples).
            Sln  = s .* ((.-(rp .* (Kc0 .- a^2 .* zp .^ 2))) ./ (sqrt(Kc0) .* Σ))
            Snm  = s .* (ζ ./ sqrt(Kc0)) .* vm .* vn
            Snmb = conj.(Snm)
            Slmb = s .* (.-(ζ ./ sqrt(Kc0))) .* vl .* vmb
            Smmb = s .* ((im * a .* zp .* (Kc0 .+ rp .^ 2)) ./ (sqrt(Kc0) .* Σ))
            Amnn   = vn .^ 2
            Amnmb  = vn .* vmb
            Ammbmb = vmb .^ 2

            # Newman-Penrose spin coefficients (perturbed; opposite sign convention vs Teukolsky paper)
            ρ  = 1 ./ ζ
            βc = .-(zp ./ (2 .* ζbar .* sqrt(2) .* sinθp))
            πc = .-((im * a .* sinθp) ./ (ζ .^ 2 .* sqrt(2)))
            τc = (im * a .* sinθp) ./ (sqrt(2) .* Σ)
            μc = Δ ./ (2 .* ζ .^ 2 .* ζbar)
            γc = (a^2 .- rp .+ im * a .* (-1 .+ rp) .* zp) ./ (2 .* ζ .^ 2 .* ζbar)
            αc = .-(((.-(rp .* zp)) .- im * a .* (-2 .+ zp .^ 2)) ./ (2 .* ζ .^ 2 .* sqrt(2) .* sinθp))

            Scdγndc  = (.-Sln) .* 2 .* real.(γc) .- Smmb .* ((.-μc) .+ conj.(μc))
            Scdγmbdc = (.-Sln) .* (πc .+ conj.(τc)) .- Snmb .* conj.(ρ) .- Slmb .* ((.-conj.(γc)) .+ γc .- μc) .- Smmb .* ((.-αc) .+ conj.(βc))
            Adnn   = (Scdγndc .* vn .- Sln .* 2 .* real.(γc) .* vn .- 2 .* real.(Snmb .* ((conj.(αc) .+ βc) .* vn .- μc .* vm)))
            Admbmb = (Scdγmbdc .* vmb .- Snmb .* ((.-πc) .* vl) .- Slmb .* (conj.(τc) .* vn .- (conj.(γc) .- γc) .* vmb) .+ Smmb .* ((.-((.-αc) .+ conj.(βc))) .* vmb))
            Adnmb  = (Scdγndc .* vmb .+ Scdγmbdc .* vn .- Sln .* (conj.(τc) .* vn .- (conj.(γc) .- γc) .* vmb) .-
                      Snmb .* (conj.(ρ) .* vn .- μc .* vl .- (conj.(αc) .- βc) .* vmb) .- Snm .* ((.-((.-αc) .+ conj.(βc))) .* vmb) .-
                      Snmb .* ((.-conj.(πc)) .* vmb .- πc .* vm) .- Slmb .* (2 .* real.(γc) .* vn) .+
                      Smmb .* ((αc .+ conj.(βc)) .* vn .- conj.(μc) .* vmb)) ./ 2

            # functions from Eq. (B3) in 2303.16798 (perturbed Kp,Δ,Σ,ζ)
            Stφn  = (-im) .* Kp ./ (2 .* Σ) .* Sln .+ (a * ω .* sinθp .- m ./ sinθp) ./ (sqrt(2) .* Σ) .* (ζ .* Snmb .- ζbar .* Snm)
            Stφmb = (-im) .* Kp .* (1 ./ Δ .* Snmb .+ 1 ./ (2 .* Σ) .* Slmb) .+ (a * ω .* sinθp .- m ./ sinθp) ./ (sqrt(2) .* ζ) .* Smmb
            Srn   = Δ ./ (2 .* Σ) .* Sln
            Srmb  = (.-Snmb) .+ Δ ./ (2 .* Σ) .* Slmb
            Sθn   = (.-(Snmb .* ζ .+ Snm .* ζbar)) ./ (sqrt(2) .* Σ)
            Sθmb  = Smmb ./ (sqrt(2) .* ζ)

            # FULL source C-coefficients (.wl L1354-1362; z1p divided by GEODESIC sinθ0p)
            C0nn   = Amnn .+ Adnn .+ Stφn .* vn
            Crnn   = Amnn .* s .* r1p .+ Srn .* vn
            Cθnn   = (.-Amnn) .* s .* z1p ./ sinθ0p .+ Sθn .* vn
            C0nmb  = Amnmb .+ Adnmb .+ (Stφn .* vmb .+ Stφmb .* vn) ./ 2
            Crnmb  = Amnmb .* s .* r1p .+ (Srn .* vmb .+ Srmb .* vn) ./ 2
            Cθnmb  = (.-Amnmb) .* s .* z1p ./ sinθ0p .+ (Sθn .* vmb .+ Sθmb .* vn) ./ 2
            C0mbmb = Ammbmb .+ Admbmb .+ Stφmb .* vmb
            Crmbmb = Ammbmb .* s .* r1p .+ Srmb .* vmb
            Cθmbmb = (.-Ammbmb) .* s .* z1p ./ sinθ0p .+ Sθmb .* vmb

            A0 = C0nn .* fnn0 .+ Crnn .* dfnn0dr .+ Cθnn .* dfnn0dθ .+ C0nmb .* fnmb0 .+ Crnmb .* dfnmb0dr .+
                 Cθnmb .* dfnmb0dθ .+ C0mbmb .* fmbmb0 .+ Crmbmb .* dfmbmb0dr .+ Cθmbmb .* dfmbmb0dθ
            A1 = C0nmb .* fnmb1 .+ Crnmb .* dfnmb1dr .+ Cθnmb .* dfnmb1dθ .+ C0mbmb .* fmbmb1 .+ Crmbmb .* dfmbmb1dr .+ Cθmbmb .* dfmbmb1dθ
            A2 = C0mbmb .* fmbmb2 .+ Crmbmb .* dfmbmb2dr .+ Cθmbmb .* dfmbmb2dθ
            B1 = (.-Srn) .* vn .* fnn0 .- (Srn .* vmb .+ Srmb .* vn) ./ 2 .* fnmb0 .- Srmb .* vmb .* fmbmb0
            B2 = (.-(Srn .* vmb .+ Srmb .* vn)) ./ 2 .* fnmb1 .- Srmb .* vmb .* fmbmb1
            B3 = (.-Srmb) .* vmb .* fmbmb2

            # quadrant sum with R Taylor-expanded about r0p (RInrp + s·r1p·dRInrp, …)
            integP = Σ .* ((A0 .* (RInrp .+ s .* r1p .* dRInrp) .- (A1 .+ B1) .* (dRInrp .+ s .* r1p .* ddRInrp) .+
                            (A2 .+ B2) .* (ddRInrp .+ s .* r1p .* dddRInrp) .- B3 .* dddRInrp) .* exprq .* expθq .* exp1)
            sumP += sum(integP)
            integM = Σ .* ((A0 .* (RUprp .+ s .* r1p .* dRUprp) .- (A1 .+ B1) .* (dRUprp .+ s .* r1p .* ddRUprp) .+
                            (A2 .+ B2) .* (ddRUprp .+ s .* r1p .* dddRUprp) .- B3 .* dddRUprp) .* exprq .* expθq .* exp1)
            sumM += sum(integM)
        end
    end

    # invariant Wronskian (.wl L1373, at the last geodesic r0p)
    W = (RInrp * dRUprp - dRInrp * RUprp) / (r0p_last^2 - 2 * r0p_last + a^2)
    pre = -2pi / (Γ * W * stepsr * stepsθ)
    CPlus = pre * sumP; CMinus = pre * sumM
    return (l = l, m = m, n = n, k = k, omega = ω, α = α, CPlus = CPlus, CMinus = CMinus,
            Spi2 = SWSH(pi / 2, 0.0),
            EnergyI = abs2(CPlus) / (4pi * ω^2), EnergyH = α * abs2(CMinus) / (4pi * ω^2),
            AngMomI = abs2(CPlus) * m / (4pi * ω^3), AngMomH = α * abs2(CMinus) * m / (4pi * ω^3),
            stepsr = stepsr, stepsθ = stepsθ)
end

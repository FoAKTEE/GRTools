# Linear-in-spin gravitational-wave flux, ECCENTRIC-EQUATORIAL s=-2 — the wr Mino-phase midpoint sum
# generalizing circular.jl to the radial libration with two velocity branches (Ur≷0). Directly reuses
# the circular helper basis (_angular_circular/_FabEq/_alpha_circ) + the Radial/AmplitudesFluxes
# machinery; reduces to circular_spin_flux as e→0. Geodesic orbit + trajectory from Spacetime
# (geodesic_trajectory_equatorial); ∂_ωF/∂_ωW via central FD (recompute λ,ν,R,S,angular at ω±h). M=1.

using Spacetime: geodesic_trajectory_equatorial

"Radial carrier bundle at frequency ω': λ, ν, R (TeukolskyRadial) + angular (S,L2S,L1L2S)."
function _eqflux_bundle(l, m, a, ω; Nmax)
    s = -2
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, a * ω)
    ν = RenormalizedAngularMomentum(s, l, m, a, ω, λ)
    R = TeukolskyRadial(s, l, m, a, ω; lambda=λ, nu=ν, Nmax=Nmax)
    ang = _angular_circular(s, l, m, a, ω; Nmax=Nmax)
    return (λ=λ, R=R, ang=ang)
end

"FabEq projections (Fnn,Fnmb,Fmbmb,Gn,Gmb) for one branch (In/Up) at radius rp, frequency ω."
function _eqflux_fab(Rbranch, rp, a, m, ω, λ, ang)
    Δ = rp^2 - 2rp + a^2; K = (rp^2 + a^2) * ω - a * m; dKdr = 2 * rp * ω; dΔdr = 2 * rp - 2
    Rv = Rbranch(rp); dRv = teukolsky_radial_deriv(Rbranch, rp); d2Rv = _d2R_ode(Rv, dRv, -2, m, a, rp, ω, λ)
    DR = dRv - im * K / Δ * Rv
    dDRdr = d2Rv - im * K / Δ * dRv - im * (dKdr * Δ - K * dΔdr) / Δ^2 * Rv
    DDR = dDRdr - im * K / Δ * DR
    return _FabEq(rp, a, Rv, DR, DDR, ang.S, ang.L2S, ang.L1L2S)
end

"Invariant Wronskian W = (RIn dRUp' − dRIn' RUp)/Δ at radius rp."
function _eqflux_W(R, rp, a)
    Δ = rp^2 - 2rp + a^2
    RIn = R.In(rp); dRIn = teukolsky_radial_deriv(R.In, rp)
    RUp = R.Up(rp); dRUp = teukolsky_radial_deriv(R.Up, rp)
    return (RIn * dRUp - dRIn * RUp) / Δ
end

"""
    equatorial_spin_flux(l, m, n, a, p, e; Nmax=80, h=3e-4) -> NamedTuple

Linear-in-spin GW flux correction (dE¹,dLz¹ at ℐ/ℋ) for an ECCENTRIC-EQUATORIAL s=-2 spinning
secondary (radial harmonic `n`). The wr midpoint sum (stepsr matched to the upstream heuristic) over
the two velocity branches; ∂_ωF/∂_ωW via central FD. Returns
`(; dE1I,dE1H,dL1I,dL1H, CPlus0,CMinus0, omega, stepsr)`.
"""
function equatorial_spin_flux(l::Integer, m::Integer, n::Integer, a::Real, p::Real, e::Real;
        Nmax::Integer=80, h::Real=3e-4)
    s = -2
    En0 = geodesic_energy(a, p, e, 1); Lz0 = geodesic_angular_momentum(a, p, e, 1)
    Q0 = carter_constant(a, p, e, 1); K0 = Q0 + (Lz0 - a * En0)^2; sK = sqrt(K0)
    fr = kerr_frequencies(a, p, e, 1; Time="Mino")          # Υr, Υθ, Υφ, Υt (geodesic)
    Uτ = proper_frequency_factor(a, p, e, 1)
    Ωr = fr.Ur / fr.Ut; Ωφ = fr.Uph / fr.Ut
    ω = m * Ωφ + n * Ωr
    ω1 = 3 * Uτ * ω / (2 * sK * fr.Ut)
    Υt1 = -3 * Uτ / (2 * sK)
    traj = geodesic_trajectory_equatorial(a, p, e)
    # stepsr matched to the upstream 2^(5+⌈log2(Δphase')⌉) heuristic (FD of the trajectory deltas)
    dphase(wr) = (t1 = traj(wr + 1e-6); t0 = traj(wr - 1e-6);
                  ω * (t1.Δtr - t0.Δtr) / 2e-6 - m * (t1.Δφr - t0.Δφr) / 2e-6 + n)
    stepterm(wr) = (d = abs(dphase(wr)); d < 1e-12 ? 0.0 : 2.0^(5 + ceil(log2(d))))
    stepsr = Int(max(stepterm(pi), stepterm(0.0), 32.0))                             # Max[...,32]
    b0 = _eqflux_bundle(l, m, a, ω; Nmax=Nmax)
    bp = _eqflux_bundle(l, m, a, ω + h; Nmax=Nmax)
    bm = _eqflux_bundle(l, m, a, ω - h; Nmax=Nmax)
    vmb = -im * (a * En0 - Lz0)
    contractF(vn, F) = vn^2 * F.Fnn + vn * vmb * F.Fnmb + vmb^2 * F.Fmbmb
    CT = Complex{typeof(En0 + Lz0 + float(ω) + float(a))}
    sumP0 = zero(CT); sumP1 = zero(CT); sumM0 = zero(CT); sumM1 = zero(CT)
    W = zero(CT); dWdω = zero(CT)
    npts = stepsr ÷ 2
    for i in 1:npts
        wr = (2i - 1) * pi / stepsr
        tj = traj(wr); rp = tj.r; Ur = tj.Ur
        vnp = -((rp^2 + a^2) * En0 - a * Lz0 + Ur) / (rp^2 - 2rp + a^2)
        vnm = -((rp^2 + a^2) * En0 - a * Lz0 - Ur) / (rp^2 - 2rp + a^2)
        expr  = exp(im * (ω * tj.Δtr - m * tj.Δφr + n * wr))
        expr1 = im * (ω1 * tj.Δtr + ω * (-3 / (2 * sK)) * tj.Δτr)
        FinP = _eqflux_fab(b0.R.In, rp, a, m, ω, b0.λ, b0.ang)
        FupP = _eqflux_fab(b0.R.Up, rp, a, m, ω, b0.λ, b0.ang)
        # ∂_ω F_ab COMPONENTS (velocity-independent) via central FD of the In/Up FabEq at rp
        dFab(Rfp, λp, angp, Rfm, λm, angm) = let
            Fp = _eqflux_fab(Rfp, rp, a, m, ω + h, λp, angp); Fm = _eqflux_fab(Rfm, rp, a, m, ω - h, λm, angm)
            (Fnn = (Fp.Fnn - Fm.Fnn) / (2h), Fnmb = (Fp.Fnmb - Fm.Fnmb) / (2h), Fmbmb = (Fp.Fmbmb - Fm.Fmbmb) / (2h))
        end
        dFin = dFab(bp.R.In, bp.λ, bp.ang, bm.R.In, bm.λ, bm.ang)
        dFup = dFab(bp.R.Up, bp.λ, bp.ang, bm.R.Up, bm.λ, bm.ang)
        # geodesic
        sumP0 += contractF(vnp, FinP) * expr + contractF(vnm, FinP) / expr
        sumM0 += contractF(vnp, FupP) * expr + contractF(vnm, FupP) / expr
        # linear-in-spin: G_a source + phase shift (expr1) + ∂_ωF·ω1 (contracted PER velocity branch)
        GPp = (vnp * FinP.Gn + vmb * FinP.Gmb) / sK; GPm = (vnm * FinP.Gn + vmb * FinP.Gmb) / sK
        GUp = (vnp * FupP.Gn + vmb * FupP.Gmb) / sK; GUm = (vnm * FupP.Gn + vmb * FupP.Gmb) / sK
        sumP1 += (GPp + contractF(vnp, FinP) * expr1 + contractF(vnp, dFin) * ω1) * expr +
                 (GPm + contractF(vnm, FinP) * (-expr1) + contractF(vnm, dFin) * ω1) / expr
        sumM1 += (GUp + contractF(vnp, FupP) * expr1 + contractF(vnp, dFup) * ω1) * expr +
                 (GUm + contractF(vnm, FupP) * (-expr1) + contractF(vnm, dFup) * ω1) / expr
        if i == 1
            rp1 = rp
            W = _eqflux_W(b0.R, rp1, a)
            dWdω = (_eqflux_W(bp.R, rp1, a) - _eqflux_W(bm.R, rp1, a)) / (2h)
        end
    end
    Upst = fr.Ut
    CPlus0  = -2pi * sumP0 / (Upst * W * stepsr); CMinus0 = -2pi * sumM0 / (Upst * W * stepsr)
    W1 = dWdω * ω1
    pref = En0 / (2 * sK) - Υt1 / Upst - W1 / W
    CPlus1  = -2pi * sumP1 / (Upst * W * stepsr) + pref * CPlus0
    CMinus1 = -2pi * sumM1 / (Upst * W * stepsr) + pref * CMinus0
    α = _alpha_circ(a, m, ω, b0.λ)
    dαdω = (_alpha_circ(a, m, ω + h, bp.λ) - _alpha_circ(a, m, ω - h, bm.λ)) / (2h)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4pi * ω^2)
    dE1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4pi * ω^3)
    dL1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4pi * ω^3)
    return (dE1I=dE1I, dE1H=dE1H, dL1I=dL1I, dL1H=dL1H, CPlus0=CPlus0, CMinus0=CMinus0,
            omega=ω, stepsr=stepsr)
end

"""
    TeukolskySpinModeEquatorialCorrectionAnalyticalPhases(l, m, n, a, p, e; Nmax=80, h=3e-4) -> NamedTuple

Eccentric-EQUATORIAL s=-2 spinning-secondary mode (BHPToolkit name). Returns the geodesic amplitudes/
fluxes + the linear-in-spin flux correction dE¹/dLz¹ via the wr Mino sum: `(; l,m,n,a,p,e,omega,
Amplitudes=(I,H), Fluxes=Energy/AngularMomentum→(I,H), FluxesCorrection=Energy/AngularMomentum→(I,H))`.
"""
function TeukolskySpinModeEquatorialCorrectionAnalyticalPhases(l::Integer, m::Integer, n::Integer,
        a::Real, p::Real, e::Real; Nmax::Integer=80, h::Real=3e-4)
    sf = equatorial_spin_flux(l, m, n, a, p, e; Nmax=Nmax, h=h)
    λ = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * sf.omega)
    eg = energy_flux(-2, l, m, a, sf.omega, λ, sf.CPlus0, sf.CMinus0)
    jg = angular_momentum_flux(-2, l, m, a, sf.omega, λ, sf.CPlus0, sf.CMinus0)
    return (l = l, m = m, n = n, a = float(a), p = float(p), e = float(e), omega = sf.omega,
            Amplitudes = (I = sf.CPlus0, H = sf.CMinus0),
            Fluxes = (Energy = (I = eg.I, H = eg.H), AngularMomentum = (I = jg.I, H = jg.H)),
            FluxesCorrection = (Energy = (I = sf.dE1I, H = sf.dE1H),
                                AngularMomentum = (I = sf.dL1I, H = sf.dL1H)))
end

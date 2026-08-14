# Linear-in-spin gravitational-wave flux, circular-equatorial s=-2 — the FabEq path. dE¹ = ω¹Z⁰ +
# ω⁰Z¹; Z¹ from the linear amplitude correction C¹ = -2π Σ¹/(Υt W) + (E/(2√K) - Υt¹/Υt -
# ∂_ωW·ω¹/W)·C⁰, with Σ¹ = (vn Gn + vmb Gmb)/√K + (Σ vᵃvᵇ ∂_ω F_ab)·ω¹. The ∂_ω F_ab, ∂_ω W, ∂_ω α
# are central finite differences (velocities fixed); G_a is the spin-dipole source. Orbit linear
# parts (ω¹, Υt¹) from Spacetime.SpinningSecondary; angular/radial from Harmonics/Radial. M=1.
#
# This is the LEAF of the SpinFluxes namespace: the shared `_alpha_circ`/`_FabEq`/`_angular_circular`
# helper basis is reused by eccentric/spherical/inclined files. Flat-included into Perturbation,
# following the Radial/AmplitudesFluxes pattern (no helper-name collisions).

using Spacetime.SpinningSecondary: spin_linear_frequencies, spin_mino_time_correction

"angular building blocks at θ=π/2 (S, L2S = ℒ†_2 S, L1L2S = ℒ†_1 ℒ†_2 S)."
function _angular_circular(s::Integer, l::Integer, m::Integer, a::Real, omega::Real; Nmax::Integer=80)
    theta0 = oftype(float(omega), pi) / 2; sin0 = sin(theta0); cos0 = cos(theta0)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms = Nmax > 60 ? :auto : Nmax)
    S0  = Sc(theta0, zero(theta0))
    dS0 = swsh_derivative(Sc, 1, 0, theta0, zero(theta0))
    d2S0 = swsh_derivative(Sc, 2, 0, theta0, zero(theta0))
    L1 = -m / sin0 + a * omega * sin0 + cos0 / sin0
    L2 = -m / sin0 + a * omega * sin0 + 2 * cos0 / sin0
    L2p = m * cos0 / sin0^2 + a * omega * cos0 - 2 / sin0^2
    L2S = dS0 + L2 * S0
    L1Sp = d2S0 + L1 * dS0
    L1L2S = L1Sp + L2p * S0 + L2 * dS0 + L1 * L2 * S0
    return (S = S0, L2S = L2S, L1L2S = L1L2S)
end

"FabEq: (R, 𝒟R, 𝒟²R, S, L2S, L1L2S) → (Fnn,Fnmb,Fmbmb,Gn,Gmb) — source projections + spin-dipole Gₐ."
function _FabEq(r, a, R, DR, DDR, S, L2S, L1L2S)
    Fnn   = -1//2 * (L1L2S - 2 * im * a * L2S / r) * R
    Fnmb  = L2S * DR - 2 / r * L2S * R
    Fmbmb = -1//2 * S * (DDR - 2 / r * DR)
    Gn    = 1//2 * im * DR * (a * L2S + im * L1L2S * r)
    Gmb   = -1//2 * L2S * (DR - DDR * r)
    return (Fnn = Fnn, Fnmb = Fnmb, Fmbmb = Fmbmb, Gn = Gn, Gmb = Gmb)
end

"""
    _fabW(l, m, a, p, omega; Nmax=80) -> (Fin, Fup, W)

The FabEq source projections (In/Up) and the radial Wronskian W at frequency `omega` for a
circular-equatorial s=-2 orbit at r0=p. Recomputes λ(ω), ν(ω), the radial R/𝒟R/𝒟²R and the angular
S/L2S/L1L2S — so finite-differencing this in ω gives the ∂_ω F_ab / ∂_ω W.
"""
function _fabW(l::Integer, m::Integer, a::Real, p::Real, omega::Real; Nmax::Integer=80)
    s = -2; r0 = float(p); Δ = r0^2 - 2r0 + a^2
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
    ν = RenormalizedAngularMomentum(s, l, m, a, omega, λ)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=λ, nu=ν, Nmax=Nmax)
    ang = _angular_circular(s, l, m, a, omega; Nmax=Nmax)
    K = (r0^2 + a^2) * omega - a * m; dKdr = 2 * r0 * omega; dΔdr = 2 * r0 - 2
    function Dops(Rv, dRv, d2Rv)
        DR = dRv - im * K / Δ * Rv
        dDRdr = d2Rv - im * K / Δ * dRv - im * (dKdr * Δ - K * dΔdr) / Δ^2 * Rv
        DDR = dDRdr - im * K / Δ * DR
        return DR, DDR
    end
    RIn = R.In(r0); dRIn = teukolsky_radial_deriv(R.In, r0); d2RIn = _d2R_ode(RIn, dRIn, s, m, a, r0, omega, λ)
    RUp = R.Up(r0); dRUp = teukolsky_radial_deriv(R.Up, r0); d2RUp = _d2R_ode(RUp, dRUp, s, m, a, r0, omega, λ)
    DRIn, DDRIn = Dops(RIn, dRIn, d2RIn)
    DRUp, DDRUp = Dops(RUp, dRUp, d2RUp)
    Fin = _FabEq(r0, a, RIn, DRIn, DDRIn, ang.S, ang.L2S, ang.L1L2S)
    Fup = _FabEq(r0, a, RUp, DRUp, DDRUp, ang.S, ang.L2S, ang.L1L2S)
    W = (RIn * dRUp - dRIn * RUp) / Δ
    return (Fin = Fin, Fup = Fup, W = W)
end

"velocity projections (Kinnersley tetrad, circular equatorial)."
_vn(a, r0, En, Lz) = -((r0^2 + a^2) * En - a * Lz) / (r0^2 - 2r0 + a^2)
_vmb(a, En, Lz) = -im * (a * En - Lz)

"horizon absorption factor α; C² = generalized Teukolsky-Starobinsky constant."
function _alpha_circ(a, m, omega, λ)
    rplus = 1 + sqrt(1 - a^2); P = omega - m * a / (2 * rplus); ϵ = sqrt(1 - a^2) / (4 * rplus)
    C2 = ((λ + 2)^2 + 4a * omega * (m - a * omega)) * (λ^2 + 36a * omega * (m - a * omega)) -
         (2λ + 3) * (48a * omega * (m - 2a * omega)) + 144 * omega^2 * (1 - a^2)
    return 256 * (2 * rplus)^5 * P * (P^2 + 4ϵ^2) * (P^2 + 16ϵ^2) * omega^3 / C2
end

"""
    _circular_amplitudes_fabeq(l, m, a, p, En, Lz, Upst, omega, lambda; Nmax=80) -> (CPlus0, CMinus0, W)

Geodesic asymptotic amplitudes C^{+,0}/C^{-,0} (= Z^inf/Z^H) via the FabEq path; validated == the
point-particle source Zinf/ZH.
"""
function _circular_amplitudes_fabeq(l::Integer, m::Integer, a::Real, p::Real,
        En::Real, Lz::Real, Upst::Real, omega::Real, lambda::Real; Nmax::Integer=80)
    r0 = float(p); vn = _vn(a, r0, En, Lz); vmb = _vmb(a, En, Lz)
    f = _fabW(l, m, a, p, omega; Nmax=Nmax)
    sumP0 = vn^2 * f.Fin.Fnn + vn * vmb * f.Fin.Fnmb + vmb^2 * f.Fin.Fmbmb
    sumM0 = vn^2 * f.Fup.Fnn + vn * vmb * f.Fup.Fnmb + vmb^2 * f.Fup.Fmbmb
    return (CPlus0 = -2 * pi * sumP0 / (Upst * f.W), CMinus0 = -2 * pi * sumM0 / (Upst * f.W), W = f.W)
end

"""
    circular_spin_flux(l, m, a, p; Nmax=80, h=3e-4) -> NamedTuple

Linear-in-spin gravitational-wave flux correction (dE¹, dLz¹ at infinity ℐ / horizon ℋ) for a
circular-equatorial s=-2 spinning secondary. The ∂_ω F_ab, ∂_ω W, ∂_ω α are central finite
differences (velocities fixed); G_a is the spin-dipole source (FabEq Gn,Gmb). Returns
`(; dE1I, dE1H, dL1I, dL1H, CPlus0, CMinus0, omega, lambda)`.
"""
function circular_spin_flux(l::Integer, m::Integer, a::Real, p::Real; Nmax::Integer=80, h::Real=3e-4)
    s = -2; r0 = float(p)
    En = geodesic_energy(a, p, 0.0, 1.0); Lz = geodesic_angular_momentum(a, p, 0.0, 1.0)
    Q = carter_constant(a, p, 0.0, 1.0); K0 = Q + (Lz - a * En)^2; sK = sqrt(K0)
    Ω = kerr_frequencies(a, p, 0.0, 1.0).Omegaph
    Upst = kerr_frequencies(a, p, 0.0, 1.0; Time="Mino").Ut
    ω = m * Ω
    ω1 = m * spin_linear_frequencies(a, p, 0.0, 1.0)[3]    # GW frequency shift = m δΩφ
    Υt1 = spin_mino_time_correction(a, p, 0.0, 1.0)        # -3Υτ/(2√K0)
    λ = SpinWeightedSpheroidalEigenvalue(s, l, m, a * ω)
    vn = _vn(a, r0, En, Lz); vmb = _vmb(a, En, Lz)
    contract(F) = vn^2 * F.Fnn + vn * vmb * F.Fnmb + vmb^2 * F.Fmbmb
    f0 = _fabW(l, m, a, p, ω; Nmax=Nmax)
    fp = _fabW(l, m, a, p, ω + h; Nmax=Nmax)
    fm = _fabW(l, m, a, p, ω - h; Nmax=Nmax)
    sumP0 = contract(f0.Fin);  sumM0 = contract(f0.Fup)
    dsumPdω = (contract(fp.Fin) - contract(fm.Fin)) / (2h)        # ∂_ω ΣvᵃvᵇF (fixed v)
    dsumMdω = (contract(fp.Fup) - contract(fm.Fup)) / (2h)
    GP = (vn * f0.Fin.Gn + vmb * f0.Fin.Gmb) / sK                 # spin-dipole source contraction
    GM = (vn * f0.Fup.Gn + vmb * f0.Fup.Gmb) / sK
    sumP1 = GP + dsumPdω * ω1;  sumM1 = GM + dsumMdω * ω1
    W = f0.W; dWdω = (fp.W - fm.W) / (2h)
    CPlus0  = -2 * pi * sumP0 / (Upst * W);  CMinus0 = -2 * pi * sumM0 / (Upst * W)
    pref = En / (2 * sK) - Υt1 / Upst - dWdω * ω1 / W
    CPlus1  = -2 * pi * sumP1 / (Upst * W) + pref * CPlus0
    CMinus1 = -2 * pi * sumM1 / (Upst * W) + pref * CMinus0
    α  = _alpha_circ(a, m, ω, λ)
    λp = SpinWeightedSpheroidalEigenvalue(s, l, m, a * (ω + h)); λm = SpinWeightedSpheroidalEigenvalue(s, l, m, a * (ω - h))
    dαdω = (_alpha_circ(a, m, ω + h, λp) - _alpha_circ(a, m, ω - h, λm)) / (2h)
    reP = 2 * real(CPlus1 * conj(CPlus0)); reM = 2 * real(CMinus1 * conj(CMinus0))
    dE1I = (reP - 2 * abs2(CPlus0) * ω1 / ω) / (4 * pi * ω^2)
    dE1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 2 * abs2(CMinus0) * ω1 / ω) / (4 * pi * ω^2)
    dL1I = (reP - 3 * abs2(CPlus0) * ω1 / ω) * m / (4 * pi * ω^3)
    dL1H = α * (dαdω * ω1 / α * abs2(CMinus0) + reM - 3 * abs2(CMinus0) * ω1 / ω) * m / (4 * pi * ω^3)
    return (dE1I = dE1I, dE1H = dE1H, dL1I = dL1I, dL1H = dL1H,
            CPlus0 = CPlus0, CMinus0 = CMinus0, omega = ω, lambda = λ)
end

"""
    TeukolskySpinModeCircularCorrectionAnalytical(l, m, a, p; Nmax=80, h=3e-4) -> NamedTuple

Circular-equatorial s=-2 spinning-secondary mode (BHPToolkit name). Returns the geodesic amplitudes/
fluxes + the linear-in-spin flux correction dE¹/dLz¹: `(; l,m,a,p,omega, Amplitudes=(I,H),
Fluxes=Energy/AngularMomentum→(I,H), FluxesCorrection=Energy/AngularMomentum→(I,H))`.
"""
function TeukolskySpinModeCircularCorrectionAnalytical(l::Integer, m::Integer, a::Real, p::Real;
        Nmax::Integer=80, h::Real=3e-4)
    sf = circular_spin_flux(l, m, a, p; Nmax=Nmax, h=h)
    eg = energy_flux(-2, l, m, a, sf.omega, sf.lambda, sf.CPlus0, sf.CMinus0)
    jg = angular_momentum_flux(-2, l, m, a, sf.omega, sf.lambda, sf.CPlus0, sf.CMinus0)
    return (l = l, m = m, a = float(a), p = float(p), omega = sf.omega,
            Amplitudes = (I = sf.CPlus0, H = sf.CMinus0),
            Fluxes = (Energy = (I = eg.I, H = eg.H), AngularMomentum = (I = jg.I, H = jg.H)),
            FluxesCorrection = (Energy = (I = sf.dE1I, H = sf.dE1H),
                                AngularMomentum = (I = sf.dL1I, H = sf.dL1H)))
end

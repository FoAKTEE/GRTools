# Relativistic Lindblad/corotation resonance torque (Hirata 1010.0758/0759 GW-overlap link) + the
# reduced-ODE forcing source 𝔉_m — the leaf that CLOSES the Modes submodule (arXiv:2601.19123).
# The Lindblad torque is computed by REUSING the Perturbation engine (MST amplitudes / point-particle
# Z's) — the same Teukolsky overlap the legacy `DiscModes.lindblad_torque` used, now routed through
# the GrTools `Perturbation` runtime dep and `Harmonics` SWSH eigenvalue. GENERAL Kerr (all spins);
# the equatorial geodesic/wave factors reduce to Schwarzschild at a=0. The resonance radius r1 is the
# Lindblad root already solved by `Resonances.lindblad_locations` (AHA — the same condition
# m(Ω−Ω_p)=±κ drives the Wronskian, the WKB turning points, and the torque). From
# DiscModes/lindblad_torque.jl. Type-generic in the disc background; the MST overlap is Float64.

using ..Modes: RelativisticDisc, disc_OmegaK, kappar2, delta_r, csound, risco_disc
using ..Modes: lindblad_residual, lindblad_locations
import Perturbation: RenormalizedAngularMomentum, TeukolskyPointParticleModeCircular,
                     TeukolskyPointParticleModeEccentric
import Perturbation: _mst_amplitudes
import Harmonics: SpinWeightedSpheroidalEigenvalue

# Kerr equatorial geodesic/wave factors (M=1; prograde a>0, retrograde a<0, ẋ>0). Reduce to
# Schwarzschild at a=0; identities: _kOmega == disc_OmegaK, _kkappa² == kappar2, _kDelta == delta_r.
_kOmega(a, r) = 1 / (r^1.5 + a)
_kkappa(a, r) = _kOmega(a, r) * sqrt(1 - 6/r + 8a*r^-1.5 - 3a^2/r^2)
_kwt(a, r)    = (1 + a*r^-1.5) / sqrt(1 - 3/r + 2a*r^-1.5)                 # u^t (≡ metric)
_kDelta(a, r) = r^2 - 2r + a^2
_kZ(a, r)     = _kkappa(a, r) * _kwt(a, r) * (r^2 / _kDelta(a, r))         # impedance Z = κ u^t r²/Δ
_kgD(m, a, r, kind) = m*_kOmega(a, r) + (kind === :ILR ? -1 : 1)*_kkappa(a, r)   # r0-indep part of D
_kDprime(m, a, r, kind; h=1e-6) = (_kgD(m, a, r+h, kind) - _kgD(m, a, r-h, kind)) / (2h)

"""
    _c13_alpha(l, m, a, w; Nm=300) -> (c13, alpha)

Kerr `c₁₃ = B_ref/B_trans` (the In reflection/transmission ratio) and the horizon-flux factor `α`
for `(ℓ,m)` at spin `a`, frequency `w`, via the `Perturbation` MST amplitudes + the `Harmonics`
spheroidal eigenvalue / `RenormalizedAngularMomentum` ν-solve (Hirata ERRATA #1/#8 prefactor).
"""
function _c13_alpha(l, m, a, w; Nm=300)
    lam_c = SpinWeightedSpheroidalEigenvalue(-2, l, m, a*w)
    nu    = RenormalizedAngularMomentum(-2, l, m, a, w, lam_c)
    amp   = _mst_amplitudes(-2, m, a, w, nu, lam_c; Nmax=Nm)
    c13   = amp.InRef / amp.InTrans; lam = real(lam_c)
    rp = 1 + sqrt(1 - a^2); Omh = a/(2*rp); pw = w - m*Omh             # ϖ = ω − mΩ_H
    Gam = 4*sqrt(1 - a^2)/(rp^2 + a^2); chi = a*w
    Csq = ((lam+2)^2 + 4*m*chi - 4*chi^2)*(lam^2 + 36*m*chi - 36*chi^2) +
          48*(2*lam+3)*chi*(2*chi - m) + 144*w^2*(1 - a^2)
    pref = 8192*rp^5*pw*w^3/Csq
    (c13 = c13, alpha = pref*(pw^2 + (Gam/2)^2)*(pw^2 + Gam^2))
end

"""
    _test_Zout(l, m, n, a, r1; eR=1e-2, Nm=300) -> (Zout, ZH)

Eccentric `n=±1` test amplitudes `(Z^out, Z^H)` at the resonance radius `r1`, Richardson-extrapolated
in `e²` (kills the a<0 small-e n=0 leakage; an identity for a≥0). Built on
`Perturbation.TeukolskyPointParticleModeEccentric`.
"""
function _test_Zout(l, m, n, a, r1; eR=1e-2, Nm=300)
    t1 = TeukolskyPointParticleModeEccentric(l, m, n, a, r1, eR;  s=-2, Nmax=Nm, N=64)
    t2 = TeukolskyPointParticleModeEccentric(l, m, n, a, r1, 2eR; s=-2, Nmax=Nm, N=64)
    (4*(t1.Zinf/eR) - (t2.Zinf/(2eR)))/3, t1.ZH/eR
end

"""
    lindblad_torque(d::RelativisticDisc, r_p, m; kind=:ILR, lmax=24, Nm=300) -> N

Relativistic Lindblad resonant torque `N` (ILR/OLR, mode `m`, perturber at `r_p`) via the Hirata
GW-overlap, REUSING the `Perturbation` Teukolsky engine. General Kerr (spin `d.a`). The resonance
radius is the Lindblad root from [`lindblad_locations`](@ref) (the same m(Ω−Ω_p)=±κ condition).
Verified vs Hirata tab:ilrS (a=0) + fig:kerrm1 m=1 ILR (a=±0.9).
"""
function lindblad_torque(d::RelativisticDisc, r_p, m; kind::Symbol=:ILR, lmax=24, Nm=300)
    a = Float64(d.a); rp = Float64(r_p); w = m*_kOmega(a, rp)
    rILR, rOLR = lindblad_locations(d, r_p, m)             # AHA: reuse the Resonances root-finder
    r1 = Float64(kind === :ILR ? rILR : rOLR)
    isnan(r1) && return NaN
    n = kind === :ILR ? -1 : 1; s2pi = sqrt(2π); Pm = 0.0im; lastgood = 0.0im
    for l in max(m, 2):lmax
        ca = _c13_alpha(l, m, a, w; Nm=Nm)
        pert = TeukolskyPointParticleModeCircular(-2, l, m, a, rp; Nmax=Nm)
        Zout_t, ZH_t = _test_Zout(l, m, n, a, r1; Nm=Nm)
        if kind === :ILR
            term = (pert.ZH/s2pi)*(ca.c13*conj(Zout_t/s2pi) + ca.alpha*conj(ZH_t/s2pi))     # eq:P-ILR
        else
            term = (pert.Zinf/s2pi)*(conj(Zout_t/s2pi) - conj(ca.c13)*conj(ZH_t/s2pi))       # eq:P-OLR
        end
        (isnan(term) || isinf(term)) && break
        Pm += term; abs(term) < 1e-15*abs(Pm) && break; lastgood = Pm
    end
    Pm = -lastgood/w^2
    Zimp = _kZ(a, r1); Sm = 2im/(m*_kOmega(a, rp)*Zimp*r1)*Pm
    sgn = kind === :ILR ? -1 : 1
    sgn * π*m*_kwt(a, r1)*Zimp/abs(_kDprime(m, a, r1, kind))*abs2(Sm)
end

# --- 𝔉_m reduced-ODE forcing source --------------------------------------------------------------
# The m-th azimuthal Fourier component of the perturber tidal potential, in the reduced-ODE
# normalization 𝔥''+A_m𝔥'+B_m𝔥=𝔉_m (eq:reduceeq). The DEFAULT potential component is the Newtonian
# razor-thin softened point-mass m-potential, whose m-mode amplitude is the Laplace coefficient
# ½ b_{1/2}^{(m)}(β) (Goldreich-Tremaine 1979). The resonant amplification lives in the homogeneous
# operator B_m (which → 0 at the Lindblad turning points lindblad_residual=0); the source itself is
# regular there. The RIGOROUS height-averaged source operator (the unpublished master-PDE reduction)
# is carried as [HOLE: H-DISC-FULLOP] — this is the honest WKB-consistent construction, validated by
# the turning-point structure + the corotation singular exponent already grounded in ReducedODE.

"""
    laplace_b(s, m, β; N=512) -> b_s^{(m)}(β)

Laplace coefficient `b_s^{(m)}(β) = (2/π) ∫₀^π cos(mφ) (1 − 2β cosφ + β²)^{-s} dφ` by the composite
trapezoidal rule (the integrand is smooth + periodic ⇒ spectral convergence). `0 ≤ β < 1`. The
analytic kernel of the Newtonian m-mode perturber potential.
"""
function laplace_b(s, m::Integer, β::T; N::Integer=512) where {T<:Real}
    acc = zero(T); h = T(π)/N
    for k in 0:N
        φ = k*h
        w = (k == 0 || k == N) ? one(T)/2 : one(T)
        acc += w * cos(m*φ) * (1 - 2β*cos(φ) + β^2)^(-s)
    end
    2*h*acc/T(π)
end

"""
    forcing_Fm(d::RelativisticDisc, r, r_p, m; Φm=nothing) -> 𝔉_m(r)

Reduced-ODE forcing source 𝔉_m(r) (eq:reduceeq RHS): the perturber m-mode tidal potential divided by
the local sound speed squared (the thin-disc density-wave normalization). `Φm` is the m-mode
potential-amplitude callable `Φm(r)`; the default is the Newtonian softened point-mass component
`−½ b_{1/2}^{(m)}(β)` with `β = min(r,r_p)/max(r,r_p)`. The source is regular at the Lindblad
resonances (the singular response sits in B_m); near corotation it carries the (r−r_p) structure that
fixes the validated corotation exponent. The RIGOROUS height-averaged source is [HOLE: H-DISC-FULLOP].
"""
function forcing_Fm(d::RelativisticDisc{T}, r, r_p, m; Φm=nothing) where {T}
    rr = T(r); rp = T(r_p)
    pot = Φm === nothing ? -laplace_b(one(T)/2, m, min(rr, rp)/max(rr, rp))/2 : T(Φm(rr))
    pot / csound(d, rr)^2
end

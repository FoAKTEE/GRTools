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
    c13_alpha(l, m, a, w; Nm=300) -> (; c13, c14, aleph, alpha, lambda, Csq)

Hirata (1010.0759) amplitude-convention bundle of the `(ℓ,m)` Teukolsky In solution at spin `a`,
frequency `w`, built on the `Perturbation` MST amplitudes and the `Harmonics` spheroidal
eigenvalue / `RenormalizedAngularMomentum` ν-solve:

  * `c13 = B_ref/B_trans`, `c14 = B_inc/B_trans` — Hirata's R1 = c₁₃R₃ + c₁₄R₄ connection
    coefficients (eq:cmatrix), reconstructed from the unit-transmission In amplitudes;
  * `aleph = 2iω c₁₄`;
  * `alpha` — the horizon-flux factor 8192 r₊⁵ ϖ ω³ (ϖ²+(Γ/2)²)(ϖ²+Γ²)/|C|² in the
    Teukolsky–Press convention Γ = 2√(1−a²)/(r₊²+a²) = 2κ_H, ε_TP = Γ/2 (ERRATA #1/#8: the
    printed (ϖ²+Γ²)(ϖ²+4Γ²) bracket is replaced by (ϖ²+(Γ/2)²)(ϖ²+Γ²)). NOTE: the legacy
    `DiscModes._c13_alpha` and the legacy arXiv:1010.0759 figure scripts carried Γ = 2κ_H·2,
    i.e. exactly the printed (erroneous) bracket; that fork FAILS the closure identity below by
    up to 9.1e-5 and is corrected here (it agrees with `Perturbation`'s independent horizon
    factor ε = √(1−a²)/(4r₊), a fold candidate);
  * `lambda` — the (real) spin-weighted spheroidal eigenvalue, `Csq` the Teukolsky–Starobinsky |C|².

These satisfy the closure identity eq:cplex1  `|c₁₃|² + α = (2ω)⁶ |ℵ|² / |C|²`  to machine precision
(the gate that pins the three convention factors; see `reproduce/arXiv-1010.0759/check_c13_alpha.jl`).
"""
function c13_alpha(l, m, a, w; Nm=300)
    lam_c = SpinWeightedSpheroidalEigenvalue(-2, l, m, a*w)
    nu    = RenormalizedAngularMomentum(-2, l, m, a, w, lam_c)
    amp   = _mst_amplitudes(-2, m, a, w, nu, lam_c; Nmax=Nm)
    c13   = amp.InRef / amp.InTrans; c14 = amp.InInc / amp.InTrans; lam = real(lam_c)
    rp = 1 + sqrt(1 - a^2); Omh = a/(2*rp); pw = w - m*Omh             # ϖ = ω − mΩ_H
    Gam = 2*sqrt(1 - a^2)/(rp^2 + a^2); chi = a*w      # Γ = 2κ_H (r₊²+a² ≡ 2r₊); ε_TP = Γ/2
    Csq = ((lam+2)^2 + 4*m*chi - 4*chi^2)*(lam^2 + 36*m*chi - 36*chi^2) +
          48*(2*lam+3)*chi*(2*chi - m) + 144*w^2*(1 - a^2)
    pref = 8192*rp^5*pw*w^3/Csq
    (c13 = c13, c14 = c14, aleph = 2im*w*c14, alpha = pref*(pw^2 + (Gam/2)^2)*(pw^2 + Gam^2),
     lambda = lam, Csq = Csq)
end

"Internal alias kept for the characterization tests; see [`c13_alpha`](@ref)."
const _c13_alpha = c13_alpha

"""
    _test_Zout(l, m, n, a, r1; eR=1e-2, Nm=300) -> (Zout, ZH)

Eccentric `n=±1` test amplitudes `(Z^out, Z^H)` at the resonance radius `r1`, per unit `e`:
`Z^out` is the linear-in-`e` amplitude `Z^∞/e` Richardson-extrapolated in `e²` from `e = eR, 2eR`;
`Z^H` is `Z^H/e` at `e = eR`. Built on `Perturbation.TeukolskyPointParticleModeEccentric`.
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

"""
    laplace_b_deriv(s, m, β; N=512) -> b_s^{(m)′}(β)

β-derivative of the Laplace coefficient [`laplace_b`](@ref), differentiated under the integral sign:

    b_s^{(m)′}(β) = (4s/π) ∫₀^π cos(mφ) (cos φ − β) (1 − 2β cos φ + β²)^{-s-1} dφ,

evaluated with the same composite trapezoidal rule (smooth periodic integrand ⇒ spectral
convergence). `0 ≤ β < 1`. Together with `laplace_b` it is the analytic kernel of the Hirata
Newtonian resonance amplitude S^{(m)} ([`lindblad_strength`](@ref)); `b_{1/2}^{(1)′}(0) = 1` is the
indirect-term coefficient of the m=1 outer resonance (1010.0759 eq:smolr).
"""
function laplace_b_deriv(s, m::Integer, β::T; N::Integer=512) where {T<:Real}
    acc = zero(T); h = T(π)/N
    for k in 0:N
        φ = k*h
        w = (k == 0 || k == N) ? one(T)/2 : one(T)
        acc += w * cos(m*φ) * (cos(φ) - β) * (1 - 2β*cos(φ) + β^2)^(-s-1)
    end
    4*s*h*acc/T(π)
end

"""
    epicyclic_impedance(d::RelativisticDisc, r) -> Z

Specific epicyclic impedance of a circular equatorial Kerr geodesic at `r`,

    Z(r) = κ_r(r) u^t(r) r² / Δ_r(r)                                    (Hirata 1010.0759 l.954),

the factor that converts the resonant overlap S^{(m)} into the normalized torque
N = ∓π m u^t Z |S^{(m)}|²/|D′| ([`lindblad_strength`](@ref), [`lindblad_torque`](@ref)). At a=0 it
reduces to (r−2)^{-1} √((r−6)/(r(r−3))). Identical (verified to 5e-9, the legacy value being
finite-difference limited) to the metric form √(C₀₂₀/C₀₀₂) used by the legacy arXiv:2510.02433
torque library — that fork is COLLAPSED onto this closed form.
"""
epicyclic_impedance(d::RelativisticDisc, r) = _kZ(d.a, r)

"""
    lindblad_detuning(d, r, r_p, m; kind=:ILR) -> D(r)

Signed Lindblad detuning `D(r) = m(Ω^k_r − Ω^k_{r_p}) ∓ κ_r` (upper sign `:ILR`, lower `:OLR`).
Its roots are the resonance radii returned by [`lindblad_locations`](@ref); its square is the
[`lindblad_residual`](@ref) (whose zeros drive the Wronskian).
"""
lindblad_detuning(d::RelativisticDisc, r, r_p, m; kind::Symbol=:ILR) =
    _kgD(m, d.a, r, kind) - m*disc_OmegaK(d, r_p)

"""
    lindblad_detuning_deriv(d, r, m; kind=:ILR, h=nothing) -> D′(r)

Radial derivative of [`lindblad_detuning`](@ref) (independent of `r_p`), central finite difference
with the SCALE-AWARE step `h = 1e-6·max(1, r)` (a fixed absolute step loses all significance at
r ≫ 10⁴, where the legacy fixed-h evaluation carries ~2e-4 relative error).
`|D′(r₁)|` is the resonance width factor of the normalized torque N.
"""
lindblad_detuning_deriv(d::RelativisticDisc, r, m; kind::Symbol=:ILR, h=nothing) =
    _kDprime(m, d.a, r, kind; h = h === nothing ? 1e-6*max(one(float(r)), r) : h)

# --- Hirata Newtonian resonance amplitude S^{(m)} -------------------------------------------------
# The ℓ-sum of the GW overlap eq:P resums (RT-A2 generating-function chain) into a closed Laplace
# bracket in ς = r₁/r_p. ILR/OLR brackets (1010.0759 l.904 / eq:smolr, the latter with the ς⁺
# argument in b′ — the published typo #2 corrected):
#   ILR: ς b′(ς) + 2m b(ς)
#   OLR: −ς b′(ς) + 2m b(ς) − 4 ς^{-2} δ_{m1}     (b^{(1)′}(0)=1 indirect term)
# Quadrature resolution for the bracket: the Laplace integrand (1−2ςcosφ+ς²)^{-1/2} develops a
# |2 sin(φ/2)|^{-1} near-singularity as ς→1 (high m: the resonance crowds onto the orbit) and
# oscillates m times, so the trapezoidal node count must grow like max(8m, 64/|1−ς|) — at the default
# N=512 the m=160 coefficient is already 11% wrong (measured), which is a ~2% error on the total
# Lindblad torque at small r_s.
_bracket_nodes(m::Integer, ς::Real) = clamp(max(512, 8m, ceil(Int, 64/max(abs(1 - ς), 1e-4))), 512, 1 << 16)

"Laplace bracket of the Hirata resonance amplitude at ς (see above); `kind` ∈ (:ILR, :OLR)."
function _resonance_bracket(m::Integer, ς::T, kind::Symbol) where {T<:Real}
    N  = _bracket_nodes(m, ς)
    b  = laplace_b(one(T)/2, m, ς; N = N)
    bp = laplace_b_deriv(one(T)/2, m, ς; N = N)
    if kind === :ILR
        return ς*bp + 2m*b
    elseif kind === :OLR
        return -ς*bp + 2m*b - (m == 1 ? 4/ς^2 : zero(T))
    else
        throw(ArgumentError("kind must be :ILR or :OLR, got $(kind)"))
    end
end

"""
    lindblad_strength(d::RelativisticDisc, r_p, m; kind=:ILR) -> N

SIGNED normalized Lindblad torque strength `N = T/(2π r₁ q² Σ)` of mode `m` for a perturber on a
circular equatorial orbit at `r_p`, in the HYBRID scheme of arXiv:2510.02433: EXACT relativistic Kerr
kinematics (resonance radius r₁ from [`lindblad_locations`](@ref), u^t, `epicyclic_impedance`, D′)
with the Hirata Newtonian-limit overlap amplitude evaluated at the relativistic geometry ς = r₁/r_p,

    N = ∓ π m u^t(r₁) Z(r₁) |S^{(m)}|² / |D′(r₁)|,   |S^{(m)}|² = (ς/4r_p) [bracket]²

(upper sign / `:ILR` ⇒ N < 0, the secondary is torqued inward; `:OLR` ⇒ N > 0). Returns `NaN` when
the resonance is not bracketed in the disc (e.g. inside the ISCO). The FULL relativistic amplitude —
the Teukolsky ℓ-sum this closed form approximates — is [`lindblad_torque`](@ref); their ratio is the
"relativistic enhancement" of the arXiv:2510.02433 slope analysis.
"""
function lindblad_strength(d::RelativisticDisc, r_p, m::Integer; kind::Symbol=:ILR)
    rILR, rOLR = lindblad_locations(d, r_p, m)
    r1 = kind === :ILR ? rILR : rOLR
    isnan(r1) && return oftype(float(r1), NaN)
    ut = _kwt(d.a, r1); Z = _kZ(d.a, r1)
    (isnan(Z) || isnan(ut)) && return oftype(float(r1), NaN)
    Dp = abs(lindblad_detuning_deriv(d, r1, m; kind=kind))
    ς = r1 / r_p
    S2 = (ς / (4 * r_p)) * _resonance_bracket(m, ς, kind)^2
    N = π * m * ut * Z * S2 / Dp
    isfinite(N) || return oftype(float(r1), NaN)
    return (kind === :ILR ? -one(N) : one(N)) * N
end

"""
    lindblad_strength_newtonian(m; kind=:ILR, T=Float64) -> N

`r_p`-INDEPENDENT Newtonian limit of [`lindblad_strength`](@ref) (a=0, r_p→∞): u^t→1, Z→r₁^{-3/2},
|D′|→(3/2) r₁^{-5/2}(m∓1) and ς → ς_m^∓ = ((m∓1)/m)^{2/3}, so every r_p power cancels:

    N = ∓ (π/6) · m/(m∓1) · ς² · [bracket]²                        (1010.0759 RT-N1)

Published values (the paper's l.1008 misprints the m/(m∓1) factor; eq:smolr misprints the b′
argument — both corrected here): −2.36 (m=2 ILR), −7.50 (m=3 ILR), +0.19 (m=1 OLR), +8.62 (m=2 OLR).
"""
function lindblad_strength_newtonian(m::Integer; kind::Symbol=:ILR, T::Type{<:Real}=Float64)
    kind in (:ILR, :OLR) || throw(ArgumentError("kind must be :ILR or :OLR, got $(kind)"))
    δ = kind === :ILR ? -1 : 1
    m + δ > 0 || throw(ArgumentError("kind=:ILR needs m ≥ 2 (m=1 has no inner Lindblad resonance)"))
    ς = T((m + δ)//m)^(T(2)/3)
    br = _resonance_bracket(m, ς, kind)
    return T(δ) * (T(π)/6) * (m/T(m + δ)) * ς^2 * br^2
end

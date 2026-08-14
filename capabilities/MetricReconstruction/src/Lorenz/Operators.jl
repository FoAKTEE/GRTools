# Lorenz-gauge operator scalars + 1D radial jet algebra + promoted CG angular mixing.
#
# LEAF of the Lorenz submodule (design §4.4). Three concerns, consolidated from KerrLorenzCirc
# operators.jl + angular.jl into one type-generic foundation:
#   (1) GEOMETRY/TETRAD/CHANDRASEKHAR scalars — geom (Δ,Σ,ρ,ρc,Δp bundle), the l±/m± tetrad legs,
#       the tetrad/equatorial inverse metric, and the Kfn/Qfn/Dcoef operator coefficients.
#   (2) JET ALGEBRA — a "jet" is a length-(N+1) `Vector{Complex{T}}` of derivative values
#       `[f,f',…,f^(N)]` at r0; jmul/jadd/jscale/jderiv/jinv + constjet, all preserving the element
#       type so a Complex{BigFloat} jet keeps full precision. `jet_ops` builds the coefficient jets
#       and the Chandrasekhar operators/ODE-closure factory bound to (r0,a,ω,m) — the single source
#       of truth the closure inherits (AHA: the legacy closure_a0/closure_kerr/radial_a0 each carried
#       a verbatim copy of this).
#   (3) ANGULAR CG MIXING — the spheroidal→spherical coupling coefficients 𝔰ⁿ/𝔠ⁿ. KLC's hand-rolled
#       Racah `cg` is PROMOTED to `Wigner.clebsch_gordan` (a thin adapter), with the a=0 ladder
#       diagonals λ̂/λ̂₂/Λ̂.
#
# Kept SEPARATE from Radiation's 2D RTJet (AHA — different decision: this is a 1D radial Taylor tower).

import Wigner
import Spacetime: kfun, qfun, delta_kerr, sigma_kerr, zeta_kerr, kerr_metric_inv
# Canonical 1D derivative-tower jet algebra now lives in NumericsCore.jets (Move A). These names are
# brought into the flat Lorenz namespace so downstream Lorenz files (Closure.jl) and ModifiedGravity's
# qualified `MRL.jmul`/… keep resolving IDENTICALLY (same function objects); Lorenz re-exports
# jinv/constjet (Lorenz.jl:27) unchanged.
using NumericsCore: constjet, jadd, jscale, jmul, jderiv, jinv

# ── (3a) CG adapter — KLC's hand-rolled Racah `cg` promoted to Wigner.clebsch_gordan ─────────────
"""
    cg(j1,m1,j2,m2,j3,m3) -> Float64

Clebsch-Gordan coefficient ⟨j1 m1; j2 m2 | j3 m3⟩. Thin adapter promoting KLC's hand-rolled Racah
implementation to the already-grounded `Wigner.clebsch_gordan` (BigInt-exact internally).
"""
cg(j1, m1, j2, m2, j3, m3) = Wigner.clebsch_gordan(j1, m1, j2, m2, j3, m3)

# ── (3b) spheroidal→spherical CG mixing coefficients 𝔰ⁿ/𝔠ⁿ (eq:mixing-coeffs) ────────────────────
"𝔰^{s1,s2}_{ℓ,ℓ'} = ⟨s2 ℓ' m|sinθ|s1 ℓ m⟩ (single-sin, Δl=±1; Δs=2 sin² coupling for s1=∓1,s2=±1)."
function frak_s(s1, s2, l, lp, m)
    if s1==1 && s2==0;   return -sqrt(2*(2l+1)/(2lp+1))*cg(l,m,1,0,lp,m)*cg(l,-1,1,1,lp,0); end
    if s1==-1 && s2==0;  return  sqrt(2*(2l+1)/(2lp+1))*cg(l,m,1,0,lp,m)*cg(l,1,1,-1,lp,0); end
    if s1==2 && s2==1;   return -sqrt(2*(2l+1)/(2lp+1))*cg(l,m,1,0,lp,m)*cg(l,-2,1,1,lp,-1); end
    if s1==-2 && s2==-1; return  sqrt(2*(2l+1)/(2lp+1))*cg(l,m,1,0,lp,m)*cg(l,2,1,-1,lp,1); end
    if s1==-1 && s2==1;  return  sqrt(8*(2l+1)/(3*(2lp+1)))*cg(l,m,2,0,lp,m)*cg(l,1,2,-2,lp,-1); end
    if s1==1 && s2==2;   return frak_s(2,1,lp,l,m); end    # symmetry 𝔰^{s2 s1}_{ℓ' ℓ}=𝔰^{s1 s2}_{ℓ ℓ'}
    if s1==1 && s2==-1;  return frak_s(-1,1,lp,l,m); end
    if s1==-1 && s2==-2; return frak_s(-2,-1,lp,l,m); end
    if s1==0 && s2==1;   return frak_s(1,0,lp,l,m); end
    if s1==0 && s2==-1;  return frak_s(-1,0,lp,l,m); end
    error("frak_s($s1,$s2) not implemented")
end
"(𝔰²)^{s1,s2}_{ℓ,ℓ'} = ⟨s2 ℓ' m|sin²θ|s1 ℓ m⟩ (Δl up to ±2)."
function frak_s2(s1, s2, l, lp, m)
    if s1==2 && s2==0;  return sqrt(8*(2l+1)/(3*(2lp+1)))*cg(l,m,2,0,lp,m)*cg(l,-2,2,2,lp,0); end
    if s1==-2 && s2==0; return sqrt(8*(2l+1)/(3*(2lp+1)))*cg(l,m,2,0,lp,m)*cg(l,2,2,-2,lp,0); end
    if s1==0 && s2==2;  return frak_s2(2,0,lp,l,m); end
    if s1==0 && s2==-2; return frak_s2(-2,0,lp,l,m); end
    error("frak_s2($s1,$s2) not implemented")
end
"𝔠^{s,s}_{ℓ,ℓ'} = ⟨s ℓ' m|cosθ|s ℓ m⟩ (single-cos, Δl=±1)."
frak_c(s, l, lp, m) = sqrt((2l+1)/(2lp+1))*cg(l,m,1,0,lp,m)*cg(l,-s,1,0,lp,-s)
"(𝔠²)^{s,s}_{ℓ,ℓ'} = ⟨s ℓ' m|cos²θ|s ℓ m⟩ (Δl=0,±2)."
frak_c2(s, l, lp, m) = (l==lp ? 1/3 : 0.0) + 2/3*sqrt((2l+1)/(2lp+1))*cg(l,m,2,0,lp,m)*cg(l,-s,2,0,lp,-s)

# a=0 angular-ladder diagonal factors (consumed by Closure's a=0 path).
"λ̂(l) = √(l(l+1)) — the s=1 ladder eigenvalue (a=0)."
λ̂(l) = sqrt(l*(l+1))
"λ̂₂(l) = √((l-1)(l+2)) — the spin-2→1 ladder eigenvalue (a=0)."
λ̂2(l) = sqrt((l-1)*(l+2))
"Λ̂(l) = √((l-1)l(l+1)(l+2)) = λ̂·λ̂₂ — the composite spin-2 ladder eigenvalue (a=0)."
Λ̂(l) = λ̂(l)*λ̂2(l)

# ── (1) geometry / tetrad / Chandrasekhar scalars (type-generic where used) ───────────────────────
"Kerr geometry scalar bundle at (a,r,θ): (Δ,Σ,ζ̄ from Spacetime) + ρ=r+ia cosθ and Δ'=2(r−1) extras."
geom(a, r, θ) = (Δ=delta_kerr(a,r), Σ=sigma_kerr(a,r,θ), ρ=complex(r,a*cos(θ)), ρc=zeta_kerr(a,r,θ), Δp=2(r-1))

"Tetrad l± leg (eq:tetrad2), as a Complex{T} 4-vector; sgn=+1 ⇒ l, sgn=−1 ⇒ n (up to normalization)."
lpm(a, r, sgn) = Complex{typeof(float(real(a*r)))}[sgn*(r^2+a^2)/(r^2-2r+a^2), 1, 0, sgn*a/(r^2-2r+a^2)]
"Tetrad m± leg as a Complex{T} 4-vector; sgn=+1 ⇒ m, sgn=−1 ⇒ m̄."
mpm(a, θ, sgn) = (T=typeof(float(real(a*θ))); Complex{T}[sgn*im*a*sin(θ), 0, 1, sgn*im/sin(θ)])

_sym(X, Y) = 0.5*(X*transpose(Y) .+ Y*transpose(X))
"Inverse Kerr metric g^{μν} reconstructed from the l±/m± tetrad (eq:inverse-metric), type-generic."
ginv_tetrad(a, r, θ) = (g=geom(a,r,θ); (g.Δ*_sym(lpm(a,r,+1),lpm(a,r,-1)) .+ _sym(mpm(a,θ,+1),mpm(a,θ,-1)))./g.Σ)

"""
    kerr_metric_inv_eq(a, r) -> 4×4 g^{μν}

Inverse Kerr metric at the equator (θ=π/2), type-generic. Thin equatorial specialization of
`Spacetime.kerr_metric_inv` — the (t,φ) block + r,θ diagonal for u^t = −E g^{tt} + L g^{tφ}.
"""
kerr_metric_inv_eq(a::Real, r::Real) = (T=typeof(float(r*a)); kerr_metric_inv(T(a), T(r), T(π)/2))

"Chandrasekhar K(r) = ω(r²+a²) − a m (eq:DLoperators) — Spacetime.kfun."
const Kfn = kfun
"Chandrasekhar Q(θ) = m/sinθ − aω sinθ — Spacetime.qfun (arg order θ-last)."
Qfn(a, θ, ω, m) = qfun(a, ω, m, θ)
"Radial operator coefficient c = iK/Δ (the c in 𝒟=∂_r−c, 𝒟†=∂_r+c)."
Dcoef(a, r, ω, m) = im*Kfn(a,r,ω,m)/(r^2-2r+a^2)

# ── (2) jet algebra — CANONICAL HOME = NumericsCore.jets (Move A) ─────────────────────────────────
# `jmul/jadd/jscale/jderiv/jinv/constjet` are `using`-imported from NumericsCore (top of file), so
# the generic-length Leibniz derivative-tower primitives are the SINGLE canonical implementation
# (characterized bit-exact vs this former copy in core/NumericsCore/test/test_jets.jl). Only `binom`
# stays local: it is the Leibniz weight consumed directly by the physics ODE closures below
# (close_jet/close_jet_src, jet_ops) — those are Kerr/Teukolsky-flavored and remain in Lorenz.
"Binomial coefficient n choose k (exact, BigInt-safe via integer division)."
binom(n, k) = factorial(n) ÷ (factorial(k)*factorial(n-k))

"""
    jet_ops(r0, a, ω, m; N=6) -> NamedTuple

Build the jet (derivative-value Taylor) algebra about `r0` for spin `a`, frequency `ω`, azimuthal
`m`, to order `N`. Returns the coefficient jets (`Δjet,Kjet,rjet,invΔjet,cjet,Δp_over_Δ`), the element
type `CT=Complex{T}`, the Chandrasekhar operators (`opD,opDd,opDn`), the Teukolsky ODE α/β builders
(`jets_sm2,jets_sp2,jets_sm1,jets_s0`), and the closures (`close_jet,close_jet_src,cj`). Type-generic:
BigFloat `r0/a/ω` → Complex{BigFloat} throughout.
"""
function jet_ops(r0::Real, a::Real, ω::Real, m::Integer; N::Integer=6)
    T = typeof(float(promote(r0, a, ω)[1]))
    CT = Complex{T}
    r0T = T(r0); aT = T(a); ωT = T(ω)
    cj(s) = constjet(s, N, CT)

    Δjet = CT[r0T^2-2r0T+aT^2, 2r0T-2, 2, 0,0,0,0][1:N+1]
    Kjet = CT[ωT*(r0T^2+aT^2)-aT*m, 2ωT*r0T, 2ωT, 0,0,0,0][1:N+1]   # K = ω(r²+a²) − a m
    rjet = CT[r0T, 1, 0, 0,0,0,0][1:N+1]
    invΔjet = jinv(Δjet)
    cjet = jscale(im, jmul(Kjet, invΔjet))                          # c = iK/Δ
    Δp_over_Δ = jmul(jderiv(Δjet), invΔjet)

    opD(F)  = jadd(jderiv(F), jscale(-one(CT), jmul(cjet, F)))      # 𝒟  = ∂_r − c
    opDd(F) = jadd(jderiv(F), jmul(cjet, F))                        # 𝒟† = ∂_r + c
    opDn(n, F) = jadd(opD(F), jscale(CT(n), jmul(Δp_over_Δ, F)))    # 𝒟_n = 𝒟 + nΔ'/Δ

    function close_jet(J0, J1, αjet, βjet)
        F = zeros(CT, N+1); F[1] = CT(J0); F[2] = CT(J1)
        for k in 0:(N-2)
            F[k+3] = sum(binom(k,j)*(αjet[j+1]*F[(k+1-j)+1]+βjet[j+1]*F[(k-j)+1]) for j in 0:k)
        end
        F
    end
    function close_jet_src(J0, J1, αjet, βjet, γjet, Sjet)
        F = zeros(CT, N+1); F[1] = CT(J0); F[2] = CT(J1)
        for k in 0:(N-2)
            F[k+3] = sum(binom(k,j)*(αjet[j+1]*F[(k+1-j)+1]+βjet[j+1]*F[(k-j)+1]) for j in 0:k) +
                     sum(binom(k,j)*γjet[j+1]*Sjet[(k-j)+1] for j in 0:k)
        end
        F
    end

    # ODE α,β jets (spin-general; c carries a). Each Teukolsky op = a2 P''+a1 P'+a0 P ⇒ α=−a1/a2, β=−a0/a2.
    function jets_sm2(Λ)  # s=−2: Δ𝒟†_{-1}𝒟P−6iωrP−ΛP=0
        cp = jderiv(cjet); csq = jmul(cjet, cjet)
        a1 = jscale(-one(CT), jderiv(Δjet))
        a0 = jadd(jadd(jmul(Δjet, jadd(jadd(jscale(-one(CT),cp), jscale(-one(CT),csq)), jmul(Δp_over_Δ,cjet))),
                jscale(-6im*ωT, rjet)), cj(-Λ))
        (jmul(jscale(-one(CT),a1), invΔjet), jmul(jscale(-one(CT),a0), invΔjet))
    end
    function jets_sp2(Λ)  # s=+2: Δ𝒟_{-1}𝒟†P+6iωrP−ΛP=0
        cp = jderiv(cjet); csq = jmul(cjet, cjet)
        a1 = jscale(-one(CT), jderiv(Δjet))
        a0 = jadd(jadd(jmul(Δjet, jadd(jadd(cp, jscale(-one(CT),csq)), jscale(-one(CT),jmul(Δp_over_Δ,cjet)))),
                jscale(6im*ωT, rjet)), cj(-Λ))
        (jmul(jscale(-one(CT),a1), invΔjet), jmul(jscale(-one(CT),a0), invΔjet))
    end
    function jets_sm1(λ)  # s=−1: Δ𝒟†𝒟P−2iωrP−λP=0
        cp = jderiv(cjet); csq = jmul(cjet, cjet)
        a0 = jadd(jadd(jmul(Δjet, jadd(jscale(-one(CT),cp), jscale(-one(CT),csq))), jscale(-2im*ωT, rjet)), cj(-λ))
        (zeros(CT, N+1), jmul(jscale(-one(CT),a0), invΔjet))
    end
    function jets_s0(λ0)  # trace/s=0: 𝒟(Δ𝒟†h)−2iωrh−λ0h=0
        cp = jderiv(cjet); csq = jmul(cjet, cjet)
        a1 = jderiv(Δjet)
        a0 = jadd(jadd(jadd(jmul(jderiv(Δjet),cjet), jmul(Δjet,cp)), jscale(-one(CT),jmul(Δjet,csq))),
                jadd(jscale(-2im*ωT, rjet), cj(-λ0)))
        (jmul(jscale(-one(CT),a1), invΔjet), jmul(jscale(-one(CT),a0), invΔjet))
    end

    (CT=CT, T=T, N=N, Δjet=Δjet, Kjet=Kjet, rjet=rjet, invΔjet=invΔjet, cjet=cjet, Δp_over_Δ=Δp_over_Δ,
     cj=cj, opD=opD, opDd=opDd, opDn=opDn,
     close_jet=close_jet, close_jet_src=close_jet_src,
     jets_sm2=jets_sm2, jets_sp2=jets_sp2, jets_sm1=jets_sm1, jets_s0=jets_s0)
end

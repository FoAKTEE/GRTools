# =============================================================================================================
# TetradComponents.jl — the six Hertz/CCK tetrad-frame metric-perturbation components for one (s=±2) mode,
# in the IRG gauge (H_nn, H_nm, H_mm; s=−2) and the ORG gauge (H_ll, H_lm, H_mm; s=+2). Each is a chain of
# the GHP/Teukolsky ladder operators (Operators.jl) on the separated mode product R[s]⊗S[s] (Modes.jl), times
# scalar background prefactors (1/ζ̄², Δ, sinθ, …). εg ≡ 1. The "barred" variants apply the GHP overbar:
# ζ↔ζ̄, R→R̄(ω̄), S→S̄(ω̄), and the i-flipped ladder operators (Dnbar/Lnbar); for the real-ω physical modes
# used here ω̄=ω, so the overbar is plain complex conjugation of the mode jets plus the i-flip.
#
# Each public Hxxx returns the COMPLEX scalar value at (ω,l,m,r,θ) (the e^{imφ−iωt} coordinate phase is
# applied downstream). The companion Hxxx_jet keeps the full analytic (r,θ) RTJet (value + spatial
# derivatives) for the linearized-Einstein residual; its `.v` field equals the scalar Hxxx exactly.
# Background scalars Δ/Σ/ζ are delegated to Spacetime. Type-generic Complex{T} (BigFloat inherited).
# =============================================================================================================

using Spacetime: delta_kerr, sigma_kerr, zeta_kerr

# Local background-scalar aliases (M=1): ζ = r − i a cosθ, ζ̄ = conj(ζ), Δ, Σ.
@inline _ζfun(a, r, θ)    = zeta_kerr(a, r, θ)
@inline _ζbarfun(a, r, θ) = conj(zeta_kerr(a, r, θ))
@inline _Δfun(a, r)       = delta_kerr(a, r)
@inline _Σfun(a, r, θ)    = sigma_kerr(a, r, θ)

# The barred mode functions R̄[s], S̄[s] are the complex conjugates of R[s], S[s] (real-ω ⇒ ω̄=ω).
_conjjet(j::Vector) = conj.(j)

"Base separable field R[s]·S[s] (one term); `barred` ⇒ R̄·S̄ (conjugated jets)."
function _RS(s::Integer, l::Integer, m::Integer, a, ω, r, θ; bc::Symbol=:Up, λ=nothing, barred::Bool=false)
    rj = radial_jet(s, l, m, a, real(ω), r; bc=bc, λ=λ)
    aj = angular_jet(s, l, m, a, real(ω), θ)
    barred && (rj = _conjjet(rj); aj = _conjjet(aj))
    sepfield(rj, aj)
end

"Pure-S field (radial jet ≡ 1); `barred` ⇒ S̄."
function _Sonly(s::Integer, l::Integer, m::Integer, a, ω, r, θ; λ=nothing, barred::Bool=false)
    C = typeof(complex(float(a*real(ω)*r*θ)))
    aj = angular_jet(s, l, m, a, real(ω), θ)
    barred && (aj = _conjjet(aj))
    one_rj = zeros(C, NJ); one_rj[1] = one(C)
    sepfield(one_rj, aj)
end

"Pure-R field (angular jet ≡ 1); `barred` ⇒ R̄."
function _Ronly(s::Integer, l::Integer, m::Integer, a, ω, r, θ; bc::Symbol=:Up, λ=nothing, barred::Bool=false)
    C = typeof(complex(float(a*real(ω)*r*θ)))
    rj = radial_jet(s, l, m, a, real(ω), r; bc=bc, λ=λ)
    barred && (rj = _conjjet(rj))
    one_aj = zeros(C, NJ); one_aj[1] = one(C)
    sepfield(rj, one_aj)
end

"Scalar (barred-aware) single R-value / S-value used in some H prefactors."
_Rval(s,l,m,a,ω,r;bc,λ,barred) = (v = radial_jet(s,l,m,a,real(ω),r;bc=bc,λ=λ)[1]; barred ? conj(v) : v)
_Sval(s,l,m,a,ω,θ;barred)      = (v = angular_jet(s,l,m,a,real(ω),θ)[1];          barred ? conj(v) : v)

# ===== IRG components (s=−2) =================================================================================
"""
    HIRG_nn(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing)

`H_nn^IRG = −(εg/4)(1/ζ̄²)( Lnadj[1, Lnadj[2, S]] − (2 i a/ζ) sinθ · Lnadj[2, S] ) R`, s=−2.
"""
function HIRG_nn(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = -2
    ζ = _ζfun(a, r, θ); ζb = _ζbarfun(a, r, θ); sθ = sin(θ)
    Sf  = _Sonly(s, l, m, a, ω, r, θ; λ=λ, barred=barred)
    L2  = Ln(-1, 2, Sf;  a=a, θ=θ, ω=ω, m=m)
    L12 = Ln(-1, 1, L2;  a=a, θ=θ, ω=ω, m=m)
    isign = barred ? -1 : 1
    inner = L12 + (-(isign*2*im*a/(barred ? ζb : ζ))*sθ) * L2
    Rval = _Rval(s, l, m, a, ω, r; bc=bc, λ=λ, barred=barred)
    pref = -(one(a)/4) * (1/(barred ? ζ : ζb)^2)
    pref * value(inner) * Rval
end

"""
    HIRG_nm(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing)

`H_nm^IRG = −(εg/(2√2))(1/ζ̄)( Dn[0, Lnadj[2, R·S]] + (2a²cosθ sinθ/(ζζ̄)) Dn[0,R]·S − (2r/(ζζ̄)) Lnadj[2,S]·R )`.
"""
function HIRG_nm(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = -2
    ζ = _ζfun(a, r, θ); ζb = _ζbarfun(a, r, θ); sθ = sin(θ); cθ = cos(θ)
    RS  = _RS(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    L2RS = Ln(-1, 2, RS; a=a, θ=θ, ω=ω, m=m)
    t1 = Dn(-1, 0, L2RS; a=a, r=r, ω=ω, m=m, bar=barred)
    Ronly = _Ronly(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    DnR = Dn(-1, 0, Ronly; a=a, r=r, ω=ω, m=m, bar=barred)
    Sval = _Sval(s, l, m, a, ω, θ; barred=barred)
    t2 = (2*a^2*cθ*sθ/(ζ*ζb)) * value(DnR) * Sval
    Sf = _Sonly(s, l, m, a, ω, r, θ; λ=λ, barred=barred)
    L2S = Ln(-1, 2, Sf; a=a, θ=θ, ω=ω, m=m)
    Rval = _Rval(s, l, m, a, ω, r; bc=bc, λ=λ, barred=barred)
    t3 = -(2*r/(ζ*ζb)) * value(L2S) * Rval
    pref = -(one(a)/(2*sqrt(oftype(real(ζ), 2)))) * (1/(barred ? ζ : ζb))
    pref * (value(t1) + t2 + t3)
end

"""
    HIRG_mm(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing)

`H_mm^IRG = −(εg/2)( Dn[0, Dn[0, R·S]] − (2/ζ) Dn[0, R·S] )`, s=−2.
"""
function HIRG_mm(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = -2
    ζ = _ζfun(a, r, θ); ζb = _ζbarfun(a, r, θ)
    RS  = _RS(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    D1  = Dn(-1, 0, RS; a=a, r=r, ω=ω, m=m, bar=barred)
    D2  = Dn(-1, 0, D1; a=a, r=r, ω=ω, m=m, bar=barred)
    inner = D2 + (-(2/(barred ? ζb : ζ))) * D1
    -(one(a)/2) * value(inner)
end

# ===== ORG components (s=+2) =================================================================================
"""
    HORG_ll(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing)

`H_ll^ORG = −(εg/4) ζ²( Ln[1, Ln[2, S]] − (2 i a/ζ) sinθ · Ln[2, S] ) R`, s=+2.
"""
function HORG_ll(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = +2
    ζ = _ζfun(a, r, θ); ζb = _ζbarfun(a, r, θ); sθ = sin(θ)
    Sf  = _Sonly(s, l, m, a, ω, r, θ; λ=λ, barred=barred)
    L2  = Ln(+1, 2, Sf; a=a, θ=θ, ω=ω, m=m)
    L12 = Ln(+1, 1, L2; a=a, θ=θ, ω=ω, m=m)
    isign = barred ? -1 : 1
    inner = L12 + (-(isign*2*im*a/(barred ? ζb : ζ))*sθ) * L2
    Rval = _Rval(s, l, m, a, ω, r; bc=bc, λ=λ, barred=barred)
    pref = -(one(a)/4) * (barred ? ζb : ζ)^2
    pref * value(inner) * Rval
end

"""
    HORG_lm(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing)

`H_lm^ORG = (εg/(4√2)) Δ (ζ²/ζ̄)( ζ̄ Dnadj[2, ζ̄⁻¹ Ln[2,R·S]] + (a² sin2θ/Σ) Dnadj[2,R]·S − (1/ζ) Ln[2,R·S] )`.
"""
function HORG_lm(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = +2
    ζ = _ζfun(a, r, θ); ζb = _ζbarfun(a, r, θ); Δ = _Δfun(a, r); Σ = _Σfun(a, r, θ)
    zA = barred ? ζb : ζ        # barred swaps ζ↔ζ̄ everywhere
    zB = barred ? ζ  : ζb
    C = typeof(complex(float(a*real(ω)*r*θ)))
    RS  = _RS(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    L2RS = Ln(+1, 2, RS; a=a, θ=θ, ω=ω, m=m)
    zBinv_rjet = barred ? _zinv_rjet(C, a, r, θ) : _zbarinv_rjet(C, a, r, θ)
    inner1f = mul_rjet(L2RS, zBinv_rjet)
    Dn2 = Dn(+1, 2, inner1f; a=a, r=r, ω=ω, m=m, bar=barred)
    t1 = zB * value(Dn2)
    Ronly = _Ronly(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    Dn2R = Dn(+1, 2, Ronly; a=a, r=r, ω=ω, m=m, bar=barred)
    Sval = _Sval(s, l, m, a, ω, θ; barred=barred)
    t2 = (a^2*sin(2*θ)/Σ) * value(Dn2R) * Sval
    t3 = -(1/zA) * value(L2RS)
    pref = (one(a)/(4*sqrt(oftype(real(ζ), 2)))) * Δ * (zA^2/zB)
    pref * (t1 + t2 + t3)
end

"""
    HORG_mm(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing)

`H_mm^ORG = −(εg/8)(ζ²/ζ̄²)( Δ Dnadj[1, Δ Dnadj[2, R]]·S − (2/ζ) Δ² Dnadj[2,R]·S )`, s=+2.
"""
function HORG_mm(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = +2
    ζ = _ζfun(a, r, θ); ζb = _ζbarfun(a, r, θ); Δ = _Δfun(a, r)
    zA = barred ? ζb : ζ; zB = barred ? ζ : ζb
    C = typeof(complex(float(a*real(ω)*r*θ)))
    Δj = _Delta_jet(C, a, r)
    Ronly = _Ronly(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    Dn2R = Dn(+1, 2, Ronly; a=a, r=r, ω=ω, m=m, bar=barred)
    ΔDn2R = mul_rjet(Dn2R, Δj)
    Dn1 = Dn(+1, 1, ΔDn2R; a=a, r=r, ω=ω, m=m, bar=barred)
    Sval = _Sval(s, l, m, a, ω, θ; barred=barred)
    t1 = Δ * value(Dn1) * Sval
    t2 = -(2/zA) * Δ^2 * value(Dn2R) * Sval
    -(one(a)/8) * (zA^2/zB^2) * (t1 + t2)
end

# ===== RTJet (analytic (r,θ) 2-jet) forms — same chains, full jet instead of value() ========================
"RTJet of R(r) (radial only)."
_R_rtjet(s,l,m,a,ω,r,θ;bc,λ,barred=false) = sepfield_rtjet(_Ronly(s,l,m,a,ω,r,θ;bc=bc,λ=λ,barred=barred))
"RTJet of S(θ) (angular only)."
_S_rtjet(s,l,m,a,ω,r,θ;λ,barred=false)    = sepfield_rtjet(_Sonly(s,l,m,a,ω,r,θ;λ=λ,barred=barred))

"""
    HIRG_nn_jet(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing) -> RTJet

Analytic `(r,θ)` 2-jet of [`HIRG_nn`](@ref) (value + ∂ᵣ/∂_θ/∂ᵣᵣ/∂_θθ/∂ᵣ_θ) for the
linearized-Einstein residual; the `.v` field equals the scalar `HIRG_nn` exactly.
"""
function HIRG_nn_jet(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = -2
    C = typeof(complex(float(a*real(ω)*r*θ)))
    ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ); sθ = _rt_sin(C,θ)
    Sf  = _Sonly(s, l, m, a, ω, r, θ; λ=λ, barred=barred)
    L2  = Ln(-1, 2, Sf;  a=a, θ=θ, ω=ω, m=m)
    L12 = Ln(-1, 1, L2;  a=a, θ=θ, ω=ω, m=m)
    L2j = sepfield_rtjet(L2); L12j = sepfield_rtjet(L12)
    zfac = barred ? ζb : ζ
    isign = barred ? -1 : 1
    inner = L12j + (rtconst(C(-isign*2*im*a))*sθ*rtinv(zfac))*L2j
    Rj = _R_rtjet(s,l,m,a,ω,r,θ;bc=bc,λ=λ,barred=barred)
    zden = barred ? ζ : ζb
    pref = rtconst(C(-one(a)/4)) * rtinv(zden*zden)
    pref * inner * Rj
end

"""
    HIRG_nm_jet(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing) -> RTJet

Analytic `(r,θ)` 2-jet of [`HIRG_nm`](@ref); the `.v` field equals the scalar `HIRG_nm` exactly.
"""
function HIRG_nm_jet(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = -2
    C = typeof(complex(float(a*real(ω)*r*θ)))
    ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ); sθ = _rt_sin(C,θ); cθ = _rt_cos(C,θ)
    RS  = _RS(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    L2RS = Ln(-1, 2, RS; a=a, θ=θ, ω=ω, m=m)
    t1f = Dn(-1, 0, L2RS; a=a, r=r, ω=ω, m=m, bar=barred)
    t1 = sepfield_rtjet(t1f)
    Ronly = _Ronly(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    DnR = Dn(-1, 0, Ronly; a=a, r=r, ω=ω, m=m, bar=barred)
    DnRj = sepfield_rtjet(DnR)
    Sj = _S_rtjet(s,l,m,a,ω,r,θ;λ=λ,barred=barred)
    t2 = (rtconst(C(2*a^2))*cθ*sθ*rtinv(ζ*ζb))*DnRj*Sj
    Sf = _Sonly(s, l, m, a, ω, r, θ; λ=λ, barred=barred)
    L2S = Ln(-1, 2, Sf; a=a, θ=θ, ω=ω, m=m)
    L2Sj = sepfield_rtjet(L2S)
    Rj = _R_rtjet(s,l,m,a,ω,r,θ;bc=bc,λ=λ,barred=barred)
    t3 = (rtconst(C(-2))*_rt_r(C,r)*rtinv(ζ*ζb))*L2Sj*Rj
    zden = barred ? ζ : ζb
    pref = rtconst(C(-one(a)/(2*sqrt(oftype(real(C(1)), 2))))) * rtinv(zden)
    pref * (t1 + t2 + t3)
end

"""
    HIRG_mm_jet(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing) -> RTJet

Analytic `(r,θ)` 2-jet of [`HIRG_mm`](@ref); the `.v` field equals the scalar `HIRG_mm` exactly.
"""
function HIRG_mm_jet(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = -2
    C = typeof(complex(float(a*real(ω)*r*θ)))
    ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ)
    RS  = _RS(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    D1  = Dn(-1, 0, RS; a=a, r=r, ω=ω, m=m, bar=barred)
    D2  = Dn(-1, 0, D1; a=a, r=r, ω=ω, m=m, bar=barred)
    D1j = sepfield_rtjet(D1); D2j = sepfield_rtjet(D2)
    zfac = barred ? ζb : ζ
    inner = D2j + (rtconst(C(-2))*rtinv(zfac))*D1j
    rtconst(C(-one(a)/2)) * inner
end

"""
    HORG_ll_jet(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing) -> RTJet

Analytic `(r,θ)` 2-jet of [`HORG_ll`](@ref); the `.v` field equals the scalar `HORG_ll` exactly.
"""
function HORG_ll_jet(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = +2
    C = typeof(complex(float(a*real(ω)*r*θ)))
    ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ); sθ = _rt_sin(C,θ)
    Sf  = _Sonly(s, l, m, a, ω, r, θ; λ=λ, barred=barred)
    L2  = Ln(+1, 2, Sf; a=a, θ=θ, ω=ω, m=m)
    L12 = Ln(+1, 1, L2; a=a, θ=θ, ω=ω, m=m)
    L2j = sepfield_rtjet(L2); L12j = sepfield_rtjet(L12)
    zfac = barred ? ζb : ζ
    isign = barred ? -1 : 1
    inner = L12j + (rtconst(C(-isign*2*im*a))*sθ*rtinv(zfac))*L2j
    Rj = _R_rtjet(s,l,m,a,ω,r,θ;bc=bc,λ=λ,barred=barred)
    zsq = barred ? (ζb*ζb) : (ζ*ζ)
    pref = rtconst(C(-one(a)/4)) * zsq
    pref * inner * Rj
end

"""
    HORG_lm_jet(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing) -> RTJet

Analytic `(r,θ)` 2-jet of [`HORG_lm`](@ref) (the ``∂_r ζ̄^{-1} = -ζ̄^{-2}`` term kept exact);
the `.v` field equals the scalar `HORG_lm` exactly.
"""
function HORG_lm_jet(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = +2
    C = typeof(complex(float(a*real(ω)*r*θ)))
    ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ); Δ = _rt_Delta(C,a,r); Σ = _rt_Sigma(C,a,r,θ)
    zA = barred ? ζb : ζ
    zB = barred ? ζ  : ζb
    RS  = _RS(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    L2RS = Ln(+1, 2, RS; a=a, θ=θ, ω=ω, m=m)
    # Dnadj[2, ζ̄⁻¹·h] = ζ̄⁻¹·Dnadj[2,h] + (∂_r ζ̄⁻¹)·h, ∂_r ζ̄⁻¹ = −ζ̄⁻² (keeps the θ-dependence of ζ̄⁻¹ exact).
    Dn2h = Dn(+1, 2, L2RS; a=a, r=r, ω=ω, m=m, bar=barred)
    Dn2hj = sepfield_rtjet(Dn2h)
    L2RSj = sepfield_rtjet(L2RS)
    zBinvJ = rtinv(zB)
    dr_zBinv = rtconst(C(-1))*(zBinvJ*zBinvJ)
    Dn2j = zBinvJ*Dn2hj + dr_zBinv*L2RSj
    t1 = zB * Dn2j
    Ronly = _Ronly(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    Dn2R = Dn(+1, 2, Ronly; a=a, r=r, ω=ω, m=m, bar=barred)
    Dn2Rj = sepfield_rtjet(Dn2R)
    Sj = _S_rtjet(s,l,m,a,ω,r,θ;λ=λ,barred=barred)
    s2θ = _rt_sin2(C,θ)
    t2 = (rtconst(C(a^2))*s2θ*rtinv(Σ))*Dn2Rj*Sj
    t3 = (rtconst(C(-1))*rtinv(zA))*L2RSj
    pref = rtconst(C(one(a)/(4*sqrt(oftype(real(C(1)), 2))))) * Δ * (zA*zA*rtinv(zB))
    pref * (t1 + t2 + t3)
end

"""
    HORG_mm_jet(a, ω, l, m, r, θ; barred=false, bc=:Up, λ=nothing) -> RTJet

Analytic `(r,θ)` 2-jet of [`HORG_mm`](@ref); the `.v` field equals the scalar `HORG_mm` exactly.
"""
function HORG_mm_jet(a, ω, l::Integer, m::Integer, r, θ; barred::Bool=false, bc::Symbol=:Up, λ=nothing)
    s = +2
    C = typeof(complex(float(a*real(ω)*r*θ)))
    ζ = _rt_zeta(C,a,r,θ); ζb = _rt_zetabar(C,a,r,θ); Δ = _rt_Delta(C,a,r)
    zA = barred ? ζb : ζ; zB = barred ? ζ : ζb
    Δj = _Delta_jet(C, a, r)
    Ronly = _Ronly(s, l, m, a, ω, r, θ; bc=bc, λ=λ, barred=barred)
    Dn2R = Dn(+1, 2, Ronly; a=a, r=r, ω=ω, m=m, bar=barred)
    ΔDn2R = mul_rjet(Dn2R, Δj)
    Dn1 = Dn(+1, 1, ΔDn2R; a=a, r=r, ω=ω, m=m, bar=barred)
    Dn1j = sepfield_rtjet(Dn1); Dn2Rj = sepfield_rtjet(Dn2R)
    Sj = _S_rtjet(s,l,m,a,ω,r,θ;λ=λ,barred=barred)
    t1 = Δ * Dn1j * Sj
    t2 = (rtconst(C(-2))*rtinv(zA))*(Δ*Δ)*Dn2Rj*Sj
    pref = rtconst(C(-one(a)/8)) * (zA*zA*rtinv(zB*zB))
    pref * (t1 + t2)
end

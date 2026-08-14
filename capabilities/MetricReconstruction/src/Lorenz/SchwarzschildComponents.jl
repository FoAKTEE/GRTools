# =============================================================================================================
# SchwarzschildComponents — a=0 FULL METRIC RECONSTRUCTION: the 10 tetrad metric-perturbation radial modes
# assembled from the validated radial functions (SchwarzschildRadial: P_{-2},P_{-1},h,κ), converted to
# Barack-Lousto-Sago (BLS) variables h̄^i (legacy KerrLorenzCirc/components.jl).
#
# The 10 tetrad components (DDKW Eq.128), spin-weighted-spherical mode order:
#   [ h_{l+l+}, h_{l-l-}, h_{m+m+}, h_{m-m-}, ρ h_{l+m+}, ρc h_{l+m-}, ρc h_{l-m+}, ρ h_{l-m-}, ΣΔ h_{l+l-}, h ]
# each = s=2 + s=1 + s=0 contributions. At a=0 the spheroidal→spherical projection is DIAGONAL in ℓ (b=I).
#
# This is the a=0 standalone gauge-field surface AND the a→0 tetrad/BLS regression oracle for KerrComponents.
# Towers are ANALYTIC (FD-free): each base homogeneous field P satisfies its 2nd-order Teukolsky/scalar radial
# ODE, so P'' is analytic in (P,P'); `field_tower5` builds (P,P',P'',P''',P'''') by exact recurrence seeded by
# the SchwarzschildRadial exact (val,der). κ is inhomogeneous (source S=½r²h) → `kappa_tower5`. Type-generic {T}
# (BigFloat-clean when the radial reconstruct ran in BigFloat).
# =============================================================================================================

import Perturbation: TeukolskyRadial, teukolsky_radial_deriv
import Printf: @printf

export Recon, build_recon, tetrad_modes, to_BLS

const _AA0 = 0.0
const _SIGDELTA_DEBUG = Ref(false)   # diagnostic: print the 3 ΣΔ sector values (s2,s1,s0) per (l,r)
# c=iK/Δ and r-derivatives (a=0): K=ωr², Δ=r²-2r.
cc(r,ω) = (Δ=r^2-2r; im*ω*r^2/Δ)
function c_tower(r,ω)   # returns c, c', c'', c'''
    Δ=r^2-2r; Δp=2r-2; Δpp=2.0
    K=ω*r^2; Kp=2ω*r; Kpp=2ω
    u=im*K; up=im*Kp; upp=im*Kpp; uppp=zero(u)
    c   = u/Δ
    cp  = (up*Δ - u*Δp)/Δ^2
    cpp = (upp*Δ - u*Δpp)/Δ^2 - 2Δp*(up*Δ - u*Δp)/Δ^3
    w   = 1/Δ
    wp  = -Δp/Δ^2
    wpp = (2Δp^2 - Δpp*Δ)/Δ^3
    wppp= (6Δp*Δpp*Δ - 6Δp^3)/Δ^4
    cppp = uppp*w + 3*upp*wp + 3*up*wpp + u*wppp
    (c,cp,cpp,cppp)
end

# =============================================================================================================
# ANALYTIC DERIVATIVE TOWER.  Length-4 scalar-jet (g,g',g'',g''') algebra for the smooth ODE coefficients, then
# the field α,β jets (P''=αP'+βP) per spin, then `field_tower5`/`kappa_tower5` building the exact 5-tuple from
# the radial exact seed (val,der).  Type-generic over Complex{T}.
# =============================================================================================================
# j4mul/j4add/j4scale/j4inv — the length-4 tuple jet kernel. DEDUP (F2): these are now thin ALIASES to
# the canonical NumericsCore jet algebra (`using`-imported into Lorenz via Lorenz/Operators.jl:25), so
# NO jet DEFINITION lives outside NumericsCore anymore:
#   • j4add/j4scale reuse the untyped `Vector` methods, which already act on tuples elementwise
#     (`f .+ g`, `s .* f` return tuples) — bit-exact to the former hand-unrolled bodies, zero-alloc.
#   • j4mul/j4inv resolve to the allocation-free `NTuple` methods added to NumericsCore/src/jets.jl,
#     emitting the SAME Leibniz recurrence term-for-term (`@generated`, bit-exact + zero-heap; the
#     `NTuple` `jinv` deliberately keeps THIS kernel's `inv(f[1])` reciprocal convention, distinct from
#     the `Vector` `jinv`'s `one/f[1]` form — see the jets.jl section note). Characterized bit-exact vs
#     the former length-4 tuple bodies in core/NumericsCore/test/test_jets.jl and gated here in
#     test/lorenz_schwarzschild_components.jl (baked `_legacy_j4*` oracles + workload === fixture).
const j4mul   = jmul
const j4add   = jadd
const j4scale = jscale
const j4inv   = jinv

# coefficient length-4 jets at r (a=0). Δ,K polynomial; c=iK/Δ.  CT is the (complex) element type.
function _coef4(r, ω, ::Type{CT}) where {CT}
    T = real(CT)
    rT = T(r); ωT = T(ω)
    Δ  = (CT(rT^2-2rT), CT(2rT-2), CT(2), CT(0))
    K  = (CT(ωT*rT^2),  CT(2ωT*rT), CT(2ωT), CT(0))
    invΔ = j4inv(Δ)
    c  = j4scale(CT(im), j4mul(K, invΔ))
    Δp = (Δ[2], Δ[3], Δ[4], CT(0))
    Kp = (K[2], K[3], K[4], CT(0))
    (Δ=Δ, Δp=Δp, K=K, Kp=Kp, invΔ=invΔ, c=c)
end

# α,β length-4 jets for spin ∈ {-2,-1,0,+2}: P''=αP'+βP from Δ P'' = A1 P' + A0 P.  eig = Λ (s=±2) / λ (s=−1) /
# λ0 (s=0).  Returns also the inhomogeneous γ-jet for κ (s=0 with a UNIT source coeff 1/Δ, since Δκ''=A1κ'+A0κ+S).
function _alphabeta4(spin, r, ω, eig, ::Type{CT}) where {CT}
    co = _coef4(r, ω, CT)
    Δ=co.Δ; Δp=co.Δp; c=co.c; invΔ=co.invΔ; K=co.K; Kp=co.Kp
    cp = (c[2], c[3], c[4], CT(0))
    csq = j4mul(c,c)
    rj  = (CT(r), CT(1), CT(0), CT(0))
    eigj = (CT(eig), CT(0), CT(0), CT(0))
    Z4 = (CT(0),CT(0),CT(0),CT(0))
    if spin == -2
        A1 = Δp
        A0 = j4add(j4add(j4mul(Δ,cp), j4mul(Δ,csq)),
                   j4add(j4scale(CT(-1), j4mul(Δp,c)),
                         j4add(j4scale(CT(6im*ω), rj), eigj)))
    elseif spin == -1
        A1 = Z4
        A0 = j4add(j4add(j4mul(Δ,cp), j4mul(Δ,csq)),
                   j4add(j4scale(CT(2im*ω), rj), eigj))
    elseif spin == 0
        A1 = j4scale(CT(-1), Δp)
        KKoverΔ = j4mul(j4mul(K,K), invΔ)
        A0 = j4scale(CT(-1), j4add(j4add(j4scale(CT(im), Kp), KKoverΔ),
                                   j4add(j4scale(CT(-2im*ω), rj), j4scale(CT(-1),eigj))))
    elseif spin == 2
        A1 = Δp
        A0 = j4add(j4add(j4scale(CT(-1),j4mul(Δ,cp)), j4mul(Δ,csq)),
                   j4add(j4mul(Δp,c),
                         j4add(j4scale(CT(-6im*ω), rj), eigj)))
    else
        error("_alphabeta4: spin=$spin unsupported")
    end
    α = j4mul(A1, invΔ); β = j4mul(A0, invΔ)
    (α, β, invΔ)   # invΔ = the κ source-coefficient jet (1/Δ) for the inhomogeneous recurrence
end

# Build the exact (P,P',P'',P''',P'''') 5-tuple of a HOMOGENEOUS field from (P0,P1)=(val,der) + (α,β) jets.
@inline function _tower5(P0, P1, α, β)
    P2 = α[1]*P1 + β[1]*P0
    P3 = α[2]*P1 + α[1]*P2 + β[2]*P0 + β[1]*P1
    P4 = α[3]*P1 + 2α[2]*P2 + α[1]*P3 + β[3]*P0 + 2β[2]*P1 + β[1]*P2
    (P0,P1,P2,P3,P4)
end

# Analytic 5-tuple tower of a homogeneous field `f` (HomField) at r, spin `spin`, eigenvalue `eig`.
function field_tower5(f, spin, ω, eig, r)
    P0 = val(f,r); P1 = der(f,r)
    CT = typeof(complex(P0))
    α, β, _ = _alphabeta4(spin, r, ω, eig, CT)
    _tower5(CT(P0), CT(P1), α, β)
end

# Analytic 5-tuple tower of the INHOMOGENEOUS κ field (s=0 op, source S=½r²h).
function kappa_tower5(kf, hT5, ω, λ0, r)
    I1, I2 = kappa_I(kf, r)
    κ0v = val(kf.hom, r) + (kf.brs0.inn(r)*I1 + kf.brs0.up(r)*I2)
    κ1v = der(kf.hom, r) + (kf.brs0.ind(r)*I1 + kf.brs0.upd(r)*I2)
    CT = typeof(complex(κ0v))
    κ0 = CT(κ0v); κ1 = CT(κ1v)
    α, β, srcc = _alphabeta4(0, r, ω, λ0, CT)   # srcc = (1/Δ) jet ; recurrence κ''=ακ'+βκ + srcc·S
    T = real(CT); rT = T(r)
    h0,h1,h2,h3,_ = hT5
    half_r2 = (CT(rT^2/2), CT(rT), CT(1), CT(0))   # ½r² jet
    hjet = (CT(h0),CT(h1),CT(h2),CT(h3))
    Sjet = j4mul(half_r2, hjet)                    # S, S', S'', S'''
    g = j4mul(srcc, Sjet)                          # γ·S and its derivatives (γ=1/Δ)
    g0,g1,g2,_ = g
    κ2 = α[1]*κ1 + β[1]*κ0 + g0
    κ3 = α[2]*κ1 + α[1]*κ2 + β[2]*κ0 + β[1]*κ1 + g1
    κ4 = α[3]*κ1 + 2α[2]*κ2 + α[1]*κ3 + β[3]*κ0 + 2β[2]*κ1 + β[1]*κ2 + g2
    (κ0,κ1,κ2,κ3,κ4)
end

# ---- Chandrasekhar operator application on a 5-tuple tower (exact, FD-free) ----
function D_apply(t, c4, sgn)        # c4=(c,c',c'',c'''); sgn=+1→𝒟 (-c), -1→𝒟† (+c)
    f0,f1,f2,f3,f4 = t
    c,cp,cpp,cppp = c4
    s = sgn==+1 ? -one(c) : +one(c)
    g0 = f1 + s*c*f0
    g1 = f2 + s*(cp*f0 + c*f1)
    g2 = f3 + s*(cpp*f0 + 2cp*f1 + c*f2)
    g3 = f4 + s*(cppp*f0 + 3cpp*f1 + 3cp*f2 + c*f3)
    (g0,g1,g2,g3,zero(c))
end
function chain_tower(baseT, ω, ops)
    r->begin
        t = baseT(r)
        c4 = c_tower(r, ω)
        for op in reverse(ops)
            t = D_apply(t, c4, op==:D ? +1 : -1)
        end
        t
    end
end
chainval(baseT, ω, ops, r) = first(chain_tower(baseT, ω, ops)(r))

"""
    Recon{T}

a=0 reconstruction bundle at `(r0, l, m)` built by [`build_recon`](@ref): the orbit/eigenvalue data
(`ω, σ, Λ, λ, Λ̂, λ̂, λ̂2, C, B, λ0`), the `SchwarzschildRadial` `reconstruct(...)` result, and the
ANALYTIC 5-tuple operator-tower closures ``r ↦ (f, f', f'', f''', f'''')`` for
``P_{-2}, P_{-1}, h, κ, P_{+2}`` that feed [`tetrad_modes`](@ref) / [`to_BLS`](@ref).
Legacy: `KerrLorenzCirc/components.jl`.
"""
struct Recon{T<:Real}
    r0::T; l::Int; m::Int; ω::T; σ::T
    Λ::T; λ::T; Λ̂::T; λ̂::T; λ̂2::T; C::Complex{T}; B::T
    λ0::T
    rec               # the SchwarzschildRadial reconstruct(...) result (Pm2,Pm1,h,κ HomFields)
    Pm2T; Pm1T; hT; κT; Pp2T   # ANALYTIC 5-tuple tower closures r->(f,f',f'',f''',f'''')
end

# +2 Weyl jumps (a=0), same closed form as jumps_pm(sgn=+1); element type from r0.
function plus2_jumps(r0, l, m)
    bg = orbit_background(0, r0, m)
    J0,J1 = jumps_pm(bg, l, +1)
    CT = Complex{typeof(float(real(r0)))}
    (CT(J0), CT(J1))
end

# An analytic P_{+2}=Δ²R_{+2} accessor, reconstructed from the s=+2 TeukolskyRadial with its own +2 jumps.
struct Plus2Field{T<:Real}
    r0::T; Aup::Complex{T}; Ain::Complex{T}; Rp2; a::T
end
function _Δ2R(pf::Plus2Field, r; leg)
    Δ=r^2-2r+pf.a^2
    R = leg===:up ? pf.Rp2.Up(r) : pf.Rp2.In(r)
    Δ^2*R
end
function _Δ2R_d(pf::Plus2Field, r; leg)
    Δ=r^2-2r+pf.a^2; Δp=2r-2
    R  = leg===:up ? pf.Rp2.Up(r) : pf.Rp2.In(r)
    Rp = leg===:up ? teukolsky_radial_deriv(pf.Rp2.Up,r) : teukolsky_radial_deriv(pf.Rp2.In,r)
    2Δ*Δp*R + Δ^2*Rp
end
p2val(pf::Plus2Field,r) = r<pf.r0 ? pf.Ain*_Δ2R(pf,r;leg=:inn) : pf.Aup*_Δ2R(pf,r;leg=:up)
p2der(pf::Plus2Field,r) = r<pf.r0 ? pf.Ain*_Δ2R_d(pf,r;leg=:inn) : pf.Aup*_Δ2R_d(pf,r;leg=:up)

# An analytic P_{-2} field for the BigFloat path (rebuilds the s=-2 leg at the active BigFloat precision; the
# Float64 path keeps the SchwarzschildRadial HomField).
struct Minus2Field{T<:Real}
    r0::T; Aup::Complex{T}; Ain::Complex{T}; Rm2; a::T
end
m2val(f::Minus2Field,r) = r<f.r0 ? f.Ain*f.Rm2.In(r) : f.Aup*f.Rm2.Up(r)
m2der(f::Minus2Field,r) = r<f.r0 ? f.Ain*teukolsky_radial_deriv(f.Rm2.In,r) : f.Aup*teukolsky_radial_deriv(f.Rm2.Up,r)
function field_tower5(f::Minus2Field, spin, ω, eig, r)
    P0 = m2val(f,r); P1 = m2der(f,r)
    CT = typeof(complex(P0))
    α, β, _ = _alphabeta4(spin, r, ω, eig, CT)
    _tower5(CT(P0), CT(P1), α, β)
end

"""
    build_recon(r0, l; m=2) -> Recon

Build the ANALYTIC operator towers (P_{-2},P_{-1},h,κ,P_{+2}) for the a=0 metric reconstruction at (r0,l).
"""
function build_recon(r0, l; m=2)
    T = typeof(float(real(r0)))
    r0T = T(r0)
    rec = reconstruct(r0T, l; m=m)
    ω = orbit_background(0, r0T, m).ω; σ = T((-1)^(l+m))
    Λ = T((l-1)*(l+2)); λ = T(l*(l+1)); λ0 = T(l*(l+1))
    Λ̂v = sqrt(T((l-1)*l*(l+1)*(l+2))); λ̂v = sqrt(T(l*(l+1))); λ̂2v = sqrt(T((l-1)*(l+2)))
    C = Λ*(Λ+2) + 12im*ω*σ; B = λ
    Pm2T = if T <: BigFloat
        Rm2 = TeukolskyRadial(-2, l, m, T(_AA0), ω)
        J0m2, J1m2 = rec.jd.Pm2
        Mu = Rm2.Up(r0T); Mi = Rm2.In(r0T)
        Mud = teukolsky_radial_deriv(Rm2.Up, r0T); Mid = teukolsky_radial_deriv(Rm2.In, r0T)
        detm = Mu*(-Mid) - (-Mi)*Mud
        Aup_m = (J0m2*(-Mid) - (-Mi)*J1m2)/detm
        Ain_m = (Mu*J1m2 - Mud*J0m2)/detm
        m2f = Minus2Field{T}(r0T, Complex{T}(Aup_m), Complex{T}(Ain_m), Rm2, T(_AA0))
        r->field_tower5(m2f, -2, ω, Λ, r)
    else
        r->field_tower5(rec.Pm2, -2, ω, Λ, r)
    end
    Pm1T = r->field_tower5(rec.Pm1, -1, ω, λ,  r)
    hT   = r->field_tower5(rec.h,    0, ω, λ0, r)
    κT   = r->kappa_tower5(rec.κ, field_tower5(rec.h,0,ω,λ0,r), ω, λ0, r)
    Rp2 = TeukolskyRadial(2, l, m, T(_AA0), ω)
    J0p2,J1p2 = plus2_jumps(r0T, l, m)
    P2up(r)=Rp2.Up(r)*(r^2-2r)^2; P2in(r)=Rp2.In(r)*(r^2-2r)^2
    P2upd(r)=(Δ=r^2-2r;Δp=2r-2; 2Δ*Δp*Rp2.Up(r)+Δ^2*teukolsky_radial_deriv(Rp2.Up,r))
    P2ind(r)=(Δ=r^2-2r;Δp=2r-2; 2Δ*Δp*Rp2.In(r)+Δ^2*teukolsky_radial_deriv(Rp2.In,r))
    det = P2up(r0T)*(-P2ind(r0T)) - (-P2in(r0T))*P2upd(r0T)
    Aup2 = (J0p2*(-P2ind(r0T)) - (-P2in(r0T))*J1p2)/det
    Ain2 = (P2up(r0T)*J1p2 - P2upd(r0T)*J0p2)/det
    CT = Complex{T}
    pf = Plus2Field{T}(r0T, CT(Aup2), CT(Ain2), Rp2, T(_AA0))
    Pp2T = r->field_tower5(pf, 2, ω, Λ, r)
    Recon{T}(r0T,l,m,ω,σ,Λ,λ,Λ̂v,λ̂v,λ̂2v,C,B,λ0,rec,Pm2T,Pm1T,hT,κT,Pp2T)
end

# field_tower5 specialization for Plus2Field (seed from p2val/p2der).
function field_tower5(f::Plus2Field, spin, ω, eig, r)
    P0 = p2val(f,r); P1 = p2der(f,r)
    CT = typeof(complex(P0))
    α, β, _ = _alphabeta4(spin, r, ω, eig, CT)
    _tower5(CT(P0), CT(P1), α, β)
end

function scaleΔ2_tower(g4tuple_fn, C, r)
    Δ=r^2-2r; Δp=2r-2; Δpp=2.0
    A=Δ^2/C; Ap=2Δ*Δp/C; App=(2Δp^2+2Δ*Δpp)/C
    g=g4tuple_fn(r); g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; g2=length(g)>=3 ? g[3] : 0.0+0im
    (A*g0, Ap*g0+A*g1, App*g0+2Ap*g1+A*g2)
end

"""
    tetrad_modes(R::Recon, r) -> NamedTuple

The 10 tetrad radial-mode coefficients at r (a=0): lplp,lmlm,mpmp,mmmm,ρlpmp,ρclpmm,ρclmmp,ρlmmm,ΣΔlplm,h.
All operator chains use the ANALYTIC towers — FD-free.
"""
function tetrad_modes(R::Recon, r)
    ω=R.ω; σ=R.σ; Δ=r^2-2r; Λ̂=R.Λ̂; λ̂=R.λ̂; λ̂2=R.λ̂2; β=1/(-6im*ω)

    # Same-r base-tower memo (byte-equal): evaluate each base tower ONCE at this r and reuse across chains.
    _memo(bT) = let tcache = bT(r), rfix = r
        _rrq -> (_rrq == rfix ? tcache : bT(_rrq))
    end
    Pm2Tm = _memo(R.Pm2T); Pp2Tm = _memo(R.Pp2T); Pm1Tm = _memo(R.Pm1T); hTm = _memo(R.hT); κTm = _memo(R.κT)

    DDPm2     = chainval(Pm2Tm, ω, [:D,:D], r)
    DPm2      = chainval(Pm2Tm, ω, [:D], r)
    DDDPm2    = chainval(Pm2Tm, ω, [:D,:D,:D], r)
    DdDDPm2   = chainval(Pm2Tm, ω, [:Dd,:D,:D], r)
    Pp2_5 = Pp2Tm
    Pp2  = Pp2Tm(r)[1]
    DdPp2   = chainval(Pp2_5, ω, [:Dd], r)
    DdDdPp2 = chainval(Pp2_5, ω, [:Dd,:Dd], r)
    DdDdDdPp2 = chainval(Pp2_5, ω, [:Dd,:Dd,:Dd], r)
    DDdDdPp2  = chainval(Pp2_5, ω, [:D,:Dd,:Dd], r)

    DPm1   = chainval(Pm1Tm, ω, [:D], r)
    Ddm1Pm1 = (t=Pm1Tm(r); cc_=cc(r,ω); Δp=2r-2; (t[2]+cc_*t[1]) + (-1)*(Δp/Δ)*t[1])
    D2Pm1_tower = chain_tower(Pm1Tm, ω, [:D,:D])
    Pp1T = r2->begin
        Δr=r2^2-2r2; Δrp=2r2-2; g=D2Pm1_tower(r2); g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im
        (Δr/R.B*g0, (Δrp*g0+Δr*g1)/R.B)
    end
    Pp1_5 = r2->begin
        Δr=r2^2-2r2; Δrp=2r2-2; Δrpp=2.0; g=D2Pm1_tower(r2)
        g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; g2=length(g)>=3 ? g[3] : 0.0+0im
        ( Δr/R.B*g0, (Δrp*g0+Δr*g1)/R.B, (Δrpp*g0+2Δrp*g1+Δr*g2)/R.B, zero(g0), zero(g0) )
    end
    Pp1  = Pp1T(r)[1]
    DPp1 = chainval(Pp1_5, ω, [:D], r)
    DdPp1 = chainval(Pp1_5, ω, [:Dd], r)
    Dm1Pp1 = (t=Pp1_5(r); cc_=cc(r,ω); Δp=2r-2; (t[2]-cc_*t[1]) + (-1)*(Δp/Δ)*t[1])

    ht = hTm(r); h0=ht[1]; Dh=chainval(hTm,ω,[:D],r); Ddh=chainval(hTm,ω,[:Dd],r)
    κt0 = κTm(r); κ0 = κt0[1]
    κT = κTm
    DDκ   = chainval(κT, ω, [:D,:D], r)
    DdDdκ = chainval(κT, ω, [:Dd,:Dd], r)
    Dκ    = chainval(κT, ω, [:D], r)
    Ddκ   = chainval(κT, ω, [:Dd], r)
    Dh_tower = chain_tower(hTm, ω, [:D])
    rDh_tower = r2->(g=Dh_tower(r2); g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (r2*g0, g0+r2*g1, zero(g0),zero(g0),zero(g0)))
    DrDh = chainval(rDh_tower, ω, [:D], r)
    Ddh_tower = chain_tower(hTm, ω, [:Dd])
    rDdh_tower = r2->(g=Ddh_tower(r2); g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (r2*g0, g0+r2*g1, zero(g0),zero(g0),zero(g0)))
    DdrDdh = chainval(rDdh_tower, ω, [:Dd], r)
    Ddκ_tower = chain_tower(κT, ω, [:Dd])
    ΔDdκ_tower = r2->(g=Ddκ_tower(r2); Δr=r2^2-2r2; Δrp=2r2-2; g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (Δr*g0, Δrp*g0+Δr*g1, zero(g0),zero(g0),zero(g0)))
    DΔDdκ = chainval(ΔDdκ_tower, ω, [:D], r)
    Dκ_tower = chain_tower(κT, ω, [:D])
    ΔDκ_tower = r2->(g=Dκ_tower(r2); Δr=r2^2-2r2; Δrp=2r2-2; g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (Δr*g0, Δrp*g0+Δr*g1, zero(g0),zero(g0),zero(g0)))
    DdΔDκ = chainval(ΔDκ_tower, ω, [:Dd], r)
    κp = κt0[2]

    # ============================ s=2 ============================
    hs2_lplp = -1/(6ω^2*Δ^2)*Λ̂*(1+σ)*Pp2
    hs2_lmlm = -1/(6ω^2*Δ^2)*Λ̂*(σ+1)*(Pm2Tm(r)[1])
    M_mpmp = (σ*DDPm2 + DdDdPp2)
    hs2_mpmp = -1/(6ω^2)*M_mpmp
    hs2_mmmm = -1/(6ω^2)*M_mpmp*σ
    Tlpmp = ( 1/(12ω^2)*(DDPm2*Λ̂ - DdDdPp2*Λ̂)*λ̂
            - 1/(12ω^2)*(DDDPm2*Λ̂ - DDdDdPp2*Λ̂)*(λ̂*r)
            - 1/(3im*ω)*(DDDPm2*λ̂2*r^2 - 2*DDPm2*λ̂2*r + 2*DPm2*λ̂2)
            - DdPp2/(3im*ω*Δ)*λ̂2*(R.λ̂^2*r^2) )
    ρhlpmp = σ*β*Tlpmp
    Tlpmm = ( -1/(12ω^2)*(DDPm2*Λ̂ - DdDdPp2*Λ̂)*λ̂
            + 1/(12ω^2)*(DDDPm2*Λ̂ - DDdDdPp2*Λ̂)*(λ̂*r)
            + 1/(3im*ω)*(DDDPm2*λ̂2*r^2 - 2*DDPm2*λ̂2*r + 2*DPm2*λ̂2)
            + DdPp2/(3im*ω*Δ)*λ̂2*(R.λ̂^2*r^2) )
    ρchlpmm = β*Tlpmm
    Tlmmp = ( 1/(12ω^2)*(DDPm2*Λ̂ - DdDdPp2*Λ̂)*λ̂
            - 1/(12ω^2)*(DdDDPm2*Λ̂ - DdDdDdPp2*Λ̂)*(λ̂*r)
            - 1/(3im*ω)*(DdDdDdPp2*λ̂2*r^2 - 2*DdDdPp2*λ̂2*r + 2*DdPp2*λ̂2)
            - DPm2/(3im*ω*Δ)*λ̂2*(R.λ̂^2*r^2) )
    ρchlmmp = β*Tlmmp
    Tlmmm = ( -1/(12ω^2)*(DDPm2*Λ̂ - DdDdPp2*Λ̂)*λ̂
            + 1/(12ω^2)*(DdDDPm2*Λ̂ - DdDdDdPp2*Λ̂)*(λ̂*r)
            + 1/(3im*ω)*(DdDdDdPp2*λ̂2*r^2 - 2*DdDdPp2*λ̂2*r + 2*DdPp2*λ̂2)
            + DPm2/(3im*ω*Δ)*λ̂2*(R.λ̂^2*r^2) )
    ρhlmmm = σ*β*Tlmmm
    DDPm2_tower = chain_tower(Pm2Tm, ω, [:D,:D])
    DdDDPm2_tower = chain_tower(Pm2Tm, ω, [:Dd,:D,:D])
    ΔDdDDPm2_tower = r2->(g=DdDDPm2_tower(r2); Δr=r2^2-2r2; Δrp=2r2-2; g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (Δr*g0,Δrp*g0+Δr*g1,zero(g0),zero(g0),zero(g0)))
    DΔDdDDPm2 = chainval(ΔDdDDPm2_tower, ω, [:D], r)
    ΔDDDPm2_tower = r2->(g=chain_tower(Pm2Tm,ω,[:D,:D,:D])(r2); Δr=r2^2-2r2; Δrp=2r2-2; g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (Δr*g0,Δrp*g0+Δr*g1,zero(g0),zero(g0),zero(g0)))
    DdΔDDDPm2 = chainval(ΔDDDPm2_tower, ω, [:Dd], r)
    LDD_Pm2 = DΔDdDDPm2 + DdΔDDDPm2
    drDDPm2 = DDPm2_tower(r)[2]
    brace1 = 1/(24ω^2)*( LDD_Pm2*Λ̂*r^2 - 4r*Δ*drDDPm2*Λ̂ )
    DdDdPp2_tower = chain_tower(Pp2_5, ω, [:Dd,:Dd])
    DDdDdPp2_tower = chain_tower(Pp2_5, ω, [:D,:Dd,:Dd])
    ΔDDdDdPp2_tower = r2->(g=DDdDdPp2_tower(r2); Δr=r2^2-2r2; Δrp=2r2-2; g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (Δr*g0,Δrp*g0+Δr*g1,zero(g0),zero(g0),zero(g0)))
    DdΔDDdDdPp2 = chainval(ΔDDdDdPp2_tower, ω, [:Dd], r)
    ΔDdDdDdPp2_tower = r2->(g=chain_tower(Pp2_5,ω,[:Dd,:Dd,:Dd])(r2); Δr=r2^2-2r2; Δrp=2r2-2; g0=g[1]; g1=length(g)>=2 ? g[2] : 0.0+0im; (Δr*g0,Δrp*g0+Δr*g1,zero(g0),zero(g0),zero(g0)))
    DΔDdDdDdPp2 = chainval(ΔDdDdDdPp2_tower, ω, [:D], r)
    LDD_Pp2 = DΔDdDdDdPp2 + DdΔDDdDdPp2
    drDdDdPp2 = DdDdPp2_tower(r)[2]
    brace2 = -1/(24ω^2)*( LDD_Pp2*Λ̂*r^2 - 4r*Δ*drDdDdPp2*Λ̂ )
    weyl_m = 1/(3im*ω)*( DDPm2*Λ̂*(r^2)*(r) - DPm2*Λ̂*(r^2) )
    weyl_p = 1/(3im*ω)*( DdDdPp2*Λ̂*(r^2)*(r) - DdPp2*Λ̂*(r^2) )
    ΣΔhlplm_s2 = β*( brace1 + brace2 + weyl_m + weyl_p )

    # ============================ s=1 ============================
    S0s1 = -(1+σ)*λ̂
    hs1_lplp = 2/Δ*(σ*Dm1Pp1)*S0s1
    hs1_lmlm = 2/Δ*(Ddm1Pm1)*S0s1
    hs1_mpmp = (σ*2*(DPm1 + σ*DdPp1))*(-λ̂2)
    hs1_mmmm = (-2*(DPm1 + σ*DdPp1))*(+λ̂2)
    DPm1_tower  = chain_tower(Pm1Tm, ω, [:D])
    DdPp1_tower = chain_tower(Pp1_5,  ω, [:Dd])
    Pcal_tower = r2->begin
        a = DPm1_tower(r2); b = DdPp1_tower(r2)
        (a[1]+σ*b[1], a[2]+σ*b[2], a[3]+σ*b[3], zero(a[1]), zero(a[1]))
    end
    Pcal = first(Pcal_tower(r))
    DPcal  = chainval(Pcal_tower, ω, [:D], r)
    DdPcal = chainval(Pcal_tower, ω, [:Dd], r)
    Pm1_0 = Pm1Tm(r)[1]
    ρhlpmp_s1 = σ*(r*DPcal - 2*Pcal)  + (σ/Δ)*Pp1*S0s1*(-r*λ̂)
    ρchlpmm_s1 = -(r*DPcal - 2*Pcal)  + (σ/Δ)*Pp1*S0s1*(r*λ̂)
    ρchlmmp_s1 = σ*(r*DdPcal - 2*Pcal) + (1/Δ)*Pm1_0*S0s1*(-r*λ̂)
    ρhlmmm_s1  = -(r*DdPcal - 2*Pcal) + (1/Δ)*Pm1_0*S0s1*(r*λ̂)
    ΣΔhlplm_s1 = Pcal*S0s1*r^2 - 2r*(Pm1_0 + σ*Pp1)*S0s1

    # ============================ s=0 ============================
    hs0_lplp = -( (1/(im*ω))*DrDh + 4*DDκ )
    hs0_lmlm = -( -(1/(im*ω))*DdrDdh + 4*DdDdκ )
    hs0_mpmp = -4*Λ̂*κ0
    hs0_mmmm = -4*Λ̂*κ0
    ρhlpmp_s0 = (1/(2im*ω))*Dh*λ̂*r^2 + 4*Dκ*λ̂*r - 4*κ0*λ̂
    ρchlpmm_s0 = -(1/(2im*ω))*Dh*λ̂*r^2 - 4*Dκ*λ̂*r + 4*κ0*λ̂
    ρchlmmp_s0 = -(1/(2im*ω))*Ddh*λ̂*r^2 + 4*Ddκ*λ̂*r - 4*κ0*λ̂
    ρhlmmm_s0  = (1/(2im*ω))*Ddh*λ̂*r^2 - 4*Ddκ*λ̂*r + 4*κ0*λ̂
    ΣΔhlplm_s0 = -h0*(r^2 - 2r^2)*(r^2) - 2*(DΔDdκ + DdΔDκ)*(r^2) + 8r*Δ*κp
    htrace = h0
    _SIGDELTA_DEBUG[] && @printf("ΣΔ sectors l=%d r=%.3f: |s2|=%.4e |s1|=%.4e |s0|=%.4e total=%.4e\n",
        R.l, Float64(r), Float64(abs(ΣΔhlplm_s2)), Float64(abs(ΣΔhlplm_s1)), Float64(abs(ΣΔhlplm_s0)),
        Float64(abs(ΣΔhlplm_s2+ΣΔhlplm_s1+ΣΔhlplm_s0)))

    (lplp = hs2_lplp+hs1_lplp+hs0_lplp,
     lmlm = hs2_lmlm+hs1_lmlm+hs0_lmlm,
     mpmp = hs2_mpmp+hs1_mpmp+hs0_mpmp,
     mmmm = hs2_mmmm+hs1_mmmm+hs0_mmmm,
     ρlpmp = ρhlpmp+ρhlpmp_s0+ρhlpmp_s1,
     ρclpmm = ρchlpmm+ρchlpmm_s0+ρchlpmm_s1,
     ρclmmp = ρchlmmp+ρchlmmp_s0+ρchlmmp_s1,
     ρlmmm = ρhlmmm+ρhlmmm_s0+ρhlmmm_s1,
     ΣΔlplm = ΣΔhlplm_s2 + ΣΔhlplm_s1 + ΣΔhlplm_s0,
     h = htrace)
end

"""
    to_BLS(R::Recon, r) -> NTuple{10}

The 10 Barack-Lousto-Sago variables h̄^1..h̄^10 at r (a=0), from the tetrad modes.
"""
function to_BLS(R::Recon, r)
    M = tetrad_modes(R, r)
    f=1-2/r; Δ=r^2-2r; Σ=r^2; Λ̂=R.Λ̂; λ̂=R.λ̂
    hlplm = M.ΣΔlplm/(Σ*Δ)
    h1 = 0.5*r*f^2*(M.lplp + M.lmlm)
    h2 = 0.5*r*f^2*(M.lplp - M.lmlm)
    hmpmm = Σ*M.h - Δ*hlplm
    h3 = hmpmm/r
    hlpmp = M.ρlpmp/r; hlpmm = M.ρclpmm/r; hlmmp = M.ρclmmp/r; hlmmm = M.ρlmmm/r
    h4 = -0.5*λ̂*f*( hlpmp - hlpmm - hlmmp + hlmmm )
    h5 = -0.5*λ̂*f*( hlpmp - hlpmm + hlmmp - hlmmm )
    h6 = -r*f*hlplm
    h7 = (Λ̂/(2r))*(M.mpmp + M.mmmm)
    h8 = -0.5im*λ̂*f*( hlpmp + hlpmm - hlmmp - hlmmm )
    h9 = -0.5im*λ̂*f*( hlpmp + hlpmm + hlmmp + hlmmm )
    h10 = (im*Λ̂/(2r))*(M.mpmp - M.mmmm)
    (h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)
end

# =============================================================================================================
# KerrSMatrices — a≠0 closure S-matrices + the radial operator chains R^{(s)}_{XX}[l](r) that mix spheroidal-l
# modes for the lplp/mpmp tetrad components (legacy KerrLorenzCirc/components_kerr.jl, the smatrices + lplp/mpmp
# block). At a≠0 the spheroidal→spherical projection is banded: h^{(s),l'}_{XX}=Σ_l S^{(s)}_{XX}(l,l')·R^{(s)}_{XX}[l](r).
#
# `smatrices_kerr` are the SAME App-H banded S-entries the (validated) closure_kerr assembles (App-H a=0.99 lower
# 6.85e-8); here exposed as standalone (l,l')→Complex functions. `radial_chain_kerr` evaluates the closure radial
# pieces as FD-free analytic functions of r on a reconstruct_kerr result. The spin-general jet towers
# (_field_tower5K/_kappa_tower5K/_chain_towerK) reuse the SchwarzschildComponents length-4 jet algebra (j4*/_tower5).
# `_attach_recmeta` adds P_{+2} (s=+2 reconstructed) + the (l,B) metadata the component chains consume.
# =============================================================================================================

import Harmonics: SpinWeightedSpheroidalHarmonicS, expansion_coefficients

export smatrices_kerr, radial_chain_kerr, tetrad_lplp_kerr, tetrad_mpmp_kerr

# spin-general c=iK/Δ and its r-derivatives c,c',c'',c''' (analytic, FD-free).
function c_tower_kerr(r, a, ω, m)
    Δ=r^2-2r+a^2; Δp=2r-2; Δpp=2.0
    K=ω*(r^2+a^2)-a*m; Kp=2ω*r; Kpp=2ω
    u=im*K; up=im*Kp; upp=im*Kpp; uppp=zero(u)
    w=1/Δ; wp=-Δp/Δ^2; wpp=(2Δp^2-Δpp*Δ)/Δ^3; wppp=(6Δp*Δpp*Δ-6Δp^3)/Δ^4
    c   = u*w
    cp  = up*w + u*wp
    cpp = upp*w + 2up*wp + u*wpp
    cppp= uppp*w + 3upp*wp + 3up*wpp + u*wppp
    (c,cp,cpp,cppp)
end

# Apply 𝒟 (sgn=+1 ⇒ −c) or 𝒟† (sgn=−1 ⇒ +c) to a 5-tuple derivative tower (f,f',f'',f''',f'''').
function D_apply_k(t, c4, sgn)
    f0,f1,f2,f3,f4 = t
    c,cp,cpp,cppp = c4
    s = sgn==+1 ? -one(c) : +one(c)
    g0 = f1 + s*c*f0
    g1 = f2 + s*(cp*f0 + c*f1)
    g2 = f3 + s*(cpp*f0 + 2cp*f1 + c*f2)
    g3 = f4 + s*(cppp*f0 + 3cpp*f1 + 3cp*f2 + c*f3)
    (g0,g1,g2,g3,zero(c))
end

# metadata accessors stored on rec via _attach_recmeta (l and the Teukolsky-Starobinsky B). jd_l is the public
# spheroidal-l accessor (also used by LorenzWeyl); falls back with a helpful error for a bare reconstruct_kerr rec.
rec_l(rec) = rec._l
jd_B(rec) = rec._B
jd_l(rec) = hasproperty(rec, :_l) ? rec._l :
    error("pass an _attach_recmeta'd rec, or use weyl_psi4_lorenz(rec,r,θ;l=…)")

"""
    radial_chain_kerr(rec, r) -> NamedTuple

The radial-operator chains for spheroidal-l at radius r on a `reconstruct_kerr` result `rec`: the six
contributions entering h_{l+l+} and h_{m+m+} BEFORE the angular S-matrix contraction. FD-free for the leading
orders (analytic spin-general c-tower).
"""
function radial_chain_kerr(rec, r)
    jd = rec.jd; ω = jd.ω; a = jd.a; m = jd.m; σ = (-1)^(rec_l(rec)+m)
    T = typeof(real(ω)); CT = Complex{T}
    Δ = r^2 - 2r + a^2
    c4 = c_tower_kerr(r, a, ω, m)
    c, cp, cpp, cppp = c4

    Pm2  = valK(rec.Pm2, r); Pm2p = derK(rec.Pm2, r); Pm2pp = der2K(rec.Pm2, r)
    Pp2  = valK(rec.Pp2, r); Pp2p = derK(rec.Pp2, r); Pp2pp = der2K(rec.Pp2, r)

    DPm2  = Pm2p - c*Pm2
    DPm2p = Pm2pp - cp*Pm2 - c*Pm2p
    DDPm2 = DPm2p - c*DPm2
    DdPp2  = Pp2p + c*Pp2
    DdPp2p = Pp2pp + cp*Pp2 + c*Pp2p
    DdDdPp2 = DdPp2p + c*DdPp2

    s2lplp = -1/(6ω^2)*(1/Δ^2)*Pp2
    s2mpmp = -1/(6ω^2)*(σ*DDPm2 + DdDdPp2)

    B = jd_B(rec)
    Pm1  = valK(rec.Pm1, r); Pm1p = derK(rec.Pm1, r); Pm1pp = der2K(rec.Pm1, r)
    DPm1  = Pm1p - c*Pm1
    DPm1p = Pm1pp - cp*Pm1 - c*Pm1p
    DDPm1 = DPm1p - c*DPm1
    Pp1 = Δ/B*DDPm1
    Δp = 2r-2
    Pm1ppp = sm1_third(rec, r, Pm1, Pm1p, Pm1pp)
    DDPm1p = (Pm1ppp - cpp*Pm1 - 2cp*Pm1p - c*Pm1pp) - (cp*DPm1 + c*DPm1p)
    Pp1p = Δp/B*DDPm1 + Δ/B*DDPm1p
    DdPp1 = Pp1p + c*Pp1
    Dm1Pp1 = (Pp1p - c*Pp1) + (-1)*(Δp/Δ)*Pp1
    s1lplp = 2*σ*(1/Δ)*Dm1Pp1
    s1mpmp = 2*σ*(σ*DdPp1 + DPm1)

    h0  = valK(rec.h, r); h0p = derK(rec.h, r)
    Dh  = h0p - c*h0
    h0pp = sh0_second(rec, r, h0, h0p)
    Dhp  = h0pp - cp*h0 - c*h0p
    rDh  = r*Dh; rDhp = Dh + r*Dhp
    DrDh = rDhp - c*rDh
    κ0  = κvalK(rec.κ, r); κ0p = κderK(rec.κ, r)
    Dκ  = κ0p - c*κ0
    κ0pp = skappa_second(rec, r, κ0, κ0p)
    Dκp = κ0pp - cp*κ0 - c*κ0p
    DDκ = Dκp - c*Dκ
    s0lplp = -(1/(im*ω)*DrDh + 4*DDκ)
    s0mpmp_κ = -4*κ0
    s0mpmp_h = (a==0) ? zero(CT) : (a/ω)*h0

    (s2lplp=CT(s2lplp), s2mpmp=CT(s2mpmp), s1lplp=CT(s1lplp), s1mpmp=CT(s1mpmp),
     s0lplp=CT(s0lplp), s0mpmp_κ=CT(s0mpmp_κ), s0mpmp_h=CT(s0mpmp_h),
     κ0=CT(κ0), h0=CT(h0))
end

# P_{-1} third derivative from the s=-1 ODE (analytic).
function sm1_third(rec, r, P, Pp, Ppp)
    jd=rec.jd; ω=jd.ω; a=jd.a; m=jd.m; λ=jd.λ
    Δ=r^2-2r+a^2; Δp=2r-2; c,cp,cpp,_ = c_tower_kerr(r,a,ω,m)
    lhs_rest = Δp*(Ppp - cp*P - c^2*P)
    rhs = 2im*ω*(P + r*Pp) + λ*Pp
    Pppp = (rhs - lhs_rest)/Δ + (cpp*P + cp*Pp + 2c*cp*P + c^2*Pp)
    Pppp
end

# h second derivative from the s=0 ODE.
function sh0_second(rec, r, h, hp)
    jd=rec.jd; ω=jd.ω; a=jd.a; m=jd.m; λ0=jd.λ0
    Δ=r^2-2r+a^2; Δp=2r-2; K=ω*(r^2+a^2)-a*m; Kp=2ω*r
    (-(Δp*hp + (im*Kp + K^2/Δ)*h - 2im*ω*r*h - λ0*h))/Δ
end

# κ second derivative from the inhomogeneous s=0 ODE L[κ]=½(r²+a²γ_ll)h.
function skappa_second(rec, r, κ, κp)
    jd=rec.jd; ω=jd.ω; a=jd.a; m=jd.m; λ0=jd.λ0; γll=jd.γll
    Δ=r^2-2r+a^2; Δp=2r-2; K=ω*(r^2+a^2)-a*m; Kp=2ω*r
    S = 0.5*(r^2 + a^2*γll)*valK(rec.h, r)
    (S - (Δp*κp + (im*Kp + K^2/Δ)*κ - 2im*ω*r*κ - λ0*κ))/Δ
end

# ---- S-matrix scalar entries (extracted from closure_kerr; identical formulas) ----
"""
    smatrices_kerr(a, r0, m; lmax) -> NamedTuple of (l,lp)→Complex functions

The closure_kerr S-matrix scalars as standalone (l,lp)→Complex functions, for assembling the projected
components h^{l'}_{lplp}, h^{l'}_{mpmp}. b-matrices via SWSH expansion_coefficients (Float64 angular cap).
"""
function smatrices_kerr(a, r0, m; lmax=20)
    bg = orbit_background(a, r0, m)
    aw0 = Float64(bg.aω)
    lmin = max(2, abs(m))
    σof(l) = (-1)^(l+m)
    bget(s, l) = (aw0==0.0 ? Dict(l=>1.0) :
                  expansion_coefficients(SpinWeightedSpheroidalHarmonicS(s, l, m, aw0; method=:SphericalExpansion)))
    bval(bd, lp) = get(bd, lp, 0.0)
    function Smat_s2_lplp(l,lp)
        s=0.0+0im
        for (j,bj) in bget(2,l)
            t=(j==lp ? Λ̂(j) : 0.0); if aw0!=0.0; t+=-2aw0*λ̂2(j)*frak_s(1,0,j,lp,m)+aw0^2*frak_s2(2,0,j,lp,m); end
            s += bj*t
        end
        for (j,bj) in bget(-2,l)
            t=(j==lp ? Λ̂(j) : 0.0); if aw0!=0.0; t+=-2aw0*λ̂2(j)*frak_s(-1,0,j,lp,m)+aw0^2*frak_s2(-2,0,j,lp,m); end
            s += σof(l)*bj*t
        end
        s
    end
    Smat_s2_mpmp(l,lp) = bval(bget(2,l),lp)
    function Smat_s1_lplp(l,lp)
        s=0.0+0im
        for (j,bj) in bget(-1,l); t=(j==lp ? -λ̂(j) : 0.0); if aw0!=0.0; t+=aw0*frak_s(-1,0,j,lp,m); end; s+=bj*t; end
        for (j,bj) in bget(1,l);  t=(j==lp ? -λ̂(j) : 0.0); if aw0!=0.0; t+=aw0*frak_s(1,0,j,lp,m); end; s+=σof(l)*bj*t; end
        s
    end
    function Smat_s1_mpmp(l,lp)
        s=0.0+0im
        for (j,bj) in bget(1,l); t=(j==lp ? -λ̂2(j) : 0.0); if aw0!=0.0; t+=aw0*frak_s(1,2,j,lp,m); end; s+=bj*t; end
        s
    end
    Smat_s0_lplp(l,lp) = bval(bget(0,l),lp)
    function Smat_s0_mpmp_κ(l,lp)
        s=0.0+0im
        for (j,bj) in bget(0,l)
            t=(j==lp ? Λ̂(j) : 0.0); if aw0!=0.0; t+=-2aw0*λ̂(j)*frak_s(1,2,j,lp,m)+aw0^2*frak_s2(0,2,j,lp,m); end
            s += bj*t
        end
        s
    end
    function Smat_s0_mpmp_h(l,lp)
        aw0==0.0 && return 0.0+0im
        s=0.0+0im
        for (j,bj) in bget(0,l)
            t0=0.0+0im
            for k in (j-2):(j+2)
                k<lmin && continue
                A=(j==k ? Λ̂(j) : 0.0) - 2aw0*λ̂(j)*frak_s(1,2,j,k,m) + aw0^2*frak_s2(0,2,j,k,m)
                A==0.0 && continue
                t0 += A*frak_c(2,k,lp,m)
            end
            t1 = λ̂(j)*frak_s(1,2,j,lp,m) - aw0*frak_s2(0,2,j,lp,m)
            s += bj*(t0+t1)
        end
        s
    end
    (s2lplp=Smat_s2_lplp, s2mpmp=Smat_s2_mpmp, s1lplp=Smat_s1_lplp, s1mpmp=Smat_s1_mpmp,
     s0lplp=Smat_s0_lplp, s0mpmp_κ=Smat_s0_mpmp_κ, s0mpmp_h=Smat_s0_mpmp_h, lmin=lmin)
end

"""
    tetrad_lplp_kerr(a, r0, lp, m, r; recs, S, ls) -> Complex

The spherical-l' projected h_{l+l+}(r) tetrad component at a≠0:  Σ_l S^{(s)}_{lplp}(l,l')·R^{(s)}_{lplp}\[l\](r).
"""
function tetrad_lplp_kerr(a, r0, lp, m, r; recs=nothing, S=nothing, ls=nothing)
    recs, S, ls = _ensure_components_inputs(a, r0, m, lp, recs, S, ls)
    acc = zero(ComplexF64) + 0im
    for l in ls
        chain = radial_chain_kerr(recs[l], r)
        acc += S.s2lplp(l,lp)*chain.s2lplp + S.s1lplp(l,lp)*chain.s1lplp + S.s0lplp(l,lp)*chain.s0lplp
    end
    acc
end

"""
    tetrad_mpmp_kerr(a, r0, lp, m, r; recs, S, ls) -> Complex

The spherical-l' projected h_{m+m+}(r) tetrad component at a≠0.
"""
function tetrad_mpmp_kerr(a, r0, lp, m, r; recs=nothing, S=nothing, ls=nothing)
    recs, S, ls = _ensure_components_inputs(a, r0, m, lp, recs, S, ls)
    acc = zero(ComplexF64) + 0im
    for l in ls
        chain = radial_chain_kerr(recs[l], r)
        acc += S.s2mpmp(l,lp)*chain.s2mpmp + S.s1mpmp(l,lp)*chain.s1mpmp +
               S.s0mpmp_κ(l,lp)*chain.s0mpmp_κ + S.s0mpmp_h(l,lp)*chain.s0mpmp_h
    end
    acc
end

function _ensure_components_inputs(a, r0, m, lp, recs, S, ls)
    if recs === nothing || S === nothing || ls === nothing
        lmx = max(2*lp, lp+8)
        sol, aux = closure_kerr(a, r0, m; lmax=lmx, return_aux=true)
        lsv = aux.ls
        recs2 = Dict{Int,Any}()
        for l in lsv
            rec = reconstruct_kerr(a, r0, l; m=m, aux=aux, sol=sol)
            recs2[l] = _attach_recmeta(rec, l, aux)
        end
        S2 = smatrices_kerr(a, r0, m; lmax=lmx)
        return recs2, S2, lsv
    end
    (recs, S, ls)
end

# =============================================================================================================
# ANALYTIC SPIN-GENERAL DERIVATIVE TOWER.  Each homogeneous field satisfies Δ P''=A1 P'+A0 P (spin-general
# K=ω(r²+a²)−am, Δ=r²−2r+a², c=iK/Δ); P''',P'''',… by EXACT recurrence (`_field_tower5K`) seeded by the
# radial_kerr exact (valK,derK).  κ inhomogeneous (`_kappa_tower5K`, source S=½(r²+a²γ_ll)h).  Reuses the
# SchwarzschildComponents length-4 jet algebra (j4mul/j4add/j4scale/j4inv) + _tower5.  FD-free, BigFloat-clean.
# =============================================================================================================
_cK(r,a,ω,m) = c_tower_kerr(r, a, ω, m)

function _coef4K(r, a, ω, m, ::Type{CT}) where {CT}
    T = real(CT); rT=T(r); aT=T(a); ωT=T(ω); mT=T(m)
    Δ  = (CT(rT^2-2rT+aT^2), CT(2rT-2), CT(2), CT(0))
    K  = (CT(ωT*(rT^2+aT^2)-aT*mT), CT(2ωT*rT), CT(2ωT), CT(0))
    invΔ = j4inv(Δ)
    c  = j4scale(CT(im), j4mul(K, invΔ))
    Δp = (Δ[2], Δ[3], Δ[4], CT(0)); Kp = (K[2], K[3], K[4], CT(0))
    (Δ=Δ, Δp=Δp, K=K, Kp=Kp, invΔ=invΔ, c=c)
end

function _alphabeta4K(spin, r, a, ω, m, eig, ::Type{CT}) where {CT}
    co = _coef4K(r, a, ω, m, CT)
    Δ=co.Δ; Δp=co.Δp; c=co.c; invΔ=co.invΔ; K=co.K; Kp=co.Kp
    cp = (c[2], c[3], c[4], CT(0)); csq = j4mul(c,c)
    rj = (CT(r), CT(1), CT(0), CT(0)); eigj = (CT(eig), CT(0), CT(0), CT(0))
    Z4 = (CT(0),CT(0),CT(0),CT(0))
    if spin == -2
        A1 = Δp
        A0 = j4add(j4add(j4mul(Δ,cp), j4mul(Δ,csq)),
                   j4add(j4scale(CT(-1), j4mul(Δp,c)), j4add(j4scale(CT(6im*ω), rj), eigj)))
    elseif spin == -1
        A1 = Z4
        A0 = j4add(j4add(j4mul(Δ,cp), j4mul(Δ,csq)), j4add(j4scale(CT(2im*ω), rj), eigj))
    elseif spin == 0
        A1 = j4scale(CT(-1), Δp)
        KKoverΔ = j4mul(j4mul(K,K), invΔ)
        A0 = j4scale(CT(-1), j4add(j4add(j4scale(CT(im), Kp), KKoverΔ),
                                   j4add(j4scale(CT(-2im*ω), rj), j4scale(CT(-1),eigj))))
    elseif spin == 2
        A1 = Δp
        A0 = j4add(j4add(j4scale(CT(-1),j4mul(Δ,cp)), j4mul(Δ,csq)),
                   j4add(j4mul(Δp,c), j4add(j4scale(CT(-6im*ω), rj), eigj)))
    else
        error("_alphabeta4K: spin=$spin unsupported")
    end
    α = j4mul(A1, invΔ); β = j4mul(A0, invΔ)
    (α, β, invΔ)
end

function _field_tower5K(f, spin, a, ω, m, eig, r)
    P0 = valK(f,r); P1 = derK(f,r); CT = typeof(complex(P0))
    α, β, _ = _alphabeta4K(spin, r, a, ω, m, eig, CT)
    _tower5(CT(P0), CT(P1), α, β)
end
function _kappa_tower5K(kf, hT5, a, ω, m, λ0, γll, r)
    κ0v = κvalK(kf,r); κ1v = κderK(kf,r); CT = typeof(complex(κ0v))
    κ0 = CT(κ0v); κ1 = CT(κ1v)
    α, β, srcc = _alphabeta4K(0, r, a, ω, m, λ0, CT)
    T = real(CT); rT=T(r); aT=T(a); γT=T(γll)
    h0,h1,h2,h3,_ = hT5
    half = ((rT^2+aT^2*γT)/2)
    src_jet = (CT(half), CT(rT), CT(1), CT(0))
    hjet = (CT(h0),CT(h1),CT(h2),CT(h3))
    Sjet = j4mul(src_jet, hjet)
    g = j4mul(srcc, Sjet); g0,g1,g2,_ = g
    κ2 = α[1]*κ1 + β[1]*κ0 + g0
    κ3 = α[2]*κ1 + α[1]*κ2 + β[2]*κ0 + β[1]*κ1 + g1
    κ4 = α[3]*κ1 + 2α[2]*κ2 + α[1]*κ3 + β[3]*κ0 + 2β[2]*κ1 + β[1]*κ2 + g2
    (κ0,κ1,κ2,κ3,κ4)
end
function _D_applyK(t, c4, op)
    f0,f1,f2,f3,f4 = t; c,cp,cpp,cppp = c4
    s = op==:D ? -one(c) : +one(c)
    (f1 + s*c*f0,
     f2 + s*(cp*f0 + c*f1),
     f3 + s*(cpp*f0 + 2cp*f1 + c*f2),
     f4 + s*(cppp*f0 + 3cpp*f1 + 3cp*f2 + c*f3),
     zero(c))
end
function _chain_towerK(baseT, a, ω, m, ops)
    r->begin
        t = baseT(r); c4 = _cK(r,a,ω,m)
        for op in reverse(ops); t = _D_applyK(t, c4, op); end
        t
    end
end
_chainvalK(baseT, a, ω, m, ops, r) = first(_chain_towerK(baseT, a, ω, m, ops)(r))

# (Δ²R)' and (Δ²R)'' for a branch leg (:up/:inn), using R,R',R'' from the branch (analytic).
function _d_Δ2R(br, a, r, leg)
    Δ=r^2-2r+a^2; Δp=2r-2
    R  = leg===:up ? br.up(r)  : br.inn(r)
    Rp = leg===:up ? br.upd(r) : br.ind(r)
    2Δ*Δp*R + Δ^2*Rp
end
function _dd_Δ2R(br, a, r, leg)
    Δ=r^2-2r+a^2; Δp=2r-2; Δpp=2.0
    R   = leg===:up ? br.up(r)   : br.inn(r)
    Rp  = leg===:up ? br.upd(r)  : br.ind(r)
    Rpp = leg===:up ? br.updd(r) : br.indd(r)
    (2Δp^2 + 2Δ*Δpp)*R + 4Δ*Δp*Rp + Δ^2*Rpp
end

merge_namedtuple(nt, extra) = (; nt..., extra...)

# Attach P_{+2} (s=+2 reconstructed, with its own +2 jumps applied to P_{+2}=Δ²R_{+2}) + the (l,B) metadata the
# component chains consume. The BigFloat-Wronskian `big_at` path mirrors radial_kerr's jump solve for the
# high-magnitude s=+2 leg; absent ⇒ Float64 path byte-unchanged.
function _attach_recmeta(rec, l, aux)
    bg = aux.bg; ω = bg.ω; a = bg.a; m = bg.m
    Λ = aux.Λ[l]
    J0p2, J1p2 = jumps_pm(bg, l, +1)
    rspan = hasproperty(rec, :rspan) ? rec.rspan : nothing
    brp2 = branchK(2, l, m, a, ω, nothing; rspan=rspan)
    T = typeof(bg.r0)
    Δ0 = bg.Δ0
    _big_at_Δ2 = (hasproperty(brp2, :big_at) && brp2.big_at !== nothing) ? function (r)
        Ru, Rud, Ri, Rid = brp2.big_at(r)
        Δ  = BigFloat(Float64(r))^2 - 2*BigFloat(Float64(r)) + BigFloat(Float64(a))^2
        Δp = 2*BigFloat(Float64(r)) - 2
        Δ2 = Δ^2; dΔ2 = 2*Δ*Δp
        (Δ2*Ru, dΔ2*Ru + Δ2*Rud, Δ2*Ri, dΔ2*Ri + Δ2*Rid)
    end : nothing
    brp2_scaled = (up=r->brp2.up(r)*(r^2-2r+a^2)^2, upd=r->_d_Δ2R(brp2,a,r,:up),
                   inn=r->brp2.inn(r)*(r^2-2r+a^2)^2, ind=r->_d_Δ2R(brp2,a,r,:inn),
                   updd=r->_dd_Δ2R(brp2,a,r,:up), indd=r->_dd_Δ2R(brp2,a,r,:inn),
                   big_at=_big_at_Δ2)
    Pp2 = HomFieldK{T}(bg.r0, jump_ampsK(brp2_scaled, bg.r0, Complex{T}(J0p2), Complex{T}(J1p2))..., brp2_scaled)
    merge_namedtuple(rec, (Pp2=Pp2, _l=l, _B=aux.BB[l]))
end

# =============================================================================================================
# KerrComponents — the a≠0 (spin-general) TEN-component tetrad metric-perturbation surface (legacy
# KerrLorenzCirc/components_kerr.jl, the mixed/ΣΔ banded assembly + the full 10-mode projection). This is the
# headline AccretionDisc/ModifiedGravity consumer: `tetrad_modes_kerr` (the 10 tetrad modes lplp,lmlm,mpmp,mmmm,
# ρlpmp,ρclpmm,ρclmmp,ρlmmm,ΣΔlplm,h projected onto spherical-l') and `to_BLS_kerr` (the Barack-Lousto-Sago
# variables h̄¹..h̄¹⁰), plus the standalone mixed surfaces `tetrad_ρlpmp_kerr` / `tetrad_ΣΔlplm_kerr`.
#
# At a≠0 the spheroidal→spherical projection is BANDED: each tetrad mode is a MATRIX product of l-indexed radial
# derivative VECTORS against banded angular matrices (b-matrices × CG mixing 𝔰ⁿ/𝔠ⁿ + ladder diagonals) over the
# l-index space:  h^{l'}_{XX}(r) = Σ_l Σ_s S^{(s)}_{XX}(l,l') · R^{(s)}_{XX}[l](r). The radial vectors reuse the
# FD-free spin-general analytic jet towers built in KerrSMatrices (_mixed_towers/_chainvalK/_field_tower5K); the
# angular matrices reuse the Operators frak_s/frak_s2/frak_c/frak_c2 CG-mixing tables (promoted onto Wigner cg).
#
# At a→0 every aω-ladder term, every ia·cos/sin piece, every a²·Cos piece VANISHES, the b-matrices collapse to I,
# ρmat→r·I, Σmat→r²·I — so the contraction reduces EXACTLY to the validated a=0 SchwarzschildComponents tetrad
# modes (gate G-K11). The diagonal lplp/mpmp pieces reuse the SAME App-H-validated banded S-matrices as
# tetrad_lplp_kerr/tetrad_mpmp_kerr (App-H a=0.99 lower 6.85e-8); the lmlm/mmmm/h and mixed/ΣΔ are the same
# b-matrix×CG structural generalization with the a=0 diagonals as the validated reduction target.
#
# Element type is INFERRED from the reconstructed radial field (CT = typeof(complex(valK(...)))), so BigFloat
# precision threads through the banded contraction (at a=0 the projection is exactly b=I ⇒ pure-CT arithmetic).
# The angular b-matrices carry the documented Float64 SWSH spheroidal cap (aω-dependent), promoted into CT.
# The r-INDEPENDENT angular matrices and the per-r radial-vector bundles are memoized under a ReentrantLock
# (blas_guard cache pattern); results are byte-identical to per-call construction (synchronization only).
# =============================================================================================================

using LinearAlgebra: I

export tetrad_modes_kerr, to_BLS_kerr, tetrad_ρlpmp_kerr, tetrad_ΣΔlplm_kerr

# -------------------------------------------------------------------------------------------------------------
# Analytic spin-general tower bundle for one spheroidal-l (the a≠0 analogue of build_recon's towers). FD-free:
# every base tower is the exact ODE-recurrence 5-tuple seeded by the radial_kerr exact (valK,derK). Reuses the
# KerrSMatrices jet machinery (_field_tower5K/_kappa_tower5K/_chain_towerK).
# -------------------------------------------------------------------------------------------------------------
function _mixed_towers(rec, a, ω, m)
    jd = rec.jd; B = jd_B(rec); σ = (-1)^(rec_l(rec)+m)
    Λ = jd.Λ; λ = jd.λ; λ0 = jd.λ0; γll = jd.γll
    Pm2T = r->_field_tower5K(rec.Pm2, -2, a, ω, m, Λ,  r)
    Pp2T = r->_field_tower5K(rec.Pp2,  2, a, ω, m, Λ,  r)   # P_{+2}=Δ²R_{+2}: s=+2 ODE (eig=Λ)
    Pm1T = r->_field_tower5K(rec.Pm1, -1, a, ω, m, λ,  r)
    hT   = r->_field_tower5K(rec.h,    0, a, ω, m, λ0, r)
    κT   = r->_kappa_tower5K(rec.κ, _field_tower5K(rec.h,0,a,ω,m,λ0,r), a, ω, m, λ0, γll, r)
    D2Pm1_tower = _chain_towerK(Pm1T, a, ω, m, [:D,:D])
    Pp1_5 = r->begin
        Δr=r^2-2r+a^2; Δrp=2r-2; Δrpp=2*one(real(complex(Δr))); g=D2Pm1_tower(r)
        g0=g[1]; g1=g[2]; g2=g[3]
        (Δr/B*g0, (Δrp*g0+Δr*g1)/B, (Δrpp*g0+2Δrp*g1+Δr*g2)/B, zero(g0), zero(g0))
    end
    (Pm2T=Pm2T, Pp2T=Pp2T, Pm1T=Pm1T, hT=hT, κT=κT, Pp1_5=Pp1_5, B=B, σ=σ)
end

# small Δ-scaling tower helper used by the deep ΣΔ chains: (Δ·X, Δ'·X + Δ·X', 0, …)
@inline function _scaleΔ_tower(g, r, a)
    Δr=r^2-2r+a^2; Δrp=2r-2; g0=g[1]; g1=g[2]
    (Δr*g0, Δrp*g0+Δr*g1, zero(g0), zero(g0), zero(g0))
end

# Per-spheroidal-l radial chain scalars for the MIXED (ρ h_{l±m±}) + ΣΔ banded contraction (spin-general).
function _mixed_radial_vecs(rec, r, m)
    jd = rec.jd; ω = jd.ω; a = jd.a
    tw = _mixed_towers(rec, a, ω, m)
    σ = tw.σ
    Pm2T=tw.Pm2T; Pp2T=tw.Pp2T; Pm1T=tw.Pm1T; hT=tw.hT; κT=tw.κT; Pp1_5=tw.Pp1_5

    DPm2     = _chainvalK(Pm2T, a,ω,m, [:D], r)
    DDPm2    = _chainvalK(Pm2T, a,ω,m, [:D,:D], r)
    DDDPm2   = _chainvalK(Pm2T, a,ω,m, [:D,:D,:D], r)
    DdDDPm2  = _chainvalK(Pm2T, a,ω,m, [:Dd,:D,:D], r)
    DdPp2     = _chainvalK(Pp2T, a,ω,m, [:Dd], r)
    DdDdPp2   = _chainvalK(Pp2T, a,ω,m, [:Dd,:Dd], r)
    DdDdDdPp2 = _chainvalK(Pp2T, a,ω,m, [:Dd,:Dd,:Dd], r)
    DDdDdPp2  = _chainvalK(Pp2T, a,ω,m, [:D,:Dd,:Dd], r)

    DDPm2_tower = _chain_towerK(Pm2T, a,ω,m, [:D,:D])
    DdDDPm2_tower = _chain_towerK(Pm2T, a,ω,m, [:Dd,:D,:D])
    ΔDdDDPm2_tower = r2->_scaleΔ_tower(DdDDPm2_tower(r2), r2, a)
    DΔDdDDPm2 = _chainvalK(ΔDdDDPm2_tower, a,ω,m, [:D], r)
    ΔDDDPm2_tower = r2->_scaleΔ_tower(_chain_towerK(Pm2T,a,ω,m,[:D,:D,:D])(r2), r2, a)
    DdΔDDDPm2 = _chainvalK(ΔDDDPm2_tower, a,ω,m, [:Dd], r)
    LDD_Pm2 = DΔDdDDPm2 + DdΔDDDPm2
    drDDPm2 = DDPm2_tower(r)[2]
    DdDdPp2_tower = _chain_towerK(Pp2T, a,ω,m, [:Dd,:Dd])
    DDdDdPp2_tower = _chain_towerK(Pp2T, a,ω,m, [:D,:Dd,:Dd])
    ΔDDdDdPp2_tower = r2->_scaleΔ_tower(DDdDdPp2_tower(r2), r2, a)
    DdΔDDdDdPp2 = _chainvalK(ΔDDdDdPp2_tower, a,ω,m, [:Dd], r)
    ΔDdDdDdPp2_tower = r2->_scaleΔ_tower(_chain_towerK(Pp2T,a,ω,m,[:Dd,:Dd,:Dd])(r2), r2, a)
    DΔDdDdDdPp2 = _chainvalK(ΔDdDdDdPp2_tower, a,ω,m, [:D], r)
    LDD_Pp2 = DΔDdDdDdPp2 + DdΔDDdDdPp2
    drDdDdPp2 = DdDdPp2_tower(r)[2]

    Pm1_0  = Pm1T(r)[1]
    Pp1    = Pp1_5(r)[1]
    DPm1_tower  = _chain_towerK(Pm1T, a,ω,m, [:D])
    DdPp1_tower = _chain_towerK(Pp1_5, a,ω,m, [:Dd])
    calPs1_tower = r2->begin
        aa=DPm1_tower(r2); bb=DdPp1_tower(r2)
        (aa[1]+σ*bb[1], aa[2]+σ*bb[2], aa[3]+σ*bb[3], zero(aa[1]), zero(aa[1]))
    end
    calPs1  = first(calPs1_tower(r))
    DcalPs1  = _chainvalK(calPs1_tower, a,ω,m, [:D], r)
    DdcalPs1 = _chainvalK(calPs1_tower, a,ω,m, [:Dd], r)

    h0 = hT(r)[1]; Dh = _chainvalK(hT,a,ω,m,[:D],r); Ddh = _chainvalK(hT,a,ω,m,[:Dd],r)
    κ0 = κT(r)[1]; Dκ = _chainvalK(κT,a,ω,m,[:D],r); Ddκ = _chainvalK(κT,a,ω,m,[:Dd],r)
    κp = κT(r)[2]
    Ddκ_tower = _chain_towerK(κT, a,ω,m, [:Dd])
    ΔDdκ_tower = r2->_scaleΔ_tower(Ddκ_tower(r2), r2, a)
    DΔDdκ = _chainvalK(ΔDdκ_tower, a,ω,m, [:D], r)
    Dκ_tower = _chain_towerK(κT, a,ω,m, [:D])
    ΔDκ_tower = r2->_scaleΔ_tower(Dκ_tower(r2), r2, a)
    DdΔDκ = _chainvalK(ΔDκ_tower, a,ω,m, [:Dd], r)

    (σ=σ, DPm2=DPm2, DDPm2=DDPm2, DDDPm2=DDDPm2, DdDDPm2=DdDDPm2, DdPp2=DdPp2, DdDdPp2=DdDdPp2,
     DdDdDdPp2=DdDdDdPp2, DDdDdPp2=DDdDdPp2,
     LDD_Pm2=LDD_Pm2, drDDPm2=drDDPm2, LDD_Pp2=LDD_Pp2, drDdDdPp2=drDdDdPp2,
     calPs1=calPs1, DcalPs1=DcalPs1, DdcalPs1=DdcalPs1, Pm1_0=Pm1_0, Pp1=Pp1,
     h0=h0, Dh=Dh, Ddh=Ddh, κ0=κ0, Dκ=Dκ, Ddκ=Ddκ, κp=κp, DΔDdκ=DΔDdκ, DdΔDκ=DdΔDκ)
end

# Per-spheroidal-l radial chain scalars for the DIAGONAL sector (lplp,lmlm,mpmp,mmmm,h).
function _diag_radial_vecs(rec, r, m)
    jd = rec.jd; ω = jd.ω; a = jd.a
    tw = _mixed_towers(rec, a, ω, m)
    Δ = r^2-2r+a^2; Δp = 2r-2
    Pm2T=tw.Pm2T; Pp2T=tw.Pp2T; Pm1T=tw.Pm1T; hT=tw.hT; κT=tw.κT; Pp1_5=tw.Pp1_5
    Pm2_0 = Pm2T(r)[1]; Pp2_0 = Pp2T(r)[1]
    DDPm2   = _chainvalK(Pm2T, a,ω,m, [:D,:D], r)
    DdDdPp2 = _chainvalK(Pp2T, a,ω,m, [:Dd,:Dd], r)
    DPm1   = _chainvalK(Pm1T, a,ω,m, [:D], r)
    DdPp1  = _chainvalK(Pp1_5, a,ω,m, [:Dd], r)
    Pp1    = Pp1_5(r)[1]; Pm1_0 = Pm1T(r)[1]
    Dm1Pp1 = _chainvalK(Pp1_5, a,ω,m, [:D], r) + (-1)*(Δp/Δ)*Pp1
    Ddm1Pm1 = _chainvalK(Pm1T, a,ω,m, [:Dd], r) + (-1)*(Δp/Δ)*Pm1_0
    h0 = hT(r)[1]; κ0 = κT(r)[1]
    DDκ   = _chainvalK(κT, a,ω,m, [:D,:D], r)
    DdDdκ = _chainvalK(κT, a,ω,m, [:Dd,:Dd], r)
    Dh_tower = _chain_towerK(hT, a,ω,m, [:D])
    rDh_tower = r2->(g=Dh_tower(r2); g0=g[1]; g1=g[2]; (r2*g0, g0+r2*g1, zero(g0),zero(g0),zero(g0)))
    DrDh = _chainvalK(rDh_tower, a,ω,m, [:D], r)
    Ddh_tower = _chain_towerK(hT, a,ω,m, [:Dd])
    rDdh_tower = r2->(g=Ddh_tower(r2); g0=g[1]; g1=g[2]; (r2*g0, g0+r2*g1, zero(g0),zero(g0),zero(g0)))
    DdrDdh = _chainvalK(rDdh_tower, a,ω,m, [:Dd], r)
    (σ=tw.σ, Pm2_0=Pm2_0, Pp2_0=Pp2_0, DDPm2=DDPm2, DdDdPp2=DdDdPp2, DPm1=DPm1, DdPp1=DdPp1,
     Dm1Pp1=Dm1Pp1, Ddm1Pm1=Ddm1Pm1, h0=h0, κ0=κ0, DDκ=DDκ, DdDdκ=DdDdκ, DrDh=DrDh, DdrDdh=DdrDdh)
end

# -------------------------------------------------------------------------------------------------------------
# Banded angular matrices over the l-list `ls` (dense, positions 1..n). r-INDEPENDENT (depend only on ls, m, aω,
# elt type). B-matrices via SWSH expansion_coefficients (Float64 spheroidal cap, promoted into CT); SinMix/CosMix
# via the Operators frak_s/frak_s2/frak_c/frak_c2 CG-mixing tables. Memoized once per (ls,m,aω,CT) under a lock;
# the matrices are immutable once built ⇒ all later calls are cache-hit reads (byte-identical, sync only).
# -------------------------------------------------------------------------------------------------------------
const _ANGMAT_CACHE = Dict{Any,Any}()
const _ANG_LOCK = ReentrantLock()
_angular_mats_kerr(ls, m, aw0, CT) = lock(_ANG_LOCK) do
    get!(_ANGMAT_CACHE, (Tuple(ls), Int(m), aw0, CT)) do
        _build_angular_mats_kerr(ls, m, aw0, CT)
    end
end

function _build_angular_mats_kerr(ls, m, aw0, CT)
    n = length(ls)
    bget(s,l) = (aw0==0.0 ? Dict(l=>1.0) :
                 expansion_coefficients(SpinWeightedSpheroidalHarmonicS(s, l, m, aw0; method=:SphericalExpansion)))
    Bcache=Dict{Int,Matrix{CT}}(); Scache=Dict{Tuple{Int,Int},Matrix{CT}}()
    C1cache=Dict{Int,Matrix{CT}}(); C2cache=Dict{Int,Matrix{CT}}()
    clk = ReentrantLock()
    Bmat(s) = lock(clk) do; get!(Bcache, s) do
        M=zeros(CT,n,n); for (i,l) in enumerate(ls); for (lp,b) in bget(s,l); j=findfirst(==(lp),ls); j===nothing && continue; M[i,j]=b; end; end; M
    end; end
    _sinel(s1,s2,l,lp) = ((s1,s2) in ((2,0),(-2,0),(0,2),(0,-2)) ? frak_s2(s1,s2,l,lp,m) : frak_s(s1,s2,l,lp,m))
    SinM(s1,s2) = lock(clk) do; get!(Scache, (s1,s2)) do
        M=zeros(CT,n,n); for i in 1:n, j in 1:n; abs(ls[i]-ls[j])<=2 && (M[i,j]=_sinel(s1,s2,ls[i],ls[j])); end; M
    end; end
    Cos1(s) = lock(clk) do; get!(C1cache, s) do
        M=zeros(CT,n,n); for i in 1:n, j in 1:n; abs(ls[i]-ls[j])<=1 && (M[i,j]=frak_c(s,ls[i],ls[j],m)); end; M
    end; end
    Cos2(s) = lock(clk) do; get!(C2cache, s) do
        M=zeros(CT,n,n); for i in 1:n, j in 1:n; abs(ls[i]-ls[j])<=2 && (M[i,j]=frak_c2(s,ls[i],ls[j],m)); end; M
    end; end
    Diagm(f) = (M=zeros(CT,n,n); for i in 1:n; M[i,i]=CT(f(ls[i])); end; M)
    λ̂m = Diagm(λ̂); λ̂2m = Diagm(λ̂2); Λ̂m = Diagm(Λ̂)
    sign = Diagm(l->CT((-1)^(l+m)))
    idm  = Matrix{CT}(I, n, n)
    (n=n, B=Bmat, Sin=SinM, Cos1=Cos1, Cos2=Cos2, λ̂=λ̂m, λ̂2=λ̂2m, Λ̂=Λ̂m, sign=sign, id=idm)
end

# Infer the complex element type from the reconstructed radial field at r (so BigFloat threads through).
_ct_at(recs, ls, r) = typeof(complex(valK(recs[ls[1]].Pm2, r)))

# -------------------------------------------------------------------------------------------------------------
# Banded MIXED (ρ h_{l±m±}) + ΣΔ (h_{l+l-}) contraction for ALL spherical-l' at once. The memo key
# includes the complete contraction signature: different l-bands or physical parameters over the same
# reconstruction/radius are distinct computations, and arbitrary-precision radii must not be Float64-folded.
# -------------------------------------------------------------------------------------------------------------
const _MIXBAND_CACHE = IdDict{Any,Dict{Any,Any}}()
const _DIAGBAND_CACHE = IdDict{Any,Dict{Any,Any}}()
const _BAND_LOCK = ReentrantLock()
const _BAND_CACHE_LIMIT = 256
const _BAND_RECON_CACHE_LIMIT = 16

# Reconstruction bundles are immutable-after-build by contract. Replacing a recs[l]
# entry changes the signature below; callers that mutate an existing record in place
# must explicitly invalidate the bundle.
function _invalidate_band_cache!(recs)
    lock(_BAND_LOCK) do
        delete!(_MIXBAND_CACHE, recs)
        delete!(_DIAGBAND_CACHE, recs)
    end
    return nothing
end

_band_scalar_key(x) = (typeof(x), x, x isa BigFloat ? precision(x) : 0)
_band_recs_signature(recs, ls) = Tuple(objectid(recs[l]) for l in ls)
_band_cache_key(recs, ls, r, m, a, ω) =
    (Tuple(ls), _band_recs_signature(recs, ls), Int(m), _band_scalar_key(a),
     _band_scalar_key(ω), _band_scalar_key(r), precision(BigFloat), rounding(BigFloat))

function _band_cache_get(cache, recs, key)
    lock(_BAND_LOCK) do
        inner = get(cache, recs, nothing)
        inner === nothing ? nothing : get(inner, key, nothing)
    end
end

function _band_cache_put!(cache, recs, key, value)
    lock(_BAND_LOCK) do
        length(cache) >= _BAND_RECON_CACHE_LIMIT && !haskey(cache, recs) && empty!(cache)
        inner = get!(cache, recs) do
            Dict{Any,Any}()
        end
        length(inner) >= _BAND_CACHE_LIMIT && !haskey(inner, key) && empty!(inner)
        inner[key] = value
    end
    return value
end

function _mixed_banded_all(recs, ls, r, m, a, ω)
    key = _band_cache_key(recs, ls, r, m, a, ω)
    c = _band_cache_get(_MIXBAND_CACHE, recs, key)
    c !== nothing && return c
    v = _mixed_banded_all_impl(recs, ls, r, m, a, ω)
    _band_cache_put!(_MIXBAND_CACHE, recs, key, v)
end

function _mixed_banded_all_impl(recs, ls, r, m, a, ω)
    CT = _ct_at(recs, ls, r)
    aw0 = Float64(a*ω)
    A = _angular_mats_kerr(ls, m, aw0, CT)
    n = A.n
    β = 1/(-6im*ω); Δ = r^2-2r+a^2
    V = [_mixed_radial_vecs(recs[ls[i]], r, m) for i in 1:n]
    σv = CT[V[i].σ for i in 1:n]
    vget(f) = CT[getfield(V[i], f) for i in 1:n]
    DPm2=vget(:DPm2); DDPm2=vget(:DDPm2); DDDPm2=vget(:DDDPm2); DdDDPm2=vget(:DdDDPm2)
    DdPp2=vget(:DdPp2); DdDdPp2=vget(:DdDdPp2); DdDdDdPp2=vget(:DdDdDdPp2); DDdDdPp2=vget(:DDdDdPp2)
    calPs1=vget(:calPs1); DcalPs1=vget(:DcalPs1); DdcalPs1=vget(:DdcalPs1); Pm1_0=vget(:Pm1_0); Pp1=vget(:Pp1)
    h0=vget(:h0); Dh=vget(:Dh); Ddh=vget(:Ddh); κ0=vget(:κ0); Dκ=vget(:Dκ); Ddκ=vget(:Ddκ); κp=vget(:κp)
    LDD_Pm2=vget(:LDD_Pm2); drDDPm2=vget(:drDDPm2); LDD_Pp2=vget(:LDD_Pp2); drDdDdPp2=vget(:drDdDdPp2)
    DΔDdκ=vget(:DΔDdκ); DdΔDκ=vget(:DdΔDκ)

    sgn = A.sign; λ̂M=A.λ̂; λ̂2M=A.λ̂2; Λ̂M=A.Λ̂; idm=A.id
    Bm2=A.B(-2); Bp2=A.B(2); Bm1=A.B(-1); Bp1=A.B(1); B0=A.B(0)
    SSplus  = Bp2*(Λ̂M - 2aw0*λ̂2M*A.Sin(1,0) + aw0^2*A.Sin(2,0))
    SSminus = Bm2*(Λ̂M - 2aw0*λ̂2M*A.Sin(-1,0) + aw0^2*A.Sin(-2,0))
    MMp2 = Bp2*(λ̂2M - aw0*A.Sin(2,1))
    MMm2 = Bm2*(λ̂2M - aw0*A.Sin(-2,-1))
    SS0  = -(Bm1 + sgn*Bp1)*λ̂M + aw0*(Bm1*A.Sin(-1,0) + sgn*Bp1*A.Sin(1,0))

    # ===================== s=2 mixed =====================
    sinmix01 = A.Sin(0,1); cosmix11 = A.Cos1(1); sinmix1m1 = A.Sin(-1,1)
    ρmat  = r*idm + im*a*cosmix11
    cosm11_m = A.Cos1(-1); sinmix0m1 = A.Sin(0,-1)
    ρcmat = r*idm - im*a*cosm11_m

    t1 = (1/(12ω^2))*( (transpose(DDPm2)*sgn*SSplus) - (transpose(DdDdPp2)*sgn*SSminus) ) * (λ̂M - aw0*sinmix01)
    t2 = -(1/(12ω^2))*( (transpose(DDDPm2)*sgn*SSplus) - (transpose(DDdDdPp2)*sgn*SSminus) ) * ((λ̂M - aw0*sinmix01)*ρmat - im*a*sinmix01)
    t3 = -(1/(3im*ω))*( (transpose(DDDPm2)*sgn*MMp2)*ρmat*ρmat - 2*(transpose(DDPm2)*sgn*MMp2)*ρmat + 2*(transpose(DPm2)*sgn*MMp2) )
    tmpA = (λ̂M*λ̂M - 2aw0*λ̂M*sinmix01 + aw0^2*sinmix1m1)*ρmat*ρmat - 2im*a*(λ̂M*sinmix01 - aw0*sinmix1m1)*ρmat - 2a^2*sinmix1m1
    t4 = -(1/(3im*ω*Δ))*( transpose(DdPp2)*(sgn*MMm2*tmpA) )
    ρhs2lpmp = β*(t1 + t2 + t3 + t4)

    sinmix = sinmix0m1; sinmix1m1b = A.Sin(1,-1)
    u1 = -(1/(12ω^2))*( (transpose(DDPm2)*SSminus) - (transpose(DdDdPp2)*SSplus) ) * (λ̂M - aw0*sinmix)
    u2 =  (1/(12ω^2))*( (transpose(DDDPm2)*SSminus) - (transpose(DDdDdPp2)*SSplus) ) * ((λ̂M - aw0*sinmix)*ρcmat - im*a*sinmix)
    u3 =  (1/(3im*ω))*( (transpose(DDDPm2)*MMm2)*ρcmat*ρcmat - 2*(transpose(DDPm2)*MMm2)*ρcmat + 2*(transpose(DPm2)*MMm2) )
    tmpB = (λ̂M*λ̂M - 2aw0*λ̂M*sinmix + aw0^2*sinmix1m1b)*ρcmat*ρcmat - 2im*a*(λ̂M*sinmix - aw0*sinmix1m1b)*ρcmat - 2a^2*sinmix1m1b
    u4 =  (1/(3im*ω*Δ))*( (transpose(DdPp2)*MMp2)*tmpB )
    ρchs2lpmm = β*(u1 + u2 + u3 + u4)

    ρcmat_v = r*idm - im*a*cosmix11
    v1 = (1/(12ω^2))*( (transpose(DDPm2)*SSminus) - (transpose(DdDdPp2)*SSplus) ) * (λ̂M - aw0*sinmix01)
    v2 = -(1/(12ω^2))*( (transpose(DdDDPm2)*SSminus) - (transpose(DdDdDdPp2)*SSplus) ) * ((λ̂M - aw0*sinmix01)*ρcmat_v + im*a*sinmix01)
    v3 = -(1/(3im*ω))*( (transpose(DdDdDdPp2)*MMp2)*ρcmat_v*ρcmat_v - 2*(transpose(DdDdPp2)*MMp2)*ρcmat_v + 2*(transpose(DdPp2)*MMp2) )
    tmpV = (λ̂M*λ̂M - 2aw0*λ̂M*sinmix01 + aw0^2*sinmix1m1)*ρcmat_v*ρcmat_v + 2im*a*(λ̂M*sinmix01 - aw0*sinmix1m1)*ρcmat_v - 2a^2*sinmix1m1
    v4 = -(1/(3im*ω*Δ))*( transpose(DPm2)*(MMm2*tmpV) )
    ρchs2lmmp = β*(v1 + v2 + v3 + v4)

    cosmix_m1_1 = cosm11_m
    ρmat_w = r*idm + im*a*cosmix_m1_1
    sinmix1m1w = A.Sin(1,-1)
    w1 = -(1/(12ω^2))*( (transpose(DDPm2)*sgn*SSplus) - (transpose(DdDdPp2)*sgn*SSminus) ) * (λ̂M - aw0*sinmix0m1)
    w2 =  (1/(12ω^2))*( (transpose(DdDDPm2)*sgn*SSplus) - (transpose(DdDdDdPp2)*sgn*SSminus) ) * ((λ̂M - aw0*sinmix0m1)*ρmat_w + im*a*sinmix0m1)
    w3 =  (1/(3im*ω))*( (transpose(DdDdDdPp2)*sgn*MMm2)*ρmat_w*ρmat_w - 2*(transpose(DdDdPp2)*sgn*MMm2)*ρmat_w + 2*(transpose(DdPp2)*sgn*MMm2) )
    tmpW = (λ̂M*λ̂M - 2aw0*λ̂M*sinmix0m1 + aw0^2*sinmix1m1w)*ρmat_w*ρmat_w + 2im*a*(λ̂M*sinmix0m1 - aw0*sinmix1m1w)*ρmat_w - 2a^2*sinmix1m1w
    w4 =  (1/(3im*ω*Δ))*( transpose(DPm2)*(sgn*MMp2*tmpW) )
    ρhs2lmmm = β*(w1 + w2 + w3 + w4)

    # ===================== s=1 mixed =====================
    cosmix = cosmix11; sinmix = sinmix01
    line1 = (transpose(sgn*(r*DcalPs1 - 2*calPs1))*Bp1) + im*a*(transpose(sgn*DcalPs1)*Bp1*cosmix)
    tmp1  = transpose((1/Δ)*(sgn*Pp1))
    tmp2  = -r*λ̂M - im*a*λ̂M*cosmix + (a*ω*r + 2im*a)*sinmix + im*a^2*ω*sinmix*cosmix
    line2 = tmp1*SS0*tmp2
    ρhs1lpmp = vec(line1) .+ vec(line2)
    line1 = -(transpose(r*DcalPs1 - 2*calPs1)*Bm1) + im*a*(transpose(DcalPs1)*Bm1*cosm11_m)
    tmp1  = transpose((1/Δ)*(sgn*Pp1))
    tmp2  = r*λ̂M - im*a*λ̂M*cosm11_m - (a*ω*r + 2im*a)*sinmix0m1 + im*a^2*ω*sinmix0m1*cosm11_m
    line2 = tmp1*SS0*tmp2
    ρchs1lpmm = vec(line1) .+ vec(line2)
    line1 = (transpose(sgn*(r*DdcalPs1 - 2*calPs1))*Bp1) - im*a*(transpose(sgn*DdcalPs1)*Bp1*cosmix)
    tmp1  = transpose((1/Δ)*Pm1_0)
    tmp2  = -r*λ̂M + im*a*λ̂M*cosmix + (a*ω*r - 2im*a)*sinmix - im*a^2*ω*sinmix*cosmix
    line2 = tmp1*SS0*tmp2
    ρchs1lmmp = vec(line1) .+ vec(line2)
    cosmix_m = cosm11_m
    line1 = -(transpose(r*DdcalPs1 - 2*calPs1)*Bm1) - im*a*(transpose(DdcalPs1)*Bm1*cosmix_m)
    tmp1  = transpose((1/Δ)*Pm1_0)
    tmp2  = r*λ̂M + im*a*λ̂M*cosmix_m - (a*ω*r - 2im*a)*sinmix0m1 - im*a^2*ω*sinmix0m1*cosmix_m
    line2 = tmp1*SS0*tmp2
    ρhs1lmmm = vec(line1) .+ vec(line2)

    # ===================== s=0 mixed =====================
    Dhvec = transpose(Dh); Ddhvec = transpose(Ddh)
    Dκvec = transpose(Dκ); Ddκvec = transpose(Ddκ)
    hvec  = transpose(h0); κvec = transpose(κ0)
    t0   = r^2*idm + a^2*A.Cos2(1)
    sinm = sinmix01; cosm = cosmix11; MM0 = λ̂M - a*ω*sinm
    p1 = (1/(2im*ω))*Dhvec*B0*MM0*t0
    p2 = -a/ω*( r*Dhvec*B0*sinm - hvec*B0*MM0*cosm )
    tmp = r*idm + im*a*cosm
    p3 = 4*Dκvec*B0*MM0*tmp
    p4 = -4*( κvec*B0*MM0 + im*a*Dκvec*B0*sinm )
    ρhs0lpmp = vec(p1) .+ vec(p2) .+ vec(p3) .+ vec(p4)
    t0   = r^2*idm + a^2*A.Cos2(-1)
    sinm = sinmix0m1; MM0 = λ̂M - a*ω*sinm
    q1 = -(1/(2im*ω))*Dhvec*B0*MM0*t0
    q2 = a/ω*( r*Dhvec*B0*sinm + hvec*B0*MM0*cosm11_m )
    tmp = r*idm - im*a*cosm11_m
    q3 = -4*( Dκvec*B0*MM0*tmp )
    q4 = 4*( κvec*B0*MM0 + im*a*Dκvec*B0*sinm )
    ρchs0lpmm = vec(q1) .+ vec(q2) .+ vec(q3) .+ vec(q4)
    t0   = r^2*idm + a^2*A.Cos2(1)
    sinm = sinmix01; cosm = cosmix11; MM0 = λ̂M - a*ω*sinm
    x1 = -(1/(2im*ω))*Ddhvec*B0*MM0*t0
    x2 = -a/ω*( r*Ddhvec*B0*sinm - hvec*B0*MM0*cosm )
    tmp = r*idm - im*a*cosm
    x3 = 4*( Ddκvec*B0*MM0*tmp )
    x4 = -4*( κvec*B0*MM0 - im*a*Ddκvec*B0*sinm )
    ρchs0lmmp = vec(x1) .+ vec(x2) .+ vec(x3) .+ vec(x4)
    t0   = r^2*idm + a^2*A.Cos2(-1)
    sinm = sinmix0m1; cosm = cosm11_m; MM0 = λ̂M - a*ω*sinm
    y1 = (1/(2im*ω))*Ddhvec*B0*MM0*t0
    y2 = a/ω*( r*Ddhvec*B0*sinm + hvec*B0*MM0*cosm )
    tmp = r*idm + im*a*cosm
    y3 = -4*( Ddκvec*B0*MM0*tmp )
    y4 = 4*( κvec*B0*MM0 - im*a*Ddκvec*B0*sinm )
    ρhs0lmmm = vec(y1) .+ vec(y2) .+ vec(y3) .+ vec(y4)

    # ===================== ΣΔ (l+ l-) =====================
    Σmat = r^2*idm + a^2*A.Cos2(0)
    cosmix01 = A.Cos1(0)
    ρmatΣ = r*idm + im*a*cosmix01
    ρcmatΣ = r*idm - im*a*cosmix01
    tmp3 = λ̂M*(A.Sin(-1,0) - A.Sin(1,0))
    s2t1 = (1/(24ω^2))*( transpose(LDD_Pm2)*SSminus*Σmat - 4*r*Δ*transpose(drDDPm2)*SSminus - 2a^2*transpose(DDPm2)*SSminus*tmp3*cosmix01 )
    s2t2 = -(1/(24ω^2))*( transpose(LDD_Pp2)*SSplus*Σmat - 4*r*Δ*transpose(drDdDdPp2)*SSplus - 2a^2*transpose(DdDdPp2)*SSplus*tmp3*cosmix01 )
    tmpΣ = 2a^2*r*cosmix01 - im*a*(r^2*idm + 3a^2*A.Cos2(0))
    s2t3 = (1/(3im*ω))*( transpose(DDPm2)*SSminus*Σmat*ρcmatΣ - transpose(DPm2)*SSminus*ρcmatΣ*ρcmatΣ
                         - transpose(DDPm2)*MMm2*A.Sin(-1,0)*tmpΣ - 2im*a*transpose(DPm2)*MMm2*A.Sin(-1,0)*ρmatΣ )
    s2t4 = (1/(3im*ω))*( transpose(DdDdPp2)*SSplus*Σmat*ρcmatΣ - transpose(DdPp2)*SSplus*ρcmatΣ*ρcmatΣ
                         + transpose(DdDdPp2)*MMp2*A.Sin(1,0)*tmpΣ + 2im*a*transpose(DdPp2)*MMp2*A.Sin(1,0)*ρmatΣ )
    ΣΔs2 = β*(s2t1 + s2t2 + s2t3 + s2t4)
    g1 = transpose(calPs1)*SS0*Σmat
    g2 = -2*r*transpose(Pm1_0 .+ σv.*Pp1)*SS0
    g3 = 2a^2*transpose(calPs1)*(Bm1*A.Sin(-1,0) - sgn*Bp1*A.Sin(1,0))*cosmix01
    ΣΔs1 = vec(g1) .+ vec(g2) .+ vec(g3)
    Kr = ω*(r^2+a^2) - a*m
    z1 = -hvec*B0*((Kr/ω)*idm - 2*Σmat)*Σmat
    z2 = -2*transpose(DΔDdκ .+ DdΔDκ)*B0*Σmat
    z3 = 8*r*Δ*transpose(κp)*B0
    z4 = 4a^2*κvec*B0*tmp3*cosmix01
    ΣΔs0 = vec(z1) .+ vec(z2) .+ vec(z3) .+ vec(z4)

    (ρlpmp = vec(ρhs2lpmp) .+ ρhs1lpmp .+ ρhs0lpmp,
     ρclpmm= vec(ρchs2lpmm) .+ ρchs1lpmm .+ ρchs0lpmm,
     ρclmmp= vec(ρchs2lmmp) .+ ρchs1lmmp .+ ρchs0lmmp,
     ρlmmm = vec(ρhs2lmmm) .+ ρhs1lmmm .+ ρhs0lmmm,
     ΣΔ    = vec(ΣΔs2) .+ ΣΔs1 .+ ΣΔs0,
     ls = ls)
end

# Banded DIAGONAL sector (lplp,lmlm,mpmp,mmmm,h) for all spherical-l'. The lplp/mpmp pieces reuse the SAME
# App-H-validated banded S-matrices as tetrad_lplp_kerr/tetrad_mpmp_kerr; lmlm/mmmm are the signmat-partner
# contractions; h = b^{(0)}-projected trace. At a=0 (b=I, aω=0) → per-l a=0 diagonals.
function _diag_banded_all(recs, ls, r, m, a, ω)
    key = _band_cache_key(recs, ls, r, m, a, ω)
    c = _band_cache_get(_DIAGBAND_CACHE, recs, key)
    c !== nothing && return c
    v = _diag_banded_all_impl(recs, ls, r, m, a, ω)
    _band_cache_put!(_DIAGBAND_CACHE, recs, key, v)
end

V_sigvec(V, CT) = CT[V[i].σ for i in 1:length(V)]

function _diag_banded_all_impl(recs, ls, r, m, a, ω)
    CT = _ct_at(recs, ls, r); aw0 = Float64(a*ω)
    A = _angular_mats_kerr(ls, m, aw0, CT); n = A.n
    Δ = r^2-2r+a^2
    V = [_diag_radial_vecs(recs[ls[i]], r, m) for i in 1:n]
    σvec = V_sigvec(V, CT)
    vget(f) = CT[getfield(V[i], f) for i in 1:n]
    Pm2_0=vget(:Pm2_0); Pp2_0=vget(:Pp2_0); DDPm2=vget(:DDPm2); DdDdPp2=vget(:DdDdPp2)
    DPm1=vget(:DPm1); DdPp1=vget(:DdPp1); Dm1Pp1=vget(:Dm1Pp1); Ddm1Pm1=vget(:Ddm1Pm1)
    h0=vget(:h0); κ0=vget(:κ0); DDκ=vget(:DDκ); DdDdκ=vget(:DdDdκ); DrDh=vget(:DrDh); DdrDdh=vget(:DdrDdh)
    sgn=A.sign; λ̂M=A.λ̂; λ̂2M=A.λ̂2; Λ̂M=A.Λ̂
    Bm2=A.B(-2); Bp2=A.B(2); Bm1=A.B(-1); Bp1=A.B(1); B0=A.B(0)
    Ss2lplp = Bp2*(Λ̂M - 2aw0*λ̂2M*A.Sin(1,0) + aw0^2*A.Sin(2,0)) + sgn*Bm2*(Λ̂M - 2aw0*λ̂2M*A.Sin(-1,0) + aw0^2*A.Sin(-2,0))
    Ss2lmlm = sgn*Ss2lplp; Ss2mpmp = Bp2; Ss2mmmm = sgn*Bm2
    SS0 = -(Bm1 + sgn*Bp1)*λ̂M + aw0*(Bm1*A.Sin(-1,0) + sgn*Bp1*A.Sin(1,0))
    Ss1mpmp = Bp1*(-λ̂2M + aw0*A.Sin(1,2)); Ss1mmmm = Bm1*(λ̂2M - aw0*A.Sin(-1,-2))
    Sh_mp  = B0*(Λ̂M - 2aw0*λ̂M*A.Sin(1,2) + aw0^2*A.Sin(0,2))*A.Cos1(2) + B0*(λ̂M*A.Sin(1,2) - aw0*A.Sin(0,2))
    Sκ_mp  = B0*(Λ̂M - 2aw0*λ̂M*A.Sin(1,2) + aw0^2*A.Sin(0,2))
    Sh_mm  = B0*(Λ̂M - 2aw0*λ̂M*A.Sin(-1,-2) + aw0^2*A.Sin(0,-2))*A.Cos1(-2) - B0*(λ̂M*A.Sin(-1,-2) - aw0*A.Sin(0,-2))
    Sκ_mm  = B0*(Λ̂M - 2aw0*λ̂M*A.Sin(-1,-2) + aw0^2*A.Sin(0,-2))
    hs2lplp = (-1/(6ω^2*Δ^2))*(transpose(Pp2_0)*Ss2lplp)
    hs2lmlm = (-1/(6ω^2*Δ^2))*(transpose(Pm2_0)*Ss2lmlm)
    mp_chain = transpose(DdDdPp2 .+ σvec.*DDPm2)
    hs2mpmp = (-1/(6ω^2))*(mp_chain*Ss2mpmp)
    hs2mmmm = (-1/(6ω^2))*(mp_chain*Ss2mmmm)
    hs1lplp = (transpose((2/Δ)*Dm1Pp1)*sgn*SS0)
    hs1lmlm = (transpose((2/Δ)*Ddm1Pm1)*SS0)
    s1mp_chain = transpose(2 .*(σvec.*DdPp1 .+ DPm1))
    hs1mpmp = (s1mp_chain*sgn*Ss1mpmp)
    hs1mmmm = -(s1mp_chain*Ss1mmmm)
    hs0lplp = transpose(-(1/(im*ω))*DrDh)*B0 - 4*transpose(DDκ)*B0
    hs0lmlm = transpose((1/(im*ω))*DdrDdh)*B0 - 4*transpose(DdDdκ)*B0
    hs0mpmp = (a/ω)*transpose(h0)*Sh_mp - 4*transpose(κ0)*Sκ_mp
    hs0mmmm = -(a/ω)*transpose(h0)*Sh_mm - 4*transpose(κ0)*Sκ_mm
    htr = transpose(h0)*B0
    (lplp = vec(hs2lplp) .+ vec(hs1lplp) .+ vec(hs0lplp),
     lmlm = vec(hs2lmlm) .+ vec(hs1lmlm) .+ vec(hs0lmlm),
     mpmp = vec(hs2mpmp) .+ vec(hs1mpmp) .+ vec(hs0mpmp),
     mmmm = vec(hs2mmmm) .+ vec(hs1mmmm) .+ vec(hs0mmmm),
     h = vec(htr), ls = ls)
end

# -------------------------------------------------------------------------------------------------------------
# Public surfaces.
# -------------------------------------------------------------------------------------------------------------
"""
    tetrad_modes_kerr(a, r0, lp, m, r; recs, S, ls) -> NamedTuple (10 modes)

The spherical-l' projected FULL set of 10 tetrad metric-perturbation modes at a≠0 — the headline
AccretionDisc/ModifiedGravity surface. Returns `(lplp,lmlm,mpmp,mmmm,ρlpmp,ρclpmm,ρclmmp,ρlmmm,ΣΔlplm,h)`
via the genuine-Kerr banded spheroidal→spherical contraction. At a→0 (b=I) reduces EXACTLY to the validated
a=0 `tetrad_modes` (gate G-K11). The l'=1 mass-dipole (not in the radiative s=−2 spheroidal list) is pure
gauge ⇒ exact zeros. `recs`/`S`/`ls` are built internally when omitted.
"""
function tetrad_modes_kerr(a, r0, lp, m, r; recs=nothing, S=nothing, ls=nothing)
    recs, S, ls = _ensure_components_inputs(a, r0, m, lp, recs, S, ls)
    lsv = collect(ls)
    if findfirst(==(lp), lsv) === nothing
        z = zero(ComplexF64)
        return (lplp=z, lmlm=z, mpmp=z, mmmm=z, ρlpmp=z, ρclpmm=z, ρclmmp=z, ρlmmm=z, ΣΔlplm=z, h=z)
    end
    ω = recs[lsv[1]].jd.ω
    db = _diag_banded_all(recs, lsv, r, m, a, ω)
    mb = _mixed_banded_all(recs, lsv, r, m, a, ω)
    j = findfirst(==(lp), lsv)
    (lplp=db.lplp[j], lmlm=db.lmlm[j], mpmp=db.mpmp[j], mmmm=db.mmmm[j],
     ρlpmp=mb.ρlpmp[j], ρclpmm=mb.ρclpmm[j], ρclmmp=mb.ρclmmp[j], ρlmmm=mb.ρlmmm[j],
     ΣΔlplm=mb.ΣΔ[j], h=db.h[j])
end

"""
    to_BLS_kerr(a, r0, lp, m, r; recs, S, ls) -> NTuple{10}

The 10 Barack-Lousto-Sago variables h̄¹..h̄¹⁰ at r for spherical-l' (a≠0), from the projected
`tetrad_modes_kerr` via the spin-general eq:hBL-inverse mapping (f=Δ/r², Σ=r²+a²; the λ̂,Λ̂ angular factors
carry the spherical-l' index). At a→0 reduces to the validated a=0 `to_BLS` 10-tuple (gate G-K11).
"""
function to_BLS_kerr(a, r0, lp, m, r; recs=nothing, S=nothing, ls=nothing)
    M = tetrad_modes_kerr(a, r0, lp, m, r; recs=recs, S=S, ls=ls)
    Δ = r^2 - 2r + a^2; Σ = r^2 + a^2; f = Δ/r^2
    Λ̂lp = Λ̂(lp); λ̂lp = λ̂(lp)
    hlplm = M.ΣΔlplm/(Σ*Δ)
    h1 = 0.5*r*f^2*(M.lplp + M.lmlm)
    h2 = 0.5*r*f^2*(M.lplp - M.lmlm)
    hmpmm = Σ*M.h - Δ*hlplm
    h3 = hmpmm/r
    hlpmp = M.ρlpmp/r; hlpmm = M.ρclpmm/r; hlmmp = M.ρclmmp/r; hlmmm = M.ρlmmm/r
    h4 = -0.5*λ̂lp*f*( hlpmp - hlpmm - hlmmp + hlmmm )
    h5 = -0.5*λ̂lp*f*( hlpmp - hlpmm + hlmmp - hlmmm )
    h6 = -r*f*hlplm
    h7 = (Λ̂lp/(2r))*(M.mpmp + M.mmmm)
    h8 = -0.5im*λ̂lp*f*( hlpmp + hlpmm - hlmmp - hlmmm )
    h9 = -0.5im*λ̂lp*f*( hlpmp + hlpmm + hlmmp + hlmmm )
    h10 = (im*Λ̂lp/(2r))*(M.mpmp - M.mmmm)
    (h1,h2,h3,h4,h5,h6,h7,h8,h9,h10)
end

"""
    tetrad_ρlpmp_kerr(a, r0, lp, m, r; recs, S, ls) -> Complex

The spherical-l' projected ρ h_{l+m+}(r) MIXED tetrad component at a≠0 (the banded contraction's ρlpmp leaf).
At a→0 reduces EXACTLY to the validated a=0 `tetrad_modes` ρlpmp (gate G-K11).
"""
function tetrad_ρlpmp_kerr(a, r0, lp, m, r; recs=nothing, S=nothing, ls=nothing)
    recs, S, ls = _ensure_components_inputs(a, r0, m, lp, recs, S, ls)
    ω = recs[collect(ls)[1]].jd.ω
    mb = _mixed_banded_all(recs, collect(ls), r, m, a, ω)
    mb.ρlpmp[findfirst(==(lp), mb.ls)]
end

"""
    tetrad_ΣΔlplm_kerr(a, r0, lp, m, r; recs, S, ls) -> Complex

The spherical-l' projected ΣΔ h_{l+l-}(r) tetrad component at a≠0. At a→0 reduces to the validated a=0
ΣΔlplm (gate G-K11).
"""
function tetrad_ΣΔlplm_kerr(a, r0, lp, m, r; recs=nothing, S=nothing, ls=nothing)
    recs, S, ls = _ensure_components_inputs(a, r0, m, lp, recs, S, ls)
    ω = recs[collect(ls)[1]].jd.ω
    mb = _mixed_banded_all(recs, collect(ls), r, m, a, ω)
    mb.ΣΔ[findfirst(==(lp), mb.ls)]
end

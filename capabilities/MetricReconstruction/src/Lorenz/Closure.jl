# Closure — the LINEAR-SYSTEM closure determining the UNKNOWN radial-function jumps per spheroidal-l:
#   x_l = [J0(P_{-1})_l, J1(P_{-1})_l, J0(κ)_l, J1(κ)_l]
# from C0/C1 continuity of the projected components h^{(s)}_{l+l+}, h^{(s)}_{m+m+}.
#
# Consolidated from KerrLorenzCirc closure.jl: the spin-general BANDED system closure_kerr (primary
# entry) + the a=0 DECOUPLED 4×4 per-l closure_a0 (fast path / regression). The jet algebra + ODE-α/β
# builders + known jumps come from Operators (jet_ops), Jumps (jumps_pm/hj1), and the frak_*/λ̂/λ̂₂/Λ̂
# CG-mixing/ladder factors (Operators).
#
# VALIDATED: closure_kerr(0.0,6,2) == App-H LOWER (7.8e-8); closure_kerr(0.99,6,2) == App-H a099 LOWER
# (6.85e-8); a=0 regression closure_kerr(0.0,…) == closure_a0(…).
#
# BigFloat: jet_ops is type-generic, so BigFloat a,r0,ω give a Complex{BigFloat} banded system M\b (a
# generic LU). The S-matrix scalars (b-matrices, 𝔰/𝔠) are Float64 (accepted cap). Parallelism: the
# n=4·nl unit-probe columns M[:,k]=eqs(e_k)+b are independent → Threads.@threads (probe_parallel),
# with BLAS.set_num_threads(1) guarding the LAPACK banded solve region.

import Harmonics: SpinWeightedSpheroidalHarmonicS, SpinWeightedSpheroidalEigenvalue,
                  SpinWeightedSphericalHarmonicY, expansion_coefficients
import NumericsCore
import Base.Threads

"""
    closure_kerr(a, r0, m; lmax=20, N=6, probe_parallel=false) -> Dict{Int,Vector}

Spin-general Lorenz-gauge jump closure. Returns a Dict mapping spheroidal-l → [J0(P_{-1}),J1(P_{-1}),
J0(κ),J1(κ)]. `a,r0` may be BigFloat (then the banded solve is Complex{BigFloat}). At a=0 the system
DECOUPLES (b=I, aω=0) and reproduces App-H LOWER (the regression gate). `probe_parallel=true` threads
the independent unit-probe columns (bit-identical to serial), guarded by BLAS.set_num_threads(1).
"""
function closure_kerr(a::Real, r0::Real, m::Integer; lmax::Int=20, N::Int=6, probe_parallel::Bool=false,
                      bofd::Bool=true, smix::Bool=true, s0h::Bool=true, dressγ::Bool=true, return_aux::Bool=false)
    bg = orbit_background(a, r0, m)
    T = typeof(bg.r0); CT = Complex{T}
    aT = bg.a; r0T = bg.r0; Ω = bg.Ω; ut = bg.ut; ω = bg.ω; aw = bg.aω
    Δ0 = bg.Δ0; calA = bg.calA; calB = bg.calB
    sw = smix ? one(T) : zero(T)

    ops = jet_ops(r0T, aT, ω, m; N=N)
    Δjet = ops.Δjet; rjet = ops.rjet; invΔjet = ops.invΔjet
    opD = ops.opD; opDd = ops.opDd; opDn = ops.opDn
    close_jet = ops.close_jet; close_jet_src = ops.close_jet_src; cj = ops.cj

    aw0 = Float64(aw)   # angular factors at Float64 aω

    lmin = max(2, abs(m))
    ls = collect(lmin:lmax)
    function bget(s, l)
        if aw0 == 0.0 || !bofd
            d = Dict{Int,Float64}(); d[l] = 1.0; return d
        end
        expansion_coefficients(SpinWeightedSpheroidalHarmonicS(s, l, m, aw0; method=:SphericalExpansion))
    end
    bval(bd, lp) = get(bd, lp, 0.0)

    nU = 4
    nl = length(ls)
    idx(li, k) = (li-1)*nU + k

    σof = Dict{Int,T}(); ΛΛ = Dict{Int,T}(); λ1 = Dict{Int,T}(); λ00 = Dict{Int,T}(); BB = Dict{Int,CT}()
    αm2=Dict{Int,Vector{CT}}(); βm2=Dict{Int,Vector{CT}}(); αp2=Dict{Int,Vector{CT}}(); βp2=Dict{Int,Vector{CT}}()
    αm1=Dict{Int,Vector{CT}}(); βm1=Dict{Int,Vector{CT}}(); α0=Dict{Int,Vector{CT}}(); β0=Dict{Int,Vector{CT}}()
    Pm2j=Dict{Int,Vector{CT}}(); Pp2j=Dict{Int,Vector{CT}}(); hj=Dict{Int,Vector{CT}}(); γκ=Dict{Int,Vector{CT}}()
    bm2=Dict(); bp2=Dict(); bm1=Dict(); bp1=Dict(); b0=Dict(); γll=Dict{Int,T}()
    for l in ls
        σof[l] = T((-1)^(l+m))
        ΛΛ[l]  = T(real(SpinWeightedSpheroidalEigenvalue(-2, l, m, aw0)))
        λ1[l]  = T(real(SpinWeightedSpheroidalEigenvalue(-1, l, m, aw0)))
        λ00[l] = T(real(SpinWeightedSpheroidalEigenvalue(0,  l, m, aw0)))
        BB[l]  = sqrt(complex(λ1[l]^2 + 4aT*m*ω - 4aT^2*ω^2))
        αm2[l],βm2[l] = ops.jets_sm2(ΛΛ[l]); αp2[l],βp2[l] = ops.jets_sp2(ΛΛ[l])
        αm1[l],βm1[l] = ops.jets_sm1(λ1[l]); α0[l],β0[l]   = ops.jets_s0(λ00[l])
        J0m2,J1m2 = jumps_pm(bg, l, -1); J0p2,J1p2 = jumps_pm(bg, l, +1)
        Pm2j[l] = close_jet(J0m2, J1m2, αm2[l], βm2[l])
        Pp2j[l] = close_jet(J0p2, J1p2, αp2[l], βp2[l])
        hj[l]   = close_jet(zero(CT), hj1(bg, l), α0[l], β0[l])
        bm2[l]=bget(-2,l); bp2[l]=bget(2,l); bm1[l]=bget(-1,l); bp1[l]=bget(1,l); b0[l]=bget(0,l)
        γ = 0.0+0im
        if dressγ
            for (j,bj) in b0[l], (k,bk) in b0[l]
                γ += bj*bk*frak_c2(0,j,k,m)
            end
        else
            γ = frak_c2(0,l,l,m)
        end
        γll[l] = T(real(γ))
        γκ[l] = jscale(T(1)/2, jmul(jadd(jmul(rjet,rjet), cj(aT^2*γll[l])), invΔjet))
    end

    # OFF-DIAGONAL κ entries (KNOWN, ∝ h): κoff[k]=Σ_{j≠k} a²γ_{jk} h_j / (2(λ0_j−λ0_k))
    κoff = Dict{Int,Vector{CT}}()
    γmat(j,k) = (g=0.0+0im; for (jp,bjp) in b0[j], (kp,bkp) in b0[k]; g += bjp*bkp*frak_c2(0,jp,kp,m); end; g)
    for k in ls
        acc = zeros(CT, N+1)
        if aT != 0
            for j in ls
                j==k && continue
                gj = γmat(j,k); dλ = λ00[j]-λ00[k]
                (gj==0.0 || abs(dλ)<1e-12) && continue
                acc .+= (aT^2*T(real(gj))/(2*dλ)) .* hj[j]
            end
        end
        κoff[k] = acc
    end

    function radial_pieces(xvec)
        Pm1=Dict{Int,Vector{CT}}(); κ=Dict{Int,Vector{CT}}(); Pp1=Dict{Int,Vector{CT}}()
        for (li,l) in enumerate(ls)
            Pm1[l] = close_jet(xvec[idx(li,1)], xvec[idx(li,2)], αm1[l], βm1[l])
            κ[l]   = close_jet_src(xvec[idx(li,3)], xvec[idx(li,4)], α0[l], β0[l], γκ[l], hj[l])
            Pp1[l] = jscale(1/BB[l], jmul(Δjet, opD(opD(Pm1[l]))))   # TS P_{+1}=(Δ/B)𝒟𝒟P_{-1}
        end
        R = Dict{Symbol,Dict{Int,Vector{CT}}}()
        for k in (:s2lplp,:s2mpmp,:s1lplp,:s1mpmp,:s0lplp,:s0mpmp_κ,:s0mpmp_h,:s0mpmp_κoff,:s0lplp_κoff)
            R[k] = Dict{Int,Vector{CT}}()
        end
        for l in ls
            R[:s2lplp][l] = jscale(-1/(6ω^2), jmul(jmul(invΔjet,invΔjet), Pp2j[l]))
            R[:s2mpmp][l] = jscale(-1/(6ω^2), jadd(jscale(σof[l], opD(opD(Pm2j[l]))), opDd(opDd(Pp2j[l]))))
            R[:s1lplp][l] = jscale(2*σof[l], jmul(invΔjet, opDn(-1, Pp1[l])))
            R[:s1mpmp][l] = jscale(2*σof[l], jadd(jscale(σof[l], opDd(Pp1[l])), opD(Pm1[l])))
            DrDh = opD(jmul(rjet, opD(hj[l])))
            R[:s0lplp][l]   = jscale(-one(CT), jadd(jscale(1/(im*ω), DrDh), jscale(4, opD(opD(κ[l])))))
            R[:s0mpmp_κ][l] = jscale(-4, κ[l])
            R[:s0mpmp_h][l] = (aT==0 || !s0h) ? zeros(CT,N+1) : jscale(aT/ω, hj[l])
        end
        for k in ls
            R[:s0mpmp_κoff][k] = jscale(-4, κoff[k])
            R[:s0lplp_κoff][k] = jscale(-4, opD(opD(κoff[k])))
        end
        R
    end

    function Smat_s2_lplp(l,lp)
        s = 0.0+0im
        for (j,bj) in bp2[l]
            term = (j==lp ? Λ̂(j) : 0.0)
            if aw0!=0.0; term += Float64(sw)*(-2aw0*λ̂2(j)*frak_s(1,0,j,lp,m) + aw0^2*frak_s2(2,0,j,lp,m)); end
            s += bj*term
        end
        for (j,bj) in bm2[l]
            term = (j==lp ? Λ̂(j) : 0.0)
            if aw0!=0.0; term += Float64(sw)*(-2aw0*λ̂2(j)*frak_s(-1,0,j,lp,m) + aw0^2*frak_s2(-2,0,j,lp,m)); end
            s += Float64(σof[l])*bj*term
        end
        s
    end
    Smat_s2_mpmp(l,lp) = bval(bp2[l],lp)
    function Smat_s1_lplp(l,lp)
        s = 0.0+0im
        for (j,bj) in bm1[l]
            t = (j==lp ? -λ̂(j) : 0.0)
            if aw0!=0.0; t += Float64(sw)*aw0*frak_s(-1,0,j,lp,m); end
            s += bj*t
        end
        for (j,bj) in bp1[l]
            t = (j==lp ? -λ̂(j) : 0.0)
            if aw0!=0.0; t += Float64(sw)*aw0*frak_s(1,0,j,lp,m); end
            s += Float64(σof[l])*bj*t
        end
        s
    end
    function Smat_s1_mpmp(l,lp)
        s = 0.0+0im
        for (j,bj) in bp1[l]
            t = (j==lp ? -λ̂2(j) : 0.0)
            if aw0!=0.0; t += Float64(sw)*aw0*frak_s(1,2,j,lp,m); end
            s += bj*t
        end
        s
    end
    Smat_s0_lplp(l,lp) = bval(b0[l],lp)
    function Smat_s0_mpmp_κ(l,lp)
        s = 0.0+0im
        for (j,bj) in b0[l]
            t = (j==lp ? Λ̂(j) : 0.0)
            if aw0!=0.0; t += Float64(sw)*(-2aw0*λ̂(j)*frak_s(1,2,j,lp,m) + aw0^2*frak_s2(0,2,j,lp,m)); end
            s += bj*t
        end
        s
    end
    function Smat_s0_mpmp_h(l,lp)
        aw0==0.0 && return 0.0+0im
        s = 0.0+0im
        for (j,bj) in b0[l]
            t0 = 0.0+0im
            for k in (j-2):(j+2)
                k<lmin && continue
                A = (j==k ? Λ̂(j) : 0.0) - 2aw0*λ̂(j)*frak_s(1,2,j,k,m) + aw0^2*frak_s2(0,2,j,k,m)
                A==0.0 && continue
                t0 += A*frak_c(2,k,lp,m)
            end
            t1 = λ̂(j)*frak_s(1,2,j,lp,m) - aw0*frak_s2(0,2,j,lp,m)
            s += bj*(t0+t1)
        end
        s
    end

    function comp_jets(xvec)
        R = radial_pieces(xvec)
        hlplp=Dict{Int,Vector{CT}}(); hmpmp=Dict{Int,Vector{CT}}()
        for lp in ls
            acc_lp = zeros(CT,N+1); acc_mp = zeros(CT,N+1)
            for l in ls
                s2l=CT(Smat_s2_lplp(l,lp)); s1l=CT(Smat_s1_lplp(l,lp)); s0l=CT(Smat_s0_lplp(l,lp))
                acc_lp .+= s2l.*R[:s2lplp][l] .+ s1l.*R[:s1lplp][l] .+ s0l.*R[:s0lplp][l]
                acc_lp .+= s0l.*R[:s0lplp_κoff][l]
                s2m=CT(Smat_s2_mpmp(l,lp)); s1m=CT(Smat_s1_mpmp(l,lp)); s0m=CT(Smat_s0_mpmp_κ(l,lp)); s0hm=CT(Smat_s0_mpmp_h(l,lp))
                acc_mp .+= s2m.*R[:s2mpmp][l] .+ s1m.*R[:s1mpmp][l] .+ s0m.*R[:s0mpmp_κ][l] .+ s0hm.*R[:s0mpmp_h][l]
                acc_mp .+= s0m.*R[:s0mpmp_κoff][l]
            end
            hlplp[lp]=acc_lp; hmpmp[lp]=acc_mp
        end
        (hlplp, hmpmp)
    end
    function srcs(lp)
        Y0 = T(real(SpinWeightedSphericalHarmonicY(0, lp, m, π/2, 0.0)))
        Y2 = T(real(SpinWeightedSphericalHarmonicY(2, lp, m, π/2, 0.0)))
        (-16*T(π)/(Δ0*ut)*calA^2*Y0, 16*T(π)/(Δ0*ut)*calB^2*Y2)
    end
    function eqs(xvec)
        hlplp,hmpmp = comp_jets(xvec)
        out = zeros(CT, 4*nl)
        for (li,lp) in enumerate(ls)
            sl,sm = srcs(lp)
            out[(li-1)*4+1] = hlplp[lp][1]
            out[(li-1)*4+2] = hlplp[lp][2]-sl
            out[(li-1)*4+3] = hmpmp[lp][1]
            out[(li-1)*4+4] = hmpmp[lp][2]-sm
        end
        out
    end

    n = 4*nl
    b = -eqs(zeros(CT, n))
    M = zeros(CT, n, n)
    if probe_parallel
        NumericsCore.blas_single_thread() do   # serialize LAPACK inside the threaded unit-probe
            Threads.@threads for k in 1:n
                e = zeros(CT, n); e[k] = one(CT)
                M[:,k] = eqs(e) .+ b
            end
        end
    else
        for k in 1:n
            e = zeros(CT, n); e[k] = one(CT)
            M[:,k] = eqs(e) .+ b
        end
    end
    x = M \ b
    sol = Dict{Int,Vector{CT}}()
    for (li,l) in enumerate(ls)
        sol[l] = [x[idx(li,1)], x[idx(li,2)], x[idx(li,3)], x[idx(li,4)]]
    end
    if return_aux
        aux = (ls=ls, γll=γll, Λ=ΛΛ, λ=λ1, λ0=λ00, BB=BB, σof=σof, ω=ω, aω=aw, Δ0=Δ0, ut=ut,
               calA=calA, calB=calB, hj=hj, κoff=κoff, bg=bg)
        return sol, aux
    end
    sol
end

"""
    closure_a0(r0, m; ls=2:10, N=6) -> Dict{Int,Vector}

a=0 fast path: DECOUPLED per-l 4×4 solve (b=I, aω=0). Identical numerics to `closure_kerr(0.0,r0,m)`
without assembling the banded system — the a=0 regression target. Type-generic in r0.
"""
function closure_a0(r0::Real, m::Integer; ls=2:10, N::Int=6)
    bg = orbit_background(0, r0, m)
    T = typeof(bg.r0); CT = Complex{T}
    r0T = bg.r0; ω = bg.ω; ut = bg.ut; Δ0 = bg.Δ0; calA = bg.calA; calB = bg.calB
    ops = jet_ops(r0T, zero(T), ω, m; N=N)
    Δjet=ops.Δjet; rjet=ops.rjet; invΔjet=ops.invΔjet
    opD=ops.opD; opDd=ops.opDd; opDn=ops.opDn
    close_jet=ops.close_jet; close_jet_src=ops.close_jet_src

    function solve_l(l)
        σ = T((-1)^(l+m))
        Λ̂l = T(Λ̂(l)); λ̂l = T(λ̂(l)); λ̂2l = T(λ̂2(l))
        Λ = T((l-1)*(l+2)); λ1 = T(l*(l+1)); λ0 = T(l*(l+1)); B = λ1
        αm2,βm2 = ops.jets_sm2(Λ); αp2,βp2 = ops.jets_sp2(Λ); αm1,βm1 = ops.jets_sm1(λ1); α0,β0 = ops.jets_s0(λ0)
        J0m2,J1m2 = jumps_pm(bg,l,-1); J0p2,J1p2 = jumps_pm(bg,l,+1); J1h = hj1(bg,l)
        Pm2jet = close_jet(J0m2,J1m2,αm2,βm2)
        Pp2jet = close_jet(J0p2,J1p2,αp2,βp2)
        hjet   = close_jet(zero(CT), J1h, α0, β0)
        γκ = jscale(T(1)/2, jmul(jmul(rjet,rjet), invΔjet))
        function comp(x)
            Pm1j0,Pm1j1,κj0,κj1 = x
            Pm1jet = close_jet(Pm1j0,Pm1j1,αm1,βm1)
            κjet   = close_jet_src(κj0,κj1,α0,β0,γκ,hjet)
            Pp1jet = jscale(1/B, jmul(Δjet, opD(opD(Pm1jet))))
            hs2lplp = jscale(-(1+σ)*Λ̂l/(6ω^2), jmul(jmul(invΔjet,invΔjet), Pp2jet))
            hs2mpmp = jscale(-1/(6ω^2), jadd(jscale(σ, opD(opD(Pm2jet))), opDd(opDd(Pp2jet))))
            hs1lplp = jscale(σ*(-(1+σ)*λ̂l), jscale(2, jmul(invΔjet, opDn(-1, Pp1jet))))
            hs1mpmp = jscale(σ*(-λ̂2l), jscale(2, jadd(jscale(σ, opDd(Pp1jet)), opD(Pm1jet))))
            hs0lplp = jadd(jscale(-1/(im*ω), opD(jmul(rjet, opD(hjet)))), jscale(-4, opD(opD(κjet))))
            hs0mpmp = jscale(-4*Λ̂l, κjet)
            (jadd(jadd(hs0lplp,hs1lplp),hs2lplp), jadd(jadd(hs0mpmp,hs1mpmp),hs2mpmp))
        end
        Y0 = T(real(SpinWeightedSphericalHarmonicY(0,l,m,π/2,0.0))); Y2 = T(real(SpinWeightedSphericalHarmonicY(2,l,m,π/2,0.0)))
        srclplp = -16*T(π)/(Δ0*ut)*calA^2*Y0; srcmpmp = 16*T(π)/(Δ0*ut)*calB^2*Y2
        function eqs(x)
            lplp,mpmp = comp(x)
            [lplp[1], lplp[2]-srclplp, mpmp[1], mpmp[2]-srcmpmp]
        end
        b = -eqs(zeros(CT,4))
        M = zeros(CT,4,4)
        for k in 1:4; e=zeros(CT,4); e[k]=one(CT); M[:,k]=eqs(e).+b; end
        M\b
    end
    sol = Dict{Int,Vector{CT}}()
    for l in ls; sol[l] = solve_l(l); end
    sol
end

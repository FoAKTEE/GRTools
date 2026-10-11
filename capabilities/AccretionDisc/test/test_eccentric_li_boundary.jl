using Test
using LinearAlgebra
using Printf
using AccretionDisc
using AccretionDisc.ForcingModes: HarmonicMode
const LIBSW = AccretionDisc.SpiralWaves
const LIBMO = AccretionDisc.Modes
isdefined(LIBSW,:radiation_impedance_row) ||
    Base.include(LIBSW,joinpath(@__DIR__,"../src/SpiralWaves/boundary.jl"))
isdefined(LIBSW,:li_boundary) ||
    Base.include(LIBSW,joinpath(@__DIR__,"../src/SpiralWaves/li_boundary.jl"))

# Independent Newtonian first-order oracle: radial/azimuthal Euler and continuity,
# MR20 src/main.tex:70-78,160-195, using the existing pressure-supported profiles.
# y=(delta h,-i delta u_r), sigma=omega-m*Omega, L=(log c_T^2)'.
function lib_newtonian(d,r,ω,m)
    Ω,Ωr = LIBMO.Omega(d,r),LIBMO.dOmegadr(d,r)
    κ2 = 4Ω^2+2r*Ω*Ωr
    σ = ω-m*Ω
    D = κ2-σ^2
    L = LIBMO.LTinv(d,r)
    [2m*Ω/(r*σ)+L D/σ;
     σ/LIBMO.csiso(d,r)^2-m^2/(r^2*σ) -(1-d.p)/r-m*κ2/(2r*Ω*σ)]
end

function lib_context(q,m,ω; R=20.0, domain=(18.0,22.0), T=Float64, χ=0.0)
    disc = LIBSW.DiscModel(χ=T(χ),r_ref=T(R),h_ref=T(0.08),pSigma=T(1.5),
        qTemp=T(q),closure=:locally_isothermal,evolution=:fixed)
    LIBSW.ModeContext(disc,HarmonicMode(m,1,T(ω));closure=:locally_isothermal,domain=T.(domain))
end

# Analytic continuous manufactured source on the ACTUAL emitted GR LI operator.
# The polynomial velocity/pressure ratio interpolates the corrected outgoing
# state at three fixed radii. No grid-dependent source or boundary-point R probe.
function lib_reflection(c,n,outward)
    lo,hi=c.domain; width=hi-lo
    zs=map((lo,(lo+hi)/2,hi)) do r
        b=LIBSW.li_boundary(c,r;outward)
        b.vector[2]/b.vector[1]
    end
    A=zs[1]; B=4zs[2]-3zs[1]-zs[3]; D=2(zs[1]+zs[3]-2zs[2])
    phase=(outward ? 1 : -1)*sign(c.mode.omega)*4
    function exact(r)
        x=(r-lo)/width; h=exp(im*phase*x)
        z=A+B*x+D*x^2
        y=h.*[1,z]
        yp=im*phase/width.*y+[0,h*(B+2D*x)/width]
        y,yp
    end
    grid=collect(range(lo,hi;length=n+1))
    mids=(grid[1:end-1]+grid[2:end])/2
    Cs=[LIBSW.master_operator(c,r;landau=false) for r in mids]
    Ss=[last(exact(r))-C*first(exact(r)) for (r,C) in zip(mids,Cs)]
    edge=outward ? hi : lo
    row=LIBSW.li_boundary(c,edge;outward).row
    y=outward ? LIBSW.solve_midpoint(grid,Cs,Ss,(1,0),first(exact(lo))[1],row,0) :
                LIBSW.solve_midpoint(grid,Cs,Ss,row,0,(1,0),first(exact(hi))[1])
    mid=(lo+hi)/2; i=n÷2+1
    a=LIBSW.li_boundary(c,mid;outward).vector
    t=LIBSW.li_boundary(c,mid;outward=!outward).vector
    ya=y[:,i]
    reflection=abs((a[1]*ya[2]-a[2]*ya[1])/(ya[1]*t[2]-ya[2]*t[1]))
    reflection,maximum(abs,y[:,i]-first(exact(mid)))
end

# Lander R1 oracle for the analytic eigenvector derivative z'=(λ'-a'-b'z)/b: difference the
# emitted operator's own eigenvector ratio z=(λ-C11)/C12 (branch continued from the returned
# eigenvalue) on a second-order stencil pointing into the domain, instead of differencing C.
function lib_transport_fd(c,r,b)
    lo,hi=c.domain; δ=1e-4*r; s=r+2δ<=hi ? 1 : -1
    z(x)=begin
        C=LIBSW.master_operator(c,x;landau=false)
        λ=LIBSW.eigen2(C)[1:2]; λ=λ[argmin(abs.(λ.-b.eigenvalue))]
        (λ-C[1,1])/C[1,2]
    end
    zr=s*(-3z(r)+4z(r+s*δ)-z(r+2s*δ))/(2δ)
    LIBSW.master_operator(c,r;landau=false)[1,2]*zr/(b.other_eigenvalue-b.eigenvalue)
end

# Lander R1: windowed incoming/outgoing ratio of a HOMOGENEOUS solution (no manufactured source),
# radiated from a Dirichlet edge into the li_boundary edge, in the corrected and frozen local bases.
function lib_transparency(q,χ,m,ω,h,outward)
    disc=LIBSW.DiscModel(χ=χ,r_ref=20.0,h_ref=h,pSigma=1.5,qTemp=q,closure=:locally_isothermal,evolution=:fixed)
    c=LIBSW.ModeContext(disc,HarmonicMode(m,1,ω);closure=:locally_isothermal,domain=(18.0,22.0))
    lo,hi=c.domain
    be=LIBSW.li_boundary(c,outward ? hi : lo;outward)
    n=ceil(Int,abs(imag(be.eigenvalue))*(hi-lo)/0.01)
    grid=collect(range(lo,hi;length=n+1))
    Cs=[LIBSW.master_operator(c,(grid[j]+grid[j+1])/2;landau=false) for j in 1:n]
    Ss=[zeros(ComplexF64,2) for _ in 1:n]
    y=outward ? LIBSW.solve_midpoint(grid,Cs,Ss,(1,0),1,be.row,0) :
                LIBSW.solve_midpoint(grid,Cs,Ss,be.row,0,(1,0),1)
    cr(u,v)=u[1]*v[2]-u[2]*v[1]
    ra,rb=lo+(hi-lo)/4,hi-(hi-lo)/4
    ec=ef=zero(ComplexF64); wsum=0.0
    for j in 1:20:n+1
        r=grid[j]; ra<r<rb || continue
        w=sin(pi*(r-ra)/(rb-ra))^2
        b=LIBSW.li_boundary(c,r;outward); bt=LIBSW.li_boundary(c,r;outward=!outward)
        C=LIBSW.master_operator(c,r;landau=false)
        Ya=(1,(b.eigenvalue-C[1,1])/C[1,2]); Yt=(1,(b.other_eigenvalue-C[1,1])/C[1,2])
        ec+=w*cr(b.vector,y[:,j])/cr(y[:,j],bt.vector)
        ef+=w*cr(Ya,y[:,j])/cr(y[:,j],Yt)
        wsum+=w
    end
    abs(ec/wsum),abs(ef/wsum)
end

@testset "LI boundary: emitted GR operator and physical branches" begin
    for q in (0.0,0.75,1.0), χ in (-0.6,0.0,0.6), (m,ω) in ((2,0.06),(2,-0.03),(0,0.04))
        c=lib_context(q,m,ω;χ)
        @test c.r_c === nothing
        for outward in (false,true)
            r=outward ? last(c.domain) : first(c.domain)
            b=LIBSW.li_boundary(c,r;outward)
            @test all(isfinite,b.row) && abs(b.mixing)<0.1
            @test abs(sum(b.row.*b.vector)) < 1e-12*norm(collect(b.row))*norm(collect(b.vector))
            @test b.transport ≈ lib_transport_fd(c,r,b) rtol=1e-6
            # Independently difference the eigenvalue branch vs TOTAL omega.
            δω=1e-7*abs(ω)
            roots=map((-δω,δω)) do shift
                cc=lib_context(q,m,ω+shift;χ)
                pair=LIBSW.eigen2(LIBSW.master_operator(cc,r;landau=false))[1:2]
                pair[argmin(abs.(pair.-b.eigenvalue))]
            end
            slope=imag(roots[2]-roots[1])/(2δω)
            @test (outward ? 1 : -1)*slope > 0
            refined=LIBSW.li_boundary(c,r;outward,radial_step=cbrt(eps())*r/2)
            @test collect(b.row) ≈ collect(refined.row) rtol=2e-8
        end
        # Explicit rows integrate through the existing closure-aware solve API.
        g=LIBSW.build_grid(c;nL=1,nR=1,n_anc=16)
        rows=((LIBSW.li_boundary(c,first(g.r);outward=false).row,1),
              (LIBSW.li_boundary(c,last(g.r);outward=true).row,0))
        sol=LIBSW.solve_mode(c;source=r->LIBSW._zero_source(c),grid=g,boundary_rows=rows,landau=false)
        @test sol.residual<1e-10 && sol.bc_residual<1e-10
    end
    for q in (0.0,0.75,1.0), ω in (-0.001,0.001), outward in (false,true)
        c=lib_context(q,0,ω)
        b=LIBSW.li_boundary(c,outward ? 22.0 : 18.0;outward)
        @test iszero(imag(b.eigenvalue))
        @test (outward ? 1 : -1)*real(b.eigenvalue)<0
        @test b.transport ≈ lib_transport_fd(c,outward ? 22.0 : 18.0,b) rtol=1e-6
    end
    c=lib_context(0.75,2,0.06)
    @test_throws ArgumentError LIBSW.li_boundary(c,18.0;outward=false,radial_step=-1)
    @test_throws ArgumentError LIBSW.li_boundary(c,18.0;outward=false,radial_step=3)
    @test_throws DomainError LIBSW.li_boundary(c,17.0;outward=false)
    rc=20.0; ω=2LIBSW.Omega(c.disc,rc)
    cc=lib_context(0.75,2,ω)
    @test_throws DomainError LIBSW.li_boundary(cc,rc;outward=true,radial_step=0.01)
    # Lander R1: the call above stops in master_inputs (omega_d=0), before the stencil guard.
    # From regular radii whose stencil straddles r_c, the guard itself must reject.
    @test_throws "stencil crosses corotation" LIBSW.li_boundary(cc,cc.r_c+0.004;outward=true,radial_step=0.01)
    @test_throws "stencil crosses corotation" LIBSW.li_boundary(cc,cc.r_c-0.004;outward=false,radial_step=0.01)
    # Lander R1: m>=1 domains normally contain corotation. Edge rows never see the Landau layer
    # and equal those of a corotation-free context with the same disc, mode and stencil.
    for q in (0.0,0.75,1.0)
        ωc=2LIBSW.Omega(lib_context(q,2,0.06).disc,20.0)
        wide=lib_context(q,2,ωc;domain=(10.0,30.0))
        @test wide.r_c !== nothing && 15<wide.r_c<25
        for (outward,r,dom) in ((false,10.0,(10.0,15.0)),(true,30.0,(25.0,30.0)))
            narrow=lib_context(q,2,ωc;domain=dom)
            @test narrow.r_c === nothing
            b=LIBSW.li_boundary(wide,r;outward)
            @test b.row == LIBSW.li_boundary(narrow,r;outward).row
            @test !iszero(imag(b.eigenvalue)) && abs(b.mixing)<0.1
        end
    end
    @test_throws DomainError LIBSW.li_boundary([0.0 1; 0 0],zeros(2,2),zeros(2,2);outward=true)
    @test_throws ArgumentError LIBSW.li_boundary([0.0 1; -1 0].+0.01im,zeros(2,2),zeros(2,2);outward=true)
    @test_throws ArgumentError LIBSW.li_boundary([0.0 1; -1 0],zeros(2,2),fill(NaN,2,2);outward=true)
    @test_throws DimensionMismatch LIBSW.li_boundary([0.0 1; -1 0],zeros(2,2),zeros(3,3);outward=true)
    @test_throws DomainError LIBSW.li_boundary([-1.0 0; 0 1],zeros(2,2),zeros(2,2);outward=true)
    ad=LIBSW.DiscModel(c.disc.reference;closure=:adiabatic,evolution=:fixed)
    ac=LIBSW.ModeContext(ad,c.mode;closure=:adiabatic,domain=c.domain)
    @test_throws ArgumentError LIBSW.li_boundary(ac,18.0;outward=false)
    # Only a precision-inheritance unit test, no BigFloat production integration.
    setprecision(BigFloat,160) do
        bc=lib_context(3//4,0,1//25;T=BigFloat)
        b=LIBSW.li_boundary(bc,big"22";outward=true)
        @test b.eigenvalue isa Complex{BigFloat}
        @test all(precision(real(x))==160 for x in b.row)
    end
end
println("PASS numerical: emitted LI outgoing/evanescent branches; qTemp=0,3/4,1; signed omega; m=0,k=1; no corotation; precision inheritance")

@testset "LI boundary: manufactured reflection convergence" begin
    for q in (0.0,0.75,1.0), (m,ω) in ((2,0.06),(2,-0.03),(0,0.04)), outward in (false,true)
        c=lib_context(q,m,ω)
        errors=[first(lib_reflection(c,n,outward)) for n in (32,64,128,256)]
        ratios=errors[1:end-1]./errors[2:end]
        @test all(3.6 .< ratios .< 4.4)
        @test last(errors)<1e-4
        @printf("REFLECTION q=%.2f m=%d omega=%+.3f %s %.9e %.9e %.9e %.9e order=%.6f\n",
                q,m,ω,outward ? "outer" : "inner",errors...,log2(last(ratios)))
    end
end
println("PASS numerical: actual LI operator manufactured interior reflection tends to zero at second order")

@testset "LI boundary: emitted-operator homogeneous transparency" begin
    # Lander R1. The manufactured test above measures discretisation of a state built from
    # li_boundary's own vector: its exact reflection is zero by construction. Here a homogeneous
    # solution carries no source. Its interior incoming/outgoing ratio is e(r)+R_edge*exp(i*Phi(r)):
    # e is the reflection a row built from that local basis would cause at r, and the edge's own
    # reflection oscillates with Phi'~2k. The sin^2-windowed mean over the central half removes the
    # oscillation. WKB parameter: h_ref, with k*dr=0.01 fixed. Corrected O(h^2), frozen O(h).
    for q in (0.0,0.75,1.0), χ in (-0.6,0.6), (m,ω) in ((2,0.06),(2,-0.03),(0,0.04)), outward in (false,true)
        corrected=Float64[]; frozen=Float64[]
        for h in (0.04,0.02,0.01)
            ec,ef=lib_transparency(q,χ,m,ω,h,outward)
            push!(corrected,ec); push!(frozen,ef)
        end
        @test all(corrected[1:end-1]./corrected[2:end] .> 3.4)
        @test all(1.8 .< frozen[1:end-1]./frozen[2:end] .< 2.3)
        @test last(corrected) < last(frozen)/20
        @printf("TRANSPARENCY q=%.2f chi=%+.1f m=%d omega=%+.3f %s corrected %.3e %.3e %.3e order %.3f | frozen %.3e %.3e %.3e order %.3f\n",
                q,χ,m,ω,outward ? "outer" : "inner",corrected...,log2(corrected[end-1]/corrected[end]),
                frozen...,log2(frozen[end-1]/frozen[end]))
    end
end
println("PASS controlled approximation: emitted GR LI operator, homogeneous waves, chi=-0.6,0.6, qTemp=0,3/4,1: corrected-basis reflection O(h^2), frozen O(h)")

@testset "LI boundary: homogeneous thermal WKB truncation order" begin
    # Exact m=0 Newtonian acoustic family: constant sigma and kappa, theta=h^2*r^-q,
    # Sigma=r^-p, p=1-q/2. The pressure equation has exact branches
    # H=r^-q exp(± integral sqrt(D/theta) dr), D=kappa^2-sigma^2.
    # For D<0 choose phase by group velocity. For D>0 choose exterior decay.
    # The test refines wavelength/background scale, NOT grid resolution.
    for q in (0.0,0.75,1.0), σ in (-1.0,1.0), κ in (0.0,2.0), outward in (false,true)
        r=1.7; n=outward ? 1 : -1; D=κ^2-σ^2
        errors=Float64[]; frozen=Float64[]
        for h in (0.08,0.04,0.02,0.01)
            θ=h^2*r^(-q); b=D/σ
            C=[-q/r b; σ/θ -q/(2r)]
            Cω=[0.0 -κ^2/σ^2-1; 1/θ 0.0]
            Cr=[q/r^2 0.0; q*σ/(r*θ) q/(2r^2)]
            z=D<0 ? im*n*sign(σ)*sqrt(-D/θ) : -n*sqrt(D/θ)
            exact=(one(z),z/b); incoming=(one(z),-z/b)
            result=LIBSW.li_boundary(C,Cω,Cr;outward)
            # Lander R1: analytic z' versus a central difference of z=(λ-C11)/C12 of C(r±δ).
            # For m=0 the (1,1) entry is the thermal gradient L=-q/r, so a' is purely thermal.
            zfam(x)=begin
                Cx=complex.([-q/x b; σ/(h^2*x^(-q)) -q/(2x)])
                λ=LIBSW.eigen2(Cx)[1:2]; λ=λ[argmin(abs.(λ.-result.eigenvalue))]
                (λ-Cx[1,1])/Cx[1,2]
            end
            δ=1e-4*r
            @test result.transport ≈ b*(zfam(r+δ)-zfam(r-δ))/(2δ)/(result.other_eigenvalue-result.eigenvalue) rtol=1e-6 atol=1e-12
            row0=last(LIBSW.outgoing_branch(C,Cω;outward))
            push!(errors,abs(sum(result.row.*exact)/sum(result.row.*incoming)))
            push!(frozen,abs(sum(row0.*exact)/sum(row0.*incoming)))
        end
        if q==0
            @test maximum(errors)<1e-14
        else
            @test all(3.8 .< errors[1:end-1]./errors[2:end] .< 4.2)
            @test all(1.9 .< frozen[1:end-1]./frozen[2:end] .< 2.1)
            @test last(errors)<last(frozen)/100
            @printf("ASYMPTOTIC q=%.2f sigma=%+.1f kappa=%.1f %s %.9e %.9e %.9e %.9e order=%.6f\n",
                q,σ,κ,outward ? "outer" : "inner",errors...,log2(errors[end-1]/errors[end]))
        end
    end
end
println("PASS controlled approximation: exact homogeneous thermal waves, propagating and evanescent; corrected reflection O(h^2), frozen O(h)")

@testset "LI boundary: weak-field and MR20 thermal amplitude" begin
    for q in (0.0,0.75,1.0), (m,ω) in ((2,3.0),(2,-2.0),(0,2.0))
        d=LIBMO.NewtonianDisc(h0=0.08,p=1.5,q=q,beta=0.0,pole_disp=0.0)
        x=1.1; N=lib_newtonian(d,x,ω,m)
        Nr=imag.(lib_newtonian(d,x+1e-20im,ω,m))/1e-20
        Nω=imag.(lib_newtonian(d,x,ω+1e-20im,m))/1e-20
        nb=LIBSW.li_boundary(N,Nω,Nr;outward=true)
        errors=Float64[]
        for R in (1e3,1e4,1e5)
            c=lib_context(q,m,ω/R^(3/2);R,domain=(R*0.9,R*1.2))
            C=LIBSW.master_operator(c,R*x;landau=false)
            scaled=[R*C[1,1] R^(3/2)*C[1,2]; sqrt(R)*C[2,1] R*C[2,2]]
            b=LIBSW.li_boundary(c,R*x;outward=true)
            push!(errors,max(norm(scaled-N)/norm(N),abs(R*b.logderivative-nb.logderivative)/abs(nb.logderivative)))
        end
        @test all(8 .< errors[1:end-1]./errors[2:end] .< 12)
        @test last(errors)<1e-4
        @printf("WEAK_FIELD q=%.2f m=%d omega=%+.1f %.9e %.9e %.9e\n",q,m,ω,errors...)
    end
    # MR20 footnote: in tightly wound waves k^2=-D/theta and amplitude
    # D'/4D - Sigma'/2Sigma - 1/2r + 3L/4. The extra L/2 is essential.
    for q in (0.0,0.75,1.0)
        errors=Float64[]
        for h in (0.04,0.02,0.01)
            d=LIBMO.NewtonianDisc(h0=h,p=1.5,q=q,beta=0.0,pole_disp=0.0)
            r=1.7; m=3; ω=Float64(m)
            C=lib_newtonian(d,r,ω,m)
            Cr=imag.(lib_newtonian(d,r+1e-20im,ω,m))/1e-20
            Cω=imag.(lib_newtonian(d,r,ω+1e-20im,m))/1e-20
            b=LIBSW.li_boundary(C,Cω,Cr;outward=true)
            expected=LIBMO.dDdr(d,r,m)/(4LIBMO.Dfreq(d,r,m))-
                LIBMO.dSigmadr(d,r)/(2LIBMO.Sigma(d,r))-1/(2r)+3LIBMO.LTinv(d,r)/4
            @test real(LIBMO.wkbC(d,r,m)) ≈ expected rtol=1e-13
            push!(errors,abs(real(b.logderivative)-expected))
            if q>0
                @test abs(real(b.logderivative)-(expected-LIBMO.LTinv(d,r)/2))>0.1
            end
        end
        @test all(3.8 .< errors[1:end-1]./errors[2:end] .< 4.2)
        @printf("THERMAL_WKB q=%.2f %.9e %.9e %.9e order=%.6f\n",q,errors...,log2(errors[end-1]/errors[end]))
    end
end
println("PASS controlled approximation: GR emitted matrix and impedance recover Newtonian limit; MR20 amplitude including extra LTinv/2")

@testset "LI boundary: free-wave Newtonian FJ/cT2 law" begin
    for q in (0.0,0.75,1.0), ω in (-2.0,3.0)
        d=LIBMO.NewtonianDisc(h0=0.12,p=1.5,q=q,beta=0.0,pole_disp=0.0)
        m=2; errors=Float64[]; fjchange=0.0
        for n in (256,512,1024)
            grid=collect(range(1.5,2.0;length=n+1))
            Cs=[complex.(lib_newtonian(d,(grid[j]+grid[j+1])/2,ω,m)) for j in 1:n]
            r=last(grid); C=lib_newtonian(d,r,ω,m)
            b=LIBSW.li_boundary(C,imag.(lib_newtonian(d,r,ω+1e-20im,m))/1e-20,
                                imag.(lib_newtonian(d,r+1e-20im,ω,m))/1e-20;outward=true)
            y=LIBSW.solve_midpoint(grid,Cs,[zeros(ComplexF64,2) for _ in 1:n],(1,0),1,b.row,0)
            flux=map(eachindex(grid)) do j
                r=grid[j]; Ω=LIBMO.Omega(d,r); σ=ω-m*Ω
                κ2=4Ω^2+2r*Ω*LIBMO.dOmegadr(d,r)
                ur=im*y[2,j]
                up=LIBSW._newtonian_delta_u_phi(m,r,Ω,κ2,σ,y[:,j],0)
                F=pi*r^2*LIBMO.Sigma(d,r)*real(ur*conj(up))
                hp=first(lib_newtonian(d,r,ω,m)*y[:,j])
                Fh=-pi*m*r*LIBMO.Sigma(d,r)/(κ2-σ^2)*imag(conj(y[1,j])*hp)
                @test F ≈ Fh rtol=1e-12
                F
            end
            law=flux./[LIBMO.csiso(d,r)^2 for r in grid]
            push!(errors,maximum(abs.(law./law[1].-1)))
            fjchange=flux[end]/flux[1]-1
        end
        @test all(3.7 .< errors[1:end-1]./errors[2:end] .< 4.3)
        @test last(errors)<2e-5
        q>0 && @test abs(fjchange)>0.1
        @printf("FLUX_LAW q=%.2f omega=%+.1f %.9e %.9e %.9e FJ_change=%+.9e\n",q,ω,errors...,fjchange)
    end
end
println("PASS numerical: source-free Reynolds FJ equals MR20 pressure current; FJ/cT2 constant under refinement, FJ varies for qTemp!=0")

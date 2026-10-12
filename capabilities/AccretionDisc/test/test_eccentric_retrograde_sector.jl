module EccentricRetrogradeSectorTests
using Test, LinearAlgebra, Printf
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW=AccretionDisc.SpiralWaves
const MO=AccretionDisc.Modes

@testset "Paper masks and disjoint physical sectors" begin
    for e in (0.0,0.06,0.12,0.5), policy in (:fr22,:fr24,:fd25)
        o=EquatorialOrbit(0.0,1000.0,e)
        mask=SW.eccentric_mode_mask(o;policy)
        newt=MO.newtonian_mode_set(MO.EccentricOrbit(1.0,e,1.0);mmax=mask.mmax,dlmax=mask.kmax)
        @test [(s.m,s.k) for s in mask.paper]==[(s.m,s.k) for s in newt]
        @test mask.selected===mask.paper
        @test length(unique(vcat(values(mask.sectors)...)))==length(mask.physical)
        @test Set(vcat(values(mask.sectors)...))==Set(mask.physical)
        @test all(s->s.m+s.k<=0,mask.sectors.excluded)
        # lander addition (review GM-F9): the physical inventory is all M(2K+1) m>=1 slots plus the K
        # independent m=0 reality representatives
        @test length(mask.physical)==mask.mmax*(2mask.kmax+1)+mask.kmax
        if policy===:fr24
            @test length(mask.paper)==(e==0 ? 150 : e<=0.06 ? 11330 : 13355)
            e>0 && @test length(mask.sectors.excluded)+length(mask.sectors.secular)==820
            @test [(s.m,s.k) for s in mask.paper]==[(s.m,s.k) for s in MO.fr24_mode_set(MO.EccentricOrbit(1.0,e,1.0))]
        end
    end
    o=EquatorialOrbit(0.0,1000.0,0.12)
    z=SW.eccentric_mode_mask(o;mask=:physical,mmax=3,kmax=4)
    @test z.selected===z.physical
    @test [(s.m,s.k) for s in z.sectors.secular]==[(1,-1)]
    @test any(s->s.omega>0,z.sectors.excluded) # k=-m is prograde in Schwarzschild
    @test any(s->s.omega<0,z.sectors.excluded)
    @test all(s->s.m>0 || s.k>0,z.physical)
    @test [(s.m,s.k) for s in z.sectors.axisymmetric]==[(0,k) for k in 1:4] # lander addition (GM-F9)
    routed=SW.eccentric_mode_mask(o;mmax=3,kmax=4,secular_rtol=0.01)
    @test any(s->s.m==2 && s.k == -2,routed.sectors.secular)
    @test Set(routed.physical)==Set(z.physical) # routing never deletes a mode
    anchor=[NewtonianMode(2,1,1.0),NewtonianMode(2,3,1.0)]
    exact=SW.eccentric_mode_mask(o;policy=:fd25,mmax=3,kmax=4,anchor_modes=anchor)
    @test [(s.m,s.k) for s in exact.paper]==[(s.m,s.k) for s in anchor]
    @test length(exact.sectors.anchor_excluded)==length(z.paper)-2
    @test Set(exact.physical)==Set(z.physical)
    @test_throws ArgumentError SW.eccentric_mode_mask(o;anchor_modes=vcat(anchor,anchor))
    @test_throws ArgumentError SW.eccentric_mode_mask(o;policy=:unknown)
    @test_throws ArgumentError SW.eccentric_mode_mask(o;kmax=-1)
    @test_throws ArgumentError SW.eccentric_mode_mask(o;secular_rtol=-1.0)
end

@testset "Manufactured signed-frequency group velocity" begin
    for w in (-3.0,3.0), normal in (-1,1)
        cs=0.4; kap=1.0; kr=sqrt(w^2-kap^2)/cs
        C=[0.0 1.0;-kr^2 0.0]; Cw=[0.0 0.0;-2w/cs^2 0.0]
        ka,Ya,kt,Yt=SW.wkb_split(C,Cw;outward=normal==1)
        @test normal*cs^2*imag(ka)/w>0
        @test normal*cs^2*imag(kt)/w<0
        @test imag(ka)≈normal*sign(w)*kr
    end
    for normal in (-1,1)
        ka,_,_,_=SW.wkb_split([0.0 1.0;4.0 0.0],zeros(2,2);outward=normal==1)
        @test normal*real(ka)<0
    end
end

function fixture(a,p,e,m,k;closure=Adiabatic)
    d=SW.DiscModel(;χ=a,r_ref=p,h_ref=0.06,pSigma=0.5,qTemp=0.0,closure,evolution=:fixed)
    pert=SW.EccentricPerturber(EquatorialOrbit(a,p,e);mass=1e-5,softening=FixedSoftening(0.03p))
    domain=(max(1.02SW.r_isco(d),0.3p),5p)
    SW.ModeContext(d,pert,m,k;closure,domain),pert
end

function slope(c,r,λ)
    δ=1e-6*max(abs(c.mode.omega),abs(c.mode.m*SW.Omega(c.disc,r)))
    roots=map((-δ,δ)) do t
        ctx=SW.ModeContext(c.disc,HarmonicMode(c.mode.m,c.mode.k,c.mode.omega+t);closure=c.closure,domain=c.domain)
        eig=SW.eigen2(SW.master_operator(ctx,r;landau=false))[1:2]
        eig[argmin(abs.(eig.-λ))]
    end
    imag(roots[2]-roots[1])/(2δ)
end

@testset "GR signed modes: boundary slopes and flux directions" begin
    for (a,p,e) in ((0.0,30.0,0.3),(0.0,30.0,0.6),(0.9,10.0,0.3),(-0.9,10.0,0.3))
        for (m,k) in ((2,-8),(2,-4),(2,-2),(2,0)), closure in (Adiabatic,LocallyIsothermalOptionA)
            c,pert=fixture(a,p,e,m,k;closure)
            rows=SW.mode_boundary_rows(c)
            for (j,r,n) in ((1,c.domain[1],-1),(2,c.domain[2],1))
                C=SW.master_operator(c,r;landau=false)
                row=rows[j][1]; y=ComplexF64[-row[2],row[1]]
                λ=closure===Adiabatic ? dot(y,C*y)/dot(y,y) : SW.li_boundary(c,r;outward=n==1).eigenvalue
                if abs(imag(λ))>1e-10
                    @test n*slope(c,r,λ)>0
                    F=SW.amf_mode(c,r,y;hmodes=SW._zero_source(c),landau=false)
                    @test n*sign(SW.omega_d(c,r))*F>0
                    if c.mode.omega<0
                        @test n*c.mode.omega/m*F>0
                    end
                else
                    @test n*real(λ)<0
                end
            end
            @printf("BRANCH a=%+.1f p=%.0f e=%.2f m=%d k=%d closure=%s omega=%+.8e wd_edges=(%+.8e,%+.8e)\n",
                a,p,e,m,k,string(closure),c.mode.omega,SW.omega_d(c,c.domain[1]),SW.omega_d(c,c.domain[2]))
        end
    end
end

@testset "GR zero-frequency continuity on a fixed source" begin
    base,_=fixture(0.0,30.0,0.3,2,-2)
    d=base.disc; r=collect(range(20.0,40.0;length=1001))
    function src(r)
        h=1e-6*exp(-((r-30)/2)^2)
        merge(SW._zero_source(base),(htt=complex(h),htt_r=complex(-(r-30)*h/2)))
    end
    states=Dict{Float64,Any}()
    for η in (-1e-5,-1e-6,0.0,1e-6,1e-5)
        mode=HarmonicMode(2,-2,η*SW.Omega_ref(d))
        c=SW.ModeContext(d,mode;closure=Adiabatic,domain=(first(r),last(r)))
        g=SW.EccentricGrid(r,c.r_c,1:0,SW.mode_cache_key(c),(;continuity=true))
        sol=SW.solve_mode(c;grid=g,source=src,landau=false,boundary_rows=SW.mode_boundary_rows(c))
        states[η]=sol.y
        @test sol.residual<1e-9 && sol.bc_residual<1e-10
        b=SW.balance_diagnostics(sol)
        @test b.amf[1]>0>b.amf[end]
        @test b.power≈mode.omega/2*b.torque rtol=2e-12 atol=1e-30
        @test abs((b.amf[end]-b.amf[1]-b.torque)/b.torque)<2e-4
        @printf("ZERO eta=%+.2e torque=%+.10e power=%+.10e fin=%+.10e fout=%+.10e balance=%.4e\n",
            η,b.torque,b.power,b.amf[1],b.amf[end],abs((b.amf[end]-b.amf[1]-b.torque)/b.torque))
    end
    near=norm(states[1e-6]-states[-1e-6])/norm(states[0.0])
    far=norm(states[1e-5]-states[-1e-5])/norm(states[0.0])
    @test near<2e-4
    @test 9.5<far/near<10.5
    @printf("CONTINUITY near=%.8e far=%.8e ratio=%.8f\n",near,far,far/near)
    # The flux-balance residual must converge under grid refinement.
    errors=Float64[]
    for n in (501,1001,2001)
        c=SW.ModeContext(d,HarmonicMode(2,-2,0.0);closure=Adiabatic,domain=(20.0,40.0))
        rr=collect(range(c.domain...;length=n))
        g=SW.EccentricGrid(rr,c.r_c,1:0,SW.mode_cache_key(c),(;balance=true))
        sol=SW.solve_mode(c;grid=g,source=src,landau=false,boundary_rows=SW.mode_boundary_rows(c))
        b=SW.balance_diagnostics(sol)
        push!(errors,abs((b.amf[end]-b.amf[1]-b.torque)/b.torque))
    end
    @test all(3.8 .< errors[1:end-1]./errors[2:end] .< 4.2)
    @test errors[end]<5e-5
    @printf("BALANCE_REFINEMENT %.8e %.8e %.8e\n",errors...)
end
end

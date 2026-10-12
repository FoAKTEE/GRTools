module FR22RetrogradeSectorTests
using Test, AccretionDisc
using AccretionDisc.ForcingModes
const MO=AccretionDisc.Modes
include("../../../../reproduce/arXiv-2407.20398/excluded_sector.jl")

@testset "Newtonian group velocity and exterior decay" begin
    d=NewtonianDisc(h0=0.1,p=0.5,q=0.0,beta=0.0,indirect=0)
    x=loggrid(d,0.1,5.0,401)
    for l in (-1,0,2,3)
        s=NewtonianMode(2,l,1.0)
        bc=MO.group_velocity_boundary_coeffs(d,x,1,length(x),s)
        for (r,n,C) in ((first(x),-1,bc[1]),(last(x),1,bc[2]))
            w=s.omega-s.m*MO.Omega(d,r)
            @test n*MO.csiso(d,r)^2*imag(C)/w>0
        end
    end
    retro=NewtonianMode(1,-1,1.0)
    ci,co=MO.group_velocity_boundary_coeffs(d,x,1,length(x),retro)
    @test imag(ci)>0>imag(co) # no legacy m=1 inner flip in a propagating region
    uniform=NewtonianDisc(h0=0.1,p=0.0,q=0.0,beta=0.0)
    @test_throws DomainError MO.group_velocity_boundary_coeffs(uniform,[1.0,2.0],1,2,NewtonianMode(1,0,1.0))
    s=NewtonianMode(2,2,1.0)
    # A narrow evanescent domain straddling corotation: check the exponential
    # part, independently of the WKB amplitude gradient.
    xx=[0.95,1.05]
    bc=MO.group_velocity_boundary_coeffs(d,xx,1,2,s)
    for (i,n) in ((1,-1),(2,1))
        amp=MO.wkbC(d,xx[i],s)-im*MO.kwave(d,xx[i],s)
        @test n*real(bc[i]-amp)<0
    end
    cache=PotentialCache(d,2,0.1,5.0;n=401)
    @test_throws ArgumentError MO.solve_lm_mode(d,NewtonianMode(2,0,1.0),x,cache)
    @test_throws ArgumentError MO.solve_lm_mode(d,NewtonianMode(0,1,1.0),x,cache;group_velocity_bc=true)
    println("NEWTONIAN coefficient and rejection checks done");flush(stdout)
    legacy=MO.solve_lm_mode(d,s,x,cache)
    @test isequal(legacy,MO.solve_lm_mode(d,s,x,cache;group_velocity_bc=false))
    for l in (0,-1)
        mode=NewtonianMode(2,l,1.0)
        f=MO.solve_lm_mode(d,mode,x,cache;group_velocity_bc=true,lindblad_regular=true,diagnostics=true)
        @test size(f.mode)==(length(x),3)
        @test all(isfinite,f.mode)
        @test MO.corotation_index(d,x,mode)==length(x) # documented outer shooting seed
        @test MO._mode_anchor(d,x,mode,:legacy)==argmin(abs.(x.-1))
        l==0 && @test mode_power(mode,1.0)==0.0
        # lander addition (review XS-F5): for omega<=0 the opt-in drops a caller's phase-gradient step
        # and anchor (the FR22 settings carry phgrad_steps=1, anchor=:corotation)
        plain=MO.solve_lm_mode(d,mode,x,cache;group_velocity_bc=true,lindblad_regular=true)
        @test isequal(plain,MO.solve_lm_mode(d,mode,x,cache;group_velocity_bc=true,lindblad_regular=true,
                                             phgrad_steps=1,phgrad_safeguard=false,anchor=:corotation))
    end
end

@testset "Forced non-isothermal weighted flux identity" begin
    # An arbitrary local state evolved by the existing first-order RHS. No
    # shooting solution or imposed radiation boundary enters this identity.
    for q in (0.0,0.75), l in (0,-1)
        d=NewtonianDisc(h0=.06,p=.5,q=q,beta=0.0,pole_disp=0.0,indirect=0)
        s=NewtonianMode(2,l,1.0)
        cache=PotentialCache(d,2,.1,5.;n=401)
        pc=MO.RegularPotential(cache); rhs=MO.RegularRHS(d,s,1,pc)
        r=.7; y=[.04,-.01,.02,.03]; dy=rhs(r,y)
        function profile_at(r,y)
            p,p1,_=MO.potential_at(pc,r); v=rhs(r,y)
            f=reshape(ComplexF64[complex(y[1],y[2])-p,
                complex(v[1],v[2])-p1,complex(y[3],y[4])],1,3)
            NewtonianExcludedSector.diagnostic_profiles(d,s,[r],f,cache;regular=true)
        end
        center=profile_at(r,y)
        slopes=Float64[]
        for step in (1e-4,5e-5)
            lo=profile_at(r-step,y-step*dy); hi=profile_at(r+step,y+step*dy)
            push!(slopes,((hi.F[1]/hi.c2[1])-(lo.F[1]/lo.c2[1]))/(2step))
        end
        limit=(4slopes[2]-slopes[1])/3
        @test limit≈center.weighted_density[1] rtol=2e-7
        q>0 && @test !isapprox(limit,center.density[1]/center.c2[1];rtol=1e-3)
        println("FORCED_FLUX_IDENTITY q=",q," l=",l," derivative=",limit,
            " corrected_rhs=",center.weighted_density[1]," uncorrected_rhs=",center.density[1]/center.c2[1])
    end
end

@testset "Static and retrograde probe torque and outgoing flux" begin
    # Use the lander's radial resolution. At n=1200 the e=.5,l=-1 torque
    # quadrature still has a 6.6e-4 error, decreasing fourfold on refinement.
    rows=NewtonianExcludedSector.probe(;n=4000,nt=512)
    @test length(rows)==12
    for r in rows
        @test r["fin_FJ0"]>0>r["fout_FJ0"]
        @test r["torque_FJ0"]<0
        haskey(r,"reference_relative_error") && @test r["reference_relative_error"]<5e-4
        r["disc"]=="FR22" && @test r["balance_relative"]<2e-4
        @test r["weighted_balance"]<2e-4
    end
end
end

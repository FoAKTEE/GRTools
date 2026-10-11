# Scoped resolution measurement, not a full FR24 mode-bank certificate.
using Test, AccretionDisc, NumericsCore, LinearAlgebra, TOML
P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes

function configuration_response(d,mode,x,cache)
    sol=solve_lm_mode(d,mode,x,cache;P.FR22_SOLVER_SETTINGS...)
    _,_,sig=vars_from_enthalpy(d,sol,x,mode,cache;lindblad_regular=true)
    # Reuse the sole torque-density implementation for a single mode.
    regular=P.RegularPotential(cache)
    phi=[first(P.potential_at(regular,r)) for r in x]
    density=P.fr24_instantaneous_torque(x,[mode],[sig],[phi],0.)
    profile=P.fr24_torque_profile(x,density)
    net=profile.at(5.)-profile.at(.05)
    (; h=sol[:,1],net)
end

@testset "FR24 adopted grid scoped refinement" begin
    d=NewtonianDisc(h0=.06,p=.5,q=.75,beta=0.,soft=.018,indirect=0)
    orb=P.EccentricOrbit(1.,.12,1.)
    rows=Dict[]
    for (m,l) in ((1,1),(2,3),(12,1))
        mode=F.NewtonianMode(m,l,1.); plan=P.fr24_mode_domain(d,mode)
        x=loggrid(d,plan.solve_domain...,4000)
        rg=loggrid(d,plan.solve_domain...,7993)
        bank=P.potential_bank(d,[orb],rg,[mode];softening=plan.softening,Nt=128,Nphi=1024)
        fine=P.harmonic_cache(bank,1,mode)
        coarse=P.PotentialCache(fine.lg[1:2:end],fine.P[1:2:end],fine.dP[1:2:end],fine.d2P[1:2:end])
        baseline=configuration_response(d,mode,x,coarse)
        source=configuration_response(d,mode,x,fine)
        xx=loggrid(d,plan.solve_domain...,7999)
        grid=configuration_response(d,mode,xx,coarse)
        logxx=log10.(xx)
        fine_h=[local_interp(log10(r),logxx,grid.h;npts=4) for r in x]
        cache_h=norm(baseline.h-source.h)/max(norm(source.h),1e-12)
        grid_h=norm(baseline.h-fine_h)/max(norm(fine_h),1e-12)
        # Net torque in FJ0 units: an absolute bound remains meaningful for
        # very weak exterior-corotation harmonics and cancelling torques.
        cache_t=abs(baseline.net-source.net)*d.h0^3
        grid_t=abs(baseline.net-grid.net)*d.h0^3
        @test cache_t < 1e-5
        @test grid_t < 1e-5
        @test cache_h < 1e-3
        @test grid_h < 1e-3
        row=Dict{String,Any}("m"=>m,"l"=>l,"cache_h"=>cache_h,"grid_h"=>grid_h,
                 "cache_torque_FJ0"=>cache_t,"grid_torque_FJ0"=>grid_t,
                 "torque_FJ0"=>baseline.net*d.h0^3,"solve_rout"=>plan.solve_domain[2])
        if !plan.corotation_in_report
            y=loggrid(d,plan.check_domain...,4000)
            b=P.potential_bank(d,[orb],loggrid(d,plan.check_domain...,3997),[mode];
                softening=plan.softening,Nt=128,Nphi=1024)
            extended=configuration_response(d,mode,y,P.harmonic_cache(b,1,mode))
            dt=abs(baseline.net-extended.net)*d.h0^3
            @test dt<1e-5
            row["boundary_torque_FJ0"]=dt
        end
        rr=[.8,1.,1.2]
        b1=P.potential_bank(d,[orb],rr,[mode];softening=plan.softening,Nt=128,Nphi=1024)
        angular=Float64[]
        for (nt,np) in ((256,1024),(128,2048))
            b2=P.potential_bank(d,[orb],rr,[mode];softening=plan.softening,Nt=nt,Nphi=np)
            err=maximum(abs,b1.caches[1].P-b2.caches[1].P)/max(1.,maximum(abs,b2.caches[1].P))
            @test err < 1e-7
            push!(angular,err)
        end
        row["angular"]=angular; push!(rows,row); println("GRID ",row)
    end
    open(io->TOML.print(io,Dict("cases"=>rows)),ARGS[1],"w")
end

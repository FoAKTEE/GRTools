using Test, AccretionDisc, NumericsCore, LinearAlgebra

@testset "FR24 integration and instantaneous conventions" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    r=[.05,.3,.9,1.2,2.,5.]
    profile=P.fr24_torque_profile(r,2r .- 3)
    exact(s)=s^2-3s+2
    @test profile.torque ≈ exact.(r) atol=1e-14
    @test profile.at(1.) == 0
    @test profile.at(.7631428283688879) ≈ exact(.7631428283688879) atol=1e-14
    @test_throws ArgumentError profile.at(5.1)
    @test_throws ArgumentError P.fr24_torque_profile(reverse(r),2r .- 3)
    @test_throws ArgumentError P.fr24_torque_profile(r,2r .- 3;a=6.)
    modes=[F.NewtonianMode(2,2,1.),F.NewtonianMode(2,3,1.),F.NewtonianMode(3,2,1.)]
    densities=[fill(.7+.2im,6),fill(-.4+.1im,6),fill(.3-.8im,6)]
    potentials=[fill(.2,6),fill(-.7,6),fill(.5,6)]
    mean_td=sum(P.fr24_instantaneous_torque(r,modes,densities,potentials,2pi*j/32) for j in 0:31)/32
    diagonal=sum(-pi .* r .* modes[j].m .* potentials[j] .* imag.(densities[j]) for j in 1:3)
    @test mean_td ≈ diagonal rtol=2e-14
    # Independent physical real fields integrated in azimuth at a fixed time.
    t=.37; physical=zeros(6)
    for az in (2pi*j/128 for j in 0:127)
        ds=sum(real.(densities[j].*exp(im*(modes[j].m*az-modes[j].omega*t))) for j in 1:3)
        grad=sum(-modes[j].m .* potentials[j] .* sin(modes[j].m*az-modes[j].omega*t) for j in 1:3)
        physical .-= (2pi/128).*r.*ds.*grad
    end
    @test P.fr24_instantaneous_torque(r,modes,densities,potentials,t) ≈ physical rtol=2e-14
end

@testset "FR24 diagnostics from solved harmonic responses" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    # One eccentric sideband at q=0 and q=.75; the q=0 bank also contains e=0.
    # Cache 2001 -> 4001 refinement is a response check, not a production accuracy claim.
    for q in (0.,.75)
        d=NewtonianDisc(h0=.1,p=.5,q=q,beta=0.,gamma=1.,soft=.03,pole_disp=1e-6,lap_nphi=1024)
        x=loggrid(d,.05,5.,8001)
        modes=[F.NewtonianMode(2,3,1.),F.NewtonianMode(2,2,1.)]
        orbits=q==0 ? [P.EccentricOrbit(1.,.1,1.),P.EccentricOrbit(1.,0.,1.)] : [P.EccentricOrbit(1.,.1,1.)]
        bank=P.potential_bank(d,orbits,loggrid(d,.05,5.,4001),modes;
            softening=F.InstantaneousScaleHeight(.3),Nt=128,Nphi=1024)
        cache=P.harmonic_cache(bank,1,modes[1])
        opts=(;reint_steps=2,phgrad_steps=1,phgrad_safeguard=true,anchor=:corotation,
                rtol=1e-10,atol=1e-12,dense=false,diagnostics=true)
        s=solve_lm_mode(d,modes[1],x,cache;opts...)
        v=P.fr24_mode_diagnostics(d,modes[1],x,s.mode,cache)
        @test all(isfinite,v.flux) && all(isfinite,v.density)
        @test v.net == v.torque[end]-v.torque[1]
        @test v.lindblad+v.corotation ≈ v.net rtol=5e-15
        @test v.power == (3/2)*v.net
        @test v.planet_power == -v.power
        @test v.FJ0 == inv(d.h0^3)
        # Fixed masks: avoid the resonance/pole and numerical boundary cutoffs.
        # Derivative uses the existing seven-point Fornberg implementation.
        ii,io=s.indices
        # FD25 oracle convention (ECC::fd25-oracle-conventions): F_J > 0 at both edges (lander, lane w16-g41).
        @test v.flux[ii] > 0 && v.flux[io] > 0
        ids=[i for i in max(ii+8,8):min(io-8,length(x)-7) if
            .15<x[i]<3.5 && abs(modes[1].omega-2P.Omega(d,x[i]))>.2 &&
            abs(P.Dfreq(d,x[i],modes[1]))>.3*P.kappa(d,x[i])^2]
        @test length(ids)>1000
        deriv=[highorder_deriv(v.flux,x,i;npts=7) for i in ids]
        ur,_,_=vars_from_enthalpy(d,s.mode,x,modes[1],cache)
        phi=first.(P.potential_at.(Ref(cache),x))
        thermal=[P.LTinv(d,r)*(v.flux[i]-pi*2r*P.Sigma(d,r)*phi[i]/
                  (modes[1].omega-2P.Omega(d,r))*real(ur[i])) for (i,r) in enumerate(x)]
        rhs=v.density[ids]+thermal[ids]
        scale=maximum(abs,rhs)
        balance=maximum(abs,deriv-rhs)/scale
        plain=maximum(abs,deriv-v.density[ids])/scale
        @test balance < 2e-3
        q==0 ? (@test plain < 2e-3) : (@test plain > .01)
        # q = 0: F_J is conserved away from the forcing (FR24 3_circular_benchmark.tex:38; FD25 oracle), so the
        # FR24 window Lindblad torque is the wave flux through the active edges (lander, lane w16-g41;
        # measured 5.2e-5; kills T_C = Delta T + Delta F, T_C = Delta T and a coorbital window).
        wave=v.flux[io]-v.flux[ii]; identity=abs(v.lindblad-wave)/abs(v.net)
        q==0 && @test identity < 5e-4
        q==0 && println("FR24_SPLIT q=$q m=2 l=3 wave=$wave TL=$(v.lindblad) identity=$identity")
        coarse=P.PotentialCache(cache.lg[1:2:end],cache.P[1:2:end],cache.dP[1:2:end],cache.d2P[1:2:end])
        sc=solve_lm_mode(d,modes[1],x,coarse;opts...)
        vc=P.fr24_mode_diagnostics(d,modes[1],x,sc.mode,coarse)
        tchange=abs(vc.net-v.net)/abs(v.net)
        fchange=norm(sc.mode[:,1]-s.mode[:,1])/norm(s.mode[:,1])
        @test tchange < 1e-3
        @test fchange < 1e-3
        # Independent angular refinement at a compact forcing holdout grid.
        rr=[.7,.9,1.,1.1,1.3]
        b1=P.potential_bank(d,orbits[1:1],rr,modes[1:1];softening=F.InstantaneousScaleHeight(.3),Nt=128,Nphi=1024)
        b2=P.potential_bank(d,orbits[1:1],rr,modes[1:1];softening=F.InstantaneousScaleHeight(.3),Nt=256,Nphi=2048)
        angular=maximum(maximum(abs,a-b)/max(1,maximum(abs,b)) for (a,b) in
            zip((b1.caches[1].P,b1.caches[1].dP,b1.caches[1].d2P),(b2.caches[1].P,b2.caches[1].dP,b2.caches[1].d2P)))
        @test angular < 1e-7
        println("FR24_RESPONSE q=$q m=2 l=3 nmask=$(length(ids)) balance=$balance uncorrected=$plain net=$(v.net) power=$(v.power) Tc=$(v.corotation) TL=$(v.lindblad) cache_torque=$tchange cache_state=$fchange angular=$angular")
        if q==0
            cc=P.harmonic_cache(bank,2,modes[2])
            se=solve_lm_mode(d,modes[2],x,cc;opts...)
            sm=solve_m_mode(d,2,x,cc;opts...)
            @test isequal(se,sm)
            ve=P.fr24_mode_diagnostics(d,modes[2],x,se.mode,cc)
            ie,oe=se.indices; wave_e=ve.flux[oe]-ve.flux[ie]; identity_e=abs(ve.lindblad-wave_e)/abs(ve.net)
            @test ve.flux[ie] > 0 && ve.flux[oe] > 0
            @test identity_e < 5e-4  # measured 1.3e-6 (lander)
            println("FR24_SPLIT q=$q m=2 l=2 wave=$wave_e TL=$(ve.lindblad) identity=$identity_e")
            # Existing circular reconstruction, integer m, same forcing cache.
            cr,cp,cs=vars_from_enthalpy(d,sm.mode,x,2,cc)
            @test ve.flux == pi .* x.^2 .* P.Sigma.(Ref(d),x) .* real.(cr.*conj.(cp))
            @test ve.density == -pi .* x .* 2 .* first.(P.potential_at.(Ref(cc),x)) .* imag.(cs)
            # Circular Laplace cache is independent of the eccentric projector.
            native=P.PotentialCache(d,2,.05,5.;n=4001)
            sn=solve_m_mode(d,2,x,native;opts...)
            vn=P.fr24_mode_diagnostics(d,modes[2],x,sn.mode,native)
            circ=abs(ve.net-vn.net)/abs(vn.net)
            @test circ < 1e-3
            @test ve.power == ve.net
            j=searchsortedfirst(x,1.4)
            physical_flux=sum(real(cr[j]*cis(2az))*real(cp[j]*cis(2az)) for az in
                              (2pi*k/64 for k in 0:63))*(2pi/64)*x[j]^2*P.Sigma(d,x[j])
            @test ve.flux[j] ≈ physical_flux rtol=2e-13
            @test_throws ArgumentError P.fr24_mode_diagnostics(NewtonianDisc(beta=.1),modes[2],x,se.mode,cc)
            @test_throws ArgumentError P.fr24_mode_diagnostics(d,modes[2],x,se.mode,cc;halfwidth=6.)
            println("FR24_CIRCULAR bitwise_response=true flux_density_exact=true relative_net_vs_Laplace=$circ net=$(ve.net)")
        end
    end
end

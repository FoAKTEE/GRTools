using Test, AccretionDisc, LinearAlgebra

@testset "FR22 eccentric path circular limit" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    # fr22-o-circular-gate-conventions: matched FR22 disc, not MR19 defaults.
    d=NewtonianDisc(h0=.1,p=0.,q=0.,beta=0.,soft=.03,indirect=0,lap_nphi=1024)
    x=loggrid(d,.1,5.,4000)
    rg=loggrid(d,.1,5.,7993)
    modes=[F.NewtonianMode(m,l,1.) for m in (2,5,10) for l in (m-1,m,m+1)]
    orbits=[P.EccentricOrbit(1.,e,1.) for e in (0.,1e-6)]
    # Deliberately use the direct double quadrature at BOTH e values; no circular dispatch.
    bank=P.potential_bank(d,orbits,rg,modes;softening=F.FixedSoftening(.03),Nt=64,Nphi=1024)
    # Lander addition (lane w17-g43 review): the e = 1e-6 sidebands are the first-order epicyclic value, not merely
    # small. With r_p = a(1 - e cos M), phi_p = M + 2e sin M and fixed softening,
    # Phi_{m,m+k} = (e/2)(2 k m Phi_m - a dPhi_m/dr_p), k = +-1, Phi_m(r; r_p) = -b_{1/2}^{(m)}(r/r_p; eps/r_p)/r_p
    # from the circular kernel, dPhi_m/dr_p by a central difference in r_p (lander probe: 2.2e-10 to 3.2e-10
    # relative with the exact derivative). Kills zeroed sidebands at small e and an apocentre-first orbit.
    rs=rg[1:40:end]; firstorder=Float64[]
    for m in (2,5,10), k in (-1,1)
        phim(rp)=[-P.laplace_b_soft(1/2,m,r/rp,.03/rp;nphi=1024)/rp for r in rs]
        pred=(1e-6/2) .* (2k*m .* phim(1.) .- (phim(1+1e-5) .- phim(1-1e-5)) ./ 2e-5)
        got=P.harmonic_cache(bank,2,F.NewtonianMode(m,m+k,1.)).P[1:40:end]
        push!(firstorder,maximum(abs,got-pred)/maximum(abs,pred))
        @test firstorder[end] < 1e-6
    end
    println("FIRST_ORDER_SIDEBAND max_relative=",maximum(firstorder))
    # Lander addition: the m = 1 diagonal through the bank equals the circular kernel with the indirect term
    # off (potential only; the certified MR core has no m = 1 solve).
    r1=rg[1:8:end]
    b1=P.potential_bank(d,orbits,r1,[F.NewtonianMode(1,1,1.)];softening=F.FixedSoftening(.03),Nt=64,Nphi=1024)
    phi1=[P.Phi_m(d,r,1) for r in r1]
    m1=[maximum(abs,b1.caches[ie,1].P-phi1)/maximum(abs,phi1) for ie in 1:2]
    @test m1[1] < 1e-10 && m1[2] < 1e-8
    println("M1_DIAGONAL e0=",m1[1]," e1e-6=",m1[2])
    rows=Dict[]
    for m in (2,5,10)
        native=P.PotentialCache(d,m,.1,5.;n=7993)
        # The reference is the existing public circular solver. It has the legacy
        # second-order Lindblad crossing, so only lindblad_regular is inapplicable.
        reference_opts=Base.structdiff(P.FR22_SOLVER_SETTINGS,(;lindblad_regular=true))
        circular=solve_m_mode(d,m,x,native;reference_opts...)
        _,_,sigc=vars_from_enthalpy(d,circular,x,m,native)
        for (ie,orb) in enumerate(orbits)
            mode=F.NewtonianMode(m,m,orb.n)
            cache=P.harmonic_cache(bank,ie,mode)
            coarse=P.PotentialCache(cache.lg[1:2:end],cache.P[1:2:end],cache.dP[1:2:end],cache.d2P[1:2:end])
            sol=solve_lm_mode(d,mode,x,coarse;P.FR22_SOLVER_SETTINGS...)
            fine=solve_lm_mode(d,mode,x,cache;P.FR22_SOLVER_SETTINGS...)
            _,_,sig=vars_from_enthalpy(d,sol,x,mode,coarse;lindblad_regular=true)
            _,_,sigf=vars_from_enthalpy(d,fine,x,mode,cache;lindblad_regular=true)
            # Max absolute error scaled by the maximum reference amplitude; no
            # componentwise relative division at nodes of the oscillatory wave.
            h_error=maximum(abs,sol[:,1]-circular[:,1])/maximum(abs,circular[:,1])
            sigma_error=maximum(abs,sig-sigc)/maximum(abs,sigc)
            cache_h=norm(sol[:,1]-fine[:,1])/norm(fine[:,1])
            cache_sigma=norm(sig-sigf)/norm(sigf)
            phi_error=maximum(abs,cache.P-native.P)/maximum(abs,native.P)
            sideband=maximum(maximum(abs,P.harmonic_cache(bank,ie,F.NewtonianMode(m,l,1.)).P)
                             for l in (m-1,m+1))/maximum(abs,native.P)
            @test size(sol)==(4000,3)  # the regular eccentric solver actually ran
            @test mode.omega/m == orb.n  # paper omega_ml is PATTERN speed
            @test phi_error < (ie==1 ? 1e-10 : 1e-8)
            @test sideband < (ie==1 ? 1e-10 : 1e-4)
            @test h_error < 1e-3
            @test sigma_error < 1e-3
            @test cache_h < 1e-4
            @test cache_sigma < 1e-4
            row=Dict("m"=>m,"e"=>orb.e,"phi"=>phi_error,"sideband"=>sideband,
                "enthalpy"=>h_error,"density"=>sigma_error,"cache_h"=>cache_h,"cache_density"=>cache_sigma)
            push!(rows,row); println("CIRCULAR ",row)
        end
    end
    # Independent angular refinements on radii closest to the softened source.
    rr=[.7,.97,1.,1.03,1.3]
    angular=Float64[]
    base=P.potential_bank(d,orbits,rr,modes;softening=F.FixedSoftening(.03),Nt=64,Nphi=1024)
    for (nt,np) in ((128,1024),(64,2048))
        fine=P.potential_bank(d,orbits,rr,modes;softening=F.FixedSoftening(.03),Nt=nt,Nphi=np)
        err=maximum(maximum(abs,a.P-b.P)/max(1.,maximum(abs,b.P)) for (a,b) in zip(base.caches,fine.caches))
        push!(angular,err); @test err < 1e-9
    end
    # Evidence capture only (lander, lane w17-g43): TOML is not an AccretionDisc test dependency, so the
    # reproduce driver (or the caller) loads it; Pkg.test never sets CHANDRA_OUTPUT.
    if haskey(ENV,"CHANDRA_OUTPUT")
        open(ENV["CHANDRA_OUTPUT"],"w") do io
            Main.TOML.print(io,Dict("cases"=>rows,"angular"=>angular,"radii"=>4000,"cache_nodes"=>[3997,7993],"Nt"=>64,"Nphi"=>1024))
        end
    end
end

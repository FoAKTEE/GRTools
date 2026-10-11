using Test, AccretionDisc
include(joinpath(@__DIR__,"../../../../reproduce/arXiv-2506.19917/fd25_setup.jl"))

@testset "FD25 disc setup and printed units" begin
    P=AccretionDisc.Modes
    discs=fd25_disc_setup()
    @test [(d.p,d.q) for d in discs]==[(1.5,1.),(1.5,0.),(.5,0.)]
    maxerr=0.
    for d in discs
        @test d.h0==.06 && d.soft==.018 && d.beta==0 && d.indirect==0
        for r in (.05,.17,.8,1.,1.4,3.,5.)
            # Hand transcription of FD25 power laws; force balance uses their
            # differentiated pressure product, not a second Omega implementation.
            sigma=r^(-d.p); cs=.06*r^(-d.q/2); h=.06*r^((1-d.q)/2)
            pressure_prime=-(d.p+d.q)*sigma*cs^2/r
            expected=(sigma,cs,h,r^(-3)+pressure_prime/(r*sigma),-d.q/r)
            actual=(P.Sigma(d,r),P.csiso(d,r),P.csiso(d,r)/(r*P.OmegaK(d,r)),P.Omega(d,r)^2,P.LTinv(d,r))
            for (a,b) in zip(actual,expected)
                err=abs(a-b)/max(abs(b),1.); maxerr=max(maxerr,err)
                @test err<2e-14
            end
        end
    end
    # FD25 4_results.tex:61,125. q_r is MASS RATIO here, never temperature q.
    norms=Dict[]
    for mass in (1e-6,2e-6,1.)
        fj=mass^2/.06^3
        rates=P.fr24_orbital_rates([(;FJ0=fj,net=0.,power=0.)];Mp=mass,a=1.,Omega_p=1.,e=.1)
        @test rates.tau0_inv ≈ mass/.06^3 rtol=2e-15
        @test fj/rates.tau0_inv ≈ mass rtol=2e-15
        push!(norms,Dict("mass_ratio"=>mass,"FJ0"=>fj,"tau0_inv"=>rates.tau0_inv))
    end
    # Lander addition (lane w17-g43 review): the package F_J0 of the modal diagnostics (G M_p = Sigma_p = a = n = 1)
    # is h_p^-3 for every FD25 disc, so q_r^2 F_J0 is the printed q_r^2/h_p^3 (4_results.tex:61).
    F=AccretionDisc.ForcingModes
    for d in discs
        r=loggrid(d,.05,5.,64); c=P.PotentialCache(d,2,.05,5.;n=64)
        @test P.fr24_mode_diagnostics(d,F.NewtonianMode(2,2,1.),r,zeros(ComplexF64,64,2),c).FJ0 == inv(.06^3)
    end
    println("FD25_SETUP max_scaled_profile_error=",maxerr," normalizations=",norms)
    # Evidence capture only (lander, lane w17-g43): TOML is not an AccretionDisc test dependency, so the
    # reproduce driver (or the caller) loads it; Pkg.test never sets CHANDRA_OUTPUT.
    if haskey(ENV,"CHANDRA_OUTPUT")
        open(io->Main.TOML.print(io,Dict("max_profile_error"=>maxerr,"normalizations"=>norms)),ENV["CHANDRA_OUTPUT"],"w")
    end
end

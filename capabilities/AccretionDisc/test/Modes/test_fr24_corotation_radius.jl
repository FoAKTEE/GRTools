using Test, AccretionDisc, NumericsCore

@testset "FR24 (35) corotation approximation on the full paper envelope" begin
    P = AccretionDisc.Modes; F = AccretionDisc.ForcingModes
    worst = (error=0.,m=0,l=0,q=0.,p=0.,h=0.,exact=0.,approx=0.,bound=0.)
    worst_all = worst
    ninside = 0; noutside = 0; maxres = 0.; maxbound = 0.; maxbound_all = 0.
    # Violations are counted per case and asserted once after the loop (lander, lane w16-g41): a
    # regression then reports a few failures, not one per case of the 2,243,640-case envelope.
    bad_slope = 0; bad_bracket = 0; bad_bound = 0; bad_window = 0
    # Superset of circular 80-cell grid, eccentric 48-cell grid, Fig 11 and B1.
    # FR24 3_circular_benchmark.tex:13; 4_eccentric_torques.tex:11-19;
    # 5_discussion.tex:100,106; app_spurious_modes.tex:9-15.
    discs = [(h,p,q) for h in (.05,.06,.07,.08,.09,.10), p in (0.,.5,1.,1.5),
                         q in (0.,.25,.5,.75,.995,1.,1.005)]
    for (h,p,q) in discs, m in 1:175, l in max(1,m-40):m+40
        d = NewtonianDisc(h0=h,p=p,q=q,beta=0.,gamma=1.)
        mode = F.NewtonianMode(m,l,1.)
        r0 = (l/m)^(-2/3); eta = h^2*(p+q)*r0^(1-q)
        approx = P.fr24_corotation_radius(d,mode)
        exact = iszero(eta) ? r0 : bisect(r -> P.Omega(d,r)-l/m,r0*(1-eta),r0)
        err = abs(exact-approx)
        # MVT bound for g(x)=x^3+eta*x^(1-q)-1 on [1-eta,1].
        # g' >= 3(1-eta)^2 - eta*max(q-1,0)*(1-eta)^(-q).
        xa = approx/r0
        slope = 3*(1-eta)^2-eta*max(q-1,0)*(1-eta)^(-q)
        bound = r0*abs(xa^3+eta*xa^(1-q)-1)/slope
        bad_slope += !(slope > 0)
        bad_bracket += !((1-eta)^3+eta*(1-eta)^(1-q)-1 <= 5e-15)
        bad_bound += !(err <= bound+2e-14)
        maxbound_all = max(maxbound_all,bound)
        maxres = max(maxres,abs(P.Omega(d,exact)/(l/m)-1))
        row = (;error=err,m,l,q,p,h,exact,approx,bound)
        err > worst_all.error && (worst_all=row)
        if .05 <= exact <= 5.
            ninside += 1
            maxbound = max(maxbound,bound)
            err > worst.error && (worst=row)
            bad_window += !(bound < .005)
        else
            noutside += 1
        end
    end
    @test (; bad_slope, bad_bracket, bad_bound, bad_window) == (; bad_slope=0, bad_bracket=0, bad_bound=0, bad_window=0)
    @test ninside + noutside == 6*4*7*13355  # 168 discs x 13355 (m,l) pairs with l >= 1
    @test maxres < 1e-13
    @test maxbound_all < .005
    for (m,l,printed) in ((2,3,.76314),(3,2,1.31037))
        d=NewtonianDisc(h0=.06,p=0.,q=0.,beta=0.,gamma=1.)
        rc=P.fr24_corotation_radius(d,F.NewtonianMode(m,l,1.))
        @test abs(rc-printed) < 5e-6
        @test isapprox(P.Omega(d,rc),l/m;rtol=1e-14)
        println("PRINTED_RADIUS m=$m l=$l r=$rc printed=$printed error=$(abs(rc-printed))")
    end
    @test_throws ArgumentError P.fr24_corotation_radius(NewtonianDisc(),F.NewtonianMode(0,1,1.))
    println("COROTATION_INSIDE n=$ninside worst=$worst max_bound=$maxbound root_residual=$maxres")
    println("COROTATION_OUTSIDE n=$noutside ALL_WORST=$worst_all max_bound_all=$maxbound_all")
end

using Test, AccretionDisc

@testset "FR24 mode and domain policy" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    rows=Dict[]
    for (e,M,K,N) in ((0.,150,0,150),(1e-6,150,40,11330),(.06,150,40,11330),(nextfloat(.06),175,40,13355),(.12,175,40,13355))
        modes=P.fr24_mode_set(P.EccentricOrbit(1.,e,1.))
        pairs=[(s.m,s.m+s.k) for s in modes]
        expected=Set((m,l) for m in 0:176 for l in -40:216 if 1<=m<=M && l>=1 && abs(l-m)<=K)
        @test Set(pairs)==expected
        @test length(modes)==length(Set(modes))==N
        @test issorted(pairs)
        @test (1,1) in pairs && (1,0) ∉ pairs
        @test all(s.omega==s.m+s.k for s in modes)
    end
    orb=P.EccentricOrbit(1.,.1,1.)
    @test length(P.newtonian_mode_set(orb;mmax=150,dlmax=8))==2514
    @test isempty(P.newtonian_mode_set(orb;mmax=0,dlmax=8))
    @test_throws ArgumentError P.newtonian_mode_set(orb;mmax=-1,dlmax=8)
    @test_throws ArgumentError P.newtonian_mode_set(orb;mmax=1,dlmax=-1)
    # Complete enumeration, cheap geometry only: no mode-bank response claim.
    for (h,p,q) in ((.06,0.,0.),(.06,.5,.75),(.1,1.5,0.),(.1,1.5,.75))
        d=NewtonianDisc(h0=h,p=p,q=q,beta=0.,soft=.3h)
        plans=[P.fr24_mode_domain(d,s) for s in P.fr24_mode_set(orb)]
        @test all(v.report_domain==(.05,5.) for v in plans)
        @test all(v.softening isa F.InstantaneousScaleHeight && v.softening.b==.3 for v in plans)
        @test all(v.solve_domain[1]<v.corotation<v.solve_domain[2] for v in plans)
        modes=P.fr24_mode_set(orb)
        residual=maximum(abs(P.Omega(d,v.corotation)-s.omega/s.m)/(s.omega/s.m) for (s,v) in zip(modes,plans))
        @test residual < 1e-12
        @test count(v.edge_policy==:fr22_edge_extension for v in plans)==4
        ext=count(!v.corotation_in_report for v in plans)
        @test ext>0
        # Lander addition (lane w17-g43 review): exterior-corotation and extended-domain counts reproduced without
        # GRTools code (evidence/fr24-mode-set/review/enumeration_independent.py); kills a dropped or narrowed 2H rule.
        @test (ext,count(v.solve_domain[2]>5 for v in plans)) ==
              Dict((.06,0.,0.)=>(60,102),(.06,.5,.75)=>(60,88),(.1,1.5,0.)=>(58,126),(.1,1.5,.75)=>(59,102))[(h,p,q)]
        row=Dict("h"=>h,"p"=>p,"q"=>q,"modes"=>length(plans),"external_corotation"=>ext,
                 "extended_domains"=>count(v.solve_domain[2]>5 for v in plans),"root_residual"=>residual,
                 "largest_outer_boundary"=>maximum(v.solve_domain[2] for v in plans))
        push!(rows,row); println("MODE_SET ",row)
    end
    d=NewtonianDisc(h0=.1,p=0.,q=0.,beta=0.)
    for m in 6:9
        v=P.fr24_mode_domain(d,F.NewtonianMode(m,1,1.))
        @test v.solve_domain==(.05,10.) && v.check_domain==(.05,12.)
    end
    for l in (-1,0)
        @test_throws ArgumentError P.fr24_mode_domain(d,F.NewtonianMode(1,l,1.))
    end
    @test_throws ArgumentError P.fr24_mode_domain(d,F.NewtonianMode(2,3,2.))
    @test P.FR24_GRID_SETTINGS.radii==4000 && P.FR24_GRID_SETTINGS.cache_nodes>=3997
    @test P.FR24_GRID_SETTINGS.Nt>2(175+40) && P.FR24_GRID_SETTINGS.Nphi>2*175
    # Evidence capture only (lander, lane w17-g43): TOML is not an AccretionDisc test dependency, so the
    # reproduce driver (or the caller) loads it; Pkg.test never sets CHANDRA_OUTPUT.
    if haskey(ENV,"CHANDRA_OUTPUT")
        open(io->Main.TOML.print(io,Dict("cases"=>rows,"omitted_nonpositive_l"=>820)),ENV["CHANDRA_OUTPUT"],"w")
    end
end

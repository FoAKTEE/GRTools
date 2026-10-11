using Test, LinearAlgebra, Printf
using AccretionDisc
using AccretionDisc.ForcingModes
if !isdefined(@__MODULE__,:EccentricModeTestSupport)
    const ECC_MODES_HELPERS_ONLY=true
    include("test_eccentric_modes.jl")
end
const EMLI=EccentricModeTestSupport
const EMLISW=AccretionDisc.SpiralWaves
isdefined(EMLISW,:li_boundary) ||
    Base.include(EMLISW,joinpath(@__DIR__,"../src/SpiralWaves/li_boundary.jl"))

function li_response_closure(sol)
    fluids=EMLISW.reconstruct_fluid(sol)
    ratios=map(eachindex(sol.r)) do i
        r=sol.r[i]; f=fluids[i]; d=sol.context.disc
        expected=EMLISW.pressure(d,r)/EMLISW.eps_profile(d,r)*f.deps
        abs(f.dp-expected)/max(abs(f.dp),abs(expected),1e-30)
    end
    maximum(ratios)
end

# Small end-to-end gate also run in fresh processes with poisoned adiabatic
# methods, and with LI dispatch deliberately mutated to adiabatic coefficients.
function li_dispatch_checks()
    @testset "ECC LI source to reconstructed closure" begin
        for q in (0.0,0.75,1.0)
            c,p=EMLI.fixture(closure=LocallyIsothermalOptionA,q=q)
            g=EMLISW.build_grid(c;nL=1,nR=1,n_anc=32)
            sol=EMLISW.solve_mode(c,p;grid=g,bank_options=(;Nq=32))
            @test sol.context.closure==sol.closure==sol.params.closure==LocallyIsothermalOptionA
            @test sol.residual<2e-9
            # Bypass the closure dispatcher to detect wrong operator/source
            # dispatch even if reconstruction still satisfies the LI law.
            inputs=[EMLISW.master_inputs(c,(g.r[i]+g.r[i+1])/2;
                hmodes=sol.source((g.r[i]+g.r[i+1])/2)) for i in 1:length(g.r)-1]
            Cs=[[EMLISW._eccentric_li_C00(x) EMLISW._eccentric_li_C01(x);
                 EMLISW._eccentric_li_C10(x) EMLISW._eccentric_li_C11(x)] for x in inputs]
            Ss=[[EMLISW._eccentric_li_S0(x),EMLISW._eccentric_li_S1(x)] for x in inputs]
            _,li_residual=EMLISW.midpoint_residual(g.r,Cs,Ss,sol.y)
            @test li_residual<2e-9
            law=li_response_closure(sol)
            @test law<2e-11
            for (i,outward) in ((1,false),(length(g.r),true))
                row=EMLISW.li_boundary(c,g.r[i];outward).row
                @test row==sol.params.rows[outward ? 2 : 1][1]
                @test abs(sum(row.*sol.y[:,i]))<1e-12
            end
            h=sol.source(g.r[length(g.r)÷2])
            @test abs(h.htr)>1e-12 && abs(h.hrphi)>1e-12
            @printf("LI_CLOSURE,%.2f,%d,%.12e,%.12e\n",q,length(g.r),law,sol.residual)
            @printf("LI_DIRECT_RESIDUAL,%.2f,%.12e\n",q,li_residual)
        end
    end
    println("PASS LI dispatch: source to context to thermal radiation rows to common BVP to OPTION A reconstruction")
end

li_dispatch_checks()
if get(ENV,"ECC_LI_DISPATCH_ONLY","false")!="true"
    @testset "ECC LI manufactured and physical forced convergence" begin
        for q in (0.0,0.75,1.0)
            c,p=EMLI.fixture(closure=LocallyIsothermalOptionA,q=q)
            EMLI.manufactured_checks(c)
            sol=EMLI.forced_checks(c,p)
            @test li_response_closure(sol)<2e-11
            # The generic engine can reuse its operator, but a different
            # closure must never accept that cache, even when qTemp=0.
            g=EMLI.uniform_grid(c,128)
            cache=EMLISW.mode_operator_cache(c,g)
            cached=EMLISW.solve_mode(c;source=sol.source,grid=g,cache,
                                    boundary_rows=EMLISW.mode_boundary_rows(c))
            @test cached.y==sol.y
            ad=EMLISW.DiscModel(c.disc.reference;closure=Adiabatic,evolution=:fixed)
            ac=EMLISW.ModeContext(ad,c.mode;closure=Adiabatic,domain=c.domain)
            @test_throws ArgumentError EMLISW.solve_mode(ac;source=sol.source,
                grid=EMLI.uniform_grid(ac,128),cache,boundary_rows=EMLISW.mode_boundary_rows(ac))
        end
        EMLI.precision_checks(LocallyIsothermalOptionA)
    end
    println("PASS LI convergence: qTemp=0,3/4,1; second-order manufactured and forced continuous defects; closure cache isolation")
    @testset "ECC LI signed and corotating physical modes" begin
        for (m,k) in ((2,-6),(2,0),(0,1))
            c,p=EMLI.fixture(closure=LocallyIsothermalOptionA,m=m,k=k,domain=(10.0,42.0))
            g=EMLISW.build_grid(c;nL=1,nR=1,nC=65,n_anc=64)
            sol=EMLISW.solve_mode(c,p;grid=g,bank_options=(;Nq=64))
            @test sol.residual<1e-8 && sol.bc_residual<1e-10
            @test li_response_closure(sol)<2e-11
            EMLI.sampled_source_checks(sol,p,64;indices=filter(!isnothing,(length(sol.r)÷2,EMLI.nearest_corotation_node(sol))))
            @test (c.r_c===nothing)==isempty(g.central)
            c.r_c===nothing || @test c.r_c in g.r
            @printf("LI_GRID,%d,%d,%.12e,%s,%d,%.12e\n",m,k,c.mode.omega,string(c.r_c),length(g.r),sol.residual)
        end
    end
    println("PASS LI physical grids: signed omega, corotation present/absent, axisymmetric nonstatic source")
    # Lander review: the dispatch gate above, repeated in a fresh process with
    # every adiabatic and circular kernel redefined to throw.
    @testset "ECC LI adiabatic/circular isolation (fresh process)" begin
        @test EMLI.isolation_run("li")
    end
end

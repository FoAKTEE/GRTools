module EccentricRegressionTests
using Test, LinearAlgebra
include(joinpath(@__DIR__,"../../../reproduce/_emri_eccentric/circular_regression.jl"))
const CR=CircularRegression
const SW=CR.SW
const RH=CR.RH

@testset "EGD driver adaptation preserves numerical checks" begin
    @test CR.adapt(:(SW.solve_mode(bg,2))) ==
          Expr(:call,GlobalRef(CR,:solve_mode),:bg,2)
    @test CR.adapt(:(SW.build_grid(bg,2))) == :(SW.build_grid(bg,2))
    @test CR.adapt(:(draw(x))) == :(nothing)
    @test CR.adapt(:(cache_write(1,x))) == :(nothing)
    signature=:(function draw(x); savefig(x,"forbidden.png"); end)
    @test CR.adapt(signature).args[1]==signature.args[1]
    ex=CR.adapt(:(ReproHarness.conclude!(opts;paper="2608.17003",id="fig10",script="x",rtol=1e-8)))
    @test ex.args[1]==GlobalRef(CR,:conclude)
    @test count(x->x isa Expr && x.head===:parameters,ex.args)==1
    @test RH.compare_values(1.0,1.0+2e-8;rtol=1e-8,atol=0.0)[1]==false
    @test_throws ErrorException CR.run_figure("fig99","unused")
    @test_throws ErrorException CR.run_figure("fig01",joinpath(CR.PAPER,"expected"))
    mktempdir() do tmp
        alias=joinpath(tmp,"alias")
        symlink(joinpath(CR.PAPER,"expected"),alias)
        @test_throws ErrorException CR.run_figure("fig01",alias)
        @test_throws ErrorException CR.run_figure("fig01",tmp)
        @test islink(alias)
    end
end

@testset "Shared source bank: identity, interpolation and isolation" begin
    bg=SW.DiscBackground(χ=0.4,aEMRI=20.0,h0=0.08,𝔭=1.5,𝔮=0.75,b=0.3,εp=1e-3)
    frozen=SW.build_grid(bg,2;nL=1,nR=1,nC=17,n_anc=32,extent=:fiducial,r_in=8.0,r_out=45.0)
    c,p,g=CR.context(bg,2,frozen;e=0.02,mass=1e-5)
    modes=[SW.ModeContext(c.disc,p,m,k;closure=CR.Adiabatic,domain=c.domain).mode
           for (m,k) in ((2,0),(3,1))]
    builder=mode->SW.ModeContext(c.disc,p,mode;closure=CR.Adiabatic,domain=c.domain)
    grid_builder=ctx->SW.EccentricGrid(copy(g.r),ctx.r_c,g.central,SW.mode_cache_key(ctx),(;))
    samples=sort!(unique!(vcat(g.r,(g.r[1:end-1]+g.r[2:end])/2)))
    # Fixed azimuth resolution and one radial block make the independent
    # per-mode and shared sample quadratures identical.
    options=(;Nq=16,Nphi=4096,mmax=3,blocksize=1)
    bank=SW.source_azimuth_bank(c.disc,p,samples;options...)
    shared=SW.solve_modes(builder,p,modes;bank,grid_builder)
    separate=[SW.solve_mode(builder(m),p;grid=grid_builder(builder(m)),
                           radial_batch_size=32,bank_options=options) for m in modes]
    @test all(isequal(a.y,b.y) for (a,b) in zip(shared,separate))
    @test all(isequal(a.source.jets,b.source.jets) for (a,b) in zip(shared,separate))
    @test all(isequal(a.residual,b.residual) for (a,b) in zip(shared,separate))
    radial=SW.source_radial_bank(c.disc,p,c.domain,modes;Nq=16,method=:direct)
    interpolated=SW.solve_modes(builder,p,modes;bank=radial,grid_builder)
    units=SW._source_units(p,true)
    budget=maximum(abs(a-b)/(1e-8*abs(b)+1e-12*u)
        for (i,mode) in enumerate(modes) for r in samples
        for (a,b,u) in zip(values(interpolated[i].source(r)),values(shared[i].source(r)),units))
    @test budget<=1
    response=maximum(norm(a.y-b.y)/norm(b.y) for (a,b) in zip(interpolated,shared))
    @test response<1e-8
    automatic=SW.solve_modes(builder,p,reverse(modes);grid_builder,
        bank_options=(;Nq=16,method=:direct),transform=s->s.y)
    @test automatic==[s.y for s in reverse(interpolated)]
    fft=SW.solve_modes(builder,p,modes;grid_builder,bank_options=(;Nq=16))
    @test maximum(norm(a.y-b.y)/norm(b.y) for (a,b) in zip(fft,interpolated))<1e-8
    events=Tuple{Symbol,Int}[]
    custom=ctx->map(SW.mode_boundary_rows(ctx)) do (row,rhs)
        (ctx.mode.m .* row,rhs)
    end
    custom_solutions=SW.solve_modes(p,modes;bank,boundary_builder=custom,
        grid_builder=ctx->begin
            push!(events,(:grid,ctx.mode.m)); grid_builder(ctx)
        end) do mode
        push!(events,(:context,mode.m)); builder(mode)
    end
    @test events==[(:context,2),(:grid,2),(:context,3),(:grid,3)]
    @test all(sol.params.rows==custom(sol.context) for sol in custom_solutions)
    @test maximum(norm(a.y-b.y)/norm(b.y) for (a,b) in zip(custom_solutions,shared))<1e-8
    @test SW.solve_modes(builder,p,[modes[2],modes[1],modes[2]];bank,grid_builder,
                         transform=s->s.context.mode)==[modes[2],modes[1],modes[2]]
    # Distinct mode domains and noncoincident nodes exercise the reason for
    # radial-bank reuse: no common corotation grid is assumed.
    different=mode->SW.ModeContext(c.disc,p,mode;closure=CR.Adiabatic,
        domain=mode.m==2 ? (8.0,39.0) : (11.0,45.0))
    different_grid=ctx->SW.EccentricGrid(collect(range(ctx.domain...;length=28+ctx.mode.m)),
        ctx.r_c,1:0,SW.mode_cache_key(ctx),(;))
    varying=SW.solve_modes(different,p,modes;grid_builder=different_grid,
        bank_options=(;Nq=16,method=:direct))
    @test varying[1].source.bank === varying[2].source.bank
    @test varying[1].source.bank.grid.edges[[1,end]]==[8.0,45.0]
    @test length(varying[1].r)!=length(varying[2].r)
    for (mode,sol) in zip(modes,varying)
        ctx=different(mode)
        direct=SW.solve_mode(ctx,p;grid=different_grid(ctx),bank_options=(;Nq=16))
        @test norm(sol.y-direct.y)/norm(direct.y)<1e-8
    end
    @test_throws ArgumentError SW.solve_modes(builder,p,modes[1:0])
    # Lander L-G38 (mutant M1): pass the matching grid builder so that each
    # rejection below is caused by the check it names, not by an unsampled grid.
    @test_throws ArgumentError SW.solve_modes(_->c,p,[modes[2]];bank,grid_builder)
    @test_throws ArgumentError SW.solve_modes(builder,p,modes;bank,bank_options=(;Nq=4))
    other=SW.EccentricPerturber(p.orbit;mass=2p.mass,softening=p.softening)
    @test_throws ArgumentError SW.solve_modes(builder,other,modes;bank,grid_builder)
    missing=SW.source_azimuth_bank(c.disc,p,[8.0,45.0];Nq=16,mmax=3)
    @test_throws ArgumentError SW.solve_modes(builder,p,modes;bank=missing,grid_builder)
    println("SHARED_BANK bitwise_states=true bitwise_jets=true source_budget_ratio=",budget,
            " response_relative_error=",response," modes=",length(modes))
end

@testset "Adiabatic eccentric circular regression and removable limit" begin
    bg=SW.DiscBackground(χ=0.4,aEMRI=20.0,h0=0.08,𝔭=1.5,𝔮=0.75,b=0.3,εp=1e-3)
    grid=SW.build_grid(bg,2;nL=1,nR=1,nC=65,n_anc=64,extent=:fiducial,r_in=8.0,r_out=45.0)
    c,p,g=CR.context(bg,2,grid;mass=1e-5)
    # The frozen Float64 EGD drivers use :auto (the sparse backend), even
    # though the general eccentric solver deliberately defaults to :banded.
    adapted=CR.solve_mode(bg,2;grid,mEMRI=p.mass)
    explicit=CR.solve_mode(bg,2;grid,mEMRI=p.mass,linsolve=:auto)
    @test adapted.params.linsolve===:auto
    @test adapted.y==explicit.y && adapted.residual==explicit.residual
    batched=only(CR.batch_modes(identity,bg,2:2;mEMRI=p.mass,
        grid_kwargs=(;nL=1,nR=1,nC=65,n_anc=64,extent=:fiducial,r_in=8.0,r_out=45.0)))
    @test batched.params.linsolve===:auto
    @test norm(batched.y-adapted.y)/norm(adapted.y)<1e-8
    modes=[SW.ModeContext(c.disc,p,2,k;closure=CR.Adiabatic,domain=c.domain).mode for k in -2:2]
    bank=SW.source_azimuth_bank(c.disc,p,[12.0,20.0,31.0];Nq=16,mmax=2)
    jets=SW.source_modes(bank,modes)
    carrier=maximum(abs(v) for j in jets[:,3] for v in values(j))
    side=maximum(abs(v) for k in (1,2,4,5) for j in jets[:,k] for v in values(j))
    @test side/carrier<1e-13
    circ=SW.solve_mode(bg,2;grid,mEMRI=p.mass,linsolve=:banded)
    totals=Float64[]; errors=Float64[]
    baseline=nothing
    for e in (0.0,1e-3,1e-4,1e-5)
        ce,pe,ge=CR.context(bg,2,grid;e,mass=p.mass)
        sol=SW.solve_mode(ce,pe;grid=ge,bank_options=(;Nq=16))
        @test sol isa SW.EccentricModeSolution
        @test sol.closure==CR.Adiabatic && sol.params.linsolve==:banded
        @test sol.residual<1e-8 && sol.bc_residual<1e-10
        @test sol.params.rows==SW.mode_boundary_rows(ce)
        push!(errors,norm(sol.y-circ.y)/norm(circ.y))
        push!(totals,SW.total_torque(sol))
        baseline===nothing && (baseline=sol)
        println("E_TO_ZERO e=",e," state_relative_error=",last(errors)," torque=",last(totals))
    end
    @test errors[1]<1e-8
    @test 80<errors[2]/errors[3]<120
    @test errors[4]<errors[3]
    response_side=Float64[]
    for k in (-2,-1,1,2)
        cs,ps,gs=CR.context(bg,2,grid;k,mass=p.mass)
        ss=SW.solve_mode(cs,ps;grid=gs,bank_options=(;Nq=16))
        push!(response_side,norm(ss.y)/norm(baseline.y))
        @test last(response_side)<1e-12
        @test ss.params.rows==SW.mode_boundary_rows(cs)
    end
    refT=SW.total_torque(bg,circ)
    @test totals[1] ≈ refT rtol=1e-8
    F=SW.amf(baseline); Fc=SW.amf(bg,circ)
    @test norm(F-Fc)/norm(Fc)<1e-8
    println("SIDEBANDS e=0 source_relative_max=",side/carrier," response_relative_norms=",response_side)
    println("CIRCULAR torque=",totals[1]," oracle=",refT," flux_relative_error=",norm(F-Fc)/norm(Fc))
    out=get(ENV,"ECC_REGRESSION_EVIDENCE","")
    if !isempty(out)
        write(joinpath(out,"limits.json"),RH.json_write(Dict(
            "e"=>[0.0,1e-3,1e-4,1e-5],"state_relative_error"=>errors,
            "torque"=>totals,"circular_torque"=>refT,"sideband_relative_max"=>side/carrier,
            "sideband_response_relative_norms"=>response_side,
            "flux_relative_error"=>norm(F-Fc)/norm(Fc),"grid_nodes"=>length(g.r),
            "scope"=>"m=2 k=0 response; source m=2 k=-2:2; chi=.4 a=20 h=.08 p=1.5 q=.75")))
    end
end

# Lander L-G38 (mutant M8): every physics entry point the fourteen frozen drivers
# call through `SW.` must be intercepted; only these non-solving helpers may remain.
@testset "Every EGD response call is adapted (parse only)" begin
    allowed=(:DiscBackground,:torque_units,:mode_torques,:mode_lindblad_torques,:m_cut,
             :r_isco,:compute_budget,:build_grid,:synthesize,:synthesize_cartesian,:arm_phase)
    function left!(acc,ex)
        ex isa Expr || return acc
        if ex.head===:call && ex.args[1] isa Expr && ex.args[1].head===:. &&
           ex.args[1].args[1]===:SW && ex.args[1].args[2] isa QuoteNode
            push!(acc,ex.args[1].args[2].value)
        end
        foreach(x->left!(acc,x),ex.args)
        acc
    end
    files=[[joinpath(CR.PAPER,"fig"*lpad(string(i),2,'0')*".jl") for i in 1:14];
           joinpath(CR.PAPER,"_profile_assembly.jl")]
    for file in files
        remaining=left!(Set{Symbol}(),CR.adapt(Meta.parseall(read(file,String))))
        @test issubset(remaining,allowed)
    end
    for name in CR.BACKEND_CALLS
        @test CR.adapt(Expr(:call,Expr(:.,:SW,QuoteNode(name)),:x)).args[1]===GlobalRef(CR,name)
    end
end

# Lander L-G38 (root cause RC-A, obligation gr-o-e0-boundary-convention): with an
# edge inside the Landau profile the C-66 rows (bare edge operator) and the EGD
# (66)-(68) rows (Landau-shifted edge matrix) define different e=0 problems. The
# eccentric context solver's default rows reproduce the circular solve; the C-66
# rows do not. The broken test flips when the boundary convention is fixed.
@testset "e=0 edge-row convention with an edge near corotation" begin
    bg=SW.DiscBackground(χ=0.4,aEMRI=20.0,h0=0.08,𝔭=1.5,𝔮=0.75,b=0.3,εp=1e-3)
    grid=SW.build_grid(bg,2;nL=1,nR=1,nC=65,n_anc=64,extent=:fiducial,r_in=17.0,r_out=45.0)
    c,p,g=CR.context(bg,2,grid;mass=1e-5)
    circ=SW.solve_mode(bg,2;grid,mEMRI=p.mass,linsolve=:banded)
    shift=SW.landau_shift(c,first(g.r))/abs(SW.omega_d(c,first(g.r)))
    @test shift>1e-4
    solve(rows)=only(SW.solve_modes(_->c,p,[c.mode];grid_builder=_->g,boundary_builder=rows,
                                    bank_options=(;Nq=2,method=:direct)))
    egd=solve(_->nothing)
    c66=solve(SW.mode_boundary_rows)
    egd_err=norm(egd.y-circ.y)/norm(circ.y); c66_err=norm(c66.y-circ.y)/norm(circ.y)
    @test egd_err<1e-10
    @test abs(SW.total_torque(egd)/SW.total_torque(bg,circ)-1)<1e-10
    @test_broken c66_err<1e-8
    println("EDGE_ROW landau_over_omega_d=",shift," EGD_rows_state_error=",egd_err,
            " C66_rows_state_error=",c66_err)
end
end

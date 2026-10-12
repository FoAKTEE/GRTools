"""
Run the frozen EGD reductions with an adiabatic eccentric response backend.

Only reproduction call sites are adapted; no package method or golden is
overwritten. The original drivers own configurations, reductions and tolerances.
Their plotting/cache writes are disabled. See README.md for the call graph.
"""
module CircularRegression
using LinearAlgebra, Printf, SHA
# The figure drivers use the umbrella reproduction environment. Direct package
# tests use AccretionDisc's environment, without adding an umbrella dependency.
if Base.find_package("GrTools")===nothing
    using AccretionDisc, Spacetime
    const AD=AccretionDisc
    const ST=Spacetime
else
    using GrTools
    const AD=GrTools.AccretionDisc
    const ST=GrTools.Spacetime
end
using .AD.ForcingModes: Adiabatic, FixedSoftening
using .ST: EquatorialOrbit
const SW = AD.SpiralWaves
include(joinpath(@__DIR__, "../_lib/repro_harness.jl"))
const RH = ReproHarness
const PAPER = normpath(joinpath(@__DIR__, "../arXiv-2608.17003"))
const OUTPUT = Ref("")
const TRACE = Dict{String,Any}[]
cpu_seconds() = Float64(ccall(:clock, Clong, ()))/1e6

function context(bg, m, grid; e=0.0, k=0, mass=1.0)
    d = SW.DiscModel(bg; closure=Adiabatic, evolution=:fixed)
    p = SW.EccentricPerturber(EquatorialOrbit(bg.χ,bg.aEMRI,e);
                             mass, softening=FixedSoftening(SW.softening_length(bg)))
    c = SW.ModeContext(d,p,m,k; closure=Adiabatic, domain=(first(grid.r),last(grid.r)))
    # A frozen radial discretization isolates regression of the eccentric
    # source/operator/boundary/diagnostic chain from grid-selection changes.
    g = SW.EccentricGrid(copy(grid.r),c.r_c,grid.blocks.C,SW.mode_cache_key(c),
                         (;frozen_egd_grid=true))
    c,p,g
end

const SINGLE_BANK = Ref{Any}(nothing)

# EGD (66)-(68): solver.jl:572-579 closes with the Landau-shifted edge
# eigenvector. Keep this convention local to the frozen regression; the
# production eccentric solver's certified C-66 rows remain unregularized.
function egd_boundary_rows(c::SW.ModeContext{T}; landau=true) where {T}
    map(((first(c.domain),false),(last(c.domain),true))) do (r,outward)
        C=SW.master_operator(c,r;landau)
        _,_,row=SW.outgoing_branch(C,SW.omega_d(c,r);outward)
        (row,zero(Complex{T}))
    end
end

function record_mode(bg, grid, sol, cpu, wall)
    c=sol.context; m=c.mode.m
    row=Dict{String,Any}("m"=>m,"k"=>0,"e"=>0.0,"chi"=>bg.χ,"a"=>bg.aEMRI,
        "h0"=>bg.h0,"p"=>bg.𝔭,"q"=>bg.𝔮,"b"=>bg.b,"eps_p"=>bg.εp,
        "nodes"=>length(sol.r),"nL"=>grid.nL,"nR"=>grid.nR,"nC"=>grid.nC,
        "omega"=>c.mode.omega,"r_c"=>c.r_c,"cpu_seconds"=>cpu_seconds()-cpu,
        "wall_seconds"=>time()-wall,"residual"=>sol.residual,
        "bc_residual"=>sol.bc_residual,
        "path"=>"shared source_radial_bank -> ModeContext -> explicit EGD edge rows -> solve_modes")
    push!(TRACE,row)
    isempty(OUTPUT[]) || open(joinpath(OUTPUT[],"mode_costs.jsonl"),"a") do io
        println(io,replace(RH.json_write(row),"\n"=>""))
    end
    sol
end

function regression_bank(d,p,domain,modes)
    cpu=cpu_seconds(); wall=time()
    bank=SW.source_radial_bank(d,p,domain,modes;Nq=2,method=:direct)
    # Independent direct samples audit all requested modes against the source
    # budget, including the orbit and both domain edges. No tolerance is tuned.
    radii=sort!(unique!([range(domain...;length=7)...,p.orbit.p]))
    direct=SW.source_azimuth_bank(d,p,radii;Nq=2,mmax=maximum(m.m for m in modes))
    ref=SW.source_modes(direct,modes)
    jets=SW.source_modes(bank,modes,radii)
    units=SW._source_units(p,true)
    ratio=maximum(abs(jets[i][j]-ref[i][j])/(1e-8*abs(ref[i][j])+1e-12*units[j])
                  for i in eachindex(jets) for j in 1:14)
    ratio<=1 || error("shared radial source exceeds budget: $ratio")
    row=Dict("modes"=>length(modes),"radii"=>length(bank.grid.rs),
        "source_budget_ratio"=>ratio,"cpu_seconds"=>cpu_seconds()-cpu,
        "wall_seconds"=>time()-wall,"a"=>p.orbit.p,"chi"=>p.orbit.a)
    isempty(OUTPUT[]) || open(joinpath(OUTPUT[],"bank_costs.jsonl"),"a") do io
        println(io,replace(RH.json_write(row),"\n"=>""))
    end
    bank
end

function solve_mode(bg, m; grid=SW.build_grid(bg,m), threads=false,
                    mEMRI=1.0, landau=true, linsolve=:auto)
    bg isa SW.DiscBackground{Float64} || error("This regression is Float64 only")
    c,p,g=context(bg,m,grid;mass=mEMRI)
    cpu=cpu_seconds(); wall=time()
    # The profile drivers request modes individually. Keep one bank, covering
    # their full mode range and fiducial domain, across those call sites.
    key=(bg,mEMRI)
    cached=SINGLE_BANK[]
    if cached===nothing || cached.key!=key || !(c.mode in cached.bank.modes) ||
       first(g.r)<first(cached.bank.grid.edges) || last(g.r)>last(cached.bank.grid.edges)
        domain=(min(SW.r_isco(bg),first(g.r)),max(5bg.aEMRI,last(g.r)))
        modes=[SW.ModeContext(c.disc,p,n,0;closure=Adiabatic,domain).mode for n in 1:max(256,m)]
        bank=regression_bank(c.disc,p,domain,modes)
        cached=(;key,bank,p); SINGLE_BANK[]=cached
    end
    sol=only(SW.solve_modes(_->c,cached.p,[c.mode];bank=cached.bank,grid_builder=_->g,
        boundary_builder=ctx->egd_boundary_rows(ctx;landau),landau,linsolve))
    record_mode(bg,grid,sol,cpu,wall)
end

function batch_modes(transform,bg,ms;grid_kwargs=(;),mEMRI=1.0,threads=false,
                     landau=true,linsolve=:auto)
    bg isa SW.DiscBackground{Float64} || error("This regression is Float64 only")
    selected=collect(ms)
    first_grid=SW.build_grid(bg,first(selected);grid_kwargs...)
    c,p,g=context(bg,first(selected),first_grid;mass=mEMRI)
    # Appendix-C adaptive domains lie inside the fiducial interval. Explicit
    # edge overrides are shared by all modes and are included via first_grid.
    domain=(min(SW.r_isco(bg),first(first_grid.r)),max(5bg.aEMRI,last(first_grid.r)))
    modes=[SW.ModeContext(c.disc,p,m,0;closure=Adiabatic,domain).mode for m in selected]
    bank=regression_bank(c.disc,p,domain,modes)
    current=Ref((first_grid,g))
    cpu=Ref(cpu_seconds()); wall=Ref(time())
    builder=mode->begin
        grid=mode==first(modes) ? first_grid : SW.build_grid(bg,mode.m;grid_kwargs...)
        ctx,_,frozen=context(bg,mode.m,grid;mass=mEMRI)
        current[]=(grid,frozen)
        ctx
    end
    SW.solve_modes(builder,p,modes;bank,grid_builder=_->current[][2],landau,linsolve,
        boundary_builder=ctx->egd_boundary_rows(ctx;landau),
        transform=sol->begin
            record_mode(bg,current[][1],sol,cpu[],wall[])
            value=transform(sol)
            cpu[]=cpu_seconds(); wall[]=time()
            value
        end)

end

amf(bg, sol; kwargs...) = SW.amf(sol;kwargs...)
torque_density(bg, sol; kwargs...) = SW.torque_density(sol;kwargs...)
cumulative_torque(bg, sol; kwargs...) = SW.cumulative_torque(sol;kwargs...)
total_torque(bg, sol; kwargs...) = SW.total_torque(sol;kwargs...)
lindblad_torque(bg, sol; kwargs...) = SW.lindblad_torque(sol;kwargs...)

function torque_survey(bg, ms; threads=true, keep_profiles=false, verbose=false,
                       halfwidth=3bg.σε, model=:relativistic, grid_kwargs=(;), solve_kwargs=(;))
    # Fig. 5's Newtonian denominator is an explicitly separate reference model.
    model===:newtonian && return SW.torque_survey(bg,ms;threads=false,keep_profiles,
        verbose,halfwidth,model,grid_kwargs,solve_kwargs)
    model===:relativistic || error("Unknown model $model")
    wall=time(); mts=SW.ModeTorque{Float64}[]
    profiles=keep_profiles ? NamedTuple{(:r,:F_J,:density),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}[] : nothing
    batch_modes(bg,ms;grid_kwargs,solve_kwargs...) do sol
        t=time()-last(TRACE)["wall_seconds"]; m=sol.context.mode.m
        D=SW.torque_density(sol); integ=SW.cumulative_torque(sol;density=D)
        # The window depends on r_c/r only; retain its frozen EGD definition.
        window=findall(r->abs(r-sol.context.r_c)<=halfwidth,sol.r)
        isempty(window) && error("No node in the frozen corotation window")
        win=first(window):last(window)
        idx=[1,first(win),last(win),length(sol.r)]
        F=keep_profiles ? SW.amf(sol) : nothing
        edge=keep_profiles ? F[idx] :
            [SW.amf_mode(sol.context,sol.r[i],view(sol.y,:,i);hmodes=sol.source(sol.r[i]),
                         landau=sol.params.landau) for i in idx]
        push!(mts,SW.ModeTorque(Int(m),last(integ),edge[4]-edge[1],
            integ[idx[3]]-integ[idx[2]]-(edge[3]-edge[2]),edge[1],edge[4],
            length(sol.r),sol.residual,time()-t))
        keep_profiles && push!(profiles,(r=sol.r,F_J=F,density=D))
    end
    # Fixed-order sums, matching survey.jl, never a parallel reduction.
    T=0.0; TL=0.0
    for mt in mts
        T+=mt.T; TL+=mt.T_L
    end
    SW.TorqueSurvey(collect(Int,ms),mts,T,TL,profiles,time()-wall,SW.compute_budget())
end

function wake_modes(bg, ms; threads=true, verbose=false, grid_kwargs=(;),solve_kwargs=(;))
    batch_modes(bg,ms;grid_kwargs=merge((;extent=:fiducial),grid_kwargs),solve_kwargs...) do s
        (;m=s.context.mode.m,r=s.r,y=s.y,residual=s.residual,seconds=s.seconds)
    end
end

const BACKEND_CALLS=(:solve_mode,:torque_survey,:wake_modes,:amf,:torque_density,
                     :cumulative_torque,:total_torque,:lindblad_torque)

"Narrow call-site substitution; leaves the original numerical expressions intact."
function adapt(ex)
    ex isa Expr || return ex
    if ex.head===:call
        f=ex.args[1]
        f===:draw && return :(nothing)
        f===:cache_write && return :(nothing)
        f===:exit && return Expr(:call,GlobalRef(CircularRegression,:driver_exit),ex.args[2:end]...)
        f===:include && return Expr(:call,GlobalRef(CircularRegression,:include_driver),
                                    :(@__MODULE__),ex.args[2:end]...)
        if f isa Expr && f.head===:. && f.args[2] isa QuoteNode
            owner,name=f.args[1],f.args[2].value
            if owner===:SW && name in BACKEND_CALLS
                return Expr(:call,GlobalRef(CircularRegression,name),map(adapt,ex.args[2:end])...)
            elseif owner===:ReproHarness && name===:conclude!
                args=map(adapt,ex.args[2:end])
                push!(args[1].args,Expr(:kw,:harness,:ReproHarness))
                return Expr(:call,GlobalRef(CircularRegression,:conclude),args...)
            end
        end
    end
    # Do not adapt function signatures (e.g. `function draw(...)`).
    if ex.head===:function
        return Expr(:function,ex.args[1],adapt(ex.args[2]))
    elseif ex.head===:(=) && ex.args[1] isa Expr && ex.args[1].head===:call
        return Expr(:(=),ex.args[1],adapt(ex.args[2]))
    end
    Expr(ex.head,map(adapt,ex.args)...)
end

function include_driver(mod,path)
    path=realpath(path)
    startswith(path,PAPER*"/") ? Base.include(adapt,mod,path) : Base.include(mod,path)
end

struct DriverExit <: Exception
    code::Int
end
driver_exit(code)=throw(DriverExit(code))

function conclude(opts; harness, paper,id,script,rtol,atol=0.0,tolerances=Dict())
    opts.mode===:check || error("Golden regeneration is forbidden")
    vals=Dict(harness.results())
    golden=RH.expected_path(paper,id)
    _,expected=RH.read_expected(golden)
    rows=Dict{String,Any}[]
    for key in sort!(collect(union(keys(vals),keys(expected))))
        tol=get(tolerances,key,(;rtol,atol))
        present=haskey(vals,key) && haskey(expected,key)
        ok,err,note=present ? RH.compare_values(expected[key],vals[key];tol...) :
                             (false,NaN,"missing key")
        push!(rows,Dict("key"=>key,"rtol"=>tol.rtol,"atol"=>tol.atol,
            "value"=>RH.encode_value(get(vals,key,nothing)),
            "golden"=>RH.encode_value(get(expected,key,nothing)),"pass"=>ok,
            "relative_error"=>RH.encode_value(err),"note"=>note))
    end
    verdict=RH.compare_expected(golden,vals;rtol,atol,tolerances)
    doc=Dict("figure"=>id,"golden_sha256"=>bytes2hex(open(sha256,golden)),
             "driver_sha256"=>bytes2hex(open(sha256,script)),"keys"=>rows,
             "all_keys_pass"=>all(r->r["pass"],rows),"modes"=>length(TRACE))
    write(joinpath(OUTPUT[],id*".json"),RH.json_write(doc))
    all(r->r["pass"],rows) ? 0 : 1
end

function run_figure(id,out)
    id in [@sprintf("fig%02d",n) for n in 1:14] || error("Expected fig01..fig14")
    # No original-driver environment overrides: configurations are frozen.
    any(startswith(k,"GRTOOLS_FIG") for k in keys(ENV)) &&
        error("Remove GRTOOLS_FIG* overrides for the frozen regression")
    any(a->a in ("--regen","--plot-only"),ARGS) && error("Read-only regression only")
    Threads.nthreads()==1 || error("Use -t1")
    BLAS.set_num_threads(1)
    output=abspath(out)
    # Resolve symlinks before any write. The external parent must already
    # exist; each figure requires its own fresh destination.
    canonical=ispath(output) ? realpath(output) : joinpath(realpath(dirname(output)),basename(output))
    reproduction=realpath(joinpath(@__DIR__,".."))
    (canonical==reproduction || startswith(canonical,reproduction*"/")) &&
        error("Evidence must be outside the reproduction tree; frozen goldens are read only")
    isdir(output) && !isempty(readdir(output)) && error("Evidence directory must be empty for a new attempt")
    OUTPUT[]=output; mkpath(OUTPUT[]); empty!(TRACE); SINGLE_BANK[]=nothing
    original=Dict(p=>bytes2hex(open(sha256,p)) for p in
                  readdir(joinpath(PAPER,"expected");join=true))
    saved_args=copy(ARGS); empty!(ARGS); append!(ARGS,["--check","--tier","heavy"])
    wall=time(); cpu=cpu_seconds(); code=2; err=nothing
    driver=Module(gensym(:EGD))
    try
        Base.include(adapt,driver,joinpath(PAPER,id*".jl"))
    catch ex
        if ex isa LoadError && ex.error isa DriverExit
            code=ex.error.code
        elseif ex isa DriverExit
            code=ex.code
        else
            err=sprint(showerror,ex,catch_backtrace()); println(stderr,err)
        end
    finally
        empty!(ARGS); append!(ARGS,saved_args)
        if isdefined(driver,:ReproHarness)
            vals=Dict(k=>RH.encode_value(v) for (k,v) in Base.invokelatest(driver.ReproHarness.results))
            write(joinpath(OUTPUT[],id*"_registered_values.json"),RH.json_write(vals))
        end
        unchanged=all(bytes2hex(open(sha256,p))==h for (p,h) in original)
        status=Dict("figure"=>id,"exit_code"=>code,"error"=>err,
            "wall_seconds"=>time()-wall,"cpu_seconds"=>cpu_seconds()-cpu,
            "completed_modes"=>length(TRACE),"old_goldens_unmodified"=>unchanged)
        write(joinpath(OUTPUT[],id*"_status.json"),RH.json_write(status))
        unchanged || error("Frozen golden changed during run")
    end
    code
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    if length(ARGS)==2 && startswith(first(ARGS),"fig")
        exit(CircularRegression.run_figure(ARGS...))
    end
    opts=CircularRegression.RH.parse_args()
    opts.mode===:check || error("Golden regeneration is forbidden")
    isempty(opts.extra) || error("Usage: circular_regression.jl fig01 OUTPUT_DIRECTORY, or --check")
    CircularRegression.RH.require_tier(opts,:heavy)
    out=get(ENV,"ECC_REGRESSION_EVIDENCE","")
    isempty(out) && error("Set ECC_REGRESSION_EVIDENCE to the absolute evidence directory")
    codes=[CircularRegression.run_figure("fig"*lpad(string(i),2,'0'),
            joinpath(out,"fig"*lpad(string(i),2,'0'))) for i in 1:14]
    exit(maximum(codes))
end

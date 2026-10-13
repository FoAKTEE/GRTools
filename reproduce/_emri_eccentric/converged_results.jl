"""Signed-frequency continuation of first_results, with one shared bank and resumable shards.

Commands:
solve EVIDENCE BANKS CHI E M KMINUS KPLUS SHARD NSHARDS [target_M target_Kminus target_Kplus] [li]
collect|collect_with_axisymmetric EVIDENCE BANKS CHI E M KMINUS KPLUS [target_M target_Kminus target_Kplus] [li]
The preregistration and source-owner hashes are required for every resumed checkpoint.
"""
module ConvergedEccentricResults
include("first_results.jl")
using .FirstEccentricResults, AccretionDisc, Serialization, SHA, Printf
using NumericsCore: CubicSplineInterpolant
const FR=FirstEccentricResults
const SW=FR.SW
const DRIVER_SHA256=bytes2hex(open(sha256,@__FILE__))

"Native-grid gravitational totals; flux and wake reconstruction only on display radii."
function reduce_mode(sol,rq,wr)
    cp=FR.cpu_seconds()
    work=SW.gravitational_work(sol)
    td=getproperty.(work,:torque); pd=getproperty.(work,:power)
    torque=SW.total_torque(sol;density=td); power=SW.total_power(sol;density=pd)
    mi=[SW.mode_interp((;m=sol.context.mode.m,r=sol.r,y=sol.y);component=j) for j in 1:2]
    flux=[SW.amf_mode(sol.context,r,[mi[1](r),mi[2](r)];hmodes=sol.source(r),landau=sol.params.landau) for r in rq]
    profiles=hcat(flux,CubicSplineInterpolant(sol.r,td).(rq,1),CubicSplineInterpolant(sol.r,pd).(rq,1))
    density=[SW.reconstruct_fluid(sol.context,r,[mi[1](r),mi[2](r)];hmodes=sol.source(r),landau=sol.params.landau).deps/SW.eps_profile(sol.context.disc,r) for r in wr]
    mode=sol.context.mode
    record=Dict{String,Any}("m"=>mode.m,"k"=>mode.k,"omega"=>mode.omega,"torque_disc"=>torque,"power_disc"=>power,
        "nodes"=>length(sol.r),"residual"=>sol.residual,"bc_residual"=>sol.bc_residual,
        "bvp_wall_seconds"=>sol.seconds,"diagnostic_cpu_seconds"=>FR.cpu_seconds()-cp,
        "power_identity_relative"=>mode.m==0 ? 0.0 : abs(power-mode.omega/mode.m*torque)/max(abs(power),floatmin(Float64)))
    (;record,profiles,density)
end

function setup(evidence,banks,a,e,M,Kminus,Kplus,closure=Adiabatic)
    prereg=joinpath(evidence,"PREREG.json")
    isfile(prereg) || error("missing preregistration")
    id=@sprintf("a%+.1f_e%.1f_%s",a,e,closure==Adiabatic ? "ad" : "li")
    d,p,domain=FR.fixture(a,e,closure;epsilon_p=1e-6)
    for parent in (dirname(evidence),dirname(dirname(evidence)))
        previous=joinpath(parent,"gr-first-results",id,"result.json")
        if isfile(previous)
            domain=Tuple(Float64.(FR.json_read(read(previous,String))["domain"]))
            break
        end
    end
    mask=SW.eccentric_mode_mask(p.orbit;mmax=M,kmax=max(Kminus,Kplus),mask=:physical)
    selected=filter(s->s.m>0 && -Kminus<=s.k<=Kplus,mask.physical)
    axisymmetric_modes=filter(s->s.m==0,mask.physical)
    sec=Set((s.m,s.k) for s in mask.sectors.secular)
    dir=joinpath(banks,id); out=joinpath(evidence,id);mkpath(dir);mkpath(out)
    paths=filter(f->endswith(f,".jl"),readdir(joinpath(FR.REPO_ROOT,"capabilities/AccretionDisc/src/SpiralWaves");join=true))
    hashes=Dict(relpath(f,FR.REPO_ROOT)=>bytes2hex(open(sha256,f)) for f in paths)
    config=Dict("chi"=>a,"e"=>e,"closure"=>string(closure),"domain"=>collect(domain),"m_max"=>M,
        "k_minus"=>Kminus,"k_plus"=>Kplus,"Nq"=>e>=.5 && max(Kminus,Kplus)>20 ? 2048 : 1024,
        "Nphi"=>e>.3 ? 16384 : 8192,"source_owner_sha256"=>hashes,
        "prereg_sha256"=>bytes2hex(open(sha256,prereg)),
        "amendment_sha256"=>bytes2hex(open(sha256,joinpath(evidence,"STAGE1_AMENDMENT.json"))),
        "epsilon_p"=>1e-6)
    path=joinpath(dir,"configuration.json")
    if isfile(path)
        FR.json_read(read(path,String))==config || error("checkpoint provenance differs")
    else
        any(endswith(f,".jls") for f in readdir(dir)) && error("unversioned source/checkpoints cannot be reused")
        FR.writejson(path,config)
    end
    (;d,p,domain,selected,axisymmetric_modes,sec,dir,out,config,bankfile=joinpath(dir,"source_manifest.json"))
end

source_block(dir,m) = joinpath(dir,@sprintf("source_m%03d_%03d.jls",5*((m-1)÷5)+1,5*((m-1)÷5)+5))

"Rules borrowed from v4: exact-cell transfer remains uncertified at p=20,e=.2."
function radial_rule(c)
    m,k=c.mode.m,c.mode.k
    secular=(m,k)==(1,-1)
    n=abs(k)>20 ? 65536 : c.r_c===nothing ? 32768 : 16384
    cellclass=secular ? "secular" : m==0 ? "axisymmetric" : c.r_c===nothing ?
        (c.mode.omega<0 ? "negative_frequency_no_corotation" : "positive_frequency_no_corotation") :
        "corotation"
    template=secular ? "none" : m==0 ? "21" : abs(k)>20 ? "23/33 (high sideband)" :
        c.r_c===nothing ? "9 (no corotation)" : "32 (corotation)"
    (;n,cellclass,template,certified=false)
end

production_grid(c) = FR.production_grid(c;n=radial_rule(c).n)

function read_source_block(x,path,modes)
    manifest=FR.json_read(read(x.bankfile,String))
    bytes2hex(open(sha256,path))==manifest[basename(path)]["sha256"] || error("source block hash differs")
    FR.validate_bank(deserialize(path),x.p,modes,x.domain)
end

function solve_shard(x,M,Kminus,Kplus,shard,nshards)
    rq=collect(range(x.domain...;length=1600))
    wr=collect(range(x.domain[1],min(x.domain[2],2.5*20/(1-x.p.orbit.e));length=768))
    selected=filter(s->s.m<=M && -Kminus<=s.k<=Kplus && mod(s.m-1,nshards)==shard-1,x.selected)
    sort!(selected;by=s->(s.m,s.k))
    remaining=filter(s->!isfile(joinpath(x.dir,"mode_$(s.m)_$(s.k).jls")),selected)
    isempty(remaining) && return
    process_start=FR.cpu_seconds()
    budget=parse(Float64,get(ENV,"ECC_PROCESS_CPU_LIMIT","28000"))
    for path in unique(source_block(x.dir,s.m) for s in remaining)
    blockmodes=filter(s->source_block(x.dir,s.m)==path,remaining)
    bank=read_source_block(x,path,blockmodes);p=bank.perturber
    build=s->SW.ModeContext(x.d,p,s;closure=x.d.closure,domain=x.domain)
    started=Ref(0.);gridcost=Ref(0.)
    grid_builder=c->begin
        started[]=FR.cpu_seconds()
        try
            FR.cpu_seconds()-process_start < budget-300 || error("CPU budget checkpoint stop")
            grid=production_grid(c)
            gridcost[]=FR.cpu_seconds()-started[]
            grid
        catch err
            FR.writejson(joinpath(x.out,"failure_$(c.mode.m)_$(c.mode.k).json"),
                Dict("mode"=>[c.mode.m,c.mode.k],"error"=>sprint(showerror,err),"phase"=>"production_grid"))
            rethrow()
        end
    end
    SW.solve_modes(build,p,blockmodes;bank,grid_builder,transform=sol->begin
        cp=FR.cpu_seconds();mode=sol.context.mode
        bvp_cpu=cp-started[]-gridcost[]
        data=reduce_mode(sol,rq,wr)
        data.record["grid_cpu_seconds"]=gridcost[]
        data.record["bvp_boundary_cpu_seconds"]=bvp_cpu
        data.record["total_cpu_seconds"]=FR.cpu_seconds()-started[]
        data.record["driver_sha256"]=DRIVER_SHA256
        rule=radial_rule(sol.context)
        data.record["radial_class"]=rule.cellclass
        data.record["radial_template_cells"]=rule.template
        data.record["radial_certified"]=rule.certified
        data.record["radial_n"]=rule.n
        data.record["epsilon_p"]=sol.context.disc.εp
        data.record["conditional_secular"]=(mode.m,mode.k) in x.sec
        data.record["reduction_cpu_seconds"]=FR.cpu_seconds()-cp
        file=joinpath(x.dir,"mode_$(mode.m)_$(mode.k).jls")
        serialize(file*".tmp",data);mv(file*".tmp",file)
        println("MODE ",(mode.m,mode.k)," T=",data.record["torque_disc"]," P=",data.record["power_disc"]," nodes=",length(sol.r)," BVP_wall=",sol.seconds," diagnostic_cpu=",data.record["diagnostic_cpu_seconds"]);flush(stdout)
        nothing
    end)
    bank=nothing; GC.gc()
    end
end

function collect_case(x,M,Kminus,Kplus;include_axisymmetric=false)
    selected=filter(s->s.m<=M && -Kminus<=s.k<=Kplus,x.selected)
    include_axisymmetric && append!(selected,x.axisymmetric_modes)
    sort!(selected;by=s->(s.m,s.k))
    out=include_axisymmetric ? joinpath(x.out,"with_axisymmetric") : x.out
    bankdir=include_axisymmetric ? joinpath(x.dir,"with_axisymmetric") : x.dir
    mkpath(out);mkpath(bankdir)
    records=Dict[]; densities=NamedTuple[]
    rq=collect(range(x.domain...;length=1600))
    wr=collect(range(x.domain[1],min(x.domain[2],2.5*20/(1-x.p.orbit.e));length=768))
    phi=collect(range(-pi,pi;length=513));profiles=zeros(length(rq),3)
    axis_profiles=zeros(size(profiles))
    for mode in selected
        file=joinpath(x.dir,"mode_$(mode.m)_$(mode.k).jls")
        isfile(file) || error("incomplete rectangle: missing $file")
        data=deserialize(file);push!(records,data.record)
        if !data.record["conditional_secular"]
            profiles .+= data.profiles
            mode.m==0 && (axis_profiles .+= data.profiles)
            push!(densities,(context=(;mode),r=wr,y=vcat(permutedims(data.density),zeros(ComplexF64,1,length(wr)))))
        end
    end
    cost=FR.json_read(read(joinpath(x.out,"bank_cost.json"),String))
    meta=Dict("status"=>"candidate signed-frequency rectangle; convergence evaluated separately", "case"=>basename(x.out),
        "chi"=>x.p.orbit.a,"e"=>x.p.orbit.e,"p"=>20.,"closure"=>string(x.d.closure),"domain"=>collect(x.domain),
        "epsilon_ref"=>x.d.epsilon_ref,"mu"=>x.p.mass,"m_max"=>M,"k_max"=>max(Kminus,Kplus),
        "k_min"=>-Kminus,"k_plus"=>Kplus,"configuration"=>x.config,
        "Omega_phi"=>x.p.orbit.Ωφ,"Omega_r"=>x.p.orbit.Ωr,"axisymmetric_included"=>include_axisymmetric,
        "bank_cpu_seconds"=>cost["cpu_seconds"],"source_bank"=>x.bankfile,"source_bank_sha256"=>cost["sha256"],"wake_directory"=>bankdir)
    if include_axisymmetric
        axisbank=joinpath(x.dir,"source_axisymmetric.jls")
        meta["axisymmetric_k_max"]=maximum(s.k for s in x.axisymmetric_modes)
        meta["axisymmetric_source_bank"]=axisbank
        meta["axisymmetric_source_bank_sha256"]=bytes2hex(open(sha256,axisbank))
        meta["status"]="candidate finite nonsecular union; axisymmetric boundary HOLE; no convergence admission"
    end
    lm=sort(unique(vcat(collect(5:5:min(M,20)),collect(30:10:M),M)))
    lk=sort(unique(vcat([k for k in (3,6,12,20,30,40,50,60,70,80) if k<=Kminus],Kminus)))
    FR.write_products(out,bankdir,meta,records,rq,wr,phi,profiles,densities,lm,lk)
    if include_axisymmetric
        ri=round.(Int,range(1,length(rq);length=5))
        wi=round.(Int,range(1,length(wr);length=5))
        phi_i=round.(Int,range(1,length(phi);length=5))
        axisdensities=filter(s->s.context.mode.m==0,densities)
        proof=Dict{String,Any}("profile_indices_zero_based"=>ri.-1,
            "wake_indices_zero_based"=>wi.-1,"phi_indices_zero_based"=>phi_i.-1,
            "axisymmetric_profiles"=>axis_profiles[ri,:])
        for phase in ("peri","apo")
            proof["axisymmetric_wake_"*phase]=SW.synthesize_eccentric(
                axisdensities,wr[wi],phi[phi_i],meta["wake"]["t_"*phase])
        end
        FR.writejson(joinpath(out,"AXISYMMETRIC_SUPERPOSITION.json"),proof)
    end
    println("COLLECTED ",length(selected)," modes")
end

function main(args)
    closure=last(args)=="li" ? LocallyIsothermalOptionA : Adiabatic
    closure==LocallyIsothermalOptionA && (args=args[1:end-1])
    kind,evidence,banks=args[1:3]
    a,e=parse.(Float64,args[4:5]);M,Kminus,Kplus=parse.(Int,args[6:8])
    snapshots=joinpath(evidence,"driver_versions");mkpath(snapshots)
    cp(@__FILE__,joinpath(snapshots,DRIVER_SHA256*".jl");force=true)
    x=setup(evidence,banks,a,e,M,Kminus,Kplus,closure)
    if kind=="solve"
        shard,nshards=parse.(Int,args[9:10])
        bounds=length(args)>=13 ? parse.(Int,args[11:13]) : [M,Kminus,Kplus]
        solve_shard(x,bounds...,shard,nshards)
    elseif kind in ("collect","collect_with_axisymmetric")
        bounds=length(args)>=11 ? parse.(Int,args[9:11]) : [M,Kminus,Kplus]
        collect_case(x,bounds...;include_axisymmetric=kind=="collect_with_axisymmetric")
    else
        error("unknown command")
    end
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    ConvergedEccentricResults.main(ARGS)
end

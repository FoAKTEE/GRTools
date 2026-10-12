"""PRELIMINARY eccentric Kerr disc sums. Source banks and response checkpoints stay outside the repo.

Usage: first_results.jl pilot|case EVIDENCE BANKS CHI E MMAX KMAX [adiabatic|li]
Requires PREREG.json in EVIDENCE. All production arithmetic is Float64.
"""
module FirstEccentricResults
using AccretionDisc, Spacetime, LinearAlgebra, Printf, Serialization, SHA
using AccretionDisc.ForcingModes
using NumericsCore: CubicSplineInterpolant
const SW = AccretionDisc.SpiralWaves
const DRIVER_SHA256 = bytes2hex(open(sha256,@__FILE__))
const REPO_ROOT = normpath(joinpath(@__DIR__,"../.."))
const BASE_COMMIT = readchomp(`git -C $REPO_ROOT rev-parse HEAD`)
include(joinpath(@__DIR__,"../_lib/repro_harness.jl"))
using .ReproHarness: json_write, json_read, encode_value
cpu_seconds() = Float64(ccall(:clock,Clong,()))/1e6
writejson(path,x) = write(path,json_write(encode_value(x)))

function modal_sum(records,key)
    value,correction=0.0,0.0
    for row in records
        value,correction=SW._kbn_add(value,correction,row[key])
    end
    value+correction
end

function fixture(a,e,closure)
    p = 20.0
    orbit = EquatorialOrbit(a,p,e) # owner checks separatrix and bound-orbit validity
    d = SW.DiscModel(;χ=a,r_ref=p/(1-e^2),h_ref=0.06,pSigma=1.5,qTemp=0.0,
                    epsilon_ref=1e-8,closure,evolution=:fixed,εp=1e-5)
    pert = SW.EccentricPerturber(orbit;mass=1e-5,softening=FixedSoftening(0.018d.r_ref))
    domain = (1.001SW.r_isco(d),8p/(1-e))
    domain[1] < p/(1+e) < p/(1-e) <= domain[2] ||
        (e==0 && domain[1]<p<domain[2]) || error("orbit outside production domain")
    d,pert,domain
end

function production_grid(c)
    # Candidate transfer from gr-radial-convergence/rerun/production_rules.csv:
    # grid-only failures use their measured next doubling; no certification implied.
    n = abs(c.mode.k)>20 ? 65536 : 8192
    width = SW.corotation_width(c)
    SW.build_grid(c;nL=n÷64,nR=n÷64,nC=n+1,n_anc=512,
                  σ_wc_frac=width===nothing ? 0.01 : 16width/c.r_c)
end

function validate_bank(bank,pert,modes,domain)
    b=bank.perturber
    (b.orbit.a,b.orbit.p,b.orbit.e,b.mass,b.softening)==
        (pert.orbit.a,pert.orbit.p,pert.orbit.e,pert.mass,pert.softening) ||
        error("source bank perturber differs")
    all(mode in bank.modes for mode in modes) || error("source bank lacks requested modes/frequencies")
    first(bank.grid.edges)<=domain[1]<domain[2]<=last(bank.grid.edges) ||
        error("source bank does not cover response domain")
    bank
end

function summarize_mode(sol,rq,wr)
    t0,cp0 = time(),cpu_seconds()
    diagnostic = SW.balance_diagnostics(sol)
    # Display profiles do not replace the full-grid quadrature used for totals.
    profiles = hcat([CubicSplineInterpolant(sol.r,v).(rq,1)
                    for v in (diagnostic.amf,diagnostic.torque_density,diagnostic.power_density)]...)
    # Interpolate the master state once, then use the certified closure reconstruction
    # at each image radius. No near-corotation downsampling enters the totals.
    mi = [SW.mode_interp((;m=sol.context.mode.m,r=sol.r,y=sol.y);component=j) for j in 1:2]
    density = [SW.reconstruct_fluid(sol.context,r,[mi[1](r),mi[2](r)];
        hmodes=sol.source(r),landau=sol.params.landau).deps/SW.eps_profile(sol.context.disc,r) for r in wr]
    mode = sol.context.mode
    record = Dict("m"=>mode.m,"k"=>mode.k,"omega"=>mode.omega,
        "torque_disc"=>diagnostic.torque,"power_disc"=>diagnostic.power,
        "conditional_secular"=>mode.m==1 && mode.k==-1,"closure"=>string(sol.closure),
        "emitted"=>diagnostic.emitted,"nodes"=>length(sol.r),"residual"=>sol.residual,
        "bc_residual"=>sol.bc_residual,"bvp_wall_seconds"=>sol.seconds,
        "diagnostic_wall_seconds"=>time()-t0,"diagnostic_cpu_seconds"=>cpu_seconds()-cp0,
        "power_identity_relative"=>abs(diagnostic.power-mode.omega/mode.m*diagnostic.torque)/
                                   max(abs(diagnostic.power),floatmin(Float64)))
    (;record,profiles,density)
end

function rates(orbit,torque,power)
    Etotal,Ltotal = -power,-torque
    if iszero(orbit.e)
        return Dict("Edot_orbit"=>Etotal,"Ldot_orbit"=>Ltotal,"Edot_specific"=>Etotal/1e-5,
                    "Ldot_specific"=>Ltotal/1e-5,"pdot"=>nothing,"edot"=>nothing,
                    "circular_map"=>"undefined: owner rejects degenerate e=0 Jacobian")
    end
    J = Spacetime.EL_to_pe_jacobian(orbit.a,orbit.p,orbit.e)
    specific = [Etotal,Ltotal]./1e-5
    pe = J*specific
    guarded = Spacetime.ELdot_to_pedot(orbit,specific...)
    isapprox(pe[1],guarded.pdot;rtol=1e-13) && isapprox(pe[2],guarded.edot;rtol=1e-13) ||
        error("Jacobian and guarded rate mapping disagree")
    Dict("Edot_orbit"=>Etotal,"Ldot_orbit"=>Ltotal,"Edot_specific"=>specific[1],
         "Ldot_specific"=>specific[2],"pdot"=>pe[1],"edot"=>pe[2],
         "fraction_per_radial_period"=>guarded.fraction_per_radial_period,
         "jacobian"=>J,"jacobian_condition_inf"=>opnorm(J,Inf)*opnorm(inv(J),Inf))
end

function write_products(out,bankdir,meta,records,rq,wr,phi,profiles,densities,ladder_m,ladder_k)
    orbit = EquatorialOrbit(meta["chi"],20.0,meta["e"])
    rows = Dict[]
    for M in ladder_m, K in ladder_k
        subset = [r for r in records if r["m"]<=M && abs(r["k"])<=K && !r["conditional_secular"]]
        # Reuse the source owner's compensated accumulator in fixed (m,k) order.
        torque,power = modal_sum(subset,"torque_disc"),modal_sum(subset,"power_disc")
        push!(rows,merge(Dict("m_max"=>M,"k_max"=>K,"modes"=>length(subset),
            "torque_disc"=>torque,"power_disc"=>power,
            "sum_abs_torque"=>sum(abs(r["torque_disc"]) for r in subset),
            "sum_abs_power"=>sum(abs(r["power_disc"]) for r in subset)),rates(orbit,torque,power)))
    end
    sec = filter(r->r["conditional_secular"],records)
    meta["secular"] = [merge(r,rates(orbit,r["torque_disc"],r["power_disc"])) for r in sec]
    meta["ladder"] = rows
    # Detailed per-mode diagnostics are in the external checkpoint bank and CSV;
    # avoid a multi-megabyte JSON duplicate in the evidence repository.
    meta["max_residual"] = maximum(r["residual"] for r in records)
    meta["max_bc_residual"] = maximum(r["bc_residual"] for r in records)
    meta["max_power_identity_relative"] = maximum(r["power_identity_relative"] for r in records)
    writejson(joinpath(out,"result.json"),meta)
    open(joinpath(out,"ladder.csv"),"w") do io
        keys = ["m_max","k_max","modes","torque_disc","power_disc","Edot_orbit","Ldot_orbit",
                "Edot_specific","Ldot_specific","pdot","edot","sum_abs_torque","sum_abs_power"]
        println(io,join(keys,','))
        for row in rows
            println(io,join([row[k]===nothing ? "" : string(row[k]) for k in keys],','))
        end
    end
    open(joinpath(out,"modes.csv"),"w") do io
        keys=["m","k","omega","torque_disc","power_disc","conditional_secular","nodes","residual","bc_residual"]
        println(io,join(keys,','))
        for row in records
            println(io,join([string(row[k]) for k in keys],','))
        end
    end
    open(joinpath(out,"profiles.csv"),"w") do io
        println(io,"r,amf_disc,torque_density_disc,power_density_disc")
        for i in eachindex(rq)
            println(io,join([rq[i];profiles[i,:]],','))
        end
    end
    # At radial phase pi: retain the trajectory's pericentre-offset convention.
    t_apo = pi/orbit.Ωr+orbit.trajectory(pi).Δtr-orbit.trajectory(0.0).Δtr
    phi_apo = pi*orbit.Υφ/orbit.Υr+orbit.trajectory(pi).Δφr-orbit.trajectory(0.0).Δφr
    meta["wake"] = Dict("field"=>"Eulerian delta epsilon / epsilon; not rest-mass surface density",
        "t_peri"=>0.0,"t_apo"=>t_apo,"phi_apo"=>phi_apo,
        "r_peri"=>20/(1+orbit.e),"r_apo"=>20/(1-orbit.e),"conditional_secular_included"=>false)
    for (name,t) in (("peri",0.0),("apo",t_apo))
        F=SW.synthesize_eccentric(densities,wr,phi,t)
        # Binary matrix and coordinate vectors are retained outside the repository.
        open(joinpath(bankdir,"wake_$name.bin"),"w") do io
            write(io,Int64(length(wr)),Int64(length(phi)),wr,phi,F)
        end
    end
    writejson(joinpath(out,"result.json"),meta)
end

function run_case(kind,evidence,banks,a,e,M,K,closure)
    isfile(joinpath(evidence,"PREREG.json")) || error("missing preregistration")
    BLAS.set_num_threads(1)
    id = @sprintf("a%+.1f_e%.1f_%s",a,e,closure==Adiabatic ? "ad" : "li")
    kind=="pilot" && (id *= "_pilot")
    tstart,cpstart = time(),cpu_seconds()
    d,pert,domain = fixture(a,e,closure)
    original_domain=domain
    amendment_path=joinpath(evidence,"DOMAIN_AMENDMENT.json")
    amendment=isfile(amendment_path) && kind=="case" ?
        get(json_read(read(amendment_path,String))["cases"],id,nothing) : nothing
    if amendment!==nothing
        domain=Tuple(Float64.(amendment["domain"]))
        original_domain[1]<=domain[1]<20/(1+e)<=20/(1-e)<domain[2]<=original_domain[2] ||
            error("amended domain must contain orbit and be contained in original domain")
    end
    out=joinpath(evidence,id)
    bankdir=joinpath(banks,id*(amendment===nothing ? "" : "_domain_v2"))
    mkpath(out); mkpath(bankdir)
    nq=e>=0.5 && K>20 ? 2048 : 1024
    nphi=e>.3 ? 16384 : 8192
    owner_paths=[joinpath(REPO_ROOT,"capabilities/AccretionDisc/src/SpiralWaves",f)
                 for f in ("eccentric_source.jl","eccentric_modes.jl","eccentric.jl","background.jl",
                           "diagnostics.jl","balance_diagnostics.jl","master.jl","grid.jl",
                           "eccentric_synthesis.jl")]
    push!(owner_paths,joinpath(REPO_ROOT,"capabilities/Spacetime/src/Geodesics/equatorial_map.jl"))
    configuration=Dict("driver_sha256"=>DRIVER_SHA256,"base_commit"=>BASE_COMMIT,
        "owner_sha256"=>Dict(relpath(path,REPO_ROOT)=>bytes2hex(open(sha256,path)) for path in owner_paths),
        "chi"=>a,"e"=>e,"m_max"=>M,
        "k_max"=>K,"closure"=>string(closure),"kind"=>kind,"Nq"=>nq,"Nphi_floor"=>nphi,
        "prereg_sha256"=>bytes2hex(open(sha256,joinpath(evidence,"PREREG.json"))))
    configuration["domain"]=collect(domain)
    if amendment!==nothing
        configuration["domain_amendment_sha256"]=bytes2hex(open(sha256,amendment_path))
        source_path=amendment["source_bank"]
        source_config_path=joinpath(dirname(source_path),"configuration.json")
        bytes2hex(open(sha256,source_config_path))==amendment["source_configuration_sha256"] ||
            error("import source configuration hash changed")
        origin=json_read(read(source_config_path,String))
        for key in ("base_commit","chi","e","m_max","k_max","Nq","Nphi_floor","prereg_sha256")
            origin[key]==configuration[key] || error("import source provenance differs: $key")
        end
        for (path,hash) in origin["owner_sha256"]
            configuration["owner_sha256"][path]==hash || error("import source owner changed: $path")
        end
        configuration["source_import"]=Dict("bank"=>source_path,
            "bank_sha256"=>bytes2hex(open(sha256,source_path)),
            "configuration"=>source_config_path,
            "configuration_sha256"=>amendment["source_configuration_sha256"])
    end
    configuration_path=joinpath(bankdir,"configuration.json")
    if isfile(configuration_path)
        json_read(read(configuration_path,String))==configuration ||
            error("checkpoint configuration/code differs; preserve old evidence and use a new bank directory")
    else
        any(endswith(f,".jls") for f in readdir(bankdir)) &&
            error("unversioned source/state/checkpoint cannot be reused")
        writejson(configuration_path,configuration)
    end
    candidates=[HarmonicMode(m,k,orbit_frequency(pert.orbit,m,k)) for m in 1:M for k in -K:K]
    selected=filter(m->m.omega>0,candidates)
    excluded=count(m->m.omega<=0,candidates)
    if kind=="pilot"
        pairs = e==0 ? [(1,0),(M,0)] : [(1,-1),(2,1),(M,K),(M,-K)]
        selected = filter(m->(m.m,m.k) in pairs,selected)
    end
    # Check every requested boundary before any expensive bank construction.
    for mode in selected
        SW.mode_boundary_rows(SW.ModeContext(d,pert,mode;closure,domain))
    end
    println("VALID ",id," domain=",domain," separatrix=",separatrix(a,e,1.0)," modes=",length(selected)); flush(stdout)
    bankpath=joinpath(bankdir,"source.jls")
    tb,cb=time(),cpu_seconds()
    if !isfile(bankpath) && amendment!==nothing
        cp(amendment["source_bank"],bankpath)
    end
    if isfile(bankpath)
        bank=validate_bank(deserialize(bankpath),pert,selected,domain)
        if amendment!==nothing
            bytes2hex(open(sha256,bankpath))==configuration["source_import"]["bank_sha256"] ||
                error("copied source bank hash differs")
        end
        pert=bank.perturber
    else
        bank=SW.source_radial_bank(d,pert,domain,selected;Nq=nq,Nphi=nphi,
                                  method=:fft,max_memory_bytes=8*1024^3)
        serialize(bankpath,bank)
    end
    bankcost=cpu_seconds()-cb
    println("BANK ",id," radial_nodes=",length(bank.grid.rs)," cpu_seconds=",bankcost," wall_seconds=",time()-tb); flush(stdout)
    rq=collect(range(domain...;length=1600))
    wr=collect(range(domain[1],min(domain[2],2.5*20/(1-e));length=768))
    phi=collect(range(-pi,pi;length=513))
    profiles=zeros(length(rq),3)
    densities=NamedTuple[]; records=Dict[]
    build=mode->SW.ModeContext(d,pert,mode;closure,domain)
    for mode in selected
        checkpoint=joinpath(bankdir,"mode_$(mode.m)_$(mode.k).jls")
        if isfile(checkpoint)
            data=deserialize(checkpoint)
        else
            ts,cs=time(),cpu_seconds()
            data=only(SW.solve_modes(build,pert,[mode];bank,grid_builder=production_grid,
                transform=s->begin
                    # Keep full solved states externally for independent reintegration;
                    # do not serialize the shared source bank again for every mode.
                    statepath=joinpath(bankdir,"state_$(mode.m)_$(mode.k).jls")
                    serialize(statepath,(;r=s.r,y=s.y,context=s.context,params=s.params))
                    reduced=summarize_mode(s,rq,wr)
                    reduced.record["state_file"]=statepath
                    reduced.record["state_sha256"]=bytes2hex(open(sha256,statepath))
                    reduced
                end))
            data.record["total_cpu_seconds"]=cpu_seconds()-cs
            data.record["total_wall_seconds"]=time()-ts
            serialize(checkpoint,data)
        end
        push!(records,data.record)
        if !data.record["conditional_secular"]
            profiles .+= data.profiles
            push!(densities,(context=(;mode),r=wr,y=vcat(permutedims(data.density),zeros(ComplexF64,1,length(wr)))))
        end
        println("MODE ",id," ",(mode.m,mode.k)," T=",data.record["torque_disc"],
                " cpu=",data.record["total_cpu_seconds"]," nodes=",data.record["nodes"]); flush(stdout)
        writejson(joinpath(out,"progress.json"),Dict("case"=>id,"finished"=>length(records),
            "planned"=>length(selected),"bank_cpu_seconds"=>bankcost,"cpu_seconds"=>cpu_seconds()-cpstart,
            "last_mode"=>data.record))
        GC.gc()
    end
    meta=Dict("status"=>"PRELIMINARY, sector restricted; not convergence certified","case"=>id,
        "configuration"=>configuration,
        "chi"=>a,"p"=>20.0,"e"=>e,"closure"=>string(closure),"domain"=>collect(domain),
        "original_domain"=>collect(original_domain),"domain_amendment"=>amendment,
        "r_ref"=>d.r_ref,"mu"=>pert.mass,"epsilon_ref"=>d.epsilon_ref,
        "m_max"=>M,"k_max"=>K,"excluded_static_count"=>count(m->iszero(m.omega),candidates),
        "excluded_retrograde_count"=>count(m->m.omega<0,candidates),"axisymmetric_included"=>false,
        "bank_cpu_seconds"=>bankcost,"bank_radial_nodes"=>length(bank.grid.rs),
        "source_bank"=>bankpath,"source_bank_sha256"=>bytes2hex(open(sha256,bankpath)),
        "Omega_phi"=>pert.orbit.Ωφ,"Omega_r"=>pert.orbit.Ωr,
        "cpu_seconds"=>cpu_seconds()-cpstart,"wall_seconds"=>time()-tstart,
        "separatrix"=>separatrix(a,e,1.0),"prereg_sha256"=>bytes2hex(open(sha256,joinpath(evidence,"PREREG.json"))))
    if kind=="pilot"
        meta["modal_records"]=records
        writejson(joinpath(out,"result.json"),meta)
    else
        # Smaller supplementary ladder entries preserve a measured trend when the
        # full requested ladder cannot fit the preregistered CPU cap.
        lm=sort(unique([m for m in (5,10,20,40,60,M) if m<=M]))
        lk=e==0 ? [0] : sort(unique([k for k in (e<.3 ? (3,6,12,20,K) : (5,10,20,30,K)) if k<=K]))
        write_products(out,bankdir,meta,records,rq,wr,phi,profiles,densities,lm,lk)
    end
    println("FINISHED ",id," cpu_seconds=",cpu_seconds()-cpstart); flush(stdout)
end

function main(args)
    length(args)>=7 || error("usage: pilot|case EVIDENCE BANKS CHI E MMAX KMAX [adiabatic|li]")
    kind,evidence,banks=args[1:3]
    kind in ("pilot","case") || error("unknown run kind")
    closure=length(args)>7 && args[8]=="li" ? LocallyIsothermalOptionA : Adiabatic
    run_case(kind,evidence,banks,parse(Float64,args[4]),parse(Float64,args[5]),parse(Int,args[6]),parse(Int,args[7]),closure)
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    FirstEccentricResults.main(ARGS)
end

# GR signed-frequency probes and a reduced FR24-mask weak-field inventory.
# Uses the existing source-adapted radial bank and the existing mode solver.
module GRExcludedSector
using AccretionDisc, Spacetime, LinearAlgebra, TOML, Printf, Serialization
using AccretionDisc.ForcingModes
const SW=AccretionDisc.SpiralWaves

function critical_balance(sol,b,flux)
    c=sol.context
    # Outside six Gaussian Landau widths the shift is at most exp(-18) of peak.
    width=6c.disc.σε
    win=findall(r->abs(r-c.r_c)<width,sol.r)
    lo,hi=first(win),last(win)
    cumulative=SW.cumulative_torque(sol;density=b.torque_density)
    deposit=(cumulative[hi]-cumulative[lo])-(flux[hi]-flux[lo])
    residual=abs(last(flux)-first(flux)+deposit-b.torque)/abs(b.torque)
    (;width,deposit,residual)
end

function solve_inventory(a,p,e,modes;Nq=256,n=8001,q=0.0,soft=0.03,domain=nothing,
                         bank_path=nothing,source_bank=nothing,
                         source_memory_bytes=512*1024^2,checkpoint_path=nothing,
                         grid_strategy=:adaptive,reuse_rows=Dict{String,Any}[],core_nodes=257)
    grid_strategy in (:adaptive,:core_log) || error("unknown grid strategy")
    aref=p/(1-e^2)
    d=SW.DiscModel(;χ=a,r_ref=aref,h_ref=0.06,pSigma=0.5,qTemp=q,
                   closure=LocallyIsothermalOptionA,evolution=:fixed,εp=1e-4)
    pert=SW.EccentricPerturber(EquatorialOrbit(a,p,e);mass=1e-5,softening=FixedSoftening(soft*aref))
    inventory=SW.eccentric_mode_mask(pert.orbit;mask=:physical,
        mmax=maximum(s.m for s in modes),kmax=maximum(abs(s.k) for s in modes))
    sector_by_pair=Dict((s.m,s.k)=>String(name) for (name,sector) in pairs(inventory.sectors) for s in sector)
    dom=something(domain,(max(1.02SW.r_isco(d),0.05aref),5aref))
    bank=source_bank===nothing ? SW.source_radial_bank(d,pert,dom,modes;Nq,max_memory_bytes=source_memory_bytes) : source_bank
    bank.modes==modes && bank.perturber.orbit.a==a && bank.perturber.orbit.p==p && bank.perturber.orbit.e==e ||
        error("source-bank replay parameters do not match")
    bank_path===nothing || serialize(bank_path,bank)
    println("BANK a=",a," p=",p," e=",e," modes=",length(modes)," source_nodes=",length(bank.grid.rs));flush(stdout)
    rows=deepcopy(reuse_rows)
    mode_by_pair=Dict((s.m,s.k)=>s for s in modes)
    reused=Set((r["m"],r["k"]) for r in rows)
    length(reused)==length(rows) || error("duplicate reused modes")
    for r in rows
        s=get(mode_by_pair,(r["m"],r["k"]),nothing)
        s!==nothing && r["omega"]==s.omega && r["a"]==a && r["p"]==p &&
            r["e"]==e && r["Nq"]==Nq && r["qTemp"]==q || error("reused response parameters differ")
        get!(r,"grid_strategy","adaptive")
    end
    for s in modes
        (s.m,s.k) in reused && continue
        c=SW.ModeContext(d,pert,s;closure=d.closure,domain=dom)
        source=SW.mode_source(bank,s)
        # The core-only option retains the same existing corotation block on
        # the fixed log mesh. It makes no claim of far-field phase convergence.
        extra = if grid_strategy===:adaptive
            SW.build_grid(c;nL=2,nR=2,nC=core_nodes,n_anc=256).r
        elseif c.r_c===nothing
            Float64[]
        else
            filter(r->dom[1]<=r<=dom[2],SW.corotation_block(c.r_c,.001c.r_c,4.0,core_nodes))
        end
        candidates=sort!(vcat(exp.(range(log(dom[1]),log(dom[2]);length=n)),extra))
        r=[dom[1]]
        for t in candidates
            t-last(r)>1e-12*max(abs(t),abs(last(r))) && dom[2]-t>1e-12*dom[2] && push!(r,t)
        end
        push!(r,dom[2])
        grid=SW.EccentricGrid(r,c.r_c,1:0,SW.mode_cache_key(c),(;sector_probe=true))
        sol=SW.solve_mode(c;source,grid,boundary_rows=SW.mode_boundary_rows(c))
        b=SW.balance_diagnostics(sol)
        # Common force and time normalizations cancel out of sector fractions.
        row=Dict{String,Any}("a"=>a,"p"=>p,"e"=>e,"m"=>s.m,"k"=>s.k,"omega"=>s.omega,
            "qTemp"=>q,"Nq"=>Nq,"nodes"=>length(r),"source_nodes"=>length(bank.grid.rs),
            "torque"=>b.torque,"power"=>b.power,"fin"=>first(b.amf),"fout"=>last(b.amf),
            "residual"=>sol.residual,"bc_residual"=>sol.bc_residual,
            "sector"=>sector_by_pair[(s.m,s.k)],"grid_strategy"=>String(grid_strategy),
            "base_nodes"=>n,"core_nodes"=>core_nodes)
        if s.m>0
            row["raw_flux_jump_relative"]=abs((last(b.amf)-first(b.amf)-b.torque)/b.torque)
            row["power_identity_relative"]=abs(b.power-s.omega/s.m*b.torque)/max(abs(b.power),eps(Float64)*abs(b.torque))
        else
            row["axisymmetric_torque_absolute"]=abs(b.torque)
        end
        if iszero(q)
            ad=SW.DiscModel(d.reference;closure=Adiabatic,evolution=:fixed)
            ca=SW.ModeContext(ad,s;closure=Adiabatic,domain=dom)
            fa=[SW.amf_mode(ca,r[i],sol.y[:,i];hmodes=source(r[i]),landau=sol.params.landau) for i in eachindex(r)]
            row["adiabatic_flux_jump_relative"]=abs((last(fa)-first(fa)-b.torque)/b.torque)
            if c.r_c !== nothing
                critical=critical_balance(sol,b,fa)
                row["critical_layer_torque"]=critical.deposit
                row["balance_outside_corotation_relative"]=critical.residual
                row["critical_halfwidth"]=critical.width
            end
            ga=SW.EccentricGrid(r,ca.r_c,1:0,SW.mode_cache_key(ca),(;sector_probe=true))
            sa=SW.solve_mode(ca;source,grid=ga,boundary_rows=SW.mode_boundary_rows(ca))
            ba=SW.balance_diagnostics(sa)
            row["adiabatic_solve_torque"]=ba.torque
            row["adiabatic_solve_fin"]=first(ba.amf)
            row["adiabatic_solve_fout"]=last(ba.amf)
            row["adiabatic_solve_residual"]=sa.residual
            row["adiabatic_solve_balance"]=abs(last(ba.amf)-first(ba.amf)-ba.torque)/abs(ba.torque)
            if ca.r_c !== nothing
                critical=critical_balance(sa,ba,ba.amf)
                row["adiabatic_solve_critical_layer_torque"]=critical.deposit
                row["adiabatic_solve_balance_outside_corotation"]=critical.residual
            end
        end
        push!(rows,row)
        if s.omega<0 || s.k == -s.m
            println("GRMODE ",row);flush(stdout)
        end
        if checkpoint_path!==nothing && s.k==inventory.kmax
            open(io->TOML.print(io,Dict("modes"=>rows,"kind"=>"partial","completed_m"=>s.m)),checkpoint_path,"w")
            println("GRFAMILY m=",s.m," completed_modes=",length(rows),"/",length(modes));flush(stdout)
        end
    end
    rows
end

function run(outpath,kind;Nq=256,n=8001,mmax=4,bank_prefix=nothing,weak_p=1000.0,
             grid_strategy=:adaptive)
    rows=Dict{String,Any}[]
    if kind in ("physical","balance","replay")
        cases=kind=="balance" ? ((0.0,30.0,0.3),) : ((0.0,30.0,0.3),(0.0,30.0,0.6),(0.9,10.0,0.3),(-0.9,10.0,0.3))
        for (a,p,e) in cases
            o=EquatorialOrbit(a,p,e)
            # -8 also gives negative omega for the fast-precessing retrograde-spin orbit.
            modes=[HarmonicMode(2,k,Spacetime.orbit_frequency(o,2,k)) for k in (-8,-4,-2)]
            bank=kind=="replay" ? deserialize(bank_prefix*".$a.$p.$e.bank.jls") : nothing
            append!(rows,solve_inventory(a,p,e,modes;Nq,n,source_bank=bank,
                bank_path=kind=="replay" ? nothing : outpath*".$a.$p.$e.bank.jls"))
            open(io->TOML.print(io,Dict("modes"=>rows,"kind"=>kind)),outpath,"w")
        end
    elseif kind in ("weak","weak-replay","weak-full","weak-full-replay","weak-full-resume","axisymmetric","axisymmetric-replay")
        o=EquatorialOrbit(0.0,weak_p,0.12)
        mask=SW.eccentric_mode_mask(o;policy=:fr24,mask=:physical,mmax,kmax=40)
        axisymmetric=startswith(kind,"axisymmetric")
        full=startswith(kind,"weak-full")
        resume=endswith(kind,"-resume")
        replay=endswith(kind,"-replay") || resume
        modes=axisymmetric ? mask.sectors.axisymmetric : filter(s->s.m>0,mask.physical)
        bank=replay ? deserialize(bank_prefix*".bank.jls") : nothing
        reuse_rows=resume ? TOML.parsefile(bank_prefix*".partial.toml")["modes"] : Dict{String,Any}[]
        append!(rows,solve_inventory(o.a,o.p,o.e,modes;Nq,n,q=0.75,soft=0.018,
            source_bank=bank,bank_path=replay ? nothing : outpath*".bank.jls",grid_strategy,reuse_rows,
            source_memory_bytes=full ? 8*1024^3 : 512*1024^2,
            checkpoint_path=full ? outpath*".partial.toml" : nothing))
        if axisymmetric
            totals=Dict("torque"=>sum(r["torque"] for r in rows),"power"=>sum(r["power"] for r in rows),"n"=>n,"Nq"=>Nq)
            open(io->TOML.print(io,Dict("modes"=>rows,"totals"=>totals,"kind"=>kind)),outpath,"w")
            println("AXISYMMETRIC_TOTAL ",totals)
            return
        end
        totals=Dict{String,Any}("mmax"=>mmax,"kmax"=>40,"p"=>o.p,"e"=>o.e,"Nq"=>Nq,"n"=>n)
        for key in ("torque","power")
            for sector in ("paper","excluded","secular")
                totals[sector*"_"*key]=sum(r[key] for r in rows if r["sector"]==sector)
            end
            total=sum(r[key] for r in rows)
            totals["total_"*key]=total
            nonsecular=totals["paper_"*key]+totals["excluded_"*key]
            totals["nonsecular_total_"*key]=nonsecular
            totals["excluded_share_nonsecular_"*key]=totals["excluded_"*key]/nonsecular
            totals["excluded_share_"*key]=totals["excluded_"*key]/total
            totals["excluded_plus_secular_share_"*key]=(totals["excluded_"*key]+totals["secular_"*key])/total
        end
        open(outpath,"w") do io
            TOML.print(io,Dict("modes"=>rows,"totals"=>totals,"kind"=>kind,
                "scope"=>"m=1:$mmax, |k|<=40 with FR24 integer mask and FR24 q=.75,pSigma=.5,h=.06; fixed softening .018a_ref (GR source supports only FixedSoftening). Not an anchor-matched FR24 forcing; m=0 measured separately; secular contribution conditional."))
        end
        println("GRTOTAL ",totals)
    else
        error("kind must be physical or weak")
    end
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    GRExcludedSector.run(ARGS[1],ARGS[2];Nq=length(ARGS)>2 ? parse(Int,ARGS[3]) : 256,
        n=length(ARGS)>3 ? parse(Int,ARGS[4]) : 8001,mmax=length(ARGS)>4 ? parse(Int,ARGS[5]) : 4,
        bank_prefix=length(ARGS)>5 ? ARGS[6] : nothing,
        weak_p=length(ARGS)>6 ? parse(Float64,ARGS[7]) : 1000.0,
        grid_strategy=length(ARGS)>7 ? Symbol(ARGS[8]) : :adaptive)
end

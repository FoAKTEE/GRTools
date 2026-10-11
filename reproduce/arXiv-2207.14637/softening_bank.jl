# Narrow softening experiment, never a full Newtonian production campaign.
# Sources: FR22 3_numerical_procedure.tex:72,79; FR24 4_eccentric_torques.tex:19.
using AccretionDisc, NumericsCore, LinearAlgebra, Serialization, TOML, Printf, SHA
const M = AccretionDisc.Modes
const F = AccretionDisc.ForcingModes
BLAS.set_num_threads(1)
cpu_seconds() = ccall(:clock,Clong,()) / 1e6

function setup(paper,rule,e,mmax)
    h,lo,K = paper == "fr22" ? (.1,.1,4) : (.06,.05,40)
    d = M.NewtonianDisc(h0=h,p=0.,q=0.,gamma=1.,beta=0.,soft=.3h,indirect=0,pole_disp=1e-6,lap_nphi=800)
    o = M.EccentricOrbit(1.,e,1.)
    soft = rule == "fixed" ? F.SemimajorScaleHeight(.3) :
           rule == "moving" ? F.InstantaneousScaleHeight(.3) :
           rule == "field" ? F.LocalScaleHeight(.3) : error("unknown softening rule")
    modes = [F.NewtonianMode(m,l,1.) for m in 1:mmax for l in max(1,m-K):m+K]
    (; d,o,soft,modes,lo,K,h,paper,rule,e,mmax)
end

function atomic_save(path,value)
    tmp=path*".tmp"
    serialize(tmp,value)
    mv(tmp,path;force=true)
end

function source_cache(s,rs,modes,dir;tag="source",block=64)
    mkpath(dir)
    pieces=[]
    for firsti in 1:block:length(rs)
        lasti=min(firsti+block-1,length(rs))
        path=joinpath(dir,"$tag-$firsti-$lasti.jls")
        if !isfile(path)
            c0=cpu_seconds();t0=time()
            b=M.potential_bank(s.d,[s.o],rs[firsti:lasti],modes;softening=s.soft,Nt=1024,Nphi=1024)
            atomic_save(path,(;bank=b,cpu_s=cpu_seconds()-c0,wall_s=time()-t0))
            println("source $tag $lasti/$(length(rs))");flush(stdout)
        end
        push!(pieces,deserialize(path))
    end
    caches=[M.PotentialCache{Float64}(vcat((p.bank.caches[1,j].lg for p in pieces)...),
         vcat((p.bank.caches[1,j].P for p in pieces)...),vcat((p.bank.caches[1,j].dP for p in pieces)...),
         vcat((p.bank.caches[1,j].d2P for p in pieces)...)) for j in eachindex(modes)]
    (;caches,cpu_s=sum(p.cpu_s for p in pieces))
end

function response(s,mode,cache,x)
    c0=cpu_seconds(); t0=time()
    f=M.solve_lm_mode(s.d,mode,x,cache;M.FR22_SOLVER_SETTINGS...)
    ur,up,sig=M.vars_from_enthalpy(s.d,f,x,mode,cache;lindblad_regular=true)
    # FR24Diagnostics.jl currently lacks the regular-state keyword. Reuse its
    # quadrature; use the SAME C2 potential and reconstruction as this solve.
    # Density: ref-paper/arXiv-2407.20398/src/Sections/2_linear_framework.tex:279-287.
    phi=M.RegularPotential(cache)
    # loggrid's final value is nextfloat(5.0) on this platform. A strict <=5
    # comparison drops a real integration node. The declared reporting grid
    # consists of the FIRST 4000 points, followed only by any BC extension.
    indices=1:4000
    isapprox(x[4000],5.;rtol=0,atol=8eps(5.)) || error("published 4000-node grid missing")
    xx=x[indices]; sigma=sig[indices]
    density=[-pi*r*mode.m*phi.spline(log10(r),1)*imag(sigma[i]) for (i,r) in enumerate(xx)]
    profile=M.fr24_torque_profile(xx,density)
    torque=(profile.torque[end]-profile.torque[1])*s.h^3
    (;mode,r=xx,sigma,torque,cpu_s=cpu_seconds()-c0,wall_s=time()-t0,
      rout=last(x),cache_nodes=length(cache.lg),grid_nodes=length(x))
end

function pilot(s,out)
    # Warm projection separately from timed samples, with the complete mask.
    M.potential_bank(s.d,[s.o],[1.],s.modes;softening=s.soft,Nt=1024,Nphi=1024)
    r=10.0 .^ range(log10(s.lo),log10(5.),length=17)
    c0=cpu_seconds(); t0=time()
    M.potential_bank(s.d,[s.o],r,s.modes;softening=s.soft,Nt=1024,Nphi=1024)
    source_cpu=cpu_seconds()-c0; source_wall=time()-t0
    chosen=[F.NewtonianMode(m,m,1.) for m in unique([1,20,80,s.mmax])]
    rs=10.0 .^ range(log10(s.lo),log10(5.),length=3997)
    bank=source_cache(s,rs,chosen,out;tag="pilot")
    x=M.loggrid(s.d,s.lo,5.,4000)
    response(s,chosen[1],bank.caches[1],x) # compile warmup
    rows=[]
    for (mode,cache) in zip(chosen,bank.caches)
        y=response(s,mode,cache,x)
        push!(rows,Dict("m"=>mode.m,"l"=>mode.m+mode.k,"cpu_s"=>y.cpu_s,"wall_s"=>y.wall_s,"torque"=>y.torque))
    end
    cost_per_m=maximum(row["cpu_s"]/row["m"] for row in rows if row["m"]>=20)
    projected=(source_cpu*3997/17+cost_per_m*sum(m.m for m in s.modes))/3600
    result=Dict("source_17_cpu_s"=>source_cpu,"source_17_wall_s"=>source_wall,
        "modes"=>length(s.modes),"per_mode"=>rows,"cost_per_m_s"=>cost_per_m,
        "projected_one_rule_core_h"=>projected,"projected_two_rules_core_h"=>2projected,
        "extension_and_refinement_reserve_fraction"=>.3)
    open(joinpath(out,"pilot.toml"),"w") do io;TOML.print(io,result);end
    println(result)
end

function run_bank(s,out)
    mkpath(out)
    config=(;s.paper,s.rule,s.e,s.mmax,s.K,source_nodes=3997,Nt=1024,Nphi=1024,
             grid_nodes=4000,settings=M.FR22_SOLVER_SETTINGS,mode_pairs=[(m.m,m.m+m.k) for m in s.modes],version=2)
    config_path=joinpath(out,"config.jls")
    if isfile(config_path)
        deserialize(config_path)==config || error("checkpoint configuration mismatch")
    else
        atomic_save(config_path,config)
    end
    rs=10.0 .^ range(log10(s.lo),log10(5.),length=3997)
    sources=source_cache(s,rs,s.modes,out)
    x=M.loggrid(s.d,s.lo,5.,4000)
    # Retain modes; extend near-edge and exterior OLR modes, integrate only the
    # published domain. The Lindblad landing measured 10a vs 12a for (6:9,1).
    edge=[j for (j,mode) in enumerate(s.modes) if begin
        rp=M.lindblad_locations(s.o,mode.m,mode.m+mode.k)[2]
        (6<=mode.m<=9 && mode.m+mode.k==1) || rp+2s.h*rp^1.5>=5
    end]
    extended=isempty(edge) ? nothing : source_cache(s,10.0 .^ range(log10(5.),log10(20.),length=2001),s.modes[edge],out;tag="extension-measured-v3")
    for (j,mode) in enumerate(s.modes)
        path=joinpath(out,"response-$(mode.m)-$(mode.m+mode.k).jls")
        isfile(path) && continue
        c=sources.caches[j]; xx=x
        if j in edge
            ce=extended.caches[findfirst(==(j),edge)]
            c=M.PotentialCache{Float64}(vcat(c.lg,ce.lg[2:end]),vcat(c.P,ce.P[2:end]),
                vcat(c.dP,ce.dP[2:end]),vcat(c.d2P,ce.d2P[2:end]))
            rplus=M.lindblad_locations(s.o,mode.m,mode.m+mode.k)[2]
            rout=rplus<5 ? 10. : 20.
            xx=vcat(x,10.0 .^ range(log10(5.),log10(rout),length=1001)[2:end])
        end
        try
            y=response(s,mode,c,xx)
            all(isfinite,y.sigma) && isfinite(y.torque) || error("nonfinite mode")
            atomic_save(path,y)
            @printf("response %d/%d (%d,%d) T=%.12g cpu=%.3f\n",j,length(s.modes),mode.m,mode.m+mode.k,y.torque,y.cpu_s)
        catch err
            open(joinpath(out,"failure-$(mode.m)-$(mode.m+mode.k).txt"),"w") do io
                showerror(io,err,catch_backtrace())
            end
            rethrow()
        end
        flush(stdout)
    end
    println("BANK FINISHED $(s.paper) $(s.rule) $(s.e)")
end

if abspath(PROGRAM_FILE)==@__FILE__
    stage,paper,rule=ARGS[1:3]
    e=parse(Float64,ARGS[4]);mmax=parse(Int,ARGS[5]);out=abspath(ARGS[6])
    s=setup(paper,rule,e,mmax);mkpath(out)
    if length(ARGS)>=8
        residue,divisor=parse.(Int,ARGS[7:8])
        s=merge(s,(;modes=filter(m -> mod(m.m,divisor)==residue,s.modes)))
    end
    stage=="pilot" ? pilot(s,out) : stage=="run" ? run_bank(s,out) : error("unknown stage")
end

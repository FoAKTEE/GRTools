# Usage: julia -t1 --project=capabilities/AccretionDisc lindblad_crossing_audit.jl
# regular-{native,ladder,refinement,oracle,edges} <prior-evidence-directory> <output-directory>
# The landed audit owns the 12-mode selection, disc settings and quadrature.
# Only deserialize checkpoints made by that audit (never source-mirror content).
include("linear_solver_audit.jl")
using Printf, Test
const PRIOR=abspath(ARGS[2])

function regular_observable(d,x,mode,c,f)
    ur,up,sig=vars_from_enthalpy(d,f,x,mode,c;lindblad_regular=true)
    smooth=MO.RegularPotential(c)
    td=[-pi*r*mode.m*smooth.spline(log10(r),1)*imag(sig[i]) for (i,r) in enumerate(x)]
    flux=pi .* x.^2 .* MO.Sigma.(Ref(d),x) .* real.(ur .* conj.(up))
    (;torque=trap(x,td),torque_abs=trap(x,abs.(td)),flux)
end

function regular_measured(d,x,mode,c;rtol=1e-12,integrator=NumericsCore.integrate_grid,kwargs...)
    t0=time_ns(); cpu0=cputime()
    s=integrator === NumericsCore.integrate_grid ?
        solve_lm_mode(d,mode,x,c;lindblad_regular=true,rtol,atol=rtol/100,dense=false,diagnostics=true,kwargs...) :
        MO._solve_m_mode(d,mode,x,MO.RegularPotential(c;integrator);rtol,atol=rtol/100,dense=false,diagnostics=true,kwargs...)
    cpu=cputime()-cpu0; wall=(time_ns()-t0)/1e9
    v=regular_observable(d,x,mode,c,s.mode)
    row=Dict{String,Any}("m"=>mode.m,"l"=>mode.m+mode.k,"rtol"=>rtol,
        "cpu_s"=>cpu,"wall_s"=>wall,"torque"=>v.torque,"torque_abs"=>v.torque_abs,
        "nx"=>length(x),"nr"=>length(c.lg),"boundary_r"=>x[collect(s.indices)],
        "X_before"=>collect(s.X_before),"X_after"=>collect(s.X_after))
    if haskey(s,:phase_safeguard)
        row["safeguard"]=string(s.phase_safeguard.status)
        row["fraction"]=s.phase_safeguard.fraction
    end
    @printf("regular m=%d l=%d rt=%.17g nr=%d nx=%d T=%.14g cpu=%.3f\n",mode.m,mode.m+mode.k,rtol,length(c.lg),length(x),v.torque,cpu)
    flush(stdout)
    row,s,v
end

function regular_native()
    d=disc((;h=.06,q=0.,p=0.)); x=loggrid(d,.05,5.,4000); mode=FM.NewtonianMode(40,40,1.)
    c=PotentialCache(d,40,first(x),last(x);n=4000)
    regular_measured(d,x,mode,c;rtol=1e-8,anchor=:legacy) # compile
    rows=[]
    for rt in [TOLS;nextfloat(1e-12)]
        row,s,v=regular_measured(d,x,mode,c;rtol=rt,anchor=:legacy)
        push!(rows,row); save("regular-native",Dict("rows"=>rows))
    end
end

function regular_ladder()
    rows=[]
    for (gi,g) in enumerate(GROUPS)
        d=disc(g); b=deserialize(joinpath(PRIOR,"bank-$gi.jls")); x=loggrid(d,.1,5.,4000)
        for (m,l) in g.ml
            mode=FM.NewtonianMode(m,l,1.); c=subcache(MO.harmonic_cache(b,1,mode),2)
            regular_measured(d,x,mode,c;rtol=1e-8,FR_OPTIONS...) # warm
            for rt in [TOLS;nextfloat(1e-12)]
                row,s,v=regular_measured(d,x,mode,c;rtol=rt,FR_OPTIONS...)
                merge!(row,Dict("group"=>gi,"e"=>g.e,"h"=>g.h,"q"=>g.q,"p"=>g.p))
                push!(rows,row);save("regular-ladder",Dict("rows"=>rows))
            end
        end
    end
end

function regular_refinement()
    rows=[]
    for gi in 0:3
        g=gi==0 ? (;h=.06,q=0.,p=0.,e=0.,ml=[(40,40)]) : GROUPS[gi]
        d=disc(g); lo=gi==0 ? .05 : .1
        b=gi==0 ? nothing : deserialize(joinpath(PRIOR,"bank-$gi.jls"))
        for (m,l) in g.ml
            mode=FM.NewtonianMode(m,l,1.)
            fine=gi==0 ? PotentialCache(d,m,lo,5.;n=7999) : MO.harmonic_cache(b,1,mode)
            opts=gi==0 ? (;anchor=:legacy) : FR_OPTIONS
            for (nx,stride) in ((4000,2),(7999,2),(15997,2),(7999,1))
                c=subcache(fine,stride); x=loggrid(d,lo,5.,nx)
                row,s,v=regular_measured(d,x,mode,c;rtol=1e-12,opts...)
                row["group"]=gi;push!(rows,row);save("regular-refinement",Dict("rows"=>rows))
            end
        end
    end
end

# DOPR853 is an independent 8th-order integrator, using the existing NumericsCore
# stage kernel. The schedule is fixed BEFORE evaluating the solution: uniform
# subdivisions of each prescribed radial cell, plus geometric cells resolving
# corotation. Halving max_phase and corotation_spacing is the Richardson check.
# No error estimator controls a step. The tolerance arguments are ignored here.
struct FixedDOP853
    refinement::Int
end
function (scheme::FixedDOP853)(rhs,y0,x,ico,dir;kwargs...)
    d=rhs.d; mode=rhs.mode
    a=first(x); b=last(x)
    for _ in 1:60
        mid=(a+b)/2
        mode.m*MO.Omega(d,mid)>mode.omega ? (a=mid) : (b=mid)
    end
    rc=(a+b)/2; width=d.pole_disp/abs(mode.m*MO.dOmegadr(d,rc))
    # Geometric near-corotation mesh, fixed in asinh((r-rc)/width).
    tlo=asinh((first(x)-rc)/width); thi=asinh((last(x)-rc)/width)
    dt=.08/scheme.refinement
    extra=rc .+ width .* sinh.(range(tlo,thi;length=ceil(Int,(thi-tlo)/dt)+1))
    mesh=sort!(unique!([collect(x); extra[2:end-1]]))
    start=searchsortedfirst(mesh,x[ico]); Y=zeros(length(x),length(y0));Y[ico,:]=y0
    y=copy(y0); yn=similar(y); arg=similar(y); k=[similar(y) for _ in 1:10]
    ode!(out,r,y,args)=rhs(out,r,y)
    oi=ico+dir
    for j in start:dir:(dir>0 ? length(mesh)-1 : 2)
        r=mesh[j]; endr=mesh[j+dir]
        wave=max(abs(MO.kwave(d,r,mode)),abs(MO.kwave(d,endr,mode)))
        n=max(1,ceil(Int,abs(endr-r)*wave*scheme.refinement/.12)); h=(endr-r)/n
        for s in 0:n-1
            NumericsCore.Dopr853._dormand_prince_steps!(ode!,r+s*h,y,h,nothing,k...,arg)
            NumericsCore.Dopr853._error!(yn,y,h,k[1],k[2],k[3],k[6],k[7],k[8],k[9],k[10],1.)
            y,yn=yn,y
        end
        if 1<=oi<=length(x) && endr==x[oi]
            Y[oi,:]=y; oi+=dir
        end
    end
    Y
end

function regular_oracle()
    rows=[]
    for gi in 0:3
        g=gi==0 ? (;h=.06,q=0.,p=0.,e=0.,ml=[(40,40)]) : GROUPS[gi]
        d=disc(g); lo=gi==0 ? .05 : .1; x=loggrid(d,lo,5.,4000)
        b=gi==0 ? nothing : deserialize(joinpath(PRIOR,"bank-$gi.jls"))
        for (m,l) in g.ml
            mode=FM.NewtonianMode(m,l,1.)
            c=gi==0 ? PotentialCache(d,m,first(x),last(x);n=4000) : subcache(MO.harmonic_cache(b,1,mode),2)
            opts=gi==0 ? (;anchor=:legacy) : FR_OPTIONS
            for n in (1,2,4)
                row,s,v=regular_measured(d,x,mode,c;rtol=1e-13,integrator=FixedDOP853(n),opts...)
                row["group"]=gi;row["refinement"]=n
                push!(rows,row);save("regular-oracle",Dict("rows"=>rows))
            end
        end
    end
end

function lindblad_roots(d,mode,lo,hi)
    grid=loggrid(d,lo,hi,2001); roots=Float64[]
    for i in 1:length(grid)-1
        a,b=grid[i:i+1]
        if signbit(MO.Dfreq(d,a,mode)) != signbit(MO.Dfreq(d,b,mode))
            for _ in 1:60
                mid=(a+b)/2
                signbit(MO.Dfreq(d,mid,mode))==signbit(MO.Dfreq(d,a,mode)) ? (a=mid) : (b=mid)
            end
            push!(roots,(a+b)/2)
        end
    end
    roots
end

function regular_edges()
    d=disc((;h=.1,q=0.,p=0.)); selected=[]
    # FR22 3_numerical_procedure.tex:72,79: positive l, m<=150, |l-m|<=8.
    for m in 1:150, l in max(1,m-8):m+8
        mode=FM.NewtonianMode(m,l,1.)
        roots=lindblad_roots(d,mode,.01,20.)
        near=[r for r in roots if min(abs(r-.1),abs(r-5.))<=2MO.csiso(d,r)/MO.Omega(d,r)]
        isempty(near) || push!(selected,Dict("m"=>m,"l"=>l,"roots"=>roots,"near_roots"=>near))
    end
    save("edge-selection",Dict("mmax"=>150,"dlmax"=>8,"rin"=>.1,"rout"=>5.,
        "h"=>.1,"q"=>0.,"p"=>0.,"modes"=>selected))
    println("EDGE MODES ",[(s["m"],s["l"]) for s in selected]);flush(stdout)
    modes=[FM.NewtonianMode(s["m"],s["l"],1.) for s in selected]
    cachefile=joinpath(OUT,"edge-bank.jls")
    if !isfile(cachefile)
        t0=time_ns();cpu0=cputime()
        b=MO.potential_bank(d,[MO.EccentricOrbit(1.,e,1.) for e in (.1,.25,.5)],
            loggrid(d,.1,12.,6001),modes;softening=FM.FixedSoftening(.03),Nt=1024,Nphi=1024)
        serialize(cachefile,b)
        save("edge-bank",Dict("cpu_s"=>cputime()-cpu0,"wall_s"=>(time_ns()-t0)/1e9,
            "nr"=>6001,"Nt"=>1024,"Nphi"=>1024,"eccentricities"=>[.1,.25,.5]))
    end
    b=deserialize(cachefile); rows=[]
    # Shared inner mesh isolates the change of boundary from torque quadrature.
    base=sort!(unique!([collect(loggrid(d,.1,12.,12001));5.;6.;8.;10.;12.]))
    for (ei,e) in enumerate((.1,.25,.5)), (mi,mode) in enumerate(modes)
        c=MO.harmonic_cache(b,ei,mode); rp=last(selected[mi]["roots"])
        for rout in (5.,6.,8.,10.,12.)
            x=base[base.<=rout]
            row,s,v=regular_measured(d,x,mode,c;rtol=1e-12,FR_OPTIONS...)
            inds=findall(r->rp+2MO.csiso(d,rp)/MO.Omega(d,rp)<=r<=rout-MO.csiso(d,rout)/MO.Omega(d,rout),x)
            row["rout"]=rout;row["e"]=e;row["r_olr"]=rp
            row["flux_outer"]=v.flux[end]
            row["plateau_samples"]=length(inds)
            if length(inds)>1
                fs=v.flux[inds]; row["plateau_interval"]=[x[first(inds)],x[last(inds)]]
                row["flux_plateau_spread"]=(maximum(fs)-minimum(fs))/maximum(abs,fs)
                row["flux_plateau_mean"]=sum(fs)/length(fs)
            end
            common=findall(r->r<=5.,x)
            row["common_domain_flux_l2"]=norm(v.flux[common])
            push!(rows,row);save("regular-edges",Dict("rows"=>rows))
            serialize(joinpath(OUT,"edge-$ei-$(mode.m)-$rout.jls"),(x=x,s=s,v=v))
        end
    end
end

if ARGS[1]=="regular-native"
    regular_native()
elseif ARGS[1]=="regular-ladder"
    regular_ladder()
elseif ARGS[1]=="regular-refinement"
    regular_refinement()
elseif ARGS[1]=="regular-oracle"
    regular_oracle()
elseif ARGS[1]=="regular-edges"
    regular_edges()
end

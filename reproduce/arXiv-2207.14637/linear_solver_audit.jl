# FR22 response audit. All integrations and reconstructions use the package solver.
# Torque convention: ref-paper/arXiv-2407.20398/src/Sections/2_linear_framework.tex:268-287
# (read as text; never executed).
using AccretionDisc, NumericsCore, LinearAlgebra, TOML, Serialization
const MO = AccretionDisc.Modes
const FM = AccretionDisc.ForcingModes
const OUT = abspath(last(ARGS))
mkpath(OUT)
const TOLS = [1e-8,1e-9,1e-10,1e-11,1e-12,1e-13]
const FR_OPTIONS = (; reint_steps=2, phgrad_steps=1, phgrad_safeguard=true, anchor=:corotation)
const GROUPS = [
    (;h=.1,q=0.,p=0.,e=.1,ml=[(1,1),(1,2),(2,1),(5,7)]),
    (;h=.06,q=0.,p=0.,e=.5,ml=[(1,3),(5,4),(40,40),(150,158)]),
    (;h=.08,q=1.,p=1.5,e=.25,ml=[(1,1),(2,3),(10,8),(80,84)])]
cputime() = ccall(:clock,Clong,()) / 1e6
trap(x,y) = sum(diff(x) .* (y[1:end-1] .+ y[2:end]) ./ 2)
disc(g) = NewtonianDisc(h0=g.h,q=g.q,p=g.p,gamma=1.,beta=0.,soft=.3g.h,indirect=0,pole_disp=1e-6,lap_nphi=800)
function observable(d,x,mode,cache,f)
    ur,up,sig = vars_from_enthalpy(d,f,x,mode,cache)
    phi = first.(MO.potential_at.(Ref(cache),x))
    td = -pi .* x .* mode.m .* phi .* imag.(sig)
    (;torque=trap(x,td), torque_abs=trap(x,abs.(td)), density=sig)
end
function measured(d,x,mode,cache; rtol, kwargs...)
    c0=cputime(); t0=time_ns()
    s=solve_lm_mode(d,mode,x,cache;rtol,atol=rtol/100,diagnostics=true,kwargs...)
    wall=(time_ns()-t0)/1e9; cpu=cputime()-c0
    v=observable(d,x,mode,cache,s.mode)
    row=Dict{String,Any}("m"=>mode.m,"l"=>mode.m+mode.k,"rtol"=>rtol,"atol"=>rtol/100,
        "cpu_s"=>cpu,"wall_s"=>wall,"torque"=>v.torque,"torque_abs"=>v.torque_abs,
        "indices"=>collect(s.indices),"boundary_r"=>x[collect(s.indices)],
        "X_before"=>collect(s.X_before),"X_after"=>collect(s.X_after),
        "peak_h"=>maximum(abs,s.mode[:,1]),"relative_state_change"=>collect(s.relative_state_change))
    if haskey(s,:phase_safeguard)
        a=s.phase_safeguard
        row["safeguard"]=Dict("status"=>string(a.status),"fraction"=>a.fraction,
            "inner_response_change"=>a.inner_response_change)
    end
    println("m=$(mode.m) l=$(mode.m+mode.k) rt=$rtol T=$(v.torque) cpu=$cpu X=$(s.X_after)");flush(stdout)
    row,s,v
end
function save(name, data)
    open(joinpath(OUT,name*".toml"),"w") do io; TOML.print(io,data);end
end
function bank(g,gi)
    file=joinpath(OUT,"bank-$gi.jls")
    isfile(file) && return deserialize(file)
    d=disc(g); modes=[FM.NewtonianMode(m,l,1.) for (m,l) in g.ml]
    rs=loggrid(d,.1,5.,7993)
    c0=cputime();t0=time_ns()
    b=MO.potential_bank(d,[MO.EccentricOrbit(1.,g.e,1.)],rs,modes;
        softening=FM.FixedSoftening(.3g.h),Nt=1024,Nphi=1024,method=:direct)
    serialize(file,b)
    save("bank-$gi",Dict("cpu_s"=>cputime()-c0,"wall_s"=>(time_ns()-t0)/1e9,
        "nr"=>7993,"Nt"=>1024,"Nphi"=>1024,"e"=>g.e,"h"=>g.h,"q"=>g.q,"p"=>g.p))
    b
end
function subcache(c,stride)
    inds=1:stride:length(c.lg)
    MO.PotentialCache(c.lg[inds],c.P[inds],c.dP[inds],c.d2P[inds])
end
function diagnosis()
    d=disc((;h=.06,q=0.,p=0.)); x=loggrid(d,.05,5.,4000)
    mode=FM.NewtonianMode(40,40,1.); cache=MO.PotentialCache(d,40,.05,5.;n=1000)
    solve_lm_mode(d,mode,x,cache;rtol=1e-8,atol=1e-10) # warm compilation
    rows=Dict{String,Any}[]
    for dense in (true,false), rt in 10.0 .^ (-8:-1:-13)
        row,s,v=measured(d,x,mode,cache;rtol=rt,dense,phgrad_steps=0,anchor=:legacy)
        row["dense"]=dense; push!(rows,row); save("diagnosis",Dict("rows"=>rows))
        serialize(joinpath(OUT,"diagnosis-$dense-$rt.jls"),(x=x,s=s,v=v))
    end
end
function ladder()
    rows=Dict{String,Any}[]
    for (gi,g) in enumerate(GROUPS)
        d=disc(g); b=bank(g,gi); x=loggrid(d,.1,5.,4000)
        for (m,l) in g.ml
            mode=FM.NewtonianMode(m,l,1.); c=subcache(MO.harmonic_cache(b,1,mode),2)
            solve_lm_mode(d,mode,x,c;rtol=1e-8,atol=1e-10,dense=false,FR_OPTIONS...) # warm
            for rt in TOLS
                row,s,v=measured(d,x,mode,c;rtol=rt,dense=false,FR_OPTIONS...)
                merge!(row,Dict("group"=>gi,"e"=>g.e,"h"=>g.h,"q"=>g.q,"p"=>g.p,"nr"=>3997,"nx"=>4000))
                push!(rows,row); save("ladder",Dict("rows"=>rows))
                serialize(joinpath(OUT,"mode-$gi-$m-$l-$rt.jls"),(x=x,s=s,v=v))
            end
        end
    end
end
ARGS[1]=="diagnosis" && diagnosis()
ARGS[1]=="ladder" && ladder()

function refinements()
    rows=Dict{String,Any}[]
    for (gi,g) in enumerate(GROUPS)
        d=disc(g); b=bank(g,gi)
        for (m,l) in g.ml
            mode=FM.NewtonianMode(m,l,1.); fine=MO.harmonic_cache(b,1,mode)
            ref=deserialize(joinpath(OUT,"mode-$gi-$m-$l-1.0e-13.jls"))
            for (nx,stride) in ((7999,2),(4000,1),(4000,4),(4000,8))
                x=loggrid(d,.1,5.,nx); c=subcache(fine,stride)
                row,s,v=measured(d,x,mode,c;rtol=1e-12,dense=false,FR_OPTIONS...)
                merge!(row,Dict("group"=>gi,"nx"=>nx,"nr"=>length(c.lg),
                    "torque_delta_vs_ladder"=>abs(v.torque-ref.v.torque),
                    "torque_relative_vs_ladder"=>abs(v.torque-ref.v.torque)/abs(ref.v.torque)))
                if nx==4000
                    row["response_relative_vs_ladder"]=norm(s.mode[:,1]-ref.s.mode[:,1])/norm(ref.s.mode[:,1])
                else
                    row["response_relative_vs_ladder"]=norm(s.mode[1:2:end,1]-ref.s.mode[:,1])/norm(ref.s.mode[:,1])
                end
                push!(rows,row);save("refinement",Dict("rows"=>rows))
                serialize(joinpath(OUT,"refine-$gi-$m-$l-$nx-$stride.jls"),(x=x,s=s,v=v))
            end
        end
    end
end
function anchors()
    rows=Dict{String,Any}[]; maps=Dict{String,Any}[]
    for (gi,g) in enumerate(GROUPS)
        d=disc(g);b=bank(g,gi); x=loggrid(d,.1,5.,4000); az=range(0,2pi;length=129)[1:end-1]
        a=zeros(length(x),length(az),4);bmap=zero(a)
        for (m,l) in g.ml
            mode=FM.NewtonianMode(m,l,1.); c=subcache(MO.harmonic_cache(b,1,mode),2)
            ref=deserialize(joinpath(OUT,"mode-$gi-$m-$l-1.0e-12.jls"))
            row,s,v=measured(d,x,mode,c;rtol=1e-12,dense=false,FR_OPTIONS...,anchor=:legacy)
            row["group"]=gi;row["reference_indices"]=collect(ref.s.indices)
            row["corotation_grid_radius"]=x[MO.corotation_index(d,x,mode)]
            row["pressureless_corotation"]=(m/l)^(2/3)
            row["relative_frequency_residual"]=abs(mode.m*MO.Omega(d,row["corotation_grid_radius"])-mode.omega)/mode.omega
            row["response_relative"]=norm(s.mode[:,1]-ref.s.mode[:,1])/norm(ref.s.mode[:,1])
            # Direct finite Fourier sum on the same sampled map (FR22 eq 15).
            # Diagnostic subset only: no claim of mode-truncation convergence.
            for (ti,t) in enumerate((0.,pi/2,pi,3pi/2)), j in eachindex(az)
                phase=cis(m*az[j]-l*t)
                a[:,j,ti] .+= real.(ref.v.density .* phase)
                bmap[:,j,ti] .+= real.(v.density .* phase)
            end
            push!(rows,row);save("anchors",Dict("rows"=>rows,"maps"=>maps))
        end
        push!(maps,Dict("group"=>gi,"modes"=>[collect(v) for v in g.ml],
            "relative_L2"=>[norm(bmap[:,:,t]-a[:,:,t])/norm(a[:,:,t]) for t in 1:4],
            "relative_Linf"=>[maximum(abs,bmap[:,:,t]-a[:,:,t])/maximum(abs,a[:,:,t]) for t in 1:4]))
        save("anchors",Dict("rows"=>rows,"maps"=>maps))
        serialize(joinpath(OUT,"maps-$gi.jls"),(corotation=a,legacy=bmap))
    end
end
ARGS[1]=="refinement" && refinements()
ARGS[1]=="anchors" && anchors()

function diagnosis_exact()
    d=disc((;h=.06,q=0.,p=0.));x=loggrid(d,.05,5.,4000)
    mode=FM.NewtonianMode(40,40,1.);c=MO.PotentialCache(d,40,.05,5.;n=1000)
    rows=Dict{String,Any}[]
    for dense in (true,false)
        row,s,v=measured(d,x,mode,c;rtol=1e-12,dense,phgrad_steps=0,anchor=:legacy)
        row["dense"]=dense
        near=deserialize(joinpath(OUT,"diagnosis-$dense-1.0000000000000002e-12.jls"))
        ref=deserialize(joinpath(OUT,"diagnosis-$dense-1.0e-13.jls"))
        delta=abs.(s.mode[:,1]-near.s.mode[:,1]);ii=argmax(delta)
        row["ulp_torque_change"]=v.torque-near.v.torque
        row["ulp_response_max_location"]=x[ii]
        row["ulp_response_max_change"]=delta[ii]
        row["relative_to_1e13"]=abs(v.torque-ref.v.torque)/abs(ref.v.torque)
        push!(rows,row);save("diagnosis-exact",Dict("rows"=>rows))
        serialize(joinpath(OUT,"diagnosis-exact-$dense.jls"),(x=x,s=s,v=v))
    end
end
ARGS[1]=="diagnosis-exact" && diagnosis_exact()

function diagnosis_locations()
    rows=Dict{String,Any}[]
    d=disc((;h=.06,q=0.,p=0.));mode=FM.NewtonianMode(40,40,1.)
    for dense in (true,false)
        ref=deserialize(joinpath(OUT,"diagnosis-$dense-1.0e-13.jls"));x=ref.x
        c=MO.PotentialCache(d,40,.05,5.;n=1000)
        phi=first.(MO.potential_at.(Ref(c),x))
        for rt in 10.0 .^ (-8:-1:-12)
            s=deserialize(joinpath(OUT,"diagnosis-$dense-$rt.jls"))
            delta=s.s.mode[:,1]-ref.s.mode[:,1]
            dt=-pi .* x .* mode.m .* phi .* imag.(s.v.density-ref.v.density)
            ii=argmax(abs.(dt)); near=abs.(x .- 1) .< .05
            push!(rows,Dict("dense"=>dense,"rtol"=>rt,"max_dh_radius"=>x[argmax(abs.(delta))],
                "max_dtorque_radius"=>x[ii],"max_dtorque_density"=>dt[ii],
                "torque_delta"=>trap(x,dt),"torque_delta_near_planet"=>trap(x,dt.*near),
                "relative_response"=>norm(delta)/norm(ref.s.mode[:,1]),
                "peak_reference"=>maximum(abs,ref.s.mode[:,1]),
                "peak_response"=>maximum(abs,s.s.mode[:,1])))
        end
    end
    save("diagnosis-locations",Dict("rows"=>rows))
    display(rows)
end
ARGS[1]=="diagnosis-locations" && diagnosis_locations()

function edge_check()
    g=(;h=.1,q=0.,p=0.,e=.5);d=disc(g);mode=FM.NewtonianMode(9,1,1.)
    rg=loggrid(d,.1,8.,3997)
    b=MO.potential_bank(d,[MO.EccentricOrbit(1.,g.e,1.)],rg,[mode];
        softening=FM.FixedSoftening(.3g.h),Nt=1024,Nphi=1024)
    c=MO.harmonic_cache(b,1,mode);x=loggrid(d,.1,5.,4000)
    extended=vcat(x,exp.(range(log(5.),log(8.);length=1001))[2:end])
    rows=Dict{String,Any}[];sols=[]
    for xx in (x,extended)
        row,s,v=measured(d,xx,mode,c;rtol=1e-12,dense=false,FR_OPTIONS...)
        row["rout"]=last(xx);push!(rows,row);push!(sols,s.mode)
    end
    wave=findall(r -> r >= 10.0^(2/3),x)
    report=Dict("rows"=>rows,"relative_response_on_original_domain"=>norm(sols[1][:,1]-sols[2][1:4000,1])/norm(sols[2][1:4000,1]),
        "relative_response_in_outer_wave_zone"=>norm(sols[1][wave,1]-sols[2][wave,1])/norm(sols[2][wave,1]),
        "r_corotation"=>(9/1)^(2/3),"r_outer_LR"=>10.0^(2/3),
        "H_outer_LR"=>MO.csiso(d,10.0^(2/3))/MO.Omega(d,10.0^(2/3)))
    save("near-edge",report)
end
ARGS[1]=="near-edge" && edge_check()

function diagnosis_cache()
    d=disc((;h=.06,q=0.,p=0.));x=loggrid(d,.05,5.,4000)
    mode=FM.NewtonianMode(40,40,1.);c=MO.PotentialCache(d,40,.05,5.;n=4000)
    rows=Dict{String,Any}[]
    for rt in (1e-10,1e-12,nextfloat(1e-12),1e-13)
        row,s,v=measured(d,x,mode,c;rtol=rt,dense=false,phgrad_steps=0,anchor=:legacy)
        row["nr"]=4000;push!(rows,row);save("diagnosis-cache",Dict("rows"=>rows))
    end
end
ARGS[1]=="diagnosis-cache" && diagnosis_cache()

function diagnosis_conditioning()
    d=disc((;h=.06,q=0.,p=0.));x=loggrid(d,.05,5.,4000)
    mode=FM.NewtonianMode(40,40,1.);c=MO.PotentialCache(d,40,.05,5.;n=1000)
    rows=Dict{String,Any}[]
    for rt in (1e-12,nextfloat(1e-12))
        h1=MO.integrate_master(d,x,(1.0+0im,0.0+0im),mode,0,c;rtol=rt,atol=rt/100,dense=false)
        h2=MO.integrate_master(d,x,(0.0+0im,1.0+0im),mode,0,c;rtol=rt,atol=rt/100,dense=false)
        Ci,Co=MO.wkb_boundary_coeffs(d,x,1,length(x),mode)
        A=[h1[1,2]-Ci*h1[1,1] h2[1,2]-Ci*h2[1,1];
           h1[end,2]-Co*h1[end,1] h2[end,2]-Co*h2[end,1]]
        scaled=A ./ sqrt.(sum(abs2,A;dims=2))
        push!(rows,Dict("rtol"=>rt,"boundary_matrix_condition"=>cond(A),
            "row_normalized_condition"=>cond(scaled),"matrix_real"=>[collect(real.(A[i,:])) for i in 1:2],
            "matrix_imag"=>[collect(imag.(A[i,:])) for i in 1:2]))
        save("diagnosis-conditioning",Dict("rows"=>rows))
    end
end
ARGS[1]=="diagnosis-conditioning" && diagnosis_conditioning()

function active_anchor()
    g=GROUPS[1];d=disc(g);mode=FM.NewtonianMode(1,3,1.)
    x=loggrid(d,.1,5.,4000)
    b=MO.potential_bank(d,[MO.EccentricOrbit(1.,g.e,1.)],loggrid(d,.1,5.,3997),[mode];
        softening=FM.FixedSoftening(.3g.h),Nt=1024,Nphi=1024)
    c=MO.harmonic_cache(b,1,mode)
    rows=Dict{String,Any}[];sols=[];dens=[]
    for anchor in (:corotation,:legacy)
        row,s,v=measured(d,x,mode,c;rtol=1e-10,dense=false,FR_OPTIONS...,anchor)
        row["anchor"]=string(anchor);push!(rows,row);push!(sols,s);push!(dens,v.density)
    end
    maps=deserialize(joinpath(OUT,"maps-1.jls"));a=copy(maps.corotation);bmap=copy(maps.legacy)
    az=range(0,2pi;length=129)[1:end-1]
    for (ti,t) in enumerate((0.,pi/2,pi,3pi/2)), j in eachindex(az)
        phase=cis(az[j]-3t)
        a[:,j,ti] .+= real.(dens[1].*phase)
        bmap[:,j,ti] .+= real.(dens[2].*phase)
    end
    report=Dict("rows"=>rows,"response_relative"=>norm(sols[2].mode[:,1]-sols[1].mode[:,1])/norm(sols[1].mode[:,1]),
        "map_modes"=>[[1,1],[1,2],[1,3],[2,1],[5,7]],"map_other_modes_rtol"=>1e-12,"added_mode_rtol"=>1e-10,
        "relative_L2"=>[norm(bmap[:,:,t]-a[:,:,t])/norm(a[:,:,t]) for t in 1:4],
        "relative_Linf"=>[maximum(abs,bmap[:,:,t]-a[:,:,t])/maximum(abs,a[:,:,t]) for t in 1:4])
    save("active-anchor",report)
    serialize(joinpath(OUT,"active-anchor.jls"),(x=x,sols=sols,corotation=a,legacy=bmap))
end
ARGS[1]=="active-anchor" && active_anchor()

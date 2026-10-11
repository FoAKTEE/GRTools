"""Float64 convergence measurements using the existing eccentric source and BVP owners.

Run with the AccretionDisc project: convergence_study.jl source|radial EVIDENCE_DIR.
PREREG.json must already exist. Per-case files support restart without replacing evidence.
Numerical gate failures are data; the independent evidence verifier exposes them.
"""
module EccentricConvergenceStudy
using AccretionDisc, Spacetime, LinearAlgebra, Printf
using NumericsCore: local_interp
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
include(joinpath(@__DIR__, "../_lib/repro_harness.jl"))
using .ReproHarness: json_write, json_read, encode_value, decode_value

const SOURCE_MODES = ((0,1),(0,40),(2,0),(2,1),(8,-40),(8,40),(20,0),(1,-1))
cpu_seconds() = Float64(ccall(:clock, Clong, ()))/1e6
writejson(path, x) = write(path, json_write(encode_value(x)))

function fixture(a, p, e; pole=1e-3)
    rref = p/(1-e^2)
    d = SW.DiscModel(;χ=a,r_ref=rref,h_ref=0.06,pSigma=1.5,qTemp=0.0,
                     closure=Adiabatic,evolution=:fixed,εp=pole)
    pert = SW.EccentricPerturber(EquatorialOrbit(a,p,e);mass=1e-5,
                               softening=FixedSoftening(0.018rref))
    d, pert
end

"Coordinate dimensions and fixed physical scale, independent of measured coefficients."
function source_floors(d, pert)
    r = d.r_ref
    weights = (one(r),r,r^2,one(r),one(r),one(r),r)
    [1e-12*2pert.mass/r*weights[(j+1)÷2]/(iseven(j) ? r : one(r)) for j in 1:14]
end

"Each mode and component uses its own radial infinity norm and preregistered floor."
function source_error(a, b, floors)
    [maximum(abs(a[i,j][n]-b[i,j][n]) for i in axes(a,1)) /
     max(maximum(abs(b[i,j][n]) for i in axes(b,1)),floors[n])
     for j in axes(a,2), n in 1:14]
end

"""PREREG_v2 acceptance: relative 1e-8 plus the preregistered absolute floor, as in
DESIGN section 4 and test_eccentric_source_bank.jl `source_budget_ratio`. A ratio <= 1 passes.
`source_error` multiplies the floor by 1e-8 (1e-20 of the source scale); it stays a diagnostic."""
function source_budget(a, b, floors)
    [maximum(abs(a[i,j][n]-b[i,j][n]) for i in axes(a,1)) /
     (1e-8*maximum(abs(b[i,j][n]) for i in axes(b,1))+floors[n])
     for j in axes(a,2), n in 1:14]
end

"The source-bank test unit (mu/p, p-weighted), a stricter sensitivity floor for PREREG_v2."
function test_unit_floors(pert)
    p = pert.orbit.p
    weights = (one(p),p,p^2,one(p),one(p),one(p),p)
    [1e-12*pert.mass/p*weights[(j+1)÷2]/(iseven(j) ? p : one(p)) for j in 1:14]
end

"Five-point fourth-order centred difference (scalars or arrays)."
fourth_difference(fm2, fm1, fp1, fp2, h) = (fm2-8fm1+8fp1-fp2)/(12h)

const V2_BUDGETS = ("orbit_budget","orbit_budget_sensitivity","azimuth_budget",
                    "azimuth_budget_sensitivity","cache_budget","cache_budget_sensitivity")
_flat(x) = x isa AbstractVector ? reduce(vcat,_flat.(x);init=Float64[]) : [x]

"PREREG_v2 gate from the reported arrays: every budget <= 1, finest fourth-order audit <= 1e-3."
prereg_v2_pass(r) = maximum(k->maximum(_flat(r[k])),V2_BUDGETS) <= 1 &&
    maximum(_flat(last(r["derivative_error_v2"]))) <= 1e-3

"Original PREREG.json gate, retained as the strict diagnostic."
prereg_v1_pass(r) = maximum(k->maximum(_flat(r[k])),("orbit_error","azimuth_error","cache_error")) <= 1e-8 &&
    maximum(_flat(last(r["derivative_error"]))) <= 1e-3

"A zero or roundoff plateau has no measured convergence order."
observed_order(d1, d2; plateau=1e-13) =
    min(d1,d2)>plateau ? log2(d1/d2) : nothing

function timed(f, costs, label, parameters)
    wall = time(); cpu = cpu_seconds()
    result = f()
    push!(costs, merge(Dict("axis_level"=>label,"cpu_seconds"=>cpu_seconds()-cpu,
                           "wall_seconds"=>time()-wall),parameters))
    println("COST ",label," ",last(costs)); flush(stdout)
    result
end

function source_table(d, pert, radii, modes, nq, nphi, costs, label; blocksize=1)
    bank = timed(costs,label,Dict("Nq"=>nq,"Nphi_floor"=>nphi,"radii"=>length(radii))) do
        SW.source_azimuth_bank(d,pert,radii;Nq=nq,Nphi=nphi,
                              mmax=maximum(m.m for m in modes),blocksize)
    end
    table = timed(costs,label*"_sums",Dict("modes"=>length(modes))) do
        SW.source_modes(bank,modes)
    end
    table, Dict("actual_Nphi"=>bank.nphi,"strip_min"=>minimum(bank.strips),
                "workspace_bound_bytes"=>bank.workspace_bound_bytes)
end

"PREREG.json + PREREG_AMENDMENT_01 + PREREG_v2; small `nq`/`mode_list` serve the unit test only."
function source_case(a,p,e; nq=1024, mode_list=SOURCE_MODES)
    d,pert = fixture(a,p,e)
    modes = [HarmonicMode(m,k,orbit_frequency(pert.orbit,m,k)) for (m,k) in mode_list]
    radii = [p/(1+e),p,p/(1-e)]
    floors = source_floors(d,pert)
    sensitivity = test_unit_floors(pert)
    costs = Dict[]
    # This prepass invokes the source owner's strip bound, not a new estimator.
    T = typeof(d.χ)
    delta = minimum(SW.strip_halfwidth(SW._eccentric_source_event(pert,
        worldline_rates(pert.orbit,2T(pi)*i/4096)),r) for i in 0:4095 for r in radii)
    nphi = SW._eccentric_nphi(delta,maximum(first,mode_list),T)
    tables = Any[]; metadata = Dict[]
    for level in (nq,2nq,4nq)
        h,meta = source_table(d,pert,radii,modes,level,nphi,costs,"orbit_$level")
        push!(tables,h); push!(metadata,meta)
    end
    errors_q = [source_error(tables[i],tables[i+1],floors) for i in 1:2]
    budget_q = [source_budget(tables[i],tables[i+1],f) for f in (floors,sensitivity) for i in 1:2]
    # PREREG_AMENDMENT_01: isolate azimuth accuracy at fixed Nq=1024.
    phitables = Any[first(tables)]
    for level in (2,4)
        h,meta = source_table(d,pert,radii,modes,nq,level*nphi,costs,"azimuth_$level")
        push!(phitables,h); push!(metadata,meta)
    end
    errors_phi = [source_error(phitables[i],phitables[i+1],floors) for i in 1:2]
    budget_phi = [source_budget(phitables[i],phitables[i+1],f) for f in (floors,sensitivity) for i in 1:2]
    # Direct radial tables have no interpolation error. Refinement changes their
    # set of samples; repeated common radii must remain invariant independently.
    errors_cache = typeof(first(errors_q))[]
    budget_cache = typeof(first(errors_q))[]
    raw_cache = Any[]
    cachebase,_ = source_table(d,pert,radii,modes,32,nphi,costs,"cache_3")
    for level in (2,4)
        rr = sort!(unique!(vcat(collect(range(radii[1],radii[2];length=level+1)),
                                 collect(range(radii[2],radii[3];length=level+1)))))
        h,_ = source_table(d,pert,rr,modes,32,nphi,costs,"cache_$(length(rr))";
                          blocksize=level)
        indices = [findfirst(==(r),rr) for r in radii]
        push!(raw_cache,encode_value(h[indices,:]))
        push!(errors_cache,source_error(h[indices,:],cachebase,floors))
        append!(budget_cache,[source_budget(h[indices,:],cachebase,f) for f in (floors,sensitivity)])
    end
    # PREREG_v2: the original centred ladder (0.08,0.04,0.02) stays a diagnostic; the
    # five-point stencil at h/eps=0.02,0.01,0.005 reuses the samples at h and 2h.
    jet = [tables[1][2,j][2n] for j in eachindex(modes), n in 1:7]
    relative(D) = [abs(D[j,n]-jet[j,n])/max(abs(jet[j,n]),floors[2n]) for j in eachindex(modes), n in 1:7]
    value(h,i) = [h[i,j][2n-1] for j in eachindex(modes), n in 1:7]
    samples = Dict{Float64,Any}(); steps = Dict{Float64,Float64}()
    raw_derivatives = Any[]
    for fraction in (0.08,0.04,0.02,0.01,0.005)
        step = fraction*0.018d.r_ref
        h,_ = source_table(d,pert,[p-step,p+step],modes,nq,nphi,costs,"derivative_$fraction")
        samples[fraction] = h; steps[fraction] = step
        push!(raw_derivatives,Dict("step"=>step,"fraction"=>fraction,"table"=>encode_value(h)))
    end
    derivative_errors = [relative((value(samples[f],2)-value(samples[f],1))/(2steps[f]))
                         for f in (0.08,0.04,0.02)]
    derivative_errors_v2 = [relative(fourth_difference(value(samples[f2],1),value(samples[f],1),
                                                       value(samples[f],2),value(samples[f2],2),steps[f]))
                            for (f,f2) in ((0.02,0.04),(0.01,0.02),(0.005,0.01))]
    # Raw projections are retained for independent recomputation of every norm.
    result = Dict("a"=>a,"p"=>p,"e"=>e,"radii"=>radii,"modes"=>collect.(mode_list),"Nq_levels"=>[nq,2nq,4nq],
         "source_floors"=>floors,"sensitivity_floors"=>sensitivity,"Nphi_base"=>nphi,"strip_min"=>delta,
         "protocol"=>"PREREG_AMENDMENT_01+PREREG_v2",
         "orbit_error"=>[collect.(eachrow(x)) for x in errors_q],
         "azimuth_error"=>[collect.(eachrow(x)) for x in errors_phi],
         "cache_error"=>[collect.(eachrow(x)) for x in errors_cache],
         "derivative_error"=>[collect.(eachrow(x)) for x in derivative_errors],
         "derivative_error_v2"=>[collect.(eachrow(x)) for x in derivative_errors_v2],
         "orbit_budget"=>[collect.(eachrow(x)) for x in budget_q[1:2]],
         "orbit_budget_sensitivity"=>[collect.(eachrow(x)) for x in budget_q[3:4]],
         "azimuth_budget"=>[collect.(eachrow(x)) for x in budget_phi[1:2]],
         "azimuth_budget_sensitivity"=>[collect.(eachrow(x)) for x in budget_phi[3:4]],
         "cache_budget"=>[collect.(eachrow(x)) for x in budget_cache[1:2:end]],
         "cache_budget_sensitivity"=>[collect.(eachrow(x)) for x in budget_cache[2:2:end]],
         "orbit_order"=>[observed_order(errors_q[1][j,n],errors_q[2][j,n]) for j in 1:length(modes) for n in 1:14],
         "azimuth_order"=>[observed_order(errors_phi[1][j,n],errors_phi[2][j,n]) for j in 1:length(modes) for n in 1:14],
         "raw_orbit"=>[encode_value([collect(values(h[i,j])) for i in axes(h,1),j in axes(h,2)]) for h in tables],
         "raw_azimuth"=>[encode_value([collect(values(h[i,j])) for i in axes(h,1),j in axes(h,2)]) for h in phitables],
         "raw_cache"=>raw_cache,"raw_cache_base"=>encode_value(cachebase),
         "raw_derivatives"=>raw_derivatives,
         "metadata"=>metadata,"costs"=>costs)
    result["pass_v1_strict"] = prereg_v1_pass(result)
    result["pass"] = prereg_v2_pass(result)
    result
end

"Compare both response components independently on fixed common radii."
function response_error(a,b,radii)
    floors = (1e-12*1e-5/b.context.disc.r_ref,
              1e-12*1e-5/sqrt(b.context.disc.r_ref))
    weights = [(radii[min(i+1,end)]-radii[max(i-1,1)])/2 for i in eachindex(radii)]
    map(1:2) do j
        va = local_interp(radii,a.r,vec(a.y[j,:]);npts=4)
        vb = local_interp(radii,b.r,vec(b.y[j,:]);npts=4)
        sqrt(sum(weights.*abs2.(va-vb))) /
            max(sqrt(sum(weights.*abs2.(vb))),sqrt(sum(weights))*floors[j])
    end
end

"Local WKB contamination; forcing and basis gradients prevent a pure reflection interpretation."
function reflection_sample(sol,r,outward)
    c = sol.context
    h = cbrt(eps(typeof(r)))*max(abs(c.mode.omega),abs(c.mode.m*SW.Omega(c.disc,r)))
    # Independent centered frequency derivative for the existing branch selector.
    cp = SW.ModeContext(c.disc,HarmonicMode(c.mode.m,c.mode.k,c.mode.omega+h);
                        closure=c.closure,domain=c.domain)
    cm = SW.ModeContext(c.disc,HarmonicMode(c.mode.m,c.mode.k,c.mode.omega-h);
                        closure=c.closure,domain=c.domain)
    C = SW.master_operator(c,r;landau=false)
    Cw = (SW.master_operator(cp,r;landau=false)-SW.master_operator(cm,r;landau=false))/(2h)
    y = [local_interp(r,sol.r,vec(sol.y[j,:]);npts=4) for j in 1:2]
    w = SW.wkb_amplitudes(C,Cw,y;outward)
    incoming = abs(w.A_toward)*norm(w.Y_toward)
    outgoing = abs(w.A_away)*norm(w.Y_away)
    floor = hypot(1e-17/c.disc.r_ref,1e-17/sqrt(c.disc.r_ref))
    Dict("r"=>r,"incoming_fraction"=>incoming/max(incoming+outgoing,floor),
         "incoming_norm"=>incoming,"outgoing_norm"=>outgoing)
end

function radial_context(a,p,e,m,k,pole,edge)
    d,pert = fixture(a,p,e;pole)
    ilo,peri,apo = ((1.05,0.6,2.0),(1.01,0.4,3.0),(1.001,0.25,4.0))[edge]
    domain = (max(ilo*SW.r_isco(d),peri*p/(1+e)),apo*p/(1-e))
    SW.ModeContext(d,pert,m,k;closure=Adiabatic,domain),pert
end

function radial_solve(c,pert,n,nq,costs,label;nanc=512,cache=false)
    width = SW.corotation_width(c)
    fraction = width === nothing ? 0.01 : 16width/c.r_c
    nw = max(1,n÷64) # grid.jl's nL/nR are cells per WKB radian, not node counts.
    sol = timed(costs,label,Dict("n"=>n,"Nq"=>nq,"n_anc"=>nanc,
                                "phase_cells_per_radian"=>nw,
                                "pole"=>c.disc.εp,
                                "domain"=>collect(c.domain))) do
        grid = SW.build_grid(c;nL=nw,nR=nw,nC=n+1,n_anc=nanc,σ_wc_frac=fraction)
        ops = cache ? SW.mode_operator_cache(c,grid) : nothing
        SW.solve_mode(c,pert;grid,cache=ops,bank_options=(;Nq=nq,blocksize=1),
                      radial_batch_size=16)
    end
    last(costs)["nodes"] = length(sol.r)
    last(costs)["bvp_wall_seconds"] = sol.seconds
    last(costs)["central_spacing_over_pole_width"] = width===nothing ? nothing :
        maximum(sol.r[i+1]-sol.r[i] for i in 1:length(sol.r)-1 if
                abs((sol.r[i+1]+sol.r[i])/2-c.r_c)<=width)/width
    # Bounded runs retain finished refinement levels even if a later level times out.
    checkpoint = get(ENV,"ECC_CONVERGENCE_CHECKPOINT_DIR","")
    if !isempty(checkpoint)
        writejson(joinpath(checkpoint,label*".json"),Dict("r"=>sol.r,"y"=>encode_value(sol.y),
                  "a"=>pert.orbit.a,"p"=>pert.orbit.p,"e"=>pert.orbit.e,
                  "m"=>c.mode.m,"k"=>c.mode.k,"r_c"=>c.r_c,
                  "residual"=>sol.residual,"bc_residual"=>sol.bc_residual,
                  "cost"=>last(costs),"pole_width"=>width))
    end
    sol
end

"Vary only the causal prescription on a fixed grid, reusing its exact source table."
function causal_audit(c,pert,base,costs,label)
    c.r_c===nothing && return [base,base,base]
    solutions = [base]
    for pole in (5e-4,2.5e-4)
        cc,_ = radial_context(pert.orbit.a,pert.orbit.p,pert.orbit.e,c.mode.m,c.mode.k,pole,1)
        cc.domain==c.domain || error("causal audit requires the registered baseline domain")
        grid = SW.EccentricGrid(copy(base.r),cc.r_c,1:0,SW.mode_cache_key(cc),(;fixed_causal_grid=true))
        sol = timed(costs,"$(label)_$pole",Dict("pole"=>pole,"nodes"=>length(grid.r),
                                              "source_reused"=>true)) do
            SW.solve_mode(cc;source=base.source,grid,boundary_rows=SW.mode_boundary_rows(cc))
        end
        push!(solutions,sol)
    end
    solutions
end

"Radial gate (PREREG.json, unchanged) from the reported budget and source probe; secular (1,-1) is conditional."
radial_pass(r) = !r["conditional"] && sum(values(r["budget"]))<=1e-3 &&
    r["source_projection_error"]!==nothing && r["source_projection_error"]<=1e-8

function radial_case(a,p,e,m,k; nq=1024, source_reference=nothing, levels=(128,256,512))
    costs = Dict[]
    c,pert = radial_context(a,p,e,m,k,1e-3,1)
    projection_error = nothing
    if source_reference!==nothing
        ref = source_reference
        j = findfirst(==([m,k]),ref["modes"])
        j===nothing && error("mode absent from source reference")
        h,_ = source_table(c.disc,pert,ref["radii"],[c.mode],nq,ref["Nphi_base"],costs,"radial_source_probe")
        baseline = decode_value(ref["raw_orbit"][1])[:,j:j]
        projection_error = maximum(source_error(h,baseline,source_floors(c.disc,pert)))
        println("SOURCE_PROBE_ERROR ",projection_error); flush(stdout)
        checkpoint = get(ENV,"ECC_CONVERGENCE_CHECKPOINT_DIR","")
        if !isempty(checkpoint)
            writejson(joinpath(checkpoint,"source_probe.json"),Dict("Nq"=>nq,
                      "error"=>projection_error,"candidate"=>encode_value(h),
                      "reference"=>encode_value(baseline),"floors"=>source_floors(c.disc,pert)))
        end
    end
    radii = collect(range(c.domain...;length=257))
    if c.r_c!==nothing
        width = SW.corotation_width(c)
        append!(radii,[c.r_c+16width*sinh(x)/sinh(4) for x in range(-4,4;length=513)
                       if first(radii)<=c.r_c+16width*sinh(x)/sinh(4)<=last(radii)])
        sort!(unique!(radii))
    end
    finest = last(levels)
    grids = Any[]
    for n in levels
        push!(grids,radial_solve(c,pert,n,nq,costs,"grid_$n"))
        if n==first(levels) && c.r_c!==nothing
            probes = causal_audit(c,pert,first(grids),costs,"coarse_causal")
            checkpoint = get(ENV,"ECC_CONVERGENCE_CHECKPOINT_DIR","")
            if !isempty(checkpoint)
                writejson(joinpath(checkpoint,"coarse_causal.json"),
                          Dict("comparison_radii"=>radii,"p"=>p,"e"=>e,
                               "raw"=>[Dict("r"=>s.r,"y"=>encode_value(s.y)) for s in probes],
                               "costs"=>costs[end-1:end],
                               "errors"=>[response_error(probes[i],probes[i+1],radii) for i in 1:2],
                               "scope"=>"conditional coarse-grid causal sensitivity; not a limit certificate"))
            end
        end
    end
    grid_error = [response_error(grids[i],grids[i+1],radii) for i in 1:2]
    # Reuse exact source values for the operator-only cache audit. Grid, source
    # and domain refinements use the end-to-end solve_mode(c,pert) entry point.
    g = SW.EccentricGrid(copy(last(grids).r),c.r_c,1:0,SW.mode_cache_key(c),(;audit=true))
    cached = timed(costs,"operator_cache",Dict("nodes"=>length(g.r))) do
        SW.solve_mode(c;source=last(grids).source,grid=g,cache=SW.mode_operator_cache(c,g),
                      boundary_rows=SW.mode_boundary_rows(c))
    end
    cache_error = response_error(cached,last(grids),radii)
    ancillary = [radial_solve(c,pert,finest,nq,costs,"ancillary_$nanc";nanc) for nanc in (128,256)]
    push!(ancillary,last(grids))
    ancillary_error = [response_error(ancillary[i],ancillary[i+1],radii) for i in 1:2]
    poles = causal_audit(c,pert,last(grids),costs,"fine_causal")
    pole_error = [response_error(poles[i],poles[i+1],radii) for i in 1:2]
    pole_mesh = c.r_c===nothing ? last(poles) :
        radial_solve(last(poles).context,pert,2finest,nq,costs,"lowest_pole_mesh_$(2finest)")
    pole_mesh_error = response_error(last(poles),pole_mesh,radii)
    edges = Any[last(grids)]
    for edge in (2,3)
        cc,pp = radial_context(a,p,e,m,k,1e-3,edge)
        # Fixed cells per WKB radian: actual node count grows with domain width.
        push!(edges,radial_solve(cc,pp,finest,nq,costs,"edge_$edge"))
    end
    edge_error = [response_error(edges[i],edges[i+1],radii) for i in 1:2]
    sourcefine = radial_solve(c,pert,finest,2nq,costs,"response_source_$(2nq)")
    source_change = response_error(last(grids),sourcefine,radii)
    reflection = [Dict("edge"=>i,"samples"=>[
        reflection_sample(s,first(s.r),false),reflection_sample(s,last(s.r),true),
        reflection_sample(s,radii[17],false),reflection_sample(s,radii[end-16],true)])
        for (i,s) in enumerate(edges)]
    diagnostics = [Dict("axis_level"=>label,"torque"=>SW.total_torque(s),
                        "power"=>SW.total_power(s),"residual"=>s.residual,
                        "bc_residual"=>s.bc_residual,"nodes"=>length(s.r),
                        "corotation"=>s.context.r_c,"pole_width"=>SW.corotation_width(s.context))
                   for (label,s) in vcat(collect(zip(["grid_$n" for n in levels],grids)),
                                          collect(zip(("pole_1","pole_2","pole_3"),poles)),
                                          collect(zip(("edge_1","edge_2","edge_3"),edges)))]
    budget = Dict("grid"=>maximum(last(grid_error)),"source"=>maximum(source_change),
                  "cache"=>maximum(cache_error),"ancillary"=>maximum(last(ancillary_error)),
                  "pole"=>maximum(last(pole_error)),"pole_mesh"=>maximum(pole_mesh_error),
                  "edge"=>maximum(last(edge_error)))
    raw_axes = Dict(label=>[Dict("r"=>s.r,"y"=>encode_value(s.y)) for s in solutions]
                   for (label,solutions) in (("grid",grids),("ancillary",ancillary),
                       ("pole",poles),("edge",edges),("cache",[cached,last(grids)]),
                       ("source",[last(grids),sourcefine]),("pole_mesh",[last(poles),pole_mesh])))
    result = Dict("a"=>a,"p"=>p,"e"=>e,"m"=>m,"k"=>k,"Nq"=>nq,
         "grid_error"=>grid_error,"grid_order"=>[observed_order(grid_error[1][j],grid_error[2][j]) for j in 1:2],
         "cache_error"=>cache_error,"ancillary_error"=>ancillary_error,
         "pole_error"=>pole_error,"edge_error"=>edge_error,"source_error"=>source_change,
         "pole_mesh_error"=>pole_mesh_error,
         "reflection"=>reflection,"diagnostics"=>diagnostics,"budget"=>budget,
         "raw_axes"=>raw_axes,
         "source_projection_error"=>projection_error,
         "comparison_radii"=>radii,
         "costs"=>costs,"conditional"=>m==1 && k==-1)
    result["pass"] = radial_pass(result)
    result
end

function radial_study(out, selected, nq)
    fixtures = vcat([(0.0,30.0,e) for e in (0.05,0.3,0.6)],
                    [(a,p,0.3) for p in (10.0,100.0) for a in (-0.9,0.9)])
    cases = [(a,p,e,m,k) for (a,p,e) in fixtures for (m,k) in ((0,1),(2,1),(8,40),(8,-40),(1,-1))]
    for i in (selected===nothing ? eachindex(cases) : [selected])
        path = joinpath(out,@sprintf("case_%02d.json",i))
        isfile(path) && continue
        checkpoint = joinpath(out,@sprintf("case_%02d_stages",i))
        mkpath(checkpoint)
        ENV["ECC_CONVERGENCE_CHECKPOINT_DIR"] = checkpoint
        println("RADIAL_CASE ",i," ",cases[i]); flush(stdout)
        result = try
            a,p,e,m,k = cases[i]
            reference_dir = joinpath(dirname(out),"gr-quadrature-convergence")
            references = [json_read(read(joinpath(reference_dir,f),String)) for f in readdir(reference_dir)
                          if startswith(f,"case_") && endswith(f,".json")]
            j = findfirst(r->r["a"]==a && r["p"]==p && r["e"]==e,references)
            ref = j===nothing ? nothing : references[j]
            radial_case(a,p,e,m,k;nq,source_reference=ref)
        catch err
            Dict("case"=>collect(cases[i]),"pass"=>false,"error"=>sprint(showerror,err,catch_backtrace()))
        end
        writejson(path,result)
        println("RADIAL_GATE ",i," ",result["pass"]); flush(stdout)
        GC.gc()
    end
end

function main(args)
    length(args)>=2 || error("usage: convergence_study.jl source|radial EVIDENCE_DIR [case-index] [Nq]")
    kind,out = args[1:2]
    isfile(joinpath(out,"PREREG.json")) || error("write the preregistration before running")
    BLAS.set_num_threads(1)
    if kind == "source"
        cases = [(a,p,e) for e in (0.05,0.3,0.6) for p in (10.0,30.0,100.0) for a in (-0.9,0.0,0.9)]
        selected = length(args)>=3 ? [parse(Int,args[3])] : eachindex(cases)
        for i in selected
            path = joinpath(out,@sprintf("case_%02d.json",i))
            isfile(path) && continue
            a,p,e = cases[i]
            println("CASE ",i," ",cases[i]); flush(stdout)
            result = try
                source_case(a,p,e)
            catch err
                # Retain failures for the reviewer; never erase a failed cell.
                Dict("a"=>a,"p"=>p,"e"=>e,"pass"=>false,
                     "error"=>sprint(showerror,err,catch_backtrace()))
            end
            writejson(path,result)
            println("SOURCE_GATE ",i," ",result["pass"]," strict_v1=",get(result,"pass_v1_strict",false)); flush(stdout)
            GC.gc()
        end
    elseif kind == "radial"
        radial_study(out,length(args)>=3 ? parse(Int,args[3]) : nothing,
                     length(args)>=4 ? parse(Int,args[4]) : 1024)
    else
        error("unknown study $kind")
    end
end

end # module

abspath(PROGRAM_FILE)==(@__FILE__) && EccentricConvergenceStudy.main(ARGS)

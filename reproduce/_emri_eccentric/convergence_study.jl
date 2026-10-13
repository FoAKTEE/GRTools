"""Float64 convergence measurements using the existing eccentric source and BVP owners.

Run with the AccretionDisc project: convergence_study.jl source|radial EVIDENCE_DIR.
The study's preregistration must already exist. Per-case files support restart without replacing evidence.
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

function fixture(a, p, e; pole=1e-3, closure=Adiabatic, qTemp=0.0)
    rref = p/(1-e^2)
    d = SW.DiscModel(;χ=a,r_ref=rref,h_ref=0.06,pSigma=1.5,qTemp,
                     closure,evolution=:fixed,εp=pole)
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
function response_error(a,b,radii; reference=b)
    floors = (1e-12*1e-5/b.context.disc.r_ref,
              1e-12*1e-5/sqrt(b.context.disc.r_ref))
    weights = [(radii[min(i+1,end)]-radii[max(i-1,1)])/2 for i in eachindex(radii)]
    map(1:2) do j
        va = local_interp(radii,a.r,vec(a.y[j,:]);npts=4)
        vb = local_interp(radii,b.r,vec(b.y[j,:]);npts=4)
        vr = reference===b ? vb : local_interp(radii,reference.r,vec(reference.y[j,:]);npts=4)
        sqrt(sum(weights.*abs2.(va-vb))) /
            max(sqrt(sum(weights.*abs2.(vr))),sqrt(sum(weights))*floors[j])
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

function radial_context(a,p,e,m,k,pole,edge; closure=Adiabatic, qTemp=0.0)
    d,pert = fixture(a,p,e;pole,closure,qTemp)
    ilo,peri,apo = ((1.05,0.6,2.0),(1.01,0.4,3.0),(1.001,0.25,4.0),
                    (1.001,0.0,6.0),(1.001,0.0,8.0))[edge]
    domain = (max(ilo*SW.r_isco(d),peri*p/(1+e)),apo*p/(1-e))
    SW.ModeContext(d,pert,m,k;closure,domain),pert
end

function radial_solve(c,pert,n,nq,costs,label;nanc=512,cache=false,source=nothing,nphi=nothing,
                      grid_context=c)
    width = SW.corotation_width(c)
    mesh_width = SW.corotation_width(grid_context)
    fraction = mesh_width === nothing ? 0.01 : 16mesh_width/grid_context.r_c
    nw = max(1,n÷64) # grid.jl's nL/nR are cells per WKB radian, not node counts.
    sol = timed(costs,label,Dict("n"=>n,"Nq"=>nq,"n_anc"=>nanc,
                                "phase_cells_per_radian"=>nw,"shared_source"=>source!==nothing,
                                "pole"=>c.disc.εp,
                                "domain"=>collect(c.domain))) do
        grid = SW.build_grid(grid_context;nL=nw,nR=nw,nC=n+1,n_anc=nanc,σ_wc_frac=fraction)
        if grid_context !== c
            grid_context.domain == c.domain || throw(ArgumentError("pole meshes need the same domain"))
            grid = SW.EccentricGrid(grid.r,c.r_c,1:0,SW.mode_cache_key(c),(;fixed_causal_grid=true))
        end
        ops = cache ? SW.mode_operator_cache(c,grid) : nothing
        source === nothing ?
            SW.solve_mode(c,pert;grid,cache=ops,bank_options=(;Nq=nq,Nphi=nphi,blocksize=1),
                          radial_batch_size=16) :
            SW.solve_mode(c;source,grid,cache=ops,boundary_rows=SW.mode_boundary_rows(c))
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
function causal_audit(c,pert,base,costs,label; divisors=(2,4))
    c.r_c===nothing && return [base,base,base]
    solutions = [base]
    for divisor in divisors
        pole = c.disc.εp/divisor
        d,_ = fixture(pert.orbit.a,pert.orbit.p,pert.orbit.e;pole,
                      closure=c.closure,qTemp=c.disc.qTemp)
        cc = SW.ModeContext(d,pert,c.mode;closure=c.closure,domain=c.domain)
        grid = SW.EccentricGrid(copy(base.r),cc.r_c,1:0,SW.mode_cache_key(cc),(;fixed_causal_grid=true))
        sol = timed(costs,"$(label)_$pole",Dict("pole"=>pole,"nodes"=>length(grid.r),
                                              "source_reused"=>true)) do
            SW.solve_mode(cc;source=base.source,grid,boundary_rows=SW.mode_boundary_rows(cc))
        end
        push!(solutions,sol)
    end
    solutions
end

"Relative observable change with no cancellation-dependent acceptance floor."
relative_change(a,b,reference) = iszero(a-b) ? zero(abs(a-b)) : abs(a-b)/abs(reference)

"""Maximum discrete-row residual divided by its largest individual term.

`matrices` supplies `(Cmid, Smid)` pairs, possibly lazily. The four matrix-times-state
terms are exactly those of `SW.midpoint_system`, plus the right-hand side `dr*S`.
This measures the solved algebraic equations without dividing a small difference
of adjacent states by the mesh spacing. Boundary rows retain their separate gate.
"""
function midpoint_row_backward_error(r, y, matrices)
    size(y)==(2,length(r)) || throw(DimensionMismatch("expected two states per node"))
    length(matrices)==length(r)-1 || throw(DimensionMismatch("expected one matrix/source per cell"))
    worst = zero(real(eltype(y)))
    for (i,(C,S)) in enumerate(matrices)
        dr = r[i+1]-r[i]; h = dr/2
        for a in 1:2
            terms = ntuple(4) do j
                c = (j+1)÷2
                delta = a==c ? one(eltype(y)) : zero(eltype(y))
                isodd(j) ? -(delta+h*C[a,c])*y[c,i] : (delta-h*C[a,c])*y[c,i+1]
            end
            rhs = dr*S[a]
            residual = abs(sum(terms)-rhs)
            scale = max(maximum(abs,terms),abs(rhs))
            worst = max(worst,iszero(scale) ? residual : residual/scale)
        end
    end
    worst
end

"Row integrity of an existing solution; no solve, source or precision is changed."
function solution_row_backward_error(s)
    matrices = (SW.master_matrices(s.context,(s.r[i]+s.r[i+1])/2;
        hmodes=s.source((s.r[i]+s.r[i+1])/2),landau=s.params.landau)
        for i in 1:length(s.r)-1)
    midpoint_row_backward_error(s.r,s.y,matrices)
end

"Store response and physical observables so independent verifiers can recompute budgets."
function solution_record(s; comparison_radii=nothing, row_residual=false)
    work = SW.gravitational_work(s)
    torque = SW.total_torque(s;density=getproperty.(work,:torque))
    power = SW.total_power(s;density=getproperty.(work,:power))
    r = comparison_radii===nothing ? s.r : comparison_radii
    y = comparison_radii===nothing ? s.y :
        permutedims(hcat([local_interp(r,s.r,vec(s.y[j,:]);npts=4) for j in 1:2]...))
    record = Dict("r"=>r,"y"=>encode_value(y),"torque"=>torque,"power"=>power,
         "comparison_samples"=>comparison_radii!==nothing,
         "residual"=>s.residual,"bc_residual"=>s.bc_residual,"nodes"=>length(s.r),
         "domain"=>collect(s.context.domain),
         "corotation"=>s.context.r_c,"pole_width"=>SW.corotation_width(s.context),
         "pole"=>s.context.disc.εp)
    if s.context.closure == LocallyIsothermalOptionA
        record["closure"] = string(s.closure)
        record["qTemp"] = s.context.disc.qTemp
    end
    row_residual && (record["row_backward_error"] = solution_row_backward_error(s))
    record
end

const RADIAL_AXES = ("grid","source","cache","ancillary","pole","pole_mesh","edge")

"""Verified second-order grid error on at least four production-domain levels.

All differences share the finest reference norm. Both consecutive orders of
both state components and both observables must qualify before any grid relief
is used. Identically zero axisymmetric torque needs no error estimate. An
unqualified ladder retains its full last differences; it never certifies an
extrapolated production value or grants relief to another refinement axis.
"""
function grid_richardson(records, levels; domain, rref, m)
    length(records)==length(levels)>=4 || throw(ArgumentError("need at least four grid levels"))
    all(>(0),levels) && all(levels[2:end].==2 .* levels[1:end-1]) ||
        throw(ArgumentError("grid levels must double"))
    all(r->get(r,"domain",nothing)==collect(domain),records) ||
        throw(ArgumentError("all grid levels must use the production domain"))
    radii=Float64.(last(records)["r"])
    all(r->r["r"]==radii,records) || throw(ArgumentError("grid samples must agree"))
    mock(r)=(;r=radii,y=decode_value(r["y"]),context=(;disc=(;r_ref=rref)))
    solutions=mock.(records)
    profiles=[response_error(solutions[j],solutions[j+1],radii;reference=last(solutions))
              for j in 1:length(records)-1]
    differences=Dict("component_$j"=>getindex.(profiles,j) for j in 1:2)
    for q in ("torque","power")
        differences[q]=[relative_change(records[j][q],records[j+1][q],last(records)[q])
                        for j in 1:length(records)-1]
    end
    orders=Dict(k=>[observed_order(d[j],d[j+1]) for j in 1:length(d)-1]
                for (k,d) in differences)
    zero_torque=m==0 && all(r->iszero(r["torque"]),records)
    integrity=all(r->0<=get(r,"row_backward_error",Inf)<=1e-9 &&
                      0<=get(r,"bc_residual",Inf)<=1e-9,records)
    qualified=integrity && all(k->(k=="torque" && zero_torque) ||
        all(p->p!==nothing && isfinite(p) && abs(p-2)<=0.1,orders[k]),keys(orders))
    estimates=Dict(k=>qualified && !(k=="torque" && zero_torque) ?
        last(d)/(2^last(orders[k])-1) : last(d) for (k,d) in differences)
    Dict("qualified"=>qualified,"levels"=>collect(levels),"domain"=>collect(domain),
         "orders"=>orders,"differences"=>differences,"estimates"=>estimates,
         "integrity"=>integrity,"structural_zero_torque"=>zero_torque)
end

"""v4 joint gate: retain every v3 control, changing only a qualified grid term.

`independent_grid` contains separately measured grid estimates (full differences
when their order is unverified). Its entries are two profile components and
relative torque/power errors with the independent grid reference. The corner
certificate was measured on the production domain. Causal, pole-mesh and final
ISCO differences remain raw and the larger independent/corner term is charged.
LI controls use the production observable scale. A preregistered absolute power
reference may replace that scale near a cancellation; it never changes profiles.
"""
function radial_v4_gate(independent, joint, certificate, independent_grid, grid_reference, corner_reference;
                        power_reference=nothing)
    ind=deepcopy(independent)
    li=get(ind,"closure",nothing)=="LocallyIsothermalOptionA"
    power_reference===nothing || (li && isfinite(power_reference) && power_reference>0) ||
        throw(ArgumentError("an absolute power reference must be positive, finite and LI"))
    reference=power_reference===nothing ? corner_reference :
        merge(corner_reference,Dict("power"=>power_reference))
    profiles=deepcopy(joint["profile"])
    terms=deepcopy(joint["absolute_observable_terms"])
    for j in 1:2
        key="component_$j"
        factor=certificate["qualified"] ? 2^last(certificate["orders"][key])-1 : 1
        profiles["grid"][j]=max(profiles["grid"][j]/factor,certificate["estimates"][key])
    end
    if certificate["qualified"]
        ind["budget"]["grid"]=maximum(independent_grid["profile"])
        for q in ("torque","power")
            ind["observable_budget"][q]["grid"]=independent_grid[q]
            terms[q]["grid"]=max(independent_grid[q]*abs(grid_reference[q]),
                certificate["estimates"][q]*abs(corner_reference[q]))
        end
    end
    if li
        for q in ("torque","power"), axis in RADIAL_AXES
            absolute=ind["observable_budget"][q][axis]*abs(grid_reference[q])
            ind["observable_budget"][q][axis]=relative_change(absolute,zero(absolute),reference[q])
        end
    end
    budgets=Dict(q=>relative_change(sum(values(terms[q])),0.,reference[q])
                 for q in ("torque","power"))
    pass=radial_pass(ind) && get(joint,"integrity",false) &&
        isempty(get(joint,"owner_grid_rejections",["missing"])) &&
        get(joint,"causal_covered",false) &&
        all(q->all(k->haskey(terms[q],k) && 0<=terms[q][k],(RADIAL_AXES...,"isco")),("torque","power")) &&
        all(k->haskey(profiles,k) && length(profiles[k])==2 &&
            all(x->0<=x<=1e-3,profiles[k]),(RADIAL_AXES...,"isco")) &&
        all(v->0<=v<=1e-3,values(budgets)) &&
        all(k->!isempty(joint["source_probe"][k]) && all(x->0<=x<=1,joint["source_probe"][k]),("primary","sensitivity")) &&
        certificate["integrity"]
    result=Dict("pass"=>pass,"independent_pass"=>radial_pass(ind),"profile"=>profiles,
         "observable_budget"=>budgets,"absolute_observable_terms"=>terms,
         "grid_relief"=>certificate["qualified"])
    if li
        result["independent_observable_budget"]=ind["observable_budget"]
        result["observable_reference"]=Dict(q=>abs(reference[q]) for q in ("torque","power"))
    end
    result
end

"PREREG_radial_v2: profile gates per axis, summed uncertainty for EACH physical observable."
function radial_pass(r)
    diagnostics = r["diagnostics"]
    v3 = get(r,"protocol","")=="PREREG_radial_v3"
    residual_key = v3 ? "row_backward_error" : "residual"
    required_nq = v3 && get(r,"e",0)==0.6 && abs(get(r,"k",0))==40 ? 2048 : 1024
    base = findfirst(d->get(d,"axis_level","")=="grid_1",diagnostics)
    # An expanded domain can acquire a corotation absent from the baseline.
    # Its causal error was not measured by the baseline pole ladder.
    base!==nothing && (diagnostics[base]["corotation"]!==nothing ||
                       all(d->d["corotation"]===nothing,diagnostics)) &&
        !r["conditional"] && r["Nq"]>=required_nq &&
        all(k -> r[k]!==nothing && !isempty(_flat(r[k])) && all(x->0<=x<=1,_flat(r[k])),
            ("source_budget","source_budget_sensitivity")) &&
        all(k -> 0<=get(r["budget"],k,Inf)<=1e-3,RADIAL_AXES) &&
        all(q -> all(k->0<=get(r["observable_budget"][q],k,Inf)<=1e-3,RADIAL_AXES) &&
                 sum(r["observable_budget"][q][k] for k in RADIAL_AXES)<=1e-3,("torque","power")) &&
        all(d -> 0<=get(d,residual_key,Inf)<=1e-9 && 0<=d["bc_residual"]<=1e-9,diagnostics)
end

function radial_case(a,p,e,m,k; nq=1024, source_reference=nothing, levels=(128,256,512),
                     source=nothing, source_fine=nothing, source_probe=nothing, nphi=nothing,
                     compact=false, edge_levels=(1,2,3), pole=1e-3, source_control_nq=2nq,
                     protocol="PREREG_radial_v2", closure=Adiabatic, qTemp=0.0)
    protocol in ("PREREG_radial_v2","PREREG_radial_v3") || throw(ArgumentError("unknown radial protocol"))
    v3 = protocol=="PREREG_radial_v3"
    source_control_nq isa Integer && source_control_nq>0 && source_control_nq!=nq ||
        throw(ArgumentError("source control Nq must be a distinct positive integer"))
    required_nphi = e==0.6 ? 16384 : 8192
    nphi = something(nphi,required_nphi)
    nphi>=required_nphi || throw(ArgumentError("radial source Nphi floor is below the qualified resolution"))
    costs = Dict[]
    c,pert = radial_context(a,p,e,m,k,pole,first(edge_levels);closure,qTemp)
    projection_error = nothing
    projection_budgets = nothing; sensitivity_budgets = nothing
    if source_reference!==nothing
        ref = source_reference
        j = findfirst(==([m,k]),ref["modes"])
        j===nothing && error("mode absent from source reference")
        h,_ = source_table(c.disc,pert,ref["radii"],[c.mode],nq,ref["Nphi_base"],costs,"radial_source_probe")
        baseline = decode_value(ref["raw_orbit"][1])[:,j:j]
        projection_error = maximum(source_error(h,baseline,source_floors(c.disc,pert)))
        projection_budgets = vec(source_budget(h,baseline,source_floors(c.disc,pert)))
        sensitivity_budgets = vec(source_budget(h,baseline,test_unit_floors(pert)))
        println("SOURCE_PROBE_ERROR ",projection_error); flush(stdout)
        checkpoint = get(ENV,"ECC_CONVERGENCE_CHECKPOINT_DIR","")
        if !isempty(checkpoint)
            writejson(joinpath(checkpoint,"source_probe.json"),Dict("Nq"=>nq,
                      "error"=>projection_error,"candidate"=>encode_value(h),
                      "reference"=>encode_value(baseline),"floors"=>source_floors(c.disc,pert)))
        end
    end
    if source_probe !== nothing
        projection_error = source_probe["strict"]
        projection_budgets = source_probe["primary"]
        sensitivity_budgets = source_probe["sensitivity"]
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
        push!(grids,radial_solve(c,pert,n,nq,costs,"grid_$n";source,nphi))
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
    grid_error = [response_error(grids[i],grids[i+1],radii) for i in 1:length(grids)-1]
    # Reuse the identical source for the operator-only cache audit. A supplied
    # qualified bank also permits reuse across the independently changed meshes.
    g = SW.EccentricGrid(copy(last(grids).r),c.r_c,1:0,SW.mode_cache_key(c),(;audit=true))
    cached = timed(costs,"operator_cache",Dict("nodes"=>length(g.r))) do
        SW.solve_mode(c;source=last(grids).source,grid=g,cache=SW.mode_operator_cache(c,g),
                      boundary_rows=SW.mode_boundary_rows(c))
    end
    cache_error = response_error(cached,last(grids),radii)
    ancillary = [radial_solve(c,pert,finest,nq,costs,"ancillary_$nanc";nanc,source,nphi) for nanc in (128,256)]
    push!(ancillary,last(grids))
    ancillary_error = [response_error(ancillary[i],ancillary[i+1],radii) for i in 1:2]
    pole_diagnostic = causal_audit(c,pert,last(grids),costs,"fine_causal";
                                  divisors=v3 ? (2,4,8) : (2,4))
    poles = pole_diagnostic[1:3]
    pole_error = [response_error(poles[i],poles[i+1],radii) for i in 1:2]
    pole_mesh = c.r_c===nothing ? last(poles) :
        radial_solve(last(poles).context,pert,finest,nq,costs,"lowest_pole_mesh_$finest";source,nphi)
    pole_mesh_error = response_error(last(poles),pole_mesh,radii)
    edges = Any[last(grids)]
    for edge in edge_levels[2:end]
        cc,pp = radial_context(a,p,e,m,k,pole,edge;closure,qTemp)
        # Fixed cells per WKB radian: actual node count grows with domain width.
        push!(edges,radial_solve(cc,pp,finest,nq,costs,"edge_$edge";source,nphi))
    end
    edge_error = [response_error(edges[i],edges[i+1],radii) for i in 1:2]
    # source_fine supplies the comparison bank; its actual resolution may be lower.
    sourceother = radial_solve(c,pert,finest,source_control_nq,costs,
        "response_source_$source_control_nq";source=source_fine,nphi)
    source_pair = source_control_nq<nq ? [sourceother,last(grids)] : [last(grids),sourceother]
    source_change = response_error(source_pair...,radii)
    reflection = [Dict("edge"=>i,"samples"=>[
        reflection_sample(s,first(s.r),false),reflection_sample(s,last(s.r),true),
        reflection_sample(s,radii[17],false),reflection_sample(s,radii[end-16],true)])
        for (i,s) in enumerate(edges)]
    budget = Dict("grid"=>maximum(last(grid_error)),"source"=>maximum(source_change),
                  "cache"=>maximum(cache_error),"ancillary"=>maximum(last(ancillary_error)),
                  "pole"=>maximum(last(pole_error)),"pole_mesh"=>maximum(pole_mesh_error),
                  "edge"=>maximum(last(edge_error)))
    axes = (("grid",grids),("ancillary",ancillary),("pole",poles),("edge",edges),
            ("cache",[cached,last(grids)]),("source",source_pair),
            ("pole_mesh",[last(poles),pole_mesh]))
    records = IdDict()
    for (_,solutions) in axes, s in solutions
        get!(records,s) do
            solution_record(s;comparison_radii=compact ? radii : nothing,row_residual=v3)
        end
    end
    raw_axes = Dict(label=>[records[s] for s in solutions] for (label,solutions) in axes)
    diagnostics = [merge(Dict("axis_level"=>label*"_$i"),
                    Dict(k=>v for (k,v) in records[s] if k ∉ ("r","y")))
                   for (label,solutions) in axes for (i,s) in enumerate(solutions)]
    reference = records[last(grids)]
    observable_budget = Dict(q=>Dict(label=>relative_change(solutions[end-1][q],
                              solutions[end][q],reference[q])
                              for (label,solutions) in raw_axes) for q in ("torque","power"))
    result = Dict("a"=>a,"p"=>p,"e"=>e,"m"=>m,"k"=>k,"Nq"=>nq,"Nphi_floor"=>nphi,
         "edge_levels"=>collect(edge_levels),"baseline_pole"=>pole,
         "source_control_Nq"=>sort([nq,source_control_nq]),
         "grid_error"=>grid_error,"grid_order"=>[[observed_order(grid_error[i][j],grid_error[i+1][j]) for j in 1:2] for i in 1:length(grid_error)-1],
         "cache_error"=>cache_error,"ancillary_error"=>ancillary_error,
         "pole_error"=>pole_error,"edge_error"=>edge_error,"source_error"=>source_change,
         "pole_mesh_error"=>pole_mesh_error,
         "reflection"=>reflection,"diagnostics"=>diagnostics,"budget"=>budget,
         "raw_axes"=>raw_axes,
         "source_projection_error"=>projection_error,
         "source_budget"=>projection_budgets,"source_budget_sensitivity"=>sensitivity_budgets,
         "observable_budget"=>observable_budget,"protocol"=>protocol,
         "legacy_profile_sum"=>sum(values(budget)),
         "comparison_radii"=>radii,
         "costs"=>costs,"conditional"=>m==1 && k==-1)
    if v3 && c.r_c!==nothing
        fourth = solution_record(last(pole_diagnostic);comparison_radii=radii,row_residual=true)
        result["causal_diagnostic"] = vcat(raw_axes["pole"],[fourth])
    end
    if c.closure == LocallyIsothermalOptionA
        result["closure"] = string(c.closure)
        result["qTemp"] = c.disc.qTemp
    end
    result["pass"] = radial_pass(result)
    result
end

"LI OPTION A branch: retain the complete closure-dispatched operator and boundary rows."
li_radial_case(a,p,e,m,k; qTemp, protocol="PREREG_radial_v3", kwargs...) =
    radial_case(a,p,e,m,k;kwargs...,closure=LocallyIsothermalOptionA,qTemp,protocol)

"Last LI grid cycle: fixed qTemp-aware ladder from the preregistered prior finest n."
function li_radial_levels(previous_n, qTemp; cell)
    qTemp in (0,1) || throw(ArgumentError("LI radial registration covers qTemp=0,1"))
    n=previous_n*(qTemp==1 ? (cell==11 ? 4 : 2) : (cell==11 ? 2 : 1))
    (n÷8,n÷4,n÷2,n)
end

function radial_study(out, selected, nq)
    nq>=1024 || throw(ArgumentError("PREREG_radial_v2 requires Nq >= 1024"))
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
    prereg = kind == "radial" ? "PREREG_radial_v2.json" : "PREREG.json"
    isfile(joinpath(out,prereg)) || error("write $prereg before running")
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

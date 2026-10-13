using Test
include(joinpath(@__DIR__,"../../../reproduce/_emri_eccentric/convergence_study.jl"))
using .EccentricConvergenceStudy
const ECS = EccentricConvergenceStudy

@testset "Convergence norms and direct source cache contract" begin
    d,p = ECS.fixture(0.0,30.0,0.3)
    floors = ECS.source_floors(d,p)
    @test floors[1] == 1e-12*2p.mass/d.r_ref
    @test floors[2] == floors[1]/d.r_ref
    @test floors[5] ≈ floors[1]*d.r_ref^2
    z = NamedTuple{ECS.SW.ECCENTRIC_SOURCE_JETS}(ntuple(_->0.0+0im,14))
    a = fill(z,2,1); b = copy(a)
    b[2,1] = merge(z,(htt=2floors[1],))
    @test ECS.source_error(a,b,floors)[1,1] == 1
    @test all(iszero,ECS.source_error(a,a,floors))
    @test ECS.observed_order(0.04,0.01) ≈ 2
    @test ECS.observed_order(0.0,0.0) === nothing
    @test ECS.observed_order(0.01,0.04) ≈ -2
    doc = ECS.json_read(ECS.json_write(ECS.encode_value(Dict("modes"=>ECS.SOURCE_MODES))))
    @test doc["modes"][1]==[0,1]
    mode = ECS.HarmonicMode(0,1,ECS.orbit_frequency(p.orbit,0,1))
    b1 = ECS.SW.source_azimuth_bank(d,p,[45.0,46.0];Nq=32,mmax=0)
    b2 = ECS.SW.source_azimuth_bank(d,p,[45.0,45.5,46.0];Nq=32,mmax=0,
                                   Nphi=maximum(b1.nphi),blocksize=1)
    s1 = ECS.SW.mode_source(b1,mode); s2 = ECS.SW.mode_source(b2,mode)
    @test s1(45.0) == s2(45.0)
    @test s1(46.0) == s2(46.0)
    @test_throws ArgumentError s1(45.5)
    @test eltype(b1.q) === Float64
end

@testset "Discrete row integrity uses the assembled terms and fails closed" begin
    # Cross-check the streaming measurement against the independent sparse assembler.
    r = [1.0,1.0+1e-10,1.25,2.0]
    C = [[2.0+im -3.0; 5.0 1.0-im] for _ in 1:3]
    S = [[0.5+0im, -0.3im] for _ in 1:3]
    y = ComplexF64[1 1+1e-10 2 3; 2 2-1e-10 1 -1]
    A,b = ECS.SW.midpoint_system(r,C,S,(1,0),0im,(0,1),0im)
    u = vec(y)
    expected = maximum(2:length(b)-1) do i
        terms = A[i,:].*u
        abs(sum(terms)-b[i])/max(maximum(abs,terms),abs(b[i]))
    end
    @test ECS.midpoint_row_backward_error(r,y,zip(C,S)) ≈ expected rtol=1e-14
    @test ECS.midpoint_row_backward_error(r,7y,zip(C,7 .* S)) ≈ expected rtol=1e-14
    z = zeros(ComplexF64,2,4)
    @test ECS.midpoint_row_backward_error(r,z,zip(C,[zeros(ComplexF64,2) for _ in 1:3]))==0
    bad = copy(y); bad[1,2]=NaN
    @test isnan(ECS.midpoint_row_backward_error(r,bad,zip(C,S)))
    @test_throws DimensionMismatch ECS.midpoint_row_backward_error(r,y[:,1:3],zip(C,S))
    @test_throws DimensionMismatch ECS.midpoint_row_backward_error(r,y,zip(C[1:2],S[1:2]))
    # A representational derivative residual is large even for a row solved to roundoff.
    tiny = [1.0,1.0+1e-10]
    cm = [zeros(ComplexF64,2,2)]; sm = [ones(ComplexF64,2)]
    smooth = ComplexF64[1 1+(tiny[2]-tiny[1]); 2 2+(tiny[2]-tiny[1])]
    @test ECS.midpoint_row_backward_error(tiny,smooth,zip(cm,sm))<1e-15
end

@testset "v3 gates row error independently of derivative diagnostics and source resolution" begin
    axes7 = ECS.RADIAL_AXES
    r = Dict{String,Any}("protocol"=>"PREREG_radial_v3","conditional"=>false,"Nq"=>1024,
        "source_budget"=>[0.5],"source_budget_sensitivity"=>[0.5],
        "budget"=>Dict(k=>1e-4 for k in axes7),
        "observable_budget"=>Dict(q=>Dict(k=>1e-4 for k in axes7) for q in ("torque","power")),
        "diagnostics"=>[Dict{String,Any}("axis_level"=>"grid_1","corotation"=>10.0,
            "residual"=>1e-4,"row_backward_error"=>1e-12,"bc_residual"=>1e-14)])
    @test ECS.radial_pass(r)
    for value in (1.01e-9,NaN,Inf,-1.0)
        rr=deepcopy(r); rr["diagnostics"][1]["row_backward_error"]=value
        @test !ECS.radial_pass(rr)
    end
    rr=deepcopy(r); delete!(rr["diagnostics"][1],"row_backward_error")
    @test !ECS.radial_pass(rr)
    rr=deepcopy(r); rr["protocol"]="PREREG_radial_v2"
    @test !ECS.radial_pass(rr)
    rr=deepcopy(r); rr["e"]=0.6; rr["k"]=40
    @test !ECS.radial_pass(rr)
    rr["Nq"]=2048
    @test ECS.radial_pass(rr)
    rr["k"]=-40; rr["Nq"]=1024
    @test !ECS.radial_pass(rr)
end

@testset "Row backward error keeps the right-hand side in its scale (lander L-G46)" begin
    # A small state under a large source: dr*S dominates every interior row, so the largest
    # term entering the row is the right-hand side, not a matrix-times-state term.
    r = [1.0,1.25,2.0]
    C = [[2.0+im -3.0; 5.0 1.0-im] for _ in 1:2]
    S = [[40.0+0im, -25.0im] for _ in 1:2]
    y = 1e-3 .* ComplexF64[1 2 3; 2 1 -1]
    A,b = ECS.SW.midpoint_system(r,C,S,(1,0),0im,(0,1),0im)
    u = vec(y)
    rows = 2:length(b)-1
    @test all(i -> abs(b[i]) > 10maximum(abs,A[i,:].*u), rows)
    expected = maximum(rows) do i
        terms = A[i,:].*u
        abs(sum(terms)-b[i])/max(maximum(abs,terms),abs(b[i]))
    end
    @test ECS.midpoint_row_backward_error(r,y,zip(C,S)) ≈ expected rtol=1e-14
end

@testset "Complete small diagnostic path with an independent source reference" begin
    d,p = ECS.fixture(0.0,30.0,0.05)
    mode = ECS.HarmonicMode(0,1,ECS.orbit_frequency(p.orbit,0,1))
    radii = [30/1.05,30.0,30/0.95]
    bank = ECS.SW.source_azimuth_bank(d,p,radii;Nq=4,mmax=0,blocksize=1)
    h = ECS.SW.source_modes(bank,[mode])
    ref = Dict("modes"=>[[0,1]],"radii"=>radii,"Nphi_base"=>maximum(bank.nphi),
               "raw_orbit"=>[ECS.encode_value([collect(values(x)) for x in h])])
    # Rebuild the reference at the common Nphi floor used by the diagnostic probe.
    bank = ECS.SW.source_azimuth_bank(d,p,radii;Nq=4,mmax=0,blocksize=1,Nphi=ref["Nphi_base"])
    h = ECS.SW.source_modes(bank,[mode])
    ref["raw_orbit"] = [ECS.encode_value([collect(values(x)) for x in h])]
    result = ECS.radial_case(0.0,30.0,0.05,0,1;nq=4,source_reference=ref,levels=(16,32,64))
    @test result["source_projection_error"]==0
    @test result["source_control_Nq"]==[4,8]
    @test all(iszero,result["cache_error"])
    @test length(result["grid_error"])==2
    @test all(isfinite,values(result["budget"]))
end

@testset "Response comparison and end-to-end diagnostic slice" begin
    d,p = ECS.fixture(0.0,30.0,0.05)
    c = ECS.SW.ModeContext(d,p,0,1;closure=ECS.Adiabatic,domain=(45.0,46.0))
    r1 = collect(range(45.0,46.0;length=5))
    r2 = collect(range(45.0,46.0;length=9))
    f(r) = (r-45)^2+im*(r-45)
    a = (;r=r1,y=permutedims(hcat(f.(r1),2f.(r1))),context=c)
    b = (;r=r2,y=permutedims(hcat(f.(r2),2f.(r2))),context=c)
    @test maximum(ECS.response_error(a,b,r2))<1e-14
    costs = Dict[]
    sol = ECS.radial_solve(c,p,8,32,costs,"test";nanc=32)
    @test sol.residual<1e-10
    @test sol.bc_residual<1e-12
    @test all(isfinite,sol.y)
    @test last(costs)["nodes"]==length(sol.r)
    @test maximum(ECS.response_error(sol,sol,sol.r))==0
    ref = ECS.reflection_sample(sol,first(sol.r),false)
    @test ref["incoming_fraction"]<1e-6
end

@testset "Causal prescription changes the operator while retaining the source" begin
    c,p = ECS.radial_context(0.0,30.0,0.3,2,1,4e-3,2)
    @test c.r_c!==nothing
    costs = Dict[]
    base = ECS.radial_solve(c,p,16,4,costs,"causal_test";nanc=32)
    sols = ECS.causal_audit(c,p,base,costs,"causal_test")
    @test [s.context.disc.εp for s in sols]==[4e-3,2e-3,1e-3]
    @test all(s->s.context.domain==c.domain,sols)
    @test all(s->s.source===base.source && s.r==base.r,sols)
    @test all(s->all(isfinite,s.y),sols)
    @test maximum(abs,sols[end].y-base.y)>0
end

@testset "PREREG_v2 source budget, sensitivity unit and fourth-order audit" begin
    d,p = ECS.fixture(0.0,30.0,0.3)
    floors = ECS.source_floors(d,p)
    z = NamedTuple{ECS.SW.ECCENTRIC_SOURCE_JETS}(ntuple(_->0.0+0im,14))
    ref = fill(z,2,1); new = copy(ref)
    # A change equal to the absolute floor on a zero reference is exactly at the budget.
    new[1,1] = merge(z,(htt=floors[1]+0im,))
    @test ECS.source_budget(new,ref,floors)[1,1] == 1
    @test all(iszero,ECS.source_budget(ref,ref,floors))
    # Sub-floor harmonic: the strict diagnostic fails, the DESIGN budget passes (no 1e-20 demand).
    ref[2,1] = merge(z,(hzz=0.5floors[9]+0im,)); new = copy(ref)
    new[2,1] = merge(z,(hzz=0.51floors[9]+0im,))
    @test ECS.source_error(new,ref,floors)[1,9] ≈ 0.01
    @test ECS.source_budget(new,ref,floors)[1,9] ≈ 0.01/(1+0.5e-8)
    # Relative part: a 2e-8 relative change of a large reference exceeds the budget.
    big = 1e12floors[5]
    ref[2,1] = merge(z,(hphiphi=big+0im,)); new = copy(ref)
    new[2,1] = merge(z,(hphiphi=big*(1+2e-8)+0im,))
    @test 1.9 < ECS.source_budget(new,ref,floors)[1,5] < 2.0
    # The budget uses the radial maximum of the reference, not the value at the worst radius.
    ref = fill(z,2,1); ref[1,1] = merge(z,(htr=1e4floors[11]+0im,)); new = copy(ref)
    new[2,1] = merge(z,(htr=1e-4floors[11]+0im,))
    @test ECS.source_budget(new,ref,floors)[1,11] ≈ 1e-4/(1e-4+1)
    # ... and of the reference, not of the candidate.
    ref = fill(z,2,1); ref[1,1] = merge(z,(htr=1e12floors[11]+0im,)); new = fill(z,2,1)
    new[1,1] = merge(z,(htr=0.5e12floors[11]+0im,))
    @test ECS.source_budget(new,ref,floors)[1,11] ≈ 0.5e12/(1e4+1)
    # Sensitivity unit = the committed source-bank test unit, smaller by 2(1-e^2) = 2p/r_ref.
    s = ECS.test_unit_floors(p)
    q = p.orbit.p
    @test s[1] == 1e-12*p.mass/q
    @test s[2] == s[1]/q && s[3] ≈ s[1]*q && s[5] ≈ s[1]*q^2 && s[13] ≈ s[1]*q
    # Ratio 2(q/r_ref)^(1-w) per component (w = power of r_ref in its floor): 2(1-e^2) for htt,
    # 2 for htphi, 2/(1-e^2) for hphiphi, one more factor q/r_ref for radial jets.
    w = (0,1,2,0,0,0,1)
    @test all(floors[j]/s[j] ≈ 2(q/d.r_ref)^(1-w[(j+1)÷2]+(iseven(j) ? 1 : 0)) for j in 1:14)
    @test floors[1]/s[1] ≈ 2(1-0.3^2)
    # Five-point stencil: exact on quartics, fourth order on a smooth oscillation.
    f(x) = 3x^4-2x^3+x-7
    @test isapprox(ECS.fourth_difference(f(1.0-0.2),f(1.0-0.1),f(1.0+0.1),f(1.0+0.2),0.1),
                   12-6+1;atol=1e-12)
    g(x) = sin(7x); x0 = 0.3
    err(h) = abs(ECS.fourth_difference(g(x0-2h),g(x0-h),g(x0+h),g(x0+2h),h)-7cos(7x0))
    @test 15 < err(0.02)/err(0.01) < 17
    A = [g(x0-0.02) g(x0-0.01); g(x0-0.04) g(x0-0.02)]
    @test ECS.fourth_difference(A, A, 2A, 2A, 0.1) ≈ (A-8A+16A-2A)/1.2 rtol=1e-14
end

@testset "source_case end to end: reported gates recompute from its raw tables" begin
    res = ECS.source_case(0.0,30.0,0.3;nq=8,mode_list=((0,1),(2,1)))
    @test res["Nq_levels"] == [8,16,32]
    @test res["protocol"] == "PREREG_AMENDMENT_01+PREREG_v2"
    floors, sens = res["source_floors"], res["sensitivity_floors"]
    raw(x) = ECS.decode_value(x)
    orbit = raw.(res["raw_orbit"]); azim = raw.(res["raw_azimuth"])
    same(rows, M) = all(rows[j] == M[j,:] for j in axes(M,1))
    for i in 1:2
        @test same(res["orbit_budget"][i], ECS.source_budget(orbit[i],orbit[i+1],floors))
        @test same(res["orbit_budget_sensitivity"][i], ECS.source_budget(orbit[i],orbit[i+1],sens))
        @test same(res["orbit_error"][i], ECS.source_error(orbit[i],orbit[i+1],floors))
        @test same(res["azimuth_budget"][i], ECS.source_budget(azim[i],azim[i+1],floors))
        @test same(res["azimuth_budget_sensitivity"][i], ECS.source_budget(azim[i],azim[i+1],sens))
    end
    # Exact direct samples: common radii of nested tables are identical.
    @test all(all(iszero,r) for b in res["cache_budget"] for r in b)
    # Derivative audits from the stored samples and the analytic jet at p of the Nq=8 table.
    tabs = Dict(x["fraction"] => (x["step"], raw(x["table"])) for x in res["raw_derivatives"])
    @test sort(collect(keys(tabs))) == [0.005,0.01,0.02,0.04,0.08]
    @test all(tabs[f][1] == f*0.018*(30.0/(1-0.3^2)) for f in keys(tabs))
    keys14 = ECS.SW.ECCENTRIC_SOURCE_JETS
    v(f,i,j,n) = tabs[f][2][i,j][String(keys14[2n-1])]
    jet(j,n) = orbit[1][2,j][2n]
    den(j,n) = max(abs(jet(j,n)),floors[2n])
    for (l,(f,f2)) in enumerate(((0.02,0.04),(0.01,0.02),(0.005,0.01)))
        h = tabs[f][1]
        for j in 1:2, n in 1:7
            d4 = (v(f2,1,j,n)-8v(f,1,j,n)+8v(f,2,j,n)-v(f2,2,j,n))/(12h)
            @test res["derivative_error_v2"][l][j][n] ≈ abs(d4-jet(j,n))/den(j,n) rtol=1e-12
        end
    end
    for (l,f) in enumerate((0.08,0.04,0.02)), j in 1:2, n in 1:7
        d2 = (v(f,2,j,n)-v(f,1,j,n))/(2tabs[f][1])
        @test res["derivative_error"][l][j][n] == abs(d2-jet(j,n))/den(j,n)
    end
    @test res["pass"] == ECS.prereg_v2_pass(res)
    @test res["pass_v1_strict"] == ECS.prereg_v1_pass(res)
    @test length(ECS._flat(res["orbit_budget"])) == 2*2*14 && length(ECS._flat(res["cache_budget"])) == 2*2*14
    @test all(haskey(res,k) for k in ECS.V2_BUDGETS)
end

@testset "PREREG_v2 and strict gates are pure functions of the reported arrays" begin
    lv(x) = [[fill(x,14) for _ in 1:2] for _ in 1:2]           # two levels of 2 modes x 14
    dv(xs) = [[fill(x,7) for _ in 1:2] for x in xs]            # three derivative levels
    base() = merge(Dict{String,Any}(k=>lv(0.5) for k in ECS.V2_BUDGETS),
                   Dict{String,Any}("orbit_error"=>lv(0.0),"azimuth_error"=>lv(0.0),"cache_error"=>lv(0.0),
                                    "derivative_error_v2"=>dv((0.5,0.1,1e-4)),
                                    "derivative_error"=>dv((0.5,0.1,1e-4))))
    r = base(); @test ECS.prereg_v2_pass(r) && ECS.prereg_v1_pass(r)
    r = base(); r["orbit_budget"][2][1][3] = 1.0; @test ECS.prereg_v2_pass(r)
    six = ("orbit_budget","orbit_budget_sensitivity","azimuth_budget","azimuth_budget_sensitivity",
           "cache_budget","cache_budget_sensitivity")
    @test Set(ECS.V2_BUDGETS) == Set(six)
    for k in six                       # each axis and both floors are gated
        r = base(); r[k][end][end][end] = 1+1e-9; @test !ECS.prereg_v2_pass(r)
    end
    r = base(); r["derivative_error_v2"][end][2][7] = 2e-3; @test !ECS.prereg_v2_pass(r)
    r = base(); r["derivative_error_v2"][1][1][1] = 10.0; @test ECS.prereg_v2_pass(r)   # finest level only
    r = base(); r["derivative_error_v2"][end][1][1] = NaN; @test !ECS.prereg_v2_pass(r)
    for k in ("orbit_error","azimuth_error","cache_error")
        r = base(); r[k][1][2][14] = 2e-8; @test !ECS.prereg_v1_pass(r) && ECS.prereg_v2_pass(r)
    end
    r = base(); r["derivative_error"][end][1][1] = 2e-3; @test !ECS.prereg_v1_pass(r)
end

@testset "Response profile norm: trapezoid-weighted relative L2 against the reference" begin
    d,p = ECS.fixture(0.0,30.0,0.05)
    c = ECS.SW.ModeContext(d,p,0,1;closure=ECS.Adiabatic,domain=(45.0,46.0))
    radii = [45.0,45.1,45.15,45.5,46.0]                    # nonuniform trapezoid weights
    base(r) = (r-44)^2+0im
    dif(r) = 1e-3*(r-45)+0im
    prof(f,g) = (;r=radii,y=permutedims(hcat(f.(radii),g.(radii))),context=c)
    b = prof(base,r->2base(r)); a = prof(r->base(r)+dif(r),r->2base(r))
    w = [(radii[min(i+1,end)]-radii[max(i-1,1)])/2 for i in eachindex(radii)]
    l2(f) = sqrt(sum(w.*abs2.(f.(radii))))
    err = ECS.response_error(a,b,radii)
    @test err[1] ≈ l2(dif)/l2(base) rtol=1e-10
    @test err[2] == 0
    # normalised by the reference (second argument), not by the candidate
    @test ECS.response_error(b,a,radii)[1] ≈ l2(dif)/l2(r->base(r)+dif(r)) rtol=1e-10
    # a zero reference falls back to the preregistered response floor
    z = (;r=radii,y=zeros(ComplexF64,2,length(radii)),context=c)
    floor1 = 1e-12*1e-5/d.r_ref
    @test ECS.response_error(a,z,radii)[1] ≈ l2(r->base(r)+dif(r))/(sqrt(sum(w))*floor1) rtol=1e-10
end

@testset "Radial v2 gates keep profile and observable budgets separate" begin
    axes7 = ("grid","source","cache","ancillary","pole","pole_mesh","edge")
    base() = Dict{String,Any}("conditional"=>false,"Nq"=>1024,"source_projection_error"=>1.0,
        "source_budget"=>[0.5],"source_budget_sensitivity"=>[0.5],
        "budget"=>Dict(k=>2e-4 for k in axes7),
        "observable_budget"=>Dict(q=>Dict(k=>1e-4 for k in axes7) for q in ("torque","power")),
        "diagnostics"=>[Dict{String,Any}("axis_level"=>"grid_1","corotation"=>nothing,
                                        "residual"=>1e-12,"bc_residual"=>1e-14)])
    @test ECS.radial_pass(base()) # profile maxima SUM is 1.4e-3: it is diagnostic only
    for key in ("source_budget","source_budget_sensitivity")
        r=base(); r[key][1]=1.01; @test !ECS.radial_pass(r)
        r=base(); r[key][1]=NaN; @test !ECS.radial_pass(r)
        r=base(); r[key]=nothing; @test !ECS.radial_pass(r)
    end
    for key in axes7
        r=base(); r["budget"][key]=1.01e-3; @test !ECS.radial_pass(r)
        r=base(); delete!(r["budget"],key); @test !ECS.radial_pass(r)
        for q in ("torque","power")
            r=base(); r["observable_budget"][q][key]=5e-4; @test !ECS.radial_pass(r)
            r=base(); delete!(r["observable_budget"][q],key); @test !ECS.radial_pass(r)
            r=base(); r["observable_budget"][q][key]=-1e-4; @test !ECS.radial_pass(r)
        end
    end
    r=base(); r["conditional"]=true; @test !ECS.radial_pass(r)
    r=base(); r["Nq"]=512; @test !ECS.radial_pass(r)
    @test_throws ArgumentError ECS.radial_study("unused",1,512)
    @test_throws ArgumentError ECS.radial_case(0.0,30.0,0.05,0,1;nphi=4096)
    @test_throws ArgumentError ECS.radial_case(0.0,30.0,0.6,0,1;nphi=8192)
    r=base(); r["diagnostics"][1]["residual"]=2e-9; @test !ECS.radial_pass(r)
    r=base(); r["diagnostics"][1]["bc_residual"]=NaN; @test !ECS.radial_pass(r)
    r=base(); empty!(r["diagnostics"]); @test !ECS.radial_pass(r)
    r=base(); edge=copy(r["diagnostics"][1]);edge["axis_level"]="edge_2";edge["corotation"]=10.0
    push!(r["diagnostics"],edge); @test !ECS.radial_pass(r)
    r["diagnostics"][1]["corotation"]=10.0; @test ECS.radial_pass(r)
    @test ECS.relative_change(0.0,0.0,0.0)==0
    @test isinf(ECS.relative_change(1.0,0.0,0.0))
    @test ECS.relative_change(-3.0,2.0,-10.0)==0.5
end

@testset "Radial shared source retains C-66 rows and measured observables" begin
    c,p=ECS.radial_context(0.0,30.0,0.05,0,1,1e-3,1)
    narrow=ECS.SW.ModeContext(c.disc,p,0,1;closure=ECS.Adiabatic,domain=(45.0,46.0))
    bank=ECS.SW.source_radial_bank(c.disc,p,narrow.domain,[c.mode];Nq=8)
    source=ECS.SW.mode_source(bank,c.mode)
    costs=Dict[]
    shared=ECS.radial_solve(narrow,p,8,8,costs,"shared";source,nanc=32)
    direct=ECS.radial_solve(narrow,p,8,8,costs,"direct";nanc=32)
    @test maximum(ECS.response_error(shared,direct,shared.r))<1e-8
    @test shared.source===source
    @test last(costs)["Nq"]==8
    doc=ECS.solution_record(shared)
    @test doc["power"]==ECS.SW.total_power(shared)
    @test doc["torque"]==0
    @test doc["bc_residual"]<1e-9
    sampled=ECS.solution_record(shared;comparison_radii=[45.0,45.5,46.0])
    @test sampled["r"]==[45.0,45.5,46.0]
    @test sampled["nodes"]==length(shared.r)
    @test sampled["comparison_samples"] && !doc["comparison_samples"]
    @test sampled["power"]==doc["power"] && sampled["torque"]==doc["torque"]
    values=ECS.decode_value(sampled["y"])
    @test values[1,:] ≈ ECS.local_interp(sampled["r"],shared.r,vec(shared.y[1,:]);npts=4)
    @test values[2,:] ≈ ECS.local_interp(sampled["r"],shared.r,vec(shared.y[2,:]);npts=4)
end


@testset "Coarser source control uses the qualified baseline as profile reference" begin
    result = ECS.radial_case(0.0,30.0,0.05,0,1;
        nq=8,source_control_nq=4,levels=(8,16,32),compact=true)
    @test result["source_control_Nq"]==[4,8]
    @test last(result["raw_axes"]["source"]) === last(result["raw_axes"]["grid"])
    control = only(filter(x->x["axis_level"]=="response_source_4",result["costs"]))
    @test control["Nq"]==4
    @test all(isfinite,result["source_error"])
    @test !result["pass"] # Small diagnostic source is not production-qualified.
    for invalid in (0,-4,8,4.5)
        @test_throws ArgumentError ECS.radial_case(0.0,30.0,0.05,0,1;
            nq=8,source_control_nq=invalid)
    end
end

@testset "Reported observable budget recomputes from the stored records (lander L-G39)" begin
    # A corotation mode: nC = n+1 refines the grid even at these tiny levels, and the
    # causal ladder changes the operator, so every three-level axis carries distinct records.
    result = ECS.radial_case(0.0,30.0,0.3,2,1;nq=8,levels=(8,16,32),compact=true)
    rax = result["raw_axes"]; ref = last(rax["grid"])
    @test all(k->length(rax[k])==3,("grid","ancillary","pole","edge"))
    for q in ("torque","power"), (axis,records) in rax
        x,y = records[end-1][q],records[end][q]
        @test result["observable_budget"][q][axis] == (x==y ? 0.0 : abs(x-y)/abs(ref[q]))
    end
    @test all(q->result["observable_budget"][q]["grid"]>0 && result["observable_budget"][q]["pole"]>0,
              ("torque","power"))
    @test first(rax["grid"])["power"] != ref["power"]
end

@testset "Expanded production-domain edge levels (lander L-G39)" begin
    for (edge,apo) in ((4,6.0),(5,8.0))
        c,_ = ECS.radial_context(0.0,30.0,0.3,2,1,1e-3,edge)
        @test c.domain == (1.001*ECS.SW.r_isco(c.disc),apo*30.0/(1-0.3))
    end
end

@testset "v3 records four causal points while budgeting the registered three-point ladder" begin
    c,p=ECS.radial_context(0.0,30.0,0.3,2,1,4e-6,1)
    largest,_=ECS.radial_context(0.0,30.0,0.3,2,1,4e-6,3)
    domain=(largest.domain[1],largest.domain[2])
    banks=[ECS.SW.source_radial_bank(c.disc,p,domain,[c.mode];Nq=nq,Nphi=8192,method=:fft)
           for nq in (8,16)]
    sources=[ECS.SW.mode_source(b,c.mode) for b in banks]
    r=ECS.radial_case(0.0,30.0,0.3,2,1;nq=8,levels=(8,16,32),compact=true,
        source=sources[1],source_fine=sources[2],pole=4e-6,protocol="PREREG_radial_v3")
    @test [x["pole"] for x in r["causal_diagnostic"]]==[4e-6,2e-6,1e-6,5e-7]
    @test [x["pole"] for x in r["raw_axes"]["pole"]]==[4e-6,2e-6,1e-6]
    @test all(x->isfinite(x["row_backward_error"]) && haskey(x,"residual"),r["diagnostics"])
    @test r["protocol"]=="PREREG_radial_v3"
    @test !r["pass"] # This deliberately small source is not production qualified.
    for q in ("power","torque")
        a,b=r["raw_axes"]["pole"][2:3]
        @test r["observable_budget"][q]["pole"] == ECS.relative_change(a[q],b[q],r["raw_axes"]["grid"][end][q])
    end
    @test only(filter(x->startswith(x["axis_level"],"lowest_pole_mesh"),r["costs"]))["n"]==32
end

@testset "v4 production-domain verified-order grid certificate" begin
    domain=[6.,80.]; radii=[6.,20.,80.]; levels=[8,16,32,64]
    function ladder(order; m=2)
        [Dict("domain"=>domain,"r"=>radii,
              "y"=>ECS.encode_value(ComplexF64[1 2 3; 3 2 1].*(1+n^(-order))),
              "torque"=>m==0 ? 0. : 2+n^(-order),"power"=>3+2n^(-order),
              "row_backward_error"=>1e-13,"bc_residual"=>0.) for n in levels]
    end
    good=ladder(2.)
    cert=ECS.grid_richardson(good,levels;domain,rref=30.,m=2)
    @test cert["qualified"]
    @test all(v->all(p->p≈2,v),values(cert["orders"]))
    @test all(k->cert["estimates"][k]≈cert["differences"][k][end]/3,keys(cert["estimates"]))
    @test !ECS.grid_richardson(ladder(1.),levels;domain,rref=30.,m=2)["qualified"]
    @test ECS.grid_richardson(ladder(2.;m=0),levels;domain,rref=30.,m=0)["qualified"]
    for key in ("power","torque")
        bad=deepcopy(good)
        for (j,n) in enumerate(levels); bad[j][key]=1+n^(-1.); end
        @test !ECS.grid_richardson(bad,levels;domain,rref=30.,m=2)["qualified"]
    end
    for component in 1:2
        bad=deepcopy(good)
        for (j,n) in enumerate(levels)
            y=ECS.decode_value(bad[j]["y"]); y[component,:].=1+n^(-1.)
            bad[j]["y"]=ECS.encode_value(y)
        end
        @test !ECS.grid_richardson(bad,levels;domain,rref=30.,m=2)["qualified"]
    end
    bad=deepcopy(good); bad[1]["domain"]=[6.006,40.]
    @test_throws ArgumentError ECS.grid_richardson(bad,levels;domain,rref=30.,m=2)
    @test_throws ArgumentError ECS.grid_richardson(good[2:end],levels[2:end];domain,rref=30.,m=2)
    @test_throws ArgumentError ECS.grid_richardson(good,[8,16,32,128];domain,rref=30.,m=2)
    # Owner 2026-10-12: a solve-domain label cannot certify an unsampled tail.
    for samples in ([6.,20.,60.], [6.001,20.,80.], [6.,90.,80.])
        truncated=deepcopy(good)
        for rec in truncated; rec["r"]=samples; end
        @test_throws ArgumentError ECS.grid_richardson(truncated,levels;domain,rref=30.,m=2)
    end
    # The LI corner windows start at the finest ISCO cutoff, charged separately on the isco axis.
    cutoff=deepcopy(good); for rec in cutoff; rec["r"]=[6*1.000001,20.,80.]; end
    @test ECS.grid_richardson(cutoff,levels;domain,rref=30.,m=2)["qualified"]
    early=deepcopy(good); for rec in early; rec["r"]=[6.00001,20.,80.]; end
    @test_throws ArgumentError ECS.grid_richardson(early,levels;domain,rref=30.,m=2)
    bad=deepcopy(good); bad[2]["power"]=NaN
    @test !ECS.grid_richardson(bad,levels;domain,rref=30.,m=2)["qualified"]
    bad=deepcopy(good); bad[2]["row_backward_error"]=2e-9
    @test !ECS.grid_richardson(bad,levels;domain,rref=30.,m=2)["qualified"]
    # One good final order does not rescue the earlier first-order ratio.
    bad=deepcopy(good); bad[1]["power"]=bad[2]["power"]+2*(bad[2]["power"]-bad[3]["power"])
    @test !ECS.grid_richardson(bad,levels;domain,rref=30.,m=2)["qualified"]
    # Lander L-G50 (test-only addition): the owner's |p - 2| <= 0.1 window is pinned at both edges.
    @test ECS.grid_richardson(ladder(1.92),levels;domain,rref=30.,m=2)["qualified"]
    @test ECS.grid_richardson(ladder(2.08),levels;domain,rref=30.,m=2)["qualified"]
    @test !ECS.grid_richardson(ladder(1.85),levels;domain,rref=30.,m=2)["qualified"]
    @test !ECS.grid_richardson(ladder(2.15),levels;domain,rref=30.,m=2)["qualified"]
end

@testset "v4 keeps raw causal, ISCO, source and integrity gates" begin
    axes7=ECS.RADIAL_AXES
    ind=Dict{String,Any}("protocol"=>"PREREG_radial_v3","conditional"=>false,"Nq"=>1024,
        "source_budget"=>[0.5],"source_budget_sensitivity"=>[0.5],
        "budget"=>Dict(k=>1e-4 for k in axes7),
        "observable_budget"=>Dict(q=>Dict(k=>1e-4 for k in axes7) for q in ("torque","power")),
        "diagnostics"=>[Dict("axis_level"=>"grid_1","corotation"=>10.,
            "row_backward_error"=>1e-12,"bc_residual"=>0.)])
    ind["budget"]["grid"]=2.4e-3
    joint=Dict{String,Any}("profile"=>Dict(k=>[1e-4,1e-4] for k in (axes7...,"isco")),
        "absolute_observable_terms"=>Dict(q=>Dict(k=>1e-4 for k in (axes7...,"isco")) for q in ("torque","power")),
        "integrity"=>true,"owner_grid_rejections"=>[],"causal_covered"=>true,
        "source_probe"=>Dict("primary"=>[.5],"sensitivity"=>[.5]))
    joint["profile"]["grid"]=[2.4e-3,1e-4]
    cert=Dict("qualified"=>true,"integrity"=>true,"estimates"=>Dict("component_1"=>8e-4,"component_2"=>1e-4,"power"=>1e-4,"torque"=>1e-4))
    cert["orders"]=Dict("component_1"=>[2.,2.],"component_2"=>[2.,2.])
    ig=Dict("profile"=>[8e-4,1e-4],"power"=>1e-4,"torque"=>1e-4)
    ref=Dict("power"=>1.,"torque"=>1.)
    @test ECS.radial_v4_gate(ind,joint,cert,ig,ref,ref)["pass"]
    c=deepcopy(cert);c["qualified"]=false
    @test !ECS.radial_v4_gate(ind,joint,c,ig,ref,ref)["pass"]
    for axis in ("pole","pole_mesh","isco")
        j=deepcopy(joint);j["absolute_observable_terms"]["power"][axis]=4e-4
        @test !ECS.radial_v4_gate(ind,j,c,ig,ref,ref)["pass"]
        @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    end
    for key in ("integrity","causal_covered")
        j=deepcopy(joint);j[key]=false
        @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    end
    j=deepcopy(joint);j["owner_grid_rejections"]=["failed control"]
    @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    j=deepcopy(joint);j["source_probe"]["sensitivity"]=[1.01]
    @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    j=deepcopy(joint);j["profile"]["isco"]=[1.1e-3,0.]
    @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    j=deepcopy(joint);delete!(j["absolute_observable_terms"]["power"],"isco")
    @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    c=deepcopy(cert);c["integrity"]=false
    @test !ECS.radial_v4_gate(ind,joint,c,ig,ref,ref)["pass"]
    # Lander L-G50 (test-only addition): either response component alone fails the joint profile gate.
    for axis in ("edge","isco"), comp in 1:2
        j=deepcopy(joint);j["profile"][axis][comp]=1.1e-3
        @test !ECS.radial_v4_gate(ind,j,cert,ig,ref,ref)["pass"]
    end
    # W-CELL11 (test-only addition): the D1 production-domain estimate is charged in the joint grid profile even when
    # the corner pair alone passes (the term that kept adiabatic cell 11 a HOLE at n = 4096; survivor MC2).
    for comp in 1:2
        c=deepcopy(cert);c["estimates"]["component_$comp"]=1.1e-3
        @test !ECS.radial_v4_gate(ind,joint,c,ig,ref,ref)["pass"]
    end
    # A nearly cancelling base-domain observable must not tighten LI's production budget.
    small=Dict("power"=>1e-3,"torque"=>1e-3)
    i=deepcopy(ind); i["observable_budget"]["power"]["edge"]=0.03
    @test !ECS.radial_v4_gate(i,joint,cert,ig,small,ref)["pass"] # adiabatic unchanged
    i["closure"]="LocallyIsothermalOptionA"
    result=ECS.radial_v4_gate(i,joint,cert,ig,small,ref)
    @test result["pass"]
    @test result["independent_observable_budget"]["power"]["edge"]≈3e-5
    @test i["observable_budget"]["power"]["edge"]==0.03 # pure gate
    # The cancellation exception changes both independent and joint power budgets only.
    cancelled=Dict("power"=>0.01,"torque"=>1.)
    @test !ECS.radial_v4_gate(i,joint,cert,ig,small,cancelled)["pass"]
    bounded=ECS.radial_v4_gate(i,joint,cert,ig,small,cancelled;power_reference=1.)
    @test bounded["pass"]
    @test bounded["profile"]==result["profile"]
    @test bounded["observable_reference"]["power"]==1.
    for bad in (0.,-1.,Inf,NaN)
        @test_throws ArgumentError ECS.radial_v4_gate(i,joint,cert,ig,small,ref;power_reference=bad)
    end
    @test_throws ArgumentError ECS.radial_v4_gate(ind,joint,cert,ig,small,ref;power_reference=1.)
    for axis in ("pole","pole_mesh")
        j=deepcopy(joint); j["profile"][axis]=[1.1e-3,0.]
        @test !ECS.radial_v4_gate(i,j,cert,ig,small,ref;power_reference=1.)["pass"]
    end
end

@testset "Preregistered LI resolution has no post-run escalation" begin
    @test ECS.li_radial_levels(4096,0;cell=6)==(512,1024,2048,4096)
    @test ECS.li_radial_levels(4096,1;cell=6)==(1024,2048,4096,8192)
    @test ECS.li_radial_levels(4096,0;cell=11)==(1024,2048,4096,8192)
    @test ECS.li_radial_levels(4096,1;cell=11)==(2048,4096,8192,16384)
    @test_throws ArgumentError ECS.li_radial_levels(4096,0.5;cell=6)
end

@testset "LI radial controls preserve closure, temperature gradients and source reuse" begin
    for q in (0.0,1.0)
        c,p=ECS.radial_context(0.0,30.0,0.3,2,1,4e-6,3;
            closure=ECS.LocallyIsothermalOptionA,qTemp=q)
        @test c.closure == c.disc.closure == ECS.LocallyIsothermalOptionA
        @test c.disc.qTemp == q
        rows=ECS.SW.mode_boundary_rows(c)
        for (j,outward) in enumerate((false,true))
            r=c.domain[j]
            @test rows[j][1] == ECS.SW.li_boundary(c,r;outward).row
        end
        d0,_=ECS.fixture(0.0,30.0,0.3)
        # Reuse an adiabatic-built bank on a real LI response, including causal controls.
        outer=8p.orbit.p/(1-p.orbit.e)
        bank=ECS.SW.source_radial_bank(d0,p,(c.domain[1],outer),[c.mode];Nq=8,Nphi=8192)
        source=ECS.SW.mode_source(bank,c.mode)
        sol=ECS.radial_solve(c,p,16,8,Dict[],"li_fixture";source)
        @test sol.closure == sol.params.closure == ECS.LocallyIsothermalOptionA
        @test sol.source === source
        narrow,_=ECS.radial_context(0.0,30.0,0.3,2,1,1e-6,3;
            closure=ECS.LocallyIsothermalOptionA,qTemp=q)
        mesh=ECS.radial_solve(narrow,p,16,8,Dict[],"equal_n_mesh";source,grid_context=c)
        @test mesh.r==sol.r
        @test mesh.context.disc.εp==1e-6
        @test mesh.source===source
        @test mesh.closure==sol.closure
        for (j,node) in enumerate((1,length(sol.r)))
            terms=rows[j][1].*sol.y[:,node]
            @test abs(sum(terms))/maximum(abs,terms)<1e-9
        end
        @test abs(ECS.SW.li_boundary(c,last(c.domain);outward=true).mixing)>1e-8
        for s in ECS.causal_audit(c,p,sol,Dict[],"li_pole")
            @test s.context.disc.qTemp == q
            @test s.closure == ECS.LocallyIsothermalOptionA
            @test s.source === source
            @test ECS.solution_row_backward_error(s)<1e-9
            @test s.bc_residual<1e-9
        end
        r=ECS.li_radial_case(0.0,30.0,0.3,2,1;qTemp=q,nq=8,closure=ECS.Adiabatic,
            source_control_nq=16,levels=(8,16,32,64),compact=true,
            source,source_fine=source,pole=4e-6,edge_levels=(3,4,5))
        @test r["protocol"] == "PREREG_radial_v3"
        @test r["closure"] == "LocallyIsothermalOptionA"
        @test r["qTemp"] == q
        @test !r["pass"] # Tiny Nq and no source qualification cannot be certified.
        @test all(x->x["closure"]==r["closure"] && x["qTemp"]==q,r["diagnostics"])
        @test all(x->x["qTemp"]==q,r["causal_diagnostic"])
        @test only(filter(x->startswith(x["axis_level"],"lowest_pole_mesh"),r["costs"]))["n"]==64
        for obs in ("power","torque"), (axis,records) in r["raw_axes"]
            @test r["observable_budget"][obs][axis] == ECS.relative_change(
                records[end-1][obs],records[end][obs],r["raw_axes"]["grid"][end][obs])
        end
    end
    ad,_=ECS.radial_context(0.0,30.0,0.3,2,1,1e-3,3)
    @test ad.closure == ECS.Adiabatic && ad.disc.qTemp == 0.0
    # Lander L-G54: the fixture's own defaults stay adiabatic (survivor M8 of round 1).
    da,_=ECS.fixture(0.0,30.0,0.3)
    @test da.closure == ECS.Adiabatic && da.qTemp == 0.0
end

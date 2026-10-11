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
    c,p = ECS.radial_context(0.0,30.0,0.3,2,1,1e-3,1)
    @test c.r_c!==nothing
    costs = Dict[]
    base = ECS.radial_solve(c,p,16,4,costs,"causal_test";nanc=32)
    sols = ECS.causal_audit(c,p,base,costs,"causal_test")
    @test [s.context.disc.εp for s in sols]==[1e-3,5e-4,2.5e-4]
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

@testset "Radial gate is a pure function of the reported budget and source probe" begin
    axes7 = ("grid","source","cache","ancillary","pole","pole_mesh","edge")
    base() = Dict{String,Any}("conditional"=>false,"source_projection_error"=>1e-9,
                              "budget"=>Dict(k=>1e-4 for k in axes7))
    @test ECS.radial_pass(base())                                  # sum 7e-4
    r = base(); r["budget"]["edge"] = 5e-4; @test !ECS.radial_pass(r) # the SUM of axes is gated
    r = base(); r["source_projection_error"] = 2e-8; @test !ECS.radial_pass(r)
    r = base(); r["source_projection_error"] = nothing; @test !ECS.radial_pass(r)
    r = base(); r["conditional"] = true; @test !ECS.radial_pass(r)
    r = base(); r["budget"]["pole"] = NaN; @test !ECS.radial_pass(r)
end

isdefined(@__MODULE__, :eccentric_legacy_sample) || include("eccentric_test_support.jl")
using LinearAlgebra
using AccretionDisc.SpiralWaves: DiscModel, ModeContext, reconstruct_fluid, linear_stress,
    EccentricMasterInputs, EccentricOperatorCache, mode_cache_key, mode_operator_cache,
    master_operator
import AccretionDisc.SpiralWaves as ESW

@testset "ECC OPTION A runtime" begin
    for T in (Float64,BigFloat), q in (0,3//4,1), χ in (-3//5,2//5)
        d = DiscModel(χ=T(χ),r_ref=T(20),h_ref=T(3)/50,pSigma=T(3)/2,qTemp=T(q),
                      closure=:locally_isothermal,evolution=:fixed)
        c = ModeContext(d,HarmonicMode(3,-2,-T(3)/100);closure=:locally_isothermal)
        r=T(13); h=eccentric_test_source(r); y=(Complex{T}(1,2)/100,Complex{T}(-2,1)/50)
        x = master_inputs(c,r;hmodes=h,y,landau=false)
        C,S = master_matrices(c,r;hmodes=h,landau=false)
        f = reconstruct_fluid(c,r,y;hmodes=h,landau=false)
        @test c.closure == LocallyIsothermalOptionA
        @test eltype(C) == eltype(S) == Complex{T}
        @test all(isfinite,C) && all(isfinite,S)
        @test f.dp ≈ (x.p/x.eps)*f.deps rtol=10000eps(T)
        α=Complex{T}(2,3)
        C2,S2=master_matrices(c,r;hmodes=map(v -> α*v,h),landau=false)
        @test C2 == C
        @test S2 ≈ α*S rtol=10000eps(T)
        for field in (:htr,:hrphi)
            odd=merge(ESW._zero_source(c),NamedTuple{(field,)}((h[field],)))
            _,So=master_matrices(c,r;hmodes=odd,landau=false)
            @test maximum(abs,So)>0
        end
        # Independent perfect-fluid product rule in coordinate components.
        z=zero(T); u=[x.W/x.enu,z,(x.omt+x.enu*x.v/x.epsi)*x.W/x.enu,z]
        g=[-x.enu^2+x.epsi^2*x.omt^2 z -x.epsi^2*x.omt z;
           z x.emu^2 z z; -x.epsi^2*x.omt z x.epsi^2 z; z z z one(T)]
        H=[h.htt h.htr h.htphi z; h.htr h.hrr h.hrphi z;
           h.htphi h.hrphi h.hphiphi z; z z z h.hzz]
        du=[f.dut,f.dur,f.duphi_up,complex(z)]
        gi=inv(g)
        @test abs(2sum(u .* (g*du))+sum(u .* (H*u))) < 10000eps(T)
        expected=(f.deps+f.dp)*u*transpose(u)+(x.eps+x.p)*(du*transpose(u)+u*transpose(du))+
                 f.dp*gi-x.p*gi*H*gi
        base=(x.eps+x.p)*u*transpose(u)+x.p*gi
        current=expected*g+base*H
        stress=linear_stress(c,r,y;hmodes=h,landau=false)
        @test eltype(stress.dT)==Complex{T}
        @test stress.dT ≈ expected rtol=10000eps(T)
        @test collect(stress.dJt) ≈ current[1:3,1] rtol=10000eps(T)
        @test collect(stress.dJphi) ≈ current[1:3,3] rtol=10000eps(T)
        @test f.du_phi ≈ sum(g[3,:].*du)+sum(H[3,:].*u) rtol=10000eps(T)

        ad=DiscModel(d.reference;closure=:adiabatic,evolution=:fixed)
        ac=ModeContext(ad,c.mode;closure=:adiabatic)
        A,B=master_matrices(ac,r;hmodes=h,landau=false)
        if q==0
            @test C ≈ A rtol=10000eps(T)
            @test S ≈ B rtol=10000eps(T)
            af=reconstruct_fluid(ac,r,y;hmodes=h,landau=false)
            @test collect(values(f)) ≈ collect(values(af)) rtol=10000eps(T)
        else
            @test maximum(abs,C-A)>T(1e-8)
        end
        @test_throws ArgumentError linear_stress(ac,r,y;hmodes=h)
        # OPTION A does not consume an adiabatic-index input.
        nt=NamedTuple{fieldnames(typeof(x))}(Tuple(getfield(x,n) for n in fieldnames(typeof(x))))
        altered=merge(nt,(Gamma1=T(7),))
        @test ESW._context_operator(Val(c.closure),altered) == C
    end
end
println("[PASS] OPTION A operator/source/reconstruction, complex radial jets, precision and independent stress product rule")

@testset "ECC closure-bound operator cache and numerical solutions" begin
    domain=(12.0,16.0); rows=(((1,0),0),((1,0),0))
    d=DiscModel(χ=0.4,r_ref=20,h_ref=0.08,pSigma=1.5,qTemp=0.75,
                closure=:locally_isothermal,evolution=:fixed)
    c=ModeContext(d,HarmonicMode(3,-2,-0.03);closure=:locally_isothermal,domain)
    g=build_grid(c;nL=2,nR=2,n_anc=32)
    cache=mode_operator_cache(c,g)
    @test cache isa EccentricOperatorCache
    @test cache.key[1].closure == LocallyIsothermalOptionA
    @test_throws ArgumentError solve_mode(c;source=eccentric_test_source,grid=g)
    sol=solve_mode(c;source=eccentric_test_source,grid=g,boundary_rows=rows)
    cached=solve_mode(c;source=eccentric_test_source,grid=g,boundary_rows=rows,cache)
    @test cached.y == sol.y
    @test sol.closure == sol.params.closure == LocallyIsothermalOptionA
    @test sol.residual < 1e-10 && sol.bc_residual < 1e-12
    @test maximum(abs,sol.y)>0
    scaled=solve_mode(c;source=r -> map(z -> (2+3im)*z,eccentric_test_source(r)),
                      grid=g,boundary_rows=rows,cache)
    @test scaled.y ≈ (2+3im)*sol.y rtol=1e-11
    ad=DiscModel(d.reference;closure=:adiabatic,evolution=:fixed)
    ac=ModeContext(ad,c.mode;closure=:adiabatic,domain)
    ag=build_grid(ac;nL=2,nR=2,n_anc=32)
    @test mode_cache_key(c) != mode_cache_key(ac)
    @test length(Dict(mode_cache_key(c)=>:li,mode_cache_key(ac)=>:adiabatic))==2
    @test_throws ArgumentError solve_mode(ac;source=eccentric_test_source,grid=ag,boundary_rows=rows,cache)
    @test_throws ArgumentError solve_mode(c;source=eccentric_test_source,grid=g,boundary_rows=rows,cache,landau=false)
    fine=build_grid(c;nL=3,nR=3,n_anc=32)
    @test_throws ArgumentError solve_mode(c;source=eccentric_test_source,grid=fine,boundary_rows=rows,cache)
    # A dynamic axisymmetric LI mode also solves without a fabricated pole.
    axis=ModeContext(d,HarmonicMode(0,1,0.03);closure=:locally_isothermal,domain)
    gs=build_grid(axis;nL=1,nR=1,n_anc=32)
    axial=solve_mode(axis;source=eccentric_test_source,grid=gs,boundary_rows=rows)
    @test axis.r_c === nothing && axial.residual < 1e-10
    setprecision(BigFloat,192) do
        bd=DiscModel(χ=big"0.4",r_ref=big"20",h_ref=big"0.08",pSigma=big"1.5",qTemp=big"0.75",
                     closure=:locally_isothermal,evolution=:fixed)
        bc=ModeContext(bd,HarmonicMode(3,-2,big"-0.03");closure=:locally_isothermal,domain=BigFloat.(domain))
        bg=build_grid(bc;nL=1,nR=1,n_anc=16)
        bcache=mode_operator_cache(bc,bg)
        bs=solve_mode(bc;source=eccentric_test_source,grid=bg,boundary_rows=rows,cache=bcache)
        @test eltype(bs.y)==Complex{BigFloat} && bs.residual<big"1e-45"
        setprecision(BigFloat,256) do
            @test_throws ArgumentError solve_mode(bc;source=eccentric_test_source,grid=bg,boundary_rows=rows,cache=bcache)
        end
    end
end
println("[PASS] closure-tagged solutions/cache, cache reuse parity, cross-closure/grid/precision rejection, nonzero numerical LI solutions")

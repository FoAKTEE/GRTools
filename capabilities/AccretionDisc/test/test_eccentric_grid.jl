isdefined(@__MODULE__, :eccentric_legacy_sample) || include("eccentric_test_support.jl")
using AccretionDisc.SpiralWaves: DiscModel, ModeContext, EccentricGrid,
    EccentricModeSolution, corotation_width, mode_cache_key
import AccretionDisc.SpiralWaves as ESW

@testset "ECC mode-specific physical grids" begin
    d = DiscModel(χ=0.4,r_ref=20,h_ref=0.08,closure=:adiabatic,evolution=:fixed,εp=1e-3)
    domain = (8.0,35.0)
    for root in (nothing,8.0,8.0001,20.0,34.9999,35.0)
        ω = root === nothing ? -0.02 : 3Omega(d,root)
        c = ModeContext(d,HarmonicMode(3,-1,ω);closure=:adiabatic,domain)
        g = build_grid(c;nL=2,nR=2,nC=65,n_anc=64)
        @test g.r[1] == domain[1] && g.r[end] == domain[2]
        @test all(diff(g.r).>0)
        @test all(x -> domain[1]<=x<=domain[2],g.r)
        if root === nothing
            @test g.r_c === nothing && isempty(g.central)
            @test corotation_width(c) === nothing
            @test g.r == build_grid(c;nL=2,nR=2,nC=129,n_anc=64).r
        else
            @test g.r_c ≈ root rtol=1e-13
            @test !isempty(g.central) && c.r_c in g.r
            @test corotation_width(c)>0
        end
    end
    # A mode whose corotation is outside the requested domain has no core.
    c = ModeContext(d,HarmonicMode(3,1,3Omega(d,40.0));closure=:adiabatic,domain)
    @test c.r_c === nothing
    @test isempty(build_grid(c;nL=2,nR=2,nC=65,n_anc=64).central)
    axis = ModeContext(d,HarmonicMode(0,1,0.03);closure=:adiabatic,domain)
    @test build_grid(axis;nL=2,nR=2,n_anc=64).r_c === nothing
    @test_throws ArgumentError build_grid(c;nC=2)
    @test_throws ArgumentError build_grid(c;s_c=0)
    edge = ModeContext(d,HarmonicMode(3,0,3Omega(d,domain[1]));closure=:adiabatic,domain)
    ge = build_grid(edge;nL=2,nR=2,nC=65,n_anc=64)
    @test_throws ArgumentError solve_mode(edge;source=eccentric_test_source,grid=ge)
    explicit = solve_mode(edge;source=eccentric_test_source,grid=ge,
                          boundary_rows=(((1,0),0),((1,0),0)))
    @test explicit.residual < 1e-9

    # Cauchy-layer quadrature: a grid convergence check, not a torque certificate.
    c = ModeContext(d,HarmonicMode(3,0,3Omega(d,20.0));closure=:adiabatic,domain)
    width = corotation_width(c)
    errors = Float64[]
    for n in (33,65,129)
        g = build_grid(c;nL=2,nR=2,nC=n,n_anc=64,σ_wc_frac=0.01)
        r = g.r[g.central]; rc=c.r_c
        val = sum((r[i+1]-r[i])*width/(((r[i+1]+r[i])/2-rc)^2+width^2) for i in 1:length(r)-1)
        exact = atan((last(r)-rc)/width)-atan((first(r)-rc)/width)
        push!(errors,abs(val-exact))
    end
    println("Landau Cauchy-layer quadrature errors nC=33,65,129: ",errors)
    @test errors[1]/errors[2]>3 && errors[2]/errors[3]>3

    for ω in (-0.03,3Omega(d,20.0))
        c = ModeContext(d,HarmonicMode(3,0,ω);closure=:adiabatic,domain)
        g = build_grid(c;nL=2,nR=2,nC=129,n_anc=64)
        sol = solve_mode(c;source=eccentric_test_source,grid=g)
        @test sol isa EccentricModeSolution && sol.closure == Adiabatic
        @test sol.residual < 1e-9 && sol.bc_residual < 1e-12
        @test maximum(abs,sol.y)>0
        wrong = ModeContext(d,HarmonicMode(4,0,ω);closure=:adiabatic,domain)
        @test_throws ArgumentError solve_mode(wrong;source=eccentric_test_source,grid=g)
    end
    setprecision(BigFloat,192) do
        b = DiscModel(χ=big"0.4",r_ref=big"20",h_ref=big"0.08",closure=:adiabatic,evolution=:fixed)
        c = ModeContext(b,HarmonicMode(2,1,big"-0.03");closure=:adiabatic,domain=(big"12",big"16"))
        g = build_grid(c;nL=1,nR=1,nC=17,n_anc=16)
        sol = solve_mode(c;source=eccentric_test_source,grid=g)
        @test eltype(g.r) == BigFloat && eltype(sol.y) == Complex{BigFloat}
        @test sol.residual < big"1e-45"
    end
end
println("[PASS] in-domain/no-root/edge grids, physical boundaries, Landau-layer refinement, shared solver")

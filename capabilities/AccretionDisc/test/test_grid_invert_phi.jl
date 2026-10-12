using Test, AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
using NumericsCore: CubicSplineInterpolant, gl_nodes_weights
const SW = AccretionDisc.SpiralWaves

@testset "Cell-count inversion converges or throws" begin
    x,w=gl_nodes_weights(Float64,4)
    s=CubicSplineInterpolant([1.,2.,3.,4.],ones(4))
    @test SW._invert_phi(s,1.,2.,0.,1.,0.25,x,w,:R) === 1.25
    # An inconsistent cumulative count has no root in the segment. A small
    # coordinate bracket must not waive the estimator's residual check.
    @test_throws ErrorException SW._invert_phi(s,1.,2.,0.,2.,1.5,x,w,:R)
    # L-G46 invert_phi_probe.jl: cell12, final left ancillary segment, levels17/18.
    d=SW.DiscModel(;χ=0.,r_ref=30/(1-.6^2),h_ref=.06,pSigma=1.5,qTemp=0.,
        closure=Adiabatic,evolution=:fixed,εp=4e-6)
    p=SW.EccentricPerturber(EquatorialOrbit(0.,30.,.6);mass=1e-5,
        softening=FixedSoftening(.018d.r_ref))
    c=SW.ModeContext(d,p,2,1;closure=Adiabatic,domain=(max(1.001SW.r_isco(d),.25*30/1.6),4*30/.4))
    n=8192; sigma=(16SW.corotation_width(c)/c.r_c)*c.r_c
    spacing,slope=SW._central_edge_spacing(sigma,4.,n+1)
    edge=c.r_c-sigma
    reg=SW._wkb_region(c,:L,c.domain[1],edge,edge,spacing,slope,128,10.)
    N,nW=SW._solve_nW(reg,n÷64); S,phi=SW._phi_nodes(reg,nW)
    a,b=reg.r_anc[end-1:end]; pa,pb=phi[end-1:end]
    nodes=Float64[]
    for j in (17,18)
        try
            r=SW._invert_phi(S,a,b,pa,pb,Float64(j),reg.glx,reg.glw,:L)
            residual=abs(SW._gl_density_integral(S,a,r,reg.glx,reg.glw)-abs(j-pa))
            # The old level18 misses by ~9 cells; representability is ~1e-8 cells.
            @test residual < 1e-5
            push!(nodes,r)
        catch err
            @test err isa ErrorException
            @test occursin("grid: cell-count inversion did not converge",sprint(showerror,err))
        end
    end
    @test length(nodes)<2 || nodes[2]<nodes[1]
    # Lander L-G50 (test-only addition): the first 20 levels of the same segment. At the committed code every
    # level of the segment converges (measured for all 3974); levels 1-14 need the bisection fallback after
    # step 200 and many need the moving-node derivative term of the residual certificate, so this pins both.
    rs=[SW._invert_phi(S,a,b,pa,pb,Float64(j),reg.glx,reg.glw,:L) for j in 1:20]
    @test all(diff(rs) .< 0)
    @test all(j->abs(SW._gl_density_integral(S,a,rs[j],reg.glx,reg.glw)-abs(j-pa))<1e-5,1:20)
    # Non-finite target cannot be accepted by a bracket-width or iteration-cap exit.
    @test_throws DomainError SW._invert_phi(s,1.,2.,0.,1.,NaN,x,w,:R)
end

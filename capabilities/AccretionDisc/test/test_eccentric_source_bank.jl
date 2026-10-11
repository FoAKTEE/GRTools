using Test, LinearAlgebra, Printf
using AccretionDisc, Spacetime, NumericsCore
using AccretionDisc.ForcingModes
import AccretionDisc.SpiralWaves as ESW
const SOURCE_TEST_STAGE = get(ENV,"SOURCE_TEST_STAGE","all")
# Registration is intentionally external to this lane; direct inclusion works
# both before and after the module owner registers the production file.
const SOURCE_LEGACY_BG = ESW.DiscBackground(χ=-0.7,aEMRI=20.0,h0=0.06,b=0.3)
const SOURCE_LEGACY_BEFORE = [ESW.singular_field_mode(SOURCE_LEGACY_BG,m,r;N=1024)
                            for m in (0,1,7), r in (12.0,20.0,30.0)]
if !isdefined(ESW,:SourceAzimuthBank)
    Base.include(ESW,joinpath(@__DIR__,"..","src","SpiralWaves","eccentric_source.jl"))
end
include(joinpath(@__DIR__,"..","derivation","src","eccentric_source.jl"))
const R3_SOURCE = EccentricSource
const SOURCE_KEYS = ESW.ECCENTRIC_SOURCE_JETS
const SOURCE_PAIRS = ((1,1),(1,3),(3,3),(2,2),(4,4),(1,2),(2,3))

function source_fixture(a=0.5,p=20.0,e=0.3;soft=0.36)
    o = EquatorialOrbit(a,p,e)
    d = ESW.DiscModel(χ=a,r_ref=p,h_ref=0.06,closure=:adiabatic,evolution=:fixed)
    s = ESW.EccentricPerturber(o;softening=FixedSoftening(soft))
    d,s
end

function r3_source_inputs(p,q)
    s = q isa Real ? worldline_rates(p.orbit,q) : q
    gbl = Spacetime.kerr_metric(p.orbit.a,s.r,typeof(s.r)(pi)/2)
    g = Matrix(gbl[[1,2,4,3],[1,2,4,3]])
    g[4,4] = one(s.r)
    s,g,[s.ut,s.ur,s.uphi,zero(s.r)]
end

# Independent composite Gauss-Legendre evaluation of the full R3 tensor.
# This neither uses the runtime distance nor samples the FFT lattice.
function r3_azimuth_oracle(p,m,r,q;panels=16,order=48)
    s,g,u = r3_source_inputs(p,q)
    x,w = gl_nodes_weights(Float64,order)
    out = zeros(ComplexF64,14)
    for panel in 0:panels-1, j in eachindex(x)
        psi = (2panel+1+x[j])*pi/panels
        f = R3_SOURCE.field_jet(g,u,s.r,r,psi,p.softening.eps,p.mass)
        factor = w[j]*cis(-m*psi)/(2panels)
        for (c,(i,k)) in enumerate(SOURCE_PAIRS)
            out[2c-1] += factor*f.h[i,k]
            out[2c] += factor*f.h_r[i,k]
        end
    end
    NamedTuple{SOURCE_KEYS}(Tuple(out))
end

# Raw relative diagnostic, retaining tiny-coefficient cancellation information.
# Acceptance separately uses the design relative-plus-absolute budget.
function source_error(a,b,p,r)
    scales = (1.0, p.orbit.p, p.orbit.p^2, 1.0, 1.0, 1.0, p.orbit.p)
    maximum(eachindex(SOURCE_KEYS)) do i
        unit = p.mass*scales[(i+1)÷2]/p.orbit.p/(iseven(i) ? p.orbit.p : 1)
        abs(a[i]-b[i])/max(abs(b[i]),1e-12*unit)
    end
end

# DESIGN §4: relative error 1e-8 plus an absolute 1e-12 source-scale
# floor. Multiplying the floor by 1e-8 would instead demand 1e-20 absolute.
function source_budget_ratio(a,b,p)
    scales = (1.0,p.orbit.p,p.orbit.p^2,1.0,1.0,1.0,p.orbit.p)
    maximum(eachindex(SOURCE_KEYS)) do i
        unit=p.mass*scales[(i+1)÷2]/p.orbit.p/(iseven(i) ? p.orbit.p : 1)
        abs(a[i]-b[i])/(1e-8*abs(b[i])+1e-12*unit)
    end
end

if SOURCE_TEST_STAGE in ("all","azimuth")
@testset "ECC gr-azimuth-projector" begin
    @testset "R3 tensor, distance and radial jets" begin
        # Exact rational metric/velocity and dyadic zero-angle distance point.
        g = Matrix(Diagonal(Rational{Int}[-1,1,4,1]))
        u = Rational{Int}[3//2,1//2,1//2,0]
        K = ESW._eccentric_source_event(2//1,1//1,1//1,-1//1,0//1,4//1,1//1,u[1:3]...)
        tensor = 2*(g+2*(g*u)*transpose(g*u))
        @test K.tensor == Tuple(tensor[i,j] for (i,j) in SOURCE_PAIRS)
        @test eltype(K.tensor) == Rational{Int}
        exact = R3_SOURCE.distance_jet(g,u,2//1,3//1,0,1//1)
        runtime = ESW._eccentric_distance(K,3//1,0)
        @test exact.s2 == runtime.s2 == 9//4
        @test exact.s2_r == runtime.s2_r == 5//2
        for spin in (-0.7,0.0,0.9), e in (0.1,0.6), q in (0.7,2pi-0.7)
            d,p = source_fixture(spin,20.0,e)
            s,g,u = r3_source_inputs(p,q)
            r,psi = s.r+0.2,0.031
            h = ESW.singular_field(p,r,psi,q)
            f = R3_SOURCE.field_jet(g,u,s.r,r,psi,p.softening.eps,p.mass)
            for (c,(i,j)) in enumerate(SOURCE_PAIRS)
                @test h[2c-1] ≈ f.h[i,j] rtol=3e-13
                @test h[2c] ≈ f.h_r[i,j] rtol=3e-12
            end
            @test abs(h.htr)>0 && abs(h.hrphi)>0
            @test abs(h.htt-ESW.singular_field(p,r,-psi,q).htt)>1e-6
            step=1e-4*r
            fp=ESW.singular_field(p,r+step,psi,q)
            fm=ESW.singular_field(p,r-step,psi,q)
            fph=ESW.singular_field(p,r+step/2,psi,q)
            fmh=ESW.singular_field(p,r-step/2,psi,q)
            for c in 1:7
                fd=(4*(fph[2c-1]-fmh[2c-1])/step-(fp[2c-1]-fm[2c-1])/(2step))/3
                @test fd ≈ h[2c] rtol=2e-8 atol=1e-11
            end
        end
    end
    @testset "Strip adaptive FFT, quadrature and two doublings" begin
        println("PARITY,spin,e,block,delta,Nphi,quad_error,double_error,quadruple_error,budget_ratio")
        for spin in (-0.7,0.5), e in (0.3,0.6)
            d,p=source_fixture(spin,20.0,e)
            radii=[p.orbit.p/(1+e),20.0,38.0,90.0]
            b=ESW.source_azimuth_bank(d,p,radii;Nq=16,mmax=12,blocksize=2)
            b2=ESW.source_azimuth_bank(d,p,radii;Nq=16,mmax=12,blocksize=2,Nphi=2maximum(b.nphi))
            b4=ESW.source_azimuth_bank(d,p,radii;Nq=16,mmax=12,blocksize=2,Nphi=4maximum(b.nphi))
            for (ib,delta) in enumerate(b.strips)
                err,err2,err4,budget=0.0,0.0,0.0,0.0
                for ir in (2ib-1):(2ib), iq in (3,15), m in (0,1,3,12)
                    h=ESW.singular_field_mode(b,m,ir,iq)
                    oracle=r3_azimuth_oracle(p,m,radii[ir],b.q[iq])
                    err=max(err,source_error(h,oracle,p,radii[ir]))
                    err2=max(err2,source_error(h,ESW.singular_field_mode(b2,m,ir,iq),p,radii[ir]))
                    err4=max(err4,source_error(ESW.singular_field_mode(b2,m,ir,iq),ESW.singular_field_mode(b4,m,ir,iq),p,radii[ir]))
                    direct=ESW.singular_field_mode(p,m,radii[ir],b.q[iq];N=2maximum(b.nphi))
                    budget=max(budget,source_budget_ratio(h,oracle,p),source_budget_ratio(h,direct,p),
                        source_budget_ratio(h,ESW.singular_field_mode(b2,m,ir,iq),p),
                        source_budget_ratio(ESW.singular_field_mode(b2,m,ir,iq),ESW.singular_field_mode(b4,m,ir,iq),p))
                    @test source_budget_ratio(h,direct,p)<=1
                end
                @test budget<=1
                @test ispow2(b.nphi[ib])
                @printf("PARITY,%.1f,%.1f,%d,%.8g,%d,%.6g,%.6g,%.6g,%.6g\n",spin,e,ib,delta,b.nphi[ib],err,err2,err4,budget)
            end
            # Conservative strip has strictly positive real squared distance.
            for K in b.events, r in radii, x in range(0,2pi;length=19)
                delta=ESW.strip_halfwidth(K,r)
                @test real(ESW._eccentric_distance(K,r,x+im*delta).s2)>0
            end
            @test Base.summarysize(b) < b.workspace_bound_bytes
            @test b.workspace_bound_bytes <= 512*1024^2
        end
    end
    @testset "Circular strip and unchanged legacy bits" begin
        for spin in (-0.7,0.0,0.9), r in (12.0,20.0,30.0)
            d,p=source_fixture(spin,20.0,0.0)
            K=ESW._eccentric_source_event(p,worldline_rates(p.orbit,0.0))
            bg=ESW.DiscBackground(χ=spin,aEMRI=20.0,h0=0.06,b=0.3)
            @test 0<ESW.strip_halfwidth(K,r)<=ESW.strip_halfwidth(bg,r)
        end
        after=[ESW.singular_field_mode(SOURCE_LEGACY_BG,m,r;N=1024)
               for m in (0,1,7),r in (12.0,20.0,30.0)]
        @test isequal(SOURCE_LEGACY_BEFORE,after)
    end
    @testset "Memory cap, input rejection, precision inheritance" begin
        d,p=source_fixture()
        @test_throws ArgumentError ESW.source_azimuth_bank(d,p,[20.0];Nq=32,mmax=3,max_memory_bytes=1000)
        @test_throws ArgumentError ESW.source_azimuth_bank(d,p,[20.0];Nq=32,mmax=3,max_memory_bytes=30000)
        @test_throws ArgumentError ESW.source_azimuth_bank(d,p,[20.0];Nphi=1000)
        @test_throws ArgumentError ESW.source_azimuth_bank(d,p,Float64[])
        @test_throws ArgumentError ESW.source_azimuth_bank(d,p,[NaN])
        @test_throws ArgumentError ESW.source_azimuth_bank(d,p,[20.0];Nq=typemax(Int),mmax=typemax(Int))
        pu=ESW.EccentricPerturber(p.orbit;softening=InstantaneousScaleHeight(0.3))
        @test_throws ArgumentError ESW.source_azimuth_bank(d,pu,[20.0];Nq=2,mmax=1)
        setprecision(BigFloat,128) do
            db,pb=source_fixture(big"0.5",big"20",big"0.3";soft=big"0.36")
            b=ESW.source_azimuth_bank(db,pb,BigFloat[60];Nq=2,mmax=1)
            @test eltype(b.kernels)==Complex{BigFloat}
            @test precision(real(b.kernels[1]))==128
            @test precision(BigFloat)==128
        end
    end
end
println("PASS R3_EXACT rational tensor contraction and dyadic distance/derivative")
println("PASS R3_FLOAT64 seven components, both radial directions, three spins, derivatives")
println("PASS AZIMUTH_PARITY FFT vs independent R3 Gauss quadrature and Nphi,2Nphi,4Nphi")
println("PASS CIRCULAR_STRIP lower bound vs exact circular zeros; legacy bit identity")
println("PASS MEMORY cap/overflow/invalid input rejection; BigFloat128 tiny inheritance unit")
end

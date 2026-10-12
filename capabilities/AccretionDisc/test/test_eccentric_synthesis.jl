module EccentricSynthesisTests
using Test, AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves

struct ShiftedPhi <: AbstractVector{Float64}
    values::Vector{Float64}
end
Base.size(v::ShiftedPhi)=size(v.values)
Base.axes(v::ShiftedPhi)=(0:length(v.values)-1,)
Base.getindex(v::ShiftedPhi,i::Int)=v.values[i+1]
Base.iterate(v::ShiftedPhi,state=1)=iterate(v.values,state)
Base.IteratorSize(::Type{ShiftedPhi})=Base.HasLength()

@testset "Eccentric reality pairs and inertial time" begin
    r = [10.0,12.0,15.0,19.0]
    phi = [0.0,0.3,1.4,3.0]
    rq = [9.0,10.0,14.0,19.0,20.0]
    modes = [HarmonicMode(2,1,0.11),HarmonicMode(2,-2,-0.03),HarmonicMode(0,1,0.04)]
    amps = [1.0+2im,-0.4+0.7im,0.3-0.2im]
    sols = [(context=(;mode),r=r,y=repeat(reshape([a,2a],2,1),1,length(r)))
            for (mode,a) in zip(modes,amps)]
    for t in (0.0,11.0), component in (1,2)
        actual = SW.synthesize_eccentric(sols,rq,phi,t;component)
        expected = [10<=x<=19 ? sum(2component*real(a*cis(m.m*p-m.omega*t))
                    for (m,a) in zip(modes,amps)) : 0.0 for x in rq,p in phi]
        @test actual ≈ expected rtol=2e-14
    end
    @test_throws ArgumentError SW.synthesize_eccentric([],rq,phi,0.0)
    @test_throws ArgumentError SW.synthesize_eccentric(sols,rq,phi,0.0;component=3)
    bad = (context=(mode=HarmonicMode(0,0,0.0),),r=r,y=sols[1].y)
    @test_throws ArgumentError SW.synthesize_eccentric([bad],rq,phi,0.0)
    @test SW.synthesize_eccentric(sols,rq,ShiftedPhi(phi),11.0) ==
          SW.synthesize_eccentric(sols,rq,phi,11.0)
end

@testset "Circular co-rotating synthesis parity" begin
    orbit = EquatorialOrbit(0.4,20.0,0.0)
    r = [10.0,15.0,20.0,30.0]
    phi = collect(range(-pi,pi;length=19))
    rq = [9.0,10.0,16.0,20.0,30.0,31.0]
    circular = [(m=m,r=r,y=[complex(x*m,x/3) for j in 1:2,x in r]) for m in 1:3]
    eccentric = [(context=(mode=HarmonicMode(s.m,0,orbit_frequency(orbit,s.m,0)),),
                   r=s.r,y=s.y) for s in circular]
    for t in (0.0,27.0), component in (1,2)
        @test SW.synthesize_eccentric(eccentric,rq,phi,t;component) ≈
              SW.synthesize(circular,rq,phi .- orbit.Ωφ*t;component) rtol=3e-14
    end
end

include(joinpath(@__DIR__,"../../../reproduce/_emri_eccentric/first_results.jl"))
@testset "First-results sign, mass conversion and domain" begin
    FR = FirstEccentricResults
    for a in (-0.9,0.0,0.9), e in (0.2,0.5)
        d,p,domain = FR.fixture(a,e,Adiabatic)
        @test domain[1]>=SW.r_isco(d)
        @test domain[1]<20/(1+e)<20/(1-e)<domain[2]
        target=[-2e-12,-3e-14]
        J=Spacetime.EL_to_pe_jacobian(a,20.0,e)
        specific=J\target
        row=FR.rates(p.orbit,-p.mass*specific[2],-p.mass*specific[1])
        @test row["pdot"]≈target[1] rtol=1e-12
        @test row["edot"]≈target[2] rtol=1e-12
        @test row["Edot_specific"]≈specific[1]
        @test row["Ldot_specific"]≈specific[2]
    end
    @test FR.rates(EquatorialOrbit(0.0,20.0,0.0),1e-16,1e-18)["edot"]===nothing
    @test FR.modal_sum([Dict("T"=>v) for v in (1e16,1.0,-1e16)],"T")==1.0
    d,p,domain=FR.fixture(0.0,0.2,Adiabatic)
    c=SW.ModeContext(d,p,2,1;closure=Adiabatic,domain)
    fake=(perturber=p,grid=(edges=collect(domain),),modes=[c.mode])
    @test FR.validate_bank(fake,p,[c.mode],(domain[1]+.1,domain[2]-.1))===fake
    @test_throws ErrorException FR.validate_bank(fake,p,[c.mode],(domain[1]-.1,domain[2]))
    @test_throws ErrorException FR.validate_bank(fake,p,[HarmonicMode(3,1,.2)],domain)
    other=SW.EccentricPerturber(p.orbit;mass=2p.mass,softening=p.softening)
    @test_throws ErrorException FR.validate_bank(fake,other,[c.mode],domain)
    r=collect(range(domain...;length=61))
    g=SW.EccentricGrid(r,c.r_c,1:0,SW.mode_cache_key(c),(;))
    sol=SW.solve_mode(c;source=x->SW._zero_source(c),grid=g,
                      boundary_rows=SW.mode_boundary_rows(c))
    reduced=FR.summarize_mode(sol,r,r)
    @test iszero(reduced.record["torque_disc"])
    @test iszero(reduced.record["power_disc"])
    @test all(iszero,reduced.profiles)
    @test all(iszero,reduced.density)
end

# Lander addition (ECC::gr-first-results): the products writer is otherwise untested. A bound geodesic
# is reflection-symmetric about pericentre, so apocentre is exactly half a radial period after it, at
# azimuth pi*Omega_phi/Omega_r: independent of the trajectory's oscillating parts.
read_wake(path) = open(path) do io
    nr,nphi = read(io,Int64),read(io,Int64)
    r,phi = read!(io,zeros(nr)),read!(io,zeros(nphi))
    (;r,phi,field=read!(io,zeros(nr,nphi)))
end
@testset "First-results products: apocentre phase and conditional-sector exclusion" begin
    FR = FirstEccentricResults
    orbit = EquatorialOrbit(0.0,20.0,0.2)
    modes = [HarmonicMode(m,k,orbit_frequency(orbit,m,k)) for (m,k) in ((1,-1),(2,1),(3,-2))]
    torque = [5e-16,-2e-15,7e-16]
    records = [Dict{String,Any}("m"=>md.m,"k"=>md.k,"omega"=>md.omega,"torque_disc"=>t,
                   "power_disc"=>md.omega/md.m*t,"conditional_secular"=>(md.m,md.k)==(1,-1),
                   "nodes"=>10,"residual"=>0.0,"bc_residual"=>0.0,"power_identity_relative"=>0.0)
               for (md,t) in zip(modes,torque)]
    wr = collect(range(15.0,26.0;length=12)); phi = collect(range(-pi,pi;length=33))
    density = [complex(1+r/30,r/40) for r in wr]
    densities = [(context=(;mode=md),r=wr,y=vcat(transpose(density),zeros(ComplexF64,1,length(wr))))
                 for md in modes[2:3]]
    mktempdir() do dir
        out,bank = mkpath(joinpath(dir,"out")),mkpath(joinpath(dir,"bank"))
        meta = Dict{String,Any}("chi"=>0.0,"e"=>0.2)
        FR.write_products(out,bank,meta,records,wr,wr,phi,zeros(length(wr),3),densities,[3],[2])
        wake = meta["wake"]
        @test wake["t_peri"] == 0
        @test wake["t_apo"] ≈ pi/orbit.Ωr rtol=1e-11
        @test wake["phi_apo"] ≈ pi*orbit.Ωφ/orbit.Ωr rtol=1e-11
        @test (wake["r_peri"],wake["r_apo"]) == (20/1.2,20/0.8)
        @test orbit.trajectory(pi).r ≈ wake["r_apo"] rtol=1e-12
        row = only(meta["ladder"])
        @test row["modes"] == 2
        @test row["torque_disc"] == torque[2]+torque[3]
        @test row["Ldot_orbit"] == -row["torque_disc"] && row["Edot_orbit"] == -row["power_disc"]
        @test only(meta["secular"])["torque_disc"] == torque[1]
        peri,apo = read_wake(joinpath(bank,"wake_peri.bin")),read_wake(joinpath(bank,"wake_apo.bin"))
        @test (peri.r,peri.phi) == (wr,phi)
        @test peri.field == SW.synthesize_eccentric(densities,wr,phi,0.0)
        @test apo.field ≈ SW.synthesize_eccentric(densities,wr,phi,pi/orbit.Ωr) rtol=1e-11
        @test !isapprox(apo.field,peri.field;rtol=1e-3)
    end
end
end

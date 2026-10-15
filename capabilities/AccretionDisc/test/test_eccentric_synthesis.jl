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
        meta = Dict{String,Any}("chi"=>0.0,"e"=>0.2,"mu"=>1e-5)
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

# Lander addition (L-G52, lane w20-g52 stage 1). The element maps take the perturber mass from the caller
# (obligation gr-first-o-driver-mass), and the converged-results driver keeps its radial rule, provenance
# guards and native reduction fixed.
@testset "First-results element maps use the supplied perturber mass" begin
    FR = FirstEccentricResults
    orbit = EquatorialOrbit(0.9,20.0,0.2)
    a,b = FR.rates(orbit,2e-14,3e-16;mass=1e-5),FR.rates(orbit,2e-14,3e-16;mass=2e-5)
    @test b["pdot"] ≈ a["pdot"]/2 rtol=1e-14
    @test b["edot"] ≈ a["edot"]/2 rtol=1e-14
    @test b["Ldot_specific"] ≈ a["Ldot_specific"]/2 rtol=1e-14
    @test (b["Ldot_orbit"],b["Edot_orbit"]) == (a["Ldot_orbit"],a["Edot_orbit"]) == (-2e-14,-3e-16)
    for bad in (0.0,-1e-5,Inf,NaN)
        @test_throws ArgumentError FR.rates(orbit,1.0,1.0;mass=bad)
    end
    @test FR.rates(EquatorialOrbit(0.9,20.0,0.0),1e-16,1e-18;mass=2e-5)["Ldot_specific"] ≈ -5e-12 rtol=1e-14
    modes = [HarmonicMode(m,k,orbit_frequency(orbit,m,k)) for (m,k) in ((1,-1),(2,1))]
    records = [Dict{String,Any}("m"=>md.m,"k"=>md.k,"omega"=>md.omega,"torque_disc"=>t,
                   "power_disc"=>md.omega/md.m*t,"conditional_secular"=>(md.m,md.k)==(1,-1),
                   "nodes"=>10,"residual"=>0.0,"bc_residual"=>0.0,"power_identity_relative"=>0.0)
               for (md,t) in zip(modes,(5e-16,-2e-15))]
    wr = collect(range(15.0,26.0;length=8)); phi = collect(range(-pi,pi;length=9))
    dens = [(context=(;mode=modes[2]),r=wr,y=vcat(transpose(complex.(wr./30)),zeros(ComplexF64,1,length(wr))))]
    mktempdir() do dir
        out,bank = mkpath(joinpath(dir,"out")),mkpath(joinpath(dir,"bank"))
        meta = Dict{String,Any}("chi"=>0.9,"e"=>0.2,"mu"=>2e-5)
        FR.write_products(out,bank,meta,records,wr,wr,phi,zeros(length(wr),3),dens,[2],[1])
        for (row,rec) in ((only(meta["ladder"]),records[2]),(only(meta["secular"]),records[1]))
            want = FR.rates(orbit,rec["torque_disc"],rec["power_disc"];mass=2e-5)
            @test row["pdot"] == want["pdot"] && row["edot"] == want["edot"]
            @test row["pdot"] ≈ FR.rates(orbit,rec["torque_disc"],rec["power_disc"];mass=1e-5)["pdot"]/2 rtol=1e-14
        end
    end
end

include(joinpath(@__DIR__,"../../../reproduce/_emri_eccentric/converged_results.jl"))
@testset "Converged-results driver: radial rule, provenance and native reduction" begin
    CR = ConvergedEccentricResults
    d,p,domain = CR.FR.fixture(0.9,0.2,Adiabatic;epsilon_p=1e-6)
    @test d.εp == 1e-6 && CR.FR.fixture(0.9,0.2,Adiabatic)[1].εp == 1e-5
    ctx(m,k) = SW.ModeContext(d,p,m,k;closure=Adiabatic,domain)
    for ((m,k),n,class) in (((2,1),16384,"corotation"),((20,-21),65536,"corotation"),((40,-20),16384,"corotation"),
                            ((10,-15),32768,"negative_frequency_no_corotation"),
                            ((30,-40),65536,"negative_frequency_no_corotation"),
                            ((8,-9),32768,"positive_frequency_no_corotation"),
                            ((23,-26),65536,"positive_frequency_no_corotation"),
                            ((1,-1),16384,"secular"),((0,1),32768,"axisymmetric"))
        rule = CR.radial_rule(ctx(m,k))
        @test (rule.n,rule.cellclass,rule.certified) == (n,class,false)
    end
    c = ctx(2,1); g = CR.production_grid(c)
    @test length(g.r) == length(CR.FR.production_grid(c;n=16384).r) != length(CR.FR.production_grid(c).r)
    mktempdir() do root
        evidence = mkpath(joinpath(root,"converged","stage1")); banks = mkpath(joinpath(root,"banks"))
        @test_throws ErrorException CR.setup(evidence,banks,.9,.2,3,4,2)
        write(joinpath(evidence,"PREREG.json"),"{}"); write(joinpath(evidence,"STAGE1_AMENDMENT.json"),"{}")
        @test CR.setup(evidence,banks,.9,.2,3,4,2).domain == domain
        rm(joinpath(banks,"a+0.9_e0.2_ad");recursive=true)
        grand = mkpath(joinpath(root,"gr-first-results","a+0.9_e0.2_ad"))
        CR.FR.writejson(joinpath(grand,"result.json"),Dict("domain"=>[3.0,150.0]))
        x = CR.setup(evidence,banks,.9,.2,3,4,2)
        @test x.domain == (3.0,150.0) && x.config["domain"] == [3.0,150.0]
        @test x.config["epsilon_p"] == 1e-6 && x.d.εp == 1e-6
        @test sort([(s.m,s.k) for s in x.selected]) == [(m,k) for m in 1:3 for k in -4:2]
        @test all(s->s.m==0,x.axisymmetric_modes) && !isempty(x.axisymmetric_modes)
        @test (1,-1) in x.sec
        @test CR.setup(evidence,banks,.9,.2,3,4,2).config == x.config
        @test_throws ErrorException CR.setup(evidence,banks,.9,.2,3,5,2)
        rm(joinpath(banks,"a+0.9_e0.2_ad");recursive=true)
        parent = mkpath(joinpath(root,"converged","gr-first-results","a+0.9_e0.2_ad"))
        CR.FR.writejson(joinpath(parent,"result.json"),Dict("domain"=>[4.0,120.0]))
        @test CR.setup(evidence,banks,.9,.2,3,4,2).domain == (4.0,120.0)
        write(joinpath(evidence,"STAGE1_AMENDMENT.json"),"{\"amended\":1}")
        @test_throws ErrorException CR.setup(evidence,banks,.9,.2,3,4,2)
        stale = mkpath(joinpath(root,"stale","a+0.9_e0.2_ad")); write(joinpath(stale,"mode_1_0.jls"),"unversioned")
        @test_throws ErrorException CR.setup(evidence,joinpath(root,"stale"),.9,.2,3,4,2)
    end
    c = ctx(2,-3)
    @test c.r_c === nothing
    r = collect(range(domain...;length=1025))
    g = SW.EccentricGrid(r,c.r_c,1:0,SW.mode_cache_key(c),(;))
    bump(x) = exp(-(x-20)^2/8)
    src = x -> merge(SW._zero_source(c),(htt=complex(bump(x)),htt_r=complex(-(x-20)/4*bump(x)),hphiphi=complex(x^2*bump(x)/4)))
    sol = SW.solve_mode(c;source=src,grid=g,boundary_rows=SW.mode_boundary_rows(c))
    ref,new = CR.FR.summarize_mode(sol,r,r),CR.reduce_mode(sol,r,r)
    @test !iszero(ref.record["torque_disc"]) && !iszero(ref.record["power_disc"])
    @test new.record["torque_disc"] ≈ ref.record["torque_disc"] rtol=1e-12
    @test new.record["power_disc"] ≈ ref.record["power_disc"] rtol=1e-12
    @test new.record["power_disc"] ≈ c.mode.omega/c.mode.m*new.record["torque_disc"] rtol=1e-9
    @test new.profiles ≈ ref.profiles rtol=1e-9
    @test new.density == ref.density
    @test new.record["residual"] == sol.residual && new.record["nodes"] == length(r)
end

@testset "Ridge source bank construction, resume and provenance" begin
    CR=ConvergedEccentricResults
    ridge=Dict("window"=>Dict("negative_floor"=>40,"positive_floor"=>20,
                              "padding"=>8,"k_apo"=>-.354,"k_peri"=>.562),
               "named_exclusions"=>[[1,22],[1,24]])
    @test CR.mode_window(60,40,42,ridge)==(-40,42)
    @test CR.mode_window(80,40,53,ridge)==(-40,53)
    @test CR.mode_window(100,44,65,ridge)==(-44,65)
    @test CR.mode_window(60,40,20,nothing)==(-40,20)
    mktempdir() do root
        evidence=mkpath(joinpath(root,"evidence"));banks=joinpath(root,"banks")
        write(joinpath(evidence,"PREREG.json"),"{}")
        write(joinpath(evidence,"STAGE1_AMENDMENT.json"),"{}")
        CR.FR.writejson(joinpath(evidence,"PREREG_ridge_v2.json"),ridge)
        @test_throws ErrorException CR.setup(evidence,banks,.9,.2,60,40,20)
        x=CR.setup(evidence,banks,.9,.2,80,40,53)
        @test length(x.selected)==5892
        @test count(s->s.m<=60,x.selected)==4110
        @test count(s->s.m<=60 && s.k>20,x.selected)==450
        @test all(s->s.k<=20,filter(s->s.m==1,x.selected))
        @test (1,-1) in Set((s.m,s.k) for s in x.selected)
        # A small physical bank exercises the same builder at inexpensive quadrature.
        x=merge(x,(selected=filter(s->s.m==2 && s.k in (0,1),x.selected),))
        x.config["Nq"]=16;x.config["Nphi"]=64
        CR.build_bank(x)
        file=CR.source_block(x.dir,2)
        bank=CR.read_source_block(x,file,x.selected)
        direct=SW.source_radial_bank(x.d,x.p,x.domain,x.selected;Nq=16,Nphi=64,method=:fft,max_memory_bytes=8*1024^3)
        @test bank.jets==direct.jets
        @test bank.grid.rs==direct.grid.rs
        firsthash=bytes2hex(open(CR.sha256,file))
        CR.build_bank(x)
        @test bytes2hex(open(CR.sha256,file))==firsthash
        # Reuse a strict subset and compute the missing column; compare the full bank.
        otherdir=mkpath(joinpath(root,"other"));otherout=mkpath(joinpath(root,"otherout"))
        reuse=Dict("manifest"=>x.bankfile,"manifest_sha256"=>bytes2hex(open(CR.sha256,x.bankfile)))
        y=merge(x,(dir=otherdir,out=otherout,bankfile=joinpath(otherdir,"source_manifest.json"),
            selected=[x.selected;HarmonicMode(2,2,orbit_frequency(x.p.orbit,2,2))],
            ridge=merge(ridge,Dict("source_reuse"=>reuse))))
        CR.build_bank(y)
        combined=CR.read_source_block(y,CR.source_block(y.dir,2),y.selected)
        @test combined.jets[:,1:2]==bank.jets
        third=SW.source_radial_bank(x.d,x.p,x.domain,y.selected[3:3];Nq=16,Nphi=64,method=:fft,max_memory_bytes=8*1024^3)
        @test combined.jets[:,3:3]==third.jets
        # Positive extent can exceed the negative extent in the ridge window.
        # The final products row must include the entire requested positive band.
        for (i,s) in enumerate(x.selected)
            record=Dict{String,Any}("m"=>s.m,"k"=>s.k,"omega"=>s.omega,
                "torque_disc"=>i*1e-16,"power_disc"=>s.omega/s.m*i*1e-16,
                "conditional_secular"=>false,"nodes"=>4,"residual"=>0.,
                "bc_residual"=>0.,"power_identity_relative"=>0.)
            CR.serialize(joinpath(x.dir,"mode_$(s.m)_$(s.k).jls"),
                (;record,profiles=zeros(1600,3),density=zeros(ComplexF64,768)))
        end
        CR.collect_case(x,2,0,1)
        lastrow=CR.FR.json_read(read(joinpath(x.out,"result.json"),String))["ladder"][end]
        @test lastrow["k_max"]==1 && lastrow["modes"]==2
        @test lastrow["torque_disc"]≈3e-16 rtol=1e-15
        open(file,"a") do io;write(io,"changed");end
        @test_throws ErrorException CR.build_bank(x)
    end
end

# Lander L-G58 additions (lane w22-g58): each assertion below kills a mutant that survived the lane's tests.
@testset "Ridge registration, named exclusions and bank-reuse guards" begin
    CR=ConvergedEccentricResults
    ridge=Dict("window"=>Dict("negative_floor"=>40,"positive_floor"=>24,"padding"=>8,"k_apo"=>-.354,"k_peri"=>.562),
               "named_exclusions"=>[[1,22],[1,24]])
    mktempdir() do root
        evidence=mkpath(joinpath(root,"evidence"));banks=joinpath(root,"banks")
        write(joinpath(evidence,"PREREG.json"),"{}");write(joinpath(evidence,"STAGE1_AMENDMENT.json"),"{}")
        regfile=joinpath(evidence,"PREREG_ridge_v2.json");CR.FR.writejson(regfile,ridge)
        x=CR.setup(evidence,banks,.9,.2,2,40,24)
        pairs=Set((s.m,s.k) for s in x.selected)
        # named exclusions leave the selection; their neighbours stay
        @test !((1,22) in pairs) && !((1,24) in pairs) && (1,23) in pairs && (2,22) in pairs && (2,24) in pairs
        @test length(x.selected)==2*65-2
        # the registration file is part of the checkpoint provenance
        @test x.config["ridge_prereg_sha256"]==bytes2hex(open(CR.sha256,regfile))
        ridge["window"]["padding"]=9;CR.FR.writejson(regfile,ridge)
        @test_throws ErrorException CR.setup(evidence,banks,.9,.2,2,40,24)
        ridge["window"]["padding"]=8;CR.FR.writejson(regfile,ridge)
        x=CR.setup(evidence,banks,.9,.2,2,40,24)
        x=merge(x,(selected=filter(s->s.m==2 && s.k in (0,1),x.selected),))
        x.config["Nq"]=16;x.config["Nphi"]=64
        CR.build_bank(x)
        # products carry the registration and the UNCONVERGED status
        for (i,s) in enumerate(x.selected)
            record=Dict{String,Any}("m"=>s.m,"k"=>s.k,"omega"=>s.omega,"torque_disc"=>i*1e-16,"power_disc"=>s.omega/s.m*i*1e-16,
                "conditional_secular"=>false,"nodes"=>4,"residual"=>0.,"bc_residual"=>0.,"power_identity_relative"=>0.)
            CR.serialize(joinpath(x.dir,"mode_$(s.m)_$(s.k).jls"),(;record,profiles=zeros(1600,3),density=zeros(ComplexF64,768)))
        end
        CR.collect_case(x,2,0,1)
        meta=CR.FR.json_read(read(joinpath(x.out,"result.json"),String))
        @test meta["ridge_registration"]["named_exclusions"]==[[1,22],[1,24]] && startswith(meta["status"],"UNCONVERGED")
        extra=HarmonicMode(2,2,orbit_frequency(x.p.orbit,2,2))
        reuse=Dict("manifest"=>x.bankfile,"manifest_sha256"=>bytes2hex(open(CR.sha256,x.bankfile)))
        reusing(name,r;domain=x.domain)=begin
            dir=mkpath(joinpath(root,name));out=mkpath(joinpath(root,name*"_out"))
            merge(x,(dir=dir,out=out,bankfile=joinpath(dir,"source_manifest.json"),selected=[x.selected;extra],domain,
                     ridge=merge(ridge,Dict("source_reuse"=>r))))
        end
        # a reuse manifest whose recorded hash differs is refused before any block is read
        @test_throws ErrorException CR.build_bank(reusing("badmanifest",merge(reuse,Dict("manifest_sha256"=>"0"^64))))
        # a reused block that no longer matches its manifest entry is refused
        copydir=joinpath(root,"copy");cp(x.dir,copydir)
        open(io->write(io,"changed"),CR.source_block(copydir,2),"a")
        copied=Dict("manifest"=>joinpath(copydir,"source_manifest.json"),
                    "manifest_sha256"=>bytes2hex(open(CR.sha256,joinpath(copydir,"source_manifest.json"))))
        @test_throws ErrorException CR.build_bank(reusing("badblock",copied))
        # columns reused from a bank on another radial grid are not merged with fresh ones
        wide=merge(x,(dir=mkpath(joinpath(root,"wide")),out=mkpath(joinpath(root,"wide_out")),
                      bankfile=joinpath(root,"wide","source_manifest.json"),domain=(x.domain[1],x.domain[2]+10)))
        CR.build_bank(wide)
        widereuse=Dict("manifest"=>wide.bankfile,"manifest_sha256"=>bytes2hex(open(CR.sha256,wide.bankfile)))
        @test_throws ErrorException CR.build_bank(reusing("grid",widereuse))
    end
end

end

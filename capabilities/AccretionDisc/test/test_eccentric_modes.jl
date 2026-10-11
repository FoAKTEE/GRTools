module EccentricModeTestSupport
using Test, LinearAlgebra, Printf
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
for (name,file) in ((:SourceAzimuthBank,"eccentric_source.jl"),
                    (:radiation_impedance_row,"boundary.jl"),
                    (:EccentricModeSource,"eccentric_modes.jl"))
    isdefined(SW,name) || Base.include(SW,joinpath(@__DIR__,"../src/SpiralWaves",file))
end

function fixture(;closure=Adiabatic,q=0.75,a=0.4,p=20.0,e=0.2,m=2,k=1,
                 domain=(28.0,32.0),soft=2.0,εp=1e-3)
    d=SW.DiscModel(;χ=a,r_ref=p,h_ref=0.08,pSigma=1.5,qTemp=q,closure,evolution=:fixed,εp)
    pert=SW.EccentricPerturber(EquatorialOrbit(a,p,e);mass=1e-5,softening=FixedSoftening(soft))
    c=SW.ModeContext(d,pert,m,k;closure,domain)
    c,pert
end

function uniform_grid(c,n)
    r=collect(range(c.domain...;length=n+1))
    SW.EccentricGrid(r,c.r_c,1:0,SW.mode_cache_key(c),(;manufactured=true))
end

function exact_state(c,r)
    lo,hi=c.domain; w=hi-lo; x=(r-lo)/w
    f=x*(1-x); df=(1-2x)/w
    z=exp((1+2im)*x)
    y=[f*z,(0.2-0.1im)*f*(1+x)]
    yp=[(df+(1+2im)*f/w)*z,(0.2-0.1im)*(df*(1+x)+f/w)]
    y,yp
end

# Independent continuous forcing y'-C*y on the actual emitted operator.
# The two free source jets below manufacture that forcing through the context
# assembly. They are algebraic jets, not a claimed physical metric solution.
function manufactured_source(c,r)
    h=merge(SW._zero_source(c),(htr=(1+2im)/r^2,hrphi=(2-im)/r))
    C,S=SW.master_matrices(c,r;hmodes=h,landau=false)
    z=SW._zero_source(c)
    a=last(SW.master_matrices(c,r;hmodes=merge(z,(htt_r=1+0im,)),landau=false))
    b=last(SW.master_matrices(c,r;hmodes=merge(z,(hzz=1+0im,)),landau=false))
    y,yp=exact_state(c,r)
    v=hcat(a,b)\(yp-C*y-S)
    merge(h,(htt_r=v[1],hzz=v[2]))
end

# True continuous defect from centered derivatives at nodes. This differs from
# the algebraic midpoint residual, which should already be at roundoff.
function continuous_defect(sol)
    c=sol.context; r=sol.r; y=sol.y
    maximum(2:length(r)-1) do i
        hm,hp=r[i]-r[i-1],r[i+1]-r[i]
        yp=-hp/(hm*(hm+hp))*y[:,i-1]+(hp-hm)/(hm*hp)*y[:,i]+hm/(hp*(hm+hp))*y[:,i+1]
        C,S=SW.master_matrices(c,r[i];hmodes=sol.source(r[i]),landau=sol.params.landau)
        norm(yp-C*y[:,i]-S,Inf)
    end
end

function manufactured_checks(c)
    errors=Float64[]; defects=Float64[]
    for n in (32,64,128,256)
        g=uniform_grid(c,n)
        sol=SW.solve_mode(c;grid=g,source=r->manufactured_source(c,r),landau=false,
                          boundary_rows=SW.mode_boundary_rows(c))
        err=maximum(norm(sol.y[:,i]-first(exact_state(c,r)),Inf) for (i,r) in enumerate(g.r))
        push!(errors,err); push!(defects,continuous_defect(sol))
        @test sol.residual<2e-10 && sol.bc_residual<1e-11
        @printf("MANUFACTURED,%s,%.2f,%d,%.12e,%.12e,%.12e\n",string(c.closure),c.disc.qTemp,n,err,last(defects),sol.residual)
    end
    @test all(3.6 .< errors[1:end-1]./errors[2:end] .< 4.4)
    @test all(3.3 .< defects[1:end-1]./defects[2:end] .< 4.7)
    println("PASS manufactured: ",c.closure," q=",c.disc.qTemp," second-order error and continuous defect; banded residual")
end

# Linux process CPU clock: excludes cgroup scheduling delays. No multi-thread
# work is launched; the evidence runner also measures total child CPU time.
cpu_seconds()=Float64(ccall(:clock,Clong,()))/1e6

# Lander review (mutants MA4/MA5/MA6/ML2): the batched table must equal the
# source owner's projection from an independently built one-radius bank, and
# reconstruct_fluid(sol) must use that source and the solve's Landau setting.
function sampled_source_checks(sol,p,Nq;indices=(length(sol.r)÷2,))
    c=sol.context
    for i in indices
        for r in (sol.r[i],(sol.r[i]+sol.r[i+1])/2)
            h=SW.source_modes(SW.source_azimuth_bank(c.disc,p,[r];Nq,mmax=c.mode.m),(c.mode,))[1,1]
            @test all(isapprox(sol.source(r)[n],h[n];rtol=1e-12,atol=1e-14*abs(h.htt)) for n in keys(h))
        end
        r=sol.r[i]; h=sol.source(r)
        f=SW.reconstruct_fluid(sol)[i]
        ref=SW.reconstruct_fluid(c,r,sol.y[:,i];hmodes=h,landau=sol.params.landau)
        @test all(isapprox(f[n],ref[n];rtol=1e-12) for n in keys(ref))
        # Both checks are sensitive at this node: the source and the Landau
        # shift each change some reconstructed amplitude.
        change(g)=maximum(n->abs(f[n]-g[n]),keys(f))/maximum(n->abs(f[n]),keys(f))
        @test change(SW.reconstruct_fluid(c,r,sol.y[:,i];hmodes=SW._zero_source(c),landau=sol.params.landau))>1e-6
        if c.r_c!==nothing && sol.params.landau
            @test change(SW.reconstruct_fluid(c,r,sol.y[:,i];hmodes=h,landau=false))>1e-8
        end
    end
end

nearest_corotation_node(sol)=sol.context.r_c===nothing ? nothing :
    argmin(i->sol.r[i]==sol.context.r_c ? Inf : abs(sol.r[i]-sol.context.r_c),1:length(sol.r)-1)

# Dispatch isolation redefines kernels to throw, so it runs in a fresh process.
function isolation_run(kind)
    script=joinpath(@__DIR__,"eccentric_modes_isolation.jl")
    cmd=`$(Base.julia_cmd()) -t1 --startup-file=no --project=$(Base.active_project()) $script $kind`
    out=IOBuffer()
    proc=run(pipeline(ignorestatus(cmd);stdout=out,stderr=out))
    print(String(take!(out)))
    proc.exitcode==0
end

function forced_checks(c,p;ns=(32,64,128),Nq=32)
    sols=Any[]; defects=Float64[]
    for n in ns
        g=uniform_grid(c,n)
        t=cpu_seconds(); wall=time()
        s=SW.mode_source(c,p,g;bank_options=(;Nq),radial_batch_size=32)
        sourcecpu=cpu_seconds()-t; sourcewall=time()-wall
        rows=SW.mode_boundary_rows(c)
        # Warm the concrete solver path, then report an uncached solve including
        # emitted operators, source assembly, boundary selection and banded LU.
        SW.solve_mode(c;grid=g,source=s,boundary_rows=rows)
        t=cpu_seconds(); wall=time()
        sol=SW.solve_mode(c;grid=g,source=s,boundary_rows=SW.mode_boundary_rows(c))
        elapsed=time()-wall; cpu=cpu_seconds()-t
        push!(sols,sol); push!(defects,continuous_defect(sol))
        @test sol.residual<2e-9 && sol.bc_residual<1e-11
        @test sol.params.linsolve==:banded && sol.closure==c.closure
        fluid=SW.reconstruct_fluid(sol)
        @test all(f->all(isfinite,values(f)),fluid)
        @test maximum(abs,sol.y)>1e-12
        @printf("FORCED,%s,%.2f,%d,%.12e,%.12e,%.12e\n",string(c.closure),c.disc.qTemp,n,last(defects),sol.residual,sol.bc_residual)
        @printf("COST,%s,%.2f,%.2f,%.1f,%.2f,%d,%d,%d,%d,%.9f,%.9f,%.9f,%.9f\n",
                string(c.closure),c.disc.qTemp,p.orbit.a,p.orbit.p,p.orbit.e,c.mode.m,c.mode.k,n+1,Nq,cpu,elapsed,sourcecpu,sourcewall)
    end
    changes=[norm(sols[i].y-sols[i+1].y[:,1:2:end])/norm(sols[i+1].y[:,1:2:end]) for i in 1:length(sols)-1]
    @test 3.3 < changes[1]/changes[2] < 4.7
    @test all(3.2 .< defects[1:end-1]./defects[2:end] .< 4.8)
    @printf("REFINEMENT,%s,%.2f,%.12e,%.12e,%.8f\n",string(c.closure),c.disc.qTemp,changes...,log2(changes[1]/changes[2]))
    sampled_source_checks(sols[end],p,Nq)
    # Independent temporal projection refinement; the fixed Nq BVP refinement
    # above must not be dominated by source quadrature error.
    g=uniform_grid(c,last(ns))
    fine=SW.solve_mode(c,p;grid=g,bank_options=(;Nq=2Nq))
    source_error=norm(fine.y-sols[end].y)/norm(fine.y)
    @test source_error<1e-7
    @printf("SOURCE_REFINEMENT,%s,%.2f,%d,%d,%.12e\n",string(c.closure),c.disc.qTemp,Nq,2Nq,source_error)
    println("PASS forced: ",c.closure," q=",c.disc.qTemp," source bank/context/BC/BVP/reconstruction and refinement")
    sols[end]
end

function precision_checks(closure)
    setprecision(BigFloat,160) do
        d=SW.DiscModel(;χ=big"0.4",r_ref=big"20",h_ref=big"0.08",pSigma=big"1.5",
                       qTemp=big"0.75",closure,evolution=:fixed)
        c=SW.ModeContext(d,HarmonicMode(2,1,big"0.04");closure,domain=(big"28",big"32"))
        rows=SW.mode_boundary_rows(c)
        @test all(x->x isa Complex{BigFloat} && precision(real(x))==160,first(rows)[1])
        # Two-cell manufactured unit, never a BigFloat production source run.
        sol=SW.solve_mode(c;grid=uniform_grid(c,2),source=r->manufactured_source(c,r),
                          landau=false,boundary_rows=rows)
        @test eltype(sol.y)==Complex{BigFloat} && precision(real(sol.y[1]))==160
        @test SW.reconstruct_fluid(sol)[2].deps isa Complex{BigFloat}
    end
    println("PASS precision: ",closure," inherited BigFloat160 two-cell unit only")
end

function adiabatic_checks()
    @testset "ECC adiabatic manufactured and forced modes" begin
        c,p=fixture()
        @test c.r_c===nothing
        manufactured_checks(c)
        sol=forced_checks(c,p)
        h=sol.source(30.0)
        @test abs(h.htr)>1e-12 && abs(h.hrphi)>1e-12
        @test abs(imag(h.htr))>1e-12 && abs(imag(h.hrphi))>1e-12
        inp=SW.master_inputs(c,30.0;hmodes=h)
        @test inp.htr===h.htr && inp.hrphi===h.hrphi
        @test_throws ArgumentError sol.source(30.001)
        @test_throws ArgumentError SW.mode_source(c,p,uniform_grid(c,8);radial_batch_size=0)
        wrong,_=fixture(k=2)
        @test_throws ArgumentError SW.mode_source(wrong,p,uniform_grid(c,8))
        @test_throws ArgumentError SW.mode_source(c,last(fixture(e=0.3)),uniform_grid(c,8))
        bank=SW.source_azimuth_bank(c.disc,p,[28.0,30.0,32.0];Nq=16,mmax=2)
        tab=SW.mode_source(bank,c.mode)
        @test tab(30.0)==SW.source_modes(bank,(c.mode,))[2,1]
        precision_checks(Adiabatic)
    end
    @testset "ECC physical mode-specific grids and C-66 boundaries" begin
        for k in (-6,0,1)
            c,p=fixture(k=k,domain=(10.0,42.0))
            g=SW.build_grid(c;nL=1,nR=1,nC=65,n_anc=64)
            @test g.r_c==c.r_c
            @test (c.r_c===nothing)==isempty(g.central)
            c.r_c===nothing || @test c.r_c in g.r
            sol=SW.solve_mode(c,p;grid=g,bank_options=(;Nq=64))
            @test sol.residual<1e-8 && sol.bc_residual<1e-10
            sampled_source_checks(sol,p,64;indices=filter(!isnothing,(length(sol.r)÷2,nearest_corotation_node(sol))))
            for (r,outward,row) in ((first(c.domain),false,sol.params.rows[1][1]),
                                    (last(c.domain),true,sol.params.rows[2][1]))
                C=SW.master_operator(c,r;landau=false)
                pair=SW.eigen2(C)
                v=(-row[2],row[1]); λ=dot(collect(v),C*collect(v))/dot(collect(v),collect(v))
                δω=1e-6*abs(c.mode.omega)
                roots=map((-δω,δω)) do delta
                    cx=SW.ModeContext(c.disc,HarmonicMode(c.mode.m,k,c.mode.omega+delta);
                                      closure=Adiabatic,domain=c.domain)
                    values=SW.eigen2(SW.master_operator(cx,r;landau=false))[1:2]
                    values[argmin(abs.(values.-λ))]
                end
                slope=imag(roots[2]-roots[1])/(2δω)
                @test (outward ? 1 : -1)*slope>0
            end
            @printf("GRID,%d,%.12e,%s,%d,%.12e\n",k,c.mode.omega,string(c.r_c),length(g.r),sol.residual)
        end
    end
    println("PASS physical: C-66 arbitrary total frequency, signed omega, corotation present/absent, adaptive build_grid")
    @testset "ECC circular gate through eccentric source and BVP" begin
        for a in (-0.4,0.4), m in (2,3)
            bg=SW.DiscBackground(χ=a,aEMRI=20.0,h0=0.08,𝔭=1.5,𝔮=0.75,b=0.3,εp=1e-3)
            d=SW.DiscModel(bg;closure=Adiabatic,evolution=:fixed)
            p=SW.EccentricPerturber(EquatorialOrbit(a,20.0,0.0);mass=1e-5,
                                    softening=FixedSoftening(SW.softening_length(bg)))
            cg=SW.build_grid(bg,m;nL=1,nR=1,nC=65,n_anc=64,extent=:fiducial,r_in=8.0,r_out=45.0)
            c=SW.ModeContext(d,p,m,0;closure=Adiabatic,domain=(first(cg.r),last(cg.r)))
            eg=SW.EccentricGrid(copy(cg.r),c.r_c,1:0,SW.mode_cache_key(c),(;circular_gate=true))
            circ=SW.solve_mode(bg,m;grid=cg,mEMRI=p.mass,linsolve=:banded)
            ecc=SW.solve_mode(c,p;grid=eg,bank_options=(;Nq=2))
            err=norm(ecc.y-circ.y)/norm(circ.y)
            @test err<2e-8
            @printf("CIRCULAR,%.2f,%d,%d,%.12e,%.12e\n",a,m,length(eg.r),err,ecc.residual)
            # Independent reconstruction route (lander review, MA4/MA5): the
            # circular route-R2 covariant delta u_phi at every node.
            du=SW.delta_u_phi(bg,circ)
            duerr=norm([f.du_phi for f in SW.reconstruct_fluid(ecc)]-du)/norm(du)
            @test duerr<2e-8
            @printf("RECONSTRUCTION,%.2f,%d,%.12e\n",a,m,duerr)
            if a==0.4 && m==2
                errors=Float64[]
                for e in (1e-3,1e-4)
                    pe=SW.EccentricPerturber(EquatorialOrbit(a,20.0,e);mass=p.mass,softening=p.softening)
                    ce=SW.ModeContext(d,pe,m,0;closure=Adiabatic,domain=c.domain)
                    ge=SW.EccentricGrid(copy(cg.r),ce.r_c,1:0,SW.mode_cache_key(ce),(;circular_gate=true))
                    se=SW.solve_mode(ce,pe;grid=ge,bank_options=(;Nq=16))
                    push!(errors,norm(se.y-circ.y)/norm(circ.y))
                end
                @test 80<errors[1]/errors[2]<120
                @printf("E_TO_ZERO,%.12e,%.12e,%.8f\n",errors...,errors[1]/errors[2])
            end
        end
    end
    println("PASS circular: e=0 eccentric bank/context/C-66/BVP vs circular G26 solver; e to zero quadratic approach")
    # Lander review (mutant MA2): a circular-kernel dispatch at e=0 reproduces
    # the gate above exactly, so the same solves are repeated in a fresh process
    # with every circular source/master/selector method redefined to throw.
    @testset "ECC e=0 circular-kernel isolation (fresh process)" begin
        @test isolation_run("circular")
    end
end
end # module EccentricModeTestSupport

if !isdefined(@__MODULE__,:ECC_MODES_HELPERS_ONLY)
    EccentricModeTestSupport.adiabatic_checks()
end

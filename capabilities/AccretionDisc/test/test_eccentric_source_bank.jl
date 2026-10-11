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

function source_measure(f)
    start=ccall(:clock,Clong,())
    t=@timed f()
    cpu=(ccall(:clock,Clong,())-start)/1e6 # POSIX CLOCKS_PER_SEC on the test host
    (;value=t.value,wall=t.time,cpu,allocated=t.bytes)
end

# FR22 physical-time double integral, with E-ECC-2's n/(2pi^2).
# Absolute phi is integrated, not a pre-factorised Laplace coefficient.
function fr22_direct_source(e,radii,modes;Nt=256,Nphi=512,n=0.37)
    result=zeros(ComplexF64,length(radii),length(modes),2)
    dt=2pi/(n*Nt)
    dphi=2pi/Nphi
    weight=n/(2pi^2)*dt*dphi
    for j in 0:Nt-1
        M=2pi*j/Nt
        s=AccretionDisc.Modes.planet_orbit(M,e;a=1.0,n)
        for l in 0:Nphi-1
            phi=2pi*l/Nphi
            c=cos(phi-s.phi_p)
            for (ir,r) in enumerate(radii)
                sep=sqrt((r-s.r_p)^2+2r*s.r_p*(1-c)+0.018^2)
                potential=-1/sep
                derivative=(r-s.r_p*c)/sep^3
                for (imode,(m,k)) in enumerate(modes)
                    phase=cis((m+k)*M-m*phi)
                    result[ir,imode,1]+=weight*potential*phase
                    result[ir,imode,2]+=weight*derivative*phase
                end
            end
        end
    end
    result
end

if SOURCE_TEST_STAGE in ("all","harmonics")
include(joinpath(@__DIR__,"..","..","Spacetime","test","equatorial_quadrature_oracle.jl"))
using .EquatorialQuadratureOracle
@testset "ECC gr-source-harmonics" begin
    @testset "C-HARM owner and independent coordinate-time gate" begin
        for spin in (-0.7,0.5), e in (0.3,0.6)
            d,p=source_fixture(spin,20.0,e)
            rs=[9.0,70.0]
            b=ESW.source_azimuth_bank(d,p,rs;Nq=256,mmax=3,blocksize=2)
            pairs=((0,2),(1,-1),(2,-3),(3,2))
            modes=[HarmonicMode(m,k,orbit_frequency(p.orbit,m,k)) for (m,k) in pairs]
            h=ESW.source_modes(b,modes)
            maxbudget=0.0
            for (j,mode) in enumerate(modes), (ir,r) in enumerate(rs)
                H(q,s)=collect(values(ESW.singular_field_mode(p,mode.m,r,q;N=2048)))
                reference=orbit_harmonic(H,p.orbit,mode.m,mode.k;N=128)
                oracle=NamedTuple{SOURCE_KEYS}(Tuple(reference))
                maxbudget=max(maxbudget,source_budget_ratio(h[ir,j],oracle,p))
                @test source_budget_ratio(h[ir,j],oracle,p)<=1
            end
            @printf("CHARM,%.1f,%.1f,%.9g\n",spin,e,maxbudget)
            @test abs(sum(b.time_weights)-1)<1e-13
            # Demonstrate sensitivity to dropping dt/dlambda from the measure.
            mode=modes[4]
            unweighted=sum(b.events[i].tensor[1]*b.kernels[mode.m+1,i,1,1]*
                cis(mode.k*b.q[i]+mode.omega*b.delta_t[i]-mode.m*b.delta_phi[i])
                for i in eachindex(b.q))/length(b.q)
            @test abs(unweighted-h[1,4].htt)>1e-5abs(h[1,4].htt)
            @test abs(imag(h[1,4].htr))>1e-8
            nod=ESW.source_modes(b,modes;derivatives=false)
            @test all(nod[i].htt==h[i].htt && nod[i].hrphi==h[i].hrphi for i in eachindex(h))
            @test size(ESW.source_modes(b,HarmonicMode{Float64}[]))==(2,0)
            @test_throws ArgumentError ESW.source_modes(b,[HarmonicMode(4,0,4p.orbit.Ωφ)])
            @test_throws ArgumentError ESW.source_modes(b,[HarmonicMode(1,0,1.0)])
        end
        # Consume the Spacetime owner's Darwin-angle quadrature oracle, whose
        # time/azimuth integrals do not use trajectory deviations or Mino FFTs.
        d,p=source_fixture(-0.7,20.0,0.3)
        qo=quadrature_orbit(p.orbit.a,p.orbit.p,p.orbit.e;n=64)
        m,k,r=2,-1,9.0
        omega=(m*qo.Υφ+k*qo.Υr)/qo.Γ
        Tr=2pi*qo.Γ/qo.Υr
        function coordinate_integrand(chi)
            rp=qo.radius(chi)
            rates=qo.rates(chi)
            dl=qo.dlambda(chi)
            t=integrate_oracle(z->qo.rates(z).t*qo.dlambda(z),0.0,chi,qo.nodes,qo.weights)
            phi=integrate_oracle(z->qo.rates(z).phi*qo.dlambda(z),0.0,chi,qo.nodes,qo.weights)
            Ur=p.orbit.p*p.orbit.e*sin(chi)/(1+p.orbit.e*cos(chi))^2/dl
            state=(r=rp,ut=rates.t/rp^2,ur=Ur/rp^2,uphi=rates.phi/rp^2)
            h=r3_azimuth_oracle(p,m,r,state;panels=8,order=24)
            collect(values(h))*cis(omega*t-m*phi)*rates.t*dl/Tr
        end
        direct=integrate_oracle(coordinate_integrand,0.0,pi,qo.nodes,qo.weights)+
               integrate_oracle(coordinate_integrand,pi,2pi,qo.nodes,qo.weights)
        b=ESW.source_azimuth_bank(d,p,[r];Nq=256,mmax=m)
        h=only(ESW.source_modes(b,[HarmonicMode(m,k,orbit_frequency(p.orbit,m,k))]))
        err=source_budget_ratio(h,NamedTuple{SOURCE_KEYS}(Tuple(direct)),p)
        @test err<=1
        @printf("COORDINATE_TIME,spin=-0.7,e=0.3,budget_ratio=%.9g\n",err)
    end
    @testset "Nq refinement and measured projection cost" begin
        println("COST,spin,p,e,Nq,nr,nm,maxNphi,bank_wall_s,bank_cpu_s,sum_wall_s,sum_cpu_s,allocated_bytes,retained_bytes,workspace_bound_bytes")
        println("NQ_REFINEMENT,spin,e,Nq_to_2Nq_budget,2Nq_to_4Nq_budget")
        # Warm the same methods before timing, with no BigFloat production.
        dw,pw=source_fixture()
        bw=ESW.source_azimuth_bank(dw,pw,[80.0];Nq=2,mmax=6)
        ESW.source_modes(bw,[HarmonicMode(1,0,pw.orbit.Ωφ)])
        for (spin,e) in ((-0.7,0.3),(0.9,0.6))
            d,p=source_fixture(spin,20.0,e)
            rs=[p.orbit.p/(1+e),20.0,38.0,70.0]
            modes=[HarmonicMode(m,k,orbit_frequency(p.orbit,m,k)) for m in 0:6 for k in -4:4]
            series=[]
            # The recorded high-e pilot resolves neither 512 nor 1024;
            # retain that failure and certify two doublings from 2048.
            nqs = e>=0.5 ? (2048,4096,8192) : (512,1024,2048)
            for nq in nqs
                bt=source_measure(()->ESW.source_azimuth_bank(d,p,rs;Nq=nq,mmax=6,blocksize=2))
                b=bt.value
                ht=source_measure(()->ESW.source_modes(b,modes))
                push!(series,ht.value)
                @printf("COST,%.1f,20,%.1f,%d,4,%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%d,%d\n",
                    spin,e,nq,length(modes),maximum(b.nphi),bt.wall,bt.cpu,ht.wall,ht.cpu,
                    bt.allocated+ht.allocated,Base.summarysize(b),b.workspace_bound_bytes)
                flush(stdout)
            end
            err1=maximum(source_budget_ratio(series[1][i],series[2][i],p) for i in eachindex(series[1]))
            err2=maximum(source_budget_ratio(series[2][i],series[3][i],p) for i in eachindex(series[2]))
            @test err1<=1
            @test err2<=1
            @printf("NQ_REFINEMENT,%.1f,%.1f,%.9g,%.9g\n",spin,e,err1,err2)
        end
    end
    @testset "Circular B1 through eccentric pipeline and e to zero" begin
        println("CIRCULAR,spin,e,max_diagonal_budget,max_sideband_source_units")
        for spin in (-0.7,0.0,0.9)
            diag_errors=Float64[]
            side_errors=Float64[]
            for e in (0.0,1e-3,1e-4)
                d,p=source_fixture(spin,20.0,e)
                rs=[18.0,20.0,25.0]
                modes=[HarmonicMode(m,k,orbit_frequency(p.orbit,m,k)) for m in (1,3) for k in (-1,0,1)]
                b=ESW.source_azimuth_bank(d,p,rs;Nq=64,mmax=3)
                h=ESW.source_modes(b,modes)
                bg=ESW.DiscBackground(χ=spin,aEMRI=20.0,h0=0.06,b=0.3)
                maxbudget,maxside,maxdiag=0.0,0.0,0.0
                for (j,mode) in enumerate(modes), (ir,r) in enumerate(rs)
                    if mode.k==0
                        old=ESW.singular_field_mode(bg,mode.m,r)
                        reference=NamedTuple{SOURCE_KEYS}(ntuple(i->haskey(old,SOURCE_KEYS[i]) ? old[SOURCE_KEYS[i]] : 0.0,14))
                        maxbudget=max(maxbudget,source_budget_ratio(h[ir,j],reference,p))
                        maxdiag=max(maxdiag,abs(h[ir,j].htt-reference.htt)*20)
                    else
                        maxside=max(maxside,abs(h[ir,j].htt)*20,abs(h[ir,j].htr)*20)
                    end
                end
                if e==0
                    @test maxbudget<=1
                    @test maxside<1e-12
                end
                push!(diag_errors,maxdiag)
                push!(side_errors,maxside)
                @printf("CIRCULAR,%.1f,%.1g,%.9g,%.9g\n",spin,e,maxbudget,maxside)
            end
            @test diag_errors[3]<0.02diag_errors[2]
            @test side_errors[3]<0.11side_errors[2]
        end
    end
    @testset "FR22 direct quadrature weak-field sequence, E-ECC-2" begin
        println("FR_WEAK,spin,e,p,max_field_error,max_derivative_error,quadrature_error")
        pairs=((1,0),(2,-1),(2,1),(3,0))
        rr=[0.4,1.8]
        for e in (0.2,0.5)
            ref=fr22_direct_source(e,rr,pairs)
            refined=fr22_direct_source(e,rr,pairs;Nt=512,Nphi=1024,n=0.73)
            qerr=maximum(abs.(ref-refined))
            @test qerr<1e-10
            @test maximum(abs.(imag.(refined)))<1e-12
            errors=[]
            for porb in (1e3,1e4,1e5)
                semimajor=porb/(1-e^2)
                d,p=source_fixture(0.5,porb,e;soft=0.018semimajor)
                modes=[HarmonicMode(m,k,orbit_frequency(p.orbit,m,k)) for (m,k) in pairs]
                b=ESW.source_azimuth_bank(d,p,semimajor.*rr;Nq=256,mmax=3)
                h=ESW.source_modes(b,modes)
                err=ntuple(2) do derivative
                    maximum(eachindex(rr)) do ir
                        maximum(eachindex(modes)) do j
                            gr=derivative==1 ? -semimajor*real(h[ir,j].htt) : -semimajor^2*real(h[ir,j].htt_r)
                            abs(gr-real(refined[ir,j,derivative]))/max(abs(refined[ir,j,derivative]),1e-8)
                        end
                    end
                end
                push!(errors,err)
                @printf("FR_WEAK,0.5,%.1f,%.0f,%.9g,%.9g,%.9g\n",e,porb,err[1],err[2],qerr)
            end
            for c in 1:2
                @test errors[3][c]<0.15errors[2][c]<0.15^2*errors[1][c]
                @test errors[3][c]<2e-3
            end
        end
    end
    @testset "Compensated real and imaginary reductions" begin
        # Lander test gap (HM-NAIVE mutant): a left-to-right sum loses the unit
        # term exactly, since 1e16+1 rounds to 1e16. The bank's real and
        # imaginary reductions must both recover it through source_modes.
        _,p=source_fixture()
        K=ESW.EccentricSourceEvent(20.0,0.36,1.0,1.0,0.0,0.0,ntuple(_->1.0,7))
        kernels=Array{ComplexF64}(undef,1,3,1,2)
        for derivative in 1:2
            kernels[1,:,1,derivative]=[complex(1e16,-1e16),complex(1.0,1.0),complex(-1e16,1e16)]
        end
        b=ESW.SourceAzimuthBank(p,[20.0],zeros(3),ones(3),zeros(3),zeros(3),[K,K,K],
                                kernels,[256],[0.1],1,0)
        h=only(ESW.source_modes(b,[HarmonicMode(0,0,0.0)]))
        @test all(==(complex(1.0,1.0)),values(h))
    end
    @testset "Harmonic precision inheritance unit" begin
        setprecision(BigFloat,128) do
            d,p=source_fixture(big"0.5",big"20",big"0.3";soft=big"0.36")
            b=ESW.source_azimuth_bank(d,p,BigFloat[60];Nq=2,mmax=1)
            h=only(ESW.source_modes(b,[HarmonicMode(1,0,p.orbit.Ωφ)]))
            @test h.htr isa Complex{BigFloat}
            @test precision(real(h.htt))==128
        end
    end
end
println("PASS C_HARM independent Spacetime orbit_harmonic and Darwin coordinate-time quadrature")
println("PASS NQ_REFINEMENT two consecutive doublings; cost per orbit/radial block/mode set recorded")
println("PASS CIRCULAR_B1 exact e=0 and e->0 through eccentric source, including vanishing sidebands")
println("PASS FR_WEAK_FIELD direct double Fourier integral; E-ECC-2; field and radial derivative limit")
println("PASS HARMONIC_PRECISION BigFloat128 tiny inheritance unit; production Float64")
end

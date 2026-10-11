module EccentricDiagnosticsTests
using Test, LinearAlgebra, Random, Printf
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
for (name,file) in ((:SourceAzimuthBank,"eccentric_source.jl"),
                    (:li_boundary,"li_boundary.jl"),
                    (:EccentricModeSource,"eccentric_modes.jl"))
    isdefined(SW,name) || Base.include(SW,joinpath(@__DIR__,"../src/SpiralWaves",file))
end
include("eccentric_test_support.jl")

# Integrate actual real fields in azimuth and one temporal phase. This uses no
# complex bilinear identity or modal weight, and also covers m=0 and omega=0.
function real_field_work(mode,stress,h,measure;N=32)
    torque=0.0; power=0.0
    for j in 0:N-1, l in 0:N-1
        phi=2pi*j/N
        phase=mode.m*phi-(iszero(mode.omega) ? 0.0 : sign(mode.omega)*2pi*l/N)
        z=cis(phase)
        s=2real.(stress*z)
        dhphi=2real.(im*mode.m*h*z)
        dht=2real.(-im*mode.omega*h*z)
        torque+=sum(s.*dhphi)/2
        power-=sum(s.*dht)/2
    end
    (;torque=measure*2pi*torque/N^2,power=measure*2pi*power/N^2)
end

function context(;closure=Adiabatic,m=2,k=1,omega=0.021,q=0.75)
    d=SW.DiscModel(;χ=0.4,r_ref=30.0,h_ref=0.08,pSigma=1.5,qTemp=q,
                   closure,evolution=:fixed)
    SW.ModeContext(d,HarmonicMode(m,k,omega);closure,domain=(10.0,70.0))
end

function tensor_checks(closure;axisymmetric=false)
    rng=MersenneTwister(2801)
    for omega in (-0.021,0.021,0.0)
        axisymmetric && iszero(omega) && continue
        m=axisymmetric ? 0 : 3
        mode=HarmonicMode(m,1,omega)
        for _ in 1:8
            h=randn(rng,ComplexF64,4,4); h=h+transpose(h)
            s=randn(rng,ComplexF64,4,4); s=s+transpose(s)
            direct=SW.gravitational_work(mode,s,h;measure=3.7)
            timed=real_field_work(mode,s,h,3.7)
            @test direct.torque ≈ timed.torque rtol=2e-12 atol=1e-11
            @test direct.power ≈ timed.power rtol=2e-12 atol=1e-13
            if axisymmetric
                @test iszero(direct.torque)
                @test !iszero(direct.power)
            else
                @test direct.power ≈ (omega/m)*direct.torque rtol=1e-13 atol=1e-14
            end
        end
        c=context(;closure,m,omega)
        for r in (19.0,38.0), landau in (false,true), random_metric in (false,true)
            h=random_metric ? (;(n=>1e-4*randn(rng,ComplexF64) for n in keys(SW._zero_source(c)))...) :
                              eccentric_test_source(r)
            y=(0.003+0.001im,-0.002+0.004im)
            s=SW.gravitational_stress(c,r,y;hmodes=h,landau)
            direct=SW.gravitational_work(c,r,y;hmodes=h,landau)
            timed=real_field_work(mode,s.stress,s.metric,s.measure)
            @test direct.torque ≈ timed.torque rtol=5e-11 atol=1e-14
            @test direct.power ≈ timed.power rtol=5e-11 atol=1e-16
            @test SW.torque_density_mode(c,r,y;hmodes=h,landau)==direct.torque
            @test SW.power_density_mode(c,r,y;hmodes=h,landau)==direct.power
            if random_metric && !iszero(omega)
                raw=SW.gravitational_work(c.mode,s.eulerian,s.metric;measure=s.measure)
                @test abs(raw.power-direct.power)>1e-16
            end
            if closure==Adiabatic && !axisymmetric && !landau
                inp=SW.master_inputs(c,r;hmodes=h,landau,y)
                nt=(;(n=>getfield(inp,n) for n in fieldnames(typeof(inp)))...)
                bars=(;(Symbol(string(n),"c")=>conj(getproperty(nt,n)) for n in
                       (:htt,:htphi,:hphiphi,:hrr,:hzz,:htr,:hrphi,:y0,:y1))...)
                oracle=4pi*s.measure*real(SW._eccentric_adiabatic_torque_mode(merge(nt,bars)))
                @test direct.torque ≈ oracle rtol=2e-11 atol=1e-13
            end
            # Lander review R1 (w10-g28 landing): the whole adiabatic stress, not only its averaged
            # contractions. stress-W:h must equal the R3-emitted delta_xi T of G26 (40)-(42) in every
            # component (W of G26 (42) built here), W:h must do no averaged work (the docstring's
            # equivalence argument, executed), and the context AMF must equal the emitted (72a)
            # kernel. Mutants M09 (sign of -p g^-1 h g^-1, zero averaged work) and M11 (pressure
            # current dropped from the context AMF) survived the default selections without these.
            if closure==Adiabatic && !landau
                inp=SW.master_inputs(c,r;hmodes=h,landau,y)
                geo=SW._diagnostic_geometry(c,r,h)
                q=geo.gi+geo.u*transpose(geo.u)
                W=(a,b,g,d)->((inp.eps+inp.p)*geo.u[a]*geo.u[b]*geo.u[g]*geo.u[d]+
                            inp.p*(geo.gi[a,b]*geo.gi[g,d]-geo.gi[a,g]*geo.gi[b,d]-geo.gi[a,d]*geo.gi[b,g])-
                            inp.Gamma1*inp.p*q[a,b]*q[g,d])/2
                Wh=[sum(W(a,b,g,d)*s.metric[g,d] for g in 1:4, d in 1:4) for a in 1:4, b in 1:4]
                scale=maximum(abs,s.stress)
                for (name,(a,b)) in pairs((tt=(1,1),tr=(1,2),tphi=(1,3),rr=(2,2),rphi=(2,3),phiphi=(3,3),zz=(4,4)))
                    emitted=getfield(SW,Symbol("_eccentric_adiabatic_dxiT_",name))(inp)
                    @test s.stress[a,b]-Wh[a,b] ≈ emitted rtol=1e-10 atol=1e-12*scale
                    @test s.stress[b,a]-Wh[b,a] ≈ emitted rtol=1e-10 atol=1e-12*scale
                end
                @test all(iszero,s.stress[4,1:3]) && all(iszero,s.stress[1:3,4])
                contact=real_field_work(mode,Wh,s.metric,s.measure)
                bound=1e-12*s.measure*maximum(abs,Wh)*maximum(abs,s.metric)*(mode.m+abs(mode.omega))*16
                @test abs(contact.torque)<=bound && abs(contact.power)<=bound
                nt=(;(n=>getfield(inp,n) for n in fieldnames(typeof(inp)))...)
                bars=(;(Symbol(string(n),"c")=>conj(getproperty(nt,n)) for n in
                       (:htt,:htphi,:hphiphi,:hrr,:hzz,:htr,:hrphi,:y0,:y1))...)
                amf_oracle=4pi*s.measure*real(SW._eccentric_adiabatic_amf_mode(merge(nt,bars)))
                @test SW.amf_mode(c,r,y;hmodes=h,landau) ≈ amf_oracle rtol=2e-11 atol=1e-13
            end
        end
    end
    # An elementary Newtonian density/potential pair fixes the sign into the
    # disc and the half-amplitude map to FR24 v1 (2_linear_framework.tex:268-273).
    mode=HarmonicMode(axisymmetric ? 0 : 2,1,0.03)
    s=zeros(ComplexF64,4,4); h=copy(s)
    s[1,1]=(1+2im)/2; h[1,1]=-0.7 # h_tt=-2 Phi_hat=-Phi_FR
    v=SW.gravitational_work(mode,s,h;measure=5.0)
    @test v.torque ≈ -pi*5*mode.m*0.7*2
    @test v.power ≈ -pi*5*mode.omega*0.7*2
    println("PASS ",closure,axisymmetric ? " axisymmetric" : " C-POWER",
            ": direct tensor, signed/zero frequency, independent real-field work, sign into disc")
end

function adiabatic_checks()
    @testset "Adiabatic direct work and legacy" begin
        tensor_checks(Adiabatic)
        include("test_eccentric_legacy.jl")
        bg=SW.DiscBackground(χ=0.4,aEMRI=20.0,h0=0.06,𝔭=1.5,𝔮=0.75)
        values=Float64[]
        for r in (9.0,17.0,23.0,61.0)
            h=SW.singular_field_mode(bg,3,r;N=128); y=(0.3+0.1im,-0.2+0.4im)
            push!(values,SW.amf_mode(bg,3,r,y;hmodes=h,rc=20.0),
                  SW.torque_density_mode(bg,3,r,y;hmodes=h,rc=20.0))
        end
        @test reinterpret(UInt64,values)==UInt64[0x409a85b7ef0cd74c,0x4085292221778c6b,
            0x40c5897e6acf4ab8,0x40b7ec04227d6e83,0xc0cae29f3e97a33f,0x40b97541fd8fb6bf,
            0xc0a0f7a3a3ff55fe,0x405393bd7efe8a20]
        forced_check(Adiabatic)
    end
    println("PASS legacy diagnostic bit-identity: eight baseline Float64 flux/torque values")
end

function forced_check(closure; a=0.9,e=0.2,m=2,k=1,n=128,Nq=64,domain=(12.0,75.0),adaptive=false)
    d=SW.DiscModel(;χ=a,r_ref=30.0,h_ref=0.08,pSigma=1.5,qTemp=0.75,
                   closure,evolution=:fixed,εp=1e-3)
    pert=SW.EccentricPerturber(EquatorialOrbit(a,30.0,e);mass=1e-5,softening=FixedSoftening(0.72))
    c=SW.ModeContext(d,pert,m,k;closure,domain)
    grid=if adaptive
        SW.build_grid(c;nL=1,nR=1,nC=65,n_anc=64)
    else
        r=collect(range(domain...;length=n+1))
        SW.EccentricGrid(r,c.r_c,1:0,SW.mode_cache_key(c),(;diagnostic_sanity=true))
    end
    r=grid.r
    sol=SW.solve_mode(c,pert;grid,bank_options=(;Nq))
    @test sol.residual<1e-8 && sol.bc_residual<1e-10
    @test all(f->all(isfinite,values(f)),SW.reconstruct_fluid(sol))
    work=SW.gravitational_work(sol)
    density=getproperty.(work,:torque); power=getproperty.(work,:power)
    torque=SW.total_torque(sol;density)
    P=SW.total_power(sol;density=power)
    @test all(isfinite,(torque,P)) && !iszero(P)
    @test SW.torque_density(sol)==density && SW.power_density(sol)==power
    @test SW.total_torque(sol)==torque && SW.total_power(sol)==P
    @test last(SW.cumulative_torque(sol;density=ones(Int,length(r)))) ≈ last(r)-first(r)
    @test last(SW.cumulative_power(sol;density=copy(r))) ≈ (last(r)^2-first(r)^2)/2
    @test_throws DimensionMismatch SW.cumulative_torque(sol;density=zeros(length(r)-1))
    @test all(isfinite,SW.amf(sol))
    if closure==Adiabatic
        @test SW.lindblad_torque(sol)==last(SW.amf(sol))-first(SW.amf(sol))
    else
        @test_throws ArgumentError SW.lindblad_torque(sol)
        @test getproperty.(SW.li_flux_decomposition(sol),:wave)==SW.amf(sol)
    end
    if m==0
        @test iszero(torque) && all(iszero,density)
    else
        @test P ≈ c.mode.omega/m*torque rtol=2e-12
    end
    i=cld(length(r),2)
    stress=SW.gravitational_stress(c,r[i],sol.y[:,i];hmodes=sol.source(r[i]))
    timed=real_field_work(c.mode,stress.stress,stress.metric,stress.measure)
    @test power[i] ≈ timed.power rtol=1e-11 atol=1e-24
    @printf("SANITY,%s,%.2f,30,%.3f,%d,%d,%.12e,%.12e,%.12e,%d,%d,%.3e,%s\n",
            string(closure),a,e,m,k,c.mode.omega,torque,P,length(r)-1,Nq,sol.residual,adaptive ? "adaptive" : "uniform")
    println("PASS forced ",closure," (",m,",",k,"): end-to-end solve/reconstruction, direct work and quadrature")
    sol
end

function circular_checks(closure)
    @testset "Summed m-power circular continuity $closure" begin
        bg=SW.DiscBackground(χ=0.9,aEMRI=30.0,h0=0.08,𝔭=1.5,𝔮=0.0,b=0.3,εp=1e-3)
        d=SW.DiscModel(bg;closure,evolution=:fixed)
        reference=0.0; totals=zeros(3)
        for m in (2,3)
            cg=SW.build_grid(bg,m;nL=1,nR=1,nC=65,n_anc=64,r_in=12.0,r_out=75.0)
            circ=SW.solve_mode(bg,m;grid=cg,mEMRI=1e-5,linsolve=:banded)
            reference+=SW.Omega_emri(bg)*SW.total_torque(bg,circ)
            for (j,e) in enumerate((0.0,1e-3,1e-4)), k in (-1,0,1)
                pert=SW.EccentricPerturber(EquatorialOrbit(0.9,30.0,e);mass=1e-5,
                                          softening=FixedSoftening(SW.softening_length(bg)))
                c=SW.ModeContext(d,pert,m,k;closure,domain=(first(cg.r),last(cg.r)))
                grid=SW.EccentricGrid(copy(cg.r),c.r_c,1:0,SW.mode_cache_key(c),(;continuity=true))
                # qTemp=0 closure overlap also requires identical boundary
                # conditions. LI's default includes a further WKB correction.
                ad=SW.DiscModel(bg;closure=Adiabatic,evolution=:fixed)
                ca=SW.ModeContext(ad,c.mode;closure=Adiabatic,domain=c.domain)
                rows=SW.mode_boundary_rows(ca)
                sol=SW.solve_mode(c,pert;grid,bank_options=(;Nq=16),boundary_rows=rows)
                totals[j]+=SW.total_power(sol)
                if e==0 && k==0
                    @test norm(sol.y-circ.y)/norm(circ.y)<2e-7
                    @test SW.torque_density(sol) ≈ SW.torque_density(bg,circ) rtol=2e-7
                    # Frozen _r2_amf treats the Doppler frequency as real under
                    # conjugation. Compare the physical formula where Landau
                    # is disabled; keep the legacy shifted output unchanged.
                    if closure==Adiabatic
                        i=1; r=first(sol.r)
                        hc=SW.singular_field_mode(bg,m,r;mEMRI=1e-5)
                        @test SW.amf_mode(c,r,sol.y[:,i];hmodes=sol.source(r),landau=false) ≈
                              SW.amf_mode(bg,m,r,circ.y[:,i];hmodes=hc,landau=false,rc=c.r_c) rtol=2e-7
                    end
                end
            end
        end
        errors=abs.(totals./reference.-1)
        @test errors[1]<2e-7
        @test errors[3]<errors[2] && 50<errors[2]/errors[3]<150
        @printf("CONTINUITY,%s,%.12e,%.12e,%.12e,%.12e,%.6e,%.6e,%.6e\n",
                string(closure),reference,totals...,errors...)
    end
    println("PASS ",closure,": e-to-zero m=2:3 k=-1:1 summed direct power against circular 2608; qTemp=0 matched closure")
end

function sanity_checks(closure)
    @testset "Exploratory GR eccentric physics table $closure" begin
        for (a,e) in ((0.9,0.2),(0.0,0.2),(0.9,0.05)), (m,k) in ((2,0),(2,1),(0,1))
            forced_check(closure;a,e,m,k,adaptive=true)
        end
    end
    println("PASS sanity cells are finite, satisfy direct work/torque and retain axisymmetric power; no mesh, tail or exterior convergence claim")
end

end # module

let selection=get(ENV,"ECC_DIAGNOSTICS_TEST","all")
    selection in ("all","adiabatic") && EccentricDiagnosticsTests.adiabatic_checks()
    selection=="circular" && EccentricDiagnosticsTests.circular_checks(EccentricDiagnosticsTests.Adiabatic)
    selection=="sanity" && EccentricDiagnosticsTests.sanity_checks(EccentricDiagnosticsTests.Adiabatic)
end

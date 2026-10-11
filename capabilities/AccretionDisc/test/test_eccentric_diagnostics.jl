module EccentricDiagnosticsTests
using Test, LinearAlgebra, Random, Printf
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
for (name,file) in ((:SourceAzimuthBank,"eccentric_source.jl"),
                    (:li_boundary,"li_boundary.jl"),
                    (:li_flux_decomposition,"li_diagnostics.jl"),
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
        # Lander review R2 (finding of lane w11-g30): at e = 0 the context AMF equals the frozen
        # circular AMF THROUGH corotation with the Landau shift on, at the original regulator
        # epsilon_p = 1e-5: the conjugated factor is continued in omega_d (R3 formal bar), not
        # literally conjugated (up to 28% apart inside the layer before the fix).
        bc=SW.DiscBackground(χ=0.9,aEMRI=30.0,h0=0.08,𝔭=1.5,𝔮=0.0,b=0.3)
        @test bc.εp==1e-5
        cc=SW.ModeContext(SW.DiscModel(bc;closure=Adiabatic,evolution=:fixed),
                          HarmonicMode(2,0,2SW.Omega_emri(bc));closure=Adiabatic,domain=(12.0,75.0))
        shifts=Float64[]
        for δ in (-1e-3,-1e-4,-3e-5,3e-5,1e-4,1e-3,1e-2)
            r=cc.r_c+δ; y=(0.3+0.1im,-0.2+0.4im)
            h=SW.singular_field_mode(bc,2,r;N=128); he=merge(map(complex,h),(htr=0im,hrphi=0im))
            inp=SW.master_inputs(cc,r;hmodes=he,landau=true,y)
            push!(shifts,imag(inp.omega_d)/abs(real(inp.omega_d)))
            @test SW.amf_mode(cc,r,y;hmodes=he,landau=true) ≈
                  SW.amf_mode(bc,2,r,y;hmodes=h,landau=true,rc=cc.r_c) rtol=1e-12
        end
        @test maximum(shifts)>1   # the sampled radii reach inside the Landau layer
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
                        # Lander review R2: the solved profiles also agree through corotation with
                        # the Landau shift on (the context AMF continues omega_d as the frozen kernel).
                        Fc=SW.amf(bg,circ); Fe=SW.amf(sol)
                        @test maximum(abs.(Fe-Fc))<=1e-6*maximum(abs,Fc)
                        @printf("AMF_PROFILE,m=%d,%.3e\n",m,maximum(abs.(Fe-Fc))/maximum(abs,Fc))
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

function axisymmetric_checks(closure)
    @testset "Dynamic axisymmetric direct power $closure" begin
        tensor_checks(closure;axisymmetric=true)
        for mode in (HarmonicMode(0,0,0.0),HarmonicMode(0,0,0.1),HarmonicMode(0,-1,-0.1))
            @test_throws ArgumentError SW.reality_pair_weight(mode)
            @test_throws ArgumentError SW.gravitational_work(mode,zeros(4,4),zeros(4,4))
        end
        @test SW.reality_pair_weight(HarmonicMode(0,1,0.1))==4pi
        @test SW.reality_pair_weight(HarmonicMode(2,-4,-0.1))==4pi
        @test SW.reality_pair_weight(HarmonicMode(2,-4,0.0))==4pi
        @test SW._fourier_weight(Float64,0)==2pi
        forced_check(closure;m=0,k=1)
    end
    println("PASS ",closure,": T_0k=0; dynamic weight 4pi; static (0,0) excluded; no division by m or omega")
end

function li_flux_checks()
    c=context(closure=LocallyIsothermalOptionA)
    for r in (19.0,38.0)
        h=eccentric_test_source(r); y=(0.003+0.001im,-0.002+0.004im)
        f=SW.reconstruct_fluid(c,r,y;hmodes=h,landau=false)
        geo=SW._diagnostic_geometry(c,r,h)
        du=[f.dut,f.dur,f.duphi_up,0im]
        @test f.du_phi ≈ (geo.g*du+geo.metric*geo.u)[3] rtol=1e-12
        # Independently differentiate the exact perfect-fluid stress jet.
        enthalpy=SW.eps_profile(c.disc,r)+SW.pressure(c.disc,r)
        eulerian=(f.deps+f.dp)*(geo.u*transpose(geo.u))+
            enthalpy*(du*transpose(geo.u)+geo.u*transpose(du))+
            f.dp*geo.gi-SW.pressure(c.disc,r)*geo.gi*geo.metric*geo.gi
        @test SW.linear_stress(c,r,y;hmodes=h,landau=false).dT ≈ eulerian rtol=2e-12
        parts=SW.li_flux_decomposition(c,r,y;hmodes=h,landau=false)
        @test parts.quadratic-parts.mean_transport ≈ parts.wave rtol=1e-12
        @test !iszero(parts.mean_transport)
        @test parts.wave==SW.amf_mode(c,r,y;hmodes=h,landau=false)
        no_cross=SW.reynolds_flux(c.mode;measure=geo.measure,enthalpy,dur=f.dur,du_phi=(geo.g*du)[3])
        @test abs(no_cross-parts.wave)>1e-8*abs(parts.wave)
        # Extract the quadratic term of the actual nonlinear mixed stress and
        # subtract the background angular momentum transported by the mean.
        δ=1e-3; full=0.0; mean=0.0; u_phi=(geo.g*geo.u)[3]
        for j in 0:63, sign in (-1,1)
            z=cis(2pi*j/64); λ=sign*δ
            g=geo.g+λ*2real.(geo.metric*z)
            u=geo.u+λ*2real.(du*z)
            e=SW.eps_profile(c.disc,r)+λ*2real(f.deps*z)
            p=SW.pressure(c.disc,r)+λ*2real(f.dp*z)
            stress=(e+p)*(u*transpose(u))+p*inv(g)
            full+=sqrt(-det(g))*(stress*g)[2,3]
            mean+=sqrt(-det(g))*(e+p)*u[2]*u_phi
        end
        full*=2pi/(128δ^2); mean*=2pi/(128δ^2)
        @test full ≈ parts.quadratic rtol=2e-6
        @test mean ≈ parts.mean_transport rtol=2e-6
        @test full-mean ≈ parts.wave rtol=2e-6
    end
    # Direct convention recovery, with arbitrary independent FR amplitudes.
    mode=HarmonicMode(2,1,0.03)
    for r in (0.7,1.3,8.0)
        ur=0.3+0.2im; up=-0.1+0.4im; Sigma=r^(-1.5)
        F=SW.reynolds_flux(mode;measure=r,enthalpy=Sigma,dur=ur/2,du_phi=r*up/2)
        @test F ≈ pi*r^2*Sigma*real(ur*conj(up)) rtol=5e-15
    end
    # Lander review R2: one conjugation convention for every closure. At qTemp = 0 the OPTION A and
    # adiabatic reconstructions coincide (C-LI-0), so through corotation with the Landau shift on the
    # LI wave flux plus the adiabatic pressure current u_phi<dp du^r> equals the adiabatic context AMF.
    ca=context(;closure=Adiabatic,q=0.0); cl=context(;closure=LocallyIsothermalOptionA,q=0.0)
    @test ca.r_c==cl.r_c
    shifts=Float64[]
    for δ in (-1e-3,-1e-4,-3e-5,3e-5,1e-4,1e-3)
        r=ca.r_c+δ; h=eccentric_test_source(r); y=(0.003+0.001im,-0.002+0.004im)
        inp=SW.master_inputs(cl,r;hmodes=h,landau=true,y)
        push!(shifts,imag(inp.omega_d)/abs(real(inp.omega_d)))
        fa=SW.reconstruct_fluid(ca,r,y;hmodes=h,landau=true)
        pc=4pi*SW.sqrt_minus_g(ca.disc,r)*real(SW.rotation_profile(ca.disc,r).uφ*fa.dp*conj(fa.dur))
        @test SW.li_flux_decomposition(cl,r,y;hmodes=h,landau=true).wave+pc ≈
              SW.amf_mode(ca,r,y;hmodes=h,landau=true) rtol=1e-10
    end
    @test maximum(shifts)>1
    println("PASS LI Reynolds: emitted stress jets, exact mixed-stress expansion, metric cross terms and mean-transport subtraction; FR24 pi r^2 Sigma Re(ur conj(uphi))")
end

# Independent coordinate-current divergence, retaining the perturbed measure.
# It tests zero-bath equations, separately from the algebraic C-POWER identity.
function current_density(c,r,y,h)
    geo=SW._diagnostic_geometry(c,r,h)
    stress=SW.gravitational_stress(c,r,y;hmodes=h,landau=false)
    p=SW.pressure(c.disc,r); e=SW.eps_profile(c.disc,r)
    background=(e+p)*(geo.u*transpose(geo.u))+p*geo.gi
    (;current=geo.measure*(stress.stress*geo.g+background*geo.metric),
      background,metric=geo.metric,measure=geo.measure)
end

function zero_bath_checks(;axisymmetric=false)
    for q in (0.0,0.75,1.0), omega in (-0.021,0.021)
        c=context(;closure=LocallyIsothermalOptionA,m=axisymmetric ? 0 : 2,omega,q)
        r=38.0; y=[0.003+0.001im,-0.002+0.004im]; h=eccentric_test_source(r)
        C,S=SW.master_matrices(c,r;hmodes=h,landau=false); yp=C*y+S
        j=current_density(c,r,y,h); errors=Float64[]
        for δ in (0.02,0.01,0.005)
            jp=current_density(c,r+δ,y+δ*yp,eccentric_test_source(r+δ))
            jm=current_density(c,r-δ,y-δ*yp,eccentric_test_source(r-δ))
            dj=(jp.current-jm.current)/(2δ)
            terms=ComplexF64[]; scales=Float64[]
            for (a,w) in ((1,-im*omega),(3,im*c.mode.m))
                td=-im*omega*j.current[1,a]; pd=im*c.mode.m*j.current[3,a]
                metricwork=j.measure*sum(j.background .* (w*j.metric))/2
                push!(terms,td+pd+dj[2,a]-metricwork)
                push!(scales,abs(td)+abs(pd)+abs(dj[2,a])+abs(metricwork))
            end
            push!(errors,maximum(abs.(terms)./scales))
        end
        @test errors[end]<2e-7
        @test 3.5<errors[1]/errors[2]<4.5 && 3.5<errors[2]/errors[3]<4.5
        @printf("ZERO_BATH,m=%d,q=%.2f,omega=%+.3f,%.6e,%.6e,%.6e\n",c.mode.m,q,omega,errors...)
    end
    println("PASS OPTION A zero bath: independent densitized E/J-current divergence with metric work; qTemp=0,0.75,1; no constant-FJ assumption")
end

function newtonian_flux_checks()
    MO=AccretionDisc.Modes
    for q in (0.0,0.75,1.0)
        d=MO.NewtonianDisc(h0=0.12,p=1.5,q=q,beta=0.0,pole_disp=0.0)
        m=2; mode=HarmonicMode(m,0,Float64(m)); errors=Float64[]; fjchange=0.0
        for n in (128,256,512)
            rs=collect(range(1.5,2.0;length=n+1))
            # Reuse the certified zero-cooling Newtonian operator, ih=0.
            Cs=map(1:n) do i
                r=(rs[i]+rs[i+1])/2
                col0=MO.master_rhs(d,r,(1.,0.,0.,0.),m,0)
                col1=MO.master_rhs(d,r,(0.,0.,1.,0.),m,0)
                ComplexF64[0 1; col0[3]+im*col0[4] col1[3]+im*col1[4]]
            end
            y=SW.solve_midpoint(rs,Cs,[zeros(ComplexF64,2) for _ in 1:n],(1,0),1,
                                 (-MO.wkbC(d,last(rs),m),1),0)
            F=map(eachindex(rs)) do i
                r=rs[i]; Ω=MO.Omega(d,r); σ=m-m*Ω; D=MO.Dfreq(d,r,m)
                # MR20 src/main.tex:925–931, read as text; physical FR amplitudes.
                ur=im/D*(σ*(y[2,i]-MO.LTinv(d,r)*y[1,i])-2m*Ω/r*y[1,i])
                up=SW._newtonian_delta_u_phi(m,r,Ω,MO.kappa(d,r)^2,σ,(y[1,i],-im*ur),0)
                flux=SW.reynolds_flux(mode;measure=r,enthalpy=MO.Sigma(d,r),dur=ur/2,du_phi=r*up/2)
                wronskian=pi*m*r*MO.Sigma(d,r)/D*imag(y[1,i]*conj(y[2,i]))
                @test flux ≈ wronskian rtol=2e-12
                flux
            end
            law=F./MO.csiso.(Ref(d),rs).^2
            push!(errors,maximum(abs.(law./first(law).-1)))
            fjchange=last(F)/first(F)-1
        end
        @test all(3.7 .< errors[1:end-1]./errors[2:end] .< 4.3)
        @test errors[end]<5e-5
        q>0 && @test abs(fjchange)>0.1
        @printf("NEWTONIAN_FJ_CT2,q=%.2f,%.6e,%.6e,%.6e,FJ_change=%+.6e\n",q,errors...,fjchange)
    end
    println("PASS MR20 source-free law: Reynolds FJ equals pressure Wronskian; FJ/cT2 conserved at second order, FJ varies for qTemp!=0")
end

function weak_field_flux_checks()
    MO=AccretionDisc.Modes
    for q in (0.0,0.75,1.0)
        errors=Float64[]; x=1.7; m=2; omega=3.0
        dn=MO.NewtonianDisc(h0=0.08,p=1.5,q=q,beta=0.0,pole_disp=0.0)
        yn=(0.002+0.001im,0.003-0.002im)
        up=SW._newtonian_delta_u_phi(m,x,MO.Omega(dn,x),MO.kappa(dn,x)^2,
                                    omega-m*MO.Omega(dn,x),yn,0)
        target=pi*x^2*MO.Sigma(dn,x)*real((2im*yn[2])*conj(2up))
        for R in (2e3,2e4,2e5)
            d=SW.DiscModel(;χ=0.9,r_ref=R,h_ref=0.08,pSigma=1.5,qTemp=q,
                           closure=LocallyIsothermalOptionA,evolution=:fixed)
            c=SW.ModeContext(d,HarmonicMode(m,1,omega/R^(3/2));
                             closure=LocallyIsothermalOptionA,domain=(R,2R))
            y=(yn[1]/R,yn[2]/sqrt(R))
            F=SW.amf_mode(c,x*R,y;hmodes=SW._zero_source(c),landau=false)/R
            push!(errors,abs(F/target-1))
        end
        @test all(8 .< errors[1:end-1]./errors[2:end] .< 12)
        @test last(errors)<1e-4
        @printf("WEAK_FIELD_FLUX,q=%.2f,%.6e,%.6e,%.6e\n",q,errors...)
    end
    println("PASS emitted GR LI Reynolds flux approaches the FR24 formula at first order in M/r; half amplitudes and covariant r factor retained")
end

function li_checks()
    @testset "LI stress, Reynolds flux and independent balance" begin
        tensor_checks(LocallyIsothermalOptionA)
        li_flux_checks()
        newtonian_flux_checks()
        weak_field_flux_checks()
        zero_bath_checks()
        forced_check(LocallyIsothermalOptionA)
    end
end

function sanity_checks(closure)
    @testset "Exploratory GR eccentric physics table $closure" begin
        for (a,e) in ((0.9,0.2),(0.0,0.2),(0.9,0.05)), (m,k) in ((2,0),(2,1),(0,1))
            forced_check(closure;a,e,m,k,adaptive=true)
        end
    end
    println("PASS sanity cells are finite, satisfy direct work/torque and retain axisymmetric power; no mesh, tail or exterior convergence claim")
end

function li_axisymmetric_checks()
    axisymmetric_checks(LocallyIsothermalOptionA)
    @testset "LI axisymmetric zero bath" begin
        zero_bath_checks(;axisymmetric=true)
    end
end

end # module

let selection=get(ENV,"ECC_DIAGNOSTICS_TEST","all")
    selection in ("all","adiabatic") && EccentricDiagnosticsTests.adiabatic_checks()
    selection in ("all","axisymmetric") && EccentricDiagnosticsTests.axisymmetric_checks(EccentricDiagnosticsTests.Adiabatic)
    selection=="circular" && EccentricDiagnosticsTests.circular_checks(EccentricDiagnosticsTests.Adiabatic)
    selection in ("all","li") && EccentricDiagnosticsTests.li_checks()
    selection in ("all","li_axisymmetric") && EccentricDiagnosticsTests.li_axisymmetric_checks()
    selection=="li_circular" && EccentricDiagnosticsTests.circular_checks(EccentricDiagnosticsTests.LocallyIsothermalOptionA)
    selection=="sanity" && EccentricDiagnosticsTests.sanity_checks(EccentricDiagnosticsTests.Adiabatic)
    selection=="li_sanity" && EccentricDiagnosticsTests.sanity_checks(EccentricDiagnosticsTests.LocallyIsothermalOptionA)
end

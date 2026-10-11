using Test, AccretionDisc, NumericsCore, LinearAlgebra

@testset "regular Lindblad crossing" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    @testset "C2 forcing and master-equation identity" begin
        # Phi(log10 r) is a cubic: the global spline and both chain-rule
        # derivatives have exact independent values, including at its knots.
        x=10.0 .^ range(-1,1;length=30)
        phi(r)=log10(r)^3-2log10(r)+1
        phi1(r)=(3log10(r)^2-2)/(r*log(10.))
        phi2(r)=6log10(r)/(r*log(10.))^2-phi1(r)/r
        c=PotentialCache(log10.(x),phi.(x),phi1.(x),phi2.(x))
        rc=P.RegularPotential(c)
        for r in [x; .37; 1.9]
            @test collect(P.potential_at(rc,r)) ≈ [phi(r),phi1(r),phi2(r)] rtol=1e-11 atol=1e-11
        end
        for q in (0.,.75,1.5), (m,l) in ((1,2),(2,1),(3,5)), pole in (0.,1e-6,.001)
            d=NewtonianDisc(h0=.08,p=.5,q=q,beta=0.,gamma=1.4,pole_disp=pole)
            mode=F.NewtonianMode(m,l,.8); r=.61
            y=[.7,.1,.2,-.3]; rhs=P.RegularRHS(d,mode,1,rc); dy=rhs(r,y)
            O=P.Omega(d,r); O1=P.dOmegadr(d,r); S=P.Sigma(d,r)
            w=mode.omega-m*O+im*pole; w1=-m*O1; D=P.Dfreq(d,r,mode)
            L=P.LTinv(d,r); L1=P.dLTinvdr(d,r)
            a=2m*O/(r*w)+L; b=-im*D/(w*r*S)
            a1=2m*(O1/r-O/r^2)/w-2m*O/(r*w)*w1/w+L1
            b1=-im/(w*r*S)*(P.dDdr(d,r,mode)-D*(w1/w+(1-d.p)/r))
            psi=complex(y[1],y[2]); flux=complex(y[3],y[4])
            psi1=complex(dy[1],dy[2]); flux1=complex(dy[3],dy[4])
            h2=a1*psi+b1*flux-L1*phi(r)-L*phi1(r)+a*psi1+b*flux1-phi2(r)
            h=psi-phi(r); h1=psi1-phi1(r)
            # Use the SAME smooth forcing in the independently landed master.
            coeff=P.master_equation_cool(r,d.gamma,P.csadi(d,r),O,O1,S,P.dSigmadr(d,r),
                L,L1,P.LSinv(d,r),P.dLSinvdr(d,r),real(w),0.,0.,D,P.dDdr(d,r,mode),
                P.Nr2(d,r),P.dNr2dr(d,r),pole,phi(r),phi1(r),phi2(r),m)
            expected=complex(coeff[5],coeff[6])-complex(coeff[1],coeff[2])*h1-complex(coeff[3],coeff[4])*h
            @test h2 ≈ expected rtol=2e-11 atol=2e-10
        end
    end
    @testset "exact D=0, harmonic frequencies and m=1 without ILR" begin
        d=NewtonianDisc(h0=.1,p=0.,q=0.,beta=0.,gamma=1.,soft=.03,lap_nphi=128)
        for (m,l) in ((2,1),(3,2),(1,2))
            mode=F.NewtonianMode(m,l,1.)
            roots=m==1 ? [((m+1)/l)^(2/3)] : [((m-1)/l)^(2/3),((m+1)/l)^(2/3)]
            x=sort!(unique!([collect(loggrid(d,.2,4.,701)); roots]))
            c=PotentialCache(d,m,first(x),last(x);n=501)
            for r in roots
                @test abs(P.Dfreq(d,r,mode)) <= 100eps()*P.kappa(d,r)^2
                @test all(isfinite,P.RegularRHS(d,mode,1,P.RegularPotential(c))(r,[1.,.1,.2,.3]))
            end
            if m==1
                @test all(P.Dfreq(d,r,mode)>0 for r in range(.01,roots[1]-.001;length=500))
                @test length(roots)==1
            end
            opts=(;lindblad_regular=true,rtol=1e-10,atol=1e-12,dense=false)
            s=solve_lm_mode(d,mode,x,c;opts...)
            @test size(s)==(length(x),3)
            @test all(isfinite,s)
            @test all(v->all(isfinite,v),vars_from_enthalpy(d,s,x,mode,c;lindblad_regular=true))
        end
        x=loggrid(d,.2,4.,101); c=PotentialCache(d,2,.2,4.;n=201); mode=F.NewtonianMode(2,2,1.)
        @test isequal(solve_lm_mode(d,mode,x,c),solve_lm_mode(d,mode,x,c;lindblad_regular=false))
        @test_throws ArgumentError solve_lm_mode(NewtonianDisc(beta=.1),mode,x,c;lindblad_regular=true)
        @test_throws ArgumentError vars_from_enthalpy(d,zeros(ComplexF64,101,2),x,mode,c;lindblad_regular=true)
        coarse=[.2,1.]; coarse_mode=F.NewtonianMode(1,2,1.)
        @test_throws ArgumentError solve_lm_mode(d,coarse_mode,coarse,c;lindblad_regular=true)
    end
end

# Lander review additions (ECC::fr22-linear-solver landing of lane w16-g40), test-only except L5's named settings.
@testset "regular Lindblad crossing (lander review)" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    rel(a,b)=norm(a-b)/norm(b)
    @testset "L1 the regular system satisfies FR22 eqs (13)-(15) at pole_disp=0" begin
        # 2_framework.tex:126-138 (read as text): continuity, radial and azimuthal momentum with
        # dP = c_iso^2 dSigma, dh = dP/Sigma and P = c_iso^2 Sigma; independent of master_equation_cool.
        x=10.0 .^ range(-1,1;length=30)
        phi(r)=log10(r)^3-2log10(r)+1
        phi1(r)=(3log10(r)^2-2)/(r*log(10.))
        phi2(r)=6log10(r)/(r*log(10.))^2-phi1(r)/r
        rc=P.RegularPotential(PotentialCache(log10.(x),phi.(x),phi1.(x),phi2.(x)))
        for q in (0.,.75), (m,l) in ((1,2),(2,1),(3,5)), r in (.37,.61,1.9)
            d=NewtonianDisc(h0=.08,p=.5,q=q,beta=0.,gamma=1.4,pole_disp=0.)
            mode=F.NewtonianMode(m,l,.8); y=[.7,.1,.2,-.3]; dy=P.RegularRHS(d,mode,1,rc)(r,y)
            psi=complex(y[1],y[2]); flux=complex(y[3],y[4]); psi1=complex(dy[1],dy[2]); flux1=complex(dy[3],dy[4])
            O=P.Omega(d,r); S=P.Sigma(d,r); S1=P.dSigmadr(d,r); c2=P.csiso(d,r)^2; k2=P.kappa(d,r)^2
            w=mode.omega-m*O; h=psi-phi(r); h1=psi1-phi1(r)
            ur=flux/(r*S); uphi=(m*psi/r-im*k2*ur/(2O))/w
            dSig=S*h/c2; dP=S*h; dP1=S1*h+S*h1; P1=c2*S*(P.LTinv(d,r)+S1/S)
            r13=-im*w*dSig+flux1/r+im*m*S*uphi/r
            r14=-im*w*ur-2O*uphi+dP1/S-P1*dSig/S^2+phi1(r)
            r15=-im*w*uphi+k2*ur/(2O)+im*m/r*(dP/S+phi(r))
            @test abs(r13) <= 1e-12*(abs(w*dSig)+abs(flux1)/r+m*S*abs(uphi)/r)
            @test abs(r14) <= 1e-12*(abs(w*ur)+2O*abs(uphi)+abs(dP1)/S+abs(P1*dSig)/S^2+abs(phi1(r)))
            @test abs(r15) <= 1e-12*(abs(w*uphi)+k2*abs(ur)/(2O)+m/r*(abs(dP/S)+abs(phi(r))))
        end
    end
    @testset "L2 regular path equals the master-ODE path where D != 0" begin
        d=NewtonianDisc(h0=.08,q=1.,p=1.5,beta=0.,gamma=1.4,soft=.1,pole_disp=1e-6,lap_nphi=128)
        mode=F.NewtonianMode(2,3,1.); x=loggrid(d,.55,.95,600); c=PotentialCache(d,2,.5,1.;n=801)
        @test all(r->P.Dfreq(d,r,mode)>0,x)
        opts=(;anchor=:corotation,rtol=1e-12,atol=1e-14,dense=false)
        o=solve_lm_mode(d,mode,x,c;opts...); n=solve_lm_mode(d,mode,x,c;opts...,lindblad_regular=true)
        vo=vars_from_enthalpy(d,o,x,mode,c); vn=vars_from_enthalpy(d,n,x,mode,c;lindblad_regular=true)
        @test rel(n[:,1],o[:,1]) <= 1e-6
        @test rel(n[:,2],o[:,2]) <= 1e-6
        @test rel(n[:,3],x.*P.Sigma.(Ref(d),x).*vo[1]) <= 1e-6
        for j in 1:3
            @test rel(vn[j],vo[j]) <= 1e-6
        end
    end
    @testset "L3 the C2 forcing spline reproduces every cache node" begin
        lg=collect(range(-1,1;length=61)); f(t)=sin(3t)+t^2/2
        rc=P.RegularPotential(PotentialCache(lg,f.(lg),zeros(61),zeros(61)))
        @test maximum(abs(P.potential_at(rc,10^t)[1]-f(t)) for t in lg) <= 1e-13
        # derivatives come from the same spline (not-a-knot slopes are O(h^3), h = 1/30: 1.3e-4 measured)
        @test maximum(abs(rc.spline(t,1;deriv_order=1)-3cos(3t)-t) for t in lg) <= 1e-3
    end
    d=NewtonianDisc(h0=.1,p=0.,q=0.,beta=0.,gamma=1.,soft=.03,lap_nphi=128)
    mode=F.NewtonianMode(2,1,1.); x=loggrid(d,.2,4.,401); c=PotentialCache(d,2,.2,4.;n=401)
    @testset "L4 the regular path forwards rtol, atol and dense; reconstruction rejects beta != 0" begin
        s(;kw...)=solve_lm_mode(d,mode,x,c;lindblad_regular=true,kw...)
        a=s(rtol=1e-10,atol=1e-12,dense=false)
        @test isequal(a,P._solve_m_mode(d,mode,x,P.RegularPotential(c);rtol=1e-10,atol=1e-12,dense=false))
        @test !isequal(a,s(dense=false))
        @test !isequal(a,s(rtol=1e-10,atol=1e-12,dense=true))
        @test_throws ArgumentError vars_from_enthalpy(NewtonianDisc(beta=.1),a,x,mode,c;lindblad_regular=true)
    end
    @testset "L5 named FR22 settings" begin
        S=P.FR22_SOLVER_SETTINGS
        @test S == (; reint_steps=2, phgrad_steps=1, phgrad_safeguard=true, anchor=:corotation,
                    rtol=1e-10, atol=1e-12, dense=false, lindblad_regular=true)
        @test isequal(solve_lm_mode(d,mode,x,c;S...),
                      solve_lm_mode(d,mode,x,c;reint_steps=2,phgrad_steps=1,phgrad_safeguard=true,
                                    anchor=:corotation,rtol=1e-10,atol=1e-12,dense=false,lindblad_regular=true))
    end
end

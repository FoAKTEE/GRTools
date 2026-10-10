using Test, Spacetime, NumericsCore, LinearAlgebra
using .EquatorialQuadratureOracle

@testset "C-ORBIT and C-HARM: coordinate-time eccentric projection" begin
    @testset "Exact rational first-integral and secular-phase identities" begin
        for a in (-3//5,0//1,9//10), r in (7//1,11//1), E in (9//10,19//20), L in (3//1,4//1)
            a,r,E,L=big.((a,r,E,L))
            Δ=r^2-2r+a^2; P=E*(r^2+a^2)-a*L
            R=P^2-Δ*(r^2+(L-a*E)^2)
            tλ=(r^2+a^2)*P/Δ+a*(L-a*E); φλ=a*P/Δ+L-a*E
            gtt=-(1-2/r); gtφ=-2a/r; gφφ=r^2+a^2+2a^2/r
            @test gtt*tλ^2+2gtφ*tλ*φλ+gφφ*φλ^2+r^2/Δ*R == -r^4
            @test -(gtt*tλ+gtφ*φλ)/r^2 == E
            @test (gtφ*tλ+gφφ*φλ)/r^2 == L
        end
        for m in -2:2, k in -2:2
            Γ,Υr,Υφ=113//1,7//3,11//3
            ω=m*Υφ/Γ+k*Υr/Γ
            @test ω*Γ/Υr-m*Υφ/Υr == k
        end
    end
    @test_throws ArgumentError EquatorialOrbit(1.,10.,0.2)
    @test_throws ArgumentError EquatorialOrbit(0.5,2.,0.2)
    @test_throws ArgumentError EquatorialOrbit(0.5,10.,1.)
    # Retrograde guard: the a=-0.9, e=0.3 separatrix is 9.5536 (reflected chart), not 2.605.
    @test_throws ArgumentError EquatorialOrbit(-0.9,9.0,0.3)
    @test EquatorialOrbit(-0.9,9.6,0.3).L > 0
    @test EquatorialOrbit(0,10,0).p isa Float64
    for (T,bits) in ((Float64,53),(BigFloat,128),(BigFloat,256))
        setprecision(BigFloat,bits) do
            tol = T===Float64 ? T(2e-10) : T(10)^(-bits÷4+4)
            # Compare equal callback counts: N full-period nodes versus N/2 paired
            # nodes. The 256-bit budget requires a refined full-period rule.
            N = bits == 53 ? 80 : bits == 128 ? 128 : 320
            paired_N = N ÷ 2
            normerr=zero(T); uniterr=zero(T); pairerr=zero(T); directerr=zero(T)
            for as in ("0.9","0","-0.7")
                a,p,e = parse.(T,(as,"10","0.3"))
                o = EquatorialOrbit(a,p,e)
                oracle=quadrature_orbit(a,p,e)
                @test isapprox(o.Γ,oracle.Γ;rtol=tol)
                @test isapprox(o.Υφ,oracle.Υφ;rtol=tol)
                for q in (T(pi)/3,4T(pi)/3,zero(T),T(pi))
                    s=worldline_rates(o,q)
                    @test s.r isa T
                    @test s.dt_dlambda > 0
                    @test s.dtau_dlambda == s.r^2
                    g = kerr_metric(a,s.r,T(pi)/2)
                    u = T[s.ut,s.ur,0,s.uphi]
                    residual=abs(dot(u,g*u)+1)
                    normerr=max(normerr,residual)
                    @test residual < tol
                    if 0<q<T(pi)
                        @test s.Ur > 0
                    elseif T(pi)<q<2T(pi)
                        @test s.Ur < 0
                    end
                end
                # Unit time average and true zero harmonics of a constant time signal.
                for k in (0,1,3)
                    v=orbit_harmonic((q,s)->one(T),o,0,k;N=N)
                    err=abs(v-(k==0 ? 1 : 0))
                    uniterr=max(uniterr,err)
                    @test err < tol
                    @test v isa Complex{T}
                end
                # A complex field with an odd radial-velocity term defeats unjustified
                # cosine-only / real-only pairing. Both branches must actually be sampled.
                H(q,s)=(1+im*one(T)/3)*s.r/o.p+(one(T)/5+im)*s.Ur/o.p
                for (m,k) in ((2,-1),(0,2),(-2,1),(1,-1))
                    full=orbit_harmonic(H,o,m,k;N=N)
                    paired=orbit_harmonic(H,o,m,k;N=paired_N,paired=true)
                    pairerr=max(pairerr,abs(full-paired))
                    @test abs(full-paired) < tol
                    @test orbit_frequency(o,m,k)==m*o.Ωφ+k*o.Ωr
                    # Reality-pair convention, without assuming H is real.
                    conjugate=orbit_harmonic((q,s)->conj(H(q,s)),o,-m,-k;N=paired_N,paired=true)
                    @test abs(conjugate-conj(paired)) < tol
                end
                # Direct coordinate-time integral via independent Darwin-angle quadrature:
                # ∫ H exp(i(ωt-mφ)) dt/T_r. No Mino-phase deviations enter this oracle.
                m,k=2,-1
                ω=(m*oracle.Υφ+k*oracle.Υr)/oracle.Γ
                Tr=2T(pi)*oracle.Γ/oracle.Υr
                function direct(χ)
                    r=oracle.radius(χ)
                    dl=oracle.dlambda(χ)
                    Ur=o.p*o.e*sin(χ)/(1+o.e*cos(χ))^2/dl
                    t=integrate_oracle(z->oracle.rates(z).t*oracle.dlambda(z),zero(T),χ,oracle.nodes,oracle.weights)
                    φ=integrate_oracle(z->oracle.rates(z).phi*oracle.dlambda(z),zero(T),χ,oracle.nodes,oracle.weights)
                    return ((1+im*one(T)/3)*r/o.p+(one(T)/5+im)*Ur/o.p)*cis(ω*t-m*φ)*oracle.rates(χ).t*dl
                end
                # Integrate both halves separately, with an independent node grid.
                reference=(integrate_oracle(direct,zero(T),T(pi),oracle.nodes,oracle.weights)+
                           integrate_oracle(direct,T(pi),2T(pi),oracle.nodes,oracle.weights))/Tr
                projected=orbit_harmonic(H,o,m,k;N=paired_N,paired=true)
                directerr=max(directerr,abs(projected-reference))
                @test abs(projected-reference) < tol
                println("HARM orbit bits=$bits a=$a done"); flush(stdout)
            end
            # Circular checks stay on the eccentric projector's path.
            circ=EquatorialOrbit(T(1)/2,T(10),zero(T))
            @test abs(orbit_harmonic((q,s)->one(T),circ,2,0;N=N)-1) < tol
            @test abs(orbit_harmonic((q,s)->one(T),circ,2,1;N=N)) < tol
            @test_throws ArgumentError orbit_harmonic((q,s)->1,circ,0,0;N=0)
            println("HARM bits=$bits normalization=$normerr unit=$uniterr paired=$pairerr direct_time=$directerr"); flush(stdout)
        end
    end
    # Newtonian limit is a sequence, not equality at a finite radius.
    errors=[abs(EquatorialOrbit(0.,p,0.2).Ωr/EquatorialOrbit(0.,p,0.2).Ωφ-1) for p in (1e3,1e4,1e5)]
    @test errors[3]<errors[2]<errors[1]
    @test errors[3]<4e-5
end

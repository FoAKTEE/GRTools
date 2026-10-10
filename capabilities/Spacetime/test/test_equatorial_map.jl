using Test, Spacetime, LinearAlgebra, NumericsCore

@testset "C-MAP constrained inverse and bounded conditioning" begin
    for T in (Float64,BigFloat)
        setprecision(BigFloat,256) do
            worst=zero(T); condition=zero(T)
            for a in T.((-0.9,-0.5,-0.1,0,0.1,0.5,0.9)), p in T.((12,20,50,250)),
                e in T.((0.05,0.1,0.5)), x in T.((-1,1))
                J=EL_to_pe_jacobian(a,p,e;x)
                E,L,_=constants_of_motion(a,p,e,x)
                # Differentiate at higher precision so cancellation in the
                # constants evaluator does not dominate the derivative oracle.
                F=setprecision(BigFloat,512) do
                    ab,pb,eb,xb=BigFloat.((a,p,e,x))
                    h=big"1e-40"
                    C(p,e)=collect(Tuple(constants_of_motion(ab,p,e,xb))[1:2])
                    T.(hcat((C(pb+pb*h,eb)-C(pb-pb*h,eb))/(2pb*h),
                            (C(pb,eb+h)-C(pb,eb-h))/(2h)))
                end
                err=opnorm(J*F-I,Inf)
                worst=max(worst,err)
                @test err < (T===Float64 ? T(1e-8) : big"1e-60")
                chain=T[1 0;0 1;-2a*(L-a*E) 2*(L-a*E)]
                full=orbital_params_jacobian(a,p,e,x)*chain
                @test maximum(abs,full[3,:]) < 100eps(T)
                @test J == full[1:2,:]
                κ=opnorm(J,Inf)*opnorm(inv(J),Inf)
                condition=max(condition,κ)
                @test κ < T(1e7)
            end
            println("C-MAP $T inverse residual=$worst maximum condition_inf=$condition; grid: 168 signed-spin/orientation cases")
        end
    end
    @test_throws ArgumentError EL_to_pe_jacobian(-0.9,9.0,0.3)
    @test_throws ArgumentError EL_to_pe_jacobian(0.9,10.0,0.0)
    @test_throws ArgumentError EL_to_pe_jacobian(0.9,10.0,0.1;x=0.5)
    @test_throws ArgumentError EL_to_pe_jacobian(0.9,10.0,0.1;max_condition=1)
    @test_throws ArgumentError EL_to_pe_jacobian(-0.9,separatrix(-0.9,0.3,1)+1e-10,0.3)
end
println("PASS C-MAP constrained K chain, xdot=0, inverse derivatives, bounded sampled conditioning and stop guards")

@testset "C-MAP reaction rates and modal rewrite" begin
    for a in (-0.9,0.0,0.9)
        o=EquatorialOrbit(a,20.0,0.2)
        J=EL_to_pe_jacobian(o.a,o.p,o.e)
        alpha=-J[2,2]/J[2,1]
        modes=[(m=1,k=-2,torque=2e-15),(m=2,k=1,torque=-1e-15),(m=3,k=0,torque=4e-15)]
        power0=3e-17; mu=0.1
        Ldot=-sum(v.torque for v in modes)/mu
        Edot=-(sum(orbit_frequency(o,v.m,v.k)*v.torque/v.m for v in modes)+power0)/mu
        numerator=eccentric_mode_numerator(o,alpha,modes;axisymmetric_power=power0,mu)
        r=ELdot_to_pedot(o,Edot,Ldot;eccentric_numerator=numerator)
        @test [r.pdot,r.edot] ≈ J*[Edot,Ldot] rtol=2e-13
        @test r.e2dot == 2o.e*r.edot
        @test numerator ≈ Edot-alpha*Ldot rtol=2e-13
        @test eccentric_mode_numerator(o,alpha,();axisymmetric_power=power0,mu) == -power0/mu
        @test_throws ArgumentError ELdot_to_pedot(o,Edot,Ldot;max_fraction_per_radial_period=0)
        @test_throws ArgumentError eccentric_mode_numerator(o,alpha,[(m=0,k=1,torque=1.0)])
    end
    # A resolved sideband survives the rewrite after direct Float64 subtraction
    # has rounded it away. Compare with the same finite inputs at high precision.
    ph=1/32; rr=1/64
    modes=[(m=1,k=0,torque=1.0),(m=1,k=1,torque=1e-18)]
    stable=eccentric_mode_numerator(ph,rr,ph,modes;axisymmetric_power=1e-20)
    direct=-sum((ph+v.k*rr)*v.torque for v in modes)-1e-20+ph*sum(v.torque for v in modes)
    oracle=-(BigFloat(rr)*BigFloat(1e-18)+BigFloat(1e-20))
    @test iszero(direct)
    @test stable ≈ oracle rtol=1e-15
    setprecision(BigFloat,256) do
        o=EquatorialOrbit(big"-0.9",big"20",big"0.2")
        J=EL_to_pe_jacobian(o.a,o.p,o.e)
        Edot,Ldot=big"-1e-30",big"-2e-29"
        r=ELdot_to_pedot(o,Edot,Ldot)
        @test all(v -> v isa BigFloat && precision(v)==256,Tuple(r))
        @test [r.pdot,r.edot] ≈ J*[Edot,Ldot] rtol=big"1e-65"
        @test_throws ArgumentError ELdot_to_pedot(EquatorialOrbit(o.a,o.p,zero(o.e)),Edot,Ldot)
    end
end
println("PASS C-MAP mode-sum reaction sign, signed sidebands, m=0 power, and rate conversion")

# These grids are deterministic polynomial-identity certificates, not random
# samples. Degree bounds and denominator exclusions are in CMAP_CERTIFICATE.md.
@testset "C-MAP exact rational coefficient certificates" begin
    S=Rational{BigInt}.((-1,1))
    for a in Rational{BigInt}.((-1,0,1)), E in S, L in S, dE in S, dL in S
        b=jadd([L,dL,0],jscale(-a,[E,dE,0]))
        @test jmul(b,b)[2] == 2(L-a*E)*(dL-a*dE)
        polar=[2a*(L-a*E) -2*(L-a*E) 1]
        @test iszero(polar*[1 0;0 1;-2a*(L-a*E) 2*(L-a*E)])
    end
    # Multiplying the residual by m*mu yields a polynomial of degree at most
    # one in each of (m,k,Omega_phi,Omega_r,alpha,torque,power0).
    for m in (1,2), k in (-1,1), ph in S, rr in S, alpha in S, torque in S, power0 in S
        mu=big(3)//7
        direct=-(((m*ph+k*rr)/m)*torque+power0)/mu+alpha*torque/mu
        stable=eccentric_mode_numerator(ph,rr,alpha,[(;m,k,torque)];axisymmetric_power=power0,mu)
        @test m*mu*(direct-stable) == 0
    end
    # Finite-sum linearity extends the coefficient identity to every finite mode
    # list, including m=0 power. No sampled probability bound is needed.
end
println("PASS C-MAP exact_rational K/polar and mode-sum coefficient identities (48 and 128 unisolvent grid points)")

include("test_equatorial_map_circular.jl")

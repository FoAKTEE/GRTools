@testset "C-MAP exact double-root expansion" begin
    # A tensor grid larger than each cleared polynomial's degree: a<=4,p<=4,
    # E<=2,L<=2,Es<=1,Ls<=1. All 900 nodes, no rejection or random sampling.
    for a in Rational{BigInt}.(-2:2), p in Rational{BigInt}.(1:5),
        E in Rational{BigInt}.((-1,0,1)), L in Rational{BigInt}.((-1,0,1)),
        Es in Rational{BigInt}.((-1,1)), Ls in Rational{BigInt}.((-1,1))
        r=jscale(p,jinv(Rational{BigInt}[1,1,0,0]))
        r2=jmul(r,r); ej=[E,0,2Es,0]; lj=[L,0,2Ls,0]
        P=jadd(jmul(ej,jadd(r2,constjet(a^2,3,Rational{BigInt}))),jscale(-a,lj))
        delta=jadd(jadd(r2,jscale(-2,r)),constjet(a^2,3,Rational{BigInt}))
        b=jadd(lj,jscale(-a,ej))
        Rj=jadd(jmul(P,P),jscale(-1,jmul(delta,jadd(r2,jmul(b,b)))))
        b0=L-a*E
        Rr=4(E^2-1)*p^3+6p^2+2(a^2*(E^2-1)-L^2)*p+2b0^2
        Rrr=12(E^2-1)*p^2+12p+2(a^2*(E^2-1)-L^2)
        Rrrr=24(E^2-1)*p+12
        RCs=(2E*p^4+2a^2*E*p^2-4a*b0*p)*Es+(-2L*p^2+4b0*p)*Ls
        RrCs=(8E*p^3+4a^2*E*p-4a*b0)*Es+(-4L*p+4b0)*Ls
        @test Rj[3] == p^2*Rrr+2p*Rr+2RCs
        @test Rj[4] == -p^3*Rrrr-6p^2*Rrr-6p*RrCs-6p*Rr
    end
end
println("PASS C-MAP exact_rational double-root expansion (900 unisolvent grid points, 1800 zero residuals)")

# Independent circular continuation from the DOUBLE radial root, not the
# eccentric turning-point Jacobian. s=e², C=(E,L), R=R_r=0 at r=p.
# Expanding R(p/(1±e),C(p,s)) through e³ gives
# R_C C_s=-p² R_rr/2 and R_rC C_s=-p R_rr-p² R_rrr/6.
function circular_forward_map(a,p)
    E,L,_=constants_of_motion(a,p,zero(p),one(p))
    b=L-a*E
    RC=[2E*p^4+2a^2*E*p^2-4a*b*p, -2L*p^2+4b*p]
    RrC=[8E*p^3+4a^2*E*p-4a*b, -4L*p+4b]
    Rrr=12(E^2-1)*p^2+12p+2(a^2*(E^2-1)-L^2)
    Rrrr=24(E^2-1)*p+12
    A=vcat(transpose(RC),transpose(RrC))
    return A \ [zero(p) -p^2*Rrr/2; -Rrr -p*Rrr-p^2*Rrrr/6]
end

@testset "C-MAP separate removable d(e²)/dt limit" begin
    setprecision(BigFloat,256) do
        for a in (big"-0.9",big"0",big"0.9"), p in (big"12",big"20",big"50")
            F0=circular_forward_map(a,p); J0=inv(F0)
            omega=inv(p*sqrt(p)+a)
            @test F0[1,1]/F0[2,1] ≈ omega rtol=big"1e-65"
            Edot,Ldot=big"2e-14",-big"3e-13"
            limit=(J0*[Edot,Ldot])[2]
            last_error=big"Inf"
            for e in (big"1e-4",big"1e-6",big"1e-8")
                J=EL_to_pe_jacobian(a,p,e)
                Js=Diagonal([one(e),2e])*J
                err=opnorm(Js-J0,Inf)/opnorm(J0,Inf)
                @test err < 100e^2
                @test err < last_error
                last_error=err
                sdot=(Js*[Edot,Ldot])[2]
                @test abs(sdot-limit) < 100e^2*abs(limit)
            end
            # For arbitrary forcing sdot has a finite, generally nonzero limit.
            # This alone DOES NOT imply finite de/e.
            @test !iszero(limit)
            println("e2-limit a=$a p=$p relative_map_error=$last_error")
        end
    end
end
println("PASS C-MAP e=0 removable d(e^2)/dt: double-root continuation, separately checked")

@testset "C-MAP separate removable de/e limit for a circular-preserving rate law" begin
    setprecision(BigFloat,256) do
        for a in (big"-0.9",big"0",big"0.9"), p in (big"12",big"20",big"50")
            F0=circular_forward_map(a,p)
            alpha0=F0[1,1]/F0[2,1]
            h=big"1e-15"
            # Mixed partials from the independent double-root C_s(p,0), not
            # a differentiation of the eccentric code under test.
            Cps=(circular_forward_map(a,p+h)[:,2]-circular_forward_map(a,p-h)[:,2])/(2h)
            alpha_s=(Cps[1]-alpha0*Cps[2])/F0[2,1]
            B0=F0[1,2]-alpha0*F0[2,2]
            Ldot=-big"1e-12"; beta=big"7e-15"
            # Smooth manufactured forcing: E_dot=Omega_c L_dot+beta*s.
            # N(0)=0, N_s(0)=beta-alpha_s*L_dot, so gamma=N_s/(2 B0).
            gamma=(beta-alpha_s*Ldot)/(2B0)
            last_error=big"Inf"
            for e in (big"1e-4",big"1e-6",big"1e-8")
                J=EL_to_pe_jacobian(a,p,e)
                rates=J*[alpha0*Ldot+beta*e^2,Ldot]
                ratio=rates[2]/e
                err=abs((ratio-gamma)/gamma)
                @test err < 1000e^2
                @test err < last_error
                last_error=err
                @test abs(2e*rates[2]) < 4abs(gamma)*e^2
            end
            @test !iszero(gamma)
            println("de/e-limit a=$a p=$p gamma=$gamma relative_error=$last_error")
        end
    end
end
println("PASS C-MAP e=0 removable de/e: N(0)=0 and differentiable rate law required; limit N_s(0)/(2 B(0))")

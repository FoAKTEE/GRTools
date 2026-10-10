using Test, Spacetime, LinearAlgebra
using Spacetime.SpinningSecondary: spin_orbital_params_jacobian

@testset "Geodesics turning-point owner and constrained C-MAP chain" begin
    @test parentmodule(turning_point_jacobian) === Spacetime
    @test parentmodule(orbital_params_jacobian) === Spacetime
    for (T,bits) in ((Float64,53),(BigFloat,256))
        setprecision(BigFloat,bits) do
            h = eps(T)^(one(T)/3)
            worst=zero(T)
            for as in ("-0.7","0","0.9"), ps in ("10","20"), es in ("0.1","0.4"), xs in ("-1","0.6","1")
                a,p,e,x = parse.(T,(as,ps,es,xs))
                J=orbital_params_jacobian(a,p,e,x)
                @test isequal(J,spin_orbital_params_jacobian(a,p,e,x))
                # The constants owner normalizes equatorial signed spin.
                if abs(x)==1
                    E,L,_ = constants_of_motion(a,p,e,x)
                    chain = T[1 0; 0 1; -2a*(L-a*E) 2*(L-a*E)]
                    constrained = J*chain
                    # Independent forward derivatives of E,L wrt p,e, at fixed x.
                    C(p,e)=collect(Tuple(constants_of_motion(a,p,e,x))[1:2])
                    hp=p*h
                    forward=hcat((C(p+hp,e)-C(p-hp,e))/(2hp),(C(p,e+h)-C(p,e-h))/(2h))
                    err=maximum(abs,constrained[1:2,:]*forward-I)
                    worst=max(worst,err)
                    @test err < 5000eps(T)^(T(2)/3)
                    @test maximum(abs,constrained[3,:]) < 100eps(T)
                end
            end
            println("C-MAP bits=$bits inverse_residual=$worst")
        end
    end
end

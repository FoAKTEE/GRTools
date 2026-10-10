using Test, Spacetime
using .EquatorialQuadratureOracle

@testset "Equatorial full-period branch and independent equations" begin
    for (T, bits) in ((Float64, 53), (BigFloat, 128), (BigFloat, 256))
        setprecision(BigFloat, bits) do
            # All decimal inputs and phases originate in the working precision.
            tol = T === Float64 ? T(2e-10) : T(10)^(-bits÷4+4)
            h = eps(T)^(one(T)/3)
            dtol = 500eps(T)^(T(2)/3)
            max_oracle = zero(T); max_derivative = zero(T); max_radial = zero(T)
            for (as, ps, es) in (("0.9","10","0.3"), ("0","12","0.3"), ("-0.7","11","0.25"))
                a,p,e = parse.(T,(as,ps,es))
                f = geodesic_trajectory_equatorial(a,p,e)
                radial = geodesic_trajectory_radial(a,p,e,one(T))
                oracle = quadrature_orbit(a,p,e)
                freq = kerr_frequencies(a,p,e,one(T); Time="Mino")
                @test isapprox(freq.Ur,oracle.Υr; rtol=tol)
                @test isapprox(freq.Ut,oracle.Γ; rtol=tol)
                for q in T.((0,1,2,3,4,5,6)) .* (2T(pi)/7)
                    g = f(q); ref = oracle.at_phase(q)
                    for key in keys(g)
                        err = abs(getproperty(g,key)-getproperty(ref,key))/(1+abs(getproperty(ref,key)))
                        max_oracle = max(max_oracle,err)
                        @test err < tol
                    end
                    for key in keys(radial(q))
                        err = abs(getproperty(g,key)-getproperty(radial(q),key))/(1+abs(getproperty(g,key)))
                        max_radial = max(max_radial,err)
                        @test err < tol
                    end
                    # The differential equations, independently evaluated from E,L and the metric.
                    E,L,_ = constants_of_motion(abs(a),p,e,a < 0 ? -one(T) : one(T))
                    a < 0 && (L = -L)
                    r = g.r; P = E*(r^2+a^2)-a*L; Δ = r^2-2r+a^2
                    expected = (g.Ur, (r^2+a^2)*P/Δ+a*(L-a*E)-freq.Ut,
                                a*P/Δ+L-a*E-freq.Uph, r^2-oracle.Υτ)
                    for (key, value) in zip((:r,:Δtr,:Δφr,:Δτr),expected)
                        d = (getproperty(f(q+h),key)-getproperty(f(q-h),key))*freq.Ur/(2h)
                        err = abs(d-value)/(1+abs(value))
                        max_derivative=max(max_derivative,err)
                        @test err < dtol
                    end
                    q > T(pi) && @test g.Ur < 0
                end
                for q in (zero(T),T(pi),2T(pi))
                    @test abs(f(q).Ur) < tol
                    @test maximum(abs,(f(q).Δtr,f(q).Δφr,f(q).Δτr)) < tol
                end
            end
            println("BRANCH bits=$bits oracle_scaled=$max_oracle radial_scaled=$max_radial derivative_scaled=$max_derivative")
        end
    end
end

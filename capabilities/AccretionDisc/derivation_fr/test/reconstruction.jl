using Random

@testset "Bounded Num runtime reconstruction" begin
    s = FR.fr22_master_residuals()
    v = s.symbols
    for T in (Float64,BigFloat)
        setprecision(192) do
            rng = MersenneTwister(221921)
            worst = zero(T)
            for j in 1:16
                # Rationally specified inputs, numerically evaluated square
                # roots/quadrature. No symbolic quadrature is expanded here.
                r = Num(T(j <= 8 ? 3//5 : 31//10))
                hp, p = Num(T(1//10)), Num(T(rand(rng,0:6)//4))
                q = Num(T(j <= 4 ? 0//1 : rand(rng,1:6)//4))
                d = FR.certificate_disc(hp,p,q)
                m = rand(rng,2:7)
                h = Num(T(rand(rng,-17:17)//19))+im*Num(T(rand(rng,-17:17)//19))
                h1 = Num(T(rand(rng,-17:17)//19))+im*Num(T(rand(rng,-17:17)//19))
                actual = FR.Modes.vars_from_enthalpy(d,reshape([h,h1],1,2),[r],m)
                point = Dict{Num,Any}(v.r=>r,v.S=>FR.Modes.Sigma(d,r),v.O=>FR.Modes.Omega(d,r),
                    v.K=>FR.Modes.kappa(d,r)^2,v.D=>FR.Modes.Dfreq(d,r,m),
                    v.w=>FR.Modes.omega_pat(d,m)-m*FR.Modes.Omega(d,r),
                    v.c2=>FR.Modes.csiso(d,r)^2,v.L=>FR.Modes.LTinv(d,r),
                    v.h=>h,v.h1=>h1,v.P=>FR.Modes.Phi_m(d,r,m),
                    v.P1=>FR.Modes.dPhidr_m(d,r,m),v.m=>m)
                expected = Symbolics.substitute.((s.ur,s.up,s.ds),Ref(point))
                for k in 1:3
                    @test actual[k][1] isa Complex{Num}
                    a = Complex{T}(Symbolics.symbolic_to_float(real(actual[k][1]))+im*Symbolics.symbolic_to_float(imag(actual[k][1])))
                    b = Complex{T}(Symbolics.symbolic_to_float(real(expected[k]))+im*Symbolics.symbolic_to_float(imag(expected[k])))
                    err = abs(a-b)/max(one(T),abs(a),abs(b))
                    worst = max(worst,err)
                    @test err <= 512eps(T)
                end
            end
            println("[PASS] runtime Num reconstruction ",T," points=16 components=3 max_scaled_error=",worst,
                    " tolerance=",512eps(T)," seed=221921; numerical, no expanded quadrature")
        end
    end
end

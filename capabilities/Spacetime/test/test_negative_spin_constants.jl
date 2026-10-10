using Test, Spacetime
include("fixtures/negative_spin_unchanged_bits.jl")

@testset "signed-spin equatorial constants and separatrix" begin
    for (args,bits) in NEGATIVE_SPIN_UNCHANGED_BITS
        @test reinterpret.(UInt64,Tuple(constants_of_motion(args...))) == bits
    end
    setprecision(BigFloat, 256) do
        worst = zero(BigFloat)
        for a in parse.(BigFloat, ("-0.9", "-0.5", "-0.1")),
            e in parse.(BigFloat, ("0.1", "0.5")), x in (-big"1", big"1")
            p = big"12"
            E,L,Q = constants_of_motion(a,p,e,x)
            @test isequal(E, geodesic_energy(a,p,e,x))
            @test isequal(L, geodesic_angular_momentum(a,p,e,x))
            @test iszero(Q)
            for r in (p/(1+e),p/(1-e))
                R = (E*(r^2+a^2)-a*L)^2-(r^2-2r+a^2)*(r^2+(L-a*E)^2+Q)
                worst = max(worst,abs(R))
                @test abs(R) < big"1e-65"
            end
            ps = separatrix(a,e,x)
            Es,Ls,_ = constants_of_motion(a,ps,e,x)
            rs = ps/(1+e)
            dR = 4rs*Es*(Es*(rs^2+a^2)-a*Ls) - (2rs-2)*(rs^2+(Ls-a*Es)^2) - 2rs*(rs^2-2rs+a^2)
            @test abs(dR) < big"1e-60"
            @test x > 0 ? ps > 6+2e : ps < 6+2e
        end
        println("negative-spin turning-point max absolute R residual = ", worst)
        for e in (big"0.1",big"0.5"), x in (-big"1",big"1")
            zero_c = constants_of_motion(big"0",big"12",e,x)
            for a in (-big"1e-20",big"1e-20")
                c = constants_of_motion(a,big"12",e,x)
                @test maximum(abs.(Tuple(c).-Tuple(zero_c))) < big"1e-18"
                @test abs(separatrix(a,e,x)-(6+2e)) < big"1e-18"
            end
        end
    end
end
println("PASS signed-spin turning points, double-root separatrix, and two-sided a=0 continuity")

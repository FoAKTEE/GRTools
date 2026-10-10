@testset "C-NORM phase, units and pair counting" begin
    EC = R3.EccentricConventions
    certs = EC.convention_certificates()
    @test length(certs) == 20
    @test all(c.pass for c in certs)
    for pair in ((0,0),(0,-1),(-1,2))
        @test_throws DomainError EC.reality_pair_weight(pair...)
    end
    @test EC.reality_pair_weight(0,1) isa Float64
    for bits in (128,256)
        setprecision(BigFloat,bits) do
            for T in (Float64,BigFloat)
                mode = EC.ForcingModes.NewtonianMode(2,-1,T(7)/3)
                @test mode.k == -3
                @test mode.omega == -T(7)/3
                @test mode.omega isa T
                A, B = complex(T(2)/3,T(4)/7), complex(-T(5)/4,T(3)/8)
                # Independent direct quadrature includes the m=0,k>0 case.
                for (m,k) in ((2,-3),(0,1),(0,3))
                    N = 32
                    integral = zero(T)
                    for j in 0:N-1, h in 0:N-1
                        phase = 2T(pi)*(m*j-k*h)/N
                        integral += 2real(A*cis(phase))*2real(B*cis(phase))
                    end
                    integral *= 2T(pi)/N^2
                    expected = EC.reality_pair_weight(m,k,T)*real(A*conj(B))
                    @test integral ≈ expected rtol=200eps(T)
                    @test EC.reality_pair_weight(m,k,T) isa T
                end
            end
        end
    end
end

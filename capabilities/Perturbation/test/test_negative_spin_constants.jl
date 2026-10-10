using Test, Perturbation

@testset "negative-spin eccentric geodesic consumer" begin
    # Turning-point Newton solve from Schwarzschild seeds, followed by 256-bit
    # Fourier-Darwin quadrature at 64/128 nodes, gives this independent orbit
    # input to the established Teukolsky convolution. See st-negative-spin-constants
    # evidence/consumer_oracle.jl; no constants/frequency/trajectory call in the oracle.
    reference=0.00081529602028259198
    for N in (64,128)
        md=TeukolskyPointParticleModeEccentric(2,2,1,-0.9,10.0,0.1;N)
        @test abs(md.Zinf) ≈ reference rtol=2e-9
        @test md.En ≈ 0.96239431059444835 rtol=2e-13
        @test md.Lz ≈ 4.2047945557533346 rtol=2e-13
    end
end
println("PASS Perturbation negative-spin consumer against independent geodesic quadrature")

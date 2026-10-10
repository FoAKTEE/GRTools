# Modes/Forcing.jl consumer of the signed-spin eccentric equatorial constants (Spacetime
# constants.jl reflected chart): the a=-0.9, r_p=250, m=1 ILR torque against an independent
# geodesic oracle (turning-point Newton solve + Fourier-Darwin quadrature at e=.01,.02 and the
# analytic d(Omega-kappa)/dr): evidence/st-negative-spin-constants/consumer_torque_oracle.jl.
# The runtime divides by the central difference _kDprime(h=1e-6), whose 8e-9 relative error
# accounts for the observed 7.9e-9 gap; rtol=2e-7 is far below the wrong-chart shifts (>=1.6e-3).
using Test, AccretionDisc
using AccretionDisc.Modes: RelativisticDisc, lindblad_torque
@testset "negative-spin Lindblad torque consumer" begin
    value=lindblad_torque(RelativisticDisc(a=-0.9),250.0,1;kind=:ILR)
    @test value ≈ -1.2622185469017366e-5 rtol=2e-7
end
println("PASS AccretionDisc negative-spin torque against independent geodesic quadrature")

using Test, AccretionDisc

@testset "FR24 Keplerian orbital rates" begin
    M = AccretionDisc.Modes
    F = AccretionDisc.ForcingModes
    # Manufactured modes use the landed Jacobi-power owner, not a second torque calculation.
    modes = [F.NewtonianMode(2,3,0.7), F.NewtonianMode(3,2,0.7), F.NewtonianMode(4,4,0.7)]
    torques = [0.2,-0.07,0.03]
    diagnostics = [(net=t, power=F.mode_power(m,t), FJ0=0.4) for (m,t) in zip(modes,torques)]
    rates = M.fr24_orbital_rates(diagnostics; Mp=0.02,a=1.3,Omega_p=0.7,e=0.2)
    @test rates.torque ≈ 0.16
    @test rates.power ≈ 0.19833333333333333
    @test rates.L ≈ 0.02*sqrt((0.7^2*1.3^3)*1.3*(1-0.2^2))
    @test rates.tau0_inv ≈ 0.4/(0.02*0.7*1.3^2)
    @test rates.Ldot == -rates.torque
    @test rates.Edot == -rates.power
    @test rates.tauL_inv ≈ -rates.Ldot/rates.L
    @test rates.taua_inv ≈ rates.Edot/rates.E
    @test rates.tauL_inv ≈ rates.taua_inv/2-0.2^2/(1-0.2^2)*rates.taue_inv
    @test rates.Ldot/rates.L ≈ rates.adot/(2*1.3)-0.2*rates.edot/(1-0.2^2)
    @test rates.adot == -1.3*rates.taua_inv
    @test rates.edot == -0.2*rates.taue_inv
    for (m,d) in zip(modes,diagnostics)
        one_mode = M.fr24_orbital_rates([d]; e=0.2)
        @test one_mode.Edot-(m.omega/m.m)*one_mode.Ldot ≈ 0 atol=1e-16
    end
    # Exact circular convention, including outward migration and a zero sum.
    for t in (-0.3,0.0,0.3)
        circular = M.fr24_orbital_rates([(net=t,power=0.7t,FJ0=0.4)]; Omega_p=0.7)
        @test circular.taue_inv === 0.0
        @test iszero(circular.edot)
        @test circular.taua_inv ≈ 2circular.tauL_inv
    end
    @test_throws ArgumentError M.fr24_orbital_rates([(net=1.,power=2.,FJ0=1.)])
    @test_throws ArgumentError M.fr24_orbital_rates(diagnostics; e=1.)
    @test_throws ArgumentError M.fr24_orbital_rates(diagnostics; e=-0.1)
    @test_throws ArgumentError M.fr24_orbital_rates(diagnostics; Mp=0.)
    @test_throws ArgumentError M.fr24_orbital_rates([])
    @test_throws ArgumentError M.fr24_orbital_rates([diagnostics[1],merge(diagnostics[2],(;FJ0=0.5))];e=0.2)
    setprecision(192) do
        b = M.fr24_orbital_rates([(net=big"0.2",power=big"0.3",FJ0=big"0.4")];e=big"0.2")
        @test all(x -> x isa BigFloat, values(b))
        @test b.tauL_inv ≈ b.taua_inv/2-big"0.04"/big"0.96"*b.taue_inv rtol=big"1e-55"
        @test b.L ≈ sqrt(big"0.96") rtol=big"1e-55"
    end
    # Consume a real diagnostics NamedTuple, from a manufactured enthalpy response.
    d = M.NewtonianDisc(h0=0.06,p=0.,q=0.,beta=0.)
    mode = F.NewtonianMode(2,2,1.)
    r = collect(range(0.5,1.5,length=31))
    cache = M.PotentialCache(log10.(r),fill(0.1,length(r)),zeros(length(r)),zeros(length(r)))
    f = [complex(0.01x,0.002x^2) for x in r, _ in 1:2]
    diag = M.fr24_mode_diagnostics(d,mode,r,f,cache)
    actual = M.fr24_orbital_rates([diag];e=0.1)
    @test actual.torque == diag.net
    @test actual.power == diag.power
    @test actual.tau0_inv == diag.FJ0
end

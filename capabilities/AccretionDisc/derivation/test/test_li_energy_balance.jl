@testset "C-LI-HEAT energy, baryons and gauge" begin
    @test_throws ArgumentError R3.li_energy_balance_certificates(; npts=0, verbose=false)
    certs = R3.li_energy_balance_certificates(; npts=NPTS)
    @test length(certs) == 19
    for c in certs
        @test c.pass || error("$(c.name): $(c.note)")
    end
end

@testset "thermostat symbolic sensitivity completion" begin
    certs=R3.thermostat_certificates()
    @test length(certs)==15
    @test all(c->c.pass && c.modality==:symbolic,certs)
    ex=R3.thermostat_expressions()
    # Negative controls: deleting the compression heat or imposing primary
    # Gamma_eff on the density law produces a nonzero symbolic residual.
    @test !R3.certify_symbolic("omitted compression heat","negative control",
        ex.energy_rate+(ex.energy+ex.pressure)*ex.expansion-ex.rho*ex.cv*ex.advtemp;verbose=false).pass
    @test !R3.certify_symbolic("primary compression index alias","negative control",
        ex.gamma_eff-(1+ex.theta);verbose=false).pass
    @test !R3.certify_symbolic("zero-bath relation19 alias","negative control",
        R3.subs(ex.lagenergy,Dict(ex.temp_r=>0))/(ex.energy+ex.pressure)-ex.drho/ex.rho;verbose=false).pass
end

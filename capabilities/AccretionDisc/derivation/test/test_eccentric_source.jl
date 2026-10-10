@testset "C-SOURCE symbolic local model" begin
    certs = R3.source_certificates()
    @test !isempty(certs)
    @test all(c.pass for c in certs)
end

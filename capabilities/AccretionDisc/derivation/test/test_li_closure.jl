# Loaded by runtests.jl; the LI certificates are also independently runnable.
@testset "Option A locally isothermal closure" begin
    LI = R3.LocallyIsothermal
    @test isequal(R3.CLOSURE_RESIDUAL, R3.pressure_closure())
    @test isequal(LI.CLOSURE_RESIDUAL, LI.dp / LI.p - LI.deps / LI.eps)
    @test_throws MethodError R3.pressure_closure(Val(:unknown))
    mktempdir() do out
        certs = LI.closure_certificates(out; npts = NPTS)
        @test !isempty(certs)
        @test all(c.pass for c in certs)
    end
end

using Test
using Symbolics
using FRDerivation

const FR = FRDerivation
@testset "FR derivation environment" begin
    @test all(c -> c.pass, FR.certify_environment())
    @test FR.certificate_disc(0.1, 1.0, 0.0).h0 isa Float64
    setprecision(192) do
        d = FR.certificate_disc(big"0.1", big"1", big"0")
        @test FR.Modes.Omega(d, big"2") isa BigFloat
        @test precision(FR.Modes.Omega(d, big"2")) == 192
    end
    @variables x
    @test !FR.certify_symbolic("negative control", "test", x+1; verbose=false).pass
    @test_throws ErrorException FR.exact_eval(x, Dict(x => 0.5))
    @test FR.exact_eval(x, Dict(x => FR.RatB(1, 3))) == FR.RatB(1, 3)
end

@testset "FR22 Kepler" begin
    certs = FR.certify_fr22_kepler()
    @test length(certs) == 9
    @test all(c -> c.pass && c.modality == :symbolic, certs)
end

@testset "FR22 background" begin
    certs = FR.certify_fr22_background()
    @test length(certs) == 15
    @test all(c -> c.pass, certs)
    # Separate high-precision runtime check of the same identity, including
    # noninteger p,q and nonzero pressure support (not the p=q=0 subspace).
    setprecision(192) do
        for (p,q,r) in ((big"0.4",big"0.75",big"0.6"),
                        (big"1.5",big"0.2",big"3.1"))
            d = FR.certificate_disc(big"0.08",p,q)
            lhs = FR.Modes.kappa(d,r)^2
            rhs = 4FR.Modes.Omega(d,r)^2+2r*FR.Modes.Omega(d,r)*FR.Modes.dOmegadr(d,r)
            @test isapprox(lhs,rhs;rtol=64eps(BigFloat))
        end
    end
end

@testset "FR24 units and dimensional audit" begin
    certs = FR.certify_fr24_normalisation()
    @test all(c -> c.pass && c.modality == :symbolic, certs)
    @test any(c -> occursin("extra time",c.name), certs)
    @test any(c -> occursin("corrected by Omega_p",c.name), certs)
    # The omitted inverse-time factor is a failing dimensional identity.
    @variables mass length time
    @test !FR.certify_symbolic("uncorrected Fourier units", "test",
        length^2/time-length^2/time^2; verbose=false).pass
end

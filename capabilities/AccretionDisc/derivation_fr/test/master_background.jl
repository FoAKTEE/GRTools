@testset "FR22 arbitrary-q master and MR20 transport" begin
    certs = FR.certify_fr22_master()
    @test all(c -> c.pass,certs)
    @test count(c -> c.modality == :exact_rational,certs) == 1
    @test any(c -> c.eqref == "2207.14637v2::eq:thermal_length",certs)
    for label in ("eq:master_eqn","eq:dsigma_dh","eq:dur_dh","eq:duphi_dh")
        @test any(c -> c.eqref == "2207.14637v2::$label",certs)
        @test any(c -> c.eqref == "2407.20398v1::$label",certs)
    end
    s = FR.fr22_master_residuals()
    (;r,O,w,L,Lr,h,h1,m) = s.symbols
    # Detect dropped temperature derivative, wrong source sign, and a q=0-only
    # operator mistakenly used for the shared arbitrary-q certificate.
    for mutation in (Lr*h, 2s.Psi, L*h1, 2m*O/(r*w)*L*h)
        @test !FR.certify_symbolic("master negative control","test",s.residuals[8]+mutation;verbose=false).pass
    end
end

@testset "FR24 v1 exact power-law background" begin
    certs = FR.certify_fr24_background()
    @test all(c -> c.pass,certs)
    for label in ("eqn:surface_density","eqn:sound_speed","eqn:eos","eqn:temperature",
                  "eqn:aspect_ratio","eqn:fluid_frequency","unlabelled-eq-8","eqn:thermal_scale")
        @test any(c -> c.eqref == "2407.20398v1::$label",certs)
    end
    @variables r h p q
    d = FR.certificate_disc(h,p,q)
    @test !FR.certify_symbolic("kappa missing 2-q negative control","test",
        FR.Modes.kappa(d,r)^2-(1-h^2*(p+q)*r^(1-q))/r^3;verbose=false).pass
    for T in (Float64,BigFloat)
        setprecision(192) do
            for (p,q,r) in ((T(2)/5,T(3)/4,T(3)/5),(T(3)/2,T(1)/5,T(31)/10))
                d = FR.certificate_disc(T(2)/25,p,q)
                @test FR.Modes.Omega(d,r) isa T
                @test FR.Modes.kappa(d,r) isa T
                @test isapprox(FR.Modes.kappa(d,r)^2,
                    4FR.Modes.Omega(d,r)^2+2r*FR.Modes.Omega(d,r)*FR.Modes.dOmegadr(d,r);rtol=64eps(T))
            end
        end
    end
end

using Test
using Symbolics
using FRDerivation

const FR = FRDerivation
include("phi_quadrature.jl")
include("master_background.jl")
include("reconstruction.jl")
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

@testset "FR22 linearised polar fluid equations" begin
    certs = FR.certify_fr22_linearised()
    @test all(c -> c.pass && c.modality == :symbolic,certs)
    for paper in ("2207.14637v2","2407.20398v1","2506.19917v3"), label in
        ("eq:density_pert_eqn","eq:ur_eqn","eq:uphi_eqn")
        @test any(c -> c.eqref == "$paper::$label",certs)
    end
    s = FR.fr22_linearised_residuals()
    (; Omega,uphi,P_r,ds,Sigma,m,omega_lm,ur) = s.symbols
    # Against the derived Euler equation, drop one Coriolis contribution,
    # reverse the density-pressure term, or use pattern instead of total frequency.
    for mutant in (s.printed[2]+Omega*uphi,
                   s.printed[2]+2P_r*ds/Sigma^2,
                   s.printed[2]+im*(m-1)*omega_lm*ur)
        @test !FR.certify_symbolic("linearised negative control","test",
            s.derived[2]-mutant;verbose=false).pass
    end
end

@testset "FR22 pressureless resonances and exceptional modes" begin
    certs = FR.certify_fr22_lindblad()
    @test all(c -> c.pass && c.modality == :symbolic,certs)
    @test any(c -> occursin("E-ECC-8",c.name),certs)
    @test any(c -> occursin("apsidal",c.name),certs)
    @variables m l n
    # The purported inner radius cannot satisfy the OUTER sign equation.
    @test !FR.certify_symbolic("swapped inner sign","test",
        l*n-(m+1)*l*n/(m-1);verbose=false).pass
end

@testset "FR22 potential factorisation certificate" begin
    certs = FR.certify_fr22_phi_factorisation()
    @test all(c -> c.pass && c.modality == :symbolic,certs)
    @test any(c -> occursin("radial derivatives",c.name),certs)
    @variables n pi_exact b rp Gmu C
    expected = -Gmu*b*C/(2pi_exact*rp)
    @test !FR.certify_symbolic("missing frequency negative control","test",
        -Gmu/(2pi_exact^2)/n*(pi_exact*b/rp)*C-expected;verbose=false).pass
    @test !FR.certify_symbolic("Laplace normalization negative control","test",
        -Gmu*n/(2pi_exact^2)/n*(2pi_exact*b/rp)*C-expected;verbose=false).pass
end

check_phi_factorisation()

@testset "FR22 Fourier decomposition and E-ECC-2" begin
    certs = FR.certify_fr22_fourier()
    @test all(c -> c.pass && c.modality == :symbolic, certs)
    @test any(c -> occursin("E-ECC-2 printed inverse",c.name),certs)
    for paper in ("2407.20398v1","2506.19917v3"), label in ("eqn:potential_expansion","eqn:modal_repsonse")
        @test any(c -> c.eqref == "$paper::$label",certs)
    end
    @variables n pi_exact
    cc = FR.real_product_coefficients(Num(1),Num(1))
    inner = (2pi_exact)^2/n*(cc.pm+cc.mp)
    @test !FR.certify_symbolic("printed Fourier gain negative control", "test",
        inner/(2pi_exact^2)-1;verbose=false).pass
    # Exact nonunit frequency exposes the normalization hidden by n=1.
    @test FR.exact_eval(n/(2pi_exact^2)*inner,Dict(n=>FR.RatB(7,3),pi_exact=>FR.RatB(1))) == 1
    @test FR.exact_eval(inner/(2pi_exact^2),Dict(n=>FR.RatB(7,3),pi_exact=>FR.RatB(1))) == 3//7
end

@testset "FR24 and FD25 shared second-order averages" begin
    certs = FR.certify_fr24_second_order_averages()
    @test all(c -> c.pass && c.modality == :symbolic, certs)
    for paper in ("2407.20398v1","2506.19917v3"), label in ("eq:AMF-def","eq:F_jml","eq:dTdr_jml")
        @test any(c -> c.eqref == "$paper::$label",certs)
    end
    @variables ar ai br bi C S pi_exact Phi m
    A, B = ar+im*ai, br+im*bi
    p = FR.real_product_coefficients(A,B)
    integral = 2pi_exact*(p.pm*(C-im*S)+p.mp*(C+im*S))
    # Unequal l must retain time dependence after integrating only phi.
    @test !FR.certify_symbolic("cross-l must not be discarded", "test", integral;verbose=false).pass
    @test !FR.certify_symbolic("instantaneous sine sign negative control", "test",
        integral-pi_exact*(real(A*conj(B))*C-imag(A*conj(B))*S);verbose=false).pass
    @test !FR.certify_symbolic("4pi weight negative control", "test",
        2pi_exact*(p.pm+p.mp)-4pi_exact*real(A*conj(B));verbose=false).pass
    q = FR.real_product_coefficients(A,im*m*Phi)
    @test !FR.certify_symbolic("torque sign negative control", "test",
        -2pi_exact*(q.pm+q.mp)-pi_exact*m*Phi*ai;verbose=false).pass
    # Exact type inheritance in the common coefficient algebra.
    @test FR.real_product_coefficients(1.0+2.0im,3.0-4.0im).pm isa ComplexF64
    setprecision(192) do
        pbig = FR.real_product_coefficients(Complex(big"1",big"2"),Complex(big"3",big"-4"))
        @test pbig.pm isa Complex{BigFloat}
        @test precision(real(pbig.pm)) == 192
    end
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

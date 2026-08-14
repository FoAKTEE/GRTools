# Puncture.Static — the static h0 first-order field (l-mode static coefficients h0_static over the
# auto-generated coefficient table), the elliptic K/E CONSUMED from NumericsCore (AHA: the legacy
# _ellipk_agm/_ellipe_agm fork is dropped onto NumericsCore.ellipK/ellipE), and the worldline
# Schwarzschild connection CONSUMED from Spacetime.
#
# Grounding:
#   • h0_static vs BAKED literal expected values from the CAPTURED Punctures oracle CSV
#     (packages/Punctures/test/h0_reference.csv: columns i,l,p,r0,M,value at N[…,30] from the live
#     Mathematica oracle gen_h0.wl). A representative spread is baked here as literals at the LEGACY
#     tolerance (Float64 rtol 1e-11 for p≤3 well-conditioned families, 1e-3 for the high-p
#     catastrophic-cancellation families; BigFloat rtol 1e-25 — verbatim from Punctures/test/test_static.jl).
#   • elliptic-collapse identity: Puncture.ellipK/ellipE === NumericsCore.ellipK/ellipE (no local AGM
#     copy) + analytic K(0)=E(0)=π/2, E(1)=1 cross-checking the consumed NumericsCore primitive.
#   • Spacetime christoffel object-identity for the worldline connection.
#   • BigFloat precision inheritance of h0_static (parsed-BigFloat oracle, rtol 1e-25).

using Test
using SelfForce
import NumericsCore
import Spacetime

const PU = SelfForce.Puncture

# (i, l, p, r0) => N[h0[i,l,p,r0,1],30]   — baked verbatim from packages/Punctures/test/h0_reference.csv.
const _H0_REF = Dict(
    (1,0,0,10.0) => 7.07845371273137952227994420147,
    (1,2,0,10.0) => 7.16885585553027236068785485534,
    (1,4,0,8.0)  => 3.77834748841002271961443291537,
    (1,8,0,6.5)  => 2.85197188534665418215502828119,
    (3,0,0,10.0) => 6.87764557890936711234959761635,
    (3,4,0,8.0)  => 5.69960017355073503332583617513,
    (3,2,2,10.0) => -0.145379823285050582136076639554,
    (5,0,0,10.0) => 0.0,
    (5,1,0,8.0)  => 1.33526354778224487495398284479,
    (5,2,2,10.0) => 0.389241268170358005272750354206,
    (6,0,0,8.0)  => -0.0862694125405939037822413648291,
    (6,2,0,10.0) => 1.06983202418721907710514277565,
    (7,0,0,10.0) => 0.0,
    (7,1,0,10.0) => 0.0,
    (7,2,0,8.0)  => -1.25442101279923257783316905963,
    (2,1,1,10.0) => 0.747760439733015598857693232124,
    (2,2,1,8.0)  => 1.01318950274883005010943216688,
    (4,1,1,8.0)  => 5.64806105332527498186648840132,
    (9,2,0,10.0) => 0.0,
    (10,2,0,8.0) => 0.0,
)

# high-p families: catastrophic cancellation in Float64 (legacy gates these at rtol 1e-3).
const _H0_REF_HIGHP = Dict(
    (1,4,4,8.0) => -0.00085813318735468287972205256707,
    (3,8,6,10.0) => -0.0000410421988598137284779991459481,
    (2,8,5,8.0) => 0.000358693422193269673754478139122,
)

@testset "Puncture.Static" begin

    @testset "h0_static vs baked Punctures CSV oracle (Float64, p≤3 rtol 1e-11)" begin
        for ((i,l,p,r0), v) in _H0_REF
            got = PU.h0_static(i, l, p, r0, 1.0)
            @test got isa Float64
            @test isapprox(got, v; rtol=1e-11, atol=1e-13)
        end
    end

    @testset "h0_static high-p cancellation families (Float64 rtol 1e-3)" begin
        for ((i,l,p,r0), v) in _H0_REF_HIGHP
            @test isapprox(PU.h0_static(i, l, p, r0, 1.0), v; rtol=1e-3, atol=1e-13)
        end
    end

    @testset "h0_static spot-checks (legacy test literals, rtol 1e-14)" begin
        @test isapprox(PU.h0_static(1, 2, 0, 10.0, 1.0), 7.1688558555302723607; rtol=1e-14)
        @test isapprox(PU.h0_static(3, 4, 0, 8.0, 1.0), 5.6996001735507350333; rtol=1e-14)
    end

    @testset "unavailable (i,p) family throws" begin
        @test_throws ArgumentError PU.h0_static(1, 2, 1, 10.0, 1.0)   # wrong parity (i=1 has no p=1)
        @test_throws ArgumentError PU.h0_static(11, 2, 0, 10.0, 1.0)  # no such i
    end

    @testset "elliptic K/E CONSUMED from NumericsCore (collapse identity + analytic)" begin
        @test PU.ellipK === NumericsCore.ellipK
        @test PU.ellipE === NumericsCore.ellipE
        @test isapprox(real(PU.ellipK(0.0)), pi/2; rtol=1e-14)   # K(0)=π/2
        @test isapprox(real(PU.ellipE(0.0)), pi/2; rtol=1e-14)   # E(0)=π/2
        @test isapprox(real(PU.ellipE(1.0)), 1.0;  rtol=1e-14)   # E(1)=1
        @test isapprox(real(PU.ellipK(0.5)), 1.8540746773013719; rtol=1e-13)  # m=k² convention
    end

    @testset "worldline Schwarzschild connection CONSUMED from Spacetime (object identity)" begin
        @test PU.christoffel === Spacetime.christoffel
        @test PU.christoffel_dr === Spacetime.christoffel_dr
        @test PU.christoffel_dtheta === Spacetime.christoffel_dtheta
    end

    @testset "BigFloat precision inheritance of h0_static (rtol 1e-25)" begin
        setprecision(BigFloat, 256) do
            r0 = BigFloat("10."); M = BigFloat("1")
            expected = BigFloat("7.07845371273137952227994420147")
            got = PU.h0_static(1, 0, 0, r0, M)
            @test got isa BigFloat
            @test isapprox(got, expected; rtol=BigFloat(10)^-25, atol=BigFloat(10)^-30)
            # a high-p family BigFloat resolves the Float64 cancellation to full precision
            r0b = BigFloat("8."); exp2 = BigFloat("-0.00085813318735468287972205256707")
            @test isapprox(PU.h0_static(1, 4, 4, r0b, M), exp2; rtol=BigFloat(10)^-25, atol=BigFloat(10)^-30)
        end
    end
end

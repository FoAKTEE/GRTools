# Geodesics/constants.jl — grounded by BAKED characterization vs the KerrGeodesics oracle
# (values captured once from KerrGeoEnergy/AngularMomentum/CarterConstant, Phase7-B) at the legacy
# runtests tolerance, plus the closed-form Schwarzschild circular-equatorial limit. Spacetime's
# geodesic_energy/angular_momentum/carter_constant must equal the frozen oracle references.
using Test
using Spacetime
include(joinpath(@__DIR__, "fixtures", "spacetime_constants_ref.jl"))  # const CONSTANTS_REF

@testset "Geodesics: constants of motion" begin
    # Representative (a,p,e,x): a=0 Schwarzschild, circular e=0, equatorial x=1, generic inclined,
    # near-extremal, retrograde a<0. x = cos(inclination).
    cases = [
        (0.0,   10.0, 0.0, 1.0),  # Schwarzschild circular equatorial
        (0.0,    8.0, 0.3, 1.0),  # Schwarzschild eccentric equatorial
        (0.9,   10.0, 0.0, 1.0),  # Kerr circular equatorial prograde
        (0.9,   10.0, 0.1, 0.5),  # generic inclined eccentric
        (0.998,  6.0, 0.2, 0.8),  # near-extremal generic
        (-0.7,  12.0, 0.4, 0.6),  # retrograde generic
        (0.5,    7.0, 0.0, 0.3),  # circular inclined
    ]
    @testset "baked vs KerrGeodesics oracle (captured once): (a,p,e,x)=$c" for (i, c) in enumerate(cases)
        a, p, e, x = c
        ref = CONSTANTS_REF[i]
        @test isapprox(geodesic_energy(a, p, e, x),           ref.E;   rtol = 1e-12)
        @test isapprox(geodesic_angular_momentum(a, p, e, x), ref.Lz;  rtol = 1e-12, atol = 1e-13)
        @test isapprox(carter_constant(a, p, e, x),           ref.Q;   rtol = 1e-12, atol = 1e-13)
        com = constants_of_motion(a, p, e, x)
        @test com.E  == geodesic_energy(a, p, e, x)
        @test com.Lz == geodesic_angular_momentum(a, p, e, x)
        @test com.Q  == carter_constant(a, p, e, x)
    end

    @testset "analytic Schwarzschild circular-equatorial limit (a=0,e=0,x=1)" begin
        # E = (1-2M/r)/sqrt(1-3M/r) = (p-2)/sqrt(p(p-3)); Lz = sqrt(Mr)/sqrt(1-3M/r) = p/sqrt(p-3); Q = 0.
        for p in (6.5, 8.0, 10.0, 20.0)
            @test isapprox(geodesic_energy(0.0, p, 0.0, 1.0),           (p - 2) / sqrt(p * (p - 3)); rtol = 1e-13)
            @test isapprox(geodesic_angular_momentum(0.0, p, 0.0, 1.0), p / sqrt(p - 3);             rtol = 1e-13)
            @test isapprox(carter_constant(0.0, p, 0.0, 1.0),           0.0;                         atol = 1e-12)
        end
    end

    @testset "optional-arg aliases ignore passed E/Lz (KerrBackwards-wrapper signature)" begin
        a, p, e, x = 0.9, 10.0, 0.1, 0.5
        E  = geodesic_energy(a, p, e, x)
        Lz = geodesic_angular_momentum(a, p, e, x)
        @test geodesic_angular_momentum(a, p, e, x, E) == Lz
        @test carter_constant(a, p, e, x, E)     == carter_constant(a, p, e, x)
        @test carter_constant(a, p, e, x, E, Lz) == carter_constant(a, p, e, x)
    end

    @testset "exact polar endpoint x=0 (BHPToolkit formulas)" begin
        # Captured from the authoritative local BHPToolkit ConstantsOfMotion.m polar branches.
        circular = constants_of_motion(0.9, 10.0, 0.0, 0.0)
        @test circular.E  ≈ 0.9559116358210666 rtol = 2e-15
        @test circular.Lz === 0.0
        @test circular.Q  ≈ 14.175086162300348 rtol = 2e-15

        eccentric = constants_of_motion(0.9, 10.0, 0.2, 0.0)
        @test eccentric.E  ≈ 0.9574725410737306 rtol = 2e-15
        @test eccentric.Lz === 0.0
        @test eccentric.Q  ≈ 14.25031969387193 rtol = 2e-15

        # Independent turning-point check of the Kerr radial potential at peri/apoapsis.
        radial_potential(r, c) =
            (c.E * (r^2 + 0.9^2) - 0.9 * c.Lz)^2 -
            (r^2 - 2r + 0.9^2) * (r^2 + (c.Lz - 0.9 * c.E)^2 + c.Q)
        for r in (10.0 / 1.2, 10.0 / 0.8)
            @test abs(radial_potential(r, eccentric)) < 1e-10
        end

        # The exact branch joins the generic Schmidt expression continuously from both signs of x.
        for e in (0.0, 0.2)
            polar = constants_of_motion(0.9, 10.0, e, 0.0)
            for x in (-1e-6, 1e-6)
                nearby = constants_of_motion(0.9, 10.0, e, x)
                @test nearby.E ≈ polar.E rtol = 1e-8
                @test nearby.Q ≈ polar.Q rtol = 5e-7
                @test abs(nearby.Lz) < 4e-6
            end
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        a, p, e, x = big"0.9", big"10.0", big"0.1", big"0.5"
        @test geodesic_energy(a, p, e, x) isa BigFloat
        @test constants_of_motion(a, p, e, x).Q isa BigFloat
        polar = constants_of_motion(a, p, e, zero(BigFloat))
        @test polar.E isa BigFloat
        @test polar.Lz isa BigFloat && iszero(polar.Lz)
        @test polar.Q isa BigFloat
    end
end

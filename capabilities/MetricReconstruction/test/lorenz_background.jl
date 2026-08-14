# Grounding for Lorenz.Background — circular-equatorial Kerr orbit scalars.
#
# (1) ANALYTIC Bardeen closed forms: gate E,L,Ω,ut at a=0 (Schwarzschild) and a≠0 (Bardeen) rtol 1e-12;
#     prograde Ω=1/(r^{3/2}+a); a=0 ⇒ E=(r−2)/√(r(r−3)), Ω=r^{−3/2}, ut=1/√(1−3/r); ω=mΩ.
# (2) ONE-SHOT KLC orbit_background fixture at (a=0.99,r0=6,m=2) — calA/calB have no simple closed form
#     (captured via julia --project=packages/KerrLorenzCirc scratchpad/cap12.jl).
# (3) BigFloat type-generic check (orbit_background with BigFloat r0 → {BigFloat} fields).

using Test
import MetricReconstruction.Lorenz as L

@testset "Lorenz.Background (circular-equatorial orbit scalars)" begin

    @testset "a=0 Schwarzschild closed forms (rtol 1e-12)" begin
        for r0 in (6.0, 10.0, 8.3)
            bg = L.orbit_background(0.0, r0, 2)
            @test isapprox(bg.E,  (r0-2)/sqrt(r0*(r0-3)); rtol=1e-12)
            @test isapprox(bg.L,  r0/sqrt(r0-3);          rtol=1e-12)   # L = r²Ω/√(1−3/r) = √r²/√(r−3)
            @test isapprox(bg.Ω,  r0^(-3/2);              rtol=1e-12)
            @test isapprox(bg.ut, 1/sqrt(1 - 3/r0);       rtol=1e-12)
            @test isapprox(bg.ω,  2*bg.Ω;                 rtol=1e-14)
            @test bg.aω == 0.0
            @test isapprox(bg.Δ0, r0^2 - 2r0;             rtol=1e-14)
        end
    end

    @testset "a≠0 Bardeen closed forms (rtol 1e-12)" begin
        for (a, r0) in ((0.99, 6.0), (0.5, 8.0), (-0.7, 7.0))
            bg = L.orbit_background(a, r0, 2)
            v = 1/sqrt(r0); den = sqrt(1 - 3v^2 + 2a*v^3)
            @test isapprox(bg.Ω,  1/(r0^(3/2) + a);               rtol=1e-12)   # prograde Bardeen
            @test isapprox(bg.E,  (1 - 2v^2 + a*v^3)/den;         rtol=1e-12)
            @test isapprox(bg.L,  r0*v*(1 - 2a*v^3 + a^2*v^4)/den; rtol=1e-12)
            @test isapprox(bg.calB, bg.L - a*bg.E;               rtol=1e-12)
            @test isapprox(bg.ω,  2*bg.Ω;                        rtol=1e-14)
            @test isapprox(bg.aω, a*bg.ω;                        rtol=1e-14)
        end
    end

    @testset "KLC one-shot fixture orbit_background(0.99,6.0,2)" begin
        bg = L.orbit_background(0.99, 6.0, 2)
        @test isapprox(bg.E,    0.9213415490283645;  rtol=1e-13)
        @test isapprox(bg.L,    2.744059898821121;   rtol=1e-13)
        @test isapprox(bg.Ω,    0.06374730179252742; rtol=1e-13)
        @test isapprox(bg.ut,   1.3397370360679266;  rtol=1e-13)
        @test isapprox(bg.Δ0,   24.9801;             rtol=1e-13)
        @test isapprox(bg.calA, 1.255186461118687;   rtol=1e-13)
        @test isapprox(bg.calB, 1.8319317652830402;  rtol=1e-13)
        @test isapprox(bg.ω,    0.12749460358505485; rtol=1e-13)
        @test isapprox(bg.aω,   0.1262196575492043;  rtol=1e-13)
    end

    @testset "BigFloat inheritance" begin
        setprecision(BigFloat, 256) do
            bg = L.orbit_background(big"0.99", big"6.0", 2)
            @test bg isa L.OrbitBackground{BigFloat}
            @test bg.calA isa BigFloat && bg.calB isa BigFloat && bg.ut isa BigFloat
            # 256-bit value matches the Float64 fixture to ~Float64 precision
            @test abs(bg.calA - big"1.255186461118687") < 1e-12
        end
    end
end

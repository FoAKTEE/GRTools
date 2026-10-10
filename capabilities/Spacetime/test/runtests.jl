using Test
using Spacetime

@testset "Spacetime" begin
    @testset "smoke" begin
        @test isdefined(Main, :Spacetime)
        @test Spacetime isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module test files are included as each submodule is TDD-ported and grounded
    # against the light legacy oracle (KerrGeodesics / SpinningOrbit) at its frozen tolerance,
    # or against an analytic identity (Schwarzschild ISCO = 6M, separatrix p = 6 + 2e at a=0,
    # metric det/inverse identities, geodesic-constant consistency).
    include("test_geodesics_constants.jl")
    include("test_negative_spin_constants.jl")
    include("test_geometry.jl")
    include("test_geodesics_frequencies.jl")
    include("equatorial_quadrature_oracle.jl")
    include("test_trajectory.jl")
    include("test_equatorial_branch.jl")
    include("test_turning_point_jacobian.jl")
    include("test_eccentric_harmonics.jl")
    include("test_separatrix.jl")
    include("test_christoffel.jl")
    include("test_schwarzschild_orbit.jl")
    include("test_tortoise_kerr.jl")
    include("test_circular_kinematics.jl")
    include("test_kerr_equatorial_potentials.jl")   # EGD: Page–Thorne (ν,ψ,μ,ω̃) + d/dr + circular dt/dτ
    include("test_pn_kinematics.jl")
    include("test_spinning_secondary.jl")
    include("test_spinning_actions.jl")   # linear-in-spin action parts + J/pex gauge transforms
    include("test_spinning_forces.jl")
    include("test_npghp.jl")            # NPGHP campaign Move B: G-B1..G-B5 (spin coeffs + {p,q} + HPHC chart)
end

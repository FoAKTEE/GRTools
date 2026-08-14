using Test
using ModifiedGravity

@testset "ModifiedGravity" begin
    @testset "smoke" begin
        # Scaffold smoke test: the package loads and its module is in scope.
        @test ModifiedGravity isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module grounded tests (characterized against frozen .dat source tables, wolframscript
    # analytic kernels, the RadHzpotmode paclet oracle, and the Sourcesamplitudes replication).
    include("test_ghp_operators.jl")
    include("test_source_tables.jl")
    include("test_source_kernels.jl")
    include("test_hertz_potential.jl")
    include("test_source_amplitudes.jl")
    include("test_rhs_amplitudes.jl")
    include("test_ring_dynamics.jl")
    include("test_weyl_from_hertz.jl")
end

using Test
using TestKit

@testset "TestKit" begin
    @testset "smoke" begin
        # `using TestKit` loaded the module + exported the public harness surface.
        for f in (:oracle_records, :reset_oracle_log!, :tolerance_class, :assert_can_only_tighten,
                  :load_fixture, :characterization_snapshot, :characterization_gate,
                  :assert_dep_budget, :assert_not_loaded,
                  :strict_oracles, :oracle_root, :oracle_or_skip)
            @test isdefined(TestKit, f)
        end
        @test isdefined(TestKit, Symbol("@approx_oracle"))
        @test TestKit.TOLERANCE_CLASSES isa AbstractDict
    end

    include("hygiene.jl")        # §7.7 gates for TestKit's OWN dep budget (before heavy tests)
    include("test_oracle.jl")
    include("test_fixtures.jl")
    include("test_hygiene.jl")   # unit tests of the §7.7 checkers themselves
    include("test_oracles.jl")   # unit tests of the GRTOOLS_* oracle gate
end

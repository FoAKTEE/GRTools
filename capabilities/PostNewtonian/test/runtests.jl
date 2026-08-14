using Test
using PostNewtonian

@testset "PostNewtonian" begin
    @testset "smoke" begin
        # The capability loads and exposes its submodules (Symbolics quarantine resolves).
        @test isdefined(PostNewtonian, :Series)
        @test isdefined(PostNewtonian, :Catalog)
        @test isdefined(PostNewtonian, :MMAInputForm)
        @test isdefined(PostNewtonian, :SelfForce)
        @test isdefined(PostNewtonian, :MST)
        @test isdefined(PostNewtonian, :RadialMode)
        @test isdefined(PostNewtonian, :SourceMode)
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module characterization suites, grounded by the captured MMA oracle / .m-file oracle /
    # analytic PN identities (NEVER a heavy live legacy dep).
    include("test_series.jl")
    include("test_catalog.jl")
    include("test_mmainputform.jl")
    include("test_selfforce.jl")
    include("test_mst.jl")
    include("test_radialmode.jl")
    include("test_sourcemode.jl")
end

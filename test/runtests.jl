# Umbrella smoke test: the suite loads, every subpackage binding is a module, and
# the tier grouping is complete. Package-level physics tests live in each
# subpackage's own `test/` (run them with `Pkg.test` in that package's project).
using Test
using GrTools

@testset "GrTools umbrella" begin
    @test GrTools isa Module
    for m in (GrTools.PACKAGES.core..., GrTools.PACKAGES.capabilities...)
        @test m isa Module
    end
    @test length(GrTools.PACKAGES.core) == 7
    @test length(GrTools.PACKAGES.capabilities) == 11
    # The qualified bindings match the loaded modules.
    @test GrTools.Spacetime === GrTools.PACKAGES.capabilities[1]
    @test GrTools.NumericsCore === GrTools.PACKAGES.core[1]
end

# Umbrella smoke test: the suite loads, every subpackage binding is a module, and
# the tier grouping is complete. Package-level physics tests live in each
# subpackage's own `test/` (run them with `Pkg.test` in that package's project).
using Test
using GrTool

@testset "GrTool umbrella" begin
    @test GrTool isa Module
    for m in (GrTool.PACKAGES.core..., GrTool.PACKAGES.capabilities...)
        @test m isa Module
    end
    @test length(GrTool.PACKAGES.core) == 7
    @test length(GrTool.PACKAGES.capabilities) == 11
    # The qualified bindings match the loaded modules.
    @test GrTool.Spacetime === GrTool.PACKAGES.capabilities[1]
    @test GrTool.NumericsCore === GrTool.PACKAGES.core[1]
end

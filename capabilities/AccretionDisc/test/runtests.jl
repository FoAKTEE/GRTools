using Test
using AccretionDisc

@testset "AccretionDisc" begin
    @testset "smoke" begin
        @test isdefined(AccretionDisc, :Modes)
    end

    include("hygiene.jl")

    include("Modes/test_background.jl")
    include("Modes/test_velocity.jl")
    include("Modes/test_offequator.jl")
    include("Modes/test_masterpde.jl")
    include("Modes/test_resonances.jl")
    include("Modes/test_reducedode.jl")
    include("Modes/test_forcing.jl")
    include("Force/test_loader.jl")
    include("Force/test_assemble.jl")
    include("Force/test_source_fluid.jl")
    include("Force/test_pipeline.jl")
    include("Force/test_source_grid.jl")
    include("Force/test_bigfloat.jl")
    include("Force/test_kerr_frontend_numerics.jl")
    include("Force/test_kerr.jl")
    include("SpiralWaves/test_scaffold.jl")   # arXiv:2608.17003 spiral-wave submodule: load + dep budget
    include("SpiralWaves/test_background.jl") # EGD::sw-background: disc equilibrium (Sec. II.A)
end

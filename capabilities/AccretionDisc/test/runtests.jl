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
    include("Modes/test_pressureless_lindblad.jl")  # EGD::analytic-lindblad-hegade: arXiv:2510.03564 (18a,18b,19,20,23,24,25,C1-C2)
    include("Modes/test_thindisc.jl")   # NT/SS thin-disc structure + total Lindblad torque (2510.02433)
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
    include("SpiralWaves/test_grid.jl")       # EGD::sw-grid: adaptive corotation/WKB grid (App. C)
    include("SpiralWaves/test_master.jl")     # EGD::sw-master-matrices: master system C, S (47), (50), Landau (62-63)
    include("SpiralWaves/test_source.jl")     # EGD::sw-singular-source: singular field (52)-(60), (B1), Fourier modes
    include("SpiralWaves/test_r2_generated.jl") # EGD::sw-mode-solver: route-R2 generated δû_φ, AMF (72a), torque density (72b)
    include("SpiralWaves/test_solver.jl")     # EGD::sw-mode-solver: implicit midpoint (64) + outgoing-wave BCs (66)-(68)
end

using Test
using EMRIWaveforms

@testset "EMRIWaveforms" begin
    @testset "smoke" begin
        # `using EMRIWaveforms` resolved; submodules present.
        @test EMRIWaveforms isa Module
        @test isdefined(EMRIWaveforms, :Mappings)
        @test isdefined(EMRIWaveforms, :Trajectory)
        @test isdefined(EMRIWaveforms, :FluxModels)
        @test isdefined(EMRIWaveforms, :AmplitudeModels)
        @test isdefined(EMRIWaveforms, :ModeSelect)
        @test isdefined(EMRIWaveforms, :ModeSum)
        @test isdefined(EMRIWaveforms, :AAK)
        @test isdefined(EMRIWaveforms, :Generate)
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # per-module characterization suites are added tests-first as each module lands:
    include("test_pn5_fluxes.jl")
    include("test_mappings.jl")
    include("test_fluxmodels.jl")
    include("test_kerrampspline.jl")
    include("test_amplitudemodels.jl")
    include("test_modeselect.jl")
    include("test_modesum.jl")
    include("test_aak.jl")
    include("test_trajectory.jl")
    include("test_generate.jl")
end

using Test
using Data

@testset "Data" begin
    @testset "smoke" begin
        # `using Data` loaded the module and re-exported the dep-free schema contract.
        @test isdefined(Data, :H1Schema)
        @test isdefined(Data, :AmpGridSchema)
        @test isdefined(Data, :FluxGridSchema)
        @test Data.SCHEMA_VERSION isa Integer
        # I/O entry points exist (bodies in weakdep ext / future round).
        for f in (:read_h1, :write_mode_h5, :read_ampgrid, :read_fluxgrid,
                  :read_kerrampgrid, :read_roman_nn, :kerr_mode_table, :KerrAmpGrid,
                  :read_dat, :read_mma_inputform, :read_wl_table, :wl_to_jld2, :read_jld2,
                  :set_datadir!, :datadir, :artifact_path)
            @test isdefined(Data, f)
        end
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module test files (added tests-first as each module is ported).
    include("test_datadir.jl")
    include("test_dlm.jl")
    include("test_wl.jl")
    include("test_h1io.jl")
    include("test_ampgrid.jl")
    include("test_fluxgrid.jl")
    include("test_kerrampgrid.jl")
    include("test_jld2io.jl")
    include("test_artifacts.jl")
end

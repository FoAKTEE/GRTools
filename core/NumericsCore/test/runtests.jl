using Test
using NumericsCore

@testset "NumericsCore" begin
    @testset "smoke" begin
        @test isdefined(Main, :NumericsCore)
        @test NumericsCore isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests
    include("test_carlson.jl")
    include("test_precision.jl")
    include("test_modeindex.jl")
    include("test_units.jl")
    include("test_blas_guard.jl")
    include("test_rootfind.jl")
    include("test_jets.jl")
    include("test_continued_fraction.jl")
    include("test_dataschema.jl")
    include("test_interp.jl")
    include("test_chebyshev.jl")
    include("test_ode.jl")
    include("test_legendre.jl")
    include("test_kummer.jl")
    include("test_dg.jl")
    include("test_fd_uniform.jl")
    include("test_backend_iface.jl")
end

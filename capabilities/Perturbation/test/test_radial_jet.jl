# Gate G-C1 (Move C): the radial_jet derivative-tower relocated from MetricReconstruction.Radiation
# into Perturbation.Radial must reproduce the pre-move outputs bit-for-bit (Float64) and preserve
# the {T}-generic BigFloat-precision-keyed path (BigFloat-256). Oracle = leg_radial_jet.jl, captured
# ONCE from the pre-move MetricReconstruction.Radiation.radial_jet (see scratchpad/capture_radial_jet.jl).
#
#   • Float64 : new Perturbation.radial_jet == fixture, rel <= 1e-15 (verbatim-port bit-fidelity).
#   • BigFloat-256 : rel <= 1e-40 (inherited arbitrary precision; the cache key carries
#     precision(BigFloat), so different setprecision regimes never collide).
#   • eltype preserved : Complex{Float64} / Complex{BigFloat} follow the input eltype.
# Also exercises clear_radial_jet_caches! (the precision-keyed cache reset that moved with radial_jet).

using Test
using Perturbation
using Perturbation: radial_jet, clear_radial_jet_caches!

include(joinpath(@__DIR__, "fixtures", "leg_radial_jet.jl"))

@testset "Perturbation.radial_jet (Move C fixture identity — G-C1)" begin
    @testset "Float64 == fixture (<=1e-15 rel) + eltype" begin
        for c in RADIAL_JET_CASES
            a = parse(Float64, c.a); ω = parse(Float64, c.omega); r = parse(Float64, c.r)
            got = radial_jet(c.s, c.l, c.m, a, ω, r; bc=:Up, order=4)
            @test eltype(got) == Complex{Float64}
            @test length(got) == 5
            for k in 1:5
                @test isapprox(got[k], c.f64[k]; rtol=1e-15, atol=1e-300)
            end
        end
    end

    @testset "BigFloat-256 == fixture (<=1e-40 rel) + eltype" begin
        setprecision(BigFloat, 256) do
            clear_radial_jet_caches!()   # start from a clean precision regime
            for c in RADIAL_JET_CASES
                a = parse(BigFloat, c.a); ω = parse(BigFloat, c.omega); r = parse(BigFloat, c.r)
                got = radial_jet(c.s, c.l, c.m, a, ω, r; bc=:Up, order=4)
                @test eltype(got) == Complex{BigFloat}
                @test length(got) == 5
                for k in 1:5
                    ref = parse(BigFloat, c.bf256[k][1]) + parse(BigFloat, c.bf256[k][2]) * im
                    @test isapprox(got[k], ref; rtol=parse(BigFloat, "1e-40"), atol=parse(BigFloat, "1e-300"))
                end
            end
            clear_radial_jet_caches!()
        end
    end

    @testset "clear_radial_jet_caches! resets the memo tables" begin
        radial_jet(-2, 2, 2, 0.1, 0.3, 6.0)   # populate radial + lambda caches
        @test !isempty(Perturbation._radial_jet_cache)
        @test clear_radial_jet_caches!() === nothing
        @test isempty(Perturbation._radial_jet_cache)
        @test isempty(Perturbation._radial_jet_lambda_cache)
    end
end

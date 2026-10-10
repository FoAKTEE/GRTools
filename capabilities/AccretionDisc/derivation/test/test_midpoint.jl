# Standalone runner (from the GRTools root; no shared registration required):
# julia --startup-file=no --project=capabilities/AccretionDisc/derivation -t1 capabilities/AccretionDisc/derivation/test/test_midpoint.jl
using Test
using LinearAlgebra
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "../../../..")))
isdefined(@__MODULE__, :R3Derivation) || include(joinpath(@__DIR__, "../src/R3Derivation.jl"))
include(joinpath(@__DIR__, "../src/midpoint.jl"))
using .MidpointCertificate
using .MidpointCertificate: certify_symbolic

@testset "C-64 midpoint certificate" begin
    @testset "symbolic derivation and actual assembly" begin
        result = midpoint_certificates()
        @test length(result.certificates) == 7
        @test all(c -> c.pass && c.modality == :symbolic, result.certificates)
        # These deliberately wrong formulas must be rejected: neither sign nor
        # the endpoint-average correction can disappear in the certificate.
        @test !certify_symbolic("wrong defect sign", "negative control",
                                result.defect + result.leading; verbose = false).pass
        @test !certify_symbolic("omit endpoint averaging", "negative control",
                                result.defect - result.quadrature_only; verbose = false).pass
        @test !certify_symbolic("reverse source sign", "negative control",
                                result.source_update + result.source_update;
                                verbose = false).pass
    end

    @testset "actual solve_midpoint, uniform and variable steps" begin
        for T in (Float64, BigFloat)
            setprecision(BigFloat, 128) do
                for variable_steps in (false, true)
                    rows = midpoint_order_check(T; variable_steps)
                    @test length(rows) == 5       # initial mesh + four halvings
                    @test all(row -> row.error isa T && row.error > 0, rows)
                    @test all(row -> row.solution_type == Complex{T}, rows)
                    @test all(row -> row.precision_bits == precision(T), rows)
                    @test all(diff([row.error for row in rows]) .< 0)
                    @test all(abs(row.order - 2) < T(1) / 20 for row in rows[2:end])
                    @test abs(rows[end].order - 2) < T(1) / 200
                    @test abs(rows[end].order - 2) < abs(rows[2].order - 2)
                    @test all(row -> row.algebraic_residual < 1000 * eps(T), rows)
                    @test all(isapprox(rows[i].hmax, 2rows[i+1].hmax;
                                       rtol = 1024eps(T)) for i in 1:4)
                end
            end
        end
    end
end

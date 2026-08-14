# Lightweight numerical-contract tests for Force.KerrFrontend. The expensive Lorenz reconstruction
# is replaced by an injected polynomial evaluator so the test directly observes the radial argument
# type and verifies the Fornberg derivative on a genuinely nonuniform grid.

using Test
using AccretionDisc
using AccretionDisc.Force

@testset "Force.KerrFrontend radial precision and nonuniform derivative" begin
    seen_types = DataType[]
    function polynomial_bls(a, r_p, l, m, r; recs=nothing, S=nothing, ls=nothing)
        push!(seen_types, typeof(r))
        ntuple(j -> Complex((j + 1) * r^4 - j * r^2 + 3r, (2j - 1) * r^3 - r), 10)
    end

    rg = BigFloat[big"3.0", big"3.2", big"3.9", big"5.1", big"6.8", big"9.0"]
    mode = kerr_mode(0.5, 8.0, 2, 2; T=BigFloat, rgrid=rg, threaded=false,
                     bls_evaluator=polynomial_bls)
    @test seen_types == fill(BigFloat, length(rg))
    @test mode.r_right == rg
    # Even loader slot 1 is BLS component h̄¹ (j=1 in the evaluator). A five-node Fornberg
    # stencil differentiates this degree-four polynomial exactly, including one-sided edge stencils.
    expected(r) = Complex{BigFloat}(8r^3 - 2r + 3, 3r^2 - 1)
    @test mode.dr_right[1] ≈ expected.(rg) rtol=big"1e-65" atol=big"1e-65"

    @test_throws ArgumentError kerr_mode(0.5, 8.0, 2, 2; rgrid=[3.0], threaded=false,
                                         bls_evaluator=polynomial_bls)
    @test_throws ArgumentError kerr_mode(0.5, 8.0, 2, 2; rgrid=[3.0, 4.0, 3.5], threaded=false,
                                         bls_evaluator=polynomial_bls)
end

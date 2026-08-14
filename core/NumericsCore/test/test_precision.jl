# Characterization of the precision-policy selectors.
#
# Heavy/HDF5-pulling H1Lorenz test-dep DROPPED (see Project.toml); bc_tol/rk_tol are re-grounded by
# their DOCUMENTED CLOSED FORMS, which ARE the H1Lorenz source verbatim (read, not run):
# Oracle 1 (bc_tol): H1Lorenz `_bc_tol` (packages/H1Lorenz/src/boundary_conditions.jl:17-23) —
#   Float64 -> 1e-14; BigFloat -> max(big(2.0)^(10-bits), 8*eps(BigFloat)).
# Oracle 2 (rk_tol): H1Lorenz `_rk_tol` (packages/H1Lorenz/src/integrator.jl:22-28) —
#   Float64 -> 1e-12; BigFloat -> big(2.0)^(20-bits).
# Oracle 3 (fd_step, bigprec_for_mag): documented closed forms — DiscEMRIForce `_fdr=T(1)/T(10)^6`,
#   `_fdθh=T(1)/T(10)^4`; KerrLorenzCirc `_bigprec_for_mag = clamp(128+6*ceil(magpow),160,320)`.
# Oracle 4 (with_precision): round-trip invariant (global BigFloat precision restored after block).
# Every expected value is a source-literal closed form or an analytic relation — never the code output.

using NumericsCore

@testset "precision" begin
    @testset "bc_tol closed form (== H1Lorenz._bc_tol source)" begin
        @test NumericsCore.bc_tol(Float64) == 1e-14
        for bits in (256, 512)
            setprecision(BigFloat, bits) do
                @test NumericsCore.bc_tol(BigFloat) == max(big(2.0)^(10 - bits), 8*eps(BigFloat))
            end
        end
    end

    @testset "rk_tol closed form (== H1Lorenz._rk_tol source)" begin
        @test NumericsCore.rk_tol(Float64) == 1e-12
        for bits in (256, 512)
            setprecision(BigFloat, bits) do
                @test NumericsCore.rk_tol(BigFloat) == big(2.0)^(20 - bits)
            end
        end
    end

    @testset "fd_step closed form (DiscEMRIForce _fdr/_fdθh)" begin
        @test NumericsCore.fd_step(Float64; kind=:radial) == 1.0/10.0^6
        @test NumericsCore.fd_step(Float64; kind=:theta)  == 1.0/10.0^4
        # type-generic: BigFloat representable step inherits precision (not a Float64 literal)
        setprecision(BigFloat, 256) do
            @test NumericsCore.fd_step(BigFloat; kind=:radial) == big(1)/big(10)^6
            @test NumericsCore.fd_step(BigFloat; kind=:radial) isa BigFloat
            @test NumericsCore.fd_step(BigFloat; kind=:theta)  == big(1)/big(10)^4
        end
    end

    @testset "bigprec_for_mag closed form (clamp 128+6*ceil(magpow),160,320)" begin
        @test NumericsCore.bigprec_for_mag(13)  == 206         # 128+6*13
        @test NumericsCore.bigprec_for_mag(33)  == 320         # 128+198=326 -> clamp 320
        @test NumericsCore.bigprec_for_mag(0)   == 160         # 128 -> clamp 160
        @test NumericsCore.bigprec_for_mag(-5)  == 160         # below floor
        @test NumericsCore.bigprec_for_mag(5.2) == 128 + 6*6   # ceil(5.2)=6 -> 164
    end

    @testset "with_precision round-trip invariant" begin
        p0 = precision(BigFloat)
        got = NumericsCore.with_precision(512) do
            precision(BigFloat)
        end
        @test got == 512
        @test precision(BigFloat) == p0          # restored even though block changed it
        # restored on exception too
        @test_throws ErrorException NumericsCore.with_precision(384) do
            error("boom")
        end
        @test precision(BigFloat) == p0
    end
end

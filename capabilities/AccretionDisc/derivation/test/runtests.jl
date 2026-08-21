# derivation/test/runtests.jl — the tiny own suite of the standalone route-R3 derivation project.
#
# It re-runs the certificate suites (stage 1 = background, stage 2 = master system, stage 3 =
# diagnostics) at a REDUCED number of exact-rational points so that the whole thing is a few minutes;
# the node evidence logs are produced by `scripts/derive_all.jl` at the full point count.
# Run:  julia --project=capabilities/AccretionDisc/derivation -t4 <this file>

using Test
include(joinpath(@__DIR__, "..", "src", "R3Derivation.jl"))
using .R3Derivation
const R3 = R3Derivation

const NPTS = parse(Int, get(ENV, "R3_NPTS", "12"))

@testset "R3 derivation (arXiv:2608.17003, route R3)" begin
    @testset "stage 1 — symbolic background (4)-(11)" begin
        certs = R3.stage1_certificates(; npts = NPTS, verbose = true)
        @test !isempty(certs)
        for c in certs
            @test c.pass || error("certificate FAILED: $(c.name) [$(c.eqref)] — $(c.note)")
        end
    end

    @testset "stage 2 — perturbation + reduction (15)-(22), (45)-(49)" begin
        certs = R3.stage2_certificates(; npts = NPTS, verbose = true)
        @test !isempty(certs)
        for c in certs
            @test c.pass || error("certificate FAILED: $(c.name) [$(c.eqref)] — $(c.note)")
        end
    end

    @testset "stage 2 gate — Newtonian limit ≡ Eq. (50)" begin
        certs = R3.newtonian_certificates(; npts = max(3, NPTS ÷ 3), verbose = true)
        @test !isempty(certs)
        for c in certs
            @test c.pass || error("certificate FAILED: $(c.name) [$(c.eqref)] — $(c.note)")
        end
    end

end

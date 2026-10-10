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

    @testset "stage 3 — diagnostics (40)-(42), (69), (72)" begin
        certs = R3.stage3_certificates(; npts = NPTS, verbose = true)
        @test !isempty(certs)
        for c in certs
            @test c.pass || error("certificate FAILED: $(c.name) [$(c.eqref)] — $(c.note)")
        end
    end

    @testset "stage 4 — code emission is deterministic and round-trips exactly" begin
        here = normpath(joinpath(@__DIR__, ".."))
        mart = joinpath(here, "out", "R3_master.jls.txt")
        dart = joinpath(here, "out", "R3_diagnostics.jls.txt")
        tmp1 = mktempdir(); tmp2 = mktempdir()
        try
            m1 = R3.emit_master(joinpath(tmp1, "master_derived.jl"), mart)
            d1 = R3.emit_diagnostics(joinpath(tmp1, "r2_derived.jl"), mart, dart)
            m2 = R3.emit_master(joinpath(tmp2, "master_derived.jl"), mart)
            d2 = R3.emit_diagnostics(joinpath(tmp2, "r2_derived.jl"), mart, dart)
            @test read(m1) == read(m2)        # byte-identical across runs
            @test read(d1) == read(d2)
            # the SHIPPED sources are exactly what this emitter produces (no hand edit ever landed)
            ship = normpath(joinpath(here, "..", "src", "SpiralWaves"))
            @test read(joinpath(ship, "master_derived.jl")) == read(m1)
            @test read(joinpath(ship, "r2_derived.jl")) == read(d1)
            c = R3.codegen_selfcheck(m1, d1, mart, dart; npts = max(4, NPTS), verbose = true)
            @test c.pass || error("codegen round-trip FAILED — $(c.note)")
        finally
            rm(tmp1; recursive = true, force = true)
            rm(tmp2; recursive = true, force = true)
        end
    end
    include("test_li_closure.jl")
    include("test_li_energy_balance.jl")
    include("test_eccentric_conventions.jl")
    include("test_midpoint.jl")
end

# Module 1: datadir — injected data-root seam. GROUNDING: closed-form / round-trip — the regex
# contract + joinpath are the oracle (deterministic), so all expected values are hand-derived.

using Test
using Data

@testset "datadir" begin
    # Start from a clean slate (no process override).
    Data.set_datadir!(nothing)

    @testset "set/get round-trip identity" begin
        tmp = mktempdir()
        ret = Data.set_datadir!(tmp)
        @test ret == abspath(tmp)
        @test Data.datadir() == abspath(tmp)            # get returns what was set
        # Clearing falls back to default resolution (which throws when unconfigured).
        Data.set_datadir!(nothing)
        withenv("GRTOOLS_DATADIR" => nothing) do
            @test_throws ErrorException Data.datadir()
        end
    end

    @testset "default_datadir env resolution" begin
        Data.set_datadir!(nothing)
        tmp = mktempdir()
        withenv("GRTOOLS_DATADIR" => tmp) do
            @test Data.default_datadir() == abspath(tmp)
            @test Data.datadir() == abspath(tmp)        # override unset → env path
        end
        withenv("GRTOOLS_DATADIR" => nothing) do
            @test_throws ErrorException Data.default_datadir()
        end
    end

    @testset "datapath joinpath correctness" begin
        tmp = mktempdir()
        Data.set_datadir!(tmp)
        @test Data.datapath("a", "b.h5") == joinpath(abspath(tmp), "a", "b.h5")
        @test Data.datapath() == abspath(tmp)
        Data.set_datadir!(nothing)
    end

    @testset "available_modes regex contract" begin
        # Hand-named fake mode files → KNOWN sorted (l,m) list (the regex h1-l{l}m{m}.h5 is the oracle).
        d = mktempdir()
        for nm in ("h1-l2m0.h5", "h1-l2m1.h5", "h1-l2m2.h5", "h1-l3m1.h5", "h1-l10m10.h5",
                   # non-matching decoys (must be ignored):
                   "h1-l2m2.dat", "notamode.h5", "h1-lxmy.h5", "README.txt")
            touch(joinpath(d, nm))
        end
        @test Data.available_modes(dir = d) == [(2, 0), (2, 1), (2, 2), (3, 1), (10, 10)]
    end
end

# Unit tests of the out-of-tree oracle gate (repro-audit G3/G4). The gate is TEST infrastructure, so
# it is grounded by exercising it against itself: present ⇒ true and nothing recorded; absent ⇒ false
# and exactly one `Broken`; strict mode ⇒ false and exactly one `Fail`.

using Test
using TestKit: strict_oracles, oracle_root, oracle_or_skip

# Silent collector: records results without PRINTING them (a DefaultTestSet echoes every Fail to
# stderr the moment it is recorded, which would make this suite's own log look broken).
mutable struct Collector <: Test.AbstractTestSet
    results::Vector{Any}
end
Collector(desc::AbstractString) = Collector(Any[])
Test.record(ts::Collector, r) = (push!(ts.results, r); r)
Test.finish(ts::Collector) = ts

@testset "oracle gate" begin

    @testset "oracle_root: ENV override wins, default otherwise" begin
        var = "GRTOOLS_TEST_ORACLE_ROOT_UNIT"
        haskey(ENV, var) && delete!(ENV, var)
        @test oracle_root(var, "/capture/site/default") == "/capture/site/default"
        ENV[var] = "/elsewhere"
        @test oracle_root(var, "/capture/site/default") == "/elsewhere"
        delete!(ENV, var)
        @test oracle_root(var, "/capture/site/default") == "/capture/site/default"
    end

    @testset "strict_oracles reads ENV at call time" begin
        old = get(ENV, "GRTOOLS_STRICT_ORACLES", nothing)
        delete!(ENV, "GRTOOLS_STRICT_ORACLES")
        @test strict_oracles() == false
        for v in ("1", "true", "YES", "on")
            ENV["GRTOOLS_STRICT_ORACLES"] = v
            @test strict_oracles() == true
        end
        for v in ("0", "false", "")
            ENV["GRTOOLS_STRICT_ORACLES"] = v
            @test strict_oracles() == false
        end
        old === nothing ? delete!(ENV, "GRTOOLS_STRICT_ORACLES") : (ENV["GRTOOLS_STRICT_ORACLES"] = old)
    end

    @testset "present oracle: returns true, records nothing" begin
        inner = Collector("inner")
        Test.push_testset(inner)
        local ret
        try
            ret = oracle_or_skip(true, "always-there"; env = "GRTOOLS_X")
        finally
            Test.pop_testset()
        end
        @test ret == true
        @test isempty(inner.results)      # a present oracle records NOTHING
    end

    @testset "absent oracle: returns false, records exactly one Broken" begin
        old = get(ENV, "GRTOOLS_STRICT_ORACLES", nothing)
        delete!(ENV, "GRTOOLS_STRICT_ORACLES")
        inner = Collector("inner")
        Test.push_testset(inner)
        local ret
        try
            ret = oracle_or_skip(false, "read_ampgrid vs raw HDF5"; env = "GRTOOLS_FEWDATA_DIR")
        finally
            Test.pop_testset()
        end
        @test ret == false
        @test length(inner.results) == 1
        r = inner.results[1]
        @test r isa Test.Broken
        @test r.test_type === :skipped
        @test occursin("read_ampgrid vs raw HDF5", string(r.orig_expr))
        @test occursin("GRTOOLS_FEWDATA_DIR", string(r.orig_expr))
        old === nothing ? nothing : (ENV["GRTOOLS_STRICT_ORACLES"] = old)
    end

    @testset "strict mode: absent oracle FAILS instead of skipping" begin
        old = get(ENV, "GRTOOLS_STRICT_ORACLES", nothing)
        ENV["GRTOOLS_STRICT_ORACLES"] = "1"
        inner = Collector("inner")
        Test.push_testset(inner)
        local ret
        try
            ret = oracle_or_skip(false, "read_ampgrid vs raw HDF5"; env = "GRTOOLS_FEWDATA_DIR")
        finally
            Test.pop_testset()
        end
        @test ret == false
        @test length(inner.results) == 1
        @test inner.results[1] isa Test.Fail
        @test !any(r -> r isa Test.Broken, inner.results)
        old === nothing ? delete!(ENV, "GRTOOLS_STRICT_ORACLES") : (ENV["GRTOOLS_STRICT_ORACLES"] = old)
    end
end

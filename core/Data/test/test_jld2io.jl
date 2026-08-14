# Module: jld2io — JLD2 parse-once-cache path for WaSABI `.data` catalogs. GROUNDING: round-trip
# identity. `wl_to_jld2(src,dst)` then `read_jld2(dst)` must reproduce `Data.read_wl_data(src)`
# bit-for-bit (BigFloat digits preserved — WaSABI stores reals as BigFloat). Self-contained: uses the
# in-package `Data.read_wl_data` (green) as the oracle over a hand-authored `.data` fixture; JLD2 is a
# live light test-dep firing DataJLD2Ext — NO heavy WaSABI dep.

using Test
using Data
using JLD2   # fires DataJLD2Ext — wl_to_jld2 / read_jld2 bodies activate

# Bit-for-bit identity for a stored real: BigFloat is a mutable heap struct, so `===` is object
# identity (always false across a round-trip). True storage identity = equal value AND equal mantissa
# precision AND equal sign (covers the full IEEE-bigfloat content the cache must preserve).
_bitsame(a::BigFloat, b::BigFloat) = precision(a) == precision(b) && a == b && signbit(a) == signbit(b)
_bitsame(a, b) = a == b

# A small WL `.data` catalog with full-precision (`prec`-marked) reals + an Around[c,e] (kept as c),
# mirroring the WaSABI SF-data record shape (a list of <|...|> associations).
const _WL_DATA_FIXTURE = """
{<|"a" -> 0.99796390114032518371715537295389520378644361909354834633123951052`101.89771765497302,
   "r" -> 5.4692447091804464349691437654109517798043438795819813021`96.375962*^6,
   "inf" -> Around[1.30780472696481215560721462144676339`112.1458*^-33, 3.2e-34]|>,
 <|"a" -> 0.5`50.0, "r" -> 10.25`50.0, "inf" -> 1.5`50.0|>}
"""

@testset "jld2io" begin
    @testset "wl_to_jld2 / read_jld2 round-trip identity (BigFloat bit-for-bit)" begin
        dir = mktempdir()
        src = joinpath(dir, "fixture.data")
        dst = joinpath(dir, "fixture.jld2")
        write(src, _WL_DATA_FIXTURE)

        oracle = Data.read_wl_data(src)             # in-package green oracle (BigFloat-exact)
        @test oracle isa Vector{Dict{String,Any}}
        @test length(oracle) == 2
        @test oracle[1]["a"] isa BigFloat

        out = Data.wl_to_jld2(src, dst)
        @test out == dst
        @test isfile(dst)

        loaded = Data.read_jld2(dst)                 # default key = "data"
        # Whole-structure identity.
        @test loaded == oracle
        @test loaded isa Vector{Dict{String,Any}}
        # Bit-for-bit on the high-precision reals (=== on BigFloat compares value+precision).
        for (lr, orc) in zip(loaded, oracle)
            for k in keys(orc)
                @test _bitsame(lr[k], orc[k])
            end
        end
        # explicit key path resolves the same value.
        @test Data.read_jld2(dst; key = "data") == oracle
    end

    @testset "read_jld2 generic value round-trip" begin
        dir = mktempdir()
        p = joinpath(dir, "g.jld2")
        v = Dict("x" => big"3.14159265358979323846", "n" => [1, 2, 3])
        JLD2.jldsave(p; data = v)
        loaded = Data.read_jld2(p)
        @test loaded == v
        @test _bitsame(loaded["x"], v["x"])   # BigFloat value+precision preserved
    end
end

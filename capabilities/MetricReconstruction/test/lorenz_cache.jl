using Test
import MetricReconstruction.Lorenz as L

@testset "Lorenz band-cache identity, context, and bounds" begin
    key = (:fixed,)
    cache = IdDict{Any,Dict{Any,Any}}()
    recs1 = Dict(:x => Ref(1))
    recs2 = Dict(:x => Ref(1))
    L._band_cache_put!(cache, recs1, key, 42)
    @test L._band_cache_get(cache, recs1, key) == 42
    @test L._band_cache_get(cache, recs2, key) === nothing

    recs = Dict(2 => Ref(1))
    r = BigFloat(8); a = BigFloat("0.5"); ω = BigFloat("0.1")
    k1 = setprecision(BigFloat, 128) do
        L._band_cache_key(recs, [2], r, 2, a, ω)
    end
    k2 = setprecision(BigFloat, 256) do
        L._band_cache_key(recs, [2], r, 2, a, ω)
    end
    k3 = setrounding(BigFloat, RoundDown) do
        L._band_cache_key(recs, [2], r, 2, a, ω)
    end
    @test k1 != k2
    @test k2 != k3

    kbefore = L._band_cache_key(recs, [2], r, 2, a, ω)
    recs[2] = Ref(1)
    kafter = L._band_cache_key(recs, [2], r, 2, a, ω)
    @test kbefore != kafter

    bounded = IdDict{Any,Dict{Any,Any}}()
    owner = Dict(:owner => true)
    for i in 1:(L._BAND_CACHE_LIMIT + 1)
        L._band_cache_put!(bounded, owner, (i,), i)
    end
    @test length(bounded[owner]) <= L._BAND_CACHE_LIMIT

    owners = IdDict{Any,Dict{Any,Any}}()
    for i in 1:(L._BAND_RECON_CACHE_LIMIT + 1)
        L._band_cache_put!(owners, Dict(:owner => i), key, i)
    end
    @test length(owners) <= L._BAND_RECON_CACHE_LIMIT

    threaded = IdDict{Any,Dict{Any,Any}}()
    Threads.@threads for i in 1:64
        L._band_cache_put!(threaded, owner, (i,), i)
    end
    @test length(threaded[owner]) == 64

    L._band_cache_put!(L._MIXBAND_CACHE, recs, key, 1)
    L._band_cache_put!(L._DIAGBAND_CACHE, recs, key, 2)
    L._invalidate_band_cache!(recs)
    @test !haskey(L._MIXBAND_CACHE, recs)
    @test !haskey(L._DIAGBAND_CACHE, recs)
end

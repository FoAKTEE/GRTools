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

# G16 lock-regression gate for the OTHER KerrComponents global — the r-independent angular-matrix
# memo `_ANGMAT_CACHE`. The wave-0 repro-audit listed it as non-thread-safe; it is NOT (the `get!`
# has always run under `_ANG_LOCK`), so this PINS the property rather than fixing it. The bundle is
# built once per (ls,m,aω,CT) and read from the threaded banded contraction, so a lost update hands
# two callers DIFFERENT matrix objects for one key. aω = 0 keeps every build on the b = I shortcut
# (no SWSH expansion), so 40 keys cost ~1 s serially.
# MUTATION-TESTED: with `lock(_ANG_LOCK)` deleted from src/Lorenz/KerrComponents.jl this testset
# died with a TaskFailedException out of `Dict.get!` (concurrent rehash) at -t4. Probabilistic, as
# every race gate is.
@testset "angular-matrix memo is race-free (G16 lock gate)" begin
    ls_sets = ([2,3], [2,3,4], [2,3,4,5], [3,4], [3,4,5], [3,4,5,6], [4,5], [4,5,6])
    keys40 = [(ls, m) for ls in ls_sets for m in -2:2]        # 40 distinct cache keys
    lock(L._ANG_LOCK) do
        empty!(L._ANGMAT_CACHE)
    end
    if Threads.nthreads() < 2
        @test_skip "concurrent memo gate needs -t2"   # single-threaded run: no interleaving to expose
    else
        jobs = [i for _ in 1:5 for i in eachindex(keys40)]
        got = Vector{Any}(undef, length(jobs))
        Threads.@threads for k in eachindex(jobs)
            ls, m = keys40[jobs[k]]
            got[k] = L._angular_mats_kerr(ls, m, 0.0, ComplexF64)
        end
        canonical = Dict{Int,Any}()
        for k in eachindex(jobs)
            canonical[jobs[k]] = get(canonical, jobs[k], got[k])
        end
        @test all(k -> got[k] === canonical[jobs[k]], eachindex(jobs))  # ONE bundle per key
        @test length(L._ANGMAT_CACHE) == length(keys40)                 # no insert lost
    end
    lock(L._ANG_LOCK) do
        empty!(L._ANGMAT_CACHE)                             # leave no state for later testsets
    end
end

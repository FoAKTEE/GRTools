# Lock-regression gate for the module-global GSN asymptotic-coefficient memo
# (src/Radial/gsn_asymptotics.jl `_GSN_COEFF_CACHE`).
#
# WHY THIS FILE EXISTS: the wave-0 repro-audit (G16) listed this cache as a non-thread-safe global.
# It is NOT — `_gsn_asym_coeff` has done its get/compute/store under `_GSN_COEFF_LOCK` since the
# GSN backend was written (GSN_PLAN §3c-2) — so the audit finding is REFUTED and there is nothing
# to fix. What was missing is a gate that KEEPS it true. This cache is the one with a
# read-modify-write body (`get` → maybe recompute a LONGER series → `setindex!`), i.e. exactly the
# shape where an unlocked Dict loses entries under a concurrent rehash, and the radial solvers that
# read it run inside threaded (l,m) sweeps.
#
# The gate is differential (concurrent coefficients BITWISE identical to the serial reference) plus
# structural (one entry per (kind,mode,type,precision) key — enough distinct keys that the Dict
# rehashes while other tasks insert — and every stored series grown to the longest order asked).
# MUTATION-TESTED: with `lock(_GSN_COEFF_LOCK)` deleted from src/Radial/gsn_asymptotics.jl, the
# entry-count assertion failed (191 of 192 keys survived, one insert lost) at -t4. Probabilistic,
# as every race gate is: it detects a removed lock with high probability, it does not prove
# race-freedom.

module TestCacheThreadSafety

using Test
import Perturbation: _gsn_coeff_in_inf, _gsn_coeff_out_inf, _gsn_coeff_in_hor, _gsn_coeff_out_hor,
                     _GSN_COEFF_CACHE, _GSN_COEFF_LOCK

# base mode = the test_radial_gsn oracle mode (a=0.7, ω=0.5, λ = s=-2,l=m=2 eigenvalue); the λ fan
# only has to produce DISTINCT cache keys, so it is a tiny offset ladder (the coefficients stay
# finite and the recurrence is the same code path).
const KINDS = (_gsn_coeff_in_inf, _gsn_coeff_out_inf, _gsn_coeff_in_hor, _gsn_coeff_out_hor)
const S, M, A, OM, LAM = -2, 2, 0.7, 0.5, 1.69660940163534
const LAMS = [LAM + 1e-6 * j for j in 0:47]        # 48 λ × 4 kinds = 192 distinct keys
const ORDERS = 0:2                                  # per-key orders exercised concurrently

_clear_gsn_cache!() = lock(_GSN_COEFF_LOCK) do
    empty!(_GSN_COEFF_CACHE)
end

@testset "GSN coefficient cache is race-free (G16 lock gate)" begin
    # serial reference on the base mode, cold cache, orders ascending (the series grows in place)
    _clear_gsn_cache!()
    ref = [f(S, M, A, OM, LAM, o, Float64) for f in KINDS, o in ORDERS]
    @test length(_GSN_COEFF_CACHE) == length(KINDS)          # one entry per kind at one mode

    if Threads.nthreads() < 2
        @test_skip "concurrent cache gate needs -t2"   # single-threaded run: no interleaving to expose
    else
        _clear_gsn_cache!()
        jobs = [(ki, li, o) for _ in 1:2 for o in ORDERS
                            for li in eachindex(LAMS) for ki in eachindex(KINDS)]
        got = Vector{ComplexF64}(undef, length(jobs))
        Threads.@threads for k in eachindex(jobs)
            ki, li, o = jobs[k]
            got[k] = KINDS[ki](S, M, A, OM, LAMS[li], o, Float64)
        end
        # (a) differential: the base-mode entries must be BITWISE the serial values
        @test all(k -> jobs[k][2] != 1 || got[k] === ref[jobs[k][1], jobs[k][3] + 1], eachindex(jobs))
        # (b) structural: no insert lost while the Dict rehashed under the other tasks
        @test length(_GSN_COEFF_CACHE) == length(KINDS) * length(LAMS)
        # (c) the read-modify-write growth path left every series at the longest order asked
        @test all(v -> length(v) == length(ORDERS), values(_GSN_COEFF_CACHE))
        _clear_gsn_cache!()                                   # leave no state for later testsets
    end
end

end # module TestCacheThreadSafety

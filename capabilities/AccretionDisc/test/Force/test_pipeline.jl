# Force.Pipeline — forced master-ODE solve + F8 torque-balance / flux-crosscheck / F5 gates.
# GROUNDING (Phase-7B: fixtures only; see test_source_fluid.jl header for fixture provenance):
#   • G-D2 F8 torque balance: ∂_r J̇^r_Adv == ∂T^φ_Mat — an ANALYTIC OPERATOR IDENTITY of the
#     master ODE (d/dr[p·Im(𝔥*𝔥')] = p·Im(𝔥*𝔉), exact for real A_m/B_m); asserted at the legacy
#     frozen gate mean_rel < 1e-6 (capture-run oracle on these fixtures: 9.882e-08).
#   • G-D3 flux cross-check: the wave-part reconstruction tracks the wave-action flux (CV smaller
#     than the forcing-dominated full CV) — the legacy relational gate (capture: 0.188 < 0.747).
#   • F5: |𝔉| peaks within 1M of the OLR (legacy gate; capture: |Δ|=0.103).
#   • ForcingInterp/forced_solve unit coverage: the legacy cov_pipeline closed-form assertions.
#   • §8 parallel: threaded build_forcing == serial source_F loop BIT-IDENTICAL (race-free
#     per-slot + blas_single_thread), plus the legacy speedup report. The speedup itself is
#     wall-clock and therefore machine-load-dependent: it is REPORTED always and ASSERTED
#     (>1, the legacy CI-variance-safe gate) only under Threads.nthreads()>1 AND
#     GRTOOLS_PERF_TESTS=1; otherwise it registers as a skipped (Broken) test.
using Test
using AccretionDisc
using AccretionDisc.Force
using Printf

const FIXDIR_PL = joinpath(@__DIR__, "..", "fixtures")

@testset "Force.Pipeline (F8/F5 gates + forced solve, a=0 fixtures)" begin
    d = disc_a0(); m = 2; r_p = 8.0
    ω0 = m * disc_OmegaK(d, r_p)
    modes = load_modes(2:6, m; dir = FIXDIR_PL)
    g = modes[m].r_right

    @testset "cov: ForcingInterp call operator" begin
        rg = 11.0:1.0:14.0
        F = build_forcing(d, modes, m, r_p, ω0, rg; lmax = 6)
        @test abs(F(F.rg[2]) - F.Fg[2]) < 1e-12                       # node value
        mid = (F.rg[2] + F.rg[3]) / 2
        @test abs(F(mid) - (F.Fg[2] + F.Fg[3]) / 2) < 1e-12           # linear ⇒ exact midpoint average
        @test F(F.rg[1] - 5) == F.Fg[1]                               # endpoint clamps
        @test F(F.rg[end] + 5) == F.Fg[end]
    end

    @testset "cov: forced_solve (RK4, zero-IC source-driven)" begin
        rg = 11.0:1.0:14.0
        F = build_forcing(d, modes, m, r_p, ω0, rg; lmax = 6)
        N = 2000
        rs, hs, hps = forced_solve(d, m, r_p, F, 14.0, 12.0, N)
        @test length(rs) == N + 1 && length(hs) == N + 1 && length(hps) == N + 1
        @test isapprox(rs[1], 14.0; atol = 1e-9) && isapprox(rs[end], 12.0; atol = 1e-9)
        @test hs[1] == 0 && hps[1] == 0                               # no incoming wave
        @test abs(hs[end]) > 0 && all(isfinite, hs) && all(isfinite, hps)
    end

    @testset "G-D2 F8 torque balance (~1e-7 closure)" begin
        res = f8_torque_balance(d, m, r_p, modes, 6; N = 16000, verbose = false)
        @printf("  G-D2 F8: mean_rel=%.3e max_rel=%.3e mean_int=%.3e (capture oracle 9.882e-08)\n",
                res.mean_rel, res.max_rel, res.mean_rel_integral)
        @test res.mean_rel < 1e-6                                     # legacy frozen G-D2 gate
        @test res.npts > 0
    end

    @testset "G-D3 F8 flux cross-check" begin
        res = f8_flux_crosscheck(d, m, r_p, modes, 6; N = 12000, verbose = false)
        @printf("  G-D3 flux: cv_full=%.3e cv_wave=%.3e\n", res.cv_full, res.cv_wave)
        @test res.cv_wave < res.cv_full                               # wave-part tracks the conserved flux
        @test isfinite(res.cv_wave)
    end

    @testset "F5 source structure (|𝔉| peaks at the OLR)" begin
        res = f5_source(d, m, r_p, modes, 6; rg = 8.2:0.25:15.0, verbose = false)
        @printf("  F5: peak r=%.3f OLR=%.4f |Δ|=%.3f\n",
                Float64(res.r_peak), Float64(res.OLR), abs(Float64(res.r_peak) - Float64(res.OLR)))
        @test abs(Float64(res.r_peak) - Float64(res.OLR)) < 1.0       # legacy frozen F5 gate
    end

    @testset "§8 parallel: threaded build_forcing == serial (bit-identical) + speedup" begin
        rg = 9.0:0.05:15.0                                            # 121 pts — amortizes thread overhead
        # warm up both code paths (exclude compilation from the timing)
        build_forcing(d, modes, m, r_p, ω0, 9.0:1.0:15.0; lmax = 6)
        let s = ComplexF64[]
            for rr in 9.0:1.0:15.0
                kk = argmin(abs.(g .- rr)); push!(s, source_F(d, modes, kk, m, r_p, ω0; lmax = 6))
            end
        end
        serial = ComplexF64[]
        t_ser = @elapsed begin
            for rr in rg
                kk = argmin(abs.(g .- rr)); push!(serial, source_F(d, modes, kk, m, r_p, ω0; lmax = 6))
            end
        end
        local F
        t_par = @elapsed (F = build_forcing(d, modes, m, r_p, ω0, rg; lmax = 6))
        @printf("  parallel (nthreads=%d, %d pts): t_ser=%.3fs t_par=%.3fs speedup=%.2fx\n",
                Threads.nthreads(), length(serial), t_ser, t_par, t_ser / t_par)
        @test F.Fg == serial                                          # BIT-IDENTICAL (§8 gate)
        speedup = t_ser / t_par
        @test isfinite(speedup) && speedup > 0
        # The wall-clock speedup is LOAD-SENSITIVE — on a busy shared machine the threaded
        # run can be the slower one — so it is reported by default and asserted only when
        # GRTOOLS_PERF_TESTS=1 declares the machine quiet enough to gate on timing.
        if Threads.nthreads() > 1 && get(ENV, "GRTOOLS_PERF_TESTS", "0") == "1"
            @test speedup > 1.0                                       # legacy CI-variance-safe gate
        else
            @printf("  speedup assertion SKIPPED (measured %.2fx, nthreads=%d); set GRTOOLS_PERF_TESTS=1 to enforce\n",
                    speedup, Threads.nthreads())
            @test_skip speedup > 1.0
        end
    end
end

# Characterization of HybridC module 3 (wasabic_inspiral): the WaSABI-C 9-var forcing assembly,
# evolution equations, stop frequency and the δx-coupled trajectory driver/solution.
#
# GROUNDING (WaSABI NOT a live dep — baked WL-oracle literals from the legacy oracle16b run,
# `wolframscript -file oracle16b.wl` over the upstream paclet; reproduced verbatim from
# packages/WaSABI/test/inspiral_WaSABI-C_test.jl):
#   * per-channel `ForcingTerms["WaSABI-C"]` vector (Fω,Fφ,Fν,Fm,Fχ1,Fχ2,Fδm,Fδν,Fδχ) (rtol 1e-9,
#     atol 1e-24) — the "FCH" blocks.
#   * `BinaryInspiral["WaSABI-C"]` 9-var trajectory states + stop time at sampled t (rtol 1e-6,
#     shifts rtol 1e-4) — the "BIM" blocks.
#   * serial == threaded `int_inspiral` parity (R3): pure-function → bit-identical, before speedup.

using Test
using Hybrid.HybridC
using Hybrid.Drivers: int_inspiral, duration

# Baked fixture copy of the upstream WaSABI data (provenance: test/fixtures/wasabi/README.md).
const _WI_BASE = joinpath(@__DIR__, "fixtures", "wasabi")

@testset "HybridC.wasabic_inspiral" begin
    mdl = wasabic_model(joinpath(_WI_BASE, "InspiralModels", "sf_data"),
                        joinpath(_WI_BASE, "AmplitudeModels", "sf_amp_data"))

    @testset "forcing channel vector vs WL ForcingTerms (per channel)" begin
        cases = [
            ((0.02, 0.0, 0.1, 1, 0.5, 0.3, 0.1, 0.05, 0.2),
             (1.1270202541171997e-6, 0.02, 2.681105260055963e-12, -3.096924505082267e-11,
              -1.942598793426137e-9, -1.1385680120862304e-10, -3.1237355576828266e-10,
              2.654294207455403e-10, -1.9431350144781482e-8)),
            ((0.03, 0.0, 0.08, 1, -0.4, 0.6, -0.1, 0.1, -0.3),
             (6.139677856352448e-6, 0.03, -2.8367233238650157e-11, 3.898391289696152e-10,
              1.600167770800929e-8, -9.46215866464914e-10, 4.8375300705718776e-9,
              -4.3614621104424614e-9, 1.9991459422547118e-7)),
        ]
        for (args, wl) in cases
            F = wc_forcing(mdl, args...)
            for i in 1:9
                @test isapprox(F[i], wl[i]; rtol = 1e-9, atol = 1e-24)
            end
        end
    end

    @testset "BinaryInspiral 9-var trajectory + stop time vs WL" begin
        runs = [
            ((; ω = 0.02, ν = 0.1, χ1 = 0.5, χ2 = 0.3), 6173.477487710508,
             [(2057.825829236836, (0.022960044450056023, 43.94133361266952, 0.1000000077645276,
               0.9999999103110571, 0.4999947692437082, 0.29999969262326326,
               -8.968894101765542e-7, 7.764526338963532e-7, -0.00005230755885355894)),
              (4115.651658473672, (0.028752400638931152, 96.31477589672535, 0.10000002687713096,
               0.9999996895321818, 0.49998411183763375, 0.2999990627119567,
               -3.1046778399404796e-6, 2.687711649461948e-6, -0.0001588815809552501))]),
            ((; ω = 0.015, ν = 0.08, χ1 = -0.3, χ2 = 0.4), 13590.819820625418,
             [(4530.273273541806, (0.017049472177817535, 72.21092752472404, 0.07999999793410834,
               1.0000000283403085, -0.2999978515628354, 0.3999998356148569,
               3.542538644192545e-7, -3.2279558777226565e-7, 2.6855465247977703e-5)),
              (9060.546547083612, (0.02098665626863837, 157.14809611358552, 0.0799999925224834,
               1.0000001026493415, -0.2999931066365324, 0.39999951550616325,
               1.28311681996549e-6, -1.1683621838926618e-6, 8.616705139644089e-5))]),
        ]
        for (ics, dur_wl, samples) in runs
            s = int_inspiral(mdl, ics; reltol = 1e-11, abstol = 1e-11)
            @test isapprox(duration(s), dur_wl; rtol = 1e-6)
            for (t, st) in samples
                u = wc_state_at(s, t)
                @test isapprox(u.ω, st[1]; rtol = 1e-6)
                @test isapprox(u.φ, st[2]; rtol = 1e-6)
                @test isapprox(u.ν, st[3]; rtol = 1e-6)
                @test isapprox(u.M, st[4]; rtol = 1e-6)
                @test isapprox(u.χ1, st[5]; rtol = 1e-6, atol = 1e-8)
                @test isapprox(u.χ2, st[6]; rtol = 1e-6)
                @test isapprox(u.δm, st[7]; rtol = 1e-4, atol = 1e-10)
                @test isapprox(u.δν, st[8]; rtol = 1e-4, atol = 1e-10)
                @test isapprox(u.δχ, st[9]; rtol = 1e-4, atol = 1e-9)
            end
        end
    end

    @testset "serial == threaded int_inspiral parity (R3)" begin
        icslist = [(; ω = 0.02, ν = 0.1, χ1 = 0.5, χ2 = 0.3),
                   (; ω = 0.015, ν = 0.08, χ1 = -0.3, χ2 = 0.4),
                   (; ω = 0.025, ν = 0.12, χ1 = 0.7, χ2 = -0.2)]
        serial = [duration(int_inspiral(mdl, ic)) for ic in icslist]
        threaded = Vector{Float64}(undef, length(icslist))
        Threads.@threads for i in eachindex(icslist)
            threaded[i] = duration(int_inspiral(mdl, icslist[i]))
        end
        for i in eachindex(icslist)
            @test threaded[i] == serial[i]
        end
    end
end

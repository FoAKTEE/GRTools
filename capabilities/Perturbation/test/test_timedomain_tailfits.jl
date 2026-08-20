# DG-P4 tail-rate estimators + the RT-chart excision policy — gates.
#
# These five functions were figure-local in the legacy arXiv:1301.1591 reproduction drivers
# (emri-disk `results/dg-teukolsky-spin/scripts/{driver_common,fit_rates}.jl`); they are promoted
# here, next to `fit_tail_window`, because each one PRODUCES A PUBLISHED NUMBER (the measured tail
# exponents of Campaign A/B and the self-convergence rates of Campaign C).
#
# GATES:
#   (T1) fit_tail_envelope     : exact exponent on a synthetic oscillating power law whose deep
#                                nulls bias the plain log-log fit (the diagnosed |s|=2 κ-oscillation
#                                failure mode); NaN when fewer than four envelope peaks.
#   (T2) fit_tail_adaptive     : exact exponent on a power law that turns up into a noise floor
#                                (the estimator must stop AT the turn-up, not fit through it).
#   (T3) fit_tail_sliding      : plateau/min/max on a clean power law collapse onto the exponent.
#   (T4) fit_convergence_order : exact order on a synthetic E ∝ K^-p ladder; REUSES the single
#                                least-squares log-log kernel `fit_tail_window` (no second copy).
#   (T5) spin_excise_policy    : cap branch (a=0, a=0.9), horizon-gap branch (near-extremal), and
#                                agreement with the Spacetime RT chart maps.
#   (L)  legacy characterization: the promoted estimators reproduce the measured numbers of the
#                                legacy DG campaign to the printed precision — the self-convergence
#                                rates of `data/verifier_conv.csv` recomputed from the tracked
#                                `data/conv_s-2_a0p9_N*_tf50.csv` ladders (baked inline below).

using Test
using Perturbation
const TDS = Perturbation.TimeDomain

@testset "TimeDomain tail-rate estimators + excise policy" begin

    # ---- (T1) envelope fit through the local maxima ---------------------------------------------
    @testset "gate T1: fit_tail_envelope on an oscillating power law" begin
        τ = collect(600.0:1.0:1400.0)
        # |ψ| = A τ^-7 (½+½cos(2πτ/74))² — an EXACT τ^-7 envelope with deep nulls between peaks,
        # the shape the diagnosed κ-mediated |s|=2 oscillation imposes on the RT-chart tails.
        ψ = @. 1e18 * τ^(-7.0) * (0.5 + 0.5 * cos(2π * τ / 74))^2
        n_env, npk = TDS.fit_tail_envelope(τ, ψ, (700.0, 1300.0))
        @test npk == 8
        @test isapprox(n_env, -7.0; atol = 1e-9)              # envelope recovers the exponent exactly
        # the plain window fit is measurably pulled off by the nulls — that is why this estimator exists
        n_raw = TDS.fit_tail_window(τ, ψ, (700.0, 1300.0))
        @test abs(n_raw - (-7.0)) > 10 * abs(n_env - (-7.0))

        # a MONOTONE record has no interior maxima ⇒ NaN, never a fabricated number
        n_mono, npk_mono = TDS.fit_tail_envelope(τ, 1e18 .* τ .^ (-7.0), (700.0, 1300.0))
        @test isnan(n_mono)
        @test npk_mono == 0
        # too few samples in the window ⇒ NaN
        τshort = collect(700.0:1.0:704.0)
        n2, npk2 = TDS.fit_tail_envelope(τshort, 1e18 .* τshort .^ (-7.0), (700.0, 704.0))
        @test isnan(n2)
        @test npk2 < 4
    end

    # ---- (T2) adaptive window stops at the noise-floor turn-up ----------------------------------
    @testset "gate T2: fit_tail_adaptive stops at the floor turn-up" begin
        τ = collect(800.0:1.0:3000.0)
        pure = τ .^ (-6.0)                    # crosses 1e-19 near τ ≈ 1479

        # (a) flat junk floor: the caller's `floor` ends the window at the crossover
        n_floor = TDS.fit_tail_adaptive(τ, max.(pure, 1e-19); t0 = 800.0, floor = 2e-19)
        @test isapprox(n_floor, -6.0; atol = 1e-6)

        # (b) turn-up: the record jumps ×10 at τ = 1500 and then rises — the 1.5× rule stops there
        jump = copy(pure)
        jump[701:end] .= 10 .* pure[701] .* (1 .+ 0.001 .* (0:(length(pure) - 701)))
        n_jump = TDS.fit_tail_adaptive(τ, jump; t0 = 800.0, floor = 1e-30)
        @test isapprox(n_jump, -6.0; atol = 1e-6)

        # fitting the WHOLE record through the floor is measurably shallower
        n_all = TDS.fit_tail_window(τ, max.(pure, 1e-19), (800.0, 3000.0))
        @test n_all > n_floor + 1.0
        # fewer than four samples above t0 ⇒ NaN
        @test isnan(TDS.fit_tail_adaptive(τ, pure; t0 = 2999.0, floor = 1e-30))
    end

    # ---- (T3) sliding-window plateau -------------------------------------------------------------
    @testset "gate T3: fit_tail_sliding plateau on a clean power law" begin
        τ = collect(500.0:1.0:3200.0)
        ψ = 1e22 .* τ .^ (-8.0)      # stays above the 1e-13 floor over the whole scan
        plateau, nmin, nmax, nwin = TDS.fit_tail_sliding(τ, ψ)
        @test nwin >= 4
        @test isapprox(plateau, -8.0; atol = 1e-9)
        @test isapprox(nmin, -8.0; atol = 1e-9)
        @test isapprox(nmax, -8.0; atol = 1e-9)
        # a series that never clears the floor gives NaN, not a number
        p2, _, _, nw2 = TDS.fit_tail_sliding(τ, fill(1e-20, length(τ)))
        @test isnan(p2)
        @test nw2 == 0
    end

    # ---- (T4) convergence order -------------------------------------------------------------------
    @testset "gate T4: fit_convergence_order exact on a synthetic ladder" begin
        Ks = [8.0, 12.0, 16.0, 24.0, 32.0, 48.0]
        for p in (3.0, 5.0, 9.0)
            E = 2.5 .* Ks .^ (-p)
            @test isapprox(TDS.fit_convergence_order(Ks, E), p; atol = 1e-10)
        end
        # ordinate window: points outside (lo, hi) are dropped
        E = 2.5 .* Ks .^ (-5.0); E[end] = 1e-16          # below the default floor
        @test isapprox(TDS.fit_convergence_order(Ks, E), 5.0; atol = 1e-10)
        @test isnan(TDS.fit_convergence_order(Ks, fill(1e-18, length(Ks))))
    end

    # ---- (T5) excision policy ----------------------------------------------------------------------
    @testset "gate T5: spin_excise_policy branches + chart agreement" begin
        for a in (0.0, 0.5, 0.9)
            pol = TDS.spin_excise_policy(1.0, a)
            @test pol.Rplus == TDS.rt_Rplus(1.0, a)
            @test pol.excise == 0.05                       # cap branch: 0.7(R₊−R₋) > 0.05
            @test pol.Rplus - pol.excise > pol.Rminus      # inner edge stays inside the trapped region
        end
        @test TDS.spin_excise_policy(1.0, 0.0).Rminus == 0.0          # a=0: no inner horizon
        # near-extremal: the horizon gap closes and the policy shrinks the excision
        pol = TDS.spin_excise_policy(1.0, 0.9995)
        @test pol.excise < 0.05
        @test isapprox(pol.excise, 0.7 * (pol.Rplus - pol.Rminus); rtol = 1e-14)
        @test pol.Rplus - pol.excise > pol.Rminus
        # R₋ is the RT image of the inner BL horizon
        rm = 1.0 - sqrt(1.0 - 0.9995^2)
        @test isapprox(pol.Rminus, TDS.rt_R_of_r(rm); rtol = 1e-14)
    end

    # ---- (L) legacy characterization: the published DG campaign-C rates --------------------------
    @testset "gate L: legacy verifier_conv rates recomputed" begin
        # tracked legacy ladders emri-disk results/dg-teukolsky-spin/data/conv_s-2_a0p9_N*_tf50.csv
        # (rel-L2-in-time vs the finest-K reference, ℓ=2 probe, a=0.9, s=−2, filter off, tf=50)
        ladders = Dict(
            2 => ([16.0, 24, 32, 48, 64, 96, 128],
                  [1.42984298e-1, 7.54665297e-2, 3.29354052e-2, 5.16751446e-3, 1.06309504e-3,
                   3.28481619e-4, 1.04684970e-4],                              # hor
                  [1.68049869e-1, 6.99253614e-2, 2.22030723e-2, 2.57247691e-3, 2.83792347e-4,
                   9.14951343e-5, 3.72330871e-5],                              # mid
                  [1.13953239e-1, 2.04558543e-2, 2.96629962e-3, 8.50945989e-5, 1.13459990e-5,
                   1.96350844e-6, 6.49196614e-7]),                             # scri
            4 => ([8.0, 12, 16, 24, 32, 48, 64],
                  [7.95390716e-3, 2.03203891e-3, 5.45985096e-4, 7.86053965e-5, 1.22732977e-5,
                   2.71428583e-6, 4.90222110e-7],
                  [7.95577715e-3, 1.39154442e-3, 4.47859528e-4, 2.95849466e-5, 3.17401306e-6,
                   9.63644266e-7, 3.10445151e-7],
                  [7.09766481e-3, 7.50898879e-4, 8.63496223e-5, 2.26132852e-6, 1.19173933e-6,
                   5.48485229e-8, 1.83936433e-9]),
            6 => ([6.0, 8, 12, 16, 24, 32],
                  [1.86408383e-2, 4.12482991e-4, 5.26712801e-6, 1.66298110e-6, 8.63157214e-8,
                   1.88618841e-8],
                  [1.37064024e-2, 3.48721478e-4, 4.35634358e-6, 7.36736943e-7, 2.24945364e-8,
                   9.42142588e-9],
                  [1.06079691e-2, 7.63309359e-5, 1.51349797e-6, 6.14783778e-8, 2.46165340e-9,
                   1.61974960e-10]))
        # legacy data/verifier_conv.csv, printed to 2 decimals
        published = Dict((2, :hor) => 3.70, (2, :mid) => 4.41, (2, :scri) => 6.21,
                         (4, :hor) => 4.74, (4, :mid) => 5.15, (4, :scri) => 7.09,
                         (6, :hor) => 8.00, (6, :mid) => 8.52, (6, :scri) => 10.36)
        for (N, (Ks, Eh, Em, Es)) in ladders
            for (loc, E) in ((:hor, Eh), (:mid, Em), (:scri, Es))
                @test isapprox(TDS.fit_convergence_order(Ks, E), published[(N, loc)]; atol = 5e-3)
            end
        end
    end
end

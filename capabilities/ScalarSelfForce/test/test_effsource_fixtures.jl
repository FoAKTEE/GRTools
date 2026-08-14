# test_effsource_fixtures.jl — ESSF node `effsource-oracle-fixtures` (schedule W2-D), written
# BEFORE the fixtures exist (TDD). Fixture-INTEGRITY gates only — this node ships test assets
# + capture tooling, no runtime code. Gates (GLOBAL_DAG #5 oracle column + node contract):
#   * fixture files load and have the FROZEN grid shape/counts:
#       configs (a, r0) ∈ {(0,10), (0.5,10), (0.9,isco+0.1), (0.998,isco), (-0.5,10),
#       (-0.998,isco)}, m ∈ {0,1,2,5,10,20} (comoving PhiS_m additionally {-2,-1} for the
#       conjugate-symmetry gate), per-config point counts npts = [40,40,38,37,40,40]
#       (near-worldline |Δr| down to 1e-6 M through far field, minus-side offsets clipped
#       above the horizon for the two near-ISCO configs);
#   * every config has points on BOTH sides of the comoving variant's z = 0.995
#     hypergeometric switch (kerr-circular-comoving.c:27, argument z = (zL−1)/(zL+1),
#     zL = (ϱ² + z_c²/2)/(ϱ·√(ϱ²+z_c²)));
#   * baked orbit constants are EXACTLY the solid Spacetime values (transcription gate:
#     kerr_isco / constants_of_motion / horizons round-trip, Float64 ==);
#   * NaN sentinel structure: exactly the 5 mixed second derivatives ∂tr, ∂tθ, ∂rθ, ∂rφ,
#     ∂θφ (10 of the 20 d² slots) are NaN in every calc_m row — the C code does not provide
#     them (kerr-circular.c:709-713, kerr-circular-comoving.c:570-574); ALL other slots
#     finite (measured: no other nonfinite value anywhere in the capture);
#   * documented invariants (all bounds MEASURED on the committed capture — see the
#     effsource_configs.jl header): Ψ_{-m} = conj(Ψ_m) (comoving, m = ±1, ±2; ≤ 4.1e-15);
#     e^{−i m φ_p} phase-factor structure of PhiS_m/calc_m pairs (real- vs
#     imaginary-rotated slots, residuals ≤ 1.3e-16); the calc_m_series `+sin(m φ_p)`
#     OPPOSITE phase convention (the phase bridge, kerr-circular-comoving.c:673-674 vs
#     :634-635); cross-variant PhiS_m agreement (elliptic E/K vs comoving 3P Legendre,
#     ≤ 2.8e-10) inside the overlap window C1 = α/β ∈ [1e-8, 1e-2) — outside it the two
#     paths are complementary oracles (elliptic loses digits near the worldline AND
#     cancels catastrophically in the far field for m ≥ 2 — the W4-D BigFloat(256)
#     adjudication; comoving/3P tracks truth to ~1e-11 over the whole frozen grid,
#     see the effsource_configs.jl header); series↔direct src handoff
#     (min-over-points |src| relative difference ≤ 6.0e-4 for every (config, m));
#   * re-run reproducibility: recompiling + rerunning the capture driver reproduces the
#     COMMITTED fixture bytes identically (run twice, diff) — executed when gcc + GSL +
#     make are available (this machine), @test_skip otherwise.
using Test
using Spacetime: kerr_isco, constants_of_motion, horizons

const EFFSOURCE_DIR = joinpath(@__DIR__, "fixtures", "effsource")

include(joinpath(EFFSOURCE_DIR, "effsource_configs.jl"))
include(joinpath(EFFSOURCE_DIR, "phis4d_elliptic.jl"))
include(joinpath(EFFSOURCE_DIR, "phis4d_comoving.jl"))
include(joinpath(EFFSOURCE_DIR, "phis_m_elliptic.jl"))
include(joinpath(EFFSOURCE_DIR, "phis_m_comoving.jl"))
include(joinpath(EFFSOURCE_DIR, "calc_m_elliptic.jl"))
include(joinpath(EFFSOURCE_DIR, "calc_m_comoving.jl"))
include(joinpath(EFFSOURCE_DIR, "calc_m_series_comoving.jl"))

# ---- frozen-grid constants (must match the generator; duplicated here on purpose so the
# ---- test pins the capture, not the other way round) -----------------------------------
const EXPECTED_NPTS = [40, 40, 38, 37, 40, 40]
const EXPECTED_MS = [0, 1, 2, 5, 10, 20]
const EXPECTED_MS_PHIS_M_COMOVING = [-2, -1, 0, 1, 2, 5, 10, 20]
const EXPECTED_DPHI = [0.0, 0.5, 2.0]

# calc_m fixture column map (37 cols; see calc_m_*.jl headers)
const C_CFG, C_M, C_IPT, C_DR, C_DTH = 1, 2, 3, 4, 5
const C_PHIS = 6                      # (re, im)
const C_D1 = 8                        # 4 pairs: t, r, θ, φ
const C_D2 = 16                       # 10 pairs: tt, tr, tθ, tφ, rr, rθ, rφ, θθ, θφ, φφ
const C_SRC = 36                      # (re, im)
# d² pair indices (1-based within the 10 pairs) the C code does NOT provide:
const NAN_D2_PAIRS = [2, 3, 6, 7, 9]  # tr, tθ, rθ, rφ, θφ
const REAL_D2_PAIRS = [1, 4, 5, 8, 10] # tt, tφ, rr, θθ, φφ  (real-rotated)

# cross-variant overlap window in C1 = α/β where BOTH m-mode paths are well conditioned
# (fixture header): below 1e-8 the elliptic E/K path degrades (k' = √C1 cancellation in
# Kcomp), above ~1e-2 the comoving Legendre combination starts cancelling (total loss by
# C1 ~ 0.1 at high m). Measured worst inside the window: 2.8e-10.
const XVAR_C1_MIN, XVAR_C1_MAX = 1e-8, 1e-2
const XVAR_GATE = 3e-9
# series↔direct src handoff: measured worst-of-min 6.0e-4 (m = 20); gate with margin.
const HANDOFF_GATE = 5e-3

cfgrow(c) = EFFSOURCE_CONFIGS[c]
# field k across all fixture records (records are stored as matrix COLUMNS)
col(M, k) = view(M, k, :)

# zL/zarg of the comoving 0.995 hypergeometric switch, from baked per-config constants
function switch_zarg(c, dr, dth)
    fx = cfgrow(c)
    alpha = fx.alpha20 * dr^2 + fx.alpha02 * dth^2
    rho, zc = sqrt(alpha), sqrt(fx.beta)
    zL = (rho^2 + zc^2 / 2) / (rho * sqrt(rho^2 + zc^2))
    (zL - 1) / (zL + 1)
end

c1_of(c, dr, dth) = (fx = cfgrow(c); (fx.alpha20 * dr^2 + fx.alpha02 * dth^2) / fx.beta)

# phase-rotation residuals: z·e^{+i m φ_p} must be purely real (real-rotated slots) or
# purely imaginary (∂t, ∂φ slots) — exact trig structure of the C code, roundoff only.
function rotated(re, im, m)
    w = complex(re, im) * cis(m * EFFSOURCE_PHIP)
    w, abs(w)
end

@testset "EffectiveSource oracle fixtures (effsource-oracle-fixtures)" begin

    @testset "frozen grid: configs, points, shapes, counts" begin
        @test length(EFFSOURCE_CONFIGS) == 6
        @test EFFSOURCE_PHIP == 0.31
        @test size(EFFSOURCE_POINTS, 1) == 4
        @test size(EFFSOURCE_POINTS, 2) == sum(EXPECTED_NPTS)
        for c in 1:6
            fx = cfgrow(c)
            @test fx.idx == c
            @test fx.npts == EXPECTED_NPTS[c]
            pts = EFFSOURCE_POINTS[:, EFFSOURCE_POINTS[1, :] .== c]
            @test size(pts, 2) == fx.npts
            @test pts[2, :] == collect(1.0:fx.npts)         # ipt contiguous
            # |Δr| reaches down to 1e-6 M and out to the far field
            drs = pts[3, :]
            @test minimum(abs.(drs[drs .!= 0])) == 1e-6
            @test maximum(drs) == 10.0
            # all field points stay outside the horizon
            @test all(fx.r0 .+ drs .> fx.rp)
            # both sides of the 0.995 hypergeometric switch are populated
            zargs = [switch_zarg(c, pts[3, k], pts[4, k]) for k in 1:size(pts, 2)]
            @test any(zargs .> 0.995) && any(zargs .< 0.995)
        end
    end

    @testset "orbit constants = solid Spacetime values (transcription gate)" begin
        # frozen (a, r0-rule) pairs of the node contract
        rules = [(0.0, () -> 10.0), (0.5, () -> 10.0),
                 (0.9, () -> kerr_isco(0.9) + 0.1), (0.998, () -> kerr_isco(0.998)),
                 (-0.5, () -> 10.0), (-0.998, () -> kerr_isco(-0.998))]
        for c in 1:6
            fx = cfgrow(c)
            a, r0f = rules[c]
            @test fx.a == a
            @test fx.r0 == r0f()
            com = constants_of_motion(a, fx.r0, 0.0, 1.0)
            @test fx.E == com.E
            @test fx.Lz == com.Lz
            @test fx.risco == kerr_isco(a)
            @test fx.rp == horizons(a)[1]
            # puncture-scale constants (main.tex:445-455, M = 1)
            @test isapprox(fx.alpha20, fx.r0^2 / (fx.a^2 + fx.r0 * (fx.r0 - 2)); rtol = 4eps())
            @test fx.alpha02 == fx.r0^2
            @test isapprox(fx.beta, 4 * (fx.Lz^2 + fx.r0^2 + fx.a^2 * (fx.r0 + 2) / fx.r0);
                           rtol = 4eps())
        end
    end

    @testset "row shapes and metadata consistency" begin
        npts_tot = sum(EXPECTED_NPTS)
        for M4 in (EFFSOURCE_PHIS4D_ELLIPTIC, EFFSOURCE_PHIS4D_COMOVING)
            @test size(M4) == (6, 3 * npts_tot)
            @test sort(unique(col(M4, 5))) == EXPECTED_DPHI
        end
        @test size(EFFSOURCE_PHIS_M_ELLIPTIC) == (7, length(EXPECTED_MS) * npts_tot)
        @test size(EFFSOURCE_PHIS_M_COMOVING) ==
              (7, length(EXPECTED_MS_PHIS_M_COMOVING) * npts_tot)
        @test sort(unique(col(EFFSOURCE_PHIS_M_ELLIPTIC, 2))) == Float64.(EXPECTED_MS)
        @test sort(unique(col(EFFSOURCE_PHIS_M_COMOVING, 2))) ==
              Float64.(EXPECTED_MS_PHIS_M_COMOVING)
        for MC in (EFFSOURCE_CALC_M_ELLIPTIC, EFFSOURCE_CALC_M_COMOVING)
            @test size(MC) == (37, length(EXPECTED_MS) * npts_tot)
            @test sort(unique(col(MC, C_M))) == Float64.(EXPECTED_MS)
        end
        @test size(EFFSOURCE_CALC_M_SERIES_COMOVING) ==
              (7, length(EXPECTED_MS) * npts_tot)
        # (dr, dth) columns agree with the baked point table everywhere
        ptkey = Dict((EFFSOURCE_POINTS[1, k], EFFSOURCE_POINTS[2, k]) =>
                     (EFFSOURCE_POINTS[3, k], EFFSOURCE_POINTS[4, k])
                     for k in 1:size(EFFSOURCE_POINTS, 2))
        for M in (EFFSOURCE_PHIS_M_ELLIPTIC, EFFSOURCE_PHIS_M_COMOVING,
                  EFFSOURCE_CALC_M_ELLIPTIC, EFFSOURCE_CALC_M_COMOVING,
                  EFFSOURCE_CALC_M_SERIES_COMOVING)
            ok = all(ptkey[(M[1, k], M[3, k])] == (M[4, k], M[5, k])
                     for k in 1:size(M, 2))
            @test ok
        end
    end

    @testset "NaN sentinel structure (calc_m mixed d² slots)" begin
        for MC in (EFFSOURCE_CALC_M_ELLIPTIC, EFFSOURCE_CALC_M_COMOVING)
            nan_ok = true; fin_ok = true
            for k in 1:size(MC, 2)
                for p in NAN_D2_PAIRS
                    nan_ok &= isnan(MC[C_D2 + 2p - 2, k]) && isnan(MC[C_D2 + 2p - 1, k])
                end
                # every non-sentinel slot is finite (measured: holds for the whole capture,
                # including the ill-conditioned regimes — digit LOSS, never overflow)
                for j in [C_PHIS, C_PHIS + 1, C_D1:(C_D1 + 7)...,
                          (C_D2 + 2p - 2 for p in REAL_D2_PAIRS)...,
                          (C_D2 + 2p - 1 for p in REAL_D2_PAIRS)..., C_SRC, C_SRC + 1]
                    fin_ok &= isfinite(MC[j, k])
                end
            end
            @test nan_ok
            @test fin_ok
        end
        for M in (EFFSOURCE_PHIS_M_ELLIPTIC, EFFSOURCE_PHIS_M_COMOVING,
                  EFFSOURCE_PHIS4D_ELLIPTIC, EFFSOURCE_PHIS4D_COMOVING,
                  EFFSOURCE_CALC_M_SERIES_COMOVING)
            @test all(isfinite, M[6:end, :])
        end
    end

    @testset "e^{-i m φ_p} phase structure (PhiS_m, calc_m slots)" begin
        # PhiS_m: real-rotated
        for M in (EFFSOURCE_PHIS_M_ELLIPTIC, EFFSOURCE_PHIS_M_COMOVING)
            worst = 0.0
            for k in 1:size(M, 2)
                isfinite(M[6, k]) || continue
                w, aw = rotated(M[6, k], M[7, k], M[2, k])
                aw == 0 && continue
                worst = max(worst, abs(imag(w)) / aw)
            end
            @test worst < 1e-12
        end
        # calc_m: ∂t, ∂φ imaginary-rotated; PhiS, ∂r, ∂θ, non-NaN d², src real-rotated
        for MC in (EFFSOURCE_CALC_M_ELLIPTIC, EFFSOURCE_CALC_M_COMOVING)
            worst_re = 0.0; worst_im = 0.0
            for k in 1:size(MC, 2)
                m = MC[C_M, k]
                isfinite(MC[C_PHIS, k]) || continue
                for base in [C_PHIS, C_D1 + 2, C_D1 + 4,
                             (C_D2 + 2p - 2 for p in REAL_D2_PAIRS)..., C_SRC]
                    w, aw = rotated(MC[base, k], MC[base + 1, k], m)
                    aw == 0 && continue
                    worst_re = max(worst_re, abs(imag(w)) / aw)
                end
                for base in (C_D1, C_D1 + 6)          # ∂t, ∂φ
                    w, aw = rotated(MC[base, k], MC[base + 1, k], m)
                    aw == 0 && continue
                    worst_im = max(worst_im, abs(real(w)) / aw)
                end
            end
            @test worst_re < 1e-12
            @test worst_im < 1e-12
        end
    end

    @testset "conjugate symmetry Ψ_{-m} = conj(Ψ_m) (comoving)" begin
        M = EFFSOURCE_PHIS_M_COMOVING
        key = Dict((M[1, k], M[2, k], M[3, k]) => complex(M[6, k], M[7, k])
                   for k in 1:size(M, 2))
        worst = 0.0
        for ((c, m, i), z) in key
            m in (1.0, 2.0) || continue
            zm = key[(c, -m, i)]
            denom = max(abs(z), abs(zm))
            denom == 0 && continue
            worst = max(worst, abs(zm - conj(z)) / denom)
        end
        @test worst < 1e-12             # measured 4.1e-15
    end

    @testset "phase bridge: calc_m_series +sin(mφ_p) convention" begin
        # kerr-circular-comoving.c:673-674 emits src·(cos, +sin) — the OPPOSITE rotation
        # of every other quantity (…:634-635 emits src·(cos, −sin)). So series slots are
        # real under rotation by e^{−i m φ_p} (not e^{+i m φ_p}).
        MS = EFFSOURCE_CALC_M_SERIES_COMOVING
        worst = 0.0
        for k in 1:size(MS, 2)
            w = complex(MS[6, k], MS[7, k]) * cis(-MS[2, k] * EFFSOURCE_PHIP)
            aw = abs(w)
            aw == 0 && continue
            worst = max(worst, abs(imag(w)) / aw)
        end
        @test worst < 1e-12
        # magnitude handoff: the direct Box[Φ^P] src is destroyed by cancellation near the
        # worldline (that is WHY upstream ships calc_m_series) and the series truncates far
        # out — but for EVERY (config, m) the two representations must meet in the overlap
        # region: min-over-points relative |src| difference ≤ HANDOFF_GATE (measured
        # worst-of-min 6.0e-4, at m = 20).
        MC = EFFSOURCE_CALC_M_COMOVING
        srckey = Dict((MC[1, k], MC[2, k], MC[3, k]) =>
                      hypot(MC[C_SRC, k], MC[C_SRC + 1, k]) for k in 1:size(MC, 2))
        minrel = Dict{Tuple{Int,Int},Float64}()
        for k in 1:size(MS, 2)
            c, m, i = MS[1, k], MS[2, k], MS[3, k]
            s = hypot(MS[6, k], MS[7, k])
            sc = srckey[(c, m, i)]
            (isfinite(s) && isfinite(sc) && sc != 0) || continue
            key = (Int(c), Int(m))
            minrel[key] = min(get(minrel, key, Inf), abs(s - sc) / sc)
        end
        @test length(minrel) == 6 * length(EXPECTED_MS)
        @test maximum(values(minrel)) < HANDOFF_GATE
    end

    @testset "cross-variant: 4D Φ^P identical, PhiS_m E/K vs 3P Legendre" begin
        ME, MC = EFFSOURCE_PHIS4D_ELLIPTIC, EFFSOURCE_PHIS4D_COMOVING
        @test all(ME[6, k] === MC[6, k] for k in 1:size(ME, 2))   # same C code path
        E, C = EFFSOURCE_PHIS_M_ELLIPTIC, EFFSOURCE_PHIS_M_COMOVING
        ckey = Dict((C[1, k], C[2, k], C[3, k]) => complex(C[6, k], C[7, k])
                    for k in 1:size(C, 2))
        worst = 0.0
        nwin = 0
        for k in 1:size(E, 2)
            c, m, i = E[1, k], E[2, k], E[3, k]
            XVAR_C1_MIN <= c1_of(Int(c), E[4, k], E[5, k]) < XVAR_C1_MAX || continue
            ze = complex(E[6, k], E[7, k])
            zc = ckey[(c, m, i)]
            denom = max(abs(ze), abs(zc))
            denom == 0 && continue
            nwin += 1
            worst = max(worst, abs(ze - zc) / denom)
        end
        @test nwin > 100                # the frozen grid populates the overlap window
        @test worst < XVAR_GATE
        # calc_m PhiS slot must equal the standalone PhiS_m evaluation (same variant)
        for (PM, CM) in ((EFFSOURCE_PHIS_M_ELLIPTIC, EFFSOURCE_CALC_M_ELLIPTIC),
                         (EFFSOURCE_PHIS_M_COMOVING, EFFSOURCE_CALC_M_COMOVING))
            pkey = Dict((PM[1, k], PM[2, k], PM[3, k]) => (PM[6, k], PM[7, k])
                        for k in 1:size(PM, 2))
            ok = all(pkey[(CM[1, k], CM[2, k], CM[3, k])] ==
                     (CM[C_PHIS, k], CM[C_PHIS + 1, k]) ||
                     !isfinite(CM[C_PHIS, k])
                     for k in 1:size(CM, 2))
            @test ok
        end
    end

    @testset "re-run reproducibility (recompile + rerun == committed bytes)" begin
        gen = joinpath(EFFSOURCE_DIR, "generator")
        toolchain = Sys.which("gcc") !== nothing && Sys.which("make") !== nothing &&
                    success(`pkg-config --exists gsl`)
        # Regeneration needs the read-only upstream EffectiveSource mirror (default = the
        # capture site; override at sites with their own mirror). Committed fixtures stay
        # authoritative either way — absence skips the gate, exactly like a missing toolchain.
        effsrc = get(ENV, "GRTOOLS_EFFSOURCE_ROOT",
                     "/data/haiyangw/claude/emri-disk/ref-code/EffectiveSource")
        if !toolchain
            @test_skip "gcc/make/GSL not available — reproducibility gate skipped"
        elseif !isfile(joinpath(effsrc, "kerr-circular.c"))
            @test_skip "EffectiveSource mirror absent — reproducibility gate skipped"
        else
            mktempdir() do tmp
                out1, out2 = joinpath(tmp, "run1"), joinpath(tmp, "run2")
                for (bdir, out) in ((joinpath(tmp, "b1"), out1), (joinpath(tmp, "b2"), out2))
                    run(pipeline(`make -C $gen fixtures BUILD_DIR=$bdir OUTDIR=$out
                                  EFFSOURCE_SRC=$effsrc`,
                                 stdout = devnull, stderr = devnull))
                end
                names = filter(f -> endswith(f, ".jl"), readdir(out1))
                @test !isempty(names)
                # run-to-run byte identity (fresh compiles)
                @test all(read(joinpath(out1, f)) == read(joinpath(out2, f)) for f in names)
                # regenerated bytes == committed fixture bytes
                @test all(read(joinpath(out1, f)) == read(joinpath(EFFSOURCE_DIR, f))
                          for f in names)
            end
        end
    end
end

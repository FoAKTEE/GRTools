# test_effective_source_series.jl — ESSF node `effective-source-series` (schedule W8-B,
# GLOBAL_DAG row 37), written BEFORE the implementation (TDD). Target:
# `EffectiveSource/series.jl` + the generated `EffectiveSource/series_coeffs.jl` — the
# NEAR-WORLDLINE POWER SERIES of the m-mode effective source (arXiv:2606.04998
# main.tex:488-504) and the deterministic SWITCH that trades it against the sealed direct
# route of node `effective-source-direct`.
#
# THE SERIES (main.tex:493-504, coefficients from ref-code/EffectiveSource
# kerr-circular-comoving.c `effsource_calc_m_series`, SHA 1cb9834f):
#
#   s(Δr, Δθ) = S^[0] + S^[1]_1 Δr
#             + [ S^[2]_2 Δr² + S^[2]_0 ϱ² + (S^[2]_{2L} Δr² + S^[2]_{0L} ϱ²) log ϱ
#                 + S^[2]_4 Δr⁴/ϱ² + S^[2]_6 Δr⁶/ϱ⁴ + S^[2]_8 Δr⁸/ϱ⁶
#                 + S^[2]_10 Δr¹⁰/ϱ⁸ + S^[2]_12 Δr¹²/ϱ¹⁰ ] + O(λ³ log λ),
#
# and S^eff_m = s · e^{−im(φ_p + χ(r))} · rΣ / (2π(r²+a²) sin^{|m|}θ) — the SAME convention
# chain the sealed `effective_source_m` produces (node 36), transcribed INDEPENDENTLY below
# as `_ess_chain` from spectre CircularOrbit.cpp:322-332 + main.tex:210-262.
#
# THE C'S GROUPING vs THE PAPER'S. The C writes the H/log block as
# `H_{−m−1/2}·(X0 ϱ² + X2 Δr²) + log(ϱ²/z_c²)·(X0 ϱ² + X2 Δr²)/2`, i.e.
# `(X0 ϱ² + X2 Δr²)·(H_{−m−1/2} + log(ϱ/z_c))`. So S^[2]_{0L} = X0, S^[2]_{2L} = X2 (both
# m-INDEPENDENT) and the remainder `X·(H_{−m−1/2} − log z_c)` is absorbed into S^[2]_0 /
# S^[2]_2 — the only two coefficients the paper says carry harmonic numbers (main.tex:503).
#
# GATE CALIBRATION (PRE-RED; evidence dir effective-source-series/proto_check{1..9}.jl,
# switch_table.tsv, calibration logs; never tuned post hoc). The measurements that fix every
# constant below:
#   * the transcription reproduces the baked compiled-C fixture to 2.2095e-13 at λ ≤ 1e-3
#     (360 columns) and 2.0119e-09 over the whole 1410-column grid; the whole residual is the
#     C's OWN Float64 coefficient roundoff at the near-extremal config 4 (a = 0.998 at its
#     ISCO): configs 1, 2, 3, 5, 6 agree to ≤ 7.22e-14 at EVERY λ.
#   * the remainder is O(λ³ log λ), NOT O(λ³): the local log-log slope of |series − direct|
#     is 3 + 1/ln λ (measured 2.94 at λ ~ 1e-6 down to 2.87 at λ ~ 1e-2), so the DAG's
#     planned "slope 3.0 ± 0.2 over λ ∈ [1e-4, 1e-2]" is UNREACHABLE BY CONSTRUCTION
#     (measured raw envelope [2.5635, 3.5011] over 36 config×m fits). Gated instead on the
#     fitted exponent p of the model `rel = K λ^p |ln λ|`: measured [2.9409, 3.0136] over
#     λ ∈ [1e-7, 1e-4] for all 36, and [2.7132, 3.1289] on the DAG's own window for m ≤ 5.
#   * the direct route's actual Float64 error is a large multiple of its `eps·cancellation`
#     diagnostic — measured median 91, p90 437, worst 2310 over the 2148 calibration points
#     where the direct value still has ≥ 3 correct digits (307 in the sealed node's own box).
#     The switch carries 1e3. MEASURED ALTERNATIVES on the calibration grid (evidence
#     `calibration_switch.log` / `repair_select_rule.py`): ULPS = 300 → worst
#     achieved/optimal factor 374.3, 1000 → 42.1, 3000 → 36.3, all three with the SAME worst
#     switched relative error; 1000 was CHOSEN because it minimises the near-crossover route
#     disagreement (1.2e-3 versus 4.0e-2 at 3000). ULPS = 2310 was never run — an earlier
#     revision of this header claimed it, and that claim is withdrawn.
#   * THE ESTIMATOR (validator finding H1, repaired here). `q = |bracket|/|s|` uses the
#     SIGNED sum of the nine λ² terms, which vanishes on codimension-1 curves in (r, z); the
#     switch is therefore NOT `q·max(q, (1+|m|)λ/r₀)` alone but carries an unconditional
#     geometric floor `((1+|m|)λ/r₀)³` that no coefficient zero can defeat. Measured
#     alternatives (`repair_select_rule.py`, scored on the 3240-point grid AND on the 46
#     located bracket zeros): replacing `q` by the bracket's Σ|tᵢ| NORM is equally safe at
#     the zeros but costs three orders on legitimate points (worst factor 1.87e4 vs 42.1,
#     worst switched error 5.67e-3 vs 6.95e-4) because the λ² bracket has a PHYSICAL internal
#     cancellation the λ³ order shares (Σ|tᵢ|/|bracket| median 2.4, p99 169, max 1.5e4);
#     adding `Σ|tᵢ|·λ/r₀` on top of the floor costs 8× in the worst factor for no measured
#     safety gain. The floor alone is the active term at 46/46 zeros with ≥ 167× margin.
#   * ON THE 3240-POINT CALIBRATION GRID `test/fixtures/effsource/switch_calibration.jl`
#     (6 configs × 6 m × 9 λ × 5 directions × both Δr signs) the repaired switch's worst
#     relative error is 6.9542e-04 against an ORACLE-optimal 5.8608e-04 (ratio 1.187), while
#     on the same points the direct route ALONE reaches 6.0385e+02 and the series ALONE
#     4.0122e+00; worst pointwise achieved/optimal factor 42.106, worst near-crossover route
#     disagreement 1.1999e-03, and for λ ≤ 3.2e-4 the rule picks the series at every one of
#     the 1080 points with factor 1. These are properties OF THAT GRID: an independent
#     9072-point sweep run by the validator reports 1.23e-3 / 2389 / 4.07e-3 for the same
#     three numbers — all in the conservative direction, none a different conclusion.
using Test
using ScalarSelfForce
using ScalarSelfForce: _series_coeffs, _series_harmonic_mhalf, kerr_varphi,
                       kerr_varphi_offset, sinm_theta
using BigSpecial: digamma_big

if !@isdefined(EFFSOURCE_CONFIGS)
    include(joinpath(@__DIR__, "fixtures", "effsource", "effsource_configs.jl"))
end
if !@isdefined(EFFSOURCE_CALC_M_SERIES_COMOVING)
    include(joinpath(@__DIR__, "fixtures", "effsource", "calc_m_series_comoving.jl"))
end
if !@isdefined(EFFSOURCE_SERIES_RATIONAL)
    include(joinpath(@__DIR__, "fixtures", "effsource", "series_coeffs_rational.jl"))
end
if !@isdefined(EFFSOURCE_SWITCH_CALIBRATION)
    include(joinpath(@__DIR__, "fixtures", "effsource", "switch_calibration.jl"))
end

# The bracket-zero hunt of gate (iii) costs ~4 min of BigFloat(256) adjudication; it runs
# under `GRTOOLS_SLOW_TESTS=1` (the repo-wide convention, cf. AccretionDisc
# test_newtonian_planet.jl). Everything else, including the H1 reproducer and the ray walk
# that would have caught it, is always on.
const _ESS_SLOW = get(ENV, "GRTOOLS_SLOW_TESTS", "0") == "1"

# ---- in-test instruments (uniquely named: three test files already live in `Main`) -------
_ess_lambda(o, dr, dth) = sqrt(o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dth^2)

# a point at exactly ϱ = L on the ray with radial fraction `frac` of the ϱ metric
function _ess_ray(o, L, frac)
    dr = frac * L / sqrt(o.rho2_dr2)
    dth = sqrt(max(zero(L), one(L) - frac^2)) * L / sqrt(o.rho2_dtheta2)
    return (o.r0 + dr, -sin(dth))
end

# THE C → package convention chain. Built through the SEALED ansatz layer rather than by
# rewriting the implementation's own expression: `mode_amplitude` is the single home of the
# `r/sin^{|m|}θ` stripping and of the `e^{+imΩt}` convention (ansatz.jl `mode_term`), so
# feeding it the C's real scalar with the horizon-regular field-point rotation and the
# `Σ/(r²+a²)` mode-projection scale of eq:elliptic_eqn reproduces the chain from the sealed
# primitives only. (This is a CONVENTION cross-check, not an independent derivation of the
# chain — the matching-zone testset is what independently covers it, by comparing against the
# sealed `effective_source_m`, which reaches the same composition through a different path.)
function _ess_chain(o, m, r, z, phi_p)
    Sigma = r^2 + o.a^2 * z^2
    amp = mode_amplitude(m, cis(-m * kerr_varphi(o.a, r, phi_p)) /
                            (2 * oftype(r, pi));
                         Omega = o.Omega, ts = zero(r), r = r, z = z)
    return amp * Sigma / (r^2 + o.a^2)
end

# least-squares slope of y vs x
function _ess_slope(x, y)
    n = length(x)
    return (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2)
end

# fitted exponent p of the model  rel(λ) = K λ^p |ln λ|   (`logcorr = false` ⇒ the RAW slope)
function _ess_fit_exponent(o, p, m, es, phi_p; logcorr::Bool = true, frac = big"0.6")
    x = Float64[]
    y = Float64[]
    for e in es
        L = big(10)^e
        r, z = _ess_ray(o, L, frac)
        ser = effective_source_series_m(p, m, r, z; phi_p = phi_p)
        dir = effective_source_m(p, m, r, z; phi_p = phi_p)
        rel = Float64(abs(ser - dir) / abs(dir))
        push!(x, Float64(e))
        push!(y, log10(logcorr ? rel / abs(log(10.0^Float64(e))) : rel))
    end
    return _ess_slope(x, y)
end

_ess_rat(nu, de) = Rational{BigInt}(parse(BigInt, nu), parse(BigInt, de))

# ALLOCATION PROBES as function barriers. `@allocated` written inline in a `@testset` also
# measures the testset's OWN boxing of any variable that a soft-scope loop above it assigns
# (measured 32-48 B of caller boxing, not callee allocation); a barrier taking every input as
# an argument measures the callee only.
_ess_alloc_coeffs(p, m) = @allocated effective_source_series_coefficients(p, m)
_ess_alloc_series(p, m, r, z, pp) =
    @allocated effective_source_series_m(p, m, r, z; phi_p = pp)
_ess_alloc_switched(p, m, r, z, pp) =
    @allocated effective_source_m_switched(p, m, r, z; phi_p = pp)
_ess_alloc_bundle(p, m, r, z, pp) =
    @allocated effective_source_switched_bundle(p, m, r, z; phi_p = pp)

const _ESS_COEFF_NAMES = (:S0, :S1_1, :S2_0rat, :S2_2rat, :S2_0L, :S2_2L,
                          :S2_4, :S2_6, :S2_8, :S2_10, :S2_12)
const _ESS_PAPER_NAMES = (:S0, :S1_1, :S2_2, :S2_0, :S2_2L, :S2_0L,
                          :S2_4, :S2_6, :S2_8, :S2_10, :S2_12)

@testset "EffectiveSource: near-worldline series + switch (effective-source-series)" begin

    orbits = Dict(c.idx => CircularOrbit(c.a, c.r0) for c in EFFSOURCE_CONFIGS)
    puncs = Dict(c.idx => Puncture(orbits[c.idx]) for c in EFFSOURCE_CONFIGS)
    phi_p = EFFSOURCE_PHIP        # 0.31, the frozen-grid particle BOYER–LINDQUIST azimuth

    @testset "H_{−m−1/2}: the exact half-integer harmonic number" begin
        # The C calls GSL's ψ (`M_EULER + gsl_sf_psi(1 + x)` at x = −1/2 − m). We do NOT
        # re-implement a digamma: at half-integer argument the recurrence terminates exactly,
        #     H_{−m−1/2} = γ + ψ(1/2 − m) = 2·Σ_{k=1}^{m} 1/(2k−1) − 2 log 2 ,
        # which is EXACT at every precision (BigSpecial's `digamma_big` is an asymptotic
        # series with a ≈1e-16 floor and would cap BigFloat inheritance).
        @test _series_harmonic_mhalf(Float64, 0) ≈ -2 * log(2.0) rtol = 1e-15
        @test _series_harmonic_mhalf(Float64, 1) ≈ 2 - 2 * log(2.0) rtol = 1e-15
        @test _series_harmonic_mhalf(Float64, 2) ≈ 2 + 2 / 3 - 2 * log(2.0) rtol = 1e-15
        worst = 0.0
        for m in 0:20
            ref = Float64(MathConstants.eulergamma) +
                  real(digamma_big(complex(0.5 - m, 0.0)))
            v = _series_harmonic_mhalf(Float64, m)
            @test abs(v - ref) <= 1e-13 * max(abs(ref), 1.0)
            worst = max(worst, abs(v - ref) / max(abs(ref), 1.0))
        end
        println("effective-source-series H_{−m−1/2} vs BigSpecial digamma_big: worst = ", worst)
        setprecision(BigFloat, 256) do
            v = _series_harmonic_mhalf(BigFloat, 7)
            @test v isa BigFloat
            @test precision(v) == 256
            @test abs(v - (2 * sum(one(BigFloat) / (2k - 1) for k in 1:7) - 2 * log(BigFloat(2)))) <=
                  4 * eps(BigFloat)
        end
    end

    @testset "coefficients: EXACT-rational double transcription (198 baked rows)" begin
        # An INDEPENDENT second transcription of the same C text (fixture
        # series_coeffs_rational.jl): Python `fractions.Fraction` evaluation of the C body
        # used verbatim as Python source, with the eleven coefficients extracted by
        # SUBSTITUTION + exact Vandermonde solve instead of the shipped generator's textual
        # surgery. Zero tolerance: `Rational{BigInt}` equality.
        nok = 0
        for (ip, P) in enumerate(EFFSOURCE_SERIES_RATIONAL_POINTS)
            a = Rational{BigInt}(numerator(P.a), denominator(P.a))
            r0 = Rational{BigInt}(numerator(P.r0), denominator(P.r0))
            Om = Rational{BigInt}(numerator(P.Omega), denominator(P.Omega))
            one_ = one(Rational{BigInt})
            for m in unique(r[2] for r in EFFSOURCE_SERIES_RATIONAL)
                c = _series_coeffs(one_, a, r0, Om, one_, m)
                for (i, mm, name, nu, de) in EFFSOURCE_SERIES_RATIONAL
                    (i == ip && mm == m) || continue
                    @test getfield(c, name) == _ess_rat(nu, de)
                    nok += 1
                end
            end
            # every C term carries exactly ONE linear factor of z_c ⇒ exact homogeneity
            c1 = _series_coeffs(one_, a, r0, Om, one_, 5)
            c7 = _series_coeffs(one_, a, r0, Om, 7 * one_ // 3, 5)
            for n in _ESS_COEFF_NAMES
                @test getfield(c7, n) == (7 * one_ // 3) * getfield(c1, n)
            end
        end
        @test nok == length(EFFSOURCE_SERIES_RATIONAL)
        println("effective-source-series exact-rational double transcription: ", nok,
                " coefficients matched EXACTLY (0 tolerance)")
    end

    @testset "coefficients: the paper's m-dependence (main.tex:503-504)" begin
        # Asserted EXACTLY in Rational{BigInt} at the fixture's algebraic probe points:
        #   S^[0], S^[1]_1 ∝ 1/((2m−1)(2m+1)) and nothing else;
        #   S^[2]_0, S^[2]_2 = (degree-1 polynomial in m²)/((2m−1)(2m+1)) + H_{−m−1/2}·(const);
        #   every other coefficient is m-INDEPENDENT.
        ms = [0, 1, 2, 3, 5, 7, 10, 20]
        for P in EFFSOURCE_SERIES_RATIONAL_POINTS
            a = Rational{BigInt}(numerator(P.a), denominator(P.a))
            r0 = Rational{BigInt}(numerator(P.r0), denominator(P.r0))
            Om = Rational{BigInt}(numerator(P.Omega), denominator(P.Omega))
            one_ = one(Rational{BigInt})
            C = [_series_coeffs(one_, a, r0, Om, one_, m) for m in ms]
            for n in (:S0, :S1_1)
                v = [getfield(C[i], n) * ((2ms[i] - 1) * (2ms[i] + 1)) for i in eachindex(ms)]
                @test all(v .== v[1])
            end
            for n in (:S2_0L, :S2_2L, :S2_4, :S2_6, :S2_8, :S2_10, :S2_12)
                v = [getfield(C[i], n) for i in eachindex(ms)]
                @test all(v .== v[1])
            end
            for n in (:S2_0rat, :S2_2rat)
                v = [getfield(C[i], n) * ((2ms[i] - 1) * (2ms[i] + 1)) for i in eachindex(ms)]
                x = [Rational{BigInt}(m)^2 for m in ms]
                # exactly linear in m² (the C's highest power of m in these blocks is m²)
                slope = (v[2] - v[1]) // (x[2] - x[1])
                @test all(v[k] == v[1] + slope * (x[k] - x[1]) for k in eachindex(ms))
                @test slope != 0
            end
        end
        # the assembled (paper) coefficients: only S^[2]_0 and S^[2]_2 carry the harmonic
        # number, and their H-content is EXACTLY S^[2]_{0L} / S^[2]_{2L}
        setprecision(BigFloat, 256) do
            for c in (EFFSOURCE_CONFIGS[1], EFFSOURCE_CONFIGS[4])
                oB = CircularOrbit(big(c.a), big(c.r0))
                pB = Puncture(oB)
                zc = sqrt(oB.zc2)
                A = [effective_source_series_coefficients(pB, m) for m in ms]
                for (i, m) in enumerate(ms)
                    HL = _series_harmonic_mhalf(BigFloat, m) - log(zc)
                    raw = _series_coeffs(oB.M, oB.a, oB.r0, oB.Omega, zc, m)
                    @test abs(A[i].S2_0 - (raw.S2_0rat + raw.S2_0L * HL)) <=
                          1e-60 * abs(A[i].S2_0)
                    @test abs(A[i].S2_2 - (raw.S2_2rat + raw.S2_2L * HL)) <=
                          1e-60 * abs(A[i].S2_2)
                    for n in (:S0, :S1_1, :S2_0L, :S2_2L, :S2_4, :S2_6, :S2_8, :S2_10, :S2_12)
                        @test getfield(A[i], n) == getfield(raw, n)
                    end
                end
                # |m| reading: coefficients are even in m (the mode function's conjugation)
                @test all(getfield(effective_source_series_coefficients(pB, -5), n) ==
                          getfield(effective_source_series_coefficients(pB, 5), n)
                          for n in _ESS_PAPER_NAMES)
            end
        end
    end

    @testset "S^eff_m vs the baked effsource_calc_m_series fixture" begin
        # The fixture's src pair carries the OPPOSITE phase convention to `calc_m`
        # (s·e^{+imφ_p}, kerr-circular-comoving.c:673-674 vs :634-635 — the W2-D header's
        # "+sin(mφ_p) bridge"), so the underlying REAL scalar is s = Re(V·e^{−imφ_p}); the
        # bridge is re-derived here by asserting the imaginary part vanishes.
        fix = EFFSOURCE_CALC_M_SERIES_COMOVING
        worst = zeros(4)
        nwin = zeros(Int, 4)
        worst_tame = 0.0
        for col in 1:size(fix, 2)
            cfg = Int(fix[1, col])
            m = Int(fix[2, col])
            o = orbits[cfg]
            p = puncs[cfg]
            dr = fix[4, col]
            dth = fix[5, col]
            r = o.r0 + dr
            z = cos(pi / 2 + dth)
            V = complex(fix[6, col], fix[7, col])
            sC = V * cis(-m * phi_p)
            @test abs(imag(sC)) <= 1e-12 * max(abs(V), 1e-300)     # the phase bridge
            theirs = real(sC) * _ess_chain(o, m, r, z, phi_p)
            ours = effective_source_series_m(p, m, r, z; phi_p = phi_p)
            rel = abs(ours - theirs) / max(abs(theirs), abs(ours))
            L = _ess_lambda(o, dr, dth)
            k = L <= 1e-3 ? 1 : (L < 1e-2 ? 2 : (L < 1e-1 ? 3 : 4))
            worst[k] = max(worst[k], rel)
            nwin[k] += 1
            cfg == 4 || (worst_tame = max(worst_tame, rel))
        end
        @test nwin == [360, 144, 312, 594]
        @test worst[1] <= 1e-10                 # THE row-37 gate (measured 2.2095e-13)
        @test maximum(worst) <= 1e-8            # whole grid (measured 2.0119e-09)
        # away from the near-extremal config 4 the agreement is pure Float64 roundoff: the
        # residual there is the C's OWN coefficient conditioning at a = 0.998 on its ISCO
        @test worst_tame <= 1e-12               # measured 7.22e-14
        println("effective-source-series vs calc_m_series fixture: worst rel by λ window ",
                "[≤1e-3, 1e-3..1e-2, 1e-2..1e-1, ≥1e-1] = ", worst, " over ", nwin,
                " columns; worst excluding cfg 4 (a=0.998) = ", worst_tame)
    end

    @testset "worldline limit: series(λ → 0) == S^[0]" begin
        # The C guards its λ-dependent orders with `(dr==0&&dth==0)? 0 :`; ϱ = 0 is a genuine
        # point of the SERIES (unlike the direct route, which diverges there and throws).
        # The S^[0] comparison is to a few ulps, NOT bitwise: `_ess_chain` is deliberately
        # built from the sealed `mode_amplitude` by a different association than the
        # implementation's own chart factor (see its comment), so bitwise equality here would
        # only re-assert that the test copied the implementation. What IS exact is asserted
        # instead: the λ-dependent orders vanish identically at ϱ = 0 (so the switched
        # evaluator returns the series value bit-for-bit), and one ulp-scale step off the
        # worldline the whole residual is the O(λ) term.
        worst = 0.0
        for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 10, 20)
            o = orbits[c.idx]
            p = puncs[c.idx]
            A = effective_source_series_coefficients(p, m)
            got = effective_source_series_m(p, m, o.r0, 0.0; phi_p = phi_p)
            want = A.S0 * _ess_chain(o, m, o.r0, 0.0, phi_p)
            @test abs(got - want) <= 4 * eps(Float64) * abs(want)   # measured ≤ 2 ulps
            @test isfinite(got)
            # EXACT: on the worldline circle the switch has no direct route to compare with,
            # so it must return the series value bit-for-bit
            @test effective_source_m_switched(p, m, o.r0, 0.0; phi_p = phi_p) === got
            @test effective_source_switched_bundle(p, m, o.r0, 0.0;
                                                   phi_p = phi_p).route === :series
            worst = max(worst, abs(got - want) / abs(want))
            # the limit is approached continuously, and what remains one ulp-scale step off
            # the worldline is EXACTLY the O(λ) term (the λ² bracket is ~1e-23 relative
            # there): a stricter statement than "close to S^[0]", which the S^[1]_1 Δr term
            # itself violates at 1.9e-10 relative.
            dr0 = 1e-12
            ch = _ess_chain(o, m, o.r0 + dr0, 0.0, phi_p)
            near = effective_source_series_m(p, m, o.r0 + dr0, 0.0; phi_p = phi_p)
            @test abs(near - (A.S0 + A.S1_1 * dr0) * ch) <= 1e-13 * abs(want)
            @test abs(near - want) <= 1e-8 * abs(want)
        end
        println("effective-source-series worldline limit |series(0) − S^[0]·chain|/|S^[0]·chain| = ",
                worst)
        # the DIRECT route throws on the worldline circle and at the poles; the SERIES only
        # at the poles (its ϱ = 0 value is the S^[0] limit above)
        p1 = puncs[1]
        @test_throws DomainError effective_source_m(p1, 2, orbits[1].r0, 0.0)
        @test_throws DomainError effective_source_series_m(p1, 2, orbits[1].r0 + 0.1, 1.0)
        @test_throws DomainError effective_source_series_m(p1, 2, orbits[1].r0 + 0.1, -1.0)
        @test isfinite(effective_source_m_switched(p1, 2, orbits[1].r0, 0.0; phi_p = phi_p))
    end

    @testset "matching zone: the remainder is O(λ³ log λ)" begin
        # THE row-37 matching gate, with the direct route at BigFloat(256) so its λ^{-2}
        # roundoff cannot pollute the slope. The DAG's planned RAW slope 3.0 ± 0.2 over
        # λ ∈ [1e-4, 1e-2] is unreachable BY CONSTRUCTION: the local slope of K λ³|ln λ| is
        # 3 + 1/ln λ, i.e. 2.86…2.94 over that window (measured raw envelope [2.5635,
        # 3.5011]). The physics claim — "accurate through order λ²" — is the fitted EXPONENT.
        # DIRECTION: everything here is measured on the frac = 0.6 ray (Δr and Δθ both
        # nonzero). The exponent is DIRECTION-DEPENDENT and this is the generic case: on the
        # pure-Δθ ray (Δr = 0) every odd order vanishes by the puncture's evenness in Δθ
        # (A_ijk = 0 for odd j, main.tex:457) and the measured exponent is 4, i.e. the series
        # is one order BETTER there — the safe direction, and the reason the switch's
        # estimator is conservative on that ray.
        setprecision(BigFloat, 256) do
            ppB = big(phi_p)
            lo = Inf
            hi = -Inf
            rawlo = Inf
            rawhi = -Inf
            nfit = 0
            for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 10, 20)
                oB = CircularOrbit(big(c.a), big(c.r0))
                pB = Puncture(oB)
                pe = _ess_fit_exponent(oB, pB, m, big(-7):big(1)//2:big(-4), ppB)
                pr = _ess_fit_exponent(oB, pB, m, big(-7):big(1)//2:big(-4), ppB;
                                       logcorr = false)
                @test abs(pe - 3) <= 0.10          # measured [2.9409, 3.0136] over all 36
                lo = min(lo, pe); hi = max(hi, pe)
                rawlo = min(rawlo, pr); rawhi = max(rawhi, pr)
                nfit += 1
            end
            @test nfit == 36                       # ALL config x m, not a subset
            println("effective-source-series matching zone λ ∈ [1e-7, 1e-4]: fitted p ∈ [",
                    lo, ", ", hi, "] ; RAW log-log slope ∈ [", rawlo, ", ", rawhi, "]")
            # the DAG's own window, where the λ⁴ term is already entering for m ≳ 10
            lo2 = Inf
            hi2 = -Inf
            for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5)
                oB = CircularOrbit(big(c.a), big(c.r0))
                pB = Puncture(oB)
                pe = _ess_fit_exponent(oB, pB, m, big(-4):big(1)//2:big(-2), ppB)
                @test abs(pe - 3) <= 0.35          # measured [2.7132, 3.1289] for m ≤ 5
                lo2 = min(lo2, pe); hi2 = max(hi2, pe)
            end
            println("effective-source-series matching zone λ ∈ [1e-4, 1e-2] (m ≤ 5): ",
                    "fitted p ∈ [", lo2, ", ", hi2, "]")
            # the pure-Δθ ray is one order better (the evenness argument above)
            oB = CircularOrbit(big"0.5", big"10.0")
            pB = Puncture(oB)
            p_theta = _ess_fit_exponent(oB, pB, 2, big(-7):big(1)//2:big(-4), ppB;
                                        frac = big(0))
            @test p_theta > 3.5
            println("effective-source-series matching zone on the PURE-Δθ ray (Δr = 0): ",
                    "fitted p = ", p_theta, " (odd orders vanish by Δθ-evenness)")
        end
    end

    @testset "the switch: route choice, achieved accuracy, crossover, continuity" begin
        setprecision(BigFloat, 256) do
            ppB = big(phi_p)
            worstF = 0.0
            worst_abs = 0.0
            worst_opt = 0.0
            worst_jump = 0.0
            nser = 0
            ndir = 0
            nsmall = 0
            njump = 0
            for ic in (1, 4, 6), m in (0, 2, 20), e in -5.0:0.5:-2.0, frac in (0.0, 0.5, 1.0)
                c = EFFSOURCE_CONFIGS[ic]
                o = orbits[ic]
                p = puncs[ic]
                oB = CircularOrbit(big(c.a), big(c.r0))
                pB = Puncture(oB)
                L = 10.0^e
                r, z = _ess_ray(o, L, frac)
                truth = effective_source_m(pB, m, big(r), big(z); phi_p = ppB)
                b = effective_source_switched_bundle(p, m, r, z; phi_p = phi_p)
                @test b.route in (:series, :direct)
                @test b.source == effective_source_m_switched(p, m, r, z; phi_p = phi_p)
                ser = effective_source_series_m(p, m, r, z; phi_p = phi_p)
                dir = effective_source_m(p, m, r, z; phi_p = phi_p)
                @test b.source == (b.route === :series ? ser : dir)
                rs = Float64(abs(big(ser) - truth) / abs(truth))
                rd = Float64(abs(big(dir) - truth) / abs(truth))
                ra = Float64(abs(big(b.source) - truth) / abs(truth))
                worstF = max(worstF, ra / max(min(rs, rd), 1e-300))
                worst_abs = max(worst_abs, ra)
                worst_opt = max(worst_opt, min(rs, rd))
                # the JUMP matters only where the rule is actually deciding: within a
                # decade of the crossover, where by construction both predicted errors are
                # equal to the local accuracy limit
                if 0.1 <= b.err_series / b.err_direct <= 10
                    njump += 1
                    worst_jump = max(worst_jump,
                                     Float64(abs(big(ser) - big(dir)) / abs(truth)))
                end
                b.route === :series ? (nser += 1) : (ndir += 1)
                if L <= 3.2e-4                # deep near-zone: the series must always win
                    @test b.route === :series
                    @test ra == rs
                    nsmall += 1
                end
            end
            @test nsmall == 3 * 3 * 4 * 3
            @test nser + ndir == 3 * 3 * 7 * 3
            @test ndir > 0                     # the switch really does select both routes
            # GATE MOVED 10 → 30 by the H1 repair, with the measurement: the unconditional
            # geometric floor declines the series on a handful of m = 20, λ ≈ 1e-3 rows where
            # it happens to be anomalously accurate (rs 8.2e-8 vs rd 3.4e-6), which is the
            # price of an estimator a bracket zero cannot defeat. Measured 19.921 here, and
            # 42.106 on the 3240-point grid — UNCHANGED from before the repair there. The two
            # error gates below did not move: worst switched error 1.7430e-04 and the
            # oracle-ratio 2.040 are bit-identical pre- and post-repair.
            @test worstF <= 30.0               # measured 19.921 here (42.106 on the 3240-pt grid)
            @test worst_abs <= 5e-4            # measured 1.7430e-04 (unchanged by the repair)
            @test worst_abs <= 3 * worst_opt   # measured ratio 2.040 (unchanged by the repair)
            @test njump > 0
            @test worst_jump <= 2e-3           # measured 3.6e-4 here (1.2e-3 on the 3240-pt grid)
            println("effective-source-series switch: series ", nser, " / direct ", ndir,
                    " ; worst achieved/optimal factor = ", worstF,
                    " ; worst switched rel = ", worst_abs, " (oracle-optimal ", worst_opt,
                    ") ; worst |series − direct|/|S| at ", njump,
                    " near-crossover points = ", worst_jump)
        end
        # CROSSOVER: λ* moves inward like eps^{1/5} (the series error ~ λ³, the direct
        # roundoff ~ eps·λ^{-2}), so at BigFloat(256) the direct route is usable far closer
        # to the worldline than at Float64. Reported, and gated on the ORDERING + scaling.
        function _ess_crossover(o, p, m, T, ppT; iters = 60)
            lo = T(-20)
            hi = T(0)
            for _ in 1:iters
                mid = (lo + hi) / 2
                r, z = _ess_ray(o, T(10)^mid, T(3) / 5)
                if effective_source_switched_bundle(p, m, r, z; phi_p = ppT).route === :series
                    lo = mid
                else
                    hi = mid
                end
            end
            return Float64((lo + hi) / 2)
        end
        x64 = _ess_crossover(orbits[1], puncs[1], 2, Float64, phi_p)
        x106 = setprecision(BigFloat, 106) do
            o = CircularOrbit(big(0.5), big(10)); _ess_crossover(o, Puncture(o), 2, BigFloat,
                                                                big(phi_p); iters = 40)
        end
        x256 = setprecision(BigFloat, 256) do
            o = CircularOrbit(big(0.5), big(10)); _ess_crossover(o, Puncture(o), 2, BigFloat,
                                                                big(phi_p); iters = 30)
        end
        println("effective-source-series crossover log10 λ* (cfg 2, m = 2): Float64 ", x64,
                " | BigFloat(106) ", x106, " | BigFloat(256) ", x256)
        @test x64 > x106 > x256
        @test -3.5 < x64 < -1.0
        # eps^{1/5} scaling between precisions, within a decade
        e64 = log10(eps(Float64))
        e106 = setprecision(() -> log10(Float64(eps(BigFloat))), BigFloat, 106)
        e256 = setprecision(() -> log10(Float64(eps(BigFloat))), BigFloat, 256)
        @test abs((x64 - x106) - (e64 - e106) / 5) <= 1.0
        @test abs((x64 - x256) - (e64 - e256) / 5) <= 1.0
    end

    @testset "the switch vs the ZEROS of the signed λ² bracket (validator finding H1)" begin
        # H1: `q = |bracket|/|s|` is built from the SIGNED sum of the nine λ² terms, which
        # vanishes on codimension-1 curves in (r, z). Across such a curve q collapses by
        # 4-16 orders while the TRUE truncation error does not, so an estimator resting on q
        # alone selects `:series` arbitrarily far outside its validity and
        # `effective_source_m_switched` — the DG fixed source of nodes 38/40/42 — returns an
        # O(1)-wrong value with no diagnostic. The repair is the unconditional geometric
        # floor `((1+|m|)ϱ/r₀)³` in `effective_source_switched_bundle`.
        o = orbits[2]                                    # a = 0.5, r₀ = 10
        p = puncs[2]

        # (i) THE REPRODUCER, verbatim from the validator report. Inside the production
        #     worldtube (λ = 0.105). Pre-repair: route = :series, 93.5 % wrong.
        for (rr, zz, want_dir) in ((10.034, 0.0098310, 3.4553e-6),
                                   (10.0407963, 0.0088854, 1.2163e-5))
            b = effective_source_switched_bundle(p, 20, rr, zz; phi_p = phi_p)
            dir = effective_source_m(p, 20, rr, zz; phi_p = phi_p)
            ser = effective_source_series_m(p, 20, rr, zz; phi_p = phi_p)
            @test b.route === :direct
            @test b.source == dir
            # the series really is unusable there — the point is not benign
            @test abs(ser - dir) / abs(dir) > 0.5
            # and the switched value is within the direct route's own predicted error
            @test abs(b.source - dir) / abs(dir) <= b.err_direct
        end

        # (ii) RAY WALK across a sign change of the bracket at fixed λ, m = 20. The bracket
        #      is reconstructed from the EXPORTED coefficients (the paper's own formula), so
        #      this locates the zeros without reaching into the implementation.
        A = effective_source_series_coefficients(p, 20)
        function bracket(dr, dth)
            rho2 = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dth^2
            d2 = dr^2
            lg = log(rho2)
            return A.S2_2 * d2 + A.S2_0 * rho2 + (A.S2_2L * d2 + A.S2_0L * rho2) * lg / 2 +
                   A.S2_4 * d2^2 / rho2 + A.S2_6 * d2^3 / rho2^2 + A.S2_8 * d2^4 / rho2^3 +
                   A.S2_10 * d2^5 / rho2^4 + A.S2_12 * d2^6 / rho2^5
        end
        L = 0.1
        nzero = 0
        worst_at_zero = 0.0
        prev = nothing
        prevt = -1.0
        for t in range(-0.999, 0.999; length = 2001)
            r, z = _ess_ray(o, L, t)
            b1 = bracket(r - o.r0, -asin(z))
            if prev !== nothing && sign(b1) != sign(prev)
                lo, hi = prevt, t
                for _ in 1:60
                    mid = (lo + hi) / 2
                    rm, zm = _ess_ray(o, L, mid)
                    sign(bracket(rm - o.r0, -asin(zm))) == sign(prev) ? (lo = mid) : (hi = mid)
                end
                nzero += 1
                # AT the zero and on BOTH sides of it, within the band where q is collapsed
                for dt in (-1e-4, -1e-6, 0.0, 1e-6, 1e-4)
                    tt = (lo + hi) / 2 + dt
                    abs(tt) < 0.999 || continue
                    r, z = _ess_ray(o, L, tt)
                    bb = effective_source_switched_bundle(p, 20, r, z; phi_p = phi_p)
                    dd = effective_source_m(p, 20, r, z; phi_p = phi_p)
                    # at λ = 0.1 the direct route is trustworthy (its own predicted error is
                    # ~1e-4 of the value), so |switched − direct| is a valid misfire probe
                    @test bb.err_direct < 1e-3
                    rel = abs(bb.source - dd) / abs(dd)
                    @test rel <= max(1e-3, 10 * bb.err_direct)
                    worst_at_zero = max(worst_at_zero, rel)
                end
            end
            prev = b1
            prevt = t
        end
        @test nzero >= 1                       # the ray really does cross a bracket zero
        println("effective-source-series H1 ray walk (a = 0.5, r₀ = 10, m = 20, λ = 0.1): ",
                nzero, " bracket sign changes; worst |switched − direct|/|direct| across ",
                "them = ", worst_at_zero)

        # (iii) SLOW: the validator's full bracket-zero hunt, adjudicated at BigFloat(256).
        if _ESS_SLOW
            setprecision(BigFloat, 256) do
                nz = 0
                nmis = 0
                wF = 0.0
                wabs = 0.0
                for c in EFFSOURCE_CONFIGS, mm in (0, 1, 2, 5, 10, 20)
                    oo = orbits[c.idx]
                    pp = puncs[c.idx]
                    oB = CircularOrbit(big(c.a), big(c.r0))
                    pB = Puncture(oB)
                    Am = effective_source_series_coefficients(pp, mm)
                    brk = function (dr, dth)
                        rho2 = oo.rho2_dr2 * dr^2 + oo.rho2_dtheta2 * dth^2
                        d2 = dr^2
                        lg = log(rho2)
                        Am.S2_2 * d2 + Am.S2_0 * rho2 + (Am.S2_2L * d2 + Am.S2_0L * rho2) * lg / 2 +
                        Am.S2_4 * d2^2 / rho2 + Am.S2_6 * d2^3 / rho2^2 + Am.S2_8 * d2^4 / rho2^3 +
                        Am.S2_10 * d2^5 / rho2^4 + Am.S2_12 * d2^6 / rho2^5
                    end
                    for LL in (1e-3, 3.1623e-3, 1e-2, 3.1623e-2, 1e-1, 2e-1, 3e-1)
                        pv = nothing
                        pt = -1.0
                        for t in range(-0.99999, 0.99999; length = 2001)
                            r, z = _ess_ray(oo, LL, t)
                            v = brk(r - oo.r0, -asin(z))
                            if pv !== nothing && sign(v) != sign(pv)
                                lo, hi = pt, t
                                for _ in 1:70
                                    mid = (lo + hi) / 2
                                    rm, zm = _ess_ray(oo, LL, mid)
                                    sign(brk(rm - oo.r0, -asin(zm))) == sign(pv) ? (lo = mid) : (hi = mid)
                                end
                                r, z = _ess_ray(oo, LL, (lo + hi) / 2)
                                truth = effective_source_m(pB, mm, big(r), big(z); phi_p = big(phi_p))
                                bb = effective_source_switched_bundle(pp, mm, r, z; phi_p = phi_p)
                                rs = Float64(abs(big(effective_source_series_m(pp, mm, r, z; phi_p = phi_p)) - truth) / abs(truth))
                                rd = Float64(abs(big(effective_source_m(pp, mm, r, z; phi_p = phi_p)) - truth) / abs(truth))
                                ra = Float64(abs(big(bb.source) - truth) / abs(truth))
                                nz += 1
                                wF = max(wF, ra / max(min(rs, rd), 1e-300))
                                wabs = max(wabs, ra)
                                (bb.route === :series && rs > 1e-3 && rd < 1e-6) && (nmis += 1)
                            end
                            pv = v
                            pt = t
                        end
                    end
                end
                @test nz >= 40                 # the hunt found the zeros it is meant to find
                @test nmis == 0                # ZERO misfires (pre-repair: 46 of 46)
                @test wF <= 5                  # measured 1.0
                @test wabs <= 1e-4             # measured 8.0642e-06 (pre-repair: 0.8704)
                println("effective-source-series H1 SLOW zero hunt: ", nz, " bracket zeros, ",
                        nmis, " misfires; worst achieved/optimal = ", wF,
                        " ; worst switched rel = ", wabs)
            end
        else
            @info "effective-source-series: H1 SLOW bracket-zero hunt skipped (~4 min); set GRTOOLS_SLOW_TESTS=1"
        end
    end

    @testset "the switch replayed on the BAKED 3240-point calibration grid" begin
        # M3: the headline switch numbers are ASSERTED against the baked table, not quoted
        # from a comment. The table holds, per grid point, the two routes' relative errors
        # against a BigFloat(256) truth and the direct route's eps·cancellation — all
        # route-INDEPENDENT, so replaying them against the live rule costs only Float64.
        cal = EFFSOURCE_SWITCH_CALIBRATION
        @test size(cal) == (3, 3240)
        worstF = 0.0
        worst_abs = 0.0
        worst_opt = 0.0
        worst_jump = 0.0
        nser = 0
        ndir = 0
        nsmall = 0
        k = 0
        for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 10, 20), e in -5.0:0.5:-1.0,
            frac in (0.0, 0.25, 0.5, 0.75, 1.0), sgn in (1, -1)

            k += 1
            o = orbits[c.idx]
            p = puncs[c.idx]
            L = 10.0^e
            dr = sgn * frac * L / sqrt(o.rho2_dr2)
            dth = sqrt(max(0.0, 1 - frac^2)) * L / sqrt(o.rho2_dtheta2)
            r = o.r0 + dr
            z = -sin(dth)
            rs, rd, ec = cal[1, k], cal[2, k], cal[3, k]
            b = effective_source_switched_bundle(p, m, r, z; phi_p = phi_p)
            # the baked eps·cancellation must still be the one the live direct route reports
            @test abs(b.err_direct / 1000 - ec) <= 1e-6 * ec
            pick = b.route === :series
            ra = pick ? rs : rd
            pick ? (nser += 1) : (ndir += 1)
            worstF = max(worstF, ra / max(min(rs, rd), 1e-300))
            worst_abs = max(worst_abs, ra)
            worst_opt = max(worst_opt, min(rs, rd))
            if 0.1 <= b.err_series / b.err_direct <= 10
                worst_jump = max(worst_jump, rs + rd)
            end
            if L <= 3.2e-4
                @test pick                     # 1080 points, all series
                nsmall += 1
            end
        end
        @test k == 3240
        @test nsmall == 1440              # 4 of the 9 λ decades × 360 points each
        @test nser > 0 && ndir > 0
        # THE HEADLINE NUMBERS (measured post-repair; pre-repair they were 42.106 / 6.9542e-04
        # / 5.8608e-04 / 1.1999e-03 — the repair changed NONE of them on this grid)
        @test worstF <= 50.0                   # measured 42.106
        @test worst_abs <= 1e-3                # measured 6.9542e-04
        @test worst_opt <= 1e-3                # measured 5.8608e-04
        @test worst_abs <= 1.5 * worst_opt     # measured ratio 1.187
        @test worst_jump <= 2e-3               # measured 1.1999e-03
        println("effective-source-series calibration replay (3240 pts): series ", nser,
                " / direct ", ndir, " ; worst achieved/optimal = ", worstF,
                " ; worst switched rel = ", worst_abs, " (oracle-optimal ", worst_opt,
                ", ratio ", worst_abs / worst_opt, ") ; worst near-crossover route ",
                "disagreement = ", worst_jump)
    end

    @testset "symmetry, precision inheritance, allocations" begin
        worstc = 0.0
        for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 20),
            (dr, dth) in ((1e-3, 0.0), (1e-4, 5e-4), (0.0, 1e-3), (-2e-3, 1e-3))

            o = orbits[c.idx]
            p = puncs[c.idx]
            r = o.r0 + dr
            z = -sin(dth)
            a = effective_source_series_m(p, m, r, z; phi_p = phi_p)
            b = effective_source_series_m(p, -m, r, z; phi_p = phi_p)
            @test b == conj(a)                                     # EXACT, |m| coefficients
            worstc = max(worstc, abs(b - conj(a)) / abs(a))
            if m == 0
                @test imag(a) == 0                                 # m = 0 reality, EXACT
                @test imag(effective_source_m_switched(p, 0, r, z; phi_p = phi_p)) == 0
            end
        end
        println("effective-source-series conjugation |S_{−m} − conj(S_m)|/|S_m| = ", worstc)

        setprecision(BigFloat, 256) do
            oB = CircularOrbit(big"0.5", big"10.0")
            pB = Puncture(oB)
            r = oB.r0 + big"1e-3"
            z = -sin(big"5e-4")
            v = effective_source_series_m(pB, 2, r, z; phi_p = big(phi_p))
            @test v isa Complex{BigFloat}
            @test precision(real(v)) == 256
            A = effective_source_series_coefficients(pB, 2)
            @test all(getfield(A, n) isa BigFloat for n in _ESS_PAPER_NAMES)
            @test effective_source_switched_bundle(pB, 2, r, z;
                                                   phi_p = big(phi_p)).source isa Complex{BigFloat}
        end

        # literal coordinates, as in the sealed direct-route allocation gate: a coordinate
        # read from a soft-scope loop variable is boxed by the testset and would measure the
        # box, not the callee. 10.001/0.0005 is inside the switch's series branch and
        # 10.3/0.03 inside its direct branch, so both are exercised.
        p = puncs[2]
        _ess_alloc_coeffs(p, 5)
        _ess_alloc_series(p, 5, 10.001, 0.0005, 0.31)
        _ess_alloc_switched(p, 5, 10.001, 0.0005, 0.31)
        _ess_alloc_switched(p, 5, 10.3, 0.03, 0.31)
        _ess_alloc_bundle(p, 5, 10.001, 0.0005, 0.31)
        a1 = _ess_alloc_coeffs(p, 5)
        a2 = _ess_alloc_series(p, 5, 10.001, 0.0005, 0.31)
        a3 = _ess_alloc_switched(p, 5, 10.001, 0.0005, 0.31)
        a4 = _ess_alloc_switched(p, 5, 10.3, 0.03, 0.31)
        a5 = _ess_alloc_bundle(p, 5, 10.001, 0.0005, 0.31)
        println("effective-source-series allocations (B): coefficients ", a1, ", series ", a2,
                ", switched[series branch] ", a3, ", switched[direct branch] ", a4,
                ", bundle ", a5)
        @test a1 == 0
        @test a2 == 0
        @test a3 == 0
        @test a4 == 0
        @test a5 == 0
        @test effective_source_switched_bundle(p, 5, 10.001, 0.0005;
                                               phi_p = 0.31).route === :series
        @test effective_source_switched_bundle(p, 5, 10.3, 0.03; phi_p = 0.31).route === :direct
    end
end

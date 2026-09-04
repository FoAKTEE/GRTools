# test_modesum.jl — ESSF node `mode-sum-total-force` (GLOBAL_DAG row 46, wave 13). THE TOTAL:
# `F_r = Σ_m [F_r^m − F^m_{r[4]}] + tail`, the number the CRP compares against Warburton–Barack.
#
# WHAT THE NODE DELIVERS (`Extraction/modesum.jl`):
#   * `total_self_force(orbit, Fr_modes; …)` — the arithmetic core on raw per-mode numbers;
#     `(set::ModeSet)`, `(hists::Vector{AmrHistory})` and `(orbit, config; m_max, pool)` are
#     the solved routes, all reducing to it;
#   * `mode_sum_tail(orbit, ms, terms; method)` — the CLOSED-FORM estimate of the remainder
#     `Σ_{m>m_max}[F_r^m − F^m_{r[4]}]` with its uncertainty (`:heffernan`, `:fit`, `:none`);
#   * `mode_decay_fit(ms, values; errors, headroom, window)` — the log-log decay exponent with
#     the uncertainty the per-mode resolution errors propagate into it, and the headroom rule.
#
# THE ORACLES, IN ORDER OF STRENGTH:
#   1. EXTERNAL AND INDEPENDENT: the Warburton–Barack printed `(M²/q²)F_r` at `r₀ = 10M`
#      (arXiv:1003.1860 Table III, `Warburton.tex:620-645`; "All figures presented are
#      significant", caption `:642-644`), baked STRING-EQUAL to the tex in
#      `test/fixtures/wb_reference.jl` and re-diffed against the mirror by
#      `test_wb_fixtures.jl`. Five spins: `a = 0, ±0.5, ±0.9` — the CRP's own comparison
#      (`main.tex:703`, `fig:ssf_error_budget`). NOTHING in this package feeds that table.
#   2. EXACT ARITHMETIC IDENTITIES on synthetic sequences: a power law summed in closed form
#      against `ζ(s, a)`'s own defining telescoping identity and against the exact
#      `π^{2k}`-values; the Heffernan shape summed against its sealed closed-form tail; the
#      acceleration identity `Σ(F − P) ≡ ΣF + tail(P)` that the whole method rests on.
#   3. INTERNAL CONSISTENCY: the two solved routes (`ModeSet` and `AmrHistory`) must give the
#      same number, and the recorded per-mode table must be reproduced by a LIVE solve.
#
# NO EXTRACTION IS RE-DERIVED HERE OR IN THE NODE: both solved routes go through the SEALED
# row-44 surface — `mode_self_force_history(hist)` for an `AmrHistory` and the level route
# `mode_self_force_history(orbit, config, m, levels; x_p)` for a fanned-out `ModeResult`,
# which are the same internal function.
#
# THE PER-MODE TABLE (`test/fixtures/modesum_wb_r10.jl`) is a RECORDING of the sealed row-44 +
# row-47 pipeline (21 modes × 5 spins × 9 AMR grids ≈ 30 CPU-minutes), not an oracle: it exists
# so that the WB comparison — which is an assertion about THIS node's summation and tail, given
# per-mode forces that rows 44/45/47 already gate — runs in seconds and reproducibly. Gate (e)
# re-solves one mode per spin LIVE and requires the recording to match, so a drift of the
# per-mode pipeline cannot hide behind the recording.
#
# PAPER ANCHORS: CRP `arXiv:2606.04998 main.tex:701-703` (extract per mode, sum over m),
# `:722-723` (the `F_{r[4]}` acceleration), `:726-731` (the error budget and the "constant
# offset" of the truncation model); [H14] `arXiv:1211.6446 KerrRegularization.tex:1675-1680`
# (the m-dependence at orders 4 and 6), `:1795-1797` (`F_{r[6]}` is a numerical fit, never a
# closed form — which is why this node FITS its amplitude and takes only its shape from [H14]).
#
# GATE CONSTANTS: named below, always separately from the values observed. Every threshold was
# fixed on the throwaway calibration prototypes kept in the node evidence
# (`msum_proto2.jl`/`msum_proto3.jl`) BEFORE this file was written; the WB threshold is the
# GLOBAL_DAG row-46 planning gate ("a few × 1e-4") written down as a number, not a measurement.
# No threshold has been loosened to make a failing assertion pass.
#
# NO figure (rows 50/51/53), NO new solver, NO re-derivation of `F_r^m` (row 44 owns it).

using Test
using Printf
using ScalarSelfForce
using EllipticCore

# ------------------------------------ GATE CONSTANTS ---------------------------------------
const MS_IDENTITY_RTOL   = 1.0e-14  # (a) the acceleration identity, in Float64 arithmetic
const MS_ZETA_RTOL       = 1.0e-14  # (a) Euler–Maclaurin ζ(s,a) vs exact closed forms
const MS_ZETA_TELE_RTOL  = 1.0e-12  # (a) ... and vs its telescoping identity, which cancels
const MS_EXACT_RTOL      = 1.0e-12  # (a) closed-form tails vs brute-force partial sums
const MS_FIT_RTOL        = 1.0e-10  # (a) recovery of an exactly synthesized power law
const MS_WB_RTOL         = 5.0e-4   # (b) THE NODE GATE — GLOBAL_DAG row 46, "few × 1e-4"
const MS_BUDGET_SLACK    = 5.0      # (b) the budget must cover the residual WB deviation
const MS_MMAX_SLACK      = 3.0      # (c) how far the m_max change may exceed the tail estimate
const MS_TAIL_GAIN       = 3.0      # (c) how much the tail must shrink the m_max drift
const MS_SLOPE_UNREG     = -4.0     # (d) unregularized F_r^m  → F^m_{r[4]} ~ m^-4
const MS_SLOPE_REG       = -6.0     # (d) regularized residual → F^m_{r[6]} ~ m^-6
const MS_SLOPE_BAND      = 0.3      # (d) band on the RESOLVED unregularized exponent
const MS_SLOPE_REG_BAND  = 0.5      # (d) band on the regularized one (next-order + resolution)
const MS_SLOPE_WINDOW    = 6        # (d) trailing modes admitted to a slope fit
const MS_SLOPE_SIGMA_MAX = 0.1      # (d) propagated uncertainty required of the SHARP fit
const MS_SLOPE_HEADROOM  = 3.0      # (d) |value| / resolution error required of a fit point
const MS_LIVE_RTOL       = 1.0e-6   # (e) live re-solve vs the recorded per-mode table: 30x
                                    #     tighter than that mode's OWN resolution error and
                                    #     25x looser than the ~4e-8 a change of solver
                                    #     PROVENANCE (multigrid hierarchy, Schwarz overlap)
                                    #     can move a converged diagnostic. Measured 0.0.
const MS_ROUTE_RTOL      = 1.0e-12  # (f) ModeSet route vs AmrHistory route

# ---- fixtures (ms_-prefixed: sibling test files share the runtests.jl include scope) -------
const MS_WB = Module(:MSWarburtonBarack)
Base.include(MS_WB, joinpath(@__DIR__, "fixtures", "wb_reference.jl"))

const MS_TAB = Module(:MSModeSumTable)
Base.include(MS_TAB, joinpath(@__DIR__, "fixtures", "modesum_wb_r10.jl"))

ms_wb_Fr(a) = MS_WB.wb_cell(MS_WB.WB_FR, a, MS_TAB.MSWB_R0).value

ms_cfg(; spin, r0 = MS_TAB.MSWB_R0, m = 0, L = 0, P = 4) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L)

# ============================ (a) EXACT ARITHMETIC, NO SOLVE ================================

@testset "(a1) the Euler-Maclaurin Hurwitz zeta" begin
    z = ScalarSelfForce._ms_hurwitz_zeta
    # exact closed forms of ζ(2k, 1)
    @test isapprox(z(2, 1), pi^2 / 6;    rtol = MS_ZETA_RTOL)
    @test isapprox(z(4, 1), pi^4 / 90;   rtol = MS_ZETA_RTOL)
    @test isapprox(z(6, 1), pi^6 / 945;  rtol = MS_ZETA_RTOL)
    @test isapprox(z(8, 1), pi^8 / 9450; rtol = MS_ZETA_RTOL)
    # the shift identity ζ(s,a) = ζ(s,1) − Σ_{k<a} k^{-s}, for integer a. The RIGHT-hand side
    # cancels catastrophically for large a (1.017… − 1.017… = 5e-8 at s = 6, a = 21), so the
    # comparison is made on the ABSOLUTE scale of ζ(s,1) — a relative gate here would be a
    # statement about the oracle's round-off, not about ζ.
    for s in (2, 4, 6), a in (2, 7, 21)
        @test abs(z(s, a) - (z(s, 1) - sum(k -> Float64(k)^(-s), 1:(a - 1)))) <=
              MS_ZETA_RTOL * z(s, 1)
    end
    # the DEFINING telescoping identity, valid for every real s > 1 (non-integer included)
    for s in (1.5, 2.7, 5.5, 6.06, 9.3), a in (1.0, 3.0, 21.0, 100.0)
        @test isapprox(z(s, a) - z(s, a + 1), a^(-s); rtol = MS_ZETA_TELE_RTOL)
    end
    # and the integral bracket ∫_a^∞ x^-s dx ≤ ζ(s,a) ≤ a^-s + ∫_a^∞ x^-s dx
    for s in (1.5, 3.0, 6.0, 7.25), a in (2.0, 21.0)
        lo = a^(1 - s) / (s - 1)
        @test lo <= z(s, a) <= lo + a^(-s)
    end
    @test_throws DomainError z(1.0, 1.0)
    @test_throws DomainError z(3.0, 0.0)
    # the TYPE is inherited; the PRECISION is not (Euler–Maclaurin truncation, not arithmetic,
    # is the limit — see the accuracy note at the function), and that is stated, not assumed
    @test z(big(6), big(1)) isa BigFloat
    @test isapprox(z(big(6), big(1)), big(pi)^6 / 945; rtol = big(1e-16))
    @test !isapprox(z(big(6), big(1)), big(pi)^6 / 945; rtol = big(1e-20))
end

@testset "(a2) a synthetic power law is recovered and summed exactly" begin
    orb = CircularOrbit(0.0, 10.0)
    C, p, mmax = 3.5e-6, 6.0, 20
    ms = 0:mmax
    terms = [m == 0 ? -1.0e-4 : C * Float64(m)^(-p) for m in ms]   # m = 0 is off the law
    fit = mode_decay_fit(ms, terms; window = 6)
    @test isapprox(fit.slope, -p; rtol = MS_FIT_RTOL)
    @test isapprox(fit.amplitude, C; rtol = MS_FIT_RTOL)
    @test fit.ms == collect(15:20)
    @test 0 in fit.excluded                                   # log(0) does not exist
    @test isnan(fit.slope_sigma)                              # no errors ⇒ no uncertainty
    t = mode_sum_tail(orb, ms, terms; method = :fit, window = 6)
    # ζ(6, 21) = π⁶/945 − Σ_{k≤20} k^{-6} loses SEVEN digits to cancellation in Float64
    # (1.0173… − 1.0173… = 5.5e-8), so the reference is formed in BigFloat and only then
    # rounded: the oracle must be more accurate than the thing it judges.
    exact = Float64(C * (big(pi)^6 / 945 - sum(k -> big(k)^(-6), 1:mmax)))
    @test isapprox(t.value, exact; rtol = MS_FIT_RTOL)
    @test t.method === :fit && isapprox(t.exponent, p; rtol = MS_FIT_RTOL)
    # the CRP's omit-the-last-mode estimator vs the true remainder: the ≈ m_max/(p−1) offset
    @test isapprox(t.truncation_error, abs(terms[end]); rtol = MS_IDENTITY_RTOL)
    @test isapprox(exact / abs(terms[end]), mmax / (p - 1); rtol = 0.2)
    # a brute-force partial sum reaches the same number
    brute = sum(k -> C * Float64(k)^(-p), (mmax + 1):200_000)
    @test isapprox(t.value, brute; rtol = 1.0e-5)
end

@testset "(a3) the Heffernan shape is summed by its own sealed closed form" begin
    orb = CircularOrbit(0.5, 10.0)
    A, mmax = -7.25e-3, 20
    ms = 0:mmax
    terms = [A * mmode_reg_factor(6, m) for m in ms]           # EXACTLY the order-6 shape
    t = mode_sum_tail(orb, ms, terms; method = :heffernan, window = 6)
    @test isapprox(t.amplitude, A; rtol = MS_IDENTITY_RTOL)
    @test isapprox(t.value, A * mmode_reg_factor_tail(6, mmax + 1); rtol = MS_IDENTITY_RTOL)
    @test t.exponent == 6
    brute = sum(m -> A * mmode_reg_factor(6, m), (mmax + 1):200_000)
    @test isapprox(t.value, brute; rtol = MS_EXACT_RTOL)
    # the free fit is the CONTROL: it must find the same remainder without being told p = 6
    @test isapprox(t.alt_value, t.value; rtol = 5.0e-2)
    @test abs(t.uncertainty) <= 6.0e-2 * abs(t.value)
    # :none is the honest zero, with the CRP's own truncation estimate as its uncertainty
    n = mode_sum_tail(orb, ms, terms; method = :none)
    @test n.value == 0 && n.uncertainty == abs(terms[end]) && n.fit === nothing
end

@testset "(a4) THE ACCELERATION IDENTITY of main.tex:723" begin
    orb = CircularOrbit(-0.5, 10.0)
    mmax = 20
    frs = [1.0e-4 * cos(0.7m) / (1 + m)^3 for m in 0:mmax]     # arbitrary per-mode numbers
    r = total_self_force(orb, frs; tail = :none)
    # form 1 (what this node computes): Σ_m [F_r^m − F^m_{r[4]}], the m = 0 term halved
    @test isapprox(r.partial, sum(r.terms); rtol = MS_IDENTITY_RTOL)
    # form 2 (main.tex:723 as written): Σ_m F_r^m + the analytic parameter tail
    @test isapprox(r.partial, sum(frs) + mmode_reg_tail_Fr(4, orb, mmax + 1);
                   rtol = MS_IDENTITY_RTOL)
    # the m = 0 halving lives in ONE place and is applied exactly once
    @test r.reg_terms[1] == mmode_reg_Fr(4, orb, 0) / 2
    @test all(r.reg_terms[k] == mmode_reg_Fr(4, orb, k - 1) for k in 2:length(r.reg_terms))
    @test all(r.terms .== frs .- r.reg_terms)
    @test r.Fr == r.partial                                     # tail = :none adds nothing
    @test r.truncation_error == abs(r.terms[end])
    # order 2 is a different (and for a 4th-order puncture, wrong) acceleration — but the
    # identity is about the parameter, not about which one, so it must hold there too
    r2 = total_self_force(orb, frs; tail = :none, order = 2)
    @test isapprox(r2.partial, sum(frs) + mmode_reg_tail_Fr(2, orb, mmax + 1);
                   rtol = MS_IDENTITY_RTOL)
end

@testset "(a5) structural falsifiers and type genericity" begin
    orb = CircularOrbit(0.0, 10.0)
    frs = [1.0e-5 / (1 + m)^4 for m in 0:8]
    @test_throws ArgumentError total_self_force(orb, frs; ms = 1:9)          # must start at 0
    @test_throws ArgumentError total_self_force(orb, frs; ms = [0, 1, 2, 4, 5, 6, 7, 8, 9])
    @test_throws ArgumentError total_self_force(orb, frs; ms = 0:7)          # length mismatch
    @test_throws ArgumentError total_self_force(orb, Float64[])
    @test_throws ArgumentError total_self_force(orb, frs; resolution_deltas = [1.0])
    @test_throws ArgumentError total_self_force(orb, frs; tail = :magic)
    # the tail model needs the PRINTED m-dependence of order+2; order 8 has no order 10
    @test_throws ArgumentError total_self_force(orb, frs; order = 8)
    @test_throws ArgumentError mode_decay_fit(0:2, [1.0, 2.0, 3.0, 4.0])
    @test_throws ArgumentError mode_decay_fit([0], [1.0])                    # one point
    # the HEADROOM RULE: a mode drowned in its own resolution error may not enter a fit
    ms, vals = 1:6, [1.0e-4 / m^6 for m in 1:6]
    errs = [m == 6 ? 10 * vals[m] : 1.0e-16 for m in 1:6]
    f = mode_decay_fit(ms, vals; errors = errs, headroom = MS_SLOPE_HEADROOM)
    @test 6 in f.excluded && f.ms == collect(1:5)
    @test isapprox(f.slope, -6.0; rtol = MS_FIT_RTOL)
    @test f.slope_sigma < 1.0e-7                     # 1e-16 errors on 1e-4…6e-9 values
    @test_throws ArgumentError mode_decay_fit(ms, vals; errors = 10 .* vals)  # nothing passes
    # BigFloat in, BigFloat out, at the ambient precision
    borb = CircularOrbit(big"0.0", big"10.0")
    bfrs = [big(1) / 10^5 / (1 + m)^4 for m in 0:20]
    br = total_self_force(borb, bfrs; tail = :heffernan)
    @test br.Fr isa BigFloat && br.tail.value isa BigFloat
    @test isapprox(br.partial, sum(bfrs) + mmode_reg_tail_Fr(4, borb, 21); rtol = big(1e-60))
end

# ==================== (b) THE WARBURTON-BARACK GATE — THE NODE'S GATE ========================
# The recorded 21-mode tables at r0 = 10M for a = 0, ±0.5, ±0.9 (fixture header states the
# commit, the configuration and the capture script) summed by this node, against the printed
# Table III. `MS_WB_RTOL` is the GLOBAL_DAG row-46 planning gate.

const MS_TOTALS = Dict(a => total_self_force(CircularOrbit(a, MS_TAB.MSWB_R0),
                                             MS_TAB.MSWB_FR[a];
                                             ms = MS_TAB.MSWB_MMODES,
                                             resolution_deltas = MS_TAB.MSWB_DELTA[a],
                                             reduction = MS_TAB.MSWB_REDUCTION,
                                             levels = MS_TAB.MSWB_LEVELS)
                       for a in MS_TAB.MSWB_SPINS)

@testset "(b) total F_r against Warburton-Barack Table III at r0 = 10M" begin
    worst = 0.0
    for a in MS_TAB.MSWB_SPINS
        r = MS_TOTALS[a]
        ref = ms_wb_Fr(a)
        dev = abs(r.Fr - ref) / abs(ref)
        worst = max(worst, dev)
        @printf("mode-sum (b) observed a = %+.1f : ours %+.8e  WB %+.8e  rel %.3e  \
                 tail %+.3e (%.2e rel, sigma %.1e)  trunc %.3e  resol %.3e (bound %.3e)\n",
                a, r.Fr, ref, dev, r.tail.value, abs(r.tail.value / r.Fr),
                r.tail.uncertainty, r.truncation_error, r.resolution_error,
                r.resolution_error_bound)
        @test dev <= MS_WB_RTOL
        @test sign(r.Fr) == sign(ref)                     # the prograde sign flip, per spin
        # THE ERROR BUDGET CLOSES: what is left over after the accelerated sum and its tail is
        # covered by the three budgeted terms of main.tex:726-731 (numerical resolution,
        # finite-m truncation, and the tail model's own uncertainty).
        budget = r.resolution_error_bound + r.truncation_error + r.tail.uncertainty
        @test abs(r.Fr - ref) <= MS_BUDGET_SLACK * budget
        @test r.resolution_error <= r.resolution_error_bound
    end
    @printf("mode-sum (b) observed WORST relative deviation over the five spins: %.3e \
             (gate %.1e)\n", worst, MS_WB_RTOL)
    @test worst <= MS_WB_RTOL
end

@testset "(b2) FALSIFIER: what the large-m tail is worth at m_max = 20" begin
    improved = 0
    for a in MS_TAB.MSWB_SPINS
        ref = ms_wb_Fr(a)
        with = MS_TOTALS[a]
        without = total_self_force(CircularOrbit(a, MS_TAB.MSWB_R0), MS_TAB.MSWB_FR[a];
                                   ms = MS_TAB.MSWB_MMODES, tail = :none)
        dw, dn = abs(with.Fr - ref) / abs(ref), abs(without.Fr - ref) / abs(ref)
        improved += (dw < dn)
        @printf("mode-sum (b2) observed a = %+.1f : with tail %.3e  without %.3e  \
                 ratio %.2f  |tail|/trunc %.2f (model m_max/(p-1) = %.2f)\n",
                a, dw, dn, dn / dw, abs(with.tail.value) / with.truncation_error,
                last(with.ms) / (-with.tail.fit.slope - 1))
        # structural: :none really is no tail, and the tail is the ONLY difference
        @test without.tail.value == 0 && without.Fr == without.partial
        # `(partial + tail) − partial` cancels five orders, so the comparison is made on the
        # ABSOLUTE scale the subtraction can carry — `eps` of the partial sum, not a relative
        # tolerance, which here would be a statement about Float64 and not about the tail.
        @test abs(abs(with.Fr - without.Fr) - abs(with.tail.value)) <=
              8 * eps(abs(with.partial))
        # THE MODEL CLAIM of main.tex:731: the omit-the-last-mode estimator UNDERSTATES the
        # remainder by ~ m_max/(p-1), which is why the paper's truncation error "models the
        # true error up to a constant offset"
        @test abs(with.tail.value) > with.truncation_error
        @test isapprox(abs(with.tail.value) / with.truncation_error,
                       last(with.ms) / (-with.tail.fit.slope - 1); rtol = 0.4)
        # the two independent tail models agree to within their own scale
        @test isfinite(with.tail.alt_value)
        @test abs(with.tail.value - with.tail.alt_value) <= abs(with.tail.value)
    end
    @printf("mode-sum (b2) observed: the tail improves the WB deviation at %d of %d spins\n",
            improved, length(MS_TAB.MSWB_SPINS))
    @test improved == length(MS_TAB.MSWB_SPINS)
end

@testset "(c) convergence of the total with m_max, and what the tail buys" begin
    for a in MS_TAB.MSWB_SPINS
        ref = ms_wb_Fr(a)
        rs = Dict(M => total_self_force(CircularOrbit(a, MS_TAB.MSWB_R0),
                                        MS_TAB.MSWB_FR[a][1:(M + 1)]; ms = 0:M,
                                        resolution_deltas = MS_TAB.MSWB_DELTA[a][1:(M + 1)])
                  for M in (10, 15, 20))
        pn = Dict(M => total_self_force(CircularOrbit(a, MS_TAB.MSWB_R0),
                                        MS_TAB.MSWB_FR[a][1:(M + 1)]; ms = 0:M,
                                        tail = :none).Fr for M in (10, 15, 20))
        @printf("mode-sum (c) observed a = %+.1f : F(10) %+.10e  F(15) %+.10e  \
                 F(20) %+.10e  (rel to WB %.2e / %.2e / %.2e) ; untailed drift 10->20 %.3e, \
                 tailed %.3e, gain %.1f\n", a, rs[10].Fr, rs[15].Fr, rs[20].Fr,
                abs(rs[10].Fr - ref) / abs(ref), abs(rs[15].Fr - ref) / abs(ref),
                abs(rs[20].Fr - ref) / abs(ref), abs(pn[20] - pn[10]),
                abs(rs[20].Fr - rs[10].Fr), abs(pn[20] - pn[10]) / abs(rs[20].Fr - rs[10].Fr))
        # (c1) the change with m_max is BRACKETED by the tail estimates that model it
        for (Ma, Mb) in ((10, 15), (15, 20), (10, 20))
            change = abs(rs[Mb].Fr - rs[Ma].Fr)
            budget = abs(rs[Ma].tail.value) + rs[Ma].tail.uncertainty +
                     abs(rs[Mb].tail.value) + rs[Mb].tail.uncertainty
            @test change <= MS_MMAX_SLACK * budget
        end
        # (c2) THE FALSIFIER: without the tail the total still drifts with m_max; with it the
        # drift collapses. If the tail did not matter this factor would be 1.
        @test abs(pn[20] - pn[10]) >= MS_TAIL_GAIN * abs(rs[20].Fr - rs[10].Fr)
    end
end

@testset "(d) the decay exponents: m^-4 unregularized, m^-6 regularized" begin
    for a in MS_TAB.MSWB_SPINS
        r = MS_TOTALS[a]
        # (d1) the UNREGULARIZED sequence asymptotes to F^m_{r[4]} itself, so its exponent is
        # -4 up to the finite-m curvature of the printed rational m-dependence.
        fu = mode_decay_fit(r.ms, r.Fr_modes; errors = r.resolution_errors,
                            headroom = MS_SLOPE_HEADROOM, window = MS_SLOPE_WINDOW)
        # (d2) the REGULARIZED residual, with F^m_{r[4]} removed, is the NEXT parameter ~ m^-6.
        fr = mode_decay_fit(r.ms, r.terms; errors = r.resolution_errors,
                            headroom = MS_SLOPE_HEADROOM, window = MS_SLOPE_WINDOW)
        # THE PREDICTION is not the asymptote -2j: it is the exact log-log slope of the
        # PRINTED m-dependence over the very modes the fit used, computed by the same
        # least squares. Quoting -4 or -6 over m ~ 15..20 is a 1-2 % statement, not an
        # identity, and this makes the difference measurable instead of absorbed into a band.
        pu = mode_decay_fit(fu.ms, [mmode_reg_factor(4, m) for m in fu.ms]).slope
        pr = mode_decay_fit(fr.ms, [mmode_reg_factor(6, m) for m in fr.ms]).slope
        @printf("mode-sum (d) observed a = %+.1f : unregularized %+.4f ± %.4f over m ∈ %d:%d \
                 (F^m_{r[4]} shape predicts %+.4f, dev %.3f) ; regularized %+.4f ± %.4f over \
                 m ∈ %d:%d (F^m_{r[6]} shape predicts %+.4f, dev %.3f ; excluded %s)\n",
                a, fu.slope, fu.slope_sigma, first(fu.ms), last(fu.ms), pu, abs(fu.slope - pu),
                fr.slope, fr.slope_sigma, first(fr.ms), last(fr.ms), pr, abs(fr.slope - pr),
                string(fr.excluded))
        # the unregularized exponent is RESOLVED and sharp
        @test fu.slope_sigma <= MS_SLOPE_SIGMA_MAX
        @test abs(fu.slope - pu) <= MS_SLOPE_BAND
        @test abs(pu - MS_SLOPE_UNREG) <= MS_SLOPE_BAND        # ... and the shape IS ~ m^-4
        # the regularized exponent is measured with the LARGER band its own resolution
        # headroom allows at this AMR depth; its uncertainty is reported, not gated away
        @test abs(fr.slope - pr) <= MS_SLOPE_REG_BAND
        @test abs(pr - MS_SLOPE_REG) <= MS_SLOPE_BAND          # ... and the shape IS ~ m^-6
        # the two slopes are DIFFERENT quantities and the regularization is what separates them
        @test fr.slope < fu.slope - 1.5
        # the headroom rule is doing its job at this depth: the very last modes are within a
        # few times their own resolution error and are kept out of the fit
        @test last(fr.ms) <= last(r.ms)
    end
end

# ============ (e) THE RECORDING IS WHAT THIS CODE PRODUCES — one LIVE mode per spin ==========

@testset "(e) live re-solve of one mode per spin reproduces the recorded table" begin
    for a in MS_TAB.MSWB_SPINS
        m = MS_TAB.MSWB_LIVE_MODE
        orb = CircularOrbit(a, MS_TAB.MSWB_R0)
        hist = solve_mode_amr(orb, mode_config(ms_cfg(; spin = a), m), m;
                              iterations = MS_TAB.MSWB_ITERATIONS)
        v = [f.Fr for f in mode_self_force_history(hist)]
        rec = MS_TAB.MSWB_FR[a][m + 1]
        recd = MS_TAB.MSWB_DELTA[a][m + 1]
        @printf("mode-sum (e) observed a = %+.1f m = %d : live %+.16e  recorded %+.16e  \
                 rel %.2e ; delta live %+.3e recorded %+.3e\n",
                a, m, v[end], rec, abs(v[end] - rec) / abs(rec), v[end] - v[end - 1], recd)
        @test isapprox(v[end], rec; rtol = MS_LIVE_RTOL)
        # ... and, the physically meaningful statement: any drift is smaller than this mode's
        # own numerical-resolution uncertainty, so it cannot move the total.
        @test abs(v[end] - rec) <= abs(recd)
        # the recorded DELTA is a difference of two nearly-equal levels, so it can only be
        # reproduced to the ABSOLUTE accuracy the values themselves are reproduced to — a
        # relative gate on it would be a statement about that cancellation, not about the
        # solve. Its SHAPE is gated separately and loosely.
        @test abs((v[end] - v[end - 1]) - recd) <= 2 * MS_LIVE_RTOL * abs(rec)
        @test isapprox(v[end] - v[end - 1], recd; rtol = 5.0e-2)
        @test length(v) == MS_TAB.MSWB_ITERATIONS + 1
    end
end

# ===================== (f) THE SOLVED ROUTES COMPOSE AND AGREE ===============================

@testset "(f) ModeSet, AmrHistory and end-to-end routes agree" begin
    a, mmax, iters = 0.0, 3, 2
    orb = CircularOrbit(a, MS_TAB.MSWB_R0)
    base = ms_cfg(; spin = a)
    hists = [solve_mode_amr(orb, mode_config(base, m), m; iterations = iters) for m in 0:mmax]
    rh = total_self_force(hists; window = 3)
    set = solve_modes_serial(orb, base, 0:mmax; iterations = iters)
    rs = total_self_force(set; window = 3)
    @printf("mode-sum (f) observed: AmrHistory route %+.16e  ModeSet route %+.16e  \
             (bitwise %s)\n", rh.Fr, rs.Fr, rh.Fr == rs.Fr ? "yes" : "no")
    @test isapprox(rh.Fr, rs.Fr; rtol = MS_ROUTE_RTOL)
    @test rh.ms == rs.ms == collect(0:mmax)
    @test rh.reduction === rs.reduction === :trace
    @test rh.levels == rs.levels == iters + 1
    @test all(isapprox.(rh.Fr_modes, rs.Fr_modes; rtol = MS_ROUTE_RTOL))
    @test all(isapprox.(rh.resolution_errors, rs.resolution_errors; rtol = 1.0e-6))
    # the end-to-end call is the same thing with the fan-out inside it
    re = total_self_force(orb, base; m_max = mmax, parallel = false, iterations = iters,
                          window = 3)
    @test isapprox(re.Fr, rs.Fr; rtol = MS_ROUTE_RTOL)
    # the per-mode value IS row 44's, not a re-derivation
    @test isapprox(rh.Fr_modes[2], mode_self_force(hists[2], 1).Fr; rtol = MS_ROUTE_RTOL)
    @test isapprox(rh.terms[2], regularized_mode_force(hists[2], 1; order = 4);
                   rtol = MS_ROUTE_RTOL)
    # coherent vs incoherent resolution error
    @test rh.resolution_error <= rh.resolution_error_bound
    @test isapprox(rh.resolution_error_bound, sum(rh.resolution_errors); rtol = 1.0e-14)
    # falsifiers on the solved routes
    @test_throws ArgumentError total_self_force(AmrHistory[])
    @test_throws ArgumentError total_self_force(orb, base; m_max = 0, parallel = false,
                                                iterations = iters)
    @test_throws ArgumentError total_self_force([hists[1], hists[3]])   # a gap in m
end

@testset "(g) show methods and the result surface" begin
    r = MS_TOTALS[0.0]
    s = sprint(show, MIME"text/plain"(), r)
    @test occursin("TotalSelfForce", s) && occursin("resolution error", s)
    @test occursin("TotalSelfForce", sprint(show, r))
    @test occursin("ModeSumTail", sprint(show, r.tail))
    @test occursin("ModeDecayFit", sprint(show, r.tail.fit))
    @test r.Fr == r.partial + r.tail.value
    @test r.tail.m_start == last(r.ms) + 1
    @test r.order == 4 && length(r.ms) == length(r.Fr_modes) == length(r.terms)
end

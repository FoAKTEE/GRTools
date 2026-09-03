# test_regularization.jl — ESSF node `heffernan-regularization` (schedule W7-E, GLOBAL_DAG
# row 45), written BEFORE the implementation (TDD). The node delivers the m-mode
# regularization parameters of the RADIAL scalar self-force for circular equatorial Kerr
# geodesics — `F^m_{r[2]}`, `F^m_{r[4]}` — plus the analytic large-m tail completion that
# `mode-sum-total-force` (row 46) sums against.
#
# SOURCES (read-only mirrors; every literal below cites file:line):
#   [H14] A. Heffernan, A. Ottewill, B. Wardell, "High-order expansions of the
#         Detweiler-Whiting singular field in Kerr spacetime", PRD 89, 024030 (2014),
#         arXiv:1211.6446 — ref-paper/arxiv-1211.6446/src/KerrRegularization.tex:
#           :1675-1680  the m-dependence gather F^m_{a[2j]} (orders 2, 4, 6, 8);
#           :1689-1699  F_{r[2]} for a circular orbit;
#           :1700-1731  F_{r[4]} for a circular orbit;
#           :1610-1671  the DERIVATION (folded m >= 0 modes, l.1583-1589);
#           :1784-1798  the a = 0.6M, r0 = 10M example (m^-6 residual after subtraction).
#   [TW]  B. Wardell, J. Thornburg, ancillary notebook
#         `m-ModeKerrRegularizationParameters.nb` of arXiv:1610.09319 — the ELECTRONIC form
#         [H14] points to at :1685-1687. Its `Fmr[4]` is F^m_{r[4]} for a GENERIC equatorial
#         Kerr geodesic in (E, L, r, a, M, m); specialised to circular orbits it is an
#         INDEPENDENT SOURCE, baked as fixtures/heffernan_mmode_reg.jl.
#   [CRP] arXiv:2606.04998 main.tex:520-522 (m^-4 mode-sum convergence, "accelerated using
#         analytically derived ... regularization parameters [Heffernan_2014]") and
#         main.tex:723 ("We use the F_{r,[4]} regularization parameter ... to accelerate the
#         sum over m-modes") — the CRP's 4th-order puncture has ALREADY subtracted the [2]
#         parameter through the effective source ([H14]:1786-1789), so [4] is the leading
#         unsubtracted one.
#
# GATES (GLOBAL_DAG row 45 oracle column):
#  (1) TRANSCRIPTION RE-DIFF — the tex brackets are transcribed TWICE: the shipped form is
#      the paper's FACTORED display (src/Extraction/regularization.jl); this file carries an
#      independently hand-EXPANDED monomial table (26 / 84 monomials, coeff * a^i M^j r^k
#      s^e with s = sqrt(M r)). Compared EXACTLY in Rational{BigInt} at rational (a, r0, M)
#      with s rational, and at Float64 / BigFloat(256). Cross-SOURCE: the assembled
#      F^m_{r[4]} vs the [TW] notebook oracle.
#  (2) a -> 0 SCHWARZSCHILD LIMIT vs closed forms in a third structural form, at
#      r0 in {6, 7, 10, 20}. FINDING (this worker, 2026-09-02): the DAG row cites
#      arXiv:1204.0794 (Heffernan-Ottewill-Wardell 2012, Paper I) for the Schwarzschild
#      m-mode parameters, but that paper has NO m-mode section at all (zero occurrences of
#      "m-mode"; its sec. IV Regularization.tex:893-1846 is l-mode mode-sum only). The
#      m-mode scheme is introduced in Paper II = [H14] sec. V. The independent Schwarzschild
#      statement used here is therefore [TW]'s generic-orbit `Fmr[4]` evaluated at a = 0
#      (gate 1's oracle rows with a = 0) plus the reduced closed forms below.
#  (3) STRUCTURE — the printed m-pattern (asserted, not assumed), the m-factor re-derived
#      EXACTLY from its defining Fourier integral by the multiple-angle expansion, the m = 0
#      folding subtlety, retrograde a < 0 and continuity through a = 0, BigFloat inheritance,
#      and the closed-form tail against brute-force partial sums.
#  (4) [FUTURE] the large-m residual slope (-6 +/- 0.3 after subtracting F^m_{r[4]} from real
#      F_r^m data) needs node 44 `mode-self-force-extraction`, which does not exist yet. What
#      IS gated here is a MANUFACTURED order check: a synthetic sequence built as
#      F^m_{r[4]} + C * (order-6 m-shape) must, after subtraction of the [4] term, have
#      log-log slope -6 within 0.3 — i.e. the subtraction removes exactly the m^-4 content.
using Test
using ScalarSelfForce
using SHA: sha256

# TestKit supplies the repo-wide oracle-availability gate (oracle_root / oracle_or_skip), but
# it is a TEST-ONLY dep that is on the load path only under `Pkg.test`. This file must ALSO run
# standalone under `--project=capabilities/ScalarSelfForce` (the ledger's verification command),
# so import it when present and otherwise fall back to a shim with IDENTICAL semantics —
# `Broken(:skipped)` on absence, hard failure under `GRTOOLS_STRICT_ORACLES=1` — rather than
# silently dropping the gate. Under `Pkg.test` the real TestKit always wins.
const W7E_HAVE_TESTKIT = try
    @eval using TestKit: oracle_root, oracle_or_skip
    true
catch
    false
end
if !W7E_HAVE_TESTKIT
    w7e_strict_oracles() =
        lowercase(strip(get(ENV, "GRTOOLS_STRICT_ORACLES", ""))) in ("1", "true", "yes", "on")
    oracle_root(envvar::AbstractString, default::AbstractString) =
        String(get(ENV, envvar, default))
    function oracle_or_skip(available::Bool, what::AbstractString; env = nothing)
        available && return true
        hint = env === nothing ? "" : " (set ENV[\"$env\"])"
        msg = string(what, " — oracle absent", hint)
        if w7e_strict_oracles()
            @test msg == "no missing oracle (GRTOOLS_STRICT_ORACLES=1)"
        else
            Test.record(Test.get_testset(), Test.Broken(:skipped, msg))
        end
        return false
    end
end

# ---------------------------------------------------------------- [TW] notebook oracle
include(joinpath(@__DIR__, "fixtures", "heffernan_mmode_reg.jl"))

# ============================================================ independent transcriptions
# (a) EXPANDED monomial re-transcription of the two tex brackets. Typed by hand from
#     KerrRegularization.tex:1691-1698 ([2]) and :1702-1730 ([4]) by multiplying out the
#     printed factored groups — a DIFFERENT algebraic form of the same coefficients, so a
#     disagreement with src/ is a transcription defect on one side or the other.
#     Entries are (coeff, i, j, k, e): coeff * a^i * M^j * r0^k * s^e, s = sqrt(M r0).
#     Every entry satisfies the dimensional-homogeneity invariant i + j + k + e == 9 ([2])
#     and == 17 ([4]) — asserted below, so a mistyped exponent cannot pass silently.
const W7E_EXPANDED = Dict(
 2 => ((24,7,2,0,0),
       (-24,6,2,0,1),(48,6,1,1,1),
       (-92,5,3,1,0),(-4,5,2,2,0),(24,5,1,3,0),
       (90,4,3,1,1),(-224,4,2,2,1),(62,4,1,3,1),
       (90,3,4,2,0),(90,3,3,3,0),(-146,3,2,4,0),(38,3,1,5,0),
       (-87,2,4,2,1),(264,2,3,3,1),(-114,2,2,4,1),(12,2,1,5,1),(-3,2,0,6,1),
       (-174,1,4,4,0),(258,1,3,5,0),(-126,1,2,6,0),(18,1,1,7,0),
       (-87,0,3,5,1),(75,0,2,6,1),(-9,0,1,7,1),(-3,0,0,8,1)),
 4 => ((-23040,14,2,0,1),
       (11520,13,3,1,0),(-92160,13,2,2,0),
       (177024,12,3,1,1),(-31104,12,2,2,1),(-138240,12,1,3,1),
       (-58944,11,4,2,0),(725760,11,3,3,0),(-236736,11,2,4,0),(-92160,11,1,5,0),
       (-547136,10,4,2,1),(229952,10,3,3,1),(1037568,10,2,4,1),(-405504,10,1,5,1),
       (-23040,10,0,6,1),
       (90720,9,5,3,0),(-2220832,9,4,4,0),(1490080,9,3,5,0),(440928,9,2,6,0),
       (-279936,9,1,7,0),
       (858240,8,5,3,1),(-695232,8,4,4,1),(-3003840,8,3,5,1),(2431104,8,2,6,1),
       (-293568,8,1,7,1),(-65664,8,0,8,1),
       (-23664,7,6,4,0),(3375456,7,5,5,0),(-3596160,7,4,6,0),(-787408,7,3,7,0),
       (1497904,7,2,8,0),(-319248,7,1,9,0),(2880,7,0,10,0),
       (-689616,6,6,4,1),(983088,6,5,5,1),(4351680,6,4,6,1),(-5500416,6,3,7,1),
       (1607824,6,2,8,1),(124208,6,1,9,1),(-66528,6,0,10,1),
       (-28404,5,7,5,0),(-2654640,5,6,6,0),(4049484,5,5,7,0),(862728,5,4,8,0),
       (-3145332,5,3,9,0),(1337976,5,2,10,0),(-173508,5,1,11,0),(7056,5,0,12,0),
       (228375,4,7,5,1),(-528921,4,6,6,1),(-3472677,4,5,7,1),(5849127,4,4,8,1),
       (-2567619,4,3,9,1),(-79515,4,2,10,1),(228705,4,1,11,1),(-24195,4,0,12,1),
       (913500,3,7,7,0),(-1831644,3,6,8,0),(-1118088,3,5,9,0),(3376968,3,4,10,0),
       (-1992756,3,3,11,0),(438660,3,2,12,0),(-41760,3,1,13,0),(5520,3,0,14,0),
       (1370250,2,6,8,1),(-2613276,2,5,9,1),(1273878,2,4,10,1),(300528,2,3,11,1),
       (-358290,2,2,12,1),(57420,2,1,13,1),(4050,2,0,14,1),
       (913500,1,6,10,0),(-1705284,1,5,11,0),(1146816,1,4,12,0),(-312912,1,3,13,0),
       (25380,1,2,14,0),(-2700,1,1,15,0),(1440,1,0,16,0),
       (228375,0,5,11,1),(-423135,0,4,12,1),(261126,0,3,13,1),(-43326,0,2,14,1),
       (-12285,0,1,15,1),(3645,0,0,16,1)))
const W7E_DEGREE = Dict(2 => 9, 4 => 17)

"Bracket of tex :1691-1698 / :1702-1730 from the hand-EXPANDED monomial table."
function w7e_bracket(order, a, r0, M, s)
    tot = zero(a * M * r0 * s)
    for (c, i, j, k, e) in W7E_EXPANDED[order]
        tot += c * a^i * M^j * r0^k * (e == 1 ? s : one(s))
    end
    return tot
end

"F_{r[order]} assembled from the EXPANDED bracket and the printed prefactor."
function w7e_reg_param(order, a, r0, M = one(float(promote_type(typeof(a), typeof(r0)))))
    s = sqrt(M * r0)
    b = w7e_bracket(order, a, r0, M, s)
    u = 2 * a * s + r0 * (r0 - 3 * M)          # [2a sqrt(Mr) + r(r-3M)]
    d = a^2 + r0 * (r0 - 2 * M)                # [a^2 + r(r-2M)]        (= Delta)
    if order == 2
        return M * b / (24 * r0^4 * (a * M + r0 * s) * (u * sqrt(u)) * (d * sqrt(d)))
    else
        return M^2 * b / (1440 * r0^9 * (a * M + r0 * s) * (u^3 * sqrt(u)) * (d * sqrt(d)))
    end
end

# (b) a = 0 Schwarzschild closed forms — a THIRD structural form, obtained by cancelling the
#     a = 0 limit of the [H14] Kerr expressions analytically (only the last bracket group
#     survives; the r^5 sqrt(Mr) / r^11 sqrt(Mr) factors cancel against the prefactor).
w7e_schw2(r0, M = one(r0)) =
    -M * (r0^3 + 3*M*r0^2 - 25*M^2*r0 + 29*M^3) /
     (8 * r0^3 * (r0 - 3*M)^(3//2) * (r0 - 2*M)^(3//2))
w7e_schw4(r0, M = one(r0)) =
    M^2 * (405*r0^5 - 1365*M*r0^4 - 4814*M^2*r0^3 + 29014*M^3*r0^2 - 47015*M^4*r0 + 25375*M^5) /
     (160 * r0^4 * (r0 - 3*M)^(7//2) * (r0 - 2*M)^(3//2))

# (c) the m-dependence gather :1675-1680 re-derived EXACTLY from its DEFINING Fourier
#     integral. [H14]:1640-1660 leaves the parameters as the m-modes of the singular-field
#     remainder [2*Theta(dphi) - 1] * sin^(2j-1)(dphi/2), whose folded (m >= 0, x2) mode is
#         F^m_{a[2j]} = 2 F_{a[2j]} * (1/2pi) int_-pi^pi |sin(x/2)|^(2j-1) e^{-imx} dx
#                     = F_{a[2j]} * (2/pi) int_0^pi sin^(2j-1)(x/2) cos(m x) dx .
#     Expanding sin^(2n+1)(u) = 4^-n sum_k (-1)^(n-k) C(2n+1,k) sin((2n+1-2k)u) and using
#     int_0^pi sin(q x/2) cos(m x) dx = 2q/(q^2 - 4m^2) for ODD q gives an EXACT rational
#     multiple of 1/pi that shares no algebraic structure with the printed product form.
function w7e_mfactor_over_pi(order, m)
    n = order ÷ 2 - 1
    tot = zero(Rational{BigInt})
    for k in 0:n
        q = 2n + 1 - 2k
        tot += Rational{BigInt}((-1)^(n - k) * binomial(big(2n + 1), big(k))) *
               Rational{BigInt}(2q, q^2 - 4 * big(m)^2)
    end
    return 2 * tot // big(4)^n
end

w7e_relerr(x, y) = abs(x - y) / max(abs(y), eps(Float64) * one(abs(y)))

# measured-deviation roll-up, printed at the end (package convention: green nested testsets
# are folded in the Pkg.test summary, so the numbers behind the gates get their own line)
const W7E_WORST = Dict{String,Float64}()
w7e_note!(k, v) = (W7E_WORST[k] = max(get(W7E_WORST, k, 0.0), Float64(v)); v)

@testset "heffernan-regularization (m-mode F_r regularization parameters, W7-E)" begin

# ---------------------------------------------------------------------------------------
@testset "(1) exported surface + printed m-dependence [H14]:1675-1680" begin
    @test HEFFERNAN_MMODE_ORDERS == (2, 4, 6, 8)
    @test HEFFERNAN_CIRCULAR_FR_ORDERS == (2, 4)

    # printed pattern, ASSERTED not assumed: F^m_{a[2j]} = c_{2j} / (pi * prod_k (2m)^2-(2k-1)^2)
    for order in HEFFERNAN_MMODE_ORDERS
        j = order ÷ 2
        c = (-1)^j * 4 * factorial(big(order - 1))     # -4, 24, -480, 20160
        for m in 0:25
            den = prod(Rational{BigInt}[(2m - (2k - 1)) * (2m + (2k - 1)) for k in 1:j])
            @test mmode_reg_factor_over_pi(order, m) == c // numerator(den)
        end
    end
    # the [2] and [4] patterns spelled out literally (the two orders this node ships)
    for m in 0:12
        @test mmode_reg_factor_over_pi(2, m) == -4 // ((2m - 1) * (2m + 1))
        @test mmode_reg_factor_over_pi(4, m) ==
              24 // ((2m - 3) * (2m - 1) * (2m + 1) * (2m + 3))
    end
    # EXACT re-derivation from the defining Fourier integral (structurally different form)
    for order in HEFFERNAN_MMODE_ORDERS, m in 0:24
        @test mmode_reg_factor_over_pi(order, m) == w7e_mfactor_over_pi(order, m)
    end
    # float wrapper == exact rational / pi, in every eltype
    for T in (Float64, BigFloat), order in HEFFERNAN_MMODE_ORDERS, m in (0, 1, 3, 7, 20)
        ref = T(mmode_reg_factor_over_pi(order, m)) / T(pi)
        @test mmode_reg_factor(T, order, m) isa T
        @test w7e_relerr(mmode_reg_factor(T, order, m), ref) <= (T === Float64 ? 1e-15 : 1e-70)
    end
    @test mmode_reg_factor(4, 7) === mmode_reg_factor(Float64, 4, 7)
    # large-m asymptotics: m^4 F^m_{[4]} -> 24/(16 pi) = 3/(2 pi); m^2 F^m_{[2]} -> -1/pi
    @test isapprox(mmode_reg_factor(4, 4000) * 4000.0^4, 3 / (2 * pi); rtol = 1e-6)
    @test isapprox(mmode_reg_factor(2, 4000) * 4000.0^2, -1 / pi; rtol = 1e-6)
    # unsupported orders are refused, not silently mis-evaluated
    @test_throws ArgumentError mmode_reg_factor_over_pi(3, 4)
    @test_throws ArgumentError mmode_reg_factor(Float64, 10, 4)
    @test_throws ArgumentError mmode_reg_factor_tail_over_pi(5, 4)
    @test_throws ArgumentError mmode_reg_factor(Float64, 4, -1)
end

# ---------------------------------------------------------------------------------------
@testset "(2) transcription re-diff: factored (src) vs expanded (test)" begin
    # (2a) EXACT Rational{BigInt} equality of the two bracket transcriptions. Points are
    #      chosen with M = p^2 u and r0 = q^2 u so that s = sqrt(M r0) = p q u is rational.
    npts = 0
    for (p, q, u) in ((1//1, 3//1, 1//1), (1//1, 5//1, 1//1), (2//1, 7//1, 3//2),
                      (1//1, 10//1, 1//1), (3//2, 11//2, 5//3), (1//1, 2//1, 1//1))
        M = Rational{BigInt}(p)^2 * Rational{BigInt}(u)
        r0 = Rational{BigInt}(q)^2 * Rational{BigInt}(u)
        s = Rational{BigInt}(p) * Rational{BigInt}(q) * Rational{BigInt}(u)
        @test s * s == M * r0
        for af in (0//1, 1//2, -1//2, 9//10, -9//10, 499//500, 3//7, -13//17)
            a = Rational{BigInt}(af) * M
            for order in HEFFERNAN_CIRCULAR_FR_ORDERS
                @test ScalarSelfForce.reg_param_Fr_bracket(order, a, r0, M, s) ==
                      w7e_bracket(order, a, r0, M, s)
                npts += 1
            end
        end
    end
    @test npts == 96
    # dimensional homogeneity of the expanded table (guards mistyped exponents)
    for order in HEFFERNAN_CIRCULAR_FR_ORDERS
        @test all(i + j + k + e == W7E_DEGREE[order] for (_, i, j, k, e) in W7E_EXPANDED[order])
        @test length(W7E_EXPANDED[order]) == (order == 2 ? 26 : 84)
    end
    # (2b) assembled parameter, Float64 and BigFloat(256), including M != 1
    setprecision(BigFloat, 256) do
        for T in (Float64, BigFloat), a in (0.0, 0.5, -0.5, 0.9, -0.9, 0.998, -0.998, 0.6),
            r0 in (6, 7, 10, 20, 50), order in HEFFERNAN_CIRCULAR_FR_ORDERS
            M = one(T)
            got = reg_param_Fr(order, T(a), T(r0); M = M)
            ref = w7e_reg_param(order, T(a), T(r0), M)
            @test got isa T
            @test w7e_note!("rediff_$T", w7e_relerr(got, ref)) <= (T === Float64 ? 1e-14 : 1e-60)
        end
        # M != 1 exercises the M-exponents of both transcriptions
        for T in (Float64, BigFloat), order in HEFFERNAN_CIRCULAR_FR_ORDERS
            M = T(7) / 5
            for a in (T(0), T(9) / 20, -T(9) / 20), x in (6, 10, 20)
                r0 = T(x) * M
                got = reg_param_Fr(order, a * M, r0; M = M)
                ref = w7e_reg_param(order, a * M, r0, M)
                @test w7e_relerr(got, ref) <= (T === Float64 ? 1e-14 : 1e-60)
                # dimensional scaling F_r[order] ~ 1/M^2
                @test w7e_relerr(got * M^2, reg_param_Fr(order, a, T(x))) <=
                      (T === Float64 ? 1e-13 : 1e-60)
            end
        end
    end
end

# ---------------------------------------------------------------------------------------
@testset "(3) cross-SOURCE oracle: [TW] notebook Fmr[4]" begin
    @test length(HEFF_MMODE_FMR4) == 200
    setprecision(BigFloat, 256) do
        worst64 = 0.0
        worstbig = big(0.0)
        for ((ar, r0, m), dec) in HEFF_MMODE_FMR4
            orb64 = CircularOrbit(Float64(ar), Float64(r0))
            got64 = mmode_reg_Fr(4, orb64, m)
            ref64 = parse(Float64, dec)
            worst64 = max(worst64, Float64(w7e_relerr(got64, ref64)))
            orbb = CircularOrbit(BigFloat(numerator(ar)) / denominator(ar), BigFloat(r0))
            gotb = mmode_reg_Fr(4, orbb, m)
            refb = parse(BigFloat, dec)
            worstbig = max(worstbig, w7e_relerr(gotb, refb))
        end
        @test w7e_note!("oracle_Float64", worst64) <= 1e-14
        @test w7e_note!("oracle_BigFloat", worstbig) <= 1e-40
        # the oracle carries the FULL m-factor: coefficient x factor == oracle
        for ar in HEFF_MMODE_ORACLE_A, r0 in HEFF_MMODE_ORACLE_R0, m in HEFF_MMODE_ORACLE_MS
            orb = CircularOrbit(Float64(ar), Float64(r0))
            @test w7e_relerr(reg_param_Fr(4, orb) * mmode_reg_factor(4, m),
                             parse(Float64, HEFF_MMODE_FMR4[(ar, r0, m)])) <= 1e-14
        end
    end
    # fixture header pins
    @test HEFF_MMODE_NB_CELL_LINES == 2806:9600
    @test length(HEFF_MMODE_NB_SHA256) == 64
    @test HEFF_MMODE_ORACLE_DIGITS == 50
    @test HEFF_MMODE_NB_PATH ==
          "ref-paper/arxiv-1610.09319/src/anc/m-ModeKerrRegularizationParameters.nb"
    # DRIFT PIN against the READ-ONLY mirror itself (in-package pattern:
    # test_runtime_config.jl:43-48, test_wb_fixtures.jl:158-164). Any change to the notebook
    # invalidates the baked oracle and must force a fixture review; an absent mirror records a
    # visible Broken (GRTOOLS_STRICT_ORACLES=1 promotes it to a failure).
    let nb = joinpath(oracle_root("GRTOOLS_REFPAPER_ROOT",
                                  "/data/haiyangw/claude/emri-disk/ref-paper"),
                      "arxiv-1610.09319", "src", "anc",
                      "m-ModeKerrRegularizationParameters.nb")
        if oracle_or_skip(isfile(nb),
                          "Wardell-Thornburg m-mode reg-parameter notebook sha256 drift pin";
                          env = "GRTOOLS_REFPAPER_ROOT")
            @test bytes2hex(sha256(read(nb))) == HEFF_MMODE_NB_SHA256
            # the recorded window still bounds exactly one `Cell[BoxData[…]], "Input",` block —
            # the SAME two structural asserts the capture script makes on the STRIPPED window
            # (its first line, 2806, is a blank separator; the cell text starts on 2807)
            nblines = readlines(nb)
            @test length(nblines) >= last(HEFF_MMODE_NB_CELL_LINES)
            cellblob = strip(join(nblines[HEFF_MMODE_NB_CELL_LINES], "\n"))
            @test startswith(cellblob, "Cell[BoxData[")
            @test endswith(cellblob, "\";\"}]], \"Input\",")
            @test count(==("Cell[BoxData["), strip.(nblines[HEFF_MMODE_NB_CELL_LINES])) == 1
        end
    end
end

# ---------------------------------------------------------------------------------------
@testset "(4) a -> 0 Schwarzschild limit + retrograde continuity" begin
    setprecision(BigFloat, 256) do
        for T in (Float64, BigFloat), r0 in (6, 7, 10, 20)
            @test w7e_note!("schwarzschild_$T",
                            w7e_relerr(reg_param_Fr(2, zero(T), T(r0)), w7e_schw2(T(r0)))) <=
                  (T === Float64 ? 1e-12 : 1e-60)
            @test w7e_note!("schwarzschild_$T",
                            w7e_relerr(reg_param_Fr(4, zero(T), T(r0)), w7e_schw4(T(r0)))) <=
                  (T === Float64 ? 1e-12 : 1e-60)
        end
    end
    # the a = 0 rows of the [TW] cross-source oracle ARE the independent Schwarzschild
    # statement (see this file's header, gate (2)): they must equal the reduced forms.
    for r0 in HEFF_MMODE_ORACLE_R0, m in HEFF_MMODE_ORACLE_MS
        ref = parse(Float64, HEFF_MMODE_FMR4[(0//1, r0, m)])
        @test w7e_relerr(w7e_schw4(Float64(r0)) * mmode_reg_factor(4, m), ref) <= 1e-12
    end
    # retrograde a < 0 is a first-class input and the parameters are CONTINUOUS through 0
    for order in HEFFERNAN_CIRCULAR_FR_ORDERS, r0 in (6.0, 10.0, 20.0)
        f0 = reg_param_Fr(order, 0.0, r0)
        for h in (1e-6, 1e-8, 1e-10)
            fp, fm = reg_param_Fr(order, h, r0), reg_param_Fr(order, -h, r0)
            @test w7e_relerr(fp, f0) <= 60 * h        # C^1 in a: |df| = O(h)
            @test w7e_relerr(fm, f0) <= 60 * h
            # the a-derivative is one-sided-consistent: (f(+h) - f(-h))/2h stays bounded
            @test abs((fp - fm) / (2h)) <= 100 * abs(f0)
        end
        # prograde/retrograde asymmetry is REAL (not a mirrored value)
        @test reg_param_Fr(order, 0.9, r0) != reg_param_Fr(order, -0.9, r0)
    end
    # retrograde values reproduce the [TW] oracle just as prograde ones do (fixture grid
    # deliberately contains a < 0, adopted decision 6)
    for r0 in HEFF_MMODE_ORACLE_R0, m in (0, 2, 20), ar in (-1//2, -9//10, -499//500)
        orb = CircularOrbit(Float64(ar), Float64(r0))
        @test w7e_relerr(mmode_reg_Fr(4, orb, m),
                         parse(Float64, HEFF_MMODE_FMR4[(ar, r0, m)])) <= 1e-14
    end
end

# ---------------------------------------------------------------------------------------
@testset "(5) m = 0 folding + the analytic tail completion" begin
    # EXACT closed-form tail: sum_{m>=mstart} F^m_{[2j]} telescopes.
    for order in HEFFERNAN_MMODE_ORDERS, mstart in (0, 1, 2, 5, 11, 21)
        for N in (mstart + 1, mstart + 7, mstart + 500, mstart + 3000)
            partial = sum(mmode_reg_factor_over_pi(order, m) for m in mstart:(N - 1))
            @test partial + mmode_reg_factor_tail_over_pi(order, N) ==
                  mmode_reg_factor_tail_over_pi(order, mstart)
        end
    end
    # printed-form spot checks of the tail (derived once, pinned here)
    for m in (1, 2, 5, 21)
        @test mmode_reg_factor_tail_over_pi(2, m) == -2 // (2m - 1)
        @test mmode_reg_factor_tail_over_pi(4, m) == 4 // ((2m - 3) * (2m - 1) * (2m + 1))
    end
    # m = 0 SPECIAL CASE (operational rule, identical wording to the mmode_reg_factor_over_pi
    # docstring): the closed form is a single analytic function of m with no m = 0 branch, so
    # the x2 real-part folding is INSIDE the parameter at every m, m = 0 included. Node 46
    # (mode-sum-total-force) must use mmode_reg_Fr(order, orbit, 0)/2 wherever its m = 0 force
    # term is the UNFOLDED q grad(Psi^R_0 / r) of Osburn-Nishimura Eq. (3.10) (first branch);
    # for m >= 1 use the value as returned. Equivalently the doubly-infinite sum vanishes:
    # tail(order, 1) == -F^0_{[2j]}/2, asserted exactly below.
    for order in HEFFERNAN_MMODE_ORDERS
        @test mmode_reg_factor_tail_over_pi(order, 1) == -mmode_reg_factor_over_pi(order, 0) // 2
        @test mmode_reg_factor_tail_over_pi(order, 0) == mmode_reg_factor_over_pi(order, 0) // 2
    end
    # float wrappers and the orbit-level tail
    setprecision(BigFloat, 256) do
        for T in (Float64, BigFloat), order in HEFFERNAN_MMODE_ORDERS, mstart in (0, 1, 5, 21)
            ref = T(mmode_reg_factor_tail_over_pi(order, mstart)) / T(pi)
            @test mmode_reg_factor_tail(T, order, mstart) isa T
            @test w7e_relerr(mmode_reg_factor_tail(T, order, mstart), ref) <=
                  (T === Float64 ? 1e-15 : 1e-70)
        end
        for order in HEFFERNAN_CIRCULAR_FR_ORDERS, a in (0.0, 0.5, -0.9), r0 in (6.0, 10.0)
            orb = CircularOrbit(a, r0)
            for mstart in (1, 5, 21)
                @test w7e_relerr(mmode_reg_tail_Fr(order, orb, mstart),
                                 reg_param_Fr(order, orb) *
                                 mmode_reg_factor_tail(Float64, order, mstart)) <= 1e-14
                # the raw 4-positional method (order, a, r0, mstart; M) agrees with the orbit
                # method and, like the coefficient itself, scales as 1/M^2 (mirrors the M != 1
                # check on reg_param_Fr above; the m-factor is dimensionless)
                @test mmode_reg_tail_Fr(order, a, r0, mstart) ==
                      mmode_reg_tail_Fr(order, orb, mstart)
                M = 7 / 5
                got = mmode_reg_tail_Fr(order, a * M, r0 * M, mstart; M = M)
                @test got isa Float64
                @test w7e_relerr(got * M^2, mmode_reg_tail_Fr(order, a, r0, mstart)) <= 1e-13
            end
        end
    end
    # brute-force truncated partial sum converges onto the closed-form tail. Tolerances are
    # the MEASURED double-precision truncation of the 1/N ([2]) and 1/N^3 ([4]) remainders
    # at N = 200_000 (calibration: 2.25e-5 and 7.03e-13; the [4] value is already at its
    # roundoff floor -- unchanged at N = 2e6).
    let mstart = 5, N = 200_000
        for order in HEFFERNAN_CIRCULAR_FR_ORDERS
            brute = sum(mmode_reg_factor(Float64, order, m) for m in mstart:N)
            @test w7e_relerr(brute, mmode_reg_factor_tail(Float64, order, mstart)) <=
                  (order == 2 ? 5e-5 : 5e-12)
        end
    end
    @test mmode_reg_factor_tail(4, 7) === mmode_reg_factor_tail(Float64, 4, 7)
end

# ---------------------------------------------------------------------------------------
@testset "(6) manufactured/tautological subtraction-order placeholder [FUTURE: real F_r^m]" begin
    # Row 45's slope gate wants the -6 +/- 0.3 log-log slope of |F_r^m - F^m_{r[4]}| on REAL
    # node-44 data, which does not exist yet [FUTURE]. This is a PLACEHOLDER and is
    # TAUTOLOGICAL BY CONSTRUCTION: the sequence is built as F^m_{[4]} + C * F^m_{[6]}, so the
    # residual after subtracting the [4] term is EXACTLY C * mmode_reg_factor(6, m) and the
    # -6 slope only re-measures that the shipped order-6 m-factor behaves as m^-6. It gates
    # the plumbing (the subtraction removes all m^-4 content, the slope estimator node 46 will
    # reuse works) and NOTHING about real residual data. It must be REPLACED, not merely
    # supplemented, once mode-self-force-extraction (row 44) lands.
    orb = CircularOrbit(0.6, 10.0)                      # the [H14]:1784-1790 example orbit
    F4 = reg_param_Fr(4, orb)
    # sub-dominant order-6 admixture: kappa = 0.4 puts the [6] content at 2.2% of the [4]
    # content at the bottom of the window (measured), so the [4] term genuinely dominates.
    # ([H14]:1795 quotes a FITTED F_{r[6]} = 0.108797 for this orbit, which is ~160x F_{r[4]}
    # and would DOMINATE over m in 10..60 -- a synthetic sequence built with it does not test
    # what this gate is about, so the amplitude here is chosen, not quoted.)
    C = 0.4 * F4
    ms = 10:60
    seq = [F4 * mmode_reg_factor(4, m) + C * mmode_reg_factor(6, m) for m in ms]
    raw = [abs(s) for s in seq]
    res = [abs(seq[i] - F4 * mmode_reg_factor(4, ms[i])) for i in eachindex(ms)]
    slope(y) = let x = log.(Float64.(ms)), ly = log.(y), n = length(x),
                   sx = sum(x), sy = sum(ly)
        (n * sum(x .* ly) - sx * sy) / (n * sum(x .^ 2) - sx^2)
    end
    W7E_WORST["slope_raw"] = slope(raw)
    W7E_WORST["slope_after_subtraction"] = slope(res)
    @test isapprox(slope(raw), -4.0; atol = 0.3)        # before subtraction: m^-4 (-4.0016)
    @test isapprox(slope(res), -6.0; atol = 0.3)        # after subtraction: m^-6 (-6.0345)
    @test maximum(res ./ raw) < 0.03                    # the [4] part really dominates (0.022)
end

# ---------------------------------------------------------------------------------------
@testset "(7) type genericity, orbit adapter, domain" begin
    setprecision(BigFloat, 256) do
        orb = CircularOrbit(BigFloat(1) / 2, BigFloat(10))
        for order in HEFFERNAN_CIRCULAR_FR_ORDERS
            @test reg_param_Fr(order, orb) isa BigFloat
            @test mmode_reg_Fr(order, orb, 3) isa BigFloat
            @test mmode_reg_tail_Fr(order, orb, 3) isa BigFloat
            @test precision(reg_param_Fr(order, orb)) == 256
            # the orbit method agrees with the raw (a, r0) method
            @test reg_param_Fr(order, orb) ==
                  reg_param_Fr(order, BigFloat(1) / 2, BigFloat(10); M = BigFloat(1))
        end
    end
    setprecision(BigFloat, 512) do
        orb = CircularOrbit(BigFloat(1) / 2, BigFloat(10))
        @test precision(reg_param_Fr(4, orb)) == 512
        lo = setprecision(BigFloat, 128) do
            reg_param_Fr(4, CircularOrbit(BigFloat(1) / 2, BigFloat(10)))
        end
        @test w7e_relerr(BigFloat(lo), reg_param_Fr(4, orb)) <= 1e-35
    end
    # mixed-eltype promotion, never a silent Float64 down-cast
    @test reg_param_Fr(4, 0, 10) isa Float64
    @test reg_param_Fr(4, BigFloat(0), 10) isa BigFloat
    @test mmode_reg_Fr(4, 0.5f0, 10.0f0, 3) isa Float32
    # orders this node does NOT ship a circular closed form for are refused loudly
    @test_throws ArgumentError reg_param_Fr(6, 0.5, 10.0)
    @test_throws ArgumentError reg_param_Fr(8, CircularOrbit(0.5, 10.0))
    @test_throws ArgumentError mmode_reg_Fr(6, CircularOrbit(0.5, 10.0), 3)
    # inside the (retrograde/prograde) photon orbit there is no circular geodesic: the
    # printed bracket [2a sqrt(Mr) + r(r-3M)] goes negative and the parameter is undefined
    @test_throws DomainError reg_param_Fr(4, 0.9, 1.5)   # a=0.9 prograde photon orbit ~1.561M
    @test_throws DomainError reg_param_Fr(2, -0.9, 3.5)
    @test_throws DomainError reg_param_Fr(4, 0.0, 2.5)
end

end

# one-line visibility markers (Pkg.test folds green nested testsets)
println("heffernan-regularization gate (1) worst factored-vs-expanded transcription rel dev: ",
        "Float64 ", W7E_WORST["rediff_Float64"], ", BigFloat256 ", W7E_WORST["rediff_BigFloat"])
println("heffernan-regularization gate (2) worst vs Wardell-Thornburg notebook oracle rel dev: ",
        "Float64 ", W7E_WORST["oracle_Float64"], ", BigFloat256 ", W7E_WORST["oracle_BigFloat"])
println("heffernan-regularization gate (3) worst a=0 Schwarzschild-limit rel dev: ",
        "Float64 ", W7E_WORST["schwarzschild_Float64"],
        ", BigFloat256 ", W7E_WORST["schwarzschild_BigFloat"])
println("heffernan-regularization gate (4) manufactured/TAUTOLOGICAL placeholder ",
        "(residual == C*factor(6,m) by construction; the row-45 slope gate needs node-44 ",
        "data [FUTURE]): log-log slope raw ", W7E_WORST["slope_raw"],
        " -> after F^m_{r[4]} subtraction ", W7E_WORST["slope_after_subtraction"])

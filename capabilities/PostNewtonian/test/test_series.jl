# Characterization tests for PostNewtonian.Series (the PNSeries Num-algebra engine).
#
# GROUNDING = the captured MMA-computed oracle (results/gw-flux/reference/gen/gen_pn_series_ops.wl),
# vendored verbatim as test/fixtures/pn_series_ops_oracle.json. Each op is compared to the oracle
# STRUCTURALLY (nmin/nmax/den slot bookkeeping) and NUMERICALLY (evaluate Normal[result] at the MMA
# probe points (v,η) and match the MMA InputForm string to atol 1e-9 — the legacy tolerance). The six
# probe series mkA..mkG and probe points (PV,PET) are the legacy baked literals; NO live TeukolskyPN
# dep. ExpandLog/ExpandGamma are grounded by closed-form log/Γ identities. polyToSeries / rescale /
# power algebra additionally cross-checked against analytic Maclaurin coefficients.

using Test
using PostNewtonian.Series
import Symbolics
import JSON

const ORACLE_PATH = joinpath(@__DIR__, "fixtures", "pn_series_ops_oracle.json")
const ORACLE = JSON.parsefile(ORACLE_PATH)

Symbolics.@variables v η x ω

# The six oracle input series (legacy mkA..mkG baked literals).
mkA() = PNSeries(:η, [Symbolics.Num(1), 1 + v, 1 + v + v^2 / 2, 1 + v + v^2 / 2,
                      1 + v + v^2 / 2, 1 + v + v^2 / 2], 0, 6, 1)
mkB() = PNSeries(:η, [v, 1 + v^2, v + v^3, v^2 + v^4, v^3 + v^5, v^4 + v^6], 0, 6, 1)
mkC() = PNSeries(:η, [Symbolics.Num(1), v, ((-1 + v) * v) / 2, (v * (2 - 3v + v^2)) / 6,
                      (v * (-6 + 11v - 6v^2 + v^3)) / 24], 0, 5, 1)
mkD() = PNSeries(:η, [Symbolics.Num(1), v, v^2 / 2, v^3 / 6, v^4 / 24, v^5 / 120], -1, 5, 1)
mkF() = PNSeries(:η, [Symbolics.Num(1), v / 2, -v^2 / 8, v^3 / 16, -5v^4 / 128], 0, 5, 1)
mkG() = PNSeries(:η, [Symbolics.Num(1), Symbolics.Num(0), v], 1, 7, 2)

# MMA probe points used in gen_pn_series_ops.wl `numProbe`.
const PV  = (0.3, 1.1, -0.7)
const PET = (0.01, 0.05, -0.03)

# Evaluate Normal[s] at (vval, etaval) with COMPLEX arithmetic (half-integer powers of negative η).
function eval_normal(s::PNSeries, vval, etaval)
    acc = ComplexF64(0)
    base = ComplexF64(etaval) - ComplexF64(Symbolics.value(s.x0))
    for (i, c) in enumerate(s.coeffs)
        isequal(c, Symbolics.Num(0)) && continue
        e = (s.nmin + i - 1) / s.den
        sub = Symbolics.value(Symbolics.substitute(c, Dict(v => vval)))
        acc += ComplexF64(sub) * base^e
    end
    acc
end

# Parse an MMA numeric probe string (InputForm of N[...,30]); supports "...*I".
function parse_mma_num(str::AbstractString)
    s = replace(String(str), r"`[0-9.]*" => "")
    s = strip(s)
    if occursin("*I", s) || s == "I" || endswith(s, "I")
        s2 = replace(s, "*I" => "im")
        s2 = replace(s2, r"(?<![0-9.)])I" => "1im")
        return ComplexF64(eval(Meta.parse(s2)))
    end
    return ComplexF64(parse(Float64, s))
end

# Numeric agreement of a PNSeries result with a vector of MMA probe strings.
function probes_match(s::PNSeries, mmastrs; atol = 1e-9)
    @assert length(mmastrs) == length(PV)
    ok = true
    maxerr = 0.0
    for k in 1:length(PV)
        got = eval_normal(s, PV[k], PET[k])
        want = parse_mma_num(mmastrs[k])
        err = abs(got - want)
        maxerr = max(maxerr, err)
        ok &= (err <= atol * (1 + abs(want)))
    end
    (ok, maxerr)
end

slots(o) = (o["nmin"], o["nmax"], o["den"])
zero_sym(d) = isequal(Symbolics.simplify(Symbolics.expand(Symbolics.Num(d))), Symbolics.Num(0))

@testset "Series: PNSeries algebra vs MMA oracle" begin

    @testset "input series structure" begin
        for (name, mk) in (("A", mkA), ("B", mkB), ("C", mkC), ("D", mkD), ("F", mkF), ("G", mkG))
            o = ORACLE["inputs"][name]
            s = mk()
            @test (s.nmin, s.nmax, s.den) == slots(o)
            @test length(s.coeffs) == o["nmax"] - o["nmin"]
        end
    end

    @testset "SeriesMinOrder / MaxOrder / Length" begin
        for (name, mk) in (("A", mkA), ("D", mkD), ("G", mkG))
            o = ORACLE["orders"][name]
            s = mk()
            @test series_min_order(s) == eval(Meta.parse(replace(o["min"], "/" => "//")))
            @test series_max_order(s) == eval(Meta.parse(replace(o["max"], "/" => "//")))
            @test series_length(s)    == eval(Meta.parse(replace(o["len"], "/" => "//")))
        end
        @test series_min_order(1) == 0
    end

    @testset "arithmetic truncation vs MMA" begin
        A, B, C, D, G = mkA(), mkB(), mkC(), mkD(), mkG()
        cases = Dict(
            "A_plus_B"      => A + B,
            "A_times_B"     => A * B,
            "inv_A"         => inv(A),
            "A_pow_3"       => A^3,
            "A_pow_m2"      => A^(-2),
            "A_div_B"       => A / B,
            "A_minus_B"     => A - B,
            "scalar_plus_A" => 3 + A,
            "num_times_A"   => v^2 * A,
            "A_plus_num"    => A + v,
            "D_times_A"     => D * A,
            "D_pow_2"       => D^2,
            "C_times_C"     => C * C,
            "G_times_G"     => G * G,
            "inv_G"         => inv(G),
        )
        maxerr_all = 0.0
        for (k, res) in cases
            o = ORACLE["arith"][k]
            @test (res.nmin, res.nmax, res.den) == slots(o)
            ok, err = probes_match(res, ORACLE["arith_numeric"][k])
            maxerr_all = max(maxerr_all, err)
            @test ok
        end
        @info "Series arithmetic max numeric |got-MMA|" maxerr = maxerr_all
    end

    @testset "SeriesTake" begin
        A, B, D, G = mkA(), mkB(), mkD(), mkG()
        @test slots(ORACLE["take"]["A_take_1"]) == (series_take(A, 1).nmin, series_take(A, 1).nmax, series_take(A, 1).den)
        @test slots(ORACLE["take"]["A_take_3"]) == (series_take(A, 3).nmin, series_take(A, 3).nmax, series_take(A, 3).den)
        @test slots(ORACLE["take"]["D_take_2"]) == (series_take(D, 2).nmin, series_take(D, 2).nmax, series_take(D, 2).den)
        @test slots(ORACLE["take"]["B_take_4"]) == (series_take(B, 4).nmin, series_take(B, 4).nmax, series_take(B, 4).den)
        @test slots(ORACLE["take"]["G_take_2"]) == (series_take(G, 2).nmin, series_take(G, 2).nmax, series_take(G, 2).den)
        t3 = series_take(A, 3)
        for (i, c) in enumerate(A.coeffs[1:3])
            @test zero_sym(t3.coeffs[i] - c)
        end
        t0 = series_take(A, 0)
        @test (t0.nmin, t0.nmax) == (0, 0)
        @test isempty(t0.coeffs)
    end

    @testset "SeriesCollect" begin
        A, B = mkA(), mkB()
        ce = series_collect(A * B, :v, Symbolics.expand)
        o  = ORACLE["collect"]["A_collect_v_expand"]
        @test (ce.nmin, ce.nmax, ce.den) == slots(o)
        @test zero_sym(ce.coeffs[1] - v)
        @test zero_sym(ce.coeffs[2] - (1 + v + 2v^2))
        ci = series_collect(A, :v)
        for k in 1:length(PV)
            @test abs(eval_normal(ci, PV[k], PET[k]) - eval_normal(A, PV[k], PET[k])) < 1e-12
        end
    end

    @testset "SeriesTerms" begin
        t1 = series_terms((1 + v * η) / (1 - η), :η, 0, 3)
        @test slots(ORACLE["terms"]["expr1_3terms"]) == (t1.nmin, t1.nmax, t1.den)
        @test zero_sym(t1.coeffs[2] - (1 + v)) && zero_sym(t1.coeffs[3] - (1 + v))
        t2 = series_terms(exp(v * η) / η, :η, 0, 2)
        @test slots(ORACLE["terms"]["expr2_2terms"]) == (t2.nmin, t2.nmax, t2.den)
        @test zero_sym(t2.coeffs[2] - v)
        t3 = series_terms(sqrt(1 + v * η), :η, 0, 4)
        @test slots(ORACLE["terms"]["expr3_4terms"]) == (t3.nmin, t3.nmax, t3.den)
        # closed-form Maclaurin coeffs of sqrt(1+vη): [1, v/2, -v^2/8, v^3/16]
        @test zero_sym(t3.coeffs[2] - v / 2)
        @test zero_sym(t3.coeffs[3] + v^2 / 8)
        @test zero_sym(t3.coeffs[4] - v^3 / 16)
    end

    @testset "SeriesTerms high Laurent orders + inverse actual-lead truncation" begin
        hi = series_terms(η^20 * (2 + 3η + 4η^2), :η, 0, 3)
        @test (hi.nmin, hi.nmax, hi.den) == (20, 23, 1)
        @test all(zero_sym(hi.coeffs[i] - (2, 3, 4)[i]) for i in eachindex(hi.coeffs))

        lo = series_terms(η^(-20) * (2 + 3η + 4η^2), :η, 0, 3)
        @test (lo.nmin, lo.nmax, lo.den) == (-20, -17, 1)
        @test all(zero_sym(lo.coeffs[i] - (2, 3, 4)[i]) for i in eachindex(lo.coeffs))

        # Preserve the same mathematical series but expose two nominal leading zero slots.  The
        # reciprocal is constrained only by the actual lead η^20 and must stop at O[η]^-17.
        padded = PNSeries(:η, vcat(Symbolics.Num[0, 0], hi.coeffs), 18, 23, 1)
        pinv = inv(padded)
        @test (pinv.nmin, pinv.nmax, pinv.den) == (-20, -17, 1)
        @test zero_sym(pinv.coeffs[1] - 1 // 2)
        @test zero_sym(pinv.coeffs[2] + 3 // 4)
        @test zero_sym(pinv.coeffs[3] - 1 // 8)

        @test_throws ArgumentError series_terms(η^100, :η, 0, 2)
        @test_throws ArgumentError series_terms(Symbolics.Num(0), :η, 0, 2)
    end

    @testset "polyToSeries" begin
        p1 = polyToSeries(1 + v * η + v^2 * η^3, :η)
        @test slots(ORACLE["polyToSeries"]["p1"]) == (p1.nmin, p1.nmax, p1.den)
        @test zero_sym(p1.coeffs[3])
        @test zero_sym(p1.coeffs[4] - v^2)
        p2 = polyToSeries(v * η^2, :η)
        @test slots(ORACLE["polyToSeries"]["p2"]) == (p2.nmin, p2.nmax, p2.den)
        @test polyToSeries(Symbolics.Num(0), :η) == Symbolics.Num(0)
    end

    @testset "ChangeSeriesParameter" begin
        A, G = mkA(), mkG()
        r1 = change_series_parameter(A, :x)
        @test slots(ORACLE["change"]["A_rename"]) == (r1.nmin, r1.nmax, r1.den)
        @test r1.var == :x
        r2 = change_series_parameter(A, Val(:rescale), 2, :ω)
        @test slots(ORACLE["change"]["A_rescale"]) == (r2.nmin, r2.nmax, r2.den)
        @test zero_sym(r2.coeffs[2] - 2 * (1 + v))
        @test zero_sym(r2.coeffs[3] - 4 * (1 + v + v^2 / 2))
        r3 = change_series_parameter(A, Val(:power), :x, 2)
        @test slots(ORACLE["change"]["A_pow2"]) == (r3.nmin, r3.nmax, r3.den)
        @test zero_sym(r3.coeffs[2])
        @test zero_sym(r3.coeffs[3] - (1 + v))
        r4 = change_series_parameter(A, Val(:facpower), 3, :x, 2)
        @test slots(ORACLE["change"]["A_facpow"]) == (r4.nmin, r4.nmax, r4.den)
        @test zero_sym(r4.coeffs[3] - 3 * (1 + v))
        @test zero_sym(r4.coeffs[5] - 9 * (1 + v + v^2 / 2))
        rg = change_series_parameter(G, Val(:power), :x, 2)
        @test slots(ORACLE["change"]["G_pow2"]) == (rg.nmin, rg.nmax, rg.den)
    end

    @testset "PowerCounting" begin
        A, D, G = mkA(), mkD(), mkG()
        pa = power_counting(A, :x)
        @test (pa.nmin, pa.nmax, pa.den) == slots(ORACLE["powercounting"]["A_pc"])
        @test pa.var == :x
        @test zero_sym(pa.coeffs[2] - (1 + v) * η)
        @test zero_sym(pa.coeffs[3] - (1 + v + v^2 / 2) * η^2)
        pd = power_counting(D, :x)
        @test (pd.nmin, pd.nmax, pd.den) == slots(ORACLE["powercounting"]["D_pc"])
        @test zero_sym(pd.coeffs[1] - η^(-1))
        pg = power_counting(G, :x)
        @test (pg.nmin, pg.nmax, pg.den) == slots(ORACLE["powercounting"]["G_pc"])
    end

    @testset "Normal pass-through" begin
        @test Normal(Symbolics.Num(5)) == Symbolics.Num(5)
    end

    @testset "ExpandLog (closed-form log identities)" begin
        Symbolics.@variables a b c n
        @test zero_sym(expand_log(log(a * b * c)) - (log(a) + log(b) + log(c)))
        @test zero_sym(expand_log(log(x^n)) - n * log(x))
        @test zero_sym(expand_log(log(a^2 * b / c)) - (2 * log(a) + log(b) - log(c)))
        @test zero_sym(expand_log(2 * log(v * η) + log((1 - η)^3)) -
                       (3 * log(1 - η) + 2 * (log(v) + log(η))))
    end

    @testset "ExpandGamma (closed-form Γ recurrence)" begin
        Symbolics.@variables a b c n y
        G3 = Series.SpecialFunctions.gamma
        @test zero_sym(expand_gamma(G3(3 + x)) - x * (1 + x) * (2 + x) * G3(x))
        @test zero_sym(expand_gamma(G3(x)) - G3(x))
        @test zero_sym(expand_gamma(G3(2 + x) + G3(1 + y)) -
                       (x * (1 + x) * G3(x) + y * G3(y)))
        @test zero_sym(expand_gamma(G3(n + 1)) - n * G3(n))
        @test zero_sym(expand_gamma(v * G3(4 + x)) -
                       v * x * (1 + x) * (2 + x) * (3 + x) * G3(x))
        @test zero_sym(expand_gamma(G3(3 + x); shift = 1) -
                       (1 + x) * (2 + x) * G3(1 + x))
    end
end

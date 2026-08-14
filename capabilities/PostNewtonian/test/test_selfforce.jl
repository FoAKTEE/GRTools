# Characterization tests for PostNewtonian.SelfForce (the PN self-force accessor layer).
#
# GROUNDING:
#  - PostNewtonianData accessors / order metadata: ANALYTIC, the closed-form upstream MakeBoxes rule
#    (startOrder = nmin//den, maxOrder = (nmax-1)//den) on a hand-built series.
#  - postnewtonian_coefficient: WOLFRAMSCRIPT 14.3 oracle on the bundled Infinity-l2m2.m (baked from
#    coefficients_test.jl): SeriesCoefficient[s,7]=153088/6615; Coefficient[SC[s,9],Logy]=2930944/11025;
#    SC[s,9]/.sub1=-240693798848/81265275; (SC[s,9]/.{Logy^9->0,Logy->0})/.sub1=-254196291488/81265275.
#  - ResummedSeriesData: ANALYTIC closed forms transcribed from Kernel/ResummedSeries.m (accessors,
#    collapse-when-not-a-series, Normal == sf*Normal(sd)).
#  - postnewtonian_expansion: structural round-trip (unique query loads the same series the loader
#    yields; empty query browses the full 91-key catalog).
#  - Schwarzschild circular-orbit energy leading orders {1,0,-1/2,0,3/8}: published PN series, via the
#    a->0 Energy.m coefficients.

using Test
using PostNewtonian.SelfForce
using PostNewtonian.Series: PNSeries, Normal
using PostNewtonian.Catalog: SERIES_DATA_DIR
using PostNewtonian.MMAInputForm: parse_pn_data_file
import Symbolics

@testset "SelfForce: PostNewtonianData accessors + order metadata" begin
    series = PNSeries(:y, [1, 0, -1//2, 0, 3//8, 0, 0, 0, 0], 3, 12, 2)   # var y, nmin 3, nmax 12, den 2
    pnd = PostNewtonianData(Dict{String,Any}(
        "Name"        => "Test Quantity",
        "Description" => "A hand-built series for the struct test.",
        "Authors"     => ["A. Author", "B. Author"],
        "References"  => ["arXiv:0000.00000"],
        "Series"      => series,
    ))
    @test pnd["Name"] == "Test Quantity"
    @test pnd["Authors"] == ["A. Author", "B. Author"]
    @test pnd["Series"] === series
    @test Set(keys(pnd)) == Set(["Name", "Description", "Authors", "References", "Series"])
    # order metadata: closed-form Rationals (startOrder nmin//den, maxOrder (nmax-1)//den).
    @test pn_start_order(pnd) == 3 // 2
    @test pn_start_order(pnd) isa Rational
    @test pn_max_order(pnd) == (12 - 1) // 2     # 11//2
    @test pn_max_order(pnd) isa Rational
    # integer-order case
    series_int = PNSeries(:y, collect(0:9) .// 1, 0, 11, 1)
    pnd_int = PostNewtonianData(Dict{String,Any}("Name" => "Int", "Series" => series_int))
    @test pn_start_order(pnd_int) == 0 // 1
    @test pn_max_order(pnd_int) == 10 // 1
    str = sprint(show, pnd)
    @test occursin("Test Quantity", str) && occursin("3//2", str) && occursin("11//2", str)
end

@testset "SelfForce: postnewtonian_coefficient on l2m2 vs wolframscript oracle" begin
    path = joinpath(SERIES_DATA_DIR, "Schwarzschild", "Circular", "Flux", "Infinity-l2m2.m")
    pnd = PostNewtonianData(parse_pn_data_file(path))

    c7 = postnewtonian_coefficient(pnd, 7)
    @test Rational{BigInt}(Symbolics.value(c7)) == 153088 // 6615

    c9L1 = postnewtonian_coefficient(pnd, 9, 1)
    @test Rational{BigInt}(Symbolics.value(c9L1)) == 2930944 // 11025

    @test Rational{BigInt}(Symbolics.value(postnewtonian_coefficient(pnd, 7, 1))) == 0
    @test Rational{BigInt}(Symbolics.value(postnewtonian_coefficient(pnd, 7, 0))) == 153088 // 6615
    @test Symbolics.value(postnewtonian_coefficient(pnd, 1000)) == 0   # off-grid -> 0

    sub1 = Dict(
        Symbolics.variable(:Pi) => Symbolics.Num(big(22) // 7),
        Symbolics.variable(:EulerGamma) => Symbolics.Num(big(3) // 5),
        Symbolics.variable(:Logy) => Symbolics.Num(big(5) // 8),
        Symbolics.variable(:Log2) => Symbolics.Num(big(4) // 9),
    )
    c9 = postnewtonian_coefficient(pnd, 9)
    @test Rational{BigInt}(Symbolics.value(Symbolics.substitute(c9, sub1))) == -240693798848 // 81265275
    c9L0 = postnewtonian_coefficient(pnd, 9, 0)
    @test Rational{BigInt}(Symbolics.value(Symbolics.substitute(c9L0, sub1))) == -254196291488 // 81265275
end

@testset "SelfForce: ResummedSeriesData (closed-form ResummedSeries.m)" begin
    p = Symbolics.variable(:p)
    sf = 1 / p^5
    sd = PNSeries(:p, 0, Symbolics.Num[Symbolics.Num(1), Symbolics.Num(2), Symbolics.Num(3)], 0, 3, 1)
    rsd = resummed_series(sf, sd)
    @test rsd isa ResummedSeriesData
    @test isequal(rsd["SingularFactor"], sf)
    @test rsd["SeriesData"] === sd

    # collapse when sd is not a series: sf*sd, NOT a struct.
    res = resummed_series(sf, Symbolics.Num(7))
    @test !(res isa ResummedSeriesData)
    @test isequal(Symbolics.value(Symbolics.simplify(res - sf * 7)), 0)

    # materialize: (1/p^5)(1 + 2p + 3p^2) = p^-5 + 2 p^-4 + 3 p^-3, and == sf*Normal(sd).
    got = Normal(rsd)
    expected = p^(-5) + 2 * p^(-4) + 3 * p^(-3)
    @test isequal(Symbolics.value(Symbolics.simplify(got - expected)), 0)
    @test isequal(Symbolics.value(Symbolics.simplify(got - sf * Normal(sd))), 0)
end

@testset "SelfForce: postnewtonian_expansion round-trip + browse" begin
    # unique query -> a PostNewtonianData whose series matches the direct loader.
    q = ["Schwarzschild", "Circular", "Flux", "Infinity-l2m2", "!Spin"]
    pnd = postnewtonian_expansion(q)
    @test pnd isa PostNewtonianData
    direct = PostNewtonianData(parse_pn_data_file(
        joinpath(SERIES_DATA_DIR, "Schwarzschild", "Circular", "Flux", "Infinity-l2m2.m")))
    @test pn_series(pnd).nmin == pn_series(direct).nmin
    @test pn_series(pnd).nmax == pn_series(direct).nmax
    @test pn_max_order(pnd) == pn_max_order(direct)
    # ambiguous/empty query browses (returns the catalog Vector).
    @test postnewtonian_expansion("") isa Vector{String}
    @test length(postnewtonian_expansion("")) == 91
end

@testset "SelfForce: Schwarzschild energy leading orders {1,0,-1/2,0,3/8} (published)" begin
    pnd = PostNewtonianData(parse_pn_data_file(
        joinpath(SERIES_DATA_DIR, "Kerr", "Spherical", "Orbit", "Energy.m")))
    a = Symbolics.variable(:a)
    a0(c) = Rational{BigInt}(Symbolics.value(Symbolics.substitute(c, Dict(a => Symbolics.Num(big(0))))))
    # var p, den 2: real orders 0, 1/2, 1, 3/2, 2 are coeffs at n = 0,1,2,3,4 (num grid).
    leads = [a0(postnewtonian_coefficient(pnd, n // 2)) for n in 0:4]
    @test leads == Rational{BigInt}[1, 0, -1//2, 0, 3//8]
end

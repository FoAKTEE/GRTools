# Puncture.AlphaIntegration — terminating-2F1 polar-angle α-integration (exactint/modeint/alphaint
# + the vector/tensor combinations alphaintv1..v4, alphaintt1/t2).
#
# Grounding:
#   • CAPTURED upstream oracle: packages/Punctures/test/alphaint_reference.csv (gen_alphaint.wl,
#     N[…,30] at exact dyadic δ ∈ {1/4, 3/4}) copied verbatim into test/fixtures/puncture/. All 9
#     function families characterized at the LEGACY tolerance from Punctures/test/test_alpha_integration.jl
#     (Float64 rtol 1e-6 / atol 1e-7 gross-error guard; BigFloat rtol 1e-25 / atol 1e-28 correctness gate).
#   • Analytic cross-checks (oracle-independent): terminating-2F1 internal === BigSpecial.f21_terminating;
#     modeint m<0 reflection identity; alphaint odd-m lowering-recurrence self-consistency;
#     modeint[p,l,0] == exactint[p,l]; reg(0,·,y)=1.

using Test
using SelfForce
import BigSpecial

const PA = SelfForce.Puncture
const AREF = joinpath(@__DIR__, "fixtures", "puncture", "alphaint_reference.csv")

# δ is stored as an exact rational string ("1/4"); parse to the requested type.
function _pa_parse_delta(::Type{T}, s) where {T}
    if occursin('/', s)
        num, den = split(s, '/')
        return T(parse(Int, num)) / T(parse(Int, den))
    end
    return parse(T, s)
end

function _pa_load_reference(path)
    rows = Tuple{String,Int,Int,Int,String,String}[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        isempty(strip(line)) && continue
        f = split(line, ',')
        push!(rows, (String(f[1]), parse(Int, f[2]), parse(Int, f[3]),
                     parse(Int, f[4]), String(f[5]), String(f[6])))
    end
    return rows
end

function _pa_dispatch(name, p, l, m, δ)
    name == "exactint"   && return PA.exactint(p, l, δ)
    name == "modeint"    && return PA.modeint(p, l, m, δ)
    name == "alphaint"   && return PA.alphaint(p, l, m, δ)
    name == "alphaintv1" && return PA.alphaintv1(p, l, m, δ)
    name == "alphaintv2" && return PA.alphaintv2(p, l, m, δ)
    name == "alphaintv3" && return PA.alphaintv3(p, l, m, δ)
    name == "alphaintv4" && return PA.alphaintv4(p, l, m, δ)
    name == "alphaintt1" && return PA.alphaintt1(p, l, m, δ)
    name == "alphaintt2" && return PA.alphaintt2(p, l, m, δ)
    error("unknown reference function $name")
end

@testset "Puncture.AlphaIntegration" begin
    rows = _pa_load_reference(AREF)

    @testset "reference table is non-trivial (9 families)" begin
        @test length(rows) > 300
        @test length(unique(r[1] for r in rows)) ≥ 9
    end

    @testset "Float64 vs captured oracle (gross-error guard, rtol 1e-6)" begin
        for (name, p, l, m, ds, vs) in rows
            δ = _pa_parse_delta(Float64, ds)
            expected = parse(Float64, vs)
            got = _pa_dispatch(name, p, l, m, δ)
            @test isapprox(float(got), expected; rtol = 1e-6, atol = 1e-7)
        end
    end

    @testset "BigFloat inherits precision (rtol 1e-25)" begin
        setprecision(BigFloat, 256) do
            for (name, p, l, m, ds, vs) in rows
                δ = _pa_parse_delta(BigFloat, ds)
                expected = parse(BigFloat, vs)
                got = _pa_dispatch(name, p, l, m, δ)
                @test got isa BigFloat
                @test isapprox(got, expected;
                               rtol = BigFloat(10)^-25, atol = BigFloat(10)^-28)
            end
        end
    end

    # ---- analytic cross-checks (oracle-independent) ----

    @testset "terminating-2F1 internal === BigSpecial.f21_terminating (collapse)" begin
        for l in (1, 2, 3, 4), c in (-1.5, 0.5, 2.0, 3.5), z in (-0.03125, -0.28125)
            @test PA._f21_term(l, c, z) == BigSpecial.f21_terminating(l, c, z)
        end
    end

    @testset "negative-m reflection identity (modeint)" begin
        for p in (-3, -1, 1), l in (2, 3, 4), mm in (-1, -2), δ in (0.25, 0.6)
            lhs = PA.modeint(p, l, mm, δ)
            rhs = (-1)^mm * factorial(l + mm) / factorial(l - mm) *
                  PA.modeint(p, l, -mm, δ)
            @test isapprox(lhs, rhs; rtol = 1e-12)
        end
    end

    @testset "odd-m lowering recurrence self-consistency (alphaint)" begin
        for p in (-3, -1, 1), l in (2, 3, 4), m in (1, 3), δ in (0.25, 0.6)
            lhs = PA.alphaint(p, l, m, δ)
            rhs = 1 / (2l + 1) *
                  ((l - m + 1) * (l - m + 2) * PA.alphaint(p, l + 1, m - 1, δ) -
                   (l + m - 1) * (l + m) * PA.alphaint(p, l - 1, m - 1, δ))
            @test isapprox(lhs, rhs; rtol = 1e-12)
        end
    end

    @testset "modeint[p,l,0] == exactint[p,l] and reg(0,·,y)=1" begin
        for p in (-7, -5, -3, -1, 1), l in (2, 3, 4), δ in (0.25, 0.6)
            @test isapprox(PA.modeint(p, l, 0, δ), PA.exactint(p, l, δ); rtol = 1e-12)
        end
        @test PA._reg_ycoeffs_n(Float64, 0) == [1.0]
    end
end

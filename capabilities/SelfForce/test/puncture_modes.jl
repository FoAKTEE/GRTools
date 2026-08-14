# Puncture-mode amplitudes hP[i][l,m,r0,M,Δr] + one-sided derivatives dhP (Puncture/PunctureModes.jl).
#
# Grounding (CAPTURE-FIXTURE): the reference tables fixtures/puncture/{hP,dhP}_reference.csv are the
# legacy green Punctures oracle outputs (Mathematica pinned-NIntegrate, captured ONCE and committed
# with the legacy package; copied here verbatim). The port evaluates the same machine-generated
# integrands with QuadGK and must match at the legacy tolerance. Plus the fabrication-proof AD-identity:
# ForwardDiff/central-FD of the already-grounded hP === dhP(n=1) — derivative-of-grounded-field,
# independent of the AD path inside dhP.

using SelfForce.Puncture: hP, dhP

const _HREF = joinpath(@__DIR__, "fixtures", "puncture", "hP_reference.csv")
const _DHREF = joinpath(@__DIR__, "fixtures", "puncture", "dhP_reference.csv")

function _pm_dyadic(::Type{T}, s) where {T}
    if occursin('/', s)
        n, d = split(s, '/'); return T(parse(Int, n)) / T(parse(Int, d))
    end
    return parse(T, s)
end
function _pm_value(::Type{C}, s) where {C<:Complex}
    R = real(C)
    if occursin('|', s)
        a, b = split(s, '|'); return C(parse(R, a), parse(R, b))
    end
    return C(parse(R, s), zero(R))
end

@testset "PunctureModes (hP / dhP)" begin
    @testset "hP vs captured legacy oracle (QuadGK)" begin
        rows = Tuple{Int,Int,Int,Int,String,String}[]
        for (k, line) in enumerate(eachline(_HREF))
            (k == 1 || isempty(strip(line))) && continue
            f = split(line, ',')
            push!(rows, (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]),
                         parse(Int, f[4]), String(f[5]), String(f[6])))
        end
        @test length(rows) > 30
        @test length(unique(r[1] for r in rows)) ≥ 8
        for (i, l, m, r0i, ds, vs) in rows
            Δr = _pm_dyadic(Float64, ds)
            expected = _pm_value(ComplexF64, vs)
            got = ComplexF64(hP(i, l, m, Float64(r0i), 1.0, Δr))
            @test isapprox(got, expected; rtol = 1e-8, atol = 1e-10)
        end
    end

    @testset "hP BigFloat inherits precision (one row per i family)" begin
        rows = Tuple{Int,Int,Int,Int,String,String}[]
        for (k, line) in enumerate(eachline(_HREF))
            (k == 1 || isempty(strip(line))) && continue
            f = split(line, ',')
            push!(rows, (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]),
                         parse(Int, f[4]), String(f[5]), String(f[6])))
        end
        seen = Set{Int}()
        subset = filter(r -> (r[1] in seen ? false : (push!(seen, r[1]); true)), rows)
        setprecision(BigFloat, 160) do
            for (i, l, m, r0i, ds, vs) in subset
                Δr = _pm_dyadic(BigFloat, ds)
                expected = _pm_value(Complex{BigFloat}, vs)
                got = hP(i, l, m, BigFloat(r0i), BigFloat(1), Δr;
                         order = 16, rtol = BigFloat(10)^-18)
                @test got isa Union{BigFloat,Complex{BigFloat}}
                @test isapprox(Complex{BigFloat}(got), expected;
                               rtol = BigFloat(10)^-12, atol = BigFloat(10)^-14)
            end
        end
    end

    @testset "dhP vs captured legacy oracle (QuadGK AD)" begin
        rows = Tuple{Int,Int,Int,Int,String,Int,String}[]
        for (k, line) in enumerate(eachline(_DHREF))
            (k == 1 || isempty(strip(line))) && continue
            f = split(line, ',')
            push!(rows, (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]),
                         parse(Int, f[4]), String(f[5]), parse(Int, f[6]), String(f[7])))
        end
        @test length(rows) ≥ 16
        @test Set(r[6] for r in rows) == Set((1, 2))            # both n present
        for (i, l, m, r0i, ds, n, vs) in rows
            Δr = _pm_dyadic(Float64, ds)
            expected = _pm_value(ComplexF64, vs)
            got = ComplexF64(dhP(i, l, m, Float64(r0i), 1.0, Δr, n))
            @test isapprox(got, expected; rtol = 1e-7, atol = 1e-9)
        end
    end

    @testset "central-FD of hP === dhP(n=1) (derivative-of-grounded-field)" begin
        h = 1e-4
        for (i, l, m) in ((1, 2, 0), (3, 3, 2), (8, 2, 1))
            r0, M, Δr = 10.0, 1.0, 0.5
            fd = (hP(i, l, m, r0, M, Δr + h) - hP(i, l, m, r0, M, Δr - h)) / (2h)
            ad = dhP(i, l, m, r0, M, Δr, 1)
            @test isapprox(ComplexF64(ad), ComplexF64(fd); rtol = 1e-6, atol = 1e-8)
        end
    end

    @testset "n argument validation" begin
        @test_throws ArgumentError dhP(1, 2, 0, 10.0, 1.0, 0.5, 3)
    end
end

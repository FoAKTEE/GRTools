# Puncture.RegSum — regularized l-sum (tail fit) reg_sum/noreg_sum producing value ± uncertainty.
#
# Grounding:
#   • CAPTURED upstream oracle: packages/Punctures/test/regsum_reference.csv (gen_regsum.wl,
#     deterministic rational l-tables, 30-digit) copied into test/fixtures/puncture/. Characterized
#     at the LEGACY tolerance from Punctures/test/test_regsum.jl (value rtol 1e-11 / atol 1e-14;
#     uncertainty rtol 1e-10; BigFloat rtol 1e-13 since the oracle Fit is only ~16-digit).
#   • Analytic (oracle-independent): noreg_sum.value == sum(field) + unc == |last term|; reg_sum of
#     a pure-b2 field totals ~0 and its uncertainty matches the closed-form digamma tail
#     (partial-fraction {1/96,-1/48,1/48,-1/96}) cross-checked vs a brute-force truncated remainder;
#     zero field → zero; BigSpecial.digamma_big closed-form spot value (collapse target).

using Test
using SelfForce
import BigSpecial

const PR = SelfForce.Puncture
const REGREF = joinpath(@__DIR__, "fixtures", "puncture", "regsum_reference.csv")

# field tags must match gen_regsum.wl exactly.
_pr_field(::Val{1}, l) = 1 / (l + 1)^2
_pr_field(::Val{2}, l) = 1 / (2l + 1)^3
_pr_field(::Val{3}, l) =
    1 / ((2l - 1) * (2l + 3)) + 1 / ((2l - 3) * (2l - 1) * (2l + 3) * (2l + 5))
_pr_field(::Val{4}, l) =
    1 / (l + 1)^2 + 7 / ((2l - 1) * (2l + 3)) -
    3 / ((2l - 3) * (2l - 1) * (2l + 3) * (2l + 5))
_pr_table(::Type{T}, tag, lmax) where {T} = T[_pr_field(Val(tag), T(l)) for l in 0:lmax]

function _pr_load_ref(path)
    rows = Tuple{Int,Int,Float64,Float64,Float64,Float64}[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        isempty(strip(line)) && continue
        f = split(line, ',')
        push!(rows, (parse(Int, f[1]), parse(Int, f[2]),
                     parse(Float64, f[3]), parse(Float64, f[4]),
                     parse(Float64, f[5]), parse(Float64, f[6])))
    end
    return rows
end

@testset "Puncture.RegSum" begin
    rows = _pr_load_ref(REGREF)
    @test length(rows) == 12

    @testset "reg_sum / noreg_sum vs captured oracle (Float64)" begin
        for (tag, lmax, rv, ru, nv, nu) in rows
            tbl = _pr_table(Float64, tag, lmax)
            rs = PR.reg_sum(tbl)
            nrs = PR.noreg_sum(tbl)
            @test isapprox(rs.value, rv; rtol = 1e-11, atol = 1e-14)
            @test isapprox(rs.uncertainty, ru; rtol = 1e-10, atol = 1e-16)
            @test isapprox(nrs.value, nv; rtol = 1e-12, atol = 1e-14)
            @test isapprox(nrs.uncertainty, nu; rtol = 1e-12, atol = 1e-16)
        end
    end

    @testset "BigFloat path vs oracle (Fit is ~16-digit, rtol 1e-13)" begin
        setprecision(BigFloat, 256) do
            for (tag, lmax, rv, ru, nv, nu) in rows
                tbl = _pr_table(BigFloat, tag, lmax)
                rs = PR.reg_sum(tbl)
                nrs = PR.noreg_sum(tbl)
                @test rs.value isa BigFloat
                @test isapprox(rs.value, BigFloat(rv); rtol = BigFloat(10)^-13, atol = BigFloat(10)^-15)
                @test isapprox(nrs.value, BigFloat(nv); rtol = BigFloat(10)^-13, atol = BigFloat(10)^-15)
                @test isapprox(rs.uncertainty, BigFloat(ru); rtol = BigFloat(10)^-13, atol = BigFloat(10)^-16)
            end
        end
    end

    @testset "BigFloat precision inheritance (exact-rational noreg, rtol 1e-50)" begin
        _exact(::Val{1}, l) = 1 // (l + 1)^2
        _exact(::Val{2}, l) = 1 // (2l + 1)^3
        _exact(::Val{3}, l) =
            1 // ((2l - 1) * (2l + 3)) + 1 // ((2l - 3) * (2l - 1) * (2l + 3) * (2l + 5))
        _exact(::Val{4}, l) =
            1 // (l + 1)^2 + 7 // ((2l - 1) * (2l + 3)) -
            3 // ((2l - 3) * (2l - 1) * (2l + 3) * (2l + 5))
        setprecision(BigFloat, 256) do
            for tag in 1:4, lmax in (10, 15, 20)
                ref = BigFloat(sum(_exact(Val(tag), big(l)) for l in 0:lmax))
                got = PR.noreg_sum(_pr_table(BigFloat, tag, lmax)).value
                @test isapprox(got, ref; rtol = BigFloat(10)^-50)
            end
        end
    end

    @testset "analytic identities (oracle-independent)" begin
        tbl = _pr_table(Float64, 1, 12)
        @test PR.noreg_sum(tbl).value ≈ sum(tbl)
        @test PR.noreg_sum(tbl).uncertainty == tbl[end]

        b2(l) = 1 / ((2l - 3) * (2l - 1) * (2l + 3) * (2l + 5))
        for lmax in (10, 15, 20)
            tbl3 = _pr_table(Float64, 3, lmax)
            @test abs(PR.reg_sum(tbl3).value) < 1e-12
            brute = sum(b2(l) for l in (lmax+1):200000)
            @test isapprox(PR.reg_sum(tbl3).uncertainty, brute; rtol = 1e-7)
        end

        z = zeros(Float64, 13)
        @test PR.reg_sum(z).value == 0
        @test PR.reg_sum(z).uncertainty == 0
    end

    @testset "digamma collapse onto BigSpecial.digamma_big (closed-form spot values)" begin
        γ = 0.5772156649015328606
        @test isapprox(real(BigSpecial.digamma_big(complex(1.0))), -γ; rtol = 1e-12)
        @test isapprox(real(BigSpecial.digamma_big(complex(0.5))), -γ - 2log(2); rtol = 1e-12)
        # tail closed form on a pure-b2 field matches the digamma partial-fraction sum.
        L = 11.0
        ψ(x) = real(BigSpecial.digamma_big(complex(x)))
        tail = -0.5 * ((1/96)*ψ(L-1.5) - (1/48)*ψ(L-0.5) + (1/48)*ψ(L+1.5) - (1/96)*ψ(L+2.5))
        @test isapprox(PR._regsum_tail(L), tail; rtol = 1e-12)
    end
end

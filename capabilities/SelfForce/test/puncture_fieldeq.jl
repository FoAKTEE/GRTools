# Puncture.FieldEqDeriv — Lorenz wave-operator second radial derivative ddhret[i,l,m] from the
# Barack-Lousto field equations (purely algebraic in h_k, drh_k, dt=-imΩ), plus the even-static m=0
# routing predicate is_evenstatic_ddh / EVENSTATIC_M0_COMPONENTS.
#
# Grounding:
#   • CAPTURED upstream oracle: packages/Punctures/test/ddhret_reference.csv (gen_ddhret.wl runs the
#     GENUINE field-equation Block on one field point with SYNTHETIC complex mode amplitudes) copied
#     into test/fixtures/puncture/. Characterized real+imag across all 10 components at the LEGACY
#     tolerance from Punctures/test/test_subtraction_fieldeq.jl (rtol 1e-10 / atol 1e-12).
#   • Complex{BigFloat} type-generic path agreeing with the Float64 oracle.
#   • Enumerated even-static predicate.

using Test
using SelfForce

const PF = SelfForce.Puncture
const DDHREF = joinpath(@__DIR__, "fixtures", "puncture", "ddhret_reference.csv")

# synthetic mode amplitudes; MUST match gen_ddhret.wl.
_pf_synH(i, l, m) = (i + 1 + l - m) / 7 + im * (2 + i - l + 2m) / 5
_pf_synDH(i, l, m) = (3 + i + 2l + m) / 11 + im * (1 + i + l - m) / 4

function _pf_hvec(l, m)
    h = ComplexF64[_pf_synH(k, l, m) for k in 1:10]
    dh = ComplexF64[_pf_synDH(k, l, m) for k in 1:10]
    return h, dh
end

function _pf_load_ref(path)
    rows = Tuple{Int,Int,Int,Float64,Float64,Float64,Float64}[]
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        isempty(strip(line)) && continue
        f = split(line, ',')
        push!(rows, (parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]),
                     parse(Float64, f[4]), parse(Float64, f[5]),
                     parse(Float64, f[6]), parse(Float64, f[7])))
    end
    return rows
end

@testset "Puncture.FieldEqDeriv" begin

    @testset "field-equation operator vs captured oracle (synthetic modes)" begin
        rows = _pf_load_ref(DDHREF)
        @test length(rows) >= 60
        for (i, l, m, r, r0, re, im_) in rows
            h, dh = _pf_hvec(l, m)
            v = PF.ddhret(i, l, m, r, r0, 1.0, h, dh)
            @test isapprox(real(v), re; rtol = 1e-10, atol = 1e-12)
            @test isapprox(imag(v), im_; rtol = 1e-10, atol = 1e-12)
        end
    end

    @testset "Complex{BigFloat} path agrees with Float64 oracle" begin
        rows = _pf_load_ref(DDHREF)
        for (i, l, m, r, r0, re, im_) in rows[1:5]
            h = Complex{BigFloat}[_pf_synH(k, l, m) for k in 1:10]
            dh = Complex{BigFloat}[_pf_synDH(k, l, m) for k in 1:10]
            v = PF.ddhret(i, l, m, big(r), big(r0), big(1.0), h, dh)
            @test v isa Complex{BigFloat}
            @test isapprox(Float64(real(v)), re; rtol = 1e-10, atol = 1e-12)
        end
    end

    @testset "even-static m=0 routing predicate" begin
        @test PF.is_evenstatic_ddh(1, 2, 0)
        @test PF.is_evenstatic_ddh(7, 4, 0)
        @test PF.is_evenstatic_ddh(2, 2, 0)
        @test !PF.is_evenstatic_ddh(1, 2, 2)   # m != 0
        @test !PF.is_evenstatic_ddh(1, 3, 0)   # odd l
        @test !PF.is_evenstatic_ddh(8, 2, 0)   # odd-parity component
        @test !PF.is_evenstatic_ddh(9, 1, 0)
        @test PF.EVENSTATIC_M0_COMPONENTS == (1, 2, 3, 4, 5, 6, 7)
    end
end

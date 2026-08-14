# Grounded tests for SelfForce.SecondOrder analytic first-order fields (h1_mass, h1_spin, d2hdr2).
#
# GROUNDING (no heavy legacy dep — captured g++ oracle + analytic closed form + field-equation):
#
# (1) h1_mass: the h/dh/ddh log-bearing closed forms (ref-code/SecondOrderRicci/src/h1_M.cc) are
#     ground-truthed by the STANDALONE g++ oracle packages/SecondOrderRicci/test/oracle/
#     h1_mass_oracle.cpp (which lifts h1_M.cc VERBATIM, pure <cmath>). It was compiled with
#       g++ -O2 -o bin h1_mass_oracle.cpp
#     and run on r ∈ {3,4,6,8,20,100}; the printed columns (r h1 dh1 ddh1 h3 dh3 ddh3 h6 dh6 ddh6)
#     are BAKED as the literals below. Comparison is at rtol 1e-12 in the well-conditioned region
#     (these r give f ≥ 1/3). An INDEPENDENT central-FD check (dh≈∂_r h, ddh≈∂_r dh) catches any
#     transcription error shared by oracle and port.
#
# (2) h1_spin: fully analytic — h[8]=-8√(π/3)r⁻², h[9]=-16√(π/3)r⁻³ + exact derivatives, transcribed
#     from h1_S.cc L46-55; closed-form value + exact-derivative FD consistency.
#
# (3) d2hdr2: the genuine upstream expression sliced VERBATIM from ref-code/SecondOrderRicci/src/
#     h1.cc L59-122 into a tiny Boost-multi_array g++ TU (no HDF5/GSL), exercised for ALL i=1..10
#     with several (l,m), m≠0 (so dt=-imΩ and Ω²m² are hit). The oracle was compiled+run ONCE by
#     scratchpad/capture_d2hdr2_oracle.jl (2026-07-02) and its inputs+outputs are BAKED in
#     test/fixtures/d2hdr2_oracle.jl (full provenance header there: g++/boost versions, h1.cc
#     checksum + git state, exact invocation); no live ref-code/toolchain dependency remains.
#     PLUS an independent field-equation
#     self-consistency check: for the static spin/mass fields the trace-reversed solution satisfies
#     the field equation exactly, so d2hdr2 must reproduce ∂²_r h̄ built from the validated h1_*.

using Test
using SelfForce.SecondOrder: h1_mass, h1_spin, d2hdr2, a_il

# ── h1_mass ──────────────────────────────────────────────────────────────────────────────────────
# Captured from h1_mass_oracle.cpp (g++ -O2). Columns: r, then (h,dh,ddh) for i=1, i=3, i=6.
const _H1MASS_ORACLE = Dict(
    3.0   => (0.46456544549436307, 0.47586981485409269, -0.32975761130308473,
             1.5625957553878977, -0.056299806301602251, -2.1530476134920824,
             7.0702210066530133, -3.3984708391429361, 2.3031804302963592),
    4.0   => (0.78332133582218766, 0.1851342087690947, -0.21006530475787547,
             1.1394343708968717, -0.46999280149331263, 0.15666426716443713,
             4.6143590134415078, -1.7233069388088127, 1.0966498701510627),
    6.0   => (0.89705979010766157, -0.015132001307888358, -0.035435428585667228,
             0.52229074983151025, -0.182157430597335, 0.086624005472792739,
             2.5549913585519897, -0.59992880970250184, 0.26069541343327879),
    8.0   => (0.82661954833812656, -0.04504460697608622, -0.0031762399606958167,
             0.28161895399926135, -0.076999422425159308, 0.029796140814555349,
             1.723842987823418, -0.28588511197774263, 0.090721133600725493),
    20.0  => (0.43993613238890111, -0.018567718535769242, 0.0014894253736686348,
             0.036820587136355616, -0.0040816275545365452, 0.00066159537049798815,
             0.56501988229157452, -0.031933052828214327, 0.0036014680382750876),
    100.0 => (0.098121450104198649, -0.00095867728100506233, 1.8703759623186696e-05,
             0.0011422902356746291, -2.3968673311966676e-05, 7.4994824238029089e-07,
             0.10242364500805202, -0.0010470822547940077, 2.1421018562119486e-05),
)

@testset "h1_mass vs baked g++ oracle (well-conditioned r)" begin
    rs = sort(collect(keys(_H1MASS_ORACLE)))
    fld = h1_mass(rs)
    @test Set(keys(fld)) == Set([(1, 0, 0), (3, 0, 0), (6, 0, 0)])
    colof = Dict(1 => 1, 3 => 4, 6 => 7)            # offset of (h,dh,ddh) per i in the tuple
    for (jr, r) in enumerate(rs)
        ex = _H1MASS_ORACLE[r]
        for i in (1, 3, 6)
            c = colof[i]
            @test real(fld[(i, 0, 0)].h[jr])   ≈ ex[c]     rtol = 1e-12
            @test real(fld[(i, 0, 0)].dh[jr])  ≈ ex[c + 1] rtol = 1e-12
            @test real(fld[(i, 0, 0)].ddh[jr]) ≈ ex[c + 2] rtol = 1e-12
            @test imag(fld[(i, 0, 0)].h[jr]) == 0
        end
    end
end

@testset "h1_mass derivative consistency (dh≈∂_r h, ddh≈∂_r dh)" begin
    h = 1e-6
    for i in (1, 3, 6), r in (4.0, 6.0, 8.0, 20.0, 100.0)
        s = h1_mass([r - h, r, r + h])
        dh_num  = (real(s[(i, 0, 0)].h[3])  - real(s[(i, 0, 0)].h[1]))  / (2h)
        ddh_num = (real(s[(i, 0, 0)].dh[3]) - real(s[(i, 0, 0)].dh[1])) / (2h)
        @test real(s[(i, 0, 0)].dh[2])  ≈ dh_num  rtol = 1e-5
        @test real(s[(i, 0, 0)].ddh[2]) ≈ ddh_num rtol = 1e-5
    end
end

@testset "h1_mass BigFloat precision path" begin
    setprecision(BigFloat, 256) do
        fld = h1_mass(BigFloat, [big"8.0"])
        @test fld[(1, 0, 0)].h[1] isa Complex{BigFloat}
        f64 = h1_mass([8.0])
        @test Float64(real(fld[(1, 0, 0)].h[1])) ≈ real(f64[(1, 0, 0)].h[1]) rtol = 1e-12
    end
end

# ── h1_spin ──────────────────────────────────────────────────────────────────────────────────────
@testset "h1_spin closed-form values (analytic)" begin
    grid = [4.0, 5.0, 6.5, 9.0, 14.0, 50.0]
    fld = h1_spin(grid)
    @test Set(keys(fld)) == Set([(8, 1, 0), (9, 1, 0)])
    c = sqrt(pi / 3)
    for (j, r) in enumerate(grid)
        @test fld[(8, 1, 0)].h[j]   ≈ -8c * r^-2  rtol = 1e-14
        @test fld[(8, 1, 0)].dh[j]  ≈ 16c * r^-3  rtol = 1e-14
        @test fld[(8, 1, 0)].ddh[j] ≈ -48c * r^-4 rtol = 1e-14
        @test fld[(9, 1, 0)].h[j]   ≈ -16c * r^-3 rtol = 1e-14
        @test fld[(9, 1, 0)].dh[j]  ≈ 48c * r^-4  rtol = 1e-14
        @test fld[(9, 1, 0)].ddh[j] ≈ -192c * r^-5 rtol = 1e-14
        @test imag(fld[(8, 1, 0)].h[j]) == 0
    end
end

@testset "h1_spin exact-derivative consistency + BigFloat" begin
    h = 1e-5
    for i in (8, 9), r in (4.0, 6.0, 8.0, 10.0)
        s = h1_spin([r - h, r, r + h])
        dh_num  = (s[(i, 1, 0)].h[3]  - s[(i, 1, 0)].h[1])  / (2h)
        ddh_num = (s[(i, 1, 0)].dh[3] - s[(i, 1, 0)].dh[1]) / (2h)
        @test real(s[(i, 1, 0)].dh[2])  ≈ real(dh_num)  rtol = 1e-7
        @test real(s[(i, 1, 0)].ddh[2]) ≈ real(ddh_num) rtol = 1e-7
    end
    setprecision(BigFloat, 256) do
        fld = h1_spin(BigFloat, [big"8.0"])
        @test fld[(8, 1, 0)].h[1] isa Complex{BigFloat}
        @test fld[(8, 1, 0)].h[1] ≈ -8 * sqrt(big(pi) / 3) * big"8.0"^-2 rtol = 1e-60
    end
end

# ── d2hdr2 ───────────────────────────────────────────────────────────────────────────────────────
# Baked one-shot capture of the verbatim h1.cc:59-122 g++ oracle (provenance in the fixture header).
include(joinpath(@__DIR__, "fixtures", "d2hdr2_oracle.jl"))

@testset "d2hdr2 source slice anchor (baked at capture; guards h1.cc line drift)" begin
    # The live-file anchors of the pre-fixture test, asserted against the capture-time record.
    @test occursin("d2hdr2", D2HDR2_SLICE_FIRST_LINE)
    @test occursin("return ddhbar", D2HDR2_SLICE_LAST_LINE)
end

@testset "d2hdr2 vs baked g++ oracle (h1.cc:59-122), all i=1..10, m≠0" begin
    # Same 40 cases the live test generated (MersenneTwister(20240617)); i=1..10 × 4 (l,m) pairs.
    @test length(D2HDR2_ORACLE_CASES) == 40
    @test Set((c[1], c[2], c[3]) for c in D2HDR2_ORACLE_CASES) ==
          Set((i, l, m) for i in 1:10 for (l, m) in ((2, 0), (3, 1), (4, 2), (5, -2)))
    r0 = D2HDR2_ORACLE_R0
    for (i, l, m, r, f, fp, hb, db, oracle) in D2HDR2_ORACLE_CASES
        @test d2hdr2(i, l, m, r, f, fp, r0, hb, db) ≈ oracle rtol = 1e-12
    end
end

@testset "d2hdr2 analytic field-equation self-consistency (static spin/mass)" begin
    rs = [5.0, 6.5, 9.0, 14.0]; M = 1.0
    @testset "h1_spin (l=1,m=0)" begin
        sp = h1_spin(rs)
        for (jr, r) in enumerate(rs)
            f = 1 - 2M / r; fp = 2M / r^2
            hbar = zeros(ComplexF64, 10); dhbar = zeros(ComplexF64, 10); ddbar = zeros(ComplexF64, 10)
            for i in (8, 9)
                a = a_il(i, 1)
                hbar[i] = sp[(i, 1, 0)].h[jr] * r / a
                dhbar[i] = (sp[(i, 1, 0)].dh[jr] * r + sp[(i, 1, 0)].h[jr]) / a
                ddbar[i] = (sp[(i, 1, 0)].ddh[jr] * r + 2 * sp[(i, 1, 0)].dh[jr]) / a
            end
            for i in (8, 9)
                @test d2hdr2(i, 1, 0, r, f, fp, 8.0, hbar, dhbar) ≈ ddbar[i] rtol = 1e-10
            end
        end
    end
    @testset "h1_mass (l=0,m=0; trace-reversal swaps 3↔6)" begin
        ms = h1_mass(rs)
        for (jr, r) in enumerate(rs)
            f = 1 - 2M / r; fp = 2M / r^2
            hbar = zeros(ComplexF64, 10); dhbar = zeros(ComplexF64, 10); ddbar = zeros(ComplexF64, 10)
            for ih in (1, 3, 6)
                i = ih == 3 ? 6 : (ih == 6 ? 3 : ih)
                a = a_il(ih, 0)
                hbar[i] = ms[(ih, 0, 0)].h[jr] * r / a
                dhbar[i] = (ms[(ih, 0, 0)].dh[jr] * r + ms[(ih, 0, 0)].h[jr]) / a
                ddbar[i] = (ms[(ih, 0, 0)].ddh[jr] * r + 2 * ms[(ih, 0, 0)].dh[jr]) / a
            end
            for i in (1, 3, 6)
                @test d2hdr2(i, 0, 0, r, f, fp, 8.0, hbar, dhbar) ≈ ddbar[i] rtol = 1e-9
            end
        end
    end
end

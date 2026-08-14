# Grounded tests for Wigner.coupling (3-harmonic spin-weighted angular coupling).
#
# Promotes SecondOrderRicci.coupling to its canonical Wigner home, re-implemented over the
# in-package wigner3j (so SecondOrderRicci drops its WignerSymbols [deps], design §7.4).
#
# Grounding sources:
#  (1) GSL gsl_sf_coupling_3j ORACLE literals baked verbatim from
#      packages/SecondOrderRicci/test/coupling_test.jl ORACLE[] (15 triples, rtol=1e-13).
#      SecondOrderRicci is HEAVY (HDF5/Kernels) — never live-dep; the literals are copied.
#  (2) Independent re-derivation sgn·norm·3j_s·3j_m via the in-package wigner3j (atol=1e-15).
#  (3) WignerSymbols.wigner3j re-derivation (light TEST-ONLY dep).
#  (4) BigFloat path coupling(BigFloat,0,0,0,1,0,0,1,0,0) ≈ 1/√(4·big(π)).

using Wigner, Test
using WignerSymbols: wigner3j as ws3j

# (lc,mc,sc, la,ma,sa, lb,mb,sb) => GSL-oracle coupling value (coupling_test.jl ORACLE[])
const ORACLE = [
    ((0,0,0, 0,0,0, 0,0,0), 0.28209479177387814),
    ((1,0,0, 0,0,0, 0,0,0), 0.0),
    ((2,1,0, 1,0,0, 1,0,0), 0.0),
    ((2,0,1, 1,0,0, 1,0,0), 0.0),
    ((0,0,0, 1,0,0, 1,0,0), 0.28209479177387814),
    ((2,0,0, 1,0,0, 1,0,0), 0.25231325220201606),
    ((2,0,0, 2,0,0, 2,0,0), 0.18022375157286832),
    ((2,1,0, 1,1,0, 1,0,0), 0.21850968611841584),
    ((2,2,0, 1,1,0, 1,1,0), 0.30901936161855176),
    ((2,0,2, 2,0,1, 2,0,1), 0.22072811544182236),
    ((1,1,0, 1,1,0, 1,0,0), 0.0),
    ((2,0,1, 2,0,1, 2,0,0), 0.090111875786434162),
    ((2,1,1, 2,1,1, 2,0,0), 0.045055937893217088),
    ((3,0,0, 2,0,0, 2,0,0), 0.0),
    ((4,1,0, 2,1,0, 2,0,0), 0.22072811544182283),
]

@testset "coupling vs GSL oracle literals" begin
    for (args, expected) in ORACLE
        got = coupling(args...)
        if expected == 0.0
            @test got == 0.0
        else
            @test got ≈ expected rtol=1e-13
        end
    end
end

# Re-derivation via the in-package wigner3j (proves the slot convention: s-slot vs m-slot).
function coupling_refcheck(lc, mc, sc, la, ma, sa, lb, mb, sb)
    sgn = isodd(mc + sc) ? -1.0 : 1.0
    norm = sqrt((2la + 1) * (2lb + 1) * (2lc + 1) / (4 * pi))
    return sgn * norm *
           wigner3j(Float64, lc, la, lb, sc, -sa, -sb) *
           wigner3j(Float64, lc, la, lb, -mc, ma, mb)
end

@testset "coupling vs in-package wigner3j re-derivation" begin
    for (args, _) in ORACLE
        @test coupling(args...) ≈ coupling_refcheck(args...) atol=1e-15
    end
end

@testset "coupling vs WignerSymbols wigner3j re-derivation" begin
    for (args, _) in ORACLE
        (lc, mc, sc, la, ma, sa, lb, mb, sb) = args
        # guard slots WignerSymbols would throw on (matches GSL zero), as the legacy does
        ok = abs(sc) <= lc && abs(sa) <= la && abs(sb) <= lb &&
             abs(mc) <= lc && abs(ma) <= la && abs(mb) <= lb
        ref = if ok
            sgn = isodd(mc + sc) ? -1.0 : 1.0
            norm = sqrt((2la + 1) * (2lb + 1) * (2lc + 1) / (4 * pi))
            sgn * norm * Float64(ws3j(lc, la, lb, sc, -sa, -sb)) *
                  Float64(ws3j(lc, la, lb, -mc, ma, mb))
        else
            0.0
        end
        @test coupling(args...) ≈ ref atol=1e-15
    end
end

@testset "coupling convention: spin-slot vs m-slot distinct" begin
    a = coupling(2, 1, 0, 1, 1, 0, 1, 0, 0)   # m carried, s zero
    b = coupling(2, 0, 1, 1, 0, 1, 1, 0, 0)   # s carried, m zero
    @test a != b
end

@testset "coupling BigFloat precision path" begin
    setprecision(BigFloat, 256) do
        gb = coupling(BigFloat, 0, 0, 0, 1, 0, 0, 1, 0, 0)
        @test gb isa BigFloat
        @test gb ≈ 1 / sqrt(4 * big(pi)) rtol=1e-60
        @test Float64(gb) ≈ 0.28209479177387814 rtol=1e-13
    end
end

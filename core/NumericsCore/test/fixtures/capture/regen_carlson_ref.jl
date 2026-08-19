# ==============================================================================================================
# NEW REGENERATOR (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# `scratchpad/gen_fixtures.jl` — the one-shot capture named in the header of
# core/NumericsCore/test/fixtures/numericscore_carlson_ref.jl — is LOST (absent from emri-disk/scratchpad,
# from emri-disk git history, and from the GRTools tree).  The ORACLE is still available (legacy
# packages/KerrGeodesics) and the evaluation grid is stated explicitly in test_carlson.jl, so the capture is
# reconstructed here from those two sources and used to RE-VERIFY every baked literal.
#
# fixture : core/NumericsCore/test/fixtures/numericscore_carlson_ref.jl
# oracle  : packages/KerrGeodesics @ e1fc76aa —
#           KG._RC/_RF/_RD/_RJ, _ellipK/_ellipE/_ellipPi, _ellipFinc/_ellipEinc/_ellipPiinc,
#           _phi_reduce, _jacobiAM/_jacobiSN/_jacobiCN/_jacobiDN
# grid    : verbatim from core/NumericsCore/test/test_carlson.jl (lines 22, 41-51, 57, 64-67, 73)
# run     : julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/KerrGeodesics \
#                 core/NumericsCore/test/fixtures/capture/regen_carlson_ref.jl
# env     : GRTOOLS_LEGACY_PACKAGES_ROOT (default /data/haiyangw/claude/emri-disk/packages)
#
# Verification-only by design: it never writes a fixture.  A nonzero deviation is a FINDING —
# report it, do NOT overwrite the fixture.
# ==============================================================================================================

const LEGACY_PACKAGES_ROOT = get(ENV, "GRTOOLS_LEGACY_PACKAGES_ROOT", "/data/haiyangw/claude/emri-disk/packages")

using KerrGeodesics
const KG = KerrGeodesics

module Fix
include(joinpath(dirname(@__DIR__), "numericscore_carlson_ref.jl"))
end

const RES = Dict(:n => 0, :bad => 0, :worst => 0.0)
function chk(name, got, ref)
    d = abs(got - ref); r = iszero(ref) ? d : d / abs(ref)
    RES[:n] += 1; RES[:worst] = max(RES[:worst], r)
    d == 0 || (RES[:bad] += 1; println("  DIFF ", name, "  got=", got, "  ref=", ref, "  |Δ|=", d, "  rel=", r))
end

cx, cy, cz, cp = 1.0 + 2.0im, 2.0 - 1.0im, 3.0 + 0.5im, 0.7 + 0.3im
chk("RC_cxy",   KG._RC(cx, cy),         Fix.REF_RC_cxy)
chk("RF_cxyz",  KG._RF(cx, cy, cz),     Fix.REF_RF_cxyz)
chk("RD_cxyz",  KG._RD(cx, cy, cz),     Fix.REF_RD_cxyz)
chk("RJ_cxyzp", KG._RJ(cx, cy, cz, cp), Fix.REF_RJ_cxyzp)

chk("RC_23",   KG._RC(2.0 + 0im, 3.0 + 0im),                       Fix.REF_RC_23)
chk("RF_123",  KG._RF(1.0 + 0im, 2.0 + 0im, 3.0 + 0im),            Fix.REF_RF_123)
chk("RD_123",  KG._RD(1.0 + 0im, 2.0 + 0im, 3.0 + 0im),            Fix.REF_RD_123)
chk("RJ_1234", KG._RJ(1.0 + 0im, 2.0 + 0im, 3.0 + 0im, 4.0 + 0im), Fix.REF_RJ_1234)

for m in (0.001, 0.3, 0.7, 0.999, 2.4972)
    chk("ellipK[$m]", KG._ellipK(m), Fix.REF_ellipK[m])
    chk("ellipE[$m]", KG._ellipE(m), Fix.REF_ellipE[m])
end
for (n, m) in ((0.2, 0.3), (0.6, 0.7))
    chk("ellipPi[($n,$m)]", KG._ellipPi(n, m), Fix.REF_ellipPi[(n, m)])
end
chk("ellipPi_nC", KG._ellipPi(0.5 + 0.2im, 0.4), Fix.REF_ellipPi_nC)
for k in (0.3, 0.6, 0.9)
    chk("ellipK_mod[$k]", real(KG._ellipK(k^2)), Fix.REF_ellipK_mod[k])
    chk("ellipE_mod[$k]", real(KG._ellipE(k^2)), Fix.REF_ellipE_mod[k])
end
for phi in (0.7, 4.0, -1.3)
    chk("ellipFinc[$phi]",  KG._ellipFinc(phi, 0.3),       Fix.REF_ellipFinc[phi])
    chk("ellipEinc[$phi]",  KG._ellipEinc(phi, 0.3),       Fix.REF_ellipEinc[phi])
    chk("ellipPiinc[$phi]", KG._ellipPiinc(0.2, phi, 0.3), Fix.REF_ellipPiinc[phi])
end
let got = KG._phi_reduce(4.0)
    RES[:n] += 1
    got == Fix.REF_phi_reduce_4 || (RES[:bad] += 1; println("  DIFF phi_reduce_4 got=", got, " ref=", Fix.REF_phi_reduce_4))
end
for (u, m) in ((0.8, 0.5), (1.7, 0.2), (2.3, 0.9))
    chk("jacobiAM[($u,$m)]", KG._jacobiAM(u, m), Fix.REF_jacobiAM[(u, m)])
    chk("jacobiSN[($u,$m)]", KG._jacobiSN(u, m), Fix.REF_jacobiSN[(u, m)])
    chk("jacobiCN[($u,$m)]", KG._jacobiCN(u, m), Fix.REF_jacobiCN[(u, m)])
    chk("jacobiDN[($u,$m)]", KG._jacobiDN(u, m), Fix.REF_jacobiDN[(u, m)])
end

println("RESULT: literals checked = ", RES[:n], " ; differing = ", RES[:bad],
        " ; worst relative deviation = ", RES[:worst])

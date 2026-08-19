# ==============================================================================================================
# NEW REGENERATOR (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# The original capture script for the legacy-ScatteringPhaseMST fixture is LOST (the `leg_*.jl` headers say
# "regenerate via scratchpad gen.jl against the same legacy commit"; no `gen.jl` exists in emri-disk/scratchpad,
# in emri-disk git history, or anywhere in the GRTools tree).  The ORACLE is still available (legacy
# `packages/ScatteringPhaseMST`) and each of the four tables is a Dict keyed by its own evaluation grid, so the
# capture is reconstructed here: read the BAKED KEYS out of the committed fixture, re-evaluate the legacy oracle
# at exactly those keys, report the max abs/rel deviation per table.
#
# fixture covered
#   capabilities/Perturbation/test/fixtures/leg_scattering.jl   (REF_STAR, REF_BREFL, REF_MSTC, REF_SPN)
# oracle
#   packages/ScatteringPhaseMST @ e1fc76aa
#     starobinsky_constant(s,l,m,q,eps; P)            <- REF_STAR
#     _brefl_over_binc(s,l,m,q,eps; n=10, nUp=5)      <- REF_BREFL   (kwargs from test_scattering.jl:30-31)
#     mst_coefficient(nn,s,l,m,q,eps)                 <- REF_MSTC    (defaults, test_scattering.jl:34)
#     scattering_phase_numerical(s,l,m,q,eps; n=40, nUp=5, P)  <- REF_SPN  (test_scattering.jl:42)
#   The non-key kwargs live in the test bodies, not in the keys; the line references above are the source of
#   truth for them.
# run
#   julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/ScatteringPhaseMST \
#         capabilities/Perturbation/test/fixtures/capture/regen_leg_scattering.jl
# env
#   GRTOOLS_LEGACY_PACKAGES_ROOT  (default /data/haiyangw/claude/emri-disk/packages)
#   GRTOOLS_ROOT                  (default: the GRTools checkout containing this file)
#
# Verification-only by design: it never writes a fixture.  A nonzero max deviation is a FINDING —
# report it, do NOT overwrite the fixture.
# ==============================================================================================================

const LEGACY_PACKAGES_ROOT = get(ENV, "GRTOOLS_LEGACY_PACKAGES_ROOT", "/data/haiyangw/claude/emri-disk/packages")
const GRTOOLS_ROOT         = get(ENV, "GRTOOLS_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..")))

using ScatteringPhaseMST
const SPM = ScatteringPhaseMST

module ScatFix
include(joinpath(dirname(@__DIR__), "leg_scattering.jl"))   # REF_STAR / REF_BREFL / REF_MSTC / REF_SPN
end

worst = 0.0
nbad = 0
nchecked = 0

function check(name, tbl, f)
    global worst, nbad, nchecked
    md = 0.0; mr = 0.0; k_worst = nothing
    for (k, v) in tbl
        got = f(k)
        nchecked += 1
        d = abs(got - v)
        r = iszero(v) ? d : d / abs(v)
        if d > md; md = d; k_worst = k; end
        if r > mr; mr = r; end
        d == 0 || (nbad += 1)
    end
    worst = max(worst, mr)
    println(rpad(name, 12), " n=", lpad(length(tbl), 3),
            "  max|Δ|=", md, "  max rel=", mr, "  (worst key ", k_worst === nothing ? "nothing" : k_worst, ")")
end

check("REF_STAR",  ScatFix.REF_STAR,  k -> SPM.starobinsky_constant(k[1], k[2], k[3], k[4], k[5]; P=k[6]))
check("REF_BREFL", ScatFix.REF_BREFL, k -> SPM._brefl_over_binc(k[1], k[2], k[3], k[4], k[5]; n=10, nUp=5))
check("REF_MSTC",  ScatFix.REF_MSTC,  k -> SPM.mst_coefficient(k[1], k[2], k[3], k[4], k[5], k[6]))
check("REF_SPN",   ScatFix.REF_SPN,   k -> SPM.scattering_phase_numerical(k[1], k[2], k[3], k[4], k[5];
                                                                          n=40, nUp=5, P=k[6]))

println("RESULT: literals checked = ", nchecked, " ; differing = ", nbad,
        " ; worst relative deviation = ", worst)

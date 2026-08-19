# ==============================================================================================================
# NEW REGENERATOR (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# The original capture script for the legacy point-particle fixture is LOST (the `leg_*.jl` headers say
# "regenerate via scratchpad gen.jl against the same legacy commit"; no `gen.jl` exists in emri-disk/scratchpad,
# in emri-disk git history, or anywhere in the GRTools tree).  The ORACLE is still available (legacy
# `packages/Teukolsky`), the three Dict tables are keyed by their own evaluation grid and the two NamedTuple
# tables have their arguments spelled out in the test body, so the capture is reconstructed here: re-build every
# baked mode from the legacy constructors and report the max abs/rel deviation per field.
#
# fixture covered
#   capabilities/Perturbation/test/fixtures/leg_pointparticle.jl  (REF_PP_CIRC, REF_PP_ECC, REF_PP_SPH, REF_PP_GEN)
# oracle
#   packages/Teukolsky @ e1fc76aa
#     TeukolskyPointParticleModeCircular(s,l,m,a,p)                        <- REF_PP_CIRC
#     TeukolskyPointParticleModeEccentric(l,m,n,a,p,e; s=s)                <- REF_PP_ECC
#     TeukolskyPointParticleModeSpherical(2,2,0,0.5,8.0,0.9)               <- REF_PP_SPH  (test_pointparticle.jl:42)
#     TeukolskyPointParticleModeGeneric(2,2,0,0,0.5,8.0,0.1,0.9)           <- REF_PP_GEN  (test_pointparticle.jl:48)
# run
#   julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/Teukolsky \
#         capabilities/Perturbation/test/fixtures/capture/regen_leg_pointparticle.jl
# env
#   GRTOOLS_LEGACY_PACKAGES_ROOT  (default /data/haiyangw/claude/emri-disk/packages)
#   GRTOOLS_ROOT                  (default: the GRTools checkout containing this file)
#
# Verification-only by design: it never writes a fixture.  A nonzero max deviation is a FINDING —
# report it, do NOT overwrite the fixture.
# ==============================================================================================================

const LEGACY_PACKAGES_ROOT = get(ENV, "GRTOOLS_LEGACY_PACKAGES_ROOT", "/data/haiyangw/claude/emri-disk/packages")
const GRTOOLS_ROOT         = get(ENV, "GRTOOLS_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..")))

using Teukolsky

module PPFix
include(joinpath(dirname(@__DIR__), "leg_pointparticle.jl"))   # REF_PP_CIRC / REF_PP_ECC / REF_PP_SPH / REF_PP_GEN
end

worst = 0.0
nbad = 0
nchecked = 0

function cmp!(label, got, want)
    global worst, nbad, nchecked
    nchecked += 1
    d = abs(got - want)
    r = iszero(want) ? d : d / abs(want)
    worst = max(worst, r)
    if d != 0
        nbad += 1
        println("  DIFF ", label, "  baked=", want, "  legacy=", got, "  abs=", d, "  rel=", r)
    end
end

println("REF_PP_CIRC — TeukolskyPointParticleModeCircular")
for (key, ref) in sort(collect(PPFix.REF_PP_CIRC), by = k -> string(k))
    s, l, m, a, p = key
    md = Teukolsky.TeukolskyPointParticleModeCircular(s, l, m, a, p)
    cmp!("$key En",    md.En,                       ref.En)
    cmp!("$key Lz",    md.Lz,                       ref.Lz)
    cmp!("$key Upst",  md.Upst,                     ref.Upst)
    cmp!("$key omega", md.omega,                    ref.omega)
    cmp!("$key lam",   md.lambda,                   ref.lambda)
    cmp!("$key Zinf",  md.Zinf,                     ref.Zinf)
    cmp!("$key efI",   md.EnergyFlux.I,             ref.efI)
    cmp!("$key efH",   md.EnergyFlux.H,             ref.efH)
    cmp!("$key jfI",   md.AngularMomentumFlux.I,    ref.jfI)
end

println("REF_PP_ECC — TeukolskyPointParticleModeEccentric")
for (key, ref) in sort(collect(PPFix.REF_PP_ECC), by = k -> string(k))
    s, l, m, n, a, p, e = key
    md = Teukolsky.TeukolskyPointParticleModeEccentric(l, m, n, a, p, e; s=s)
    cmp!("$key omega", md.omega,                    ref.omega)
    cmp!("$key lam",   md.lambda,                   ref.lambda)
    cmp!("$key Zinf",  md.Zinf,                     ref.Zinf)
    cmp!("$key efI",   md.EnergyFlux.I,             ref.efI)
    cmp!("$key jfI",   md.AngularMomentumFlux.I,    ref.jfI)
end

println("REF_PP_SPH / REF_PP_GEN — Spherical + Generic")
let sph = Teukolsky.TeukolskyPointParticleModeSpherical(2, 2, 0, 0.5, 8.0, 0.9),
    gen = Teukolsky.TeukolskyPointParticleModeGeneric(2, 2, 0, 0, 0.5, 8.0, 0.1, 0.9)
    cmp!("SPH omega", sph.omega,        PPFix.REF_PP_SPH.omega)
    cmp!("SPH Zinf",  sph.Zinf,         PPFix.REF_PP_SPH.Zinf)
    cmp!("SPH efI",   sph.EnergyFlux.I, PPFix.REF_PP_SPH.efI)
    cmp!("GEN omega", gen.omega,        PPFix.REF_PP_GEN.omega)
    cmp!("GEN Zinf",  gen.Zinf,         PPFix.REF_PP_GEN.Zinf)
    cmp!("GEN efI",   gen.EnergyFlux.I, PPFix.REF_PP_GEN.efI)
end

println("RESULT: literals checked = ", nchecked, " ; differing = ", nbad,
        " ; worst relative deviation = ", worst)

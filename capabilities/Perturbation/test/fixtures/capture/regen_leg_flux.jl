# ==============================================================================================================
# NEW REGENERATOR (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# The original capture script for the legacy-Teukolsky flux fixture is LOST: every `leg_*.jl` header
# says "regenerate via scratchpad gen.jl against the same legacy commit", and no file named `gen.jl`
# exists in emri-disk/scratchpad, in emri-disk git history, or anywhere in the GRTools tree.
# The ORACLE is still available (legacy `packages/Teukolsky`), and this fixture is a Dict keyed by its
# own evaluation grid, so the capture can be reconstructed from the fixture + the test body: this script
# reads the BAKED KEYS out of the committed fixture, reads the (lambda, Z^inf, Z^H) sample list out of
# `test/test_fluxes.jl` (the keys alone do not carry them), re-evaluates the legacy oracle at exactly
# those arguments, and reports the max abs/rel deviation.
#
# fixture covered
#   capabilities/Perturbation/test/fixtures/leg_flux.jl                (REF_FLUX_LEG)
# oracle
#   packages/Teukolsky @ e1fc76aa — energy_flux(s,l,m,a,omega,lambda,ZI,ZH),
#                                   angular_momentum_flux(s,l,m,a,omega,lambda,ZI,ZH)  (src/fluxes.jl)
# run
#   julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/Teukolsky \
#         capabilities/Perturbation/test/fixtures/capture/regen_leg_flux.jl
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

# --- the baked fixture (keys + reference values) ---------------------------------------------------
module FluxFix
include(joinpath(dirname(@__DIR__), "leg_flux.jl"))   # const REF_FLUX_LEG
end

# --- the (lambda, ZI, ZH) sample list, read out of the test body (single source of truth) -----------
const TESTFILE = joinpath(@__DIR__, "..", "..", "test_fluxes.jl")
const _src = read(TESTFILE, String)
const _m = match(r"samples\s*=\s*\[(.*?)\]\n"s, _src)
_m === nothing && error("could not locate the `samples = [...]` block in $TESTFILE")
const SAMPLES = eval(Meta.parse("[" * _m.captures[1] * "]"))

worst = 0.0
nbad = 0
nchecked = 0

for (s, l, m, a, ω, λ, ZI, ZH) in SAMPLES
    key = (s, l, m, a, ω)
    haskey(FluxFix.REF_FLUX_LEG, key) || error("fixture has no baked entry for $key")
    ref = FluxFix.REF_FLUX_LEG[key]
    efn = Teukolsky.energy_flux(s, l, m, a, ω, λ, ZI, ZH)
    jfn = Teukolsky.angular_momentum_flux(s, l, m, a, ω, λ, ZI, ZH)
    for (name, got, want) in (("eI", efn.I, ref.eI), ("eH", efn.H, ref.eH),
                              ("jI", jfn.I, ref.jI), ("jH", jfn.H, ref.jH))
        global worst, nbad, nchecked
        nchecked += 1
        d = abs(got - want)
        r = iszero(want) ? d : d / abs(want)
        worst = max(worst, r)
        if d != 0
            nbad += 1
            println("  DIFF ", key, " ", name, "  baked=", want, "  legacy=", got,
                    "  abs=", d, "  rel=", r)
        end
    end
end

println("RESULT: literals checked = ", nchecked, " ; differing = ", nbad,
        " ; worst relative deviation = ", worst)

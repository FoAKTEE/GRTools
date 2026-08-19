# ==============================================================================================================
# NEW REGENERATOR (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# The original capture scripts for the spin-weighted spherical/spheroidal reference fixtures
# (`scratchpad/capture_swsh.jl` for Harmonics, `scratchpad/gen_fixtures.jl` for NumericsCore) are LOST:
# neither is present in emri-disk/scratchpad, in emri-disk git history, or anywhere in the GRTools tree.
# The ORACLE, however, is still available (legacy `packages/SpinWeightedSpheroidalHarmonics`), and every
# affected fixture is a Dict keyed by its own evaluation grid.  This script therefore reconstructs the
# capture from the fixture headers' own description: it reads the BAKED KEYS out of the committed fixture,
# re-evaluates the legacy oracle at exactly those keys, and reports the max abs/rel deviation per table.
#
# fixtures covered
#   capabilities/Harmonics/test/fixtures/spherical_harmonics_ref.jl   (FIXSPH_Y)
#   core/NumericsCore/test/fixtures/numericscore_swsh_ref.jl          (REF_Y, REF_Y_highl, REF_ydtheta)
# oracle
#   packages/SpinWeightedSpheroidalHarmonics @ e1fc76aa
#     SpinWeightedSphericalHarmonicY(s,l,m,θ,φ)  and  _y_dtheta(s,l,m,θ,1)
# run
#   julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/SpinWeightedSpheroidalHarmonics \
#         capabilities/Harmonics/test/fixtures/capture/regen_swsh_ref.jl
# env
#   GRTOOLS_LEGACY_PACKAGES_ROOT  (default /data/haiyangw/claude/emri-disk/packages)
#   GRTOOLS_ROOT                  (default: the GRTools checkout containing this file)
#
# Verification-only by design: it never writes a fixture.  A nonzero max deviation is a FINDING —
# report it, do NOT overwrite the fixture.
# ==============================================================================================================

const LEGACY_PACKAGES_ROOT = get(ENV, "GRTOOLS_LEGACY_PACKAGES_ROOT", "/data/haiyangw/claude/emri-disk/packages")
const GRTOOLS_ROOT         = get(ENV, "GRTOOLS_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..")))

using SpinWeightedSpheroidalHarmonics
const SWSH = SpinWeightedSpheroidalHarmonics

module HarmonicsFix
include(joinpath(dirname(@__DIR__), "spherical_harmonics_ref.jl"))
end
module NumericsFix
const F = joinpath(get(ENV, "GRTOOLS_ROOT",
                       normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))),
                   "core", "NumericsCore", "test", "fixtures", "numericscore_swsh_ref.jl")
include(F)
end

worst = 0.0
nbad = 0

function check(name, tbl, f)
    global worst, nbad
    md = 0.0; mr = 0.0; k_worst = nothing
    for (k, v) in tbl
        got = f(k)
        d = abs(got - v)
        r = iszero(v) ? d : d / abs(v)
        if d > md; md = d; k_worst = k; end
        mr = max(mr, r)
    end
    worst = max(worst, mr)
    md == 0 || (nbad += 1)
    println(rpad(name, 26), " n=", lpad(length(tbl), 4),
            "  max|Δ|=", md, "  max rel=", mr, "  (worst key ", k_worst, ")")
end

# Harmonics FIXSPH_Y :: Dict{(s,l,m,θ,φ)::(Int,Int,Int,Float64,Float64)} -> ComplexF64
check("Harmonics FIXSPH_Y", HarmonicsFix.FIXSPH_Y,
      k -> SWSH.SpinWeightedSphericalHarmonicY(k[1], k[2], k[3], k[4], k[5]))
# NumericsCore tables :: keys are all-Float64 tuples
check("NumericsCore REF_Y", NumericsFix.REF_Y,
      k -> SWSH.SpinWeightedSphericalHarmonicY(Int(k[1]), Int(k[2]), Int(k[3]), k[4], k[5]))
check("NumericsCore REF_Y_highl", NumericsFix.REF_Y_highl,
      k -> SWSH.SpinWeightedSphericalHarmonicY(Int(k[1]), Int(k[2]), Int(k[3]), k[4], k[5]))
check("NumericsCore REF_ydtheta", NumericsFix.REF_ydtheta,
      k -> SWSH._y_dtheta(Int(k[1]), Int(k[2]), Int(k[3]), k[4], 1))

println("\nRESULT: tables with a nonzero deviation = ", nbad, " ; worst relative deviation = ", worst)

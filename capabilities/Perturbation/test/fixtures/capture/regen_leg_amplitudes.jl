# ==============================================================================================================
# NEW REGENERATOR (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# The original capture script for the legacy-Teukolsky MST amplitude fixture is LOST (every `leg_*.jl`
# header says "regenerate via scratchpad gen.jl against the same legacy commit"; no `gen.jl` exists in
# emri-disk/scratchpad, in emri-disk git history, or anywhere in the GRTools tree).  The ORACLE is still
# available (legacy `packages/Teukolsky`) and the fixture is keyed by its own evaluation grid, so the
# capture is reconstructed here: read the BAKED KEYS out of the committed fixture, re-evaluate the legacy
# oracle at exactly those keys, report the max abs/rel deviation per table.
#
# fixture covered
#   capabilities/Perturbation/test/fixtures/leg_amplitudes.jl   (REF_AMP1, REF_AMP2, REF_AMP3)
# oracle
#   packages/Teukolsky @ e1fc76aa
#     REF_AMP1 : TeukolskyRadial(s,l,m,a,omega)  ->  In/Up .nu .lambda .amplitudes.{Incidence,Reflection}
#                (amplitudes are already normalized to unit transmission, src/teukolsky_radial.jl:54)
#     REF_AMP2 : _mst_amplitudes(s,m,a,omega,nu,lambda; Nmax=80) at the LEGACY-recomputed (nu,lambda)
#     REF_AMP3 : the (lambda,nu) pair of the s=-2,l=2,m=2,a=0.5,omega=0.2 mode
# run
#   julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/Teukolsky \
#         capabilities/Perturbation/test/fixtures/capture/regen_leg_amplitudes.jl
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

module AmpFix
include(joinpath(dirname(@__DIR__), "leg_amplitudes.jl"))   # REF_AMP1 / REF_AMP2 / REF_AMP3
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

println("REF_AMP1 — TeukolskyRadial normalized amplitudes")
for (key, ref) in sort(collect(AmpFix.REF_AMP1), by = k -> string(k[1]))
    s, l, m, a, ω = key
    R = Teukolsky.TeukolskyRadial(s, l, m, a, ω)
    cmp!("$key nu",    R.In.nu,                      ref.nu)
    cmp!("$key lam",   R.In.lambda,                  ref.lam)
    cmp!("$key inInc", R.In.amplitudes.Incidence,    ref.inInc)
    cmp!("$key inRef", R.In.amplitudes.Reflection,   ref.inRef)
    cmp!("$key upInc", R.Up.amplitudes.Incidence,    ref.upInc)
    cmp!("$key upRef", R.Up.amplitudes.Reflection,   ref.upRef)
end

println("REF_AMP2 — raw _mst_amplitudes at the legacy-recomputed (nu, lambda)")
for (key, ref) in sort(collect(AmpFix.REF_AMP2), by = k -> string(k[1]))
    s, l, m, a, ω = key
    R = Teukolsky.TeukolskyRadial(s, l, m, a, ω)
    ν = R.In.nu; λ = R.In.lambda
    cmp!("$key lam", λ, ref.lam)
    cmp!("$key nu",  ν, ref.nu)
    amp = Teukolsky._mst_amplitudes(s, m, a, ω, ν, λ; Nmax=80)
    for k in (:InTrans, :InInc, :InRef, :UpTrans, :UpInc, :UpRef, :Aplus, :Aminus)
        cmp!("$key $k", getfield(amp, k), getfield(ref, k))
    end
end

println("REF_AMP3 — (lambda, nu) of the s=-2,l=2,m=2,a=0.5,omega=0.2 mode")
let R = Teukolsky.TeukolskyRadial(-2, 2, 2, 0.5, 0.2)
    cmp!("REF_AMP3.lam64", real(R.In.lambda), AmpFix.REF_AMP3.lam64)
    cmp!("REF_AMP3.nu64",  real(R.In.nu),     AmpFix.REF_AMP3.nu64)
end

println("RESULT: literals checked = ", nchecked, " ; differing = ", nbad,
        " ; worst relative deviation = ", worst)

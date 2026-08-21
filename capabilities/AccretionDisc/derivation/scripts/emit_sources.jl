# scripts/emit_sources.jl — route R3 stage 4 (EGD::jd-codegen): emit the GENERATED sources of
# `AccretionDisc.SpiralWaves` from the serialized R3 expressions, then run the exact round-trip
# self-check on the emitted TEXT.
#
#   julia --project=capabilities/AccretionDisc/derivation \
#         capabilities/AccretionDisc/derivation/scripts/emit_sources.jl [npts] [outdir]
#
# Deterministic: the same artefacts produce byte-identical files (asserted below by re-emitting into
# a temporary directory and comparing the sha256 of both runs).

using Printf
using SHA

include(joinpath(@__DIR__, "..", "src", "R3Derivation.jl"))
using .R3Derivation
const R3 = R3Derivation

const NPTS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
const HERE = normpath(joinpath(@__DIR__, ".."))
const OUTDIR = length(ARGS) >= 2 ? ARGS[2] : normpath(joinpath(HERE, "..", "src", "SpiralWaves"))
const MASTER_ART = joinpath(HERE, "out", "R3_master.jls.txt")
const DIAG_ART = joinpath(HERE, "out", "R3_diagnostics.jls.txt")

_sha(p) = bytes2hex(SHA.sha256(read(p)))

println("=== route R3 stage 4 — Julia code emission ===")
println("julia      : ", VERSION)
println("artefacts  : R3_master.jls.txt  sha256 ", _sha(MASTER_ART))
println("             R3_diagnostics.jls.txt  sha256 ", _sha(DIAG_ART))
println("out dir    : ", OUTDIR)
println()

t0 = time()
mpath = R3.emit_master(joinpath(OUTDIR, "master_derived.jl"), MASTER_ART)
dpath = R3.emit_diagnostics(joinpath(OUTDIR, "r2_derived.jl"), MASTER_ART, DIAG_ART)
println(@sprintf("emit wall  : %.1f s", time() - t0))
println("wrote ", mpath, "  (", filesize(mpath), " bytes)  sha256 ", _sha(mpath))
println("wrote ", dpath, "  (", filesize(dpath), " bytes)  sha256 ", _sha(dpath))
println()

# --- determinism: a second emission into a scratch directory must be byte-identical -------------
tmp = mktempdir()
m2 = R3.emit_master(joinpath(tmp, "master_derived.jl"), MASTER_ART)
d2 = R3.emit_diagnostics(joinpath(tmp, "r2_derived.jl"), MASTER_ART, DIAG_ART)
det = (_sha(mpath) == _sha(m2)) && (_sha(dpath) == _sha(d2))
println("determinism: second emission byte-identical = ", det,
        "   (", _sha(m2), ", ", _sha(d2), ")")
rm(tmp; recursive = true)
det || error("codegen: emission is NOT deterministic")
println()

# --- exact round-trip of the emitted TEXT --------------------------------------------------------
println("--- round-trip self-check (emitted Julia re-parsed, interpreted over Rational{BigInt}) ---")
c = R3.codegen_selfcheck(mpath, dpath, MASTER_ART, DIAG_ART; npts = NPTS)
println()
print(R3.certificate_table([c]))
c.pass || error("codegen: round-trip self-check FAILED")
println()
println(@sprintf("=== SUMMARY === emission + round-trip OK   total wall %.1f s", time() - t0))

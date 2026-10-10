# Re-derive the default artifacts, emit with the existing printer, compare SHA256
# against the frozen artifacts AND shipped runtime files. Never writes the latter.
# Usage: julia --project=capabilities/AccretionDisc/derivation <this file> OUTDIR
using SHA
include(joinpath(@__DIR__, "..", "src", "R3Derivation.jl"))
const R3 = R3Derivation
length(ARGS) == 1 || error("usage: check_li_adiabatic_bytes.jl OUTDIR")
out = abspath(ARGS[1])
mkpath(out)
here = normpath(joinpath(@__DIR__, ".."))
ok = true
for (name, exports) in (("R3_master.jls.txt", R3.master_exports()),
                        ("R3_diagnostics.jls.txt", R3.diagnostics_exports()))
    ref = joinpath(here, "out", name)
    # Reuse only the historical header; every expression is freshly derived.
    lines = readlines(ref)
    first_manifest = findfirst(==("# R3 symbol manifest — arXiv:2608.17003, route R3 (Julia-native derivation)"), lines)
    header = join(lines[1:first_manifest-1], "\n")
    path = R3.write_exports(joinpath(out, name), exports; header)
    same = sha256(read(path)) == sha256(read(ref))
    global ok &= same
    println(same ? "[PASS] " : "[FAIL] ", "fresh adiabatic artifact ", name,
            " sha256=", bytes2hex(sha256(read(path))), " reference=", bytes2hex(sha256(read(ref))))
end
mart = joinpath(out, "R3_master.jls.txt")
dart = joinpath(out, "R3_diagnostics.jls.txt")
R3.emit_master(joinpath(out, "master_derived.jl"), mart)
R3.emit_diagnostics(joinpath(out, "r2_derived.jl"), mart, dart)
for name in ("master_derived.jl", "r2_derived.jl")
    got = read(joinpath(out, name))
    ref = read(joinpath(here, "..", "src", "SpiralWaves", name))
    same = sha256(got) == sha256(ref)
    global ok &= same
    println(same ? "[PASS] " : "[FAIL] ", "adiabatic emitted bytes ", name,
            " sha256=", bytes2hex(sha256(got)), " reference=", bytes2hex(sha256(ref)))
    if !same
        a, b = split(String(got), '\n'), split(String(ref), '\n')
        for k in eachindex(a, b)
            a[k] == b[k] && continue
            println("  differing line ", k, " emitted: ", a[k], "\n  reference: ", b[k])
        end
    end
end
c = R3.codegen_selfcheck(joinpath(out, "master_derived.jl"), joinpath(out, "r2_derived.jl"),
                         mart, dart; npts = 20)
ok && c.pass || error("adiabatic byte-identity gate failed; see differing lines above")
println("[PASS] adiabatic re-derivation and re-emission are byte-identical")

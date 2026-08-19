#!/usr/bin/env julia
# check_provenance.jl — re-verify every baked data asset against its PROVENANCE.md record.
#
# Each vendored data directory carries a PROVENANCE.md whose table rows are
#     | `<path relative to that directory>` | <bytes> | `<sha256>` | …
# This script re-hashes every recorded file and fails (exit 1) on a missing file, a size
# mismatch or a hash mismatch — so "the shipped fixtures are the ones the record names" is a
# CHECKED claim, not a comment.
#
# Usage:  julia tools/check_provenance.jl            # from the repository root
#
# Adding a data directory: give it a PROVENANCE.md in the same table format and list it below.

using SHA: sha256

const ROOT = normpath(joinpath(@__DIR__, ".."))

const PROVENANCE_DIRS = [
    "capabilities/PostNewtonian/SeriesData",
    "capabilities/ModifiedGravity/test/data",
    "capabilities/SelfForce/test/fixtures/puncture",
    "core/Wigner/test",
    "capabilities/Hybrid/assets",
]

# `| `rel/path` | 12345 | `<64 hex>` |` — trailing columns (e.g. the generator name) are ignored.
const ROW = r"^\|\s*`([^`]+)`\s*\|\s*(\d+)\s*\|\s*`([0-9a-f]{64})`"

function check_dir(dir::AbstractString)
    md = joinpath(ROOT, dir, "PROVENANCE.md")
    isfile(md) || (println("MISSING  $dir/PROVENANCE.md"); return (0, 1))
    checked = 0; bad = 0
    for line in eachline(md)
        m = match(ROW, line)
        m === nothing && continue
        rel, size, hash = m[1], parse(Int, m[2]), m[3]
        f = joinpath(ROOT, dir, rel)
        checked += 1
        if !isfile(f)
            println("NOFILE   $dir/$rel"); bad += 1
        elseif filesize(f) != size
            println("SIZE     $dir/$rel  recorded=$size actual=$(filesize(f))"); bad += 1
        else
            got = bytes2hex(open(sha256, f))
            got == hash || (println("HASH     $dir/$rel  actual=$got"); bad += 1)
        end
    end
    println(rpad(dir, 52), " recorded=", checked, bad == 0 ? "  OK" : "  $bad MISMATCH")
    return (checked, bad)
end

function main()
    total = 0; bad = 0
    for d in PROVENANCE_DIRS
        c, b = check_dir(d); total += c; bad += b
    end
    println("checked $total recorded files, $bad mismatches")
    return bad
end

exit(main() == 0 ? 0 : 1)

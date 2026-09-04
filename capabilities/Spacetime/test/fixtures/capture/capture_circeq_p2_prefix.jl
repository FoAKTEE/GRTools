# capture_circeq_p2_prefix.jl — regenerates test/fixtures/spacetime_circeq_p2_prefix_ref.jl.
#
# The fixture is a BITWISE REGRESSION PIN, not an external oracle: it records what
# Spacetime's own Schmidt closed forms returned BEFORE the equatorial p = 2 removable
# singularity was repaired, so the repair can be proved to have moved nothing else.
#
# Its oracle state is therefore the PRE-FIX tree and only that tree. Run it in a detached
# worktree at GRTools 5243f42 (the parent of the fix) to reproduce the fixture byte-identically:
#
#   WT=$(mktemp -d)
#   git -C <GRTools> worktree add --detach "$WT" 5243f42
#   cd "$WT" && JULIA_PKG_OFFLINE=true julia --project=capabilities/Spacetime -t1 \
#       capabilities/Spacetime/test/fixtures/capture/capture_circeq_p2_prefix.jl
#   git -C <GRTools> worktree remove --force "$WT"
#
# Run at HEAD it will NOT reproduce the file: the (a, p = 2, e = 0, |x| = 1) rows that the
# pre-fix code returned as NaN — and that this script therefore skips — are finite now and
# would be emitted as extra rows. That is the point of the pin, and why FIXTURES.toml records
# it as `tracked-notrerunnable`.
using Spacetime

const SPINS = (0.0, 0.5, -0.5, 0.9, -0.9, 0.998, -0.998)
const CIRC_EQ_KEYS = [(a, p, 0.0, x) for a in SPINS for p in (1.9, 2.5, 3.0, 6.0, 10.0) for x in (1.0, -1.0)]
const OTHER_BRANCH_KEYS = [
    (0.9, 10.0, 0.0, 0.5), (0.5, 7.0, 0.0, 0.3), (-0.7, 12.0, 0.4, 0.6),
    (0.9, 10.0, 0.1, 0.5), (0.998, 6.0, 0.2, 0.8), (0.0, 8.0, 0.3, 1.0),
    (0.9, 10.0, 0.0, 0.0), (0.9, 10.0, 0.2, 0.0), (0.998, 2.0, 0.0, 0.5),
    (0.998, 2.0, 0.2, 1.0), (0.9, 2.0, 0.0, 0.9),
]

# A key is baked only if the pre-fix code returned a fully finite triple there; the keys it
# could not evaluate (no circular geodesic => DomainError, or the p = 2 NaN itself) are the
# ones the fix is allowed to change and are deliberately absent.
function row(a, p, e, x)
    E  = geodesic_energy(a, p, e, x)
    Lz = geodesic_angular_momentum(a, p, e, x)
    Q  = carter_constant(a, p, e, x)
    (isfinite(E) && isfinite(Lz) && isfinite(Q)) || return nothing
    "    (a = $(repr(a)), p = $(repr(p)), e = $(repr(e)), x = $(repr(x)), " *
    "E = $(repr(E)), Lz = $(repr(Lz)), Q = $(repr(Q))),"
end
emit(keys) = filter(!isnothing, [try row(k...) catch; nothing end for k in keys])

function main()
    circ, other = emit(CIRC_EQ_KEYS), emit(OTHER_BRANCH_KEYS)
    out = joinpath(@__DIR__, "..", "spacetime_circeq_p2_prefix_ref.jl")
    open(out, "w") do io
        println(io, "# CITATION: the PRE-FIX Schmidt evaluation of Spacetime's own circular/generic closed forms")
        println(io, "# (geodesic_energy / geodesic_angular_momentum / carter_constant) captured once at GRTools")
        println(io, "# commit 5243f42 in a detached worktree, BEFORE the equatorial p = 2 removable singularity")
        println(io, "# was repaired. This is a BITWISE REGRESSION PIN, not an external oracle: its purpose is to")
        println(io, "# prove the repair moved no other (a, p, e, x). Full round-trip Float64 literals; compared")
        println(io, "# with `===` in test_geodesics_constants.jl.")
        println(io, "# Provenance: capture script test/fixtures/capture/capture_circeq_p2_prefix.jl (re-runnable")
        println(io, "# only at 5243f42 — see its header and the FIXTURES.toml row).")
        println(io, "const CIRCEQ_P2_PREFIX_PIN = [")
        println(io, "    # circular equatorial (e = 0, x = ±1), away from p = 2")
        foreach(r -> println(io, r), circ)
        println(io, "    # other branches (inclined circular, eccentric equatorial, polar), incl. three AT p = 2")
        println(io, "    # on branches the fix must NOT touch")
        foreach(r -> println(io, r), other)
        println(io, "]")
    end
    println("wrote $(length(circ)) circular-equatorial + $(length(other)) other-branch rows to $out")
end

main()

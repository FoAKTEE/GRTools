# OFFLINE GENERATOR + PROVENANCE for the baked QQNM source tables (src/QNM/qqnm_source_tables.jl).
#
# The owner's coefficient .tex files are themselves wolframscript-generated (derive-spin2-*.m under
# ref-code/QQNM-coeff/teukolsky-category/code/). This fixture is the ONE-SHOT offline bridge: it
# freezes the sha256 of each source .tex (the provenance ledger baked into the runtime module
# header) and — when the ref-code tree is present — RE-VERIFIES the current files hash-match, i.e.
# the baked Julia literals were transcribed from exactly these blobs. The runtime path
# (qqnm_source.jl) is PURE JULIA; wolframscript/tex-parsing lives only here.
#
# To regenerate after an owner table update: re-run the owner derive-*.m (wolframscript), re-hash the
# .tex, update QQNM_TEX_SHA256 here + the module header, and re-transcribe any changed literal.

using SHA

# The frozen ledger — sha256 of each owner source .tex (matches qqnm_source_tables.jl header).
const QQNM_TEX_SHA256 = Dict(
    "spin2-plus-plus-coefficients.generated.tex"        => "ea1668805fb4b0e5d6596bb86356b76f703012fc0f1b5745140d155f73e62301",
    "spin2-plus-plus-coefficients-full.generated.tex"   => "f26e3bea15a42b17a2353a5b769e4d12a4f0fa4ef37b663a6d691eb3aec880dc",
    "spin2-plus-minus-sixteen-basis.generated.tex"      => "ccffcda6d124029a39504fcc6da465455cafbe87eb1419324a3e587f9ae8c727",
    "spin2-plus-minus-four-basis.generated.tex"         => "bf562422665ba7473dabad991a8c02c04b45ab064d8486d681561d22c1d19c98",
    "spin2-minus-minus-zero-proof.generated.tex"        => "15356b9bbdffcd34280b5db844991767018d0d25c67e5203c011c0a56f6d424b",
    "spin2-minus-minus-org-cck-slot-proof.generated.tex" => "3830a58e73a6b91fdaa9249498733ebe805b2b50661bf7a8b8d96892a2bcaacd",
)

# Candidate roots for the owner coeff/ tree: ENV override first, then the walk-up
# (in-tree mirror layouts), then the recorded capture-site mirror. Absence is fine —
# the baked literals stay authoritative and the re-verification gate skips.
function _qqnm_coeff_dir()
    envroot = get(ENV, "GRTOOLS_QQNM_COEFF_ROOT", "")
    !isempty(envroot) && isdir(envroot) && return envroot
    here = @__DIR__
    for up in 1:8
        root = reduce((p,_)->dirname(p), 1:up; init=here)
        cand = joinpath(root, "ref-code", "QQNM-coeff", "teukolsky-category", "coeff")
        isdir(cand) && return cand
    end
    fallback = "/data/haiyangw/claude/emri-disk/ref-code/QQNM-coeff/teukolsky-category/coeff"
    isdir(fallback) && return fallback
    return nothing
end

"Return (checked::Int, missing::Int) after verifying each present .tex hashes to the frozen ledger."
function qqnm_verify_tex_provenance()
    dir = _qqnm_coeff_dir()
    dir === nothing && return (0, length(QQNM_TEX_SHA256))
    checked = 0; miss = 0
    for (fn, want) in QQNM_TEX_SHA256
        path = joinpath(dir, fn)
        if isfile(path)
            got = bytes2hex(sha256(read(path)))
            got == want || error("QQNM provenance MISMATCH for $fn:\n  want $want\n  got  $got")
            checked += 1
        else
            miss += 1
        end
    end
    (checked, miss)
end

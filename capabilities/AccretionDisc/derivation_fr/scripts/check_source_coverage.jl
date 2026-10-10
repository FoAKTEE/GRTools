# Usage: julia --project=.../derivation_fr scripts/check_source_coverage.jl FR24_SRC
# Read-only check against the source mirror; certificates retain source-qualified anchors.
using FRDerivation

length(ARGS) == 1 || error("supply the FR24 v1 src directory")
certs = FRDerivation.certify_fr24_dimensions()
all(c -> c.pass, certs) || error("dimensional certificates failed")
covered = Set(c.eqref for c in certs if occursin(".tex:", c.eqref))
expected = Set{String}()
pattern = r"\\begin\{(equation|align|gather)\*?\}.*?\\end\{\1\*?\}"s
for file in sort(readdir(joinpath(ARGS[1], "Sections"); join=true))
    endswith(file,".tex") || continue
    source = read(file,String)
    for m in eachmatch(pattern,source)
        line = 1+count(==('\n'),SubString(source,1,prevind(source,m.offset)))
        push!(expected, "$(basename(file)):$line")
    end
end
missing, stale = setdiff(expected,covered), setdiff(covered,expected)
isempty(missing) || error("missing source formula blocks: $missing")
isempty(stale) || error("stale formula anchors: $stale")
println("[PASS] FR24 source coverage: ", length(expected),
        " equation/align blocks, including all numbered formulas; ", length(certs), " dimensional certificates")

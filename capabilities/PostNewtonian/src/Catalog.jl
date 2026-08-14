"""
    PostNewtonian.Catalog

PN `.m` SeriesData discovery + substring browse layer over the stored series tree (stdlib `walkdir`
only — no PN-type dep). [`SERIES_DATA_DIR`](@ref) is the single source of truth every later loader
resolves series files against. Ports the upstream `FileNames`/`StringDelete`/`Fold`/
`SelectFromPNSeriesList` selector logic.
"""
module Catalog

export SERIES_DATA_DIR, pn_series_catalog, select_pn_series, pn_browse

"""
Absolute path to the bundled `SeriesData` directory (the 91-file MMA InputForm store). Catalog keys
are reported relative to this root, so it is the single source of truth every later loader/parser
resolves series files against. The store is package-OWNED (a real copy, no symlink into legacy);
the externalize-to-`Data.artifacts` move is deferred.
"""
const SERIES_DATA_DIR = normpath(joinpath(@__DIR__, "..", "SeriesData"))

"""
    pn_series_catalog(root=SERIES_DATA_DIR) -> Vector{String}

Discover all stored PN series under `root` as catalog keys: each `*.m` file path relative to `root`,
with a leading `/`, forward slashes, and the `.m` extension stripped (e.g.
`"/Schwarzschild/Circular/Flux/Infinity-l2m2"`). Upstream analogue:
`StringDelete[FileNames["*.m", dir, Infinity], {dir, ".m"}]`.
"""
function pn_series_catalog(root::AbstractString = SERIES_DATA_DIR)
    keys = String[]
    for (dir, _, files) in walkdir(root)
        for f in files
            endswith(f, ".m") || continue
            rel = relpath(joinpath(dir, f), root)
            key = "/" * replace(rel, '\\' => '/')
            push!(keys, chop(key; tail = 2))  # drop ".m"
        end
    end
    return keys
end

# Case-insensitive substring test, matching MMA `StringContainsQ[#, needle, IgnoreCase -> True]`.
_contains_ci(needle::AbstractString, hay::AbstractString) =
    occursin(lowercase(needle), lowercase(hay))

"""
    select_pn_series(list, pattern) -> Vector{String}

One filter step of the PN series browser. An empty `pattern` returns `list` unchanged. A pattern
beginning with `!` keeps entries that do NOT contain the remainder (case-insensitive); otherwise
keeps entries that DO contain `pattern` (case-insensitive). Upstream analogue:
`SelectFromPNSeriesList`.
"""
function select_pn_series(list, pattern::AbstractString)
    isempty(pattern) && return collect(list)
    if startswith(pattern, '!')
        needle = pattern[2:end]
        return filter(s -> !_contains_ci(needle, s), list)
    end
    return filter(s -> _contains_ci(pattern, s), list)
end

"""
    pn_browse(query, root=SERIES_DATA_DIR) -> Vector{String}

Apply a sequence of filter patterns to the full catalog and return the matching series names. `query`
is a single pattern string or an iterable of pattern strings, applied left to right. Upstream
analogue: `Fold[SelectFromPNSeriesList, allSeries, Flatten[{quantity},1]]`.
"""
function pn_browse(query, root::AbstractString = SERIES_DATA_DIR)
    patterns = query isa AbstractString ? (query,) : query
    return foldl(select_pn_series, patterns; init = pn_series_catalog(root))
end

end # module Catalog

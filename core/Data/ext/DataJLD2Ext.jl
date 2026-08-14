# JLD2 reader-writer BODIES for the Data spine (design §5.4 weakdep extension). Activates only when
# the loading package (e.g. Hybrid) declares JLD2 in its own [deps]. Provides the WaSABI `.data`
# parse-once-cache path: `wl_to_jld2` (parse + store) and `read_jld2` (cache load).
module DataJLD2Ext

using Data
using Data: read_wl_data
using JLD2

"""
    Data.wl_to_jld2(src, dst) -> dst

Parse a WL `.data` catalog with `Data.read_wl_data` (BigFloat-exact) and store the records under key
`"data"` in a JLD2 file at `dst`. Rewritten over the shared `read_wl_data` reader so the cache is a
bit-for-bit snapshot of the parse (de-dups WaSABI/wl_parser.jl).
"""
function Data.wl_to_jld2(src::AbstractString, dst::AbstractString)
    data = read_wl_data(src)
    JLD2.jldsave(dst; data = data)
    return dst
end

"""
    Data.read_jld2(path; key="data") -> value

Load `path[key]` from a JLD2 file (default key `"data"`). The cache-load sibling of
[`Data.wl_to_jld2`](@ref): a catalog is parsed once then loaded, not re-parsed.
"""
function Data.read_jld2(path::AbstractString; key::AbstractString = "data")
    return JLD2.load(path, key)
end

end # module DataJLD2Ext

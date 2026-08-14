# JLD2 ingest/cache path for parsed catalogs (the WaSABI `.data` round-trip). A WL `.data` catalog is
# parsed once with `read_wl_data` (BigFloat-exact) and cached as JLD2, then LOADED — not re-parsed —
# on subsequent runs. De-dups WaSABI/wl_parser.jl `wl_to_jld2`, rewritten over `Data.read_wl_data`.
# Bodies live in ext/DataJLD2Ext.jl; the stubs here error until JLD2 is loaded by the consumer.

const _NEED_JLD2 = "load JLD2 (`import JLD2`) — Data's JLD2 reader-writer lives in a weakdep extension"

"""
    wl_to_jld2(src, dst) -> dst

Convert a WL `.data` catalog to a JLD2 cache: parse `src` once with [`read_wl_data`](@ref) (full
`BigFloat` precision) and store the records under key `"data"` at `dst`. Parse-once-load-later, so the
expensive InputForm parse is amortized and the cache round-trips bit-for-bit. Method by `DataJLD2Ext`.
"""
function wl_to_jld2 end

"""
    read_jld2(path; key="data") -> value

Load a value (default key `"data"`) from a JLD2 cache produced by [`wl_to_jld2`](@ref) (or any JLD2
file). The sibling load entry so a catalog is parsed once then loaded, not re-parsed. Method by
`DataJLD2Ext`.
"""
function read_jld2 end

# Untyped `(Any,...)` fallbacks: the JLD2 ext supplies strictly-more-specific `(AbstractString,...)`
# methods, so loading JLD2 refines dispatch WITHOUT method-overwrite (illegal during ext
# precompilation). Without JLD2 these fire and point at the missing trigger.
wl_to_jld2(src, dst) = error("wl_to_jld2: " * _NEED_JLD2)
read_jld2(path; kwargs...) = error("read_jld2: " * _NEED_JLD2)

export wl_to_jld2, read_jld2

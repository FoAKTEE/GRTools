# Injected data-directory policy. Replaces the legacy hard-coded `ref-code/.../run_r8/fields`
# (DiscEMRIForce loader.jl:17 H1DIR) and the FEW `.fewdata` `@__DIR__` `..` escapes
# (ampinterp2d.jl SCHWARZ_AMP_FILE, trajectory/flux.jl, amplitude/romannet.jl) with a process-global
# root resolved explicitly or from `GRTOOLS_DATADIR` — never a `..` walk out of a package directory.

const _DATADIR = Ref{Union{Nothing,String}}(nothing)

"""
    set_datadir!(path) -> String

Set the process-global GrTools data root (where on-disk corpora — h1 fields, FEW grids, PN/WaSABI
catalogs — are resolved from). Returns the normalized absolute path. Pass `nothing` to clear the
override and fall back to [`default_datadir`](@ref).
"""
function set_datadir!(path::Union{Nothing,AbstractString})
    if path === nothing
        _DATADIR[] = nothing
        return ""
    end
    isdir(path) || throw(ArgumentError("set_datadir!: not a directory: $path"))
    _DATADIR[] = abspath(path)
    return _DATADIR[]
end

"""
    default_datadir() -> String

Resolve the data root when no explicit override is set: the `GRTOOLS_DATADIR` env var if present,
else (future) a pinned `Artifacts` root. Throws if nothing is configured — there is deliberately NO
implicit `ref-code` fallback, so a missing configuration fails loudly instead of escaping a package.
"""
function default_datadir()
    env = get(ENV, "GRTOOLS_DATADIR", "")
    isempty(env) || return abspath(env)
    # Artifacts hook: a registered "GrTools_data" artifact resolves here once pinned (none yet).
    throw(ErrorException(
        "Data: no data root configured — call set_datadir!(path) or set ENV[\"GRTOOLS_DATADIR\"]"))
end

"""
    datadir() -> String

The active GrTools data root: the [`set_datadir!`](@ref) override if set, else [`default_datadir`](@ref).
"""
datadir() = _DATADIR[] === nothing ? default_datadir() : _DATADIR[]

"""
    datapath(parts...) -> String

`joinpath` under the active [`datadir`](@ref). The canonical way for every reader to resolve an
on-disk asset path.
"""
datapath(parts::AbstractString...) = joinpath(datadir(), parts...)

# `h1-l{l}m{m}.h5` mode-file name contract (DiscEMRIForce loader.jl / SecondOrderRicci read_h1.jl).
const _H1_MODE_RE = r"^h1-l(\d+)m(\d+)\.h5$"

"""
    available_modes(; dir=datadir()) -> Vector{NTuple{2,Int}}

Discover the `(l, m)` h1 mode files in `dir` (names matching `h1-l{l}m{m}.h5`), returned sorted.
De-dups `DiscEMRIForce.available_modes`; the path is datadir-injected (no `H1DIR` escape).
"""
function available_modes(; dir::AbstractString = datadir())
    out = NTuple{2,Int}[]
    for b in readdir(dir)
        mch = match(_H1_MODE_RE, b)
        mch === nothing && continue
        push!(out, (parse(Int, mch.captures[1]), parse(Int, mch.captures[2])))
    end
    return sort(out)
end

export set_datadir!, default_datadir, datadir, datapath, available_modes

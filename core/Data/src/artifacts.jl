# Artifacts resolution + offline/HPC pre-staging (design §8.2, §10). A GrTool data corpus `name`
# resolves through an explicit OVERRIDE -> ARTIFACT -> ERROR chain: the injected `datadir`
# (`set_datadir!` / GRTOOLS_DATADIR) wins FIRST (precedence — HPC/offline + tests), else a declared +
# provisioned `Artifacts.toml` entry resolves lazily via its `Artifacts` stdlib tree-hash (the runtime
# `@artifact_str` path: `artifact_hash` -> `artifact_exists` -> `artifact_path(hash)`), else a loud
# throw via `datapath` — never a `..` escape out of a package directory.

# The shipped Artifacts.toml (carries a real bound `kerr_ecc_eq_flux` entry, provisioned locally from
# the cached `.fewdata` corpus). `_ARTIFACTS_TOML[]` is host/test-redirectable via `set_artifacts_toml!`
# so a resolver can be pointed at a freshly-bound temporary toml (round-trip grounding).
const _DEFAULT_ARTIFACTS_TOML = normpath(joinpath(@__DIR__, "..", "Artifacts.toml"))
const _ARTIFACTS_TOML = Ref{String}(_DEFAULT_ARTIFACTS_TOML)

"""
    set_artifacts_toml!(path) -> String

Point the artifact resolver at `path`; pass `nothing` to restore the package's shipped `Artifacts.toml`.
Used to resolve a temporary, freshly-`bind_artifact!`ed toml in tests and for host-injected corpora.
"""
function set_artifacts_toml!(path::Union{Nothing,AbstractString})
    _ARTIFACTS_TOML[] = path === nothing ? _DEFAULT_ARTIFACTS_TOML : String(path)
    return _ARTIFACTS_TOML[]
end

"""
    artifacts_toml() -> String

The active `Artifacts.toml` the resolver consults (the shipped one unless [`set_artifacts_toml!`](@ref)
redirected it).
"""
artifacts_toml() = _ARTIFACTS_TOML[]

# Is `name` a declared, resolvable artifact in the active Artifacts.toml? (Backed by a real bound
# `kerr_ecc_eq_flux` entry; arbitrary names still fall through to the injected datadir root.)
function _artifact_declared(name::AbstractString)
    toml = _ARTIFACTS_TOML[]
    isfile(toml) || return false
    return Artifacts.artifact_meta(name, toml) !== nothing
end

# The configured data-root override (`set_datadir!` / GRTOOLS_DATADIR), or `nothing` if none — a
# non-throwing peek at [`datadir`](@ref) so the resolver can give the override precedence without raising.
function _datadir_override()
    _DATADIR[] === nothing || return _DATADIR[]
    env = get(ENV, "GRTOOLS_DATADIR", "")
    return isempty(env) ? nothing : abspath(env)
end

# override -> artifact -> error resolution chain (the explicit precedence ladder, design §10).
function _resolve_artifact(name::AbstractString)
    # 1. Injected datadir override (`set_datadir!` / GRTOOLS_DATADIR) wins — HPC/offline + test redirect.
    root = _datadir_override()
    root === nothing || return joinpath(root, name)
    # 2. A declared + provisioned artifact resolves lazily via its tree-hash (runtime @artifact_str path).
    if _artifact_declared(name)
        h = Artifacts.artifact_hash(name, _ARTIFACTS_TOML[])
        if h !== nothing && Artifacts.artifact_exists(h)
            return Artifacts.artifact_path(h)
        end
        # Declared-but-not-provisioned (lazy, undownloaded) ⇒ fall through to the loud error below.
    end
    # 3. No override, no provisioned artifact ⇒ `datapath` throws loudly (no implicit ref-code escape).
    return datapath(name)
end

"""
    artifact_path(name) -> String

Resolve the on-disk root of a GrTool data corpus `name` through the override -> artifact -> error
chain (`_resolve_artifact`): the injected [`datadir`](@ref) override (`set_datadir!` /
`GRTOOLS_DATADIR`) FIRST (precedence), else a declared + provisioned `Artifacts.toml` entry resolved
lazily via the `Artifacts` stdlib tree-hash, else a loud throw (no `..` escape).
"""
artifact_path(name::AbstractString) = _resolve_artifact(name)

"""
    stage_datadir!(src, dst) -> String

Offline/HPC pre-staging hook (§8.2): copy a corpus `src` (a directory's contents, or a single file)
into a local data root `dst`, then `set_datadir!(dst)` so subsequent reads resolve there — for sites
without artifact download (e.g. the 4.74GB ZNAmps grid staged from shared scratch). Returns the
activated root.
"""
function stage_datadir!(src::AbstractString, dst::AbstractString)
    ispath(src) || throw(ArgumentError("stage_datadir!: source does not exist: $src"))
    mkpath(dst)
    if isdir(src)
        for entry in readdir(src)
            cp(joinpath(src, entry), joinpath(dst, entry); force = true)
        end
    else
        cp(src, joinpath(dst, basename(src)); force = true)
    end
    return set_datadir!(dst)
end

"""
    ensure_local(name) -> String

Resolve corpus `name` to a present on-disk path under the active data root ([`artifact_path`]), or
throw loudly with the staging hint. The site-local guarantee for the large local-only corpora
(`stage_datadir!` them first, or set `GRTOOLS_DATADIR`) — never a silent download or `ref-code` escape.
"""
function ensure_local(name::AbstractString)
    p = artifact_path(name)
    ispath(p) || throw(ErrorException(
        "ensure_local: corpus '$name' not found at $p — pre-stage it with " *
        "stage_datadir!(src, root) or set ENV[\"GRTOOLS_DATADIR\"]"))
    return p
end

export artifact_path, stage_datadir!, ensure_local, set_artifacts_toml!, artifacts_toml

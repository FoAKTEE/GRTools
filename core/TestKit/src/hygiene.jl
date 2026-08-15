# §7.7 binding CI hygiene gates (design GRTOOLS_DESIGN_CAPABILITY.md §7.7, gates 1/2/5) — the
# executable dep-budget + closure-load checkers shared by every GrTools capability/core package.
#
# Two checkers:
#   * `assert_dep_budget` — gate 1 (dep-budget) + gate 5 (extension-trigger, via `must_declare`)
#     as a TOML-STATIC check of the package's OWN Project.toml ([deps]/[weakdeps]/[sources]).
#     Static by design: it fails on a *declared* poison edge even before any code loads it.
#   * `assert_not_loaded` — gate 2 (closure-load): asserts the RUNTIME dependency closure did not
#     pull a poison module into `Base.loaded_modules`.
#
# Dependency: TOML (stdlib, zero cost) + Test (already a hard TestKit dep — this IS the harness).

import TOML
using Test: Test

# ---------------------------------------------------------------------------------------------------
# Name matching: exact, or trailing-`*` prefix glob ("OrdinaryDiffEq*" ⊃ OrdinaryDiffEqVerner/…Core)
# ---------------------------------------------------------------------------------------------------

"""
    _name_matches(pattern, name) -> Bool

`true` iff `name` matches `pattern`. A `pattern` ending in `*` matches by prefix (so the design's
"OrdinaryDiffEq*" poison axis catches `OrdinaryDiffEqVerner`, `OrdinaryDiffEqCore`, …); any other
pattern must match exactly.
"""
_name_matches(pattern::AbstractString, name::AbstractString) =
    endswith(pattern, "*") ? startswith(name, chop(pattern)) : pattern == name

# ---------------------------------------------------------------------------------------------------
# Gate 1 + 5: TOML-static dep budget
# ---------------------------------------------------------------------------------------------------

"""
    assert_dep_budget(project_toml_path;
                      allow=nothing, forbid=String[], must_declare=String[],
                      weakdeps_allow=nothing)

§7.7 gate 1 (dep-budget) + gate 5 (extension-trigger) as a TOML-static `@test` block over the
package's OWN `Project.toml`. Parses `[deps]`, `[weakdeps]` and `[sources]` (TOML stdlib; nothing
is loaded) and asserts:

  1. **allow** (if given): every `[deps]` name ∈ `allow` — the sanctioned closed budget. Names in
     `allow` are exact (no glob). NOTE this also enforces "X forbidden in hard deps but sanctioned
     as a weakdep" (e.g. Data's HDF5/JLD2): leave X out of `allow`, put it in `weakdeps_allow`.
  2. **forbid**: no name matching a `forbid` pattern appears in `[deps]` OR `[weakdeps]`.
     Patterns support a trailing-`*` prefix glob (`"OrdinaryDiffEq*"`). `[extras]` is deliberately
     NOT scanned — test-only oracles (e.g. Wigner's WignerSymbols) are sanctioned there.
  3. **must_declare**: every name ∈ `[deps]` — the gate-5 extension triggers (an ext fires only
     when the LOADING package declares the trigger in its own `[deps]`, design §7.2).
  4. **weakdeps_allow** (if given): every `[weakdeps]` name ∈ `weakdeps_allow`.
  5. **always**: no `[sources]` entry path contains `"packages/"` or `"ref-code/"` — the legacy
     trees are dead as load-path roots (Phase 7-B invariant).

Runs inside its own labelled `Test.@testset` so the gate is visible in the suite output. Returns
`nothing`.
"""
function assert_dep_budget(project_toml_path::AbstractString;
                           allow::Union{Nothing,Vector{String}}=nothing,
                           forbid::Vector{String}=String[],
                           must_declare::Vector{String}=String[],
                           weakdeps_allow::Union{Nothing,Vector{String}}=nothing)
    proj = TOML.parsefile(project_toml_path)
    pkgname  = get(proj, "name", "?")
    deps     = sort!(collect(keys(get(proj, "deps",     Dict{String,Any}()))))
    weakdeps = sort!(collect(keys(get(proj, "weakdeps", Dict{String,Any}()))))
    sources  = get(proj, "sources", Dict{String,Any}())

    Test.@testset "§7.7 gates 1+5: dep budget ($pkgname)" begin
        # (1) closed allow-budget over hard [deps]
        if allow !== nothing
            unsanctioned_deps = setdiff(deps, allow)
            Test.@test unsanctioned_deps == String[]
        end
        # (2) poison axes: forbidden anywhere in [deps]/[weakdeps]
        forbidden_hits = ["$pat => $name" for pat in forbid
                          for name in vcat(deps, weakdeps) if _name_matches(pat, name)]
        Test.@test forbidden_hits == String[]
        # (3) gate-5 extension triggers: must be DIRECT hard [deps]
        missing_triggers = [name for name in must_declare if !(name in deps)]
        Test.@test missing_triggers == String[]
        # (4) closed allow-budget over [weakdeps]
        if weakdeps_allow !== nothing
            unsanctioned_weakdeps = setdiff(weakdeps, weakdeps_allow)
            Test.@test unsanctioned_weakdeps == String[]
        end
        # (5) no [sources] path may point into the legacy packages/* or ref-code/* trees
        legacy_sources = ["$name => $(spec["path"])" for (name, spec) in sources
                          if spec isa AbstractDict && haskey(spec, "path") &&
                             (occursin("packages/", spec["path"]) || occursin("ref-code/", spec["path"]))]
        Test.@test legacy_sources == String[]
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------------
# Gate 2: runtime closure-load poison check
# ---------------------------------------------------------------------------------------------------

"""
    assert_not_loaded(mods::Vector{String})

§7.7 gate 2 (closure-load): asserts that NO module whose name matches an entry of `mods` (exact, or
trailing-`*` prefix glob) is present in `Base.loaded_modules`.

CAVEAT (binding usage contract): run this AFTER `using <Pkg>` — it asserts the package's runtime
dependency closure did not pull a poison module in. The test env's own extras are excluded BY
CONSTRUCTION: the caller lists only poison names that nothing in the test env legitimately loads
(e.g. Wigner's test env legitimately loads its WignerSymbols test oracle, so Wigner's list must not
contain it — its runtime elimination is asserted TOML-statically by `assert_dep_budget` instead).

Runs inside its own labelled `Test.@testset`. Returns `nothing`.
"""
function assert_not_loaded(mods::Vector{String})
    loaded = Set(string(m) for m in values(Base.loaded_modules))
    Test.@testset "§7.7 gate 2: closure-load poison" begin
        poison_loaded = ["$pat => $name" for pat in mods
                         for name in loaded if _name_matches(pat, name)]
        Test.@test poison_loaded == String[]
    end
    return nothing
end

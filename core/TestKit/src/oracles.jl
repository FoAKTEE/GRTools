# Oracle-availability gate (repro-audit G3/G4): make every "external oracle absent" branch VISIBLE.
#
# Before this module the packages guarded their mirror-backed gates with a bare `if isfile(path) …
# else @info … end`, so a suite ran green whether or not the oracle was ever read. Three rules now:
#   1. mirror roots are resolved through a NAMED `GRTOOLS_*` environment variable with the capture-site
#      path as the default (`oracle_root`), so another site redirects them without editing tests;
#   2. a missing oracle records a `Broken` result (`oracle_or_skip`), so the Test summary's `Broken`
#      column counts exactly how many gates did not run;
#   3. `ENV["GRTOOLS_STRICT_ORACLES"] = "1"` promotes every such skip to a hard `@test` failure — the
#      CI/reproduction mode in which the corpus is REQUIRED to be present.

"""
    strict_oracles() -> Bool

`true` when `ENV["GRTOOLS_STRICT_ORACLES"]` is set to one of `1` / `true` / `yes` (case-insensitive).
In strict mode [`oracle_or_skip`](@ref) turns a missing oracle into a test FAILURE instead of a
recorded `Broken`. Read at call time (not at precompile), so a test file can flip it locally.
"""
function strict_oracles()
    v = lowercase(strip(get(ENV, "GRTOOLS_STRICT_ORACLES", "")))
    return v in ("1", "true", "yes", "on")
end

"""
    oracle_root(envvar::AbstractString, default::AbstractString) -> String

Resolve the filesystem root of an out-of-tree oracle corpus: `ENV[envvar]` when set, otherwise the
capture-site `default`. This is the single sanctioned way for a test to reach a mirror outside the
package tree — relative `..` walk-ups are forbidden because they silently resolve to a non-existent
path once the package is checked out anywhere else (repro-audit §a.2).

Every `envvar` used in this repository is documented in `docs/src/topics/data-artifacts.md`.
"""
oracle_root(envvar::AbstractString, default::AbstractString) = String(get(ENV, envvar, default))

"""
    oracle_or_skip(available::Bool, what::AbstractString; env=nothing) -> Bool

Gate an oracle-backed testset. Returns `available`. When the oracle is ABSENT:

* default mode — records a `Broken` (`:skipped`) result carrying `what`, so the missing gate shows up
  in the `Broken` column of the Test summary instead of vanishing into an `@info`;
* strict mode ([`strict_oracles`](@ref)) — records a `@test` FAILURE naming the same message.

`env` names the environment variable that would supply the corpus; it is appended to the message as a
hint. Use as the condition of the `if` that wraps the oracle-backed asserts:

```julia
if oracle_or_skip(isfile(f), "read_ampgrid vs raw HDF5"; env = "GRTOOLS_FEWDATA_DIR")
    ...
end
```
"""
function oracle_or_skip(available::Bool, what::AbstractString; env = nothing)
    available && return true
    hint = env === nothing ? "" : " (set ENV[\"$env\"])"
    msg = string(what, " — oracle absent", hint)
    if strict_oracles()
        # String comparison so the Test failure report prints the full message on its "Evaluated" line.
        Test.@test msg == "no missing oracle (GRTOOLS_STRICT_ORACLES=1)"
    else
        Test.record(Test.get_testset(), Test.Broken(:skipped, msg))
    end
    return false
end

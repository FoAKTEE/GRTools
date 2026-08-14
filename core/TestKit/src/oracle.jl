# Leaf tolerance + oracle harness for the GrTool capability packages (design §5.6 / DAG §GrToolTestKit).
#
# This is the dependency LEAF of TestKit: it depends on nothing but `Test` (+ stdlib `Printf` for one
# error message), so it cold-compiles instantly. It owns three pieces:
#   * `@approx_oracle`  — an `isapprox` assert that ALSO RECORDS the (verbatim-frozen) tolerance, so a
#                         CI pass can prove every new oracle's tol was registered and never silently
#                         loosened. The macro RETURNS the boolean; wrap it in `@test` to assert.
#   * tolerance CLASSES — `:machine`/`:spectral`/`:physics`/`:gpu` named (rtol, atol) bundles for NEW
#                         tests, ordered tightest→loosest.
#   * `assert_can_only_tighten` — CI guard enforcing `class_tol <= original_tol` (tolerances may only
#                         TIGHTEN over a relocation, never loosen).

import Printf
using Test: Test

# ---------------------------------------------------------------------------------------------------
# Tolerance recording
# ---------------------------------------------------------------------------------------------------

"""
    OracleRecord

One recorded `@approx_oracle` evaluation: the `label`, the `(rtol, atol)` actually used, and whether
the `isapprox` held. The recording is what lets CI prove a new oracle's tolerance was frozen and not
loosened over time.
"""
struct OracleRecord
    label::String
    rtol::Float64
    atol::Float64
    passed::Bool
end

const _ORACLE_LOG = OracleRecord[]

"""
    oracle_records() -> Vector{OracleRecord}

The append-only log of every [`@approx_oracle`](@ref) evaluation since the last [`reset_oracle_log!`](@ref).
"""
oracle_records() = _ORACLE_LOG

"""
    reset_oracle_log!() -> Vector{OracleRecord}

Clear the [`oracle_records`](@ref) log (call at the top of a fresh test run).
"""
reset_oracle_log!() = (empty!(_ORACLE_LOG); _ORACLE_LOG)

"""
    @approx_oracle(value, ref; rtol=…, atol=…, label="…")

Evaluate `isapprox(value, ref; rtol, atol)`, RECORD the tolerance used in [`oracle_records`](@ref),
and RETURN the boolean. Wrap in `@test` to assert: `@test @approx_oracle(got, exp; rtol=1e-12)`.

The recorded tolerance is the verbatim-frozen number — the harness never re-derives it, so a CI pass
can prove no new oracle silently loosened its tol. `rtol`/`atol` default to Base's `isapprox` defaults
(rtol = √eps of the promoted type, atol = 0); pass them explicitly to freeze a specific tolerance.
"""
macro approx_oracle(args...)
    rtol_ex = nothing
    atol_ex = nothing
    label_ex = nothing
    positional = Any[]
    # Flatten `;`-grouped kwargs (a `:parameters` Expr) and inline `name=value` kwargs alike.
    function _take_kw(name, val)
        if name === :rtol
            rtol_ex = val
        elseif name === :atol
            atol_ex = val
        elseif name === :label
            label_ex = val
        else
            throw(ArgumentError("@approx_oracle: unknown keyword `$name`"))
        end
    end
    for a in args
        if a isa Expr && a.head === :parameters
            for kw in a.args
                (kw isa Expr && kw.head === :kw) ||
                    throw(ArgumentError("@approx_oracle: bad keyword in `;` group: $kw"))
                _take_kw(kw.args[1], kw.args[2])
            end
        elseif a isa Expr && (a.head === :(=) || a.head === :kw)
            _take_kw(a.args[1], a.args[2])
        else
            push!(positional, a)
        end
    end
    length(positional) == 2 ||
        throw(ArgumentError("@approx_oracle expects (value, ref; rtol, atol, label), " *
                            "got $(length(positional)) positional args"))
    value, ref = positional[1], positional[2]
    label = label_ex === nothing ? string(value, " ≈ ", ref) : label_ex
    # Build the isapprox call, passing only the kwargs the caller supplied so Base's defaults survive.
    isapprox_kws = Expr[]
    rtol_ex !== nothing && push!(isapprox_kws, Expr(:kw, :rtol, esc(rtol_ex)))
    atol_ex !== nothing && push!(isapprox_kws, Expr(:kw, :atol, esc(atol_ex)))
    return quote
        local _v = $(esc(value))
        local _r = $(esc(ref))
        local _rt = $(rtol_ex === nothing ? :(_default_rtol(_v, _r)) : esc(rtol_ex))
        local _at = $(atol_ex === nothing ? 0.0 : esc(atol_ex))
        local _ok = isapprox(_v, _r, $(isapprox_kws...))
        push!(_ORACLE_LOG, OracleRecord($(esc(label)), Float64(_rt), Float64(_at), _ok))
        _ok
    end
end

# Base's isapprox default rtol = √eps(promoted real type) when atol==0; used only for the LOG entry
# when the caller did not freeze rtol (the actual comparison still uses Base's own default).
_default_rtol(v, r) = sqrt(eps(float(real(_promote_eltype(v, r)))))
_promote_eltype(v, r) = promote_type(_eltype(v), _eltype(r))
_eltype(x::Number) = typeof(x)
_eltype(x) = eltype(x)

# ---------------------------------------------------------------------------------------------------
# Tolerance classes (NEW tests only) + the can-only-tighten guard
# ---------------------------------------------------------------------------------------------------

"""
    TOLERANCE_CLASSES :: Dict{Symbol,NamedTuple}

Named `(rtol, atol)` bundles for NEW tests, ordered tightest→loosest:

  * `:machine`  — bit-level / closed-form identities (rtol 1e-14).
  * `:spectral` — spectral/series-convergent numerics (rtol 1e-10).
  * `:physics`  — physics-level agreement, accumulated error (rtol 1e-6).
  * `:gpu`      — CPU-vs-GPU `sincos`/`fma` divergence (rtol 1e-4) — the only class wider than physics.

These are defaults for greenfield tests; a relocation must instead freeze the LEGACY tolerance verbatim
and may only TIGHTEN to a class via [`assert_can_only_tighten`](@ref).
"""
const TOLERANCE_CLASSES = Dict{Symbol,NamedTuple{(:rtol, :atol),Tuple{Float64,Float64}}}(
    :machine  => (rtol = 1e-14, atol = 0.0),
    :spectral => (rtol = 1e-10, atol = 0.0),
    :physics  => (rtol = 1e-6,  atol = 0.0),
    :gpu      => (rtol = 1e-4,  atol = 0.0),
)

"""
    tolerance_class(class::Symbol) -> NamedTuple{(:rtol,:atol)}

The default `(rtol, atol)` for a NEW-test tolerance `class` (`:machine`/`:spectral`/`:physics`/`:gpu`).
Throws on an unknown class.
"""
function tolerance_class(class::Symbol)
    haskey(TOLERANCE_CLASSES, class) ||
        throw(ArgumentError("unknown tolerance class $(repr(class)); valid: " *
                            string(sort(collect(keys(TOLERANCE_CLASSES))))))
    return TOLERANCE_CLASSES[class]
end

"""
    assert_can_only_tighten(class_tol, original_tol) -> true

CI guard: a NEW test's tolerance may only TIGHTEN the legacy/original one. Asserts
`class_tol <= original_tol` and returns `true`; throws an `AssertionError` (naming both values) if the
candidate tolerance is LOOSER. `class_tol` may be a `Symbol` (resolved via [`tolerance_class`]'s rtol)
or a number.
"""
function assert_can_only_tighten(class_tol, original_tol::Real)
    ct = class_tol isa Symbol ? tolerance_class(class_tol).rtol : float(class_tol)
    ct <= original_tol ||
        throw(AssertionError(Printf.@sprintf(
            "tolerance may only tighten: class_tol=%.3e exceeds original_tol=%.3e", ct, original_tol)))
    return true
end

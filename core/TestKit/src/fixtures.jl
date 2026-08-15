# Data-backed fixture loader + characterization-gate harness (design §5.6 / DAG §GrToolsTestKit).
#
# Built second: depends on `oracle.jl` (the tolerance/recording layer) and on `Data` (fixture IO). It
# owns two pieces:
#   * `load_fixture`    — captured-fixture loader, riding `Data.datapath` + `Data.read_dat`.
#   * characterization gate — snapshot a reference function's full-precision output on a set of inputs,
#                         then assert a candidate reproduces it to the recorded ULP (a perturbed
#                         candidate must FAIL). This is the additive-then-sweep relocation gate.

using Data: read_dat, datapath

# ---------------------------------------------------------------------------------------------------
# Fixture loader (rides the Data spine)
# ---------------------------------------------------------------------------------------------------

"""
    load_fixture(parts...; kwargs...) -> Matrix

Load a captured characterization fixture as a numeric matrix, resolving the path under the active
GrTools data root (`Data.datapath`) and reading it via `Data.read_dat`. `kwargs` pass through to the
underlying `readdlm` (e.g. `skipstart`). Configure the root first with `Data.set_datadir!(dir)` (or
`ENV["GRTOOLS_DATADIR"]`); there is deliberately no implicit fallback.
"""
load_fixture(parts...; kwargs...) = read_dat(datapath(parts...); kwargs...)

# ---------------------------------------------------------------------------------------------------
# Characterization gate
# ---------------------------------------------------------------------------------------------------

"""
    CharacterizationSnapshot

A frozen record of a reference function's full-precision outputs (`outputs`) on a fixed input set
(`inputs`), with the element type preserved (BigFloat snapshots stay BigFloat). Feed to
[`characterization_gate`](@ref) to prove a candidate reproduces it to the recorded ULP.
"""
struct CharacterizationSnapshot{I,O}
    inputs::Vector{I}
    outputs::Vector{O}
end

"""
    characterization_snapshot(f, inputs) -> CharacterizationSnapshot

Snapshot reference function `f` at full precision over `inputs` (one output per input). The outputs
are stored exactly as `f` returns them — no rounding — so a later candidate is gated against the true
reference, never against its own emission.
"""
function characterization_snapshot(f, inputs)
    iv = collect(inputs)
    ov = [f(x) for x in iv]
    return CharacterizationSnapshot(iv, ov)
end

"""
    characterization_gate(candidate, snap::CharacterizationSnapshot; ulp=0, rtol=0.0, atol=0.0) -> Bool

Return `true` iff `candidate` reproduces every snapshot output. Each pair must satisfy EITHER an
`isapprox(got, ref; rtol, atol)` (when `rtol`/`atol` given) OR an exact / `ulp`-bounded equality
(default `ulp=0` ⇒ bit-exact for floats, `==` otherwise). A candidate that perturbs any output beyond
the bound FAILS. This is the additive-then-sweep relocation gate: snapshot the legacy primitive, then
gate the relocated one before deleting the fork.
"""
function characterization_gate(candidate, snap::CharacterizationSnapshot;
                               ulp::Integer = 0, rtol::Real = 0.0, atol::Real = 0.0)
    ulp >= 0 || throw(ArgumentError("ulp budget must be nonnegative (got $ulp)"))
    for (x, ref) in zip(snap.inputs, snap.outputs)
        got = candidate(x)
        _gate_match(got, ref; ulp, rtol, atol) || return false
    end
    return true
end

# Per-output match: isapprox when a tolerance is requested, else exact/ULP equality.
function _gate_match(got, ref; ulp::Integer, rtol::Real, atol::Real)
    if rtol != 0 || atol != 0
        return isapprox(got, ref; rtol, atol)
    end
    return _ulp_equal(got, ref, ulp)
end

_ulp_equal(a, b, ulp::Integer) = a == b
# BigFloat (and any non-IEEE float) has no fixed bit layout: gate bit-exact at ulp=0, else within an
# `ulp`-scaled local `eps`. Must precede the IEEE method below (BigFloat <: AbstractFloat).
function _ulp_equal(a::T, b::T, ulp::Integer) where {T<:Union{BigFloat}}
    (isnan(a) && isnan(b)) && return true
    ulp == 0 && return a == b
    return abs(a - b) <= ulp * eps(max(abs(a), abs(b), one(T)))
end
const _IEEEFloat = Union{Float16,Float32,Float64}

function _ulp_equal(a::T, b::T, ulp::Integer) where {T<:_IEEEFloat}
    (isnan(a) && isnan(b)) && return true
    ulp == 0 && return a == b
    # Count representable floats between a and b (signbit-aware), compare to the ULP budget.
    a == b && return true
    (isfinite(a) && isfinite(b)) || return false
    return _ulp_within(_ulp_distance(a, b), ulp)
end
function _ulp_distance(a::T, b::T) where {T<:_IEEEFloat}
    # Monotonic unsigned bit ordering. The sign/magnitude fold maps both signed zeros to the
    # same key, avoids signed subtraction overflow across zero, and follows the IEEE width of T.
    ua = _ulp_key(a)
    ub = _ulp_key(b)
    return ua >= ub ? ua - ub : ub - ua
end
@inline _ulp_uint_type(::Type{Float16}) = UInt16
@inline _ulp_uint_type(::Type{Float32}) = UInt32
@inline _ulp_uint_type(::Type{Float64}) = UInt64
function _ulp_key(x::T) where {T<:_IEEEFloat}
    U = _ulp_uint_type(T)
    bits = reinterpret(U, x)
    signmask = one(U) << (8 * sizeof(U) - 1)
    magnitude = bits & (signmask - one(U))
    return iszero(bits & signmask) ? signmask + magnitude : signmask - magnitude
end
@inline _ulp_within(distance::U, budget::Integer) where {U<:Unsigned} =
    budget >= 0 && (budget > typemax(U) || distance <= U(budget))
# Complex: ULP-compare the parts independently (max of the two distances).
function _ulp_equal(a::Complex, b::Complex, ulp::Integer)
    return _ulp_equal(real(a), real(b), ulp) && _ulp_equal(imag(a), imag(b), ulp)
end

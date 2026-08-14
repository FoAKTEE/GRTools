# Typed, deterministic schema boundary for a future compiled quadratic-Ricci plan.  This file owns
# metadata only: it does not compile angular couplings, pack jets, or execute an R2 kernel.

"""Marker supertype for the normalization returned by a quadratic Ricci plan."""
abstract type R2Normalization end

"""
    OrderedRicciBilinear()

The ordered generated-kernel contribution `B(hᴬ,hᴮ)`.  The operand slots are deliberately not
interchanged: in general `B(A,B) != B(B,A)`.  Its relation to the symmetric Fréchet Hessian is

`B(A,A) = D²R[A,A]`, `D²R[A,B] = (B(A,B) + B(B,A))/2`.

Hence the cross term in `D²R[A+B,A+B]` is `B(A,B)+B(B,A)`.  This marker prevents an ordered
execution plan from being mislabeled as a symmetric second variation.
"""
struct OrderedRicciBilinear <: R2Normalization end

"""
    SecondVariationRicci()

The symmetric second variation `D²R[A,B] = [ε_Aε_B]R[g+ε_A A+ε_B B]`.  On the diagonal,
`D²R[A,A] = 2[ε²]R[g+ε A]`.  The Q0 plan is ordered and therefore rejects this marker until a
later symmetrized composition explicitly evaluates both operand orders.
"""
struct SecondVariationRicci <: R2Normalization end

"""Validated harmonic mode identifier.  Parent modes may have either sign of `m`."""
struct ModeID
    l::Int16
    m::Int16
    function ModeID(l::Integer, m::Integer)
        0 <= l <= typemax(Int16) || throw(ArgumentError("ModeID: l=$l is outside 0:$(typemax(Int16))"))
        -l <= m <= l || throw(ArgumentError("ModeID: require |m| <= l, got (l=$l,m=$m)"))
        new(Int16(l), Int16(m))
    end
end
ModeID(x::Tuple{Integer,Integer}) = ModeID(x...)
Base.Tuple(x::ModeID) = (Int(x.l), Int(x.m))
Base.isless(a::ModeID, b::ModeID) = (a.l, a.m) < (b.l, b.m)

"""Validated selected output component `(i,l,m)`; targets use the legacy nonnegative-`m` half."""
struct R2TargetID
    i::UInt8
    mode::ModeID
    function R2TargetID(i::Integer, l::Integer, m::Integer)
        mode = ModeID(l, m)
        1 <= i <= 10 || throw(ArgumentError("R2TargetID: component i=$i is outside 1:10"))
        m >= 0 || throw(ArgumentError("R2TargetID: selected targets require m >= 0"))
        ((i <= 7 && isodd(l + m)) || (i > 7 && iseven(l + m))) &&
            throw(ArgumentError("R2TargetID: component parity excludes (i=$i,l=$l,m=$m)"))
        (i in (4, 5, 8, 9) && l < 1) &&
            throw(ArgumentError("R2TargetID: component i=$i requires l >= 1"))
        (i in (7, 10) && l < 2) &&
            throw(ArgumentError("R2TargetID: component i=$i requires l >= 2"))
        new(UInt8(i), mode)
    end
end
R2TargetID(x::Tuple{Integer,Integer,Integer}) = R2TargetID(x...)
Base.Tuple(x::R2TargetID) = (Int(x.i), Int(x.mode.l), Int(x.mode.m))
Base.isless(a::R2TargetID, b::R2TargetID) =
    (a.mode.l, a.mode.m, a.i) < (b.mode.l, b.mode.m, b.i)

"""
    EffectiveBandwidth(h1lmax, dlmax, lmax)

Canonical, nonnegative bandwidth consumed by a future plan compiler.  Raw legacy sentinels are
accepted only by [`legacy_effective_bandwidth`](@ref).
"""
struct EffectiveBandwidth
    h1lmax::Int16
    dlmax::Int16
    lmax::Int16
    function EffectiveBandwidth(h1lmax::Integer, dlmax::Integer, lmax::Integer)
        0 <= h1lmax <= typemax(Int16) || throw(ArgumentError("EffectiveBandwidth: h1lmax must be nonnegative Int16"))
        0 <= dlmax <= h1lmax || throw(ArgumentError("EffectiveBandwidth: require 0 <= dlmax <= h1lmax"))
        0 <= lmax <= h1lmax || throw(ArgumentError("EffectiveBandwidth: require 0 <= lmax <= h1lmax"))
        new(Int16(h1lmax), Int16(dlmax), Int16(lmax))
    end
end

"""Adapt the exact `Ricci.cc` negative-`dlmax` sentinel into a canonical bandwidth."""
function legacy_effective_bandwidth(h1lmax::Integer, dlmax::Integer; lmax::Integer = h1lmax)
    0 <= h1lmax <= typemax(Int16) || throw(ArgumentError(
        "legacy_effective_bandwidth: h1lmax must be in the ModeID Int16 range"))
    lmax >= 0 || throw(ArgumentError("legacy_effective_bandwidth: lmax must be nonnegative"))
    effective_h1lmax = Int(h1lmax)
    effective_dlmax = if dlmax < 0
        dlmax < -effective_h1lmax && throw(ArgumentError(
            "legacy_effective_bandwidth: sentinel removes every first-order multipole"))
        effective_h1lmax += Int(dlmax)
        effective_h1lmax
    else
        dlmax > effective_h1lmax ? effective_h1lmax : Int(dlmax)
    end
    effective_lmax = lmax > effective_h1lmax ? effective_h1lmax : Int(lmax)
    return EffectiveBandwidth(effective_h1lmax, effective_dlmax, effective_lmax)
end

"""The deterministic typed form of `selection_modes(lmax)` (order: `l`, then `m`, then `i`)."""
function quadratic_target_order(lmax::Integer)
    lmax >= 0 || throw(ArgumentError("quadratic_target_order: lmax must be nonnegative"))
    lmax <= typemax(Int16) || throw(ArgumentError(
        "quadratic_target_order: lmax exceeds the ModeID Int16 range"))
    return R2TargetID[R2TargetID(x) for x in selection_modes(lmax)]
end

"""Semantic metadata that must agree between the two ordered input corpora."""
struct H1CorpusConventions
    basis::Symbol
    gauge::Symbol
    slicing::Symbol
    field_normalization::Symbol
    derivative_order::UInt8
    function H1CorpusConventions(; basis::Symbol = :barack_lousto_physical,
                                 gauge::Symbol = :lorenz,
                                 slicing::Symbol = :schwarzschild_t,
                                 field_normalization::Symbol = :physical_non_trace_reversed,
                                 derivative_order::Integer = 2)
        all(!isempty(String(x)) for x in (basis, gauge, slicing, field_normalization)) ||
            throw(ArgumentError("H1CorpusConventions: convention names must be nonempty"))
        0 <= derivative_order <= typemax(UInt8) ||
            throw(ArgumentError("H1CorpusConventions: derivative_order is outside UInt8"))
        new(basis, gauge, slicing, field_normalization, UInt8(derivative_order))
    end
end

const R2_FINGERPRINT_SCHEMA_VERSION = UInt16(1)
const R2_COUPLING_SPEC_VERSION = UInt16(1)
const _FNV128_OFFSET = UInt128(0x6c62272e07bb014262b821756295c58d)
const _FNV128_PRIME = UInt128(0x0000000001000000000000000000013b)

@inline _fp_byte(h::UInt128, b::UInt8) = xor(h, UInt128(b)) * _FNV128_PRIME
function _fp_u64(h::UInt128, x::UInt64)
    @inbounds for s in 0:8:56
        h = _fp_byte(h, UInt8((x >> s) & 0xff))
    end
    return h
end
_fp_i64(h, x::Integer) = _fp_u64(h, reinterpret(UInt64, Int64(x)))
function _fp_text(h::UInt128, x)
    bytes = codeunits(String(x))
    h = _fp_u64(h, UInt64(length(bytes)))
    @inbounds for b in bytes
        h = _fp_byte(h, b)
    end
    return h
end
_fp_float(h::UInt128, x::Float16) = _fp_u64(h, UInt64(reinterpret(UInt16, x)))
_fp_float(h::UInt128, x::Float32) = _fp_u64(h, UInt64(reinterpret(UInt32, x)))
_fp_float(h::UInt128, x::Float64) = _fp_u64(h, reinterpret(UInt64, x))
_fp_float(::UInt128, x::AbstractFloat) = throw(ArgumentError(
    "quadratic plan fingerprints currently support Float16/Float32/Float64 snapshots, got $(typeof(x))"))
_normalization_tag(::OrderedRicciBilinear) = "ordered-R2-bilinear/v1"
_normalization_tag(::SecondVariationRicci) = "symmetric-D2R/v1"

"""
Snapshot of the plan-relevant parts of an H1 corpus.  The fingerprint covers every grid value,
the canonical mode set, `r0`, numeric type, background mass, conventions, normalization, and
coupling-spec version.  Field amplitudes are intentionally excluded so a plan can be reused.  The
kernel-consumed `f` and `fp` arrays are validated against the snapshotted Schwarzschild grid/mass,
so they need not be duplicated in the fingerprint.

`radial_grid` and `modes` are owned copies.  Mutating those public vectors invalidates the stored
fingerprint; [`validate_schema`](@ref) detects that before compatibility or plan-schema creation.
"""
struct BilinearCorpusSchema{T<:AbstractFloat,N<:R2Normalization}
    radial_grid::Vector{T}
    modes::Vector{ModeID}
    r0::T
    background_mass::T
    conventions::H1CorpusConventions
    normalization::N
    coupling_spec_version::UInt16
    fingerprint::UInt128
end

Base.:(==)(a::BilinearCorpusSchema, b::BilinearCorpusSchema) =
    typeof(a) === typeof(b) && a.radial_grid == b.radial_grid && a.modes == b.modes &&
    a.r0 == b.r0 && a.background_mass == b.background_mass &&
    a.conventions == b.conventions && a.normalization == b.normalization &&
    a.coupling_spec_version == b.coupling_spec_version && a.fingerprint == b.fingerprint

function _corpus_fingerprint(::Type{T}, grid, modes, r0, mass, conventions, normalization,
                             spec_version) where {T}
    h = _fp_text(_FNV128_OFFSET, "GRTools.R2Corpus")
    h = _fp_u64(h, UInt64(R2_FINGERPRINT_SCHEMA_VERSION))
    h = _fp_text(h, string(T))
    h = _fp_float(h, r0); h = _fp_float(h, mass)
    for x in (conventions.basis, conventions.gauge, conventions.slicing,
              conventions.field_normalization)
        h = _fp_text(h, x)
    end
    h = _fp_u64(h, UInt64(conventions.derivative_order))
    h = _fp_text(h, _normalization_tag(normalization))
    h = _fp_u64(h, UInt64(spec_version))
    h = _fp_u64(h, UInt64(length(grid)))
    for x in grid; h = _fp_float(h, x); end
    h = _fp_u64(h, UInt64(length(modes)))
    for x in modes
        h = _fp_i64(h, x.l); h = _fp_i64(h, x.m)
    end
    return h
end

function _legacy_mode_set(fields)
    modes = Set{ModeID}()
    N = length(fields.r)
    for d in (fields.h, fields.dh, fields.ddh), ((i, l, m), values) in d
        1 <= i <= 10 || throw(ArgumentError("legacy_corpus_schema: component i=$i is outside 1:10"))
        length(values) == N || throw(DimensionMismatch(
            "legacy_corpus_schema: component (i=$i,l=$l,m=$m) has $(length(values)) values, expected $N"))
        push!(modes, ModeID(l, m))
    end
    isempty(modes) && throw(ArgumentError("legacy_corpus_schema: corpus has no present modes"))
    return sort!(collect(modes))
end

"""Snapshot an existing `H1RadialJetCorpus` without changing any field value or dictionary."""
function legacy_corpus_schema(fields::H1RadialJetCorpus{T};
                              background_mass::Real = one(T),
                              conventions::H1CorpusConventions = H1CorpusConventions(),
                              normalization::N = OrderedRicciBilinear(),
                              coupling_spec_version::Integer = R2_COUPLING_SPEC_VERSION) where
                             {T<:AbstractFloat,N<:R2Normalization}
    T in (Float16, Float32, Float64) || throw(ArgumentError(
        "legacy_corpus_schema: metadata snapshots currently support Float16/Float32/Float64, got $T"))
    Ngrid = length(fields.r)
    Ngrid > 0 || throw(ArgumentError("legacy_corpus_schema: radial grid is empty"))
    length(fields.f) == Ngrid && length(fields.fp) == Ngrid ||
        throw(DimensionMismatch("legacy_corpus_schema: r, f, and fp lengths must agree"))
    all(isfinite, fields.r) && isfinite(fields.r0) ||
        throw(ArgumentError("legacy_corpus_schema: grid and r0 must be finite"))
    all(>(zero(T)), fields.r) && fields.r0 > zero(T) ||
        throw(ArgumentError("legacy_corpus_schema: grid and r0 must be positive"))
    issorted(fields.r) || throw(ArgumentError(
        "legacy_corpus_schema: radial grid must be nondecreasing"))
    0 <= coupling_spec_version <= typemax(UInt16) ||
        throw(ArgumentError("legacy_corpus_schema: coupling_spec_version is outside UInt16"))
    mass = T(background_mass)
    isfinite(mass) && mass > zero(T) || throw(ArgumentError(
        "legacy_corpus_schema: background mass must be finite and positive"))
    all(isfinite, fields.f) && all(isfinite, fields.fp) || throw(ArgumentError(
        "legacy_corpus_schema: Schwarzschild f and fp arrays must be finite"))
    expected_f = @. one(T) - T(2) * mass / fields.r
    expected_fp = @. T(2) * mass / fields.r^2
    fields.f == expected_f || throw(ArgumentError(
        "legacy_corpus_schema: f does not match 1-2M/r for the declared background mass"))
    fields.fp == expected_fp || throw(ArgumentError(
        "legacy_corpus_schema: fp does not match 2M/r² for the declared background mass"))
    modes = _legacy_mode_set(fields)
    grid = copy(fields.r)
    spec = UInt16(coupling_spec_version)
    fp = _corpus_fingerprint(T, grid, modes, fields.r0, mass, conventions, normalization, spec)
    return BilinearCorpusSchema{T,N}(grid, modes, fields.r0, mass, conventions, normalization, spec, fp)
end

Base.hash(s::BilinearCorpusSchema, h::UInt) =
    hash((:GRToolsR2BilinearCorpusSchema, s.fingerprint), h)

"""Recompute a schema fingerprint, rejecting mutation of its owned snapshot vectors."""
function validate_schema(s::BilinearCorpusSchema{T}) where {T}
    issorted(s.modes) && allunique(s.modes) ||
        throw(ArgumentError("BilinearCorpusSchema: modes are not canonical and unique"))
    actual = _corpus_fingerprint(T, s.radial_grid, s.modes, s.r0, s.background_mass,
                                 s.conventions, s.normalization, s.coupling_spec_version)
    actual == s.fingerprint || throw(ArgumentError("BilinearCorpusSchema: snapshot fingerprint mismatch"))
    return true
end

"""
Reject incompatible ordered corpora.  Their mode sets may differ: distinct `A,B` mode sets are a
required bilinear use case and are guarded by their separate fingerprints.
"""
function assert_bilinear_compatible(a::BilinearCorpusSchema{TA}, b::BilinearCorpusSchema{TB}) where {TA,TB}
    validate_schema(a); validate_schema(b)
    TA === TB || throw(ArgumentError("bilinear corpus numeric types differ: $TA vs $TB"))
    a.radial_grid == b.radial_grid || throw(ArgumentError("bilinear corpus radial grids differ"))
    a.r0 == b.r0 || throw(ArgumentError("bilinear corpus r0 values differ"))
    a.background_mass == b.background_mass || throw(ArgumentError("bilinear corpus backgrounds differ"))
    a.conventions == b.conventions || throw(ArgumentError("bilinear corpus conventions differ"))
    typeof(a.normalization) === typeof(b.normalization) ||
        throw(ArgumentError("bilinear corpus normalizations differ"))
    a.coupling_spec_version == b.coupling_spec_version ||
        throw(ArgumentError("bilinear corpus coupling-spec versions differ"))
    return true
end

"""Reject a changed corpus topology/grid relative to an earlier snapshot."""
function assert_corpus_matches(snapshot::BilinearCorpusSchema, fields::H1RadialJetCorpus)
    current = legacy_corpus_schema(fields; background_mass = snapshot.background_mass,
                                   conventions = snapshot.conventions,
                                   normalization = snapshot.normalization,
                                   coupling_spec_version = snapshot.coupling_spec_version)
    current.fingerprint == snapshot.fingerprint ||
        throw(ArgumentError("H1 corpus no longer matches its quadratic-plan snapshot"))
    return true
end

"""Fingerprint-only Q0 plan contract; it contains no angular coefficients or packed field data."""
struct QuadraticPlanSchema{N<:R2Normalization}
    normalization::N
    effective_bandwidth::EffectiveBandwidth
    targets::Vector{R2TargetID}
    corpus_a_fingerprint::UInt128
    corpus_b_fingerprint::UInt128
    coupling_spec_version::UInt16
    fingerprint::UInt128
    function QuadraticPlanSchema{N}(normalization::N, effective_bandwidth,
                                    targets, corpus_a_fingerprint, corpus_b_fingerprint,
                                    coupling_spec_version, fingerprint) where
                                   {N<:R2Normalization}
        normalization isa OrderedRicciBilinear || throw(ArgumentError(
            "QuadraticPlanSchema: only OrderedRicciBilinear is valid for an ordered A/B plan"))
        new{N}(normalization, effective_bandwidth, targets, corpus_a_fingerprint,
               corpus_b_fingerprint, coupling_spec_version, fingerprint)
    end
end


Base.hash(s::QuadraticPlanSchema, h::UInt) =
    hash((:GRToolsR2QuadraticPlanSchema, s.fingerprint), h)


Base.:(==)(a::QuadraticPlanSchema, b::QuadraticPlanSchema) =
    typeof(a) === typeof(b) && a.normalization == b.normalization &&
    a.effective_bandwidth == b.effective_bandwidth && a.targets == b.targets &&
    a.corpus_a_fingerprint == b.corpus_a_fingerprint &&
    a.corpus_b_fingerprint == b.corpus_b_fingerprint &&
    a.coupling_spec_version == b.coupling_spec_version && a.fingerprint == b.fingerprint

function _plan_fingerprint(normalization, bw, targets, afp, bfp, spec)
    h = _fp_text(_FNV128_OFFSET, "GRTools.R2Plan")
    h = _fp_u64(h, UInt64(R2_FINGERPRINT_SCHEMA_VERSION))
    h = _fp_text(h, _normalization_tag(normalization))
    for x in (bw.h1lmax, bw.dlmax, bw.lmax); h = _fp_i64(h, x); end
    h = _fp_u64(h, UInt64(length(targets)))
    for t in targets
        h = _fp_i64(h, t.i); h = _fp_i64(h, t.mode.l); h = _fp_i64(h, t.mode.m)
    end
    low64(x::UInt128) = UInt64(x & UInt128(typemax(UInt64)))
    h = _fp_u64(_fp_u64(h, UInt64(afp >> 64)), low64(afp))
    h = _fp_u64(_fp_u64(h, UInt64(bfp >> 64)), low64(bfp))
    return _fp_u64(h, UInt64(spec))
end

function quadratic_plan_schema(normalization::N, a::BilinearCorpusSchema,
                               b::BilinearCorpusSchema;
                               effective_bandwidth::EffectiveBandwidth,
                               targets = quadratic_target_order(effective_bandwidth.lmax)) where
                              {N<:R2Normalization}
    normalization isa OrderedRicciBilinear || throw(ArgumentError(
        "quadratic_plan_schema: an ordered A/B plan requires OrderedRicciBilinear; " *
        "SecondVariationRicci requires a future explicit two-order symmetrized composition"))
    assert_bilinear_compatible(a, b)
    typeof(a.normalization) === N || throw(ArgumentError("plan/corpus normalization mismatch"))
    typed_targets = R2TargetID[t isa R2TargetID ? t : R2TargetID(t) for t in targets]
    sort!(typed_targets)
    allunique(typed_targets) || throw(ArgumentError("quadratic plan targets must be unique"))
    all(t -> t.mode.l <= effective_bandwidth.lmax, typed_targets) ||
        throw(ArgumentError("quadratic plan target exceeds effective lmax"))
    spec = a.coupling_spec_version
    fp = _plan_fingerprint(normalization, effective_bandwidth, typed_targets,
                           a.fingerprint, b.fingerprint, spec)
    return QuadraticPlanSchema{N}(normalization, effective_bandwidth, typed_targets,
                                  a.fingerprint, b.fingerprint, spec, fp)
end
quadratic_plan_schema(::Type{OrderedRicciBilinear}, args...; kwargs...) =
    quadratic_plan_schema(OrderedRicciBilinear(), args...; kwargs...)
quadratic_plan_schema(::Type{SecondVariationRicci}, args...; kwargs...) =
    quadratic_plan_schema(SecondVariationRicci(), args...; kwargs...)

function validate_schema(s::QuadraticPlanSchema)
    s.normalization isa OrderedRicciBilinear || throw(ArgumentError(
        "QuadraticPlanSchema: ordered plan has a non-ordered normalization"))
    issorted(s.targets) && allunique(s.targets) ||
        throw(ArgumentError("QuadraticPlanSchema: targets are not canonical and unique"))
    actual = _plan_fingerprint(s.normalization, s.effective_bandwidth, s.targets,
                               s.corpus_a_fingerprint, s.corpus_b_fingerprint,
                               s.coupling_spec_version)
    actual == s.fingerprint || throw(ArgumentError("QuadraticPlanSchema: fingerprint mismatch"))
    return true
end

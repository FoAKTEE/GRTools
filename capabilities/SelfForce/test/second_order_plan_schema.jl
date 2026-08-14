using Test
using SelfForce.SecondOrder: H1RadialJetCorpus, OrderedRicciBilinear, SecondVariationRicci,
    ModeID, R2TargetID,
    EffectiveBandwidth, H1CorpusConventions, QuadraticPlanSchema, legacy_effective_bandwidth,
    quadratic_target_order, selection_modes, legacy_corpus_schema, assert_bilinear_compatible,
    assert_corpus_matches, quadratic_plan_schema, validate_schema, read_h1

function _schema_corpus(::Type{T}, fields; reverse_dicts = false) where {T<:AbstractFloat}
    convert_dict(d) = begin
        pairs = collect(d)
        reverse_dicts && reverse!(pairs)
        Dict{NTuple{3,Int},Vector{Complex{T}}}(k => Complex{T}.(v) for (k, v) in pairs)
    end
    r = T.(fields.r)
    f = @. one(T) - T(2) / r
    fp = @. T(2) / r^2
    H1RadialJetCorpus{T}(T(fields.r0), r, f, fp,
                         convert_dict(fields.h), convert_dict(fields.dh), convert_dict(fields.ddh))
end

@testset "quadratic plan Q0 identifiers and canonical bandwidth" begin
    @test Tuple(ModeID(3, -2)) == (3, -2)
    @test Tuple(R2TargetID(1, 0, 0)) == (1, 0, 0)
    @test_throws ArgumentError ModeID(-1, 0)
    @test_throws ArgumentError ModeID(2, 3)
    @test_throws ArgumentError R2TargetID(8, 0, 0)
    @test_throws ArgumentError EffectiveBandwidth(3, -1, 2)
    @test_throws ArgumentError EffectiveBandwidth(3, 4, 2)

    @test legacy_effective_bandwidth(3, -1; lmax = 3) == EffectiveBandwidth(2, 2, 2)
    @test legacy_effective_bandwidth(UInt(3), -1; lmax = UInt(3)) == EffectiveBandwidth(2, 2, 2)
    @test legacy_effective_bandwidth(3, 9; lmax = 9) == EffectiveBandwidth(3, 3, 3)
    @test_throws ArgumentError legacy_effective_bandwidth(0, -1)
    @test_throws ArgumentError legacy_effective_bandwidth(typemin(Int), -typemax(Int); lmax = 0)
    @test_throws ArgumentError legacy_effective_bandwidth(
        UInt(typemax(Int)) + UInt(2), typemin(Int); lmax = UInt(0))
    @test Tuple.(quadratic_target_order(2)) == selection_modes(2)
    @test_throws ArgumentError quadratic_target_order(-1)
    @test_throws ArgumentError quadratic_target_order(typemax(Int16) + 1)
end

@testset "quadratic corpus fingerprints are canonical snapshots" begin
    fields = read_h1(joinpath(@__DIR__, "fixtures", "h1"); l_max = 2)
    reordered = _schema_corpus(Float64, fields; reverse_dicts = true)
    a = legacy_corpus_schema(fields)
    b = legacy_corpus_schema(reordered)

    # Dictionary insertion order and field amplitudes are outside the topology fingerprint.
    @test a.fingerprint == b.fingerprint
    @test a == b
    # Golden value is independently reproduced in fresh Julia processes; Base.hash is forbidden
    # inside the stable fingerprint. Public hash must nevertheless agree whenever schemas are equal.
    @test a.fingerprint == 0xfb702b5eb8378dcc7233ebc6cc2eb6c6
    @test hash(a) == hash(b)
    @test length(Set((a, b))) == 1
    @test length(Dict(a => :a, b => :b)) == 1
    @test a.modes == b.modes
    @test validate_schema(a)
    @test assert_corpus_matches(a, reordered)

    changed_values = _schema_corpus(Float64, fields)
    first(values(changed_values.h))[1] += 1
    @test legacy_corpus_schema(changed_values).fingerprint == a.fingerprint

    changed_grid = _schema_corpus(Float64, fields)
    changed_grid.r[1] = nextfloat(changed_grid.r[1])
    changed_grid.f[1] = 1 - 2 / changed_grid.r[1]
    changed_grid.fp[1] = 2 / changed_grid.r[1]^2
    grid_schema = legacy_corpus_schema(changed_grid)
    @test grid_schema.fingerprint != a.fingerprint
    @test_throws ArgumentError assert_bilinear_compatible(a, grid_schema)
    @test_throws ArgumentError assert_corpus_matches(a, changed_grid)

    changed_modes = _schema_corpus(Float64, fields)
    doomed = last(a.modes)
    for d in (changed_modes.h, changed_modes.dh, changed_modes.ddh)
        for key in collect(keys(d))
            (key[2], key[3]) == Tuple(doomed) && delete!(d, key)
        end
    end
    mode_schema = legacy_corpus_schema(changed_modes)
    @test mode_schema.fingerprint != a.fingerprint
    @test assert_bilinear_compatible(a, mode_schema) # distinct ordered A/B mode sets are allowed
    @test_throws ArgumentError assert_corpus_matches(a, changed_modes)

    convention_schema = legacy_corpus_schema(fields; conventions = H1CorpusConventions(gauge = :other))
    @test convention_schema.fingerprint != a.fingerprint
    @test_throws ArgumentError assert_bilinear_compatible(a, convention_schema)

    float32_schema = legacy_corpus_schema(_schema_corpus(Float32, fields))
    @test float32_schema.fingerprint != a.fingerprint
    @test_throws ArgumentError assert_bilinear_compatible(a, float32_schema)

    mutated = legacy_corpus_schema(fields)
    mutated.radial_grid[1] = nextfloat(mutated.radial_grid[1])
    @test_throws ArgumentError validate_schema(mutated)

    bad_f = _schema_corpus(Float64, fields)
    fill!(bad_f.f, NaN)
    @test_throws ArgumentError legacy_corpus_schema(bad_f)
    @test_throws ArgumentError assert_corpus_matches(a, bad_f)
    bad_fp = _schema_corpus(Float64, fields)
    fill!(bad_fp.fp, Inf)
    @test_throws ArgumentError legacy_corpus_schema(bad_fp)
    @test_throws ArgumentError legacy_corpus_schema(fields; background_mass = 2.0)
    @test_throws ArgumentError legacy_corpus_schema(fields; background_mass = 0.0)
end

@testset "quadratic plan schema fixes normalization, order, and both corpus identities" begin
    fields = read_h1(joinpath(@__DIR__, "fixtures", "h1"); l_max = 2)
    a = legacy_corpus_schema(fields)
    b = legacy_corpus_schema(_schema_corpus(Float64, fields; reverse_dicts = true))
    bw = legacy_effective_bandwidth(2, 2; lmax = 2)

    canonical = quadratic_plan_schema(OrderedRicciBilinear, a, b; effective_bandwidth = bw)
    reversed_input = quadratic_plan_schema(OrderedRicciBilinear(), a, b;
        effective_bandwidth = bw, targets = reverse(selection_modes(2)))
    @test validate_schema(canonical)
    @test canonical.targets == reversed_input.targets == quadratic_target_order(2)
    @test canonical.fingerprint == reversed_input.fingerprint
    @test canonical == reversed_input
    @test canonical.corpus_a_fingerprint == a.fingerprint
    @test canonical.corpus_b_fingerprint == b.fingerprint
    @test hash(canonical) == hash(reversed_input)
    @test length(Set((canonical, reversed_input))) == 1
    @test_throws ArgumentError quadratic_plan_schema(OrderedRicciBilinear(), a, b;
        effective_bandwidth = bw, targets = [(1, 0, 0), (1, 0, 0)])
    symmetric_a = legacy_corpus_schema(fields; normalization = SecondVariationRicci())
    @test_throws ArgumentError quadratic_plan_schema(SecondVariationRicci(), symmetric_a, symmetric_a;
        effective_bandwidth = bw)
    @test_throws ArgumentError QuadraticPlanSchema{SecondVariationRicci}(
        SecondVariationRicci(), bw, copy(canonical.targets), symmetric_a.fingerprint,
        symmetric_a.fingerprint, canonical.coupling_spec_version, canonical.fingerprint)
end

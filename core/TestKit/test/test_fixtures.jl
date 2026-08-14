# Self-tests for fixtures.jl (design §5.6). Grounded by a Data round-trip (write a known table, reload
# it, assert byte/precision-exact equality against the written literals) and by the characterization
# gate catching a >1-ULP perturbation of a closed-form function — never by a heavy legacy dep.

using Test
using TestKit
using Data: set_datadir!

@testset "fixture round-trip via Data" begin
    dir = mktempdir()
    ref = [1.0 2.0 3.0; 4.0 5.0 6.0]
    open(joinpath(dir, "fix.dat"), "w") do io
        for r in 1:size(ref, 1)
            println(io, join(ref[r, :], " "))
        end
    end
    set_datadir!(dir)
    try
        got = load_fixture("fix.dat")
        @test got == ref
    finally
        set_datadir!(nothing)
    end
end

@testset "characterization gate detects a perturbed output" begin
    f = x -> x^2 + 1.0
    inputs = [0.0, 1.0, 2.5, -3.0]
    snap = characterization_snapshot(f, inputs)
    @test snap.outputs == [f(x) for x in inputs]

    # Identical candidate passes bit-exact (ulp=0 default).
    @test characterization_gate(x -> x^2 + 1.0, snap)
    # A perturbed candidate (1 ULP off) FAILS the bit-exact gate...
    perturbed = x -> nextfloat(x^2 + 1.0)
    @test !characterization_gate(perturbed, snap)
    # ...but passes a 1-ULP-budget gate.
    @test characterization_gate(perturbed, snap; ulp = 1)
    @test_throws ArgumentError characterization_gate(perturbed, snap; ulp = -1)
    # A grossly perturbed candidate fails even with a tolerance.
    @test !characterization_gate(x -> x^2 + 1.1, snap; rtol = 1e-6)

    # BigFloat snapshots stay BigFloat and gate bit-exact / ULP-scaled.
    bsnap = characterization_snapshot(x -> x^2 + 1, BigFloat[0.5, 2.0])
    @test eltype(bsnap.outputs) == BigFloat
    @test characterization_gate(x -> x^2 + 1, bsnap)
    @test !characterization_gate(x -> x^2 + 1 + eps(BigFloat), bsnap)

    # Complex outputs gate per-part.
    csnap = characterization_snapshot(x -> x + 2im, ComplexF64[1.0 + 0im, 0.0 + 3im])
    @test characterization_gate(x -> x + 2im, csnap)
    @test !characterization_gate(x -> x + 2.0im + nextfloat(0.0) * im + 1e-3, csnap)
end

@testset "IEEE ULP distance is width-correct and overflow-safe" begin
    for T in (Float16, Float32, Float64)
        oneT = one(T)
        @test TestKit._ulp_distance(oneT, nextfloat(oneT)) == 1
        @test TestKit._ulp_distance(-oneT, prevfloat(-oneT)) == 1
        @test TestKit._ulp_distance(-zero(T), zero(T)) == 0
        @test TestKit._ulp_distance(-nextfloat(zero(T)), nextfloat(zero(T))) == 2

        snap = characterization_snapshot(identity, T[oneT, -oneT])
        @test characterization_gate(x -> x > 0 ? nextfloat(x) : prevfloat(x), snap; ulp = 1)
        @test !characterization_gate(x -> x > 0 ? nextfloat(nextfloat(x)) : prevfloat(prevfloat(x)),
                                    snap; ulp = 1)
    end

    # The cross-sign distance is larger than typemax(Int64) for Float64; computing it must not
    # overflow or turn negative, and an arbitrarily large Integer budget remains supported.
    d = TestKit._ulp_distance(-floatmax(Float64), floatmax(Float64))
    @test d isa UInt64
    @test d > UInt64(typemax(Int64))
    snap = characterization_snapshot(identity, [-floatmax(Float64)])
    @test !characterization_gate(_ -> floatmax(Float64), snap; ulp = typemax(Int64))
    @test characterization_gate(_ -> floatmax(Float64), snap; ulp = big(typemax(UInt64)) + 1)
end

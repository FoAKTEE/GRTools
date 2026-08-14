# Grounded tests for SelfForce.SecondOrder mode machinery (a_il + selection_modes).
#
# GROUNDING — fully ANALYTIC + core-characterization, no heavy legacy dep:
#
# (1) a_il is a RE-EXPORT of NumericsCore.a_il (committed-green core). We assert object IDENTITY
#     (=== the core function) AND ground the VALUES against (a) the core implementation on an
#     i∈1:10, l∈0:6 grid at Float64+BigFloat, and (b) the analytic closed form transcribed from
#     ref-code/SecondOrderRicci/src/h1.cc L31-51: base 2^(-1/2); ×1/√(l(l+1)) for i∈{4,5,8,9};
#     ×1/√((l-1)l(l+1)(l+2)) for i∈{7,10}; unscaled for i∈{1,2,3,6}.
#
# (2) selection_modes: the valid (i,l,m) target set hand-enumerated directly from the C++ parity/
#     floor rules in ref-code/SecondOrderRicci/src/Ricci.cc L86-99 (skip i≤7 when l+m odd; skip i>7
#     when l+m even; skip i∈{4,5,8,9} when l<1; skip i∈{7,10} when l<2), in driver push-order
#     (l outer, then m, then i). Expected lmax=2 (29 modes) + the l=3 block (→49 for lmax=3) are the
#     verbatim hand lists; invariants are re-checked for every emitted mode up to lmax=5.

using Test
using SelfForce.SecondOrder: a_il, selection_modes
using NumericsCore: a_il as core_a_il

@testset "a_il is the NumericsCore re-export (identity + values)" begin
    @test a_il === core_a_il
    base = 1 / sqrt(2.0)
    for i in (1, 2, 3, 6), l in (0, 1, 2, 5)
        @test a_il(i, l) == core_a_il(i, l)
        @test a_il(i, l) ≈ base rtol = 1e-15
    end
    for i in (4, 5, 8, 9), l in (1, 2, 3, 7)
        @test a_il(i, l) == core_a_il(i, l)
        @test a_il(i, l) ≈ base / sqrt(l * (l + 1)) rtol = 1e-15
    end
    for i in (7, 10), l in (2, 3, 5, 8)
        @test a_il(i, l) == core_a_il(i, l)
        @test a_il(i, l) ≈ base / sqrt((l - 1) * l * (l + 1) * (l + 2)) rtol = 1e-15
    end
end

@testset "a_il BigFloat precision path" begin
    setprecision(BigFloat, 256) do
        a = a_il(BigFloat, 7, 3)
        @test a isa BigFloat
        @test a ≈ 1 / sqrt(big(2)) / sqrt(big((3 - 1) * 3 * (3 + 1) * (3 + 2))) rtol = 1e-60
        @test a == core_a_il(BigFloat, 7, 3)
    end
end

# Hand-enumerated push-order list for lmax=2 (tuples are (i,l,m)), from Ricci.cc L86-99.
const _MODES_LMAX2 = [
    (1, 0, 0), (2, 0, 0), (3, 0, 0), (6, 0, 0),
    (8, 1, 0), (9, 1, 0),
    (1, 1, 1), (2, 1, 1), (3, 1, 1), (4, 1, 1), (5, 1, 1), (6, 1, 1),
    (1, 2, 0), (2, 2, 0), (3, 2, 0), (4, 2, 0), (5, 2, 0), (6, 2, 0), (7, 2, 0),
    (8, 2, 1), (9, 2, 1), (10, 2, 1),
    (1, 2, 2), (2, 2, 2), (3, 2, 2), (4, 2, 2), (5, 2, 2), (6, 2, 2), (7, 2, 2),
]

@testset "selection_modes exact members + count (lmax=2)" begin
    got = selection_modes(2)
    @test length(got) == 29
    @test got == _MODES_LMAX2
end

@testset "selection_modes lmax=3 (prefix stability + l=3 block)" begin
    got = selection_modes(3)
    @test length(got) == 49
    @test got[1:29] == _MODES_LMAX2          # push-order is stable under raising lmax
    l3block = [
        (8, 3, 0), (9, 3, 0), (10, 3, 0),                                      # m=0, l+m odd
        (1, 3, 1), (2, 3, 1), (3, 3, 1), (4, 3, 1), (5, 3, 1), (6, 3, 1), (7, 3, 1),  # m=1 even
        (8, 3, 2), (9, 3, 2), (10, 3, 2),                                      # m=2 odd
        (1, 3, 3), (2, 3, 3), (3, 3, 3), (4, 3, 3), (5, 3, 3), (6, 3, 3), (7, 3, 3),  # m=3 even
    ]
    @test got[30:end] == l3block
end

@testset "selection_modes parity/floor invariants (lmax=5)" begin
    for (i, l, m) in selection_modes(5)
        @test 1 <= i <= 10 && 0 <= m <= l
        i <= 7 ? (@test iseven(l + m)) : (@test isodd(l + m))
        i in (4, 5, 8, 9) && (@test l >= 1)
        i in (7, 10) && (@test l >= 2)
    end
end

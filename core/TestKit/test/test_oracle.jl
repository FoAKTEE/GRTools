# Self-tests for oracle.jl (design §5.6). The oracle layer is infrastructure, not physics, so it is
# grounded by INVARIANTS of its own behavior — the macro passes/fails correctly and records the frozen
# tolerance, the tolerance-class ordering, and the can-only-tighten guard — never by characterization
# against a physics package. Oracle = closed-form numeric identities (1.0 vs 1.0+1e-13 under a frozen
# rtol, the class-ordering invariant, and the tighten-guard catching a loosened tol).

using Test
using TestKit

@testset "@approx_oracle passes/fails correctly + records tol" begin
    reset_oracle_log!()
    # Crafted EQUAL within the frozen tol → true.
    @test @approx_oracle(1.0, 1.0 + 1e-13; rtol = 1e-12, label = "eq")
    # Crafted UNEQUAL beyond the frozen tol → false (NOT a test failure).
    @test !(@approx_oracle(1.0, 2.0; rtol = 1e-12, label = "neq"))
    # atol-only path: 1e-9 apart, atol 1e-6 → true.
    @test @approx_oracle(0.0, 1e-9; atol = 1e-6, label = "atol")

    recs = oracle_records()
    @test length(recs) == 3
    # The recorded tolerance equals the literal passed (frozen, not re-derived).
    @test recs[1].label == "eq" && recs[1].passed && recs[1].rtol == 1e-12
    @test recs[2].label == "neq" && !recs[2].passed
    @test recs[3].atol == 1e-6
    # Default (un-frozen) rtol is logged as Base's √eps default, atol 0.
    reset_oracle_log!()
    @test @approx_oracle(1.0, 1.0)
    @test oracle_records()[1].rtol ≈ sqrt(eps(Float64))
    @test oracle_records()[1].atol == 0.0
end

@testset "tolerance-class ordering invariants" begin
    m = tolerance_class(:machine).rtol
    s = tolerance_class(:spectral).rtol
    p = tolerance_class(:physics).rtol
    g = tolerance_class(:gpu).rtol
    # tightest → loosest; the binding invariant is :machine strictly tighter than :physics.
    @test m < s < p < g
    @test m < p
    @test_throws ArgumentError tolerance_class(:nonsense)
end

@testset "class_tol <= original_tol tighten-guard" begin
    # Tightening (or equal) is allowed.
    @test assert_can_only_tighten(1e-12, 1e-10)
    @test assert_can_only_tighten(1e-10, 1e-10)
    # A class can resolve through its rtol.
    @test assert_can_only_tighten(:machine, 1e-10)
    # LOOSENING is caught.
    @test_throws AssertionError assert_can_only_tighten(1e-6, 1e-10)
    @test_throws AssertionError assert_can_only_tighten(:physics, 1e-10)
end

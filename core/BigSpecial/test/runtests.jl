using Test
using BigSpecial
# NOTE: the former LIVE Teukolsky characterization oracle is now BAKED into fixture literals
# inside the per-module test files (captured once at 256-bit; bit-identical `==` preserved),
# so the test suite no longer depends on any packages/* legacy package.
import SpecialFunctions   # for direct upstream-equality assertions (Re z ≥ 1/2 branch)
import HypergeometricFunctions   # ₂F₁ oracle for the safe_2f1/zinf closed-form / pass-through checks
import NumericsCore   # self-contained kummer_1F1 cross-check oracle (distinct primitive)

@testset "BigSpecial" begin
    @testset "smoke" begin
        @test isdefined(Main, :BigSpecial)
        @test BigSpecial isa Module
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module test files are added as each primitive is TDD-ported + characterization-gated
    # against BAKED Teukolsky-oracle fixtures (gamma_big/pochhammer/2F1) and analytic identities.
    include("test_pochhammer.jl")
    include("test_biggamma.jl")
    include("test_digamma.jl")
    include("test_hyp2f1.jl")
    include("test_kummer.jl")
    include("test_legendre_p3.jl")
end

# §7.7 binding CI hygiene gates for Spacetime (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2,
# §3.2 TIER 2 [6]: Spacetime → NumericsCore only — spine root, "zero heavy deps, closure ==
# {LinearAlgebra, NumericsCore}"; Roots DROPPED — A42 vendored in NumericsCore (§4.1 Rev-2 fix,
# the S1-restored invariant)).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # spine-root budget: EXACTLY the §3.2 pin {LinearAlgebra, NumericsCore}
        allow=["LinearAlgebra", "NumericsCore"],
        forbid=["Roots", "HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "WignerSymbols", "BigSpecial", "Wigner", "Kernels", "Data",
                "SpecialFunctions", "HypergeometricFunctions", "NLsolve", "SciMLBase", "QuadGK"])
    # gate 2 (closure-load): run AFTER `using Spacetime`; the runtime closure must stay exactly
    # {LinearAlgebra, NumericsCore} — none of these may appear in Base.loaded_modules. Lists only
    # poison nothing in this test env legitimately loads (test target = Test + TestKit + Data).
    TestKit.assert_not_loaded(["Roots", "HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA",
                               "AMDGPU", "KernelAbstractions", "WignerSymbols",
                               "SpecialFunctions", "HypergeometricFunctions"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [Spacetime]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

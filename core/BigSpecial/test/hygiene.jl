# §7.7 binding CI hygiene gates for BigSpecial (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2,
# §3.2 TIER 1 [1]: BigSpecial → NumericsCore, plus the sanctioned HypGeo/SpecialFunctions closure —
# the ONLY place the heavy special-function stack is allowed below tier 3).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["NumericsCore", "HypergeometricFunctions", "SpecialFunctions"],
        forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "WignerSymbols", "Kernels", "Data", "Roots", "NLsolve",
                "SciMLBase", "QuadGK"])
    # gate 2 (closure-load): run AFTER `using BigSpecial`; SpecialFunctions/HypergeometricFunctions
    # are legitimately loaded here (sanctioned hard deps), so they are NOT in this poison list.
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions", "WignerSymbols", "Roots", "NLsolve",
                               "SciMLBase"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [BigSpecial]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

# §7.7 binding CI hygiene gates for TestKit itself (design GRTOOLS_DESIGN_CAPABILITY.md §7.7
# gates 1+2, §3.2 TIER 1 [5]: TestKit → Data (test-only harness), plus Test/Printf/TOML stdlib).
# Budget: allow == {NumericsCore, Data, Test, Printf, TOML} — `Test` is a normal [deps] because
# TestKit IS the harness (§5.6); TOML is the stdlib backing the §7.7 static dep-budget gate.
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["NumericsCore", "Data", "Test", "Printf", "TOML"],
        forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Kernels", "WignerSymbols", "Roots", "NLsolve",
                "SciMLBase", "SpecialFunctions", "HypergeometricFunctions", "QuadGK"])
    # gate 2 (closure-load): run AFTER `using TestKit`; lists only poison nothing in this
    # test env legitimately loads (Data's HDF5/JLD2 weakdep-exts are NOT triggered here).
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA",
                               "AMDGPU", "KernelAbstractions", "WignerSymbols"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [TestKit]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

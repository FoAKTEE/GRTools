# §7.7 binding CI hygiene gates for MetricReconstruction (design GRTOOLS_DESIGN_CAPABILITY.md
# §7.7 gates 1+2, §3.2 TIER 4 [9]: MetricReconstruction → Perturbation, Harmonics, Wigner,
# Spacetime, NumericsCore (+QuadGK); §4.4 "EXCLUDED: HDF5, Symbolics, KernelAbstractions").
# BigSpecial is a SANCTIONED Rev-3 amendment (terminating-2F1/Γ primitives used by the Lorenz
# submodule) on top of the §3.2 line.
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["LinearAlgebra", "Printf", "NumericsCore", "BigSpecial", "Wigner", "Spacetime",
               "Harmonics", "Perturbation", "QuadGK"],
        forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Kernels", "Data", "WignerSymbols", "Roots", "NLsolve",
                "SciMLBase"])
    # gate 2 (closure-load): run AFTER `using MetricReconstruction`. Data is EXCLUDED from this
    # list only because the TestKit harness itself rides Data in the test env — Data's absence
    # from the RUNTIME [deps] is asserted TOML-statically above.
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions", "WignerSymbols", "Roots", "NLsolve",
                               "SciMLBase"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [MetricReconstruction]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

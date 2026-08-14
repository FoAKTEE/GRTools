# §7.7 binding CI hygiene gates for Perturbation (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates
# 1+2, §3.2 TIER 3 [8]: Perturbation → Spacetime, Harmonics, BigSpecial, NumericsCore — the
# LIGHT-DEP FIREWALL: "HDF5-free, Symbolics-free, external-ODE-free" (§4.3, §7.5 binding:
# radial_numint stays on the VENDORED integrator, never OrdinaryDiffEq).
# NOTE: the pre-existing "firewall closure-load gate" testset in runtests.jl is kept UNTOUCHED;
# this file adds the TOML-static budget + the TestKit-hosted closure check on top.
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # HypGeo/SpecialFunctions: the sanctioned NOT-poison closure through BigSpecial (§4.3).
        allow=["LinearAlgebra", "Printf", "NumericsCore", "BigSpecial", "Spacetime", "Harmonics",
               "HypergeometricFunctions", "SpecialFunctions"],
        forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Kernels", "Data", "WignerSymbols", "Roots", "NLsolve",
                "SciMLBase"])
    # gate 2 (closure-load): run AFTER `using Perturbation`. Data is EXCLUDED from this list only
    # because the TestKit harness itself rides Data in the test env — Data's absence from the
    # RUNTIME [deps] is asserted TOML-statically above.
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions", "WignerSymbols", "Roots", "NLsolve",
                               "SciMLBase"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [Perturbation]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

# §7.7 binding CI hygiene gates for NumericsCore (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gate 1:
# "NumericsCore stays LinearAlgebra-only"; §3.2 TIER 0 [0]: the UNIQUE SINK; §5.1: "A CI dep-budget
# gate prevents core re-growth"). Complements the narrower GPU-axis src/deps scan already in
# test_backend_iface.jl (kept untouched).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # the sink budget: LinearAlgebra EXACTLY — any new hard dep is core re-growth.
        allow=["LinearAlgebra"],
        # forbid everything heavy + every GrTools package (nothing may point UP from the sink).
        forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU", "Metal",
                "oneAPI", "GPUArrays", "KernelAbstractions", "WignerSymbols", "Roots", "NLsolve",
                "SciMLBase", "QuadGK", "ForwardDiff", "SpecialFunctions", "HypergeometricFunctions",
                "BigSpecial", "Wigner", "Kernels", "Data", "TestKit", "Spacetime", "Harmonics",
                "Perturbation", "MetricReconstruction", "PostNewtonian", "ModifiedGravity",
                "SelfForce", "EMRIWaveforms", "Hybrid", "AccretionDisc"])
    # gate 2 (closure-load): run AFTER `using NumericsCore`; lists only poison nothing in this
    # test env legitimately loads (test target = Test + TestKit + Data, all poison-free).
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions", "WignerSymbols", "SpecialFunctions",
                               "HypergeometricFunctions", "Roots", "NLsolve", "SciMLBase",
                               "QuadGK", "ForwardDiff"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [NumericsCore]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

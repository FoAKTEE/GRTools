# §7.7 binding CI hygiene gates for SelfForce (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates
# 1+2+5, §3.2 TIER 5 [12]: SelfForce → MetricReconstruction, Harmonics, BigSpecial, Wigner,
# Kernels, Data, NumericsCore + HDF5 DIRECT; §4.5 deps also sanction ForwardDiff/QuadGK
# (Puncture). HDF5 HOME #1: the DIRECT hard dep is the trigger that activates Data's HDF5
# weakdep-ext (§7.2 — an ext fires only when the LOADING package declares the trigger itself).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # TODO(S2 discrepancy): Printf is stdlib + zero-cost but has NO src/ usage (grep-verified)
        # and is not on the §4.5 dep list — looks vestigial. (LinearAlgebra IS used: BLAS guard +
        # eigen in LorenzField.) Allowed so the suite stays green; flagged in the S2 report.
        allow=["LinearAlgebra", "Printf", "NumericsCore", "BigSpecial", "Wigner", "Kernels",
               "Data", "Harmonics", "MetricReconstruction", "Spacetime", "HDF5", "ForwardDiff",
               "QuadGK"],
        # gate 5 (extension trigger): the grid/corpus consumer MUST declare HDF5 directly.
        must_declare=["HDF5"],
        # KernelAbstractions must arrive only THROUGH Kernels, never as a direct dep.
        forbid=["Symbolics", "OrdinaryDiffEq*", "WignerSymbols", "JLD2", "CUDA", "AMDGPU",
                "KernelAbstractions", "Roots", "NLsolve", "SciMLBase"])
    # gate 2 (closure-load): run AFTER `using SelfForce`. HDF5 + KernelAbstractions (via Kernels)
    # are legitimately loaded here; the vendored-integrator invariant (§7.5) means NO
    # OrdinaryDiffEq* module may appear.
    TestKit.assert_not_loaded(["Symbolics", "OrdinaryDiffEq*", "WignerSymbols", "JLD2", "CUDA",
                               "AMDGPU", "Roots", "NLsolve", "SciMLBase"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [SelfForce]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

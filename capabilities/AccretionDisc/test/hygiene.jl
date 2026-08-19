# §7.7 binding CI hygiene gates for AccretionDisc (design GRTOOLS_DESIGN_CAPABILITY.md §7.7
# gates 1+2+5, §3.2 TIER 6 [15]: AccretionDisc → MetricReconstruction, Perturbation, Harmonics,
# Spacetime, Data, NumericsCore + HDF5 DIRECT (Force submodule only — the Data HDF5-ext trigger,
# §5.4/§7.7 gate 5) + the SparseArrays stdlib (SpiralWaves sparse implicit solver); §7.3 L3 leak fix: the a=0/a≠0 Lorenz radial physics is reached ONLY through
# MetricReconstruction.Lorenz — a SelfForce or Kernels edge is the exact leak §7.3 forbids).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["LinearAlgebra", "Printf", "SparseArrays", "NumericsCore", "Data", "Spacetime",
               "Harmonics", "Perturbation", "MetricReconstruction", "HDF5",
               "SpecialFunctions"],
        forbid=["SelfForce", "Kernels", "Symbolics", "JLD2", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "WignerSymbols", "Roots", "NLsolve", "SciMLBase",
                "BigSpecial", "Wigner"],
        must_declare=["HDF5"])
    # gate 2 (closure-load): run AFTER `using AccretionDisc`. HDF5 is legitimately loaded (hard
    # dep, the Data-ext trigger); Data rides both the runtime [deps] and the TestKit harness.
    TestKit.assert_not_loaded(["SelfForce", "Kernels", "Symbolics", "JLD2", "OrdinaryDiffEq",
                               "OrdinaryDiffEqVerner", "CUDA", "AMDGPU", "KernelAbstractions",
                               "WignerSymbols", "Roots", "NLsolve", "SciMLBase"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [AccretionDisc]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

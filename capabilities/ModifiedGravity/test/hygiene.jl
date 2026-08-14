# §7.7 binding CI hygiene gates for ModifiedGravity (design GRTOOLS_DESIGN_CAPABILITY.md §7.7
# gates 1+2, §3.2 TIER 4 [11]: ModifiedGravity → Perturbation, Harmonics, NumericsCore
# (+Printf, Data opt); HDF5-free + Symbolics-free (§4.7 — the only data surface is 12-col .dat).
# Move C (NPGHP campaign): the former MetricReconstruction dep (h44 GHP fill jet calculus) is
# DELETED — the jet algebra is now the canonical NumericsCore.jets, and radial_jet moved to
# Perturbation.Radial. The closure below is back on the pure §3.2 TIER-4 line (MR absent).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["Printf", "NumericsCore", "Harmonics", "Perturbation", "Data"],
        # MetricReconstruction added to `forbid`: Move C DELETED the MG→MR dep — this decisively
        # asserts it is gone from [deps]/[weakdeps]/[sources], not merely un-allowed.
        forbid=["MetricReconstruction",
                "HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Kernels", "WignerSymbols", "Roots", "NLsolve",
                "SciMLBase"])
    # gate 2 (closure-load): run AFTER `using ModifiedGravity` (Data is a sanctioned hard dep here).
    # MetricReconstruction is asserted NOT loaded — the runtime closure no longer pulls it (Move C).
    TestKit.assert_not_loaded(["MetricReconstruction",
                               "HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions", "WignerSymbols", "Roots", "NLsolve",
                               "SciMLBase"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [ModifiedGravity]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

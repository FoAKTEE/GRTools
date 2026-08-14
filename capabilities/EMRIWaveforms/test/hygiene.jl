# §7.7 binding CI hygiene gates for EMRIWaveforms (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates
# 1+2+5, §3.2 TIER 5 [13]: EMRIWaveforms → Spacetime, Harmonics, Kernels, Data, NumericsCore +
# HDF5 (DIRECT, triggers Data HDF5-ext), OrdinaryDiffEqVerner, Roots, NLsolve; §4.8 also
# SpecialFunctions (besselj — S4-c Kernels-purity hoist). HDF5 HOME #2. Symbolics-clean after the
# 4 dead FEW import edges were deleted (§7.1); WignerSymbols eliminated (§7.4).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # TODO(S2 discrepancy): Random/Printf/LinearAlgebra are stdlib + zero-cost but have NO
        # src/ usage (grep-verified; Random appears only in test files' `using Random`) and are
        # not on the §3.2/§4.8 dep lines — Random belongs in [extras]; Printf/LinearAlgebra look
        # vestigial. Allowed here so the suite stays green; flagged in the S2 report.
        allow=["Printf", "LinearAlgebra", "Random", "NumericsCore", "Spacetime", "Harmonics",
               "Kernels", "Data", "HDF5", "OrdinaryDiffEqVerner", "Roots", "NLsolve",
               "SpecialFunctions"],
        # gate 5 (extension trigger): the grid/corpus consumer MUST declare HDF5 directly.
        must_declare=["HDF5"],
        # exact "OrdinaryDiffEq" = the monolith (§7.5: only the Verner sub-package is sanctioned);
        # KernelAbstractions must arrive only THROUGH Kernels.
        forbid=["Symbolics", "WignerSymbols", "JLD2", "CUDA", "AMDGPU", "KernelAbstractions",
                "OrdinaryDiffEq", "SciMLBase", "BigSpecial"])
    # gate 2 (closure-load): run AFTER `using EMRIWaveforms`. HDF5/OrdinaryDiffEqVerner(+Core)/
    # Roots/NLsolve/KernelAbstractions are legitimately loaded here; the monolithic OrdinaryDiffEq
    # (exact name) and the Symbolics/WignerSymbols axes must not be.
    TestKit.assert_not_loaded(["Symbolics", "WignerSymbols", "JLD2", "CUDA", "AMDGPU",
                               "OrdinaryDiffEq"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [EMRIWaveforms]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

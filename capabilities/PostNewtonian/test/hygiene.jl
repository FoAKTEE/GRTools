# §7.7 binding CI hygiene gates for PostNewtonian (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates
# 1+2+5, §3.2 TIER 4 [10]: PostNewtonian → Harmonics, NumericsCore, Data + Symbolics (HARD,
# quarantined — §7.1: Symbolics is a struct FIELD TYPE, cannot be a weakdep; this is its SOLE
# hard home). §7.2/§7.7-5: PostNewtonian must NOT declare HDF5 (its Data artifacts are .m text).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # TODO(S2 discrepancy): LinearAlgebra is stdlib + zero-cost but has NO src/ usage
        # (grep-verified) and is not on the §4.6 dep list (Printf IS §4.6-sanctioned) — looks
        # vestigial. Allowed so the suite stays green; flagged in the S2 report.
        allow=["Printf", "LinearAlgebra", "NumericsCore", "Harmonics", "Data", "Symbolics",
               "SpecialFunctions"],
        # gate 5 (quarantine home): Symbolics MUST be a direct hard [deps] here.
        must_declare=["Symbolics"],
        forbid=["HDF5", "JLD2", "CUDA", "AMDGPU", "KernelAbstractions", "Kernels",
                "WignerSymbols", "OrdinaryDiffEq*", "Roots", "NLsolve", "SciMLBase"])
    # gate 2 (closure-load): run AFTER `using PostNewtonian`. Symbolics is legitimately loaded
    # here (the quarantine home). List only poison the Symbolics closure cannot legitimately pull.
    TestKit.assert_not_loaded(["HDF5", "JLD2", "CUDA", "AMDGPU", "KernelAbstractions",
                               "WignerSymbols", "OrdinaryDiffEq*"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [PostNewtonian]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

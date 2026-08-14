# §7.7 binding CI hygiene gates for Data (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2,
# §3.2 TIER 1 [4]: Data → NumericsCore; §5.4/§7.2: HDF5/JLD2 reader-writer BODIES live in WEAKDEP
# extensions — the exts fire only when a tier-5/6 consumer declares the trigger in its own [deps]).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # Printf: stdlib, in current [deps] (formatted writers) — sanctioned-light.
        allow=["NumericsCore", "DelimitedFiles", "Artifacts", "Printf"],
        # HDF5/JLD2 are deliberately NOT in forbid (they are the sanctioned WEAK axes); their
        # exclusion from hard [deps] is enforced by the closed allow budget above, and the weak
        # surface is pinned to exactly {HDF5, JLD2} by weakdeps_allow.
        forbid=["Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU", "KernelAbstractions",
                "WignerSymbols", "Kernels", "Roots", "NLsolve", "SciMLBase",
                "SpecialFunctions", "HypergeometricFunctions", "QuadGK"],
        weakdeps_allow=["HDF5", "JLD2"])
    # gate 2 (closure-load): run AFTER `using Data`. HDF5/JLD2 are EXCLUDED from this list — this
    # test env legitimately loads them ([extras]) to fire the weakdep-exts under round-trip test.
    TestKit.assert_not_loaded(["Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions", "WignerSymbols"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [Data]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

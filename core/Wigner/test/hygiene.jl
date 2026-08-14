# §7.7 binding CI hygiene gates for Wigner (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2,
# §3.2 TIER 1 [2]: Wigner → NumericsCore ONLY; §7.4/L4: the WignerSymbols PACKAGE is ELIMINATED
# from every runtime path — the coupling machinery is hand-rolled).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["NumericsCore"],
        # WignerSymbols is the DEFINING poison here (§7.4): runtime-eliminated. It remains a
        # SANCTIONED test-only oracle in [extras] — extras are deliberately not scanned.
        forbid=["WignerSymbols", "HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Kernels", "Data", "SpecialFunctions",
                "HypergeometricFunctions", "Roots", "NLsolve", "SciMLBase", "QuadGK"])
    # gate 2 (closure-load): run AFTER `using Wigner`. WignerSymbols is EXCLUDED from this list —
    # this test env legitimately loads it as the characterization oracle ([extras]); its runtime
    # elimination is asserted TOML-statically above instead.
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "KernelAbstractions"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [Wigner]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

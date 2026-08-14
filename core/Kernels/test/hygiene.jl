# §7.7 binding CI hygiene gates for Kernels (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2,
# §3.2 TIER 1 [3]: Kernels → NumericsCore + KernelAbstractions; GPU backends CUDA/AMDGPU are
# WEAKDEP extensions only (§5.5) — never hard deps).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        # LinearAlgebra: stdlib, in current [deps] (mul!/GEMM hosts) — sanctioned-light.
        allow=["NumericsCore", "KernelAbstractions", "LinearAlgebra"],
        # CUDA/AMDGPU are NOT in forbid (they are sanctioned WEAK); their exclusion from hard
        # [deps] is enforced by the closed allow budget above.
        forbid=["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "WignerSymbols", "Data", "Roots",
                "NLsolve", "SciMLBase", "SpecialFunctions", "HypergeometricFunctions", "QuadGK"],
        weakdeps_allow=["CUDA", "AMDGPU"])
    # gate 2 (closure-load): run AFTER `using Kernels`. KernelAbstractions is legitimately loaded
    # (hard dep); the CPU-only test env must NOT have pulled a GPU stack.
    TestKit.assert_not_loaded(["HDF5", "JLD2", "Symbolics", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                               "WignerSymbols"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [Kernels]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

# §7.7 binding CI hygiene gates for Hybrid (design GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2+5,
# §3.2 TIER 5 [14]: Hybrid → Spacetime, BigSpecial, Data, NumericsCore + OrdinaryDiffEqVerner,
# SciMLBase; "no HDF5; JLD2 DIRECT, triggers Data JLD2-ext". Top-of-DAG leaf: the ODE/SciML
# precompile floor is confined here (§4.9).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["Printf", "LinearAlgebra", "NumericsCore", "BigSpecial", "Spacetime", "Data",
               "OrdinaryDiffEqVerner", "SciMLBase", "JLD2"],
        # gate 5 (extension trigger): Hybrid MUST declare JLD2 directly (fires Data's JLD2-ext).
        must_declare=["JLD2"],
        # §7.7-5: Hybrid does NOT declare HDF5. Exact "OrdinaryDiffEq" = the monolith (§7.5:
        # WaSABI's OrdinaryDiffEq was swapped to Verner-only).
        forbid=["HDF5", "Symbolics", "WignerSymbols", "Kernels", "KernelAbstractions", "CUDA",
                "AMDGPU", "OrdinaryDiffEq", "Roots", "NLsolve"])
    # gate 2 (closure-load): run AFTER `using Hybrid`. JLD2/OrdinaryDiffEqVerner(+Core)/SciMLBase
    # are legitimately loaded here; HDF5 and the GPU/Kernels axes must not be.
    TestKit.assert_not_loaded(["HDF5", "Symbolics", "WignerSymbols", "Kernels",
                               "KernelAbstractions", "CUDA", "AMDGPU", "OrdinaryDiffEq"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [Hybrid]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

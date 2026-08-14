# §7.7 binding CI hygiene gates for ScalarSelfForce (elliptic-ssf mission: DESIGN.md §1.3 budget,
# GLOBAL_DAG node `pkg-scaffold-scalarselfforce`, schedule W1-B).
#
# Budget (closed): EXACTLY {LinearAlgebra, Printf, NumericsCore, BigSpecial, Spacetime,
# EllipticCore, Data}. Poison axes:
#   * GSL / live-C-oracle channel — the EffectiveSource C oracles (kerr-circular.c /
#     kerr-circular-comoving.c) enter ONLY as baked fixtures ("legacy oracles are NEVER live test
#     deps"); GSL*/Libdl in [deps] would be the live-ccall channel, so both are forbidden.
#   * CUDA / AMDGPU / KernelAbstractions — GPU is banned for this mission.
#   * SpecialFunctions / HypergeometricFunctions as DIRECT deps — special functions must come
#     through BigSpecial (TOML-static ban only: BigSpecial itself legitimately loads both, so the
#     closure-load gate below must not list them).
#   * external Krylov/AMG stacks (EllipticCore vendors GMRES/MG/Schwarz in-tree) and QuadGK
#     (open question Q7: default = baked quadrature fixtures).
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["LinearAlgebra", "Printf", "NumericsCore", "BigSpecial", "Spacetime",
               "EllipticCore", "Data"],
        forbid=["GSL*", "Libdl", "CUDA", "AMDGPU", "KernelAbstractions",
                "SpecialFunctions", "HypergeometricFunctions",
                "Krylov", "IterativeSolvers", "AlgebraicMultigrid",
                "HDF5", "JLD2", "QuadGK", "OrdinaryDiffEq*", "ForwardDiff",
                "Symbolics", "Roots", "NLsolve", "SciMLBase", "WignerSymbols"])
    # gate 2 (closure-load): run AFTER `using ScalarSelfForce`; lists only poison that nothing in
    # this test env legitimately loads. SpecialFunctions/HypergeometricFunctions are excluded here
    # (legitimate through BigSpecial) and Libdl is excluded (loaded by LinearAlgebra's BLAS
    # machinery) — their direct-dep bans are asserted TOML-statically above.
    TestKit.assert_not_loaded(["GSL*", "CUDA", "AMDGPU", "KernelAbstractions",
                               "Krylov", "IterativeSolvers", "AlgebraicMultigrid",
                               "HDF5", "JLD2", "QuadGK", "OrdinaryDiffEq*", "ForwardDiff",
                               "Symbolics", "WignerSymbols"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [ScalarSelfForce]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

# Binding CI hygiene gates for EllipticCore (elliptic-ssf DESIGN.md §1.3 + §6.4; GrTool
# GRTOOLS_DESIGN_CAPABILITY.md §7.7 gates 1+2). Budget: {LinearAlgebra, SparseArrays,
# NumericsCore} — zero physics: no external Krylov/AMG package (GMRES/MG/Schwarz are
# in-tree), no GPU/AD/IO axes, and no GrTool capability may ever enter this core closure.
using Test
using TestKit

ts_hygiene = @testset "hygiene (§7.7 gates)" begin
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "Project.toml");
        allow=["LinearAlgebra", "SparseArrays", "NumericsCore"],
        forbid=["Krylov", "IterativeSolvers", "AlgebraicMultigrid", "CUDA", "AMDGPU",
                "KernelAbstractions", "OrdinaryDiffEq*", "HDF5", "JLD2", "ForwardDiff",
                "SpecialFunctions", "HypergeometricFunctions", "QuadGK",
                # physics never enters core: every GrTool capability is poison here
                "AccretionDisc", "EMRIWaveforms", "Harmonics", "Hybrid",
                "MetricReconstruction", "ModifiedGravity", "Perturbation",
                "PostNewtonian", "SelfForce", "ScalarSelfForce", "Spacetime"])
    # gate 2 (closure-load): runs AFTER `using EllipticCore` (runtests.jl loads it first).
    # Nothing in this test env legitimately loads any name below (TestKit contract).
    TestKit.assert_not_loaded(["Krylov", "IterativeSolvers", "AlgebraicMultigrid",
                               "CUDA", "AMDGPU", "KernelAbstractions", "OrdinaryDiffEq*",
                               "HDF5", "JLD2", "ForwardDiff"])
end
# one-line visibility marker: green nested testsets are folded in the Pkg.test summary.
println("hygiene (§7.7 gates) [EllipticCore]: ",
        (ts_hygiene isa Test.DefaultTestSet && ts_hygiene.anynonpass) ? "FAIL" : "PASS")

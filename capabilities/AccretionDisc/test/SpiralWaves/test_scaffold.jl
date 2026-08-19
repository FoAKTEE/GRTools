# SpiralWaves/test_scaffold.jl — EGD::sw-scaffold (arXiv:2608.17003 spiral-density-wave pipeline).
# Hygiene for the scaffold: the submodule loads and is reachable as `AccretionDisc.SpiralWaves`,
# ships NO physics stubs (its export list is empty; only PAPER + sw_version are defined), and the
# SparseArrays stdlib it will use for the sparse implicit solver is a DECLARED hard dep of the
# package (gate-5 style must_declare through TestKit.assert_dep_budget) and resolvable from the
# package environment. Physics nodes add their own test files as they land.
using Test
using AccretionDisc
using TestKit

@testset "SpiralWaves scaffold: load + dep budget" begin
    @test isdefined(AccretionDisc, :SpiralWaves)
    @test AccretionDisc.SpiralWaves isa Module
    @test :SpiralWaves in names(AccretionDisc)                   # re-exported module name
    SW = AccretionDisc.SpiralWaves
    # scaffold ships no exported symbols (no stubs): names(SW) is just the module itself
    @test names(SW) == [:SpiralWaves]
    @test SW.PAPER == "arXiv:2608.17003"
    @test SW.sw_version() isa VersionNumber
    @test SW.sw_version() == v"0.1.0"
    # no placeholder physics: none of the planned node entry points exist yet
    for sym in (:disc_background, :master_equations, :emri_source, :adaptive_grid, :solve_mode,
                :torque_balance, :survey)
        @test !isdefined(SW, sym)
    end
    # SparseArrays is a declared hard dep (and resolvable from the package env)
    @test Base.identify_package(AccretionDisc, "SparseArrays") !== nothing
    @test isdefined(SW, :SparseArrays) && SW.SparseArrays isa Module
    @test isdefined(SW, :Spacetime)                             # Kerr geometry ONLY via Spacetime
    TestKit.assert_dep_budget(joinpath(@__DIR__, "..", "..", "Project.toml");
        must_declare=["SparseArrays", "Spacetime"],
        forbid=["SelfForce", "Kernels", "Symbolics", "JLD2", "OrdinaryDiffEq*", "CUDA", "AMDGPU",
                "KernelAbstractions", "Roots", "NLsolve", "SciMLBase"])
end

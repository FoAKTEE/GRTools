# SpiralWaves/test_scaffold.jl — EGD::sw-scaffold (arXiv:2608.17003 spiral-density-wave pipeline).
# Hygiene for the scaffold: the submodule loads and is reachable as `AccretionDisc.SpiralWaves`,
# ships NO physics stubs ahead of their tests (every exported name belongs to a LANDED node — the
# landed set is asserted below and grows with the DAG; planned-but-unlanded entry points must NOT
# exist), and the SparseArrays stdlib it will use for the sparse implicit solver is a DECLARED hard
# dep of the package (gate-5 style must_declare through TestKit.assert_dep_budget) and resolvable
# from the package environment. Physics nodes add their own test files as they land.
using Test
using AccretionDisc
using TestKit

@testset "SpiralWaves scaffold: load + dep budget" begin
    @test isdefined(AccretionDisc, :SpiralWaves)
    @test AccretionDisc.SpiralWaves isa Module
    @test :SpiralWaves in names(AccretionDisc)                   # re-exported module name
    SW = AccretionDisc.SpiralWaves
    # every exported symbol belongs to a LANDED node (no stubs): the landed export set is exact
    landed_background = [:DiscBackground, :rotation_profile, :Omega, :omega_d, :Omega_emri, :eps_profile,
                         :pressure, :gamma1, :dlogp_dr_over_enthalpy, :kappa2, :r_isco, :softening_length,
                         :F_J0, :F_J1, :disc_OmegaK, :kappar2, :csound, :eden0]     # EGD::sw-background
    landed_grid = [:SpiralGrid, :build_grid, :corotation_radius, :corotation_block, :smoothmin,
                   :wkb_wavenumber, :wkb_wavenumber2, :wkb_wavenumber_floored, :WKBRegion, :wkb_region,
                   :cell_count]                                                      # EGD::sw-grid
    @test Set(names(SW)) == Set(vcat(:SpiralWaves, landed_background, landed_grid))
    @test SW.PAPER == "arXiv:2608.17003"
    @test SW.sw_version() isa VersionNumber
    @test SW.sw_version() == v"0.3.0"
    # no placeholder physics: none of the UNLANDED planned node entry points exist yet
    for sym in (:master_matrices, :master_matrices_newtonian, :emri_source, :singular_field,
                :solve_mode, :ModeSolution, :torque_balance, :torque_survey)
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

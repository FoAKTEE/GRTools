using Test
using ScalarSelfForce

@testset "ScalarSelfForce" begin
    @testset "smoke: module + spine wiring" begin
        @test isdefined(Main, :ScalarSelfForce)
        @test ScalarSelfForce isa Module
        # [sources]-composition gate: every budget dep resolves and loads through the package
        # (this is the real product of the scaffold node — the relative-path wiring works,
        # including the concurrently created core/EllipticCore).
        for dep in (:LinearAlgebra, :Printf, :NumericsCore, :BigSpecial, :Spacetime,
                    :EllipticCore, :Data)
            @test getfield(ScalarSelfForce, dep) isa Module
        end
    end

    include("hygiene.jl")   # §7.7 binding CI hygiene gates (dep budget + closure-load), before heavy tests

    # Per-module test files are included as each elliptic-ssf GLOBAL_DAG node lands (TDD, one
    # commit per node): Orbits/ → MModeOperator/ → Puncture/ → EffectiveSource/ → Worldtube/ →
    # Domain/ → Solve/ → Extraction/, each grounded against its oracle column (baked
    # EffectiveSource-C fixtures, spectre-verified formulas, WB tables).
    include("test_orbits.jl")   # orbits-adapters (W2-C)
    include("test_effsource_fixtures.jl")   # effsource-oracle-fixtures (W2-D)
    include("test_mode_ansatz.jl")   # mode-ansatz-decomposition (W3-C)
    include("test_puncture_field4d.jl")   # puncture-field-4d (W3-D)
    include("test_operator_tsliced.jl")   # elliptic-operator-tsliced (W4-C)
    include("test_hyperboloidal.jl")   # hyperboloidal-transformation (W5-C)
    include("test_mmode_puncture.jl")   # mmode-puncture-integrals (W4-D)
    include("test_puncture_derivatives.jl")   # puncture-derivatives (W6-C)
    include("test_runtime_config.jl")   # runtime-config-fixture (W4-E)
    include("test_vtu.jl")   # vtu-slicing-jumps (W5-D)
    include("test_wb_fixtures.jl")   # wb-reference-fixtures (W5-F)
    include("test_fluxform.jl")   # horizon-penetrating-compactification (W6-D)
    include("test_first_order_flux_form.jl")   # first-order-flux-form (W6-E)
    include("test_equatorial_symmetry_z2.jl")   # equatorial-symmetry-z2 (W6-F)
    include("test_bayliss_turkel.jl")   # bayliss-turkel-bc (W7-D)
    include("test_effective_source_direct.jl")   # effective-source-direct (W7-C)
    include("test_effective_source_series.jl")   # effective-source-series (W8-B)
    include("test_regularization.jl")   # heffernan-regularization (W7-E)
    include("test_domain_blocks.jl")   # domain-blocks-map (W7-D)
    include("test_worldtube.jl")   # worldtube-formulation (W8-C)
    include("test_bc_glue.jl")   # dg-boundary-conditions (W8-D)
    include("test_jump_conditions.jl")   # dg-jump-conditions (W9-A)
    include("test_assemble.jl")   # fig-matrix commit 1: operator assembly (W9-C)
    include("test_driver.jl")   # solver-stack-ssf (W10-A)
    include("test_amr_driver.jl")   # amr-driver (ESSF-W11-A)
    include("test_selfforce.jl")   # mode-self-force-extraction (ESSF-W12-A)
    include("test_parallel.jl")   # mmode-parallel-driver (ESSF-W12-B)
    include("test_selfforce_levels.jl")   # fig-convergence: the level route (ESSF-W13-B)
    include("test_modesum.jl")   # mode-sum-total-force (ESSF-W13-A)
end

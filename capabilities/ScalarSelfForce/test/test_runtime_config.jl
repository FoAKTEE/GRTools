# runtime-config-fixture (ESSF node, schedule W4-E) — characterization of `crp_paper_config()`
# (src/Domain/config.jl), the Julia literal transcription of the CRP production input file
# anc/ScalarSelfForce.yaml (arXiv:2606.04998; baked byte-exact at
# test/fixtures/runtime_config/ScalarSelfForce.yaml, see PROVENANCE.md there).
#
# Four gates:
#   1. drift pin — sha256 of the baked yaml (and of the ref-paper mirror when the repo
#      checkout is present) must equal the recorded hash; any yaml change forces a review
#      of the transcription and a re-record here.
#   2. spot characterization — ~30 field values asserted against LITERALS duplicated
#      straight from the yaml text at the header-example parameters (yaml lines 9-12:
#      spin=0, orbital_radius=10, m_mode=2, inner_radius=2, outer_radius=1e4,
#      worldtube_radius=4.8, worldtube_dz2=0.25, P=4, L=0) + the jinja derivations at a
#      second parameter point (so the FORMULAS v→t = r0 − 1.5wr, t→u = r0 + 3wr, bounds,
#      L+1 pre-refinement are pinned, not just frozen numbers).
#   3. completeness — a hand-built checklist of EVERY yaml key path (parse-free; built
#      during transcription) maps each path to a struct-field accessor; its LENGTH is
#      asserted and every accessor must resolve on the default config.
#   4. sanity/composition — finiteness, positivity, r^v_* < r0 < r^u_* ordering, m = 0
#      no-outer-BC special case, CircularOrbit isco/radius/horizon consistency, and
#      BigFloat inheritance per the package type policy.

using SHA: sha256
using TestKit: oracle_root, oracle_or_skip

@testset "runtime-config-fixture (Domain/config.jl yaml-twin)" begin

    # ------------------------------------------------------------------------------------
    # 1. drift pin: content hash of the mirrored yaml
    # ------------------------------------------------------------------------------------
    # sha256 of ref-paper/arxiv-2606.04998/src/anc/ScalarSelfForce.yaml (6562 bytes),
    # recorded at transcription time (2026-07-28). If this fails, the yaml mirror or the
    # baked fixture changed: RE-REVIEW every field of crp_paper_config() before re-recording.
    YAML_SHA256 = "480aba9f8d6a7f14541b136cd2d7ba07d0a3d7cb022eba866978972f6a61b45d"

    @testset "yaml content-hash pin" begin
        fixture = joinpath(@__DIR__, "fixtures", "runtime_config", "ScalarSelfForce.yaml")
        @test isfile(fixture)
        @test bytes2hex(sha256(read(fixture))) == YAML_SHA256
        # cross-check the READ-ONLY ref-paper mirror; its root is GRTOOLS_REFPAPER_ROOT and an
        # absent mirror records a Broken (GRTOOLS_STRICT_ORACLES=1 makes it a failure).
        mirror = joinpath(oracle_root("GRTOOLS_REFPAPER_ROOT",
                                      "/data/haiyangw/claude/emri-disk/ref-paper"),
                          "arxiv-2606.04998", "src", "anc", "ScalarSelfForce.yaml")
        if oracle_or_skip(isfile(mirror), "ScalarSelfForce.yaml ref-paper mirror sha256 cross-check";
                          env = "GRTOOLS_REFPAPER_ROOT")
            @test bytes2hex(sha256(read(mirror))) == YAML_SHA256
        end
    end

    # ------------------------------------------------------------------------------------
    # 2a. spot characterization at the yaml header-example parameters (the defaults)
    # ------------------------------------------------------------------------------------
    cfg = crp_paper_config()

    @testset "spot values: Background.CircularOrbit + ScalarSelfForce (yaml 24-49)" begin
        @test cfg.orbit.black_hole_mass == 1.0                 # yaml 26: `BlackHoleMass: &bh_mass 1.`
        @test cfg.orbit.black_hole_spin == 0.0                 # yaml 27 + example `spin=0`
        @test cfg.orbit.orbital_radius == 10.0                 # yaml 28 + example `orbital_radius=10`
        @test cfg.orbit.m_mode_number == 2                     # yaml 29 + example `m_mode=2`
        # yaml 31-34: the jinja expressions VERBATIM (in IEEE Float64, 10.0 − 4.8·1.5 is
        # 2.8000000000000007, NOT 2.8 — the yaml renders and spectre parses that double)
        @test cfg.orbit.hyperboloidal_slicing_transitions ==
              (10.0 - 4.8 * 1.5, 10.0 - 4.8 * 1.5, 10.0 + 4.8 * 3, 10.0 + 4.8 * 3)
        @test cfg.orbit.penetrating_horizon === true           # yaml 35
        @test cfg.orbit.impose_equatorial_symmetry === true    # yaml 36
        # yaml 39-45: v-slicing Block0/Block5, u-slicing Block4/Block9
        @test cfg.null_slicing_blocks == ["Block0", "Block5", "Block4", "Block9"]
        @test cfg.initial_guess_is_background === true         # yaml 47: `InitialGuess: *background`
        @test cfg.randomize_initial_guess === nothing          # yaml 49: `RandomizeInitialGuess: None`
    end

    @testset "spot values: DomainCreator.AlignedLattice (yaml 51-127)" begin
        # yaml 55-61: {inner, v→t, r0−wr, r0+wr, t→u, outer} — jinja expressions verbatim
        @test cfg.domain.block_bounds_r ==
              [2.0, 10.0 - 4.8 * 1.5, 10.0 - 4.8, 10.0 + 4.8, 10.0 + 4.8 * 3, 1.0e4]
        @test cfg.domain.block_bounds_z2 == [0.0, 0.25, 1.0]   # yaml 63: `[0., worldtube_dz2, 1.]`
        # yaml 66: Linear ×4 + Inverse (outer shell); yaml 68: Linear ×2
        @test cfg.domain.distributions_r == [:Linear, :Linear, :Linear, :Linear, :Inverse]
        @test cfg.domain.distributions_z2 == [:Linear, :Linear]
        @test cfg.domain.singularity_positions_r == zeros(5)   # yaml 71
        @test cfg.domain.singularity_positions_z2 == zeros(2)  # yaml 73
        @test cfg.domain.initial_grid_points == (4, 4)         # yaml 74: `[P, P]`, P = 4
        @test cfg.domain.initial_levels == (0, 0)              # yaml 75: `[L, L]`, L = 0
        @test cfg.domain.refined_grid_points == NTuple{2,Int}[]  # yaml 76
        @test cfg.domain.blocks_to_exclude == Int[]            # yaml 77
        # yaml 79-110: 10 pre-refinement entries (middle row 5 + top row 5); source-region
        # column pre-refined to L+1, the asymmetric transition/right entries verbatim.
        @test length(cfg.domain.refined_levels) == 10
        @test cfg.domain.refined_levels[1] == BlockRefinement((0, 0), (1, 1), (0, 0))   # yaml 80-82
        @test cfg.domain.refined_levels[3] == BlockRefinement((2, 0), (3, 1), (1, 1))   # yaml 86-88 source region [L+1, L+1]
        @test cfg.domain.refined_levels[4] == BlockRefinement((3, 0), (4, 1), (0, 1))   # yaml 89-91 [L, L+1]
        @test cfg.domain.refined_levels[5] == BlockRefinement((4, 0), (5, 1), (1, 0))   # yaml 92-94 [L+1, L]
        @test cfg.domain.refined_levels[8] == BlockRefinement((2, 1), (3, 2), (1, 1))   # yaml 102-104
        @test cfg.domain.refined_levels[10] == BlockRefinement((4, 1), (5, 2), (1, 0))  # yaml 108-110
        # yaml 111-127: radial lower None; upper Sommerfeld (m ≠ 0), order 2, hyperboloidal
        @test cfg.domain.lower_radial_bc === nothing
        @test cfg.domain.angular_bc === nothing
        bc = cfg.domain.outer_radial_bc
        @test bc isa SommerfeldBC
        @test bc.black_hole_mass == 1.0 && bc.black_hole_spin == 0.0
        @test bc.orbital_radius == 10.0 && bc.m_mode_number == 2
        @test bc.hyperboloidal_slicing === true                # yaml 123
        @test bc.order == 2                                    # yaml 124
    end

    @testset "spot values: Amr + PhaseChangeAndTriggers (yaml 129-165)" begin
        @test cfg.amr.verbosity == :Quiet                      # yaml 130
        @test cfg.amr.criteria == [:RefineAtPuncture, :IncreaseResolution]  # yaml 131-133
        @test cfg.amr.isotropy == :Anisotropic                 # yaml 145
        @test cfg.amr.enforce_two_to_one_balance_in_normal_direction === true  # yaml 146
        @test cfg.amr.limit_num_grid_points == :auto           # yaml 148: `Auto`
        @test cfg.amr.limit_refinement_level == :auto          # yaml 149: `Auto`
        @test cfg.amr.error_beyond_limits === true             # yaml 150
        @test cfg.amr.allow_coarsening === true                # yaml 151
        @test cfg.amr.max_coarse_levels == :auto               # yaml 152: `Auto`
        @test cfg.amr.iterations == 8                          # yaml 153
        @test cfg.amr.trigger_every_n == 1                     # yaml 159
        @test cfg.amr.trigger_offset == 1                      # yaml 160 (skip the initial guess)
        @test cfg.amr.phase_changes == ["VisitAndReturn(EvaluateAmrCriteria)",
                                        "VisitAndReturn(AdjustDomain)",
                                        "VisitAndReturn(UpdateSections)",
                                        "VisitAndReturn(CheckDomain)"]  # yaml 162-165
    end

    @testset "spot values: Discretization + Observers (yaml 167-176)" begin
        @test cfg.dg.penalty_parameter == 1.5                  # yaml 169: C = 1.5
        @test cfg.dg.massive === true                          # yaml 170
        @test cfg.dg.quadrature == :Gauss                      # yaml 171
        @test cfg.dg.formulation == :StrongInertial            # yaml 172
        @test cfg.observers.volume_file_name == "ScalarSelfForceVolume"        # yaml 175
        @test cfg.observers.reduction_file_name == "ScalarSelfForceReductions" # yaml 176
    end

    @testset "spot values: LinearSolver stack (yaml 178-204)" begin
        @test cfg.gmres.max_iterations == 500                  # yaml 181
        @test cfg.gmres.relative_residual == 0.0               # yaml 182
        @test cfg.gmres.absolute_residual == 1.0e-12           # yaml 183
        @test cfg.gmres.verbosity == :Quiet                    # yaml 184
        @test cfg.multigrid.iterations == 1                    # yaml 187 (1 V-cycle / GMRES iter)
        @test cfg.multigrid.initial_coarse_levels == 0         # yaml 188: `{{ L }}`, L = 0
        @test cfg.multigrid.pre_smoothing === true             # yaml 189
        @test cfg.multigrid.use_bottom_solver === true         # yaml 190
        @test cfg.multigrid.post_smoothing_at_bottom === false # yaml 191
        @test cfg.multigrid.verbosity == :Silent               # yaml 192
        @test cfg.multigrid.output_volume_data === false       # yaml 193
        @test cfg.multigrid.complex_shift == 0.0               # yaml 194 (vtu slicing: no shift)
        @test cfg.schwarz.iterations == 2                      # yaml 197
        @test cfg.schwarz.max_overlap == :auto                 # yaml 198: `Auto`
        @test cfg.schwarz.verbosity == :Silent                 # yaml 199
        @test cfg.schwarz.subdomain_write_matrix_to_file === nothing  # yaml 202: `None`
        @test cfg.schwarz.subdomain_complex_shift == 0.0       # yaml 203
        @test cfg.schwarz.observe_per_core_reductions === false  # yaml 204
        @test cfg.schwarz.output_volume_data === false         # yaml 205
    end

    @testset "spot values: events + extraction + run plumbing (yaml 207-237)" begin
        @test cfg.observers.trigger == :Always                 # yaml 208
        @test cfg.observers.event_self_force == :ObserveSelfForce  # yaml 210
        @test cfg.observers.fields_subfile_name == "VolumeData"    # yaml 212
        @test cfg.observers.variables_to_observe ==
              ["MMode", "SingularField", "RadiallyCompressedCoordinates"]  # yaml 213-216
        @test cfg.observers.blocks_to_observe == :All          # yaml 217
        @test cfg.observers.interpolate_to_mesh === nothing    # yaml 218: `None`
        @test cfg.observers.coordinates_floating_point_type == :Double  # yaml 219
        @test cfg.observers.floating_point_types == [:Double]  # yaml 220
        @test cfg.observers.rcc_inner_radius == 10.0 + 4.8 * 3 # yaml 223: *t_to_u_transition
        @test cfg.observers.rcc_outer_radius == 1.0e4          # yaml 224: *outer_radius
        @test cfg.observers.rcc_compression == :Inverse        # yaml 225: *outer_shell_distribution
        @test cfg.build_matrix_subfile_name === nothing        # yaml 228: `None`
        @test cfg.build_matrix_verbosity == :Silent            # yaml 229
        @test cfg.enable_direct_solve === true                 # yaml 230
        @test cfg.element_distribution == :NumGridPoints       # yaml 233
        @test cfg.avoid_global_proc0 === false                 # yaml 236
        @test cfg.singletons == :auto                          # yaml 237: `Auto`
    end

    # ------------------------------------------------------------------------------------
    # 2b. jinja derivations pinned at a SECOND parameter point (formulas, not frozen numbers)
    # ------------------------------------------------------------------------------------
    @testset "jinja derivations at non-default parameters" begin
        c2 = crp_paper_config(spin=0.5, orbital_radius=8.0, m_mode=3, inner_radius=1.4,
                              outer_radius=500.0, worldtube_radius=2.0, worldtube_dz2=0.125,
                              P=6, L=1)
        # v→t = 8 − 1.5·2 = 5, t→u = 8 + 3·2 = 14 (yaml 31-33)
        @test c2.orbit.hyperboloidal_slicing_transitions == (5.0, 5.0, 14.0, 14.0)
        # bounds {inner, v→t, r0−wr, r0+wr, t→u, outer} (yaml 55-61)
        @test c2.domain.block_bounds_r == [1.4, 5.0, 6.0, 10.0, 14.0, 500.0]
        @test c2.domain.block_bounds_z2 == [0.0, 0.125, 1.0]
        @test c2.domain.initial_grid_points == (6, 6)
        @test c2.domain.initial_levels == (1, 1)
        # source-region pre-refinement to L+1 = 2; asymmetric neighbours verbatim (yaml 86-94)
        @test c2.domain.refined_levels[3].refinement == (2, 2)
        @test c2.domain.refined_levels[4].refinement == (1, 2)
        @test c2.domain.refined_levels[5].refinement == (2, 1)
        @test c2.multigrid.initial_coarse_levels == 1          # yaml 188: `{{ L }}`
        bc2 = c2.domain.outer_radial_bc
        @test bc2.black_hole_spin == 0.5 && bc2.orbital_radius == 8.0 && bc2.m_mode_number == 3
        @test c2.observers.rcc_inner_radius == 14.0 && c2.observers.rcc_outer_radius == 500.0
    end

    @testset "m = 0 outer-BC special case (yaml 115-125 jinja branch)" begin
        @test crp_paper_config(m_mode=0).domain.outer_radial_bc === nothing
        @test crp_paper_config(m_mode=1).domain.outer_radial_bc isa SommerfeldBC
        # paper production m-sweep: 21 separate runs of this yaml (main.tex:127,712)
        @test CRP_PAPER_MMODES == 0:20
    end

    # ------------------------------------------------------------------------------------
    # 3. completeness: every yaml scalar key path is covered by a struct field.
    # Hand-built during transcription (parse-free); list-valued keys are one entry; the 10
    # RefinedLevels entries and the 6 Sommerfeld sub-keys are enumerated individually.
    # ------------------------------------------------------------------------------------
    YAML_KEY_CHECKLIST = Pair{String,Function}[
        # Background.CircularOrbit (yaml 24-36) — 7
        "Background.CircularOrbit.BlackHoleMass"                   => c -> c.orbit.black_hole_mass,
        "Background.CircularOrbit.BlackHoleSpin"                   => c -> c.orbit.black_hole_spin,
        "Background.CircularOrbit.OrbitalRadius"                   => c -> c.orbit.orbital_radius,
        "Background.CircularOrbit.MModeNumber"                     => c -> c.orbit.m_mode_number,
        "Background.CircularOrbit.HyperboloidalSlicingTransitions" => c -> c.orbit.hyperboloidal_slicing_transitions,
        "Background.CircularOrbit.PenetratingHorizon"              => c -> c.orbit.penetrating_horizon,
        "Background.CircularOrbit.ImposeEquatorialSymmetry"        => c -> c.orbit.impose_equatorial_symmetry,
        # ScalarSelfForce (yaml 38-45) — 1
        "ScalarSelfForce.NullSlicingBlocks"                        => c -> c.null_slicing_blocks,
        # top level (yaml 47-49) — 2
        "InitialGuess"                                             => c -> c.initial_guess_is_background,
        "RandomizeInitialGuess"                                    => c -> c.randomize_initial_guess,
        # DomainCreator.AlignedLattice (yaml 51-127) — 29
        "DomainCreator.AlignedLattice.BlockBounds[r]"              => c -> c.domain.block_bounds_r,
        "DomainCreator.AlignedLattice.BlockBounds[z2]"             => c -> c.domain.block_bounds_z2,
        "DomainCreator.AlignedLattice.Distributions[r]"            => c -> c.domain.distributions_r,
        "DomainCreator.AlignedLattice.Distributions[z2]"           => c -> c.domain.distributions_z2,
        "DomainCreator.AlignedLattice.SingularityPositions[r]"     => c -> c.domain.singularity_positions_r,
        "DomainCreator.AlignedLattice.SingularityPositions[z2]"    => c -> c.domain.singularity_positions_z2,
        "DomainCreator.AlignedLattice.InitialGridPoints"           => c -> c.domain.initial_grid_points,
        "DomainCreator.AlignedLattice.InitialLevels"               => c -> c.domain.initial_levels,
        "DomainCreator.AlignedLattice.RefinedGridPoints"           => c -> c.domain.refined_grid_points,
        "DomainCreator.AlignedLattice.BlocksToExclude"             => c -> c.domain.blocks_to_exclude,
        "DomainCreator.AlignedLattice.RefinedLevels[1]"            => c -> c.domain.refined_levels[1],
        "DomainCreator.AlignedLattice.RefinedLevels[2]"            => c -> c.domain.refined_levels[2],
        "DomainCreator.AlignedLattice.RefinedLevels[3]"            => c -> c.domain.refined_levels[3],
        "DomainCreator.AlignedLattice.RefinedLevels[4]"            => c -> c.domain.refined_levels[4],
        "DomainCreator.AlignedLattice.RefinedLevels[5]"            => c -> c.domain.refined_levels[5],
        "DomainCreator.AlignedLattice.RefinedLevels[6]"            => c -> c.domain.refined_levels[6],
        "DomainCreator.AlignedLattice.RefinedLevels[7]"            => c -> c.domain.refined_levels[7],
        "DomainCreator.AlignedLattice.RefinedLevels[8]"            => c -> c.domain.refined_levels[8],
        "DomainCreator.AlignedLattice.RefinedLevels[9]"            => c -> c.domain.refined_levels[9],
        "DomainCreator.AlignedLattice.RefinedLevels[10]"           => c -> c.domain.refined_levels[10],
        "DomainCreator.AlignedLattice.BoundaryConditions[r].Lower" => c -> c.domain.lower_radial_bc,
        "DomainCreator.AlignedLattice.BoundaryConditions[r].Upper" => c -> c.domain.outer_radial_bc,
        "…BoundaryConditions[r].Upper.Sommerfeld.BlackHoleMass"    => c -> c.domain.outer_radial_bc.black_hole_mass,
        "…BoundaryConditions[r].Upper.Sommerfeld.BlackHoleSpin"    => c -> c.domain.outer_radial_bc.black_hole_spin,
        "…BoundaryConditions[r].Upper.Sommerfeld.OrbitalRadius"    => c -> c.domain.outer_radial_bc.orbital_radius,
        "…BoundaryConditions[r].Upper.Sommerfeld.MModeNumber"      => c -> c.domain.outer_radial_bc.m_mode_number,
        "…BoundaryConditions[r].Upper.Sommerfeld.HyperboloidalSlicing" => c -> c.domain.outer_radial_bc.hyperboloidal_slicing,
        "…BoundaryConditions[r].Upper.Sommerfeld.Order"            => c -> c.domain.outer_radial_bc.order,
        "DomainCreator.AlignedLattice.BoundaryConditions[z2]"      => c -> c.domain.angular_bc,
        # Amr (yaml 129-153) — 10
        "Amr.Verbosity"                                            => c -> c.amr.verbosity,
        "Amr.Criteria"                                             => c -> c.amr.criteria,
        "Amr.Policies.Isotropy"                                    => c -> c.amr.isotropy,
        "Amr.Policies.EnforceTwoToOneBalanceInNormalDirection"     => c -> c.amr.enforce_two_to_one_balance_in_normal_direction,
        "Amr.Policies.Limits.NumGridPoints"                        => c -> c.amr.limit_num_grid_points,
        "Amr.Policies.Limits.RefinementLevel"                      => c -> c.amr.limit_refinement_level,
        "Amr.Policies.Limits.ErrorBeyondLimits"                    => c -> c.amr.error_beyond_limits,
        "Amr.Policies.AllowCoarsening"                             => c -> c.amr.allow_coarsening,
        "Amr.MaxCoarseLevels"                                      => c -> c.amr.max_coarse_levels,
        "Amr.Iterations"                                           => c -> c.amr.iterations,
        # PhaseChangeAndTriggers (yaml 155-165) — 3
        "PhaseChangeAndTriggers[1].Trigger.EveryNIterations.N"      => c -> c.amr.trigger_every_n,
        "PhaseChangeAndTriggers[1].Trigger.EveryNIterations.Offset" => c -> c.amr.trigger_offset,
        "PhaseChangeAndTriggers[1].PhaseChanges"                    => c -> c.amr.phase_changes,
        # Discretization.DiscontinuousGalerkin (yaml 167-172) — 4
        "Discretization.DiscontinuousGalerkin.PenaltyParameter"    => c -> c.dg.penalty_parameter,
        "Discretization.DiscontinuousGalerkin.Massive"             => c -> c.dg.massive,
        "Discretization.DiscontinuousGalerkin.Quadrature"          => c -> c.dg.quadrature,
        "Discretization.DiscontinuousGalerkin.Formulation"         => c -> c.dg.formulation,
        # Observers (yaml 174-176) — 2
        "Observers.VolumeFileName"                                 => c -> c.observers.volume_file_name,
        "Observers.ReductionFileName"                              => c -> c.observers.reduction_file_name,
        # LinearSolver.Gmres (yaml 178-184) — 4
        "LinearSolver.Gmres.ConvergenceCriteria.MaxIterations"     => c -> c.gmres.max_iterations,
        "LinearSolver.Gmres.ConvergenceCriteria.RelativeResidual"  => c -> c.gmres.relative_residual,
        "LinearSolver.Gmres.ConvergenceCriteria.AbsoluteResidual"  => c -> c.gmres.absolute_residual,
        "LinearSolver.Gmres.Verbosity"                             => c -> c.gmres.verbosity,
        # LinearSolver.Multigrid (yaml 186-194) — 8
        "LinearSolver.Multigrid.Iterations"                        => c -> c.multigrid.iterations,
        "LinearSolver.Multigrid.InitialCoarseLevels"               => c -> c.multigrid.initial_coarse_levels,
        "LinearSolver.Multigrid.PreSmoothing"                      => c -> c.multigrid.pre_smoothing,
        "LinearSolver.Multigrid.UseBottomSolver"                   => c -> c.multigrid.use_bottom_solver,
        "LinearSolver.Multigrid.PostSmoothingAtBottom"             => c -> c.multigrid.post_smoothing_at_bottom,
        "LinearSolver.Multigrid.Verbosity"                         => c -> c.multigrid.verbosity,
        "LinearSolver.Multigrid.OutputVolumeData"                  => c -> c.multigrid.output_volume_data,
        "LinearSolver.Multigrid.ComplexShift"                      => c -> c.multigrid.complex_shift,
        # LinearSolver.SchwarzSmoother (yaml 196-205) — 7
        "LinearSolver.SchwarzSmoother.Iterations"                  => c -> c.schwarz.iterations,
        "LinearSolver.SchwarzSmoother.MaxOverlap"                  => c -> c.schwarz.max_overlap,
        "LinearSolver.SchwarzSmoother.Verbosity"                   => c -> c.schwarz.verbosity,
        "LinearSolver.SchwarzSmoother.SubdomainSolver.ExplicitInverse.WriteMatrixToFile" => c -> c.schwarz.subdomain_write_matrix_to_file,
        "LinearSolver.SchwarzSmoother.SubdomainSolver.ExplicitInverse.ComplexShift"      => c -> c.schwarz.subdomain_complex_shift,
        "LinearSolver.SchwarzSmoother.ObservePerCoreReductions"    => c -> c.schwarz.observe_per_core_reductions,
        "LinearSolver.SchwarzSmoother.OutputVolumeData"            => c -> c.schwarz.output_volume_data,
        # EventsAndTriggersAtIterations (yaml 207-220) — 8
        "EventsAndTriggersAtIterations[1].Trigger"                 => c -> c.observers.trigger,
        "EventsAndTriggersAtIterations[1].Events[1]"               => c -> c.observers.event_self_force,
        "…Events[2].ObserveFields.SubfileName"                     => c -> c.observers.fields_subfile_name,
        "…Events[2].ObserveFields.VariablesToObserve"              => c -> c.observers.variables_to_observe,
        "…Events[2].ObserveFields.BlocksToObserve"                 => c -> c.observers.blocks_to_observe,
        "…Events[2].ObserveFields.InterpolateToMesh"               => c -> c.observers.interpolate_to_mesh,
        "…Events[2].ObserveFields.CoordinatesFloatingPointType"    => c -> c.observers.coordinates_floating_point_type,
        "…Events[2].ObserveFields.FloatingPointTypes"              => c -> c.observers.floating_point_types,
        # RadiallyCompressedCoordinates (yaml 222-225) — 3
        "RadiallyCompressedCoordinates.InnerRadius"                => c -> c.observers.rcc_inner_radius,
        "RadiallyCompressedCoordinates.OuterRadius"                => c -> c.observers.rcc_outer_radius,
        "RadiallyCompressedCoordinates.Compression"                => c -> c.observers.rcc_compression,
        # BuildMatrix (yaml 227-230) — 3
        "BuildMatrix.MatrixSubfileName"                            => c -> c.build_matrix_subfile_name,
        "BuildMatrix.Verbosity"                                    => c -> c.build_matrix_verbosity,
        "BuildMatrix.EnableDirectSolve"                            => c -> c.enable_direct_solve,
        # Parallelization (yaml 232-233) — 1
        "Parallelization.ElementDistribution"                      => c -> c.element_distribution,
        # ResourceInfo (yaml 235-237) — 2
        "ResourceInfo.AvoidGlobalProc0"                            => c -> c.avoid_global_proc0,
        "ResourceInfo.Singletons"                                  => c -> c.singletons,
    ]

    @testset "completeness: yaml key-path checklist" begin
        # 7+1+2+29+10+3+4+2+4+8+7+8+3+3+1+2 = 94 key paths, hand-counted at transcription.
        @test length(YAML_KEY_CHECKLIST) == 94
        @test allunique(first.(YAML_KEY_CHECKLIST))
        for (path, accessor) in YAML_KEY_CHECKLIST
            @test accessor(cfg) !== missing   # field resolves (a missing field throws → fail)
        end
    end

    # ------------------------------------------------------------------------------------
    # 4. sanity + composition with CircularOrbit + BigFloat inheritance
    # ------------------------------------------------------------------------------------
    @testset "type sanity + domain ordering" begin
        for c in (cfg, crp_paper_config(spin=0.9, orbital_radius=10.0, inner_radius=1.3))
            r0 = c.orbit.orbital_radius
            b = c.domain.block_bounds_r
            @test all(isfinite, b) && issorted(b) && allunique(b)   # strictly increasing
            v_to_t, t_to_u = c.orbit.hyperboloidal_slicing_transitions[1],
                             c.orbit.hyperboloidal_slicing_transitions[3]
            # r^v_* < worldtube ∋ r0 < r^u_* nesting (decomp-crp §2 derived ordering)
            @test b[1] < v_to_t == b[2] < b[3] < r0 < b[4] < t_to_u == b[5] < b[6]
            @test all(isfinite, c.domain.block_bounds_z2)
            @test 0.0 < c.domain.block_bounds_z2[2] < 1.0           # 0 < dz² < 1
            @test c.dg.penalty_parameter > 0
            @test c.gmres.absolute_residual > 0 && c.gmres.relative_residual >= 0
            @test c.gmres.max_iterations > 0 && c.amr.iterations > 0
            @test c.multigrid.iterations > 0 && c.schwarz.iterations > 0
            @test all(p >= 1 for p in c.domain.initial_grid_points)
            @test all(l >= 0 for l in c.domain.initial_levels)
        end
    end

    @testset "composes with CircularOrbit (isco/radius/horizon consistency)" begin
        orb = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)
        @test orb.r_isco == 6.0                       # a = 0 Schwarzschild ISCO
        @test cfg.orbit.orbital_radius > orb.r_isco   # paper orbit r0 = 10 is stable
        # PenetratingHorizon: inner block bound at/inside r₊ (a = 0: r₊ = 2 = inner_radius)
        @test cfg.domain.block_bounds_r[1] <= orb.rp
        # worldtube brackets the particle
        @test cfg.domain.block_bounds_r[3] < orb.r0 < cfg.domain.block_bounds_r[4]
        # spinning production point: same consistency at a = 0.9 (r₊ = 1 + √(1−0.81))
        c9 = crp_paper_config(spin=0.9, orbital_radius=10.0, inner_radius=1.3)
        o9 = CircularOrbit(c9.orbit.black_hole_spin, c9.orbit.orbital_radius)
        @test c9.orbit.orbital_radius > o9.r_isco
        @test c9.domain.block_bounds_r[1] <= o9.rp
    end

    @testset "BigFloat inheritance (package type policy)" begin
        cb = crp_paper_config(spin=big"0.5", orbital_radius=big"10.0", worldtube_radius=big"4.8")
        @test cb isa RunConfig{BigFloat}
        @test cb.orbit.hyperboloidal_slicing_transitions[3] == big"10.0" + big"4.8" * 3
        ob = CircularOrbit(cb.orbit.black_hole_spin, cb.orbit.orbital_radius)
        @test ob isa CircularOrbit{BigFloat}
        @test cb.orbit.orbital_radius > ob.r_isco
    end
end

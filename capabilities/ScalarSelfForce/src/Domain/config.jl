# Domain/config.jl — runtime-config-fixture (ESSF node, schedule W4-E): the typed run
# configuration mirroring the CRP production input file anc/ScalarSelfForce.yaml of
# arXiv:2606.04998 (byte-exact baked copy: test/fixtures/runtime_config/ScalarSelfForce.yaml,
# sha256-pinned by test/test_runtime_config.jl — any yaml change forces a review here).
#
# TRANSCRIPTION, not parsing: `crp_paper_config()` is a Julia LITERAL twin of the yaml — one
# struct field per yaml key (key path in the comment on each field), the Jinja2 parameters of
# the yaml header as keyword arguments (defaults = the header example, yaml:9-12, the exact
# production point behind the paper's figures), and the Jinja2 derivations (v→t = r0 − 1.5wr,
# t→u = r0 + 3wr, block bounds, L+1 source pre-refinement, m = 0 no-outer-BC branch) computed
# exactly as the template does. NO yaml dependency by design (dep budget §7.7).
#
# Value conventions (spectre option syntax → Julia):
#   * yaml `None` ↦ `nothing`;  yaml `Auto` ↦ `:auto`;  enum options ↦ `Symbol`s verbatim
#     (`:Quiet`, `:Silent`, `:Gauss`, `:StrongInertial`, `:Linear`, `:Inverse`, …).
#   * numeric fields are type-generic `{T}` (package type policy: BigFloat inherited, never
#     hard-Float64'd); counts/levels/orders stay `Int`.

"""
    SommerfeldBC(black_hole_mass, black_hole_spin, orbital_radius, m_mode_number,
                 hyperboloidal_slicing, order)

`DomainCreator.AlignedLattice.BoundaryConditions[r].Upper.Sommerfeld` (yaml:118-124): the
order-2 Bayliss–Turkel/Sommerfeld outer radiation condition in hyperboloidal form
(main.tex:403-411). Present only for `m ≠ 0` — the yaml's Jinja2 branch (yaml:115-125)
selects `None` for `m = 0` (`RunConfig.domain.outer_radial_bc === nothing`).
"""
struct SommerfeldBC{T<:AbstractFloat}
    black_hole_mass::T        # …Upper.Sommerfeld.BlackHoleMass    (yaml:119, *bh_mass)
    black_hole_spin::T        # …Upper.Sommerfeld.BlackHoleSpin    (yaml:120, *bh_spin)
    orbital_radius::T         # …Upper.Sommerfeld.OrbitalRadius    (yaml:121, *orbital_radius)
    m_mode_number::Int        # …Upper.Sommerfeld.MModeNumber      (yaml:122, *m_mode_number)
    hyperboloidal_slicing::Bool  # …Upper.Sommerfeld.HyperboloidalSlicing (yaml:123)
    order::Int                # …Upper.Sommerfeld.Order            (yaml:124)
end

"""
    BlockRefinement(lower_corner_index, upper_corner_index, refinement)

One entry of `DomainCreator.AlignedLattice.RefinedLevels` (yaml:78-110): the (r, z²) block
patch `[lower, upper)` in lattice corner indices pre-refined to `refinement` levels.
"""
struct BlockRefinement
    lower_corner_index::NTuple{2,Int}  # …RefinedLevels[i].LowerCornerIndex
    upper_corner_index::NTuple{2,Int}  # …RefinedLevels[i].UpperCornerIndex
    refinement::NTuple{2,Int}          # …RefinedLevels[i].Refinement
end

Base.:(==)(a::BlockRefinement, b::BlockRefinement) =
    a.lower_corner_index == b.lower_corner_index &&
    a.upper_corner_index == b.upper_corner_index && a.refinement == b.refinement

"`Background.CircularOrbit` (yaml:24-36): the physics point — orbit, m-mode, vtu slicing."
struct OrbitConfig{T<:AbstractFloat}
    black_hole_mass::T        # Background.CircularOrbit.BlackHoleMass       (yaml:26, M = 1)
    black_hole_spin::T        # Background.CircularOrbit.BlackHoleSpin       (yaml:27, a)
    orbital_radius::T         # Background.CircularOrbit.OrbitalRadius       (yaml:28, r0)
    m_mode_number::Int        # Background.CircularOrbit.MModeNumber         (yaml:29)
    # yaml:30-34: (v→t, v→t, t→u, t→u) with v→t = r0 − 1.5wr, t→u = r0 + 3wr (Jinja2)
    hyperboloidal_slicing_transitions::NTuple{4,T}  # …HyperboloidalSlicingTransitions
    penetrating_horizon::Bool         # Background.CircularOrbit.PenetratingHorizon (yaml:35)
    impose_equatorial_symmetry::Bool  # …ImposeEquatorialSymmetry (yaml:36, z² variant)
end

"`DomainCreator.AlignedLattice` (yaml:51-127): 5×2 block lattice (r, z²) + boundary conditions."
struct DomainConfig{T<:AbstractFloat}
    # yaml:54-61: radial bounds {inner, v→t, r0−wr, r0+wr, t→u, R}
    block_bounds_r::Vector{T}          # …BlockBounds[r]
    block_bounds_z2::Vector{T}         # …BlockBounds[z²]: [0, dz², 1] (yaml:63)
    distributions_r::Vector{Symbol}    # …Distributions[r]: Linear ×4 + Inverse (yaml:66)
    distributions_z2::Vector{Symbol}   # …Distributions[z²]: Linear ×2 (yaml:68)
    singularity_positions_r::Vector{T} # …SingularityPositions[r] (yaml:71)
    singularity_positions_z2::Vector{T}  # …SingularityPositions[z²] (yaml:73)
    initial_grid_points::NTuple{2,Int} # …InitialGridPoints: [P, P] (yaml:74)
    initial_levels::NTuple{2,Int}      # …InitialLevels: [L, L] (yaml:75)
    refined_grid_points::Vector{NTuple{2,Int}}  # …RefinedGridPoints: [] (yaml:76)
    blocks_to_exclude::Vector{Int}     # …BlocksToExclude: [] (yaml:77)
    refined_levels::Vector{BlockRefinement}  # …RefinedLevels: 10 patches (yaml:78-110)
    lower_radial_bc::Nothing           # …BoundaryConditions[r].Lower: None (yaml:113)
    # …BoundaryConditions[r].Upper (yaml:114-125): Sommerfeld for m ≠ 0, None for m = 0
    outer_radial_bc::Union{SommerfeldBC{T},Nothing}
    angular_bc::Nothing                # …BoundaryConditions[z²]: None (yaml:127)
end

"`Amr` + `PhaseChangeAndTriggers` (yaml:129-165): the hp-AMR driver knobs."
struct AmrConfig
    verbosity::Symbol         # Amr.Verbosity (yaml:130)
    criteria::Vector{Symbol}  # Amr.Criteria (yaml:131-133): RefineAtPuncture + IncreaseResolution
    isotropy::Symbol          # Amr.Policies.Isotropy (yaml:145): Anisotropic
    enforce_two_to_one_balance_in_normal_direction::Bool  # Amr.Policies.Enforce… (yaml:146)
    limit_num_grid_points::Symbol   # Amr.Policies.Limits.NumGridPoints (yaml:148): Auto
    limit_refinement_level::Symbol  # Amr.Policies.Limits.RefinementLevel (yaml:149): Auto
    error_beyond_limits::Bool       # Amr.Policies.Limits.ErrorBeyondLimits (yaml:150)
    allow_coarsening::Bool          # Amr.Policies.AllowCoarsening (yaml:151)
    max_coarse_levels::Symbol       # Amr.MaxCoarseLevels (yaml:152): Auto
    iterations::Int                 # Amr.Iterations (yaml:153): 8
    trigger_every_n::Int      # PhaseChangeAndTriggers[1].Trigger.EveryNIterations.N (yaml:159)
    trigger_offset::Int       # …Trigger.EveryNIterations.Offset (yaml:160): skip initial guess
    phase_changes::Vector{String}  # PhaseChangeAndTriggers[1].PhaseChanges (yaml:161-165)
end

"`Discretization.DiscontinuousGalerkin` (yaml:167-172)."
struct DgConfig{T<:AbstractFloat}
    penalty_parameter::T      # …PenaltyParameter (yaml:169): C = 1.5 (main.tex:595)
    massive::Bool             # …Massive (yaml:170)
    quadrature::Symbol        # …Quadrature (yaml:171): Gauss
    formulation::Symbol       # …Formulation (yaml:172): StrongInertial
end

"`LinearSolver.Gmres` (yaml:178-184): outermost complex-valued Krylov solver."
struct GmresConfig{T<:AbstractFloat}
    max_iterations::Int       # …ConvergenceCriteria.MaxIterations (yaml:181): 500
    relative_residual::T      # …ConvergenceCriteria.RelativeResidual (yaml:182): 0
    absolute_residual::T      # …ConvergenceCriteria.AbsoluteResidual (yaml:183): 1e-12
    verbosity::Symbol         # LinearSolver.Gmres.Verbosity (yaml:184)
end

"`LinearSolver.Multigrid` (yaml:186-194): the V-cycle preconditioner."
struct MultigridConfig{T<:AbstractFloat}
    iterations::Int           # …Iterations (yaml:187): 1 V-cycle per GMRES iteration
    initial_coarse_levels::Int  # …InitialCoarseLevels (yaml:188): L (Jinja2)
    pre_smoothing::Bool       # …PreSmoothing (yaml:189)
    use_bottom_solver::Bool   # …UseBottomSolver (yaml:190): direct coarse solve
    post_smoothing_at_bottom::Bool  # …PostSmoothingAtBottom (yaml:191)
    verbosity::Symbol         # …Verbosity (yaml:192)
    output_volume_data::Bool  # …OutputVolumeData (yaml:193)
    complex_shift::T          # …ComplexShift (yaml:194): 0 under vtu slicing (main.tex:680)
end

"`LinearSolver.SchwarzSmoother` (yaml:196-205): overlapping-Schwarz smoother."
struct SchwarzConfig{T<:AbstractFloat}
    iterations::Int           # …Iterations (yaml:197): 2 smoothing passes
    max_overlap::Symbol       # …MaxOverlap (yaml:198): Auto
    verbosity::Symbol         # …Verbosity (yaml:199)
    subdomain_write_matrix_to_file::Nothing  # …SubdomainSolver.ExplicitInverse.WriteMatrixToFile (yaml:202): None
    subdomain_complex_shift::T  # …SubdomainSolver.ExplicitInverse.ComplexShift (yaml:203): 0
    observe_per_core_reductions::Bool  # …ObservePerCoreReductions (yaml:204)
    output_volume_data::Bool  # …OutputVolumeData (yaml:205)
end

"`Observers` + `EventsAndTriggersAtIterations` + `RadiallyCompressedCoordinates` (yaml:174-225)."
struct ObserverConfig{T<:AbstractFloat}
    volume_file_name::String     # Observers.VolumeFileName (yaml:175)
    reduction_file_name::String  # Observers.ReductionFileName (yaml:176)
    trigger::Symbol              # EventsAndTriggersAtIterations[1].Trigger (yaml:208): Always
    event_self_force::Symbol     # …Events[1] (yaml:210): ObserveSelfForce (F_r^m extraction)
    fields_subfile_name::String  # …Events[2].ObserveFields.SubfileName (yaml:212)
    variables_to_observe::Vector{String}  # …ObserveFields.VariablesToObserve (yaml:213-216)
    blocks_to_observe::Symbol    # …ObserveFields.BlocksToObserve (yaml:217): All
    interpolate_to_mesh::Nothing # …ObserveFields.InterpolateToMesh (yaml:218): None
    coordinates_floating_point_type::Symbol  # …CoordinatesFloatingPointType (yaml:219): Double
    floating_point_types::Vector{Symbol}     # …FloatingPointTypes (yaml:220): [Double]
    rcc_inner_radius::T          # RadiallyCompressedCoordinates.InnerRadius (yaml:223): t→u
    rcc_outer_radius::T          # RadiallyCompressedCoordinates.OuterRadius (yaml:224): R
    rcc_compression::Symbol      # RadiallyCompressedCoordinates.Compression (yaml:225): Inverse
end

"""
    RunConfig{T}

Typed twin of the full anc/ScalarSelfForce.yaml run configuration (arXiv:2606.04998): every
physics/solver parameter of the paper's production runs, one field per yaml key (key paths in
the field comments). Build with [`crp_paper_config`](@ref); one instance = ONE m-mode run
(the paper's production sweep is [`CRP_PAPER_MMODES`](@ref) separate runs).
"""
struct RunConfig{T<:AbstractFloat}
    orbit::OrbitConfig{T}             # Background.CircularOrbit (yaml:24-36)
    null_slicing_blocks::Vector{String}  # ScalarSelfForce.NullSlicingBlocks (yaml:39-45)
    initial_guess_is_background::Bool # InitialGuess: *background (yaml:47)
    randomize_initial_guess::Nothing  # RandomizeInitialGuess: None (yaml:49)
    domain::DomainConfig{T}           # DomainCreator.AlignedLattice (yaml:51-127)
    amr::AmrConfig                    # Amr + PhaseChangeAndTriggers (yaml:129-165)
    dg::DgConfig{T}                   # Discretization.DiscontinuousGalerkin (yaml:167-172)
    observers::ObserverConfig{T}      # Observers + Events… + RadiallyCompressed… (yaml:174-225)
    gmres::GmresConfig{T}             # LinearSolver.Gmres (yaml:178-184)
    multigrid::MultigridConfig{T}     # LinearSolver.Multigrid (yaml:186-194)
    schwarz::SchwarzConfig{T}         # LinearSolver.SchwarzSmoother (yaml:196-205)
    build_matrix_subfile_name::Nothing  # BuildMatrix.MatrixSubfileName: None (yaml:228)
    build_matrix_verbosity::Symbol    # BuildMatrix.Verbosity (yaml:229)
    enable_direct_solve::Bool         # BuildMatrix.EnableDirectSolve (yaml:230)
    element_distribution::Symbol      # Parallelization.ElementDistribution (yaml:233)
    avoid_global_proc0::Bool          # ResourceInfo.AvoidGlobalProc0 (yaml:236)
    singletons::Symbol                # ResourceInfo.Singletons: Auto (yaml:237)
end

"""
    CRP_PAPER_MMODES

The paper's production m-mode sweep, `0:20` — 21 separate runs of the yaml, one
`crp_paper_config(m_mode=m)` each ("20 m-modes", "up to m = 20", main.tex:127, 712;
decomp-crp §2). NOT a yaml key: the yaml holds a single `MModeNumber` per run.
"""
const CRP_PAPER_MMODES = 0:20

"""
    crp_paper_config(; spin=0.0, orbital_radius=10.0, m_mode=2, inner_radius=2.0,
                       outer_radius=1.0e4, worldtube_radius=4.8, worldtube_dz2=0.25,
                       P=4, L=0) -> RunConfig{T}

Literal transcription of anc/ScalarSelfForce.yaml. Keywords = the yaml's Jinja2 template
parameters; defaults = the production example of the yaml header (yaml:9-12). Numeric type
`T = float(promote_type(...))` of the six real parameters (BigFloat inherited); `m_mode`,
`P` (points per element per dimension) and `L` (initial refinement level) are `Int`s.

Jinja2 derivations reproduced exactly (yaml:31-33, 54-61, 79-110, 115-125, 188):
  * `v→t = orbital_radius − worldtube_radius·1.5`, `t→u = orbital_radius + worldtube_radius·3`;
  * radial block bounds `{inner, v→t, r0 − wr, r0 + wr, t→u, outer}`, angular `{0, dz², 1}`;
  * source-region blocks pre-refined to `L + 1` (the 10 RefinedLevels patches verbatim);
  * `m_mode == 0` ⇒ no outer radial BC (`outer_radial_bc === nothing`), else order-2
    hyperboloidal Sommerfeld; `Multigrid.InitialCoarseLevels = L`.
"""
function crp_paper_config(; spin::Real=0.0,            # -p spin            (yaml:10)
                            orbital_radius::Real=10.0, # -p orbital_radius  (yaml:10)
                            m_mode::Integer=2,         # -p m_mode          (yaml:10)
                            inner_radius::Real=2.0,    # -p inner_radius    (yaml:11)
                            outer_radius::Real=1.0e4,  # -p outer_radius    (yaml:11)
                            worldtube_radius::Real=4.8,   # -p worldtube_radius (yaml:12)
                            worldtube_dz2::Real=0.25,     # -p worldtube_dz2    (yaml:12)
                            P::Integer=4,              # -p P               (yaml:12)
                            L::Integer=0,              # -p L               (yaml:12)
                            equatorial_symmetry::Bool=true)  # …ImposeEquatorialSymmetry (yaml:36)
    T = float(promote_type(typeof(spin), typeof(orbital_radius), typeof(inner_radius),
                           typeof(outer_radius), typeof(worldtube_radius),
                           typeof(worldtube_dz2)))
    a, r0 = T(spin), T(orbital_radius)
    wr, dz2 = T(worldtube_radius), T(worldtube_dz2)
    rin, rout = T(inner_radius), T(outer_radius)
    m, Pn, Ln = Int(m_mode), Int(P), Int(L)

    # Jinja2 derivations (yaml:31-33); computed ONCE and aliased exactly like the yaml anchors,
    # with the yaml's own literals (1.5 and 3 are exact in binary FP, so for T = Float64 these
    # reproduce the rendered doubles bit-for-bit, e.g. 10.0 − 4.8·1.5 = 2.8000000000000007)
    v_to_t = r0 - wr * T(1.5)        # &v_to_t_transition: orbital_radius − worldtube_radius·1.5
    t_to_u = r0 + wr * 3             # &t_to_u_transition: orbital_radius + worldtube_radius·3

    orbit = OrbitConfig{T}(
        one(T),                      # Background.CircularOrbit.BlackHoleMass (yaml:26)
        a,                           # Background.CircularOrbit.BlackHoleSpin (yaml:27)
        r0,                          # Background.CircularOrbit.OrbitalRadius (yaml:28)
        m,                           # Background.CircularOrbit.MModeNumber (yaml:29)
        (v_to_t, v_to_t, t_to_u, t_to_u),  # …HyperboloidalSlicingTransitions (yaml:31-34)
        true,                        # …PenetratingHorizon (yaml:35)
        equatorial_symmetry)         # …ImposeEquatorialSymmetry (yaml:36, default True)

    # …BoundaryConditions[r].Upper (yaml:114-125): Jinja2 `{% if m_mode == 0 %} None {% else %}`
    outer_bc = m == 0 ? nothing :
        SommerfeldBC{T}(one(T), a, r0, m, true, 2)   # yaml:118-124

    domain = DomainConfig{T}(
        [rin, v_to_t, r0 - wr, r0 + wr, t_to_u, rout],  # …BlockBounds[r] (yaml:55-61)
        [zero(T), dz2, one(T)],                         # …BlockBounds[z²] (yaml:63)
        [:Linear, :Linear, :Linear, :Linear, :Inverse], # …Distributions[r] (yaml:66)
        [:Linear, :Linear],                             # …Distributions[z²] (yaml:68)
        zeros(T, 5),                                    # …SingularityPositions[r] (yaml:71)
        zeros(T, 2),                                    # …SingularityPositions[z²] (yaml:73)
        (Pn, Pn),                                       # …InitialGridPoints (yaml:74)
        (Ln, Ln),                                       # …InitialLevels (yaml:75)
        NTuple{2,Int}[],                                # …RefinedGridPoints (yaml:76)
        Int[],                                          # …BlocksToExclude (yaml:77)
        BlockRefinement[                                # …RefinedLevels (yaml:78-110)
            BlockRefinement((0, 0), (1, 1), (Ln,     Ln)),      # middle left (yaml:80-82)
            BlockRefinement((1, 0), (2, 1), (Ln,     Ln)),      # middle transition left (yaml:83-85)
            BlockRefinement((2, 0), (3, 1), (Ln + 1, Ln + 1)),  # source region (yaml:86-88)
            BlockRefinement((3, 0), (4, 1), (Ln,     Ln + 1)),  # middle transition right (yaml:89-91)
            BlockRefinement((4, 0), (5, 1), (Ln + 1, Ln)),      # middle right (yaml:92-94)
            BlockRefinement((0, 1), (1, 2), (Ln,     Ln)),      # top left (yaml:96-98)
            BlockRefinement((1, 1), (2, 2), (Ln,     Ln)),      # top transition left (yaml:99-101)
            BlockRefinement((2, 1), (3, 2), (Ln + 1, Ln + 1)),  # top middle (yaml:102-104)
            BlockRefinement((3, 1), (4, 2), (Ln,     Ln + 1)),  # top transition right (yaml:105-107)
            BlockRefinement((4, 1), (5, 2), (Ln + 1, Ln)),      # top right (yaml:108-110)
        ],
        nothing,                     # …BoundaryConditions[r].Lower: None (yaml:113)
        outer_bc,                    # …BoundaryConditions[r].Upper (yaml:114-125)
        nothing)                     # …BoundaryConditions[z²]: None (yaml:127)

    amr = AmrConfig(
        :Quiet,                                    # Amr.Verbosity (yaml:130)
        [:RefineAtPuncture, :IncreaseResolution],  # Amr.Criteria (yaml:131-133)
        :Anisotropic,                              # Amr.Policies.Isotropy (yaml:145)
        true,    # Amr.Policies.EnforceTwoToOneBalanceInNormalDirection (yaml:146)
        :auto,   # Amr.Policies.Limits.NumGridPoints: Auto (yaml:148)
        :auto,   # Amr.Policies.Limits.RefinementLevel: Auto (yaml:149)
        true,    # Amr.Policies.Limits.ErrorBeyondLimits (yaml:150)
        true,    # Amr.Policies.AllowCoarsening (yaml:151)
        :auto,   # Amr.MaxCoarseLevels: Auto (yaml:152)
        8,       # Amr.Iterations (yaml:153)
        1,       # PhaseChangeAndTriggers[1].Trigger.EveryNIterations.N (yaml:159)
        1,       # …Trigger.EveryNIterations.Offset (yaml:160)
        ["VisitAndReturn(EvaluateAmrCriteria)",    # …PhaseChanges (yaml:161-165)
         "VisitAndReturn(AdjustDomain)",
         "VisitAndReturn(UpdateSections)",
         "VisitAndReturn(CheckDomain)"])

    dg = DgConfig{T}(
        T(1.5),          # Discretization.DiscontinuousGalerkin.PenaltyParameter: 1.5 (yaml:169)
        true,            # …Massive (yaml:170)
        :Gauss,          # …Quadrature (yaml:171)
        :StrongInertial) # …Formulation (yaml:172)

    observers = ObserverConfig{T}(
        "ScalarSelfForceVolume",      # Observers.VolumeFileName (yaml:175)
        "ScalarSelfForceReductions",  # Observers.ReductionFileName (yaml:176)
        :Always,                      # EventsAndTriggersAtIterations[1].Trigger (yaml:208)
        :ObserveSelfForce,            # …Events[1] (yaml:210)
        "VolumeData",                 # …Events[2].ObserveFields.SubfileName (yaml:212)
        ["MMode", "SingularField", "RadiallyCompressedCoordinates"],  # …VariablesToObserve (yaml:213-216)
        :All,                         # …ObserveFields.BlocksToObserve (yaml:217)
        nothing,                      # …ObserveFields.InterpolateToMesh: None (yaml:218)
        :Double,                      # …CoordinatesFloatingPointType (yaml:219)
        [:Double],                    # …FloatingPointTypes (yaml:220)
        t_to_u,                       # RadiallyCompressedCoordinates.InnerRadius (yaml:223, *t_to_u_transition)
        rout,                         # RadiallyCompressedCoordinates.OuterRadius (yaml:224, *outer_radius)
        :Inverse)                     # RadiallyCompressedCoordinates.Compression (yaml:225, *outer_shell_distribution)

    gmres = GmresConfig{T}(
        500,          # LinearSolver.Gmres.ConvergenceCriteria.MaxIterations (yaml:181)
        zero(T),      # …ConvergenceCriteria.RelativeResidual (yaml:182)
        T(1.0e-12),   # …ConvergenceCriteria.AbsoluteResidual (yaml:183): the DOUBLE literal
                      # `1.e-12` exactly as spectre parses it (10^-12 is not a double)
        :Quiet)       # LinearSolver.Gmres.Verbosity (yaml:184)

    multigrid = MultigridConfig{T}(
        1,          # LinearSolver.Multigrid.Iterations (yaml:187)
        Ln,         # …InitialCoarseLevels: L (yaml:188)
        true,       # …PreSmoothing (yaml:189)
        true,       # …UseBottomSolver (yaml:190)
        false,      # …PostSmoothingAtBottom (yaml:191)
        :Silent,    # …Verbosity (yaml:192)
        false,      # …OutputVolumeData (yaml:193)
        zero(T))    # …ComplexShift (yaml:194)

    schwarz = SchwarzConfig{T}(
        2,          # LinearSolver.SchwarzSmoother.Iterations (yaml:197)
        :auto,      # …MaxOverlap: Auto (yaml:198)
        :Silent,    # …Verbosity (yaml:199)
        nothing,    # …SubdomainSolver.ExplicitInverse.WriteMatrixToFile: None (yaml:202)
        zero(T),    # …SubdomainSolver.ExplicitInverse.ComplexShift (yaml:203)
        false,      # …ObservePerCoreReductions (yaml:204)
        false)      # …OutputVolumeData (yaml:205)

    RunConfig{T}(
        orbit,
        ["Block0", "Block5", "Block4", "Block9"],  # ScalarSelfForce.NullSlicingBlocks (yaml:39-45)
        true,          # InitialGuess: *background (yaml:47)
        nothing,       # RandomizeInitialGuess: None (yaml:49)
        domain,
        amr,
        dg,
        observers,
        gmres,
        multigrid,
        schwarz,
        nothing,       # BuildMatrix.MatrixSubfileName: None (yaml:228)
        :Silent,       # BuildMatrix.Verbosity (yaml:229)
        true,          # BuildMatrix.EnableDirectSolve (yaml:230)
        :NumGridPoints,  # Parallelization.ElementDistribution (yaml:233)
        false,         # ResourceInfo.AvoidGlobalProc0 (yaml:236)
        :auto)         # ResourceInfo.Singletons: Auto (yaml:237)
end

# ===================== the observer's radially compressed coordinate =========================

"""
    radially_compressed(obs::ObserverConfig, r) -> Real

The `RadiallyCompressedCoordinates` radius of `obs` at the physical radius `r` (yaml:222-225)
— the coordinate the paper's own volume observer writes so a domain that runs to
`R = 10^4 M` can be drawn at all (`VariablesToObserve` includes
`RadiallyCompressedCoordinates`, yaml:213-216), and therefore the abscissa of
arXiv:2606.04998 Fig. 1.

TRANSCRIBED from spectre, not re-derived: `src/Domain/RadiallyCompressedCoordinates.cpp:27-52`
leaves `r ≤ r_in` untouched, sets the *compressed* outer radius to `r_in · log10(R)` ("chosen
so it increases with both the inner and the outer radius but exponentials are tamed"), and
returns the INVERSE of the `Interval` map that distributes the grid points radially,

    Interval{A = r_in, B = r_in·log10(R), a = r_in, b = R, distribution, singularity_pos = 0}

so this function is `Interval::inverse` (`src/Domain/CoordinateMaps/Interval.cpp:123-131`) for
`Distribution::Inverse`, the distribution the production yaml selects
(`Compression: *outer_shell_distribution`, the same `Inverse` as the wavezone's own radial
map):

    r_c(r) = ½ (A + B + (B − A)/(b − a) · (a + b − 2ab/r))    for r > r_in,   r_c(r) = r  below.

At the production point (`r_in = 24.4`, `R = 10^4`, hence `B = 97.6`) this is
`r_c(r) = 97.77904 − 1790.4487/r`, carrying `[24.4, 10^4]` onto `[24.4, 97.6]`.

Only `Inverse` is implemented; any other `rcc_compression` throws rather than silently
drawing a different axis. Type-generic in the config's own real type.
See [`radially_compressed_inv`](@ref) for the return trip.
"""
function radially_compressed(obs::ObserverConfig{T}, r::Real) where {T}
    obs.rcc_compression === :Inverse || throw(ArgumentError(
        "radially_compressed: only spectre's `Inverse` compression is implemented, got " *
        "$(obs.rcc_compression) (RadiallyCompressedCoordinates.Compression, yaml:225)"))
    A = obs.rcc_inner_radius
    b = obs.rcc_outer_radius
    rr = T(r)
    rr <= A && return rr
    B = A * log10(b)
    return (A + B + (B - A) / (b - A) * (A + b - 2 * A * b / rr)) / 2
end

"""
    radially_compressed_inv(obs::ObserverConfig, rc) -> Real

The physical radius whose [`radially_compressed`](@ref) coordinate is `rc` — spectre's
`Interval::operator()` for `Distribution::Inverse`
(`src/Domain/CoordinateMaps/Interval.cpp:81-89`) with the same
`{A = r_in, B = r_in·log10(R), a = r_in, b = R}` and singularity at the origin:

    r(r_c) = 2ab / (a + b − (b − a)/(B − A) · (2 r_c − B − A))   for r_c > r_in.

Below `r_in` the map is the identity, so `radially_compressed_inv(obs, radially_compressed(obs, r)) == r`
to round-off on the whole domain. It exists so a caller can lay points out UNIFORMLY in the
drawn coordinate — which is what makes a probe set on a domain spanning four decades in `r`
sample the picture rather than its outer decade.
"""
function radially_compressed_inv(obs::ObserverConfig{T}, rc::Real) where {T}
    obs.rcc_compression === :Inverse || throw(ArgumentError(
        "radially_compressed_inv: only spectre's `Inverse` compression is implemented, got " *
        "$(obs.rcc_compression) (RadiallyCompressedCoordinates.Compression, yaml:225)"))
    A = obs.rcc_inner_radius
    b = obs.rcc_outer_radius
    x = T(rc)
    x <= A && return x
    B = A * log10(b)
    return 2 * A * b / (A + b - (b - A) / (B - A) * (2 * x - B - A))
end

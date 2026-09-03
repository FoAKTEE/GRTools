# Solve/assemble.jl — operator assembly (ESSF node `fig-matrix`, GLOBAL_DAG row 49, commit 1;
# the ASSEMBLY half of the row-42 `solver-stack-ssf` contract, landed once here so the solver
# driver, the AMR driver and the reproduction figures all call the SAME composition).
#
# Everything the CRP DG operator of arXiv:2606.04998 needs is already sealed by the
# predecessor nodes; what was missing is the ONE place that composes them:
#
#   coefficients  Domain/blocks.jl   `domain_coefficients(dom, orbit, config)`        (row 39)
#   hooks         Solve/jump_glue.jl `jump_hooks(dom, wt, orbit, config, m; mesh)`    (row 40)
#   bcs           Solve/bc_glue.jl   `boundary_conditions(dom, orbit, config; mesh)`  (row 41)
#   penalty       Domain/config.jl   `config.dg.penalty_parameter` = C = 1.5 (yaml:169)
#   operator      EllipticCore       `MeshOperator(T, mesh, A, β, γ; C, hooks, bcs)`  (row 15)
#
# [`assemble_operator`](@ref) is exactly that composition and NOTHING else: zero new physics,
# zero re-derived coefficients, no solve. Its value is that the composition exists once, is
# gated BITWISE against the manual spelling, and is callable on ANY mesh built from the block
# lattice — which is what the multigrid `build_op` factory ([`operator_builder`](@ref)) and
# the hp-AMR driver need.
#
# WHAT THIS FILE DELIBERATELY DOES NOT OWN. A [`WorldtubeData`](@ref) is MESH-BOUND (its jump
# data carry one `FaceState` per mortar of one mesh and its fixed source is a nodal vector in
# that mesh's layout), so a level walk needs one per level — but REBUILDING it is row 38's
# `worldtube_data(dom, orbit, config, m; mesh = …)`, not this node's. Hence
# [`operator_builder`](@ref) and [`coarsest_operator`](@ref) take the per-level bundle as a
# CALLER-SUPPLIED factory `worldtube(mesh) -> WorldtubeData`; the CRP one is written in the
# `operator_builder` docstring and in `jump_glue.jl`'s per-level recipe. No second copy of the
# worldtube assembly is created here.
#
# THE COARSEST GRID LEVEL (what the paper's Fig. 2 plots). `LinearSolver.Multigrid.
# InitialCoarseLevels: L` (yaml:188) counts the COARSER levels only, so the production run's
# multigrid hierarchy has `L + 1` levels in total and its COARSEST level is the initial mesh
# h-coarsened `L` times — at the production `L = 0` that is the initial mesh itself (20
# elements, 320 dofs; main.tex:622-626 and the figure's element numbering 1…20).
#
# SCOPE: assembly only. No GMRES, no multigrid V-cycle, no Schwarz smoother, no AMR — nodes
# `solver-stack-ssf` and `amr-driver` own those and consume this file.

export assemble_operator, operator_builder, coarsest_operator
export linearized_rhs

"""
    assemble_operator(T, dom::BlockDomain, wt::WorldtubeData, orbit::CircularOrbit,
                      config::RunConfig, m = config.orbit.m_mode_number;
                      mesh = dom.mesh) -> MeshOperator

The CRP production DG operator: the sealed coefficient bundle, hook table and boundary-
condition table of arXiv:2606.04998 composed into one EllipticCore [`MeshOperator`](@ref).
Identical, BITWISE, to writing

    fcs = domain_coefficients(dom, orbit, config)
    MeshOperator(T, mesh, fcs.A, fcs.beta, fcs.gamma;
                 C     = config.dg.penalty_parameter,
                 hooks = jump_hooks(dom, wt, orbit, config, m; mesh = mesh),
                 bcs   = boundary_conditions(dom, orbit, config; mesh = mesh))

by hand (gated on `build_matrix`, `operator_on_zero` and `apply!`), which is the point: the
composition is written ONCE and every consumer — the single-mode driver, the multigrid
`build_op` factory, the AMR driver, the reproduction figures — calls this.

`T` is the requested solver scalar type (`ComplexF64` on the production path; the operator
scalar promotes it with the evaluated coefficients, so even a real `T` yields a complex
operator here because `β` carries `imΩ`). `m` must equal the config's `MModeNumber` and the
worldtube data's mode — the sealed guards say so, and a mismatch would silently solve a
different equation.

`mesh` may be ANY mesh built from `dom`'s block lattice: `dom.mesh`, a
[`coarsen_mesh`](@ref) level, or an AMR-refined one. The coefficient bundle is
mesh-independent, the boundary and junction tables are recomputed for `mesh` from the block
lattice (block ids survive coarsening and AMR while element ids do not), and `wt` MUST already
belong to `mesh` — rebuild it with `worldtube_data(dom, orbit, config, m; mesh = mesh)`, never
pass a stale one (the sealed `jump_hooks` guard refuses it rather than mis-applying it).

The linearized [`apply!`](@ref) of the result is the GMRES matrix — it carries the vtu junction
hooks, whose contribution is the ONLY jump term in `build_matrix` — while the worldtube jump
data are trial-field independent and reach the right-hand side through `operator_on_zero`,
which is what [`linearized_rhs`](@ref) folds.
"""
function assemble_operator(::Type{T}, dom::BlockDomain, wt::WorldtubeData,
                           orbit::CircularOrbit, config::RunConfig,
                           m::Integer = config.orbit.m_mode_number;
                           mesh::Mesh2D = dom.mesh) where {T<:Number}
    fcs = domain_coefficients(dom, orbit, config)
    return MeshOperator(T, mesh, fcs.A, fcs.beta, fcs.gamma;
                        C = config.dg.penalty_parameter,
                        hooks = jump_hooks(dom, wt, orbit, config, m; mesh = mesh),
                        bcs = boundary_conditions(dom, orbit, config; mesh = mesh))
end

"""
    operator_builder(T, dom, orbit, config, m = config.orbit.m_mode_number; worldtube)
        -> (mesh::Mesh2D -> MeshOperator)

The per-level operator FACTORY the sealed [`MultigridHierarchy`](@ref) takes as `build_op`:
`mesh -> assemble_operator(T, dom, worldtube(mesh), orbit, config, m; mesh = mesh)`.

It rediscretizes the SAME problem on whatever mesh it is handed — same coefficient closures,
same penalty `C`, boundary conditions and junction hooks rebuilt for that mesh's own faces —
which is precisely the factory contract `MultigridHierarchy` documents. Use it for every
hierarchy rebuild (multigrid levels, AMR levels) so the level operators cannot drift from the
fine one.

`worldtube(mesh) -> WorldtubeData` is REQUIRED and is the caller's, because a `WorldtubeData`
is MESH-BOUND and rebuilding it belongs to the worldtube node, not to this one. The CRP
factory is

    worldtube = mesh -> worldtube_data(dom, orbit, config, m; mesh = mesh)

(`jump_glue.jl`'s per-level recipe). A caller with an already-built bundle for a single mesh
passes `_ -> wt`; anything else is refused downstream by the sealed mesh-boundness guards
rather than mis-applied.
"""
function operator_builder(::Type{T}, dom::BlockDomain, orbit::CircularOrbit,
                          config::RunConfig, m::Integer = config.orbit.m_mode_number;
                          worldtube) where {T<:Number}
    return (mesh::Mesh2D) -> assemble_operator(T, dom, worldtube(mesh)::WorldtubeData,
                                               orbit, config, m; mesh = mesh)
end

"""
    coarsest_operator(T, dom, orbit, config, m = config.orbit.m_mode_number; worldtube,
                      maxlevels = config.multigrid.initial_coarse_levels + 1,
                      mesh = dom.mesh)
        -> (op = MeshOperator, mesh = Mesh2D, wt = WorldtubeData, nlevels = Int)

The operator at the COARSEST multigrid level — what arXiv:2606.04998 Fig. 2 plots
(main.tex:622-626, "matrix representation of the DG operator at the coarsest grid level").

`mesh` is h-coarsened with the sealed [`coarsen_mesh`](@ref) until either `maxlevels` levels
have been produced or no further progress is possible (one element per block, the irreducible
block structure of main.tex:672), and the operator is re-assembled on the result through
[`operator_builder`](@ref). `nlevels` reports how many levels that took, the finest included.

The default `maxlevels` is the yaml twin of `LinearSolver.Multigrid.InitialCoarseLevels: L`
(yaml:188), which counts the COARSER levels only — the same `L + 1` convention
[`MultigridHierarchy`](@ref) documents. At the production point (`L = 0`) it is therefore
`1`: the hierarchy has a single level and the "coarsest grid level" IS the initial mesh, the
20-element / 320-dof lattice the paper's figure numbers 1…20 (5×2 blocks with the source
column pre-refined to `(L+1, L+1)`, the `r₀+w`…`r_u` column split angularly and the
`r_u`…`R` column split radially, yaml:78-110). Pass `maxlevels = nothing` to coarsen all the
way down instead: 10 elements / 160 dofs on the 5×2 lattice, 15 / 240 on the 5×3 one.

The level meshes are exactly [`MultigridHierarchy`](@ref)'s at its default `shrink = 1` (one
coarsening step per level) and the returned operator is bitwise its `levels[end].op` — gated
in `test_assemble.jl`. This function exists so that a caller who only wants the coarse
OPERATOR (the figure, a bottom-solve study) does not pay for the per-level Schwarz
factorizations a full hierarchy builds.
"""
function coarsest_operator(::Type{T}, dom::BlockDomain, orbit::CircularOrbit,
                           config::RunConfig, m::Integer = config.orbit.m_mode_number;
                           worldtube,
                           maxlevels::Union{Nothing,Integer} =
                               config.multigrid.initial_coarse_levels + 1,
                           mesh::Mesh2D = dom.mesh) where {T<:Number}
    maxlevels === nothing || Int(maxlevels) >= 1 || throw(ArgumentError(
        "coarsest_operator: maxlevels must be ≥ 1 — it counts the TOTAL number of multigrid " *
        "levels, the finest included, while the yaml's InitialCoarseLevels counts only the " *
        "coarser ones (its twin is L + 1); got $maxlevels"))
    cur = mesh
    nlevels = 1
    while maxlevels === nothing || nlevels < Int(maxlevels)
        co = coarsen_mesh(cur)
        co === nothing && break
        cur = co[1]
        nlevels += 1
    end
    wtc = worldtube(cur)::WorldtubeData
    return (op = assemble_operator(T, dom, wtc, orbit, config, m; mesh = cur),
            mesh = cur, wt = wtc, nlevels = nlevels)
end

# ------------------------------------------------------------------------------------------
# The right-hand side. EllipticCore's `solver-stack-core` node owns the GENERAL fold
# `linearized_rhs(op, b_phys) = b_phys − operator_on_zero(op)`; when that method is present in
# the loaded EllipticCore this file EXTENDS it (one function, two methods) instead of defining
# a second function of the same name — which `using EllipticCore` would make an error anyway.
# Until it lands, the general method is defined here so the CRP driver has it. Either way the
# arithmetic is written exactly ONCE, in one place, and both spellings are gated.
if isdefined(EllipticCore, :linearized_rhs)
    import EllipticCore: linearized_rhs
else
    linearized_rhs(op::MeshOperator, b_phys::AbstractVector) = b_phys .- operator_on_zero(op)
end

"""
    linearized_rhs(op::MeshOperator, wt::WorldtubeData) -> Vector

The right-hand side of the CRP linear system: the worldtube fixed source folded with the
operator's fixed-data lift,

    b = worldtube_massive_source(wt) − operator_on_zero(op) ,

so that the LINEARIZED [`apply!`](@ref) of `op` is the matrix (main.tex:596-608; spectre
`DgOperator.hpp:1081-1130` harvests the same modification on zero data into the fixed
sources). `op` must be the operator assembled on `wt`'s mesh — [`assemble_operator`](@ref)
with the same `wt` — or the two vectors belong to different discretizations.

Skipping the `operator_on_zero` term is an O(1) wrong answer, not a small one: on the
production lattice it is of the same order as the source itself, because it IS the worldtube
jump data entering through the numerical flux.

NOTE for node `solver-stack-ssf`: the VOLUME source inside the tube that the production driver
hands to `worldtube_data` is the SWITCHED evaluator (`effective_source_m_switched`); this
function only folds whatever fixed source the bundle carries with the operator's own lift.
"""
linearized_rhs(op::MeshOperator, wt::WorldtubeData) =
    linearized_rhs(op, worldtube_massive_source(wt))

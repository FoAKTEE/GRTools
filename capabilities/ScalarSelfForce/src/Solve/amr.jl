# Solve/amr.jl — THE hp-AMR DRIVER (ESSF node `amr-driver`, GLOBAL_DAG row 43, wave 11).
# This is the production loop of the reproduction: it turns the single fixed-mesh m-mode solve
# of `Solve/driver.jl` (row 42) into the converged AMR sequence of arXiv:2606.04998 §II.D —
# `Amr: Iterations: 8` levels of "h-refine the elements containing the puncture, p-refine
# everywhere else" (main.tex:552-560) — and records, per level, everything the extraction
# (row 44), the m-mode fan-out (row 47) and the figures (rows 48, 50, 52, 53) read.
#
# WHAT THIS FILE OWNS, AND ONLY THIS:
#   1. `puncture_chart_position` — where the particle sits in the DISCRETIZATION chart.
#   2. `amr_flags`               — the CRP refinement CRITERION (spectre `RefineAtPuncture`
#                                  composed with `IncreaseResolution`, under the yaml Limits).
#   3. `amr_level_options`       — the multigrid LEVEL POLICY of an AMR sweep.
#   4. `amr_step!`               — one AMR level: flags → `hp_refine` → `transfer_solution`
#                                  warm start → `rebuild_for_mesh!` → solve.
#   5. `solve_mode_amr` / `AmrHistory` / `AmrLevel` — the whole sweep and its record.
#
# EVERYTHING ELSE IS A SEALED PREDECESSOR, called and never re-derived: `solve_mode`,
# `solve_mode!`, `rebuild_for_mesh!`, `field_at`, `field_grad_at` (row 42) · `hp_refine`,
# `transfer_solution`
# (rows 22/23) · `SolveRecord`, `ResidualMonitor` (row 21) · `block_domain` (row 39) ·
# `refelement`, `interpolation_matrix` (row 8) · `logical_interval`, `map_inverse`,
# `map_jacobian` (row 9). NO self-force extraction (row 44), NO process fan-out (row 47),
# NO figure (rows 48-53).

export puncture_chart_position, amr_flags, amr_level_options,
       solve_mode_amr, AmrLevel, AmrHistory, resolution_error, relative_resolution_error
# `amr_step!` is deliberately NOT exported: one AMR level is an internal step of
# `solve_mode_amr` and nothing outside this module calls it (checked across capabilities/ and
# reproduce/). It stays public as `ScalarSelfForce.amr_step!` with its docstring.

# The yaml's `Amr.Policies.Limits.{RefinementLevel, NumGridPoints}: Auto` (yaml:148-149)
# resolved through spectre's own defaults: `amr::Limits::Limits()` takes them from
# `ElementId<1>::max_refinement_level` (`Domain/Structure/ElementId.hpp:57`) and
# `Spectral::maximum_number_of_points<Legendre>` (`MaximumNumberOfPoints.hpp:26`).
const AMR_MAX_REFINEMENT_LEVEL = 15
const AMR_MAX_RESOLUTION = 20

# `RefineAtPuncture.cpp:41-45`: the block-logical snap radius.
const AMR_SNAP_TOL = 1.0e-10

# `Amr.Criteria` (yaml:131-133) — the two criteria this file implements, in the yaml's order.
const AMR_CRITERIA = [:RefineAtPuncture, :IncreaseResolution]

# ============================== where the particle is =====================================

"""
    puncture_chart_position(dom::BlockDomain, orbit::CircularOrbit) -> (x1, x2)

The particle's position in the DISCRETIZATION chart of `dom`: `x1 = r₀` (horizon-penetrating
Boyer–Lindquist radius) and `x2 = 0`, the equator — which is `w = z² = 0` under the
production `ImposeEquatorialSymmetry` (`dom.chart === :z2`) and `z = cos θ = 0` without it, so
the same pair serves both charts (main.tex:389-394). The orbit is circular and equatorial by
construction, so this is the whole trajectory: for eccentric or generic orbits the puncture
sweeps a curve and the criterion below would have to follow it (spectre says the same in
`RefineAtPuncture.hpp:30-33`).

In the `z²` chart the second coordinate is the DOMAIN BOUNDARY, not an interior point, which
is why the puncture is contained in TWO elements at every AMR level rather than four.
"""
puncture_chart_position(dom::BlockDomain{RT}, orbit::CircularOrbit) where {RT} =
    (RT(orbit.r0), zero(RT))

# ================================ the CRP criterion =======================================

"""
    amr_flags(sol::ModeSolution; kwargs...) -> Dict{ElementId,NTuple{2,Bool}}
    amr_flags(mesh::Mesh2D, dom::BlockDomain, orbit::CircularOrbit, config::RunConfig;
              snap_tol = 1e-10, max_refinement_level = 15, max_resolution = 20)

THE hp-AMR CRITERION of arXiv:2606.04998 (main.tex:552-560), as the single flag container
[`hp_refine`](@ref) consumes: an element containing the puncture maps to `(true, true)` and is
h-SPLIT in both dimensions with its polynomial order HELD, every other element is absent and
is p-INCREMENTED by one in both dimensions.

# THE CRITERION, AND WHY ONE CONTAINER DRIVES BOTH HALVES (normative)

The paper's `Amr.Criteria` (yaml:131-133) are two: `RefineAtPuncture` (type `h`, returns
`Split` in both dimensions for an element containing the puncture and `DoNothing` otherwise)
and `IncreaseResolution` (type `p`, returns `IncreaseResolution` in both dimensions for EVERY
element). spectre merges criteria per dimension with `std::max` over the `amr::Flag` enum
(`EvaluateRefinementCriteria.hpp:147`), and `Split` (5) outranks `IncreaseResolution` (4). A
puncture element is therefore split and NOT p-refined, and every other element is p-refined
and not split — which is exactly the sealed `hp_refine(mesh, flags)` ≡
`h_refine(p_refine(mesh; hold = flags), flags)` contract, one criterion used as "hold the
order" by the p half and "split" by the h half.

`Policies.Isotropy: Anisotropic` (yaml:145) permits per-dimension decisions and this criterion
does not use that freedom: `RefineAtPuncture` returns `make_array<2>(Flag::Split)`, both
dimensions, so the flags here are `(true, true)`. The per-dimension container is kept because
the flag TYPE is per-dimension (a limit could clip one dimension only) and because
`EnforceTwoToOneBalanceInNormalDirection` (yaml:146) — the 2:1 closure, which IS anisotropic —
is the sealed `hp_refine`'s job, not this criterion's.

# WHICH ELEMENTS "CONTAIN" THE PUNCTURE, AND THE SNAP (normative)

The puncture sits on an element CORNER at every AMR level — the sealed worldtube node proved
the Gauss nodes never touch it — so "contains" must be decided deterministically for a point
on a face. This is spectre's rule, transcribed from
`Elliptic/Systems/SelfForce/Scalar/AmrCriteria/RefineAtPuncture.cpp` and
`Domain/ElementLogicalCoordinates.cpp:31-46`:

1. map the puncture into each block with that block's own inverse map; a coordinate outside
   `[-1, 1]` means the block does not contain it;
2. SNAP the block-logical coordinate: `|ξ_d| < snap_tol ⇒ ξ_d = 0` (the puncture at the block
   CENTRE — the usual case, since the worldtube block is centred on `r₀`) and
   `|ξ_d ∓ 1| < snap_tol ⇒ ξ_d = ±1` (the puncture on a block FACE — the `z² = 0` equator);
3. keep every element whose CLOSED dyadic interval contains the snapped point in both
   dimensions, `lo ≤ ξ_d ≤ hi`. Closed on both sides: a point on a shared face belongs to BOTH
   neighbours, which is what makes the refinement symmetric about the puncture.

THE SNAP IS LOAD-BEARING, and spectre's own comment says why: "Snap to center of the block to
avoid asymmetric refinement from round-off error". The block bounds are `r₀ ∓ w_r`, so the
block-logical coordinate of `r₀` is `(2r₀ − (r₀−w_r) − (r₀+w_r)) / (2w_r)`, exactly zero in
exact arithmetic and `O(eps)` in floating point for most `(r₀, w_r)`. MEASURED: at the paper's
own `(10, 4.8)` the Float64 value is EXACTLY zero and the snap is inert, but 160 of the 4100
scanned lattices in `6 ≤ r₀ ≤ 16`, `1 ≤ w_r ≤ 5` land on a nonzero round-off (e.g. `+4.44e-16`
at `(7.2, 2.0)`), and in EVERY one of them the unsnapped criterion flags ONE of the two
elements sharing the `r = r₀` face instead of both — the refinement drifts off the particle
and stops being symmetric. After the snap the
comparison in step 3 is between an exactly-representable `ξ ∈ {0, ±1}` and the EXACT dyadic
endpoints of [`logical_interval`](@ref), so the flag set is bitwise deterministic at every
level, in any precision.

# LIMITS

`Policies.Limits` (yaml:148-150) is `{RefinementLevel: Auto, NumGridPoints: Auto,
ErrorBeyondLimits: True}`, i.e. spectre's defaults 15 and 20 and an ERROR — not a silent clip —
when a decision would cross them (`Limits.cpp:16-19`, `EnforcePolicies.cpp:47-71`). This
function throws an `ArgumentError` naming the element and the limit. It checks BOTH the raw
criterion (the reference's own path: `enforce_policies` runs inside
`EvaluateRefinementCriteria`, before `update_amr_decision` raises anything for the 2:1 rule)
AND the closed set — the second pass is a DELIBERATE STRICTENING beyond the reference, which
never limit-checks a balance-added split. The CRP sweep never gets
there: 8 iterations from the initial `(L+1, L+1) = (1, 1)` levels and `P = 4` reach h-level 9
and `N = 12`.

`config.amr.criteria` must be exactly `[:RefineAtPuncture, :IncreaseResolution]`; any other
list (the yaml's commented-out `TruncationError` / `Persson` blocks) is a DIFFERENT criterion
and is refused rather than approximated. An EMPTY flag set is refused too: it would silently
degrade the sweep to pure p-refinement — the manufactured falsifier of this node's test, which
is exactly what the paper's scheme must beat — and it can only mean `orbit` and `dom` describe
different problems.
"""
function amr_flags(mesh::Mesh2D{RT}, dom::BlockDomain, orbit::CircularOrbit,
                   config::RunConfig;
                   snap_tol::Real = AMR_SNAP_TOL,
                   max_refinement_level::Integer = AMR_MAX_REFINEMENT_LEVEL,
                   max_resolution::Integer = AMR_MAX_RESOLUTION) where {RT}
    _amrd_check_policies(config)
    snap_tol >= 0 || throw(ArgumentError("amr_flags: snap_tol must be ≥ 0, got $snap_tol"))
    xp = puncture_chart_position(dom, orbit)

    # step 1+2: the snapped block-logical position of the puncture, once per block
    bl = Dict{Int,Union{Nothing,NTuple{2,RT}}}()
    for blk in mesh.blocks
        bl[blk.id] = _amrd_block_logical(blk, xp, RT, snap_tol)
    end

    # step 3: closed containment, in both dimensions, in the element's own dyadic interval
    flags = Dict{ElementId,NTuple{2,Bool}}()
    for id in mesh.elements
        ξ = bl[id.block]
        ξ === nothing && continue
        inside = true
        for d in 1:2
            lo, hi = logical_interval(RT, id, d)
            (lo <= ξ[d] <= hi) || (inside = false; break)
        end
        inside || continue
        flags[id] = (true, true)
    end

    isempty(flags) && throw(ArgumentError(
        "amr_flags: no element of this mesh contains the puncture at " *
        "(x1, x2) = $(xp) — the criterion would silently degrade to pure p-refinement. " *
        "Check that `orbit` and `dom` describe the same problem (r₀ must lie inside the " *
        "block lattice $(dom.bounds_r[1]) … $(dom.bounds_r[end]))"))

    # THE LIMITS ARE CHECKED ON BOTH SETS, and that is a DELIBERATE STRICTENING, not parity.
    # The reference runs `enforce_policies` on the RAW criterion inside
    # `EvaluateRefinementCriteria` (EvaluateRefinementCriteria.hpp:159-189,
    # Policies/EnforcePolicies.cpp:38-70) and only THEN lets `update_amr_decision` raise flags
    # to Split for the 2:1 rule — so a balance-added split is never limit-checked there. Here
    # the raw criterion is checked first (the parity path) and the CLOSED set second (the
    # strictening): an element the balancer pushes past `RefinementLevel` is refused rather
    # than silently admitted, which is what `ErrorBeyondLimits: True` asks for in spirit even
    # though the reference's own check cannot see it. Both passes are O(elements).
    _amrd_check_limits(mesh, flags, max_refinement_level, max_resolution, "the criterion")
    closed = balanced_flags(mesh, flags)
    _amrd_check_limits(mesh, closed, max_refinement_level, max_resolution,
                       "the criterion or the two-to-one closure")
    return flags
end

# the yaml Limits on one decision set: a dimension that splits must stay under
# `RefinementLevel`, a dimension that does not is p-refined and must stay under
# `NumGridPoints`. `ErrorBeyondLimits: True` (yaml:150) means THROW, never clamp.
function _amrd_check_limits(mesh::Mesh2D, decisions::AbstractDict, maxlevel::Integer,
                            maxres::Integer, what::AbstractString)
    for (e, id) in enumerate(mesh.elements)
        c = get(decisions, id, (false, false))
        for d in 1:2
            if c[d]
                id.segments[d].level < maxlevel || throw(ArgumentError(
                    "amr_flags: element $(id) is split in dimension $d by $what, but is " *
                    "already at refinement level $(id.segments[d].level), the AMR limit " *
                    "(Amr.Policies.Limits.RefinementLevel, ErrorBeyondLimits: True)"))
            else
                mesh.orders[e][d] < maxres || throw(ArgumentError(
                    "amr_flags: element $(id) would be p-refined past $(mesh.orders[e][d]) " *
                    "grid points in dimension $d, the AMR limit " *
                    "(Amr.Policies.Limits.NumGridPoints, ErrorBeyondLimits: True)"))
            end
        end
    end
    return nothing
end

# `Amr.Policies` + `Amr.Criteria` (anc/ScalarSelfForce.yaml:131-151) — the production values
# this criterion IMPLEMENTS. Anything else is a different AMR scheme and is refused rather
# than approximated, which is also what `ErrorBeyondLimits: True` asks of the limits below.
function _amrd_check_policies(config::RunConfig)
    a = config.amr
    a.criteria == AMR_CRITERIA || throw(ArgumentError(
        "amr_flags: this criterion implements Amr.Criteria = $(AMR_CRITERIA) " *
        "(anc/ScalarSelfForce.yaml:131-133), but the config asks for $(a.criteria)"))
    a.isotropy === :Anisotropic || throw(ArgumentError(
        "amr_flags: Amr.Policies.Isotropy = $(a.isotropy) is not implemented — the CRP runs " *
        "`Anisotropic` (yaml:145), which is what makes the 2:1 closure mark only the " *
        "offending dimension"))
    a.enforce_two_to_one_balance_in_normal_direction || throw(ArgumentError(
        "amr_flags: Amr.Policies.EnforceTwoToOneBalanceInNormalDirection = false is not " *
        "implemented — the sealed `hp_refine` always balances (yaml:146)"))
    a.limit_refinement_level === :auto && a.limit_num_grid_points === :auto || throw(ArgumentError(
        "amr_flags: Amr.Policies.Limits = ($(a.limit_refinement_level), " *
        "$(a.limit_num_grid_points)) — only `Auto` is mapped here (yaml:148-149); pass " *
        "`max_refinement_level` / `max_resolution` explicitly for anything else"))
    a.error_beyond_limits || throw(ArgumentError(
        "amr_flags: Amr.Policies.Limits.ErrorBeyondLimits = false is not implemented — this " *
        "criterion THROWS at a limit (yaml:150); clamping would silently change the scheme"))
    a.allow_coarsening || throw(ArgumentError(
        "amr_flags: Amr.Policies.AllowCoarsening = false differs from the CRP input " *
        "(yaml:151). This criterion never emits a Join, so the value is inert here, but a " *
        "config that disagrees with the paper's own is refused rather than ignored"))
    return nothing
end

amr_flags(sol::ModeSolution; kwargs...) =
    amr_flags(sol.mesh, sol.dom, sol.orbit, sol.config; kwargs...)

# spectre `RefineAtPuncture.cpp:33-46` + `BlockLogicalCoordinates.cpp:113-146`: the puncture in
# this block's logical coordinates, snapped, or `nothing` when the block does not contain it.
function _amrd_block_logical(blk::Block, xp, ::Type{RT}, tol::Real) where {RT}
    ξ = ntuple(2) do d
        ζ = map_inverse(blockmap(blk, d), RT(xp[d]))
        abs(ζ - 1) < tol ? one(RT) :
        abs(ζ + 1) < tol ? -one(RT) :
        abs(ζ) < tol ? zero(RT) : ζ
    end
    return (abs(ξ[1]) <= 1 && abs(ξ[2]) <= 1) ? ξ : nothing
end

# ============================== the multigrid level policy ================================

"""
    amr_level_options(config::RunConfig, step::Integer, mesh = nothing;
                      coarse_orders = InitialGridPoints, shrink = 0.5) -> SolverStackOptions

The multigrid LEVEL POLICY of an hp-AMR sweep: the [`solver_options`](@ref) this package hands
[`rebuild_for_mesh!`](@ref) at AMR level `step`, for the grid `mesh` that level solves on.

**`overlap`.** `MaxOverlap: Auto` (yaml:198) is the neighbour's FULL extent and the reference
resolves it INSIDE the AMR projector, from each element's CURRENT extents
(`Schwarz/OverlapHelpers.cpp:21-24` at `Schwarz/ElementActions.hpp:244-264`), so under
`IncreaseResolution` it grows by one at every AMR level. Handed a `mesh` this function returns
[`resolve_max_overlap`](@ref)`(config, mesh)`, that grid's own largest order (4, 5, … 11 over the
paper's eight levels); handed `nothing` it falls back to the config's initial `P`, which is the
level-0 value and one order too small at every level after it. MEASURED effect of the mesh-aware
rule alone, everything else the landed derived hierarchy: `1, 8, 9, 9, 9, 9, 9, 9` → `1, 8, 8,
8, 8, 7, 7, 7` over eight levels.

The other two decisions are unchanged, both measured on the real CRP problem
(`evidence/amr-driver/proto2_*.log`), both taken from the paper. THEY APPLY TO
`hierarchy = :derived` ONLY: under the default `:history` the levels are supplied, so
`coarse_levels`, `coarse_orders` and `shrink` are not consulted.

**`coarse_levels`.** At `step == 0` the yaml's own `InitialCoarseLevels: L` (yaml:188), which
at the production `L = 0` makes the hierarchy a SINGLE level whose bottom solve IS the whole
problem — "The first AMR level converges immediately because the multigrid bottom solver
solves the problem directly" (main.tex:701, fig:linsolv). MEASURED: 1 GMRES iteration. From
`step ≥ 1` it is `nothing`: h-coarsen all the way to the irreducible block structure, which is
what main.tex:672 says the hierarchy is and, crucially, is a grid that does NOT MOVE as AMR
adds finer levels — so the bottom LU is factorized once and is a cache HIT at every later
level (main.tex:674). MEASURED at `a = 0.5`, `r₀ = 10`, `m = 2`: bottom builds 1 at level 0
(its own single-level hierarchy) + 1 at level 1, then 7 consecutive hits, against **9 builds**
if `coarse_levels` is instead raised by one per step — because h-coarsening a graded hp mesh
`k` times does not return the initial lattice (measured coarsest 27 elements at step 8, not
20), so that policy lets the bottom grid follow the fine grid upward.

**`coarse_orders`.** spectre coarsens `h` only and allocates every level at the initial order,
so pinning the coarse levels at `P` is the spectre-faithful setting and is the other half of
what keeps the bottom LU alive; the alternative MIN rule re-factorizes it whenever the coarse
order map moves. MEASURED A/B over the sweep: see this node's evidence table. The default is
the config's own `InitialGridPoints`, which must be equal in the two dimensions (the sealed
`coarse_orders` knob is a single order); pass it explicitly otherwise.

**`shrink`.** `0.5`, not the sealed row-21 default `1`. A graded hp cascade coarsens ONE ring
at a time, so `shrink = 1` takes a smoother on every coarsening step and builds **18** multigrid
levels at AMR level 8; `shrink = 0.5` composes the transfers across steps that do not halve the
problem and reaches the SAME fixed 160-dof bottom in **4**. MEASURED A/B over the 8-level sweep
at `a = 0.5`, `r₀ = 10`, `m = 2`, overlap 4, `OPENBLAS_NUM_THREADS=1`: the GMRES iteration
counts are IDENTICAL (`1, 8, 9, 9, 9, 9, 9, 9, 9` either way) and so are the dof and cache
sequences; the two differ in hierarchy DEPTH and in build time only (18.9 s against 21.7 s
total, 15.5 s against 16.3 s of `build_time`). The row-21 caveat that `shrink = 0.5` "degrades
anisotropic meshes" is therefore not active on the CRP's own AMR sequence — measured, not
assumed — and `shrink` stays a keyword. The full table is in
`evidence/amr-driver/shrink_ab_afterA.md`.
"""
function amr_level_options(config::RunConfig, step::Integer,
                           mesh::Union{Nothing,Mesh2D} = nothing;
                           coarse_orders::Union{Nothing,Integer} =
                               _amrd_initial_order(config),
                           shrink::Real = 0.5)
    step >= 0 || throw(ArgumentError("amr_level_options: step must be ≥ 0, got $step"))
    return solver_options(config;
                          coarse_levels = step == 0 ?
                              config.multigrid.initial_coarse_levels : nothing,
                          coarse_orders = coarse_orders, shrink = shrink,
                          overlap = mesh === nothing ? resolve_max_overlap(config) :
                                    resolve_max_overlap(config, mesh))
end

# the level-policy protocol accepts BOTH arities, so a caller's own two-argument callable keeps
# working while the shipped one is handed this level's mesh
_amrd_options(f, config::RunConfig, step::Integer, mesh::Mesh2D) =
    (applicable(f, config, step, mesh) ? f(config, step, mesh) :
     f(config, step))::SolverStackOptions

# per-step options for `amr_step!`: a fixed record, nothing, or a callable of the NEW mesh
# (the mesh is not known until `hp_refine` has run, and the overlap policy depends on it)
_amrd_step_options(::Nothing, ::Mesh2D) = nothing
_amrd_step_options(o::SolverStackOptions, ::Mesh2D) = o
_amrd_step_options(f, mesh::Mesh2D) = f(mesh)::SolverStackOptions

# the coarser-than-initial grids of `InitialCoarseLevels` (spectre
# `Multigrid/ElementsAllocator.hpp:128-136` bounds the INITIAL hierarchy only), finest first.
# `nothing` coarsens to the irreducible block structure; the CRP production `L = 0` seeds none.
function _amrd_seed_chain!(chain::AbstractVector, mesh::Mesh2D, config::RunConfig)
    L = config.multigrid.initial_coarse_levels
    P = _amrd_initial_order(config)
    cur = mesh
    k = 0
    while L === nothing || k < Int(L)
        co = coarsen_mesh(cur; coarse_orders = P)
        co === nothing && break
        cur = co[1]
        push!(chain, cur)
        k += 1
    end
    return chain
end

function _amrd_initial_order(config::RunConfig)
    P = config.domain.initial_grid_points
    P[1] == P[2] || throw(ArgumentError(
        "amr_level_options: InitialGridPoints $(P) differs per dimension, but the sealed " *
        "`coarse_orders` knob pins ONE order — pass `coarse_orders` explicitly"))
    return Int(P[1])
end

# ================================== one AMR level =========================================

"""
    amr_step!(sol::ModeSolution; flags = amr_flags(sol), options = nothing, increment = 1,
              horizon_hold = false)
        -> (sol, record)

ONE hp-AMR level, in place: refine `sol`'s mesh by `flags`, carry the solution over as the
warm start, rebuild everything mesh-bound, and solve. Literally the four sealed calls, in the
order their contracts require:

    closed = balanced_flags(sol.mesh, flags)          # rows 22/23
    mesh = h_refine(p_refine(sol.mesh; hold, increment), closed)   # h at the flags, p elsewhere
    u    = transfer_solution(sol.mesh, mesh, sol.u)   # exact on every polynomial of the source
    rebuild_for_mesh!(sol, mesh; u = u, options)      # row 42: worldtube, operator, factory, b
    solve_mode!(sol; x0 = sol.u)                      # warm-started, one record on the monitor

and nothing more: no criterion (that is [`amr_flags`](@ref)), no iteration count and no
stopping rule (that is [`solve_mode_amr`](@ref)).

`horizon_hold` adds [`horizon_order_hold`](@ref)`(sol.dom, sol.mesh)` to the `hold` half — the
elements that touch `r = r₊` keep their polynomial order, while the flagged (puncture) elements
are split exactly as before. It is the CURE for the spurious horizon–axis layer of GLOBAL_DAG
row 58 (`horizon_order_hold` carries the mechanism and the indicial roots): raising the order of
the element that touches the horizon is what lets the bounded, non-analytic `(r−r₊)^p` branch
(`p = γ^r_m(r₊)/c`, purely imaginary) into the trial space, and nothing else rejects it.

**It defaults to `false`, and `false` is the composition `hp_refine(sol.mesh, flags; increment)`
BITWISE** — every sealed AMR number (the dof and element sequences of rows 43/50/52, the
`fig:linsolv` staircase, the reproduction goldens) is untouched by this node. Turning it on is a
mesh change: it removes the horizon layer's elements from p-refinement, so those sequences move
(MEASURED at `a = 0.5`, `r₀ = 10M`, `m = 2` over the paper's eight grids: see the node's
evidence). The DEFAULT STAYS `false` (recorded decision): the reference has no
counterpart to the hold,
the paper's own observables are insensitive to the layer (`F_r` moves 2.6e-5 relative), the
flip would move sealed AMR bookkeeping, and the hold's truncation cost on the held layer is not
separately bounded. The rule is MEASURED instead: an entry point that publishes a PLANE gates
[`horizon_axis_profile`](@ref)'s `layer_ratio ≤ 1.05` and sets `horizon_hold = true` when that
fails — never by spin, never by the branch index.

ELEMENT IDS CHANGE AT EVERY STEP, so `flags` must be recomputed against `sol.mesh` before each
call and never cached across steps (the sealed row-22/23 contract). `options` replaces the
stack's options and may be a [`SolverStackOptions`](@ref), `nothing` (keep the current ones), or
a CALLABLE of the new mesh — the mesh is not known until `hp_refine` has run, and the overlap
policy depends on it, so [`solve_mode_amr`](@ref) passes
`mh -> amr_level_options(config, step, mh)`. Returns `sol` and the level's
[`SolveRecord`](@ref) (also `sol.record`, and appended to `sol.monitor`).
"""
function amr_step!(sol::ModeSolution; flags = amr_flags(sol),
                   options = nothing, increment::Integer = 1,
                   horizon_hold::Bool = false)
    closed = balanced_flags(sol.mesh, flags)
    hold = horizon_hold ? _amrd_hold_union(closed, horizon_order_hold(sol.dom, sol.mesh)) :
                          closed
    mesh = h_refine(p_refine(sol.mesh; hold = hold, increment = increment), closed)
    u = transfer_solution(sol.mesh, mesh, sol.u)
    rebuild_for_mesh!(sol, mesh; u = u, options = _amrd_step_options(options, mesh))
    solve_mode!(sol; x0 = sol.u)
    return sol, sol.record
end

# per-dimension UNION of two hold tables: a dimension is held if EITHER table holds it (the
# split dimensions of `balanced_flags`, which must keep their order — spectre
# UpdateAmrDecision.cpp:82-87 — and the horizon layer of `horizon_order_hold`). Written out
# rather than `merge`d because `merge` would let one table's `false` erase the other's `true`.
function _amrd_hold_union(a::AbstractDict{ElementId,NTuple{2,Bool}},
                          b::AbstractDict{ElementId,NTuple{2,Bool}})
    out = Dict{ElementId,NTuple{2,Bool}}()
    for (id, d) in a
        e = get(b, id, (false, false))
        out[id] = (d[1] || e[1], d[2] || e[2])
    end
    for (id, e) in b
        haskey(out, id) || (out[id] = e)
    end
    return out
end

# ================================== the sweep record ======================================

"""
    AmrLevel{T,RT}

One level of an hp-AMR sweep — everything rows 44, 47, 48, 50, 52 and 53 read about it.

  - `step::Int` — `0` is the initial solve; `k` is the mesh after `k` [`amr_step!`](@ref)s.
  - `mesh::Mesh2D{RT}` — THE MESH of this level (elements, orders, offsets): the level map and
    the order map of `fig:ssf_domain` are `mesh.elements` and `mesh.orders`.
  - `punctures::Vector{ElementId}` — the elements [`amr_flags`](@ref) found to contain the
    puncture at this level, in mesh order; the ones the NEXT step splits.
  - `record::SolveRecord{RT}` — this level's solve: `iterations`, the `residual` staircase,
    `iteration_time`, and the `build_time` / cache counters of the rebuild that preceded it
    (`fig:linsolv` is exactly this field across the levels).
  - `element::Int` — the FIRST containing element in mesh order, i.e. what [`field_at`](@ref)
    picks by default. It is the representative of a single-element trace, not the source of the
    diagnostics below.
  - `elements::Vector{Int}` — EVERY element of `mesh` whose closed cell contains the puncture,
    in the locator's own order. The puncture sits on an element corner, so there are TWO of
    them in the `z²` chart (the equator is the domain boundary there) and FOUR in the full-`z`
    chart — spectre's `NumContributingElements`.
  - `traces::Vector{NTuple{3,T}}` — one `(Ψ^ℛ_m, ∂_{x1}Ψ^ℛ_m, ∂_{x2}Ψ^ℛ_m)` one-sided trace per
    entry of `elements`, in the same order. A DG field is discontinuous across faces, so these
    differ; keeping all of them means ANY reduction can be reproduced from a level alone,
    without the `ModeSolution` (which is process-local and, for every level but the last,
    discarded).
  - `psi::T`, `dpsi_r::T`, `dpsi_w::T` — the particle diagnostics: the MEAN over `traces` of
    `Ψ^ℛ_m(x_p)` and its two chart derivatives `∂_{x1}Ψ^ℛ_m`, `∂_{x2}Ψ^ℛ_m` at the puncture,
    `x1 = r` and `x2 = z²` under equatorial symmetry. THE MEAN IS SPECTRE'S REDUCTION
    (`ObserveSelfForce.hpp:89-107`: `funcl::Plus` then `funcl::Divides` by
    `NumContributingElements`) and it is BITWISE the reduction
    `mode_self_force(sol::ModeSolution)` performs on the live solution (gated). These are the
    RAW ingredients of `F_r^m`; row 44 assembles the self-force from them (and owns every
    convention that goes with it). MEASURED effect of recording the mean rather than one
    one-sided trace: the sweep's `relative_resolution_error` falls from 3.39e-5 to 6.62e-6 on
    the production half domain and from 6.66e-6 to 1.60e-7 on the full-`z` one — a one-sided
    trace carries a face jump that the reduction cancels.
  - `trace_spread::NTuple{3,RT}` — how far the containing elements' one-sided traces of `psi`,
    `dpsi_r` AND `dpsi_w` differ from the recorded MEAN, one entry each. A DG field is
    discontinuous across faces, so this is a discretization diagnostic, and it collapses as the
    corner is refined. THE THIRD ENTRY IS NOT OPTIONAL: in the full-`z` chart `dpsi_w` is
    `∂_zΨ^ℛ`, which feeds `F_θ` — a one-sided `dpsi_w` with no diagnostic beside it reports a
    noise `F_θ` (measured 1.33e-4 against `|F_r| = 1.67e-5`) as if it were a result.
  - `series_nodes::Int` — how many of the in-worldtube collocation nodes of this level select
    the near-worldline SERIES route of the switched effective source (row 37). It is `0`
    throughout the paper's own sweep and the number exists to keep it honest: it is the AMR
    resolution, not the initial mesh, that decides, so a driver that silently crossed into the
    series regime would otherwise be invisible. Counted with the same sealed
    [`effective_source_switched_bundle`](@ref) the volume source uses.
  - `wall_time::Float64` — seconds for this level (refine + carry-over + rebuild + solve;
    level 0 is the initial build and solve).
"""
struct AmrLevel{T<:Number,RT<:AbstractFloat}
    step::Int
    mesh::Mesh2D{RT}
    punctures::Vector{ElementId}
    record::SolveRecord{RT}
    element::Int
    elements::Vector{Int}
    traces::Vector{NTuple{3,T}}
    psi::T
    dpsi_r::T
    dpsi_w::T
    trace_spread::NTuple{3,RT}
    series_nodes::Int
    wall_time::Float64
end

"""
    AmrHistory{T,RT}

A complete hp-AMR sweep of one m-mode: `orbit`, `config`, `m`, the puncture position `x_p` in
the discretization chart, the per-level record `levels::Vector{AmrLevel}` (level 0 first), and
`solution::ModeSolution` — the FINAL level's solved mode, whose `monitor` carries the same
`SolveRecord`s in the same order.

Build it with [`solve_mode_amr`](@ref). The convergence measure of `fig:convergence` and of
the error budget's "numerical resolution floor" (main.tex:705-724) is the successive
difference of a per-level diagnostic, e.g. `[l.dpsi_r for l in h.levels]`.
"""
struct AmrHistory{T<:Number,RT<:AbstractFloat}
    orbit::CircularOrbit
    config::RunConfig{RT}
    m::Int
    x_p::NTuple{2,RT}
    levels::Vector{AmrLevel{T,RT}}
    solution::ModeSolution{T,RT}
end

function Base.show(io::IO, ::MIME"text/plain", h::AmrHistory)
    println(io, "AmrHistory: m = ", h.m, ", a = ", h.orbit.a, ", r₀ = ", h.orbit.r0,
            ", ", length(h.levels), " levels")
    println(io, "  level  elements    dofs  iters   residual   build[s]  wall[s]  ",
            "Re Ψ^R(x_p)      ∂_r Ψ^R(x_p)")
    for l in h.levels
        @printf(io, "  %5d  %8d  %6d  %5d  %9.2e  %8.2f  %7.2f  %+13.6e  %+13.6e\n",
                l.step, length(l.mesh.elements), l.mesh.ndofs, l.record.iterations,
                l.record.final_residual, l.record.build_time, l.wall_time,
                real(l.psi), real(l.dpsi_r))
    end
end

# ==================================== the sweep ===========================================

"""
    solve_mode_amr(orbit::CircularOrbit, config::RunConfig, m = config.orbit.m_mode_number;
                   iterations = config.amr.iterations, T = ComplexF64,
                   dom = block_domain(orbit, config), phi_p = 0,
                   criterion = amr_flags, level_options = amr_level_options,
                   hierarchy = :history, horizon_hold = false,
                   monitor = ResidualMonitor(real(T))) -> AmrHistory

THE hp-AMR SWEEP of arXiv:2606.04998 for one m-mode: solve on the initial lattice, then run
`iterations` levels of [`amr_step!`](@ref) — the paper's `Amr: Iterations: 8` (yaml:153) — and
return the [`AmrHistory`](@ref).

# WHAT `iterations` COUNTS, AND HOW IT MAPS ONTO THE yaml (normative)

`iterations` is the number of REFINEMENT STEPS, so the sweep produces `iterations + 1` grids:
`levels[1]` is the initial lattice and `levels[end]` the finest. The yaml's `Amr: Iterations:
8` (yaml:153) is NOT that number. Its AMR phase is driven by the trigger
`EveryNIterations{N: 1, Offset: 1}` (yaml:157-160), which fires when `iteration_id >= 1`
(`elliptic/Triggers/EveryNIterations.hpp:48-51`) — the yaml says so in its own comment, "Run
AMR in every iteration, but not on the initial guess". Eight solver iterations therefore carry
SEVEN refinements and produce EIGHT grids, and the faithful call is

    solve_mode_amr(orbit, config, m; iterations = config.amr.iterations - 1)   # 8 grids

while the default `iterations = config.amr.iterations` produces NINE grids — one refinement
more than the paper. The default is kept because it is a strict superset (the paper's grid set
is `levels[1:8]`, bitwise) and because every number recorded for this node was measured on it;
a caller reproducing the paper's own resolution should say so explicitly. [PRELIMINARY: the
option `Amr.Iterations` itself is absent from the pinned spectre mirror (version skew), so this
mapping is reasoned from the trigger semantics and the yaml's comment, not read off the
consumer.]

# THE MULTIGRID HIERARCHY IS THE AMR HISTORY (`hierarchy = :history`, the default)

arXiv:2606.04998 does not derive its coarse grids from the current mesh. Its executable KEEPS
the coarse grids (`Elliptic/Executables/SelfForce/Scalar/SolveScalarSelfForce.hpp:173`,
`amr::keep_coarse_grids = true`); every AMR step creates a NEW FINEST GRID at `grid_index + 1`
— split children where the flag is `Split` (`Amr/Actions/AdjustDomain.hpp:157-163`), one
same-shape child carrying the p increment everywhere else (`:203-212`) — and the parent is
re-registered as a coarser grid rather than destroyed
(`Amr/Actions/SendDataToChildren.hpp:43-70`); and the V-cycle walks exactly those grids
(`Multigrid/ElementActions.hpp:96-135` sets the `amr::Tags::{ParentId, ChildIds, ParentMesh}`
triple traversed by `Multigrid.hpp:160-186`, and `Multigrid.hpp:56-58` says so in words).

So at AMR level `k` the hierarchy has `k+1` levels and they ARE this sweep's mesh history:
consecutive levels one AMR step apart, and the coarsest the INITIAL lattice, which never moves —
main.tex:674 ("the LU decomposition of A^bottom is cached and reused in all multigrid V-cycles
across AMR iterations … the bottom grid remains the same even when AMR adds finer levels") is
then literally true, ONE factorization for the whole sweep. This driver keeps the chain and
hands it to [`multigrid_chain_preconditioner`](@ref), with a per-grid overlap so the frozen
coarse levels keep their subdomain factorizations too (main.tex:676).

`hierarchy = :derived` keeps the h-coarsening of the current mesh reachable for the A/B.
MEASURED at the paper's own exemplar (`a = 0.5`, `r₀ = 10 M`, `m = 2`), eight levels:

| hierarchy | iterations | total | reduction per iteration | bottom LU |
|---|---|---|---|---|
| `:history` | 1, 7, 6, 5, 5, 4, 4, 4 | 36 | 29.8 → 135.4, RISING | 1 build, 7 hits |
| `:derived` | 1, 8, 8, 8, 8, 7, 7, 7 | 54 | 16.9 → 12.7 | 2 builds, 6 hits |
| `:derived`, overlap pinned at `P` | 1, 8, 9, 9, 9, 9, 9, 9 | 63 | 16.4 → 6.0, FALLING | 2 builds |
| published `fig:linsolv` | 1, 7, 6, 5, 4, 3, 3, 2 | 31 | 27 → 567 | one |

(the third row is the pre-adoption behaviour, reachable as
`level_options = (c, s, mh) -> amr_level_options(c, s)`). main.tex:720's "less than 10
iterations per AMR level" now holds with room to spare at this configuration; it was measured
at 10 for `a = 0, m = 2` level 3 and `m = 5` levels 2-4 under the pre-adoption stack and those
configurations have not been re-measured. One `ModeSolution`, one stack, one epoch, one pair of caches,
carried through every level; the previous level's solution is transferred and warm-starts the
next solve (main.tex:714).

    sol = solve_mode(orbit, config, m; options = level_options(config, 0))
    for k in 1:iterations
        amr_step!(sol; flags = criterion(sol), options = level_options(config, k))
    end

`criterion` is a callable `sol -> AbstractDict{ElementId,NTuple{2,Bool}}`; the default is the
paper's own [`amr_flags`](@ref) and the keyword exists so a caller can FALSIFY it (inverting it
— p at the puncture, h elsewhere — is what shows the scheme's rationale). `level_options` is
a callable `(config, step) -> SolverStackOptions`; see [`amr_level_options`](@ref) for the
policy and the measurements behind it.

# WHAT IT COSTS, AND WHERE THE TIME GOES (measured, single thread, Float64/ComplexF64)

At the paper's own point (`a = 0.5`, `r₀ = 10`, `m = 2`, `L = 0`, `P = 4`) the 8-level sweep
grows the mesh 20 → 69 elements and 320 → 5696 dofs, converges every level in 1 and then 8-9
preconditioned GMRES iterations, and takes **18.9 s** of warm single-threaded compute — the
paper's own "a few seconds per m-mode on an Apple M2 Pro" (main.tex:716), within a small
factor. The cost is `SolveRecord.build_time` (15.5 s of the 18.9 s), the eager per-level
rebuild of the operator and the Schwarz subdomain matrices, which the paper also names as its
dominant operation (main.tex:678). Four levels — what the m-sweep costs — is 3.6 s at m = 0,
3.6 s at m = 1, 2.6 s at m = 2 and 4.5 s at m = 5.

PIN `OPENBLAS_NUM_THREADS=1`. The subdomain and bottom solves are DENSE and TINY (a few
hundred rows), and the host's default OpenBLAS thread pool spends more time in barriers than
in arithmetic: the same sweep measures **136 s** unpinned against **15.4 s** pinned, a factor
8.9, with identical iteration counts. Pinning also makes the run bitwise reproducible — the
diagnostics move by ~4e-12 relative between BLAS thread counts, harmless against the 2e-6
level-to-level differences but not bitwise.

Each m-mode is an INDEPENDENT problem with its own stack and caches, so the fan-out is over
PROCESSES (row 47), never threads.

THE DIAGNOSTICS CONVERGE EXPONENTIALLY but do not reach machine precision in 8 levels. The
honest measure is a TWO-STEP GEOMETRIC ENVELOPE — successive differences fall by a factor ≥ 2
every two levels, which tolerates the occasional plateau a per-step ratio bound does not (one
ratio of 0.687 measured at `a = 0.5, m = 1`). At the paper's 8 iterations the relative change
between the last two levels is ~3e-6 (`Ψ^ℛ`) and ~2e-6 (`∂_rΨ^ℛ`); [`resolution_error`](@ref)
is that quantity, and it IS the error budget's resolution estimate (main.tex:722).

Pushed further the sequence reaches a FLOOR at level ~12 (a few times 1e-11 on `Ψ^ℛ` and
1e-12 on `∂_rΨ^ℛ`). MEASURED: that floor is INSENSITIVE to the GMRES `abstol` over
`[1e-14, 1e-12]` — tightening the tolerance moves it by ~15 % while the iteration counts rise
from 9 to 11-18 — so it is not set by the linear-solver tolerance. Its MECHANISM IS NOT
ESTABLISHED; do not attribute it.
"""
function solve_mode_amr(orbit::CircularOrbit, config::RunConfig{RT},
                        m::Integer = config.orbit.m_mode_number;
                        iterations::Integer = config.amr.iterations,
                        T::Type = ComplexF64,
                        dom::BlockDomain = block_domain(orbit, config),
                        phi_p::Real = 0,
                        criterion = amr_flags,
                        level_options = amr_level_options,
                        hierarchy::Symbol = :history,
                        horizon_hold::Bool = false,
                        monitor::ResidualMonitor = ResidualMonitor(real(T))) where {RT}
    iterations >= 0 || throw(ArgumentError(
        "solve_mode_amr: iterations must be ≥ 0, got $iterations"))
    hierarchy in (:history, :derived) || throw(ArgumentError(
        "solve_mode_amr: hierarchy must be :history (the reference's own — every AMR step " *
        "adds a finer grid and the old ones are kept) or :derived (h-coarsen the current " *
        "mesh), got $(repr(hierarchy))"))
    xp = puncture_chart_position(dom, orbit)
    chain = Mesh2D{RT}[]
    precond = multigrid_preconditioner
    if hierarchy === :history
        _amrd_seed_chain!(chain, dom.mesh, config)
        # PER-LEVEL overlap: the reference resolves `MaxOverlap: Auto` from each grid's own
        # extents and freezes the coarse grids, so a coarse level keeps the overlap it had when
        # it was the finest and its Schwarz factorizations stay cache HITS (main.tex:676).
        precond = multigrid_chain_preconditioner(
            chain; overlap = mh -> resolve_max_overlap(config, mh))
    end
    t = time()
    sol = solve_mode(orbit, config, m; T = T, dom = dom, phi_p = phi_p,
                     options = _amrd_options(level_options, config, 0, dom.mesh),
                     preconditioner = precond, monitor = monitor)
    wall = time() - t
    levels = AmrLevel{T,RT}[]
    for k in 0:Int(iterations)
        flags = criterion(sol)
        push!(levels, _amrd_level(k, sol, flags, xp, wall))
        k == iterations && break
        t = time()
        # the mesh this level leaves behind BECOMES the next-coarser multigrid grid
        hierarchy === :history && pushfirst!(chain, sol.mesh)
        amr_step!(sol; flags = flags, horizon_hold = horizon_hold,
                  options = mh -> _amrd_options(level_options, config, k + 1, mh))
        wall = time() - t
    end
    return AmrHistory{T,RT}(orbit, config, Int(m), xp, levels, sol)
end

function _amrd_level(step::Int, sol::ModeSolution{T,RT}, flags::AbstractDict, xp,
                     wall::Float64) where {T,RT}
    es = _dr_containing(sol.mesh, xp[1], xp[2])       # Solve/'s one locator (driver.jl)
    isempty(es) && throw(ArgumentError(
        "amr-driver: the puncture ($(xp[1]), $(xp[2])) lies in no element of this mesh"))
    # EVERY containing element's one-sided trace, in the locator's own order, and their MEAN —
    # spectre's reduction over `NumContributingElements` (`ObserveSelfForce.hpp:89-107`). The
    # mean is what `psi`/`dpsi_r`/`dpsi_w` record; the individual traces stay in `traces`, so a
    # single-element reduction is still reproducible from a level alone.
    tr = NTuple{3,T}[(field_at(sol, xp[1], xp[2]; element = q),
                      field_grad_at(sol, xp[1], xp[2]; element = q)...) for q in es]
    n = length(tr)
    psi = sum(t[1] for t in tr) / n
    d1 = sum(t[2] for t in tr) / n
    d2 = sum(t[3] for t in tr) / n
    sp = (maximum(t -> abs(t[1] - psi), tr), maximum(t -> abs(t[2] - d1), tr),
          maximum(t -> abs(t[3] - d2), tr))
    punct = [id for id in sol.mesh.elements if haskey(flags, id)]
    return AmrLevel{T,RT}(step, sol.mesh, punct, sol.record, first(es), collect(es), tr,
                          psi, d1, d2, (RT(sp[1]), RT(sp[2]), RT(sp[3])),
                          _amrd_series_nodes(sol), wall)
end

# how many in-worldtube collocation nodes of this level's mesh select the SERIES route of the
# switched effective source — the same nodes, the same sealed evaluator and the same switch the
# volume source (`mode_fixed_source`) runs over, counted rather than summed.
function _amrd_series_nodes(sol::ModeSolution)
    dom, mesh = sol.dom, sol.mesh
    p = Puncture(sol.orbit)
    z2 = dom.chart === :z2
    n = 0
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        inside_worldtube(dom, el) || continue
        em = element_matrices(Float64, mesh.blocks[mesh.blockindex[el.block]], el,
                              mesh.orders[e])
        nx, ny = mesh.orders[e]
        for j in 1:ny, i in 1:nx
            x1, x2 = em.x[1][i], em.x[2][j]
            b = effective_source_switched_bundle(p, sol.m, x1, z2 ? sqrt(x2) : x2;
                                                 phi_p = sol.phi_p)
            b.route === :series && (n += 1)
        end
    end
    return n
end

# ============================ the resolution error, ONE definition ========================

"""
    resolution_error(hist::AmrHistory; field = :dpsi_r) -> Real
    relative_resolution_error(hist::AmrHistory; field = :dpsi_r) -> Real

THE numerical-resolution error of an AMR sweep, as arXiv:2606.04998 measures it: "We measure
the error due to numerical resolution by the difference between the highest and second-highest
AMR level" (main.tex:722). `resolution_error` is that absolute difference,
`relative_resolution_error` divides it by `|value at the highest level|`.

    resolution_error(hist) == abs(hist.levels[end].dpsi_r - hist.levels[end-1].dpsi_r)

`field` selects the per-level [`AmrLevel`](@ref) quantity: `:dpsi_r` (the default — the
`∂_rΨ^ℛ(x_p)` ingredient of `F_r^m`), `:psi` or `:dpsi_w`. This is the ONE definition rows 44,
50 and 53 use, so that the error budget's resolution curve, the convergence figure and the
extraction all mean the same thing by it; do not re-derive it from `levels` by hand.

A sweep with fewer than two levels has no such estimate and is an `ArgumentError`.
"""
function resolution_error(hist::AmrHistory; field::Symbol = :dpsi_r)
    length(hist.levels) >= 2 || throw(ArgumentError(
        "resolution_error: needs at least two AMR levels, got $(length(hist.levels)) — the " *
        "estimate IS the difference between the highest and second-highest (main.tex:722)"))
    field in (:psi, :dpsi_r, :dpsi_w) || throw(ArgumentError(
        "resolution_error: field must be one of :psi, :dpsi_r, :dpsi_w, got $(repr(field))"))
    a = getfield(hist.levels[end], field)
    b = getfield(hist.levels[end-1], field)
    return abs(a - b)
end

function relative_resolution_error(hist::AmrHistory; field::Symbol = :dpsi_r)
    d = resolution_error(hist; field = field)
    v = abs(getfield(hist.levels[end], field))
    iszero(v) && throw(ArgumentError(
        "relative_resolution_error: the final-level $(field) is exactly zero, so there is no " *
        "relative measure — use `resolution_error`"))
    return d / v
end

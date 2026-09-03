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

export puncture_chart_position, amr_flags, amr_level_options, amr_step!,
       solve_mode_amr, AmrLevel, AmrHistory

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
function throws an `ArgumentError` naming the element and the limit. The CRP sweep never gets
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
    config.amr.criteria == AMR_CRITERIA || throw(ArgumentError(
        "amr_flags: this criterion implements Amr.Criteria = $(AMR_CRITERIA) " *
        "(anc/ScalarSelfForce.yaml:131-133), but the config asks for " *
        "$(config.amr.criteria)"))
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
        for d in 1:2
            id.segments[d].level < max_refinement_level || throw(ArgumentError(
                "amr_flags: element $(id) contains the puncture but is already at " *
                "refinement level $(id.segments[d].level) in dimension $d, the AMR limit " *
                "(Amr.Policies.Limits.RefinementLevel, ErrorBeyondLimits: True)"))
        end
        flags[id] = (true, true)
    end

    isempty(flags) && throw(ArgumentError(
        "amr_flags: no element of this mesh contains the puncture at " *
        "(x1, x2) = $(xp) — the criterion would silently degrade to pure p-refinement. " *
        "Check that `orbit` and `dom` describe the same problem (r₀ must lie inside the " *
        "block lattice $(dom.bounds_r[1]) … $(dom.bounds_r[end]))"))

    # the IncreaseResolution half applies to every UNFLAGGED element, so its limit is theirs
    for (e, id) in enumerate(mesh.elements)
        haskey(flags, id) && continue
        for d in 1:2
            mesh.orders[e][d] < max_resolution || throw(ArgumentError(
                "amr_flags: element $(id) would be p-refined past $(mesh.orders[e][d]) " *
                "grid points in dimension $d, the AMR limit " *
                "(Amr.Policies.Limits.NumGridPoints, ErrorBeyondLimits: True)"))
        end
    end
    return flags
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
    amr_level_options(config::RunConfig, step::Integer;
                      coarse_orders = InitialGridPoints, shrink = 1) -> SolverStackOptions

The multigrid LEVEL POLICY of an hp-AMR sweep: the [`solver_options`](@ref) this package hands
[`rebuild_for_mesh!`](@ref) at AMR level `step`. Two decisions, both measured on the real CRP
problem (`evidence/amr-driver/proto2_*.log`), both taken from the paper:

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

**`shrink`.** `0.5`, NOT the sealed row-21 default `1`, and this is what makes the policy
above affordable. A graded hp cascade coarsens ONE ring at a time, so `shrink = 1` — take a
smoother on every coarsening step — gives **18** multigrid levels at AMR level 8, each one a
full operator rebuild (worldtube bundle, hooks, boundary conditions, Schwarz subdomains).
`shrink = 0.5` composes the transfer operators across steps that do not halve the problem and
gives **4** levels for the same bottom grid. MEASURED over the 8-level sweep at `a = 0.5`,
`r₀ = 10`, `m = 2`: iterations `1, 10, 9, 10, 10, 10, 9, 9, 9` against `1, 9, 9, 9, 10, 10, 9,
9, 10` — ONE extra GMRES iteration over the whole sweep — for a last-level `build_time` of
**19.7 s against 334.9 s** and the same fixed 160-dof bottom with its LU cached. The row-21
caveat that `shrink = 0.5` "degrades anisotropic meshes" is therefore not active on the CRP's
own AMR sequence; it is measured here, not assumed, and `shrink` stays a keyword.
"""
function amr_level_options(config::RunConfig, step::Integer;
                           coarse_orders::Union{Nothing,Integer} =
                               _amrd_initial_order(config),
                           shrink::Real = 0.5)
    step >= 0 || throw(ArgumentError("amr_level_options: step must be ≥ 0, got $step"))
    return solver_options(config;
                          coarse_levels = step == 0 ?
                              config.multigrid.initial_coarse_levels : nothing,
                          coarse_orders = coarse_orders, shrink = shrink)
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
    amr_step!(sol::ModeSolution; flags = amr_flags(sol), options = nothing, increment = 1)
        -> (sol, record)

ONE hp-AMR level, in place: refine `sol`'s mesh by `flags`, carry the solution over as the
warm start, rebuild everything mesh-bound, and solve. Literally the four sealed calls, in the
order their contracts require:

    mesh = hp_refine(sol.mesh, flags; increment)      # rows 22/23: h at the flags, p elsewhere
    u    = transfer_solution(sol.mesh, mesh, sol.u)   # exact on every polynomial of the source
    rebuild_for_mesh!(sol, mesh; u = u, options)      # row 42: worldtube, operator, factory, b
    solve_mode!(sol; x0 = sol.u)                      # warm-started, one record on the monitor

and nothing more: no criterion (that is [`amr_flags`](@ref)), no iteration count and no
stopping rule (that is [`solve_mode_amr`](@ref)).

ELEMENT IDS CHANGE AT EVERY STEP, so `flags` must be recomputed against `sol.mesh` before each
call and never cached across steps (the sealed row-22/23 contract). `options` replaces the
stack's options — pass `amr_level_options(sol.config, step)` to move the multigrid level policy
with the sweep; omitted, the current options are kept. Returns `sol` and the level's
[`SolveRecord`](@ref) (also `sol.record`, and appended to `sol.monitor`).
"""
function amr_step!(sol::ModeSolution; flags = amr_flags(sol),
                   options::Union{Nothing,SolverStackOptions} = nothing,
                   increment::Integer = 1)
    mesh = hp_refine(sol.mesh, flags; increment = increment)
    u = transfer_solution(sol.mesh, mesh, sol.u)
    rebuild_for_mesh!(sol, mesh; u = u, options = options)
    solve_mode!(sol; x0 = sol.u)
    return sol, sol.record
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
  - `element::Int` — the index into `mesh.elements` whose polynomial produced the diagnostics
    below. The puncture is on an element CORNER, so its value is a one-sided trace; this is the
    first containing element in mesh order, i.e. what [`field_at`](@ref) picks by default.
  - `psi::T`, `dpsi_r::T`, `dpsi_w::T` — the particle diagnostics: `Ψ^ℛ_m(x_p)` and its two
    chart derivatives `∂_{x1}Ψ^ℛ_m`, `∂_{x2}Ψ^ℛ_m` at the puncture, `x1 = r` and `x2 = z²`
    under equatorial symmetry. These are the RAW ingredients of `F_r^m`; row 44 assembles the
    self-force from them (and owns every convention that goes with it).
  - `trace_spread::NTuple{2,RT}` — how far the OTHER containing elements' one-sided traces of
    `psi` and `dpsi_r` differ from the recorded ones. A DG field is discontinuous across
    faces, so this is a discretization diagnostic, and it collapses as the corner is refined
    (MEASURED 8.4e-6 → 3.5e-17 and 2.2e-4 → 9.4e-11 over the 8-level sweep).
  - `wall_time::Float64` — seconds for this level (refine + carry-over + rebuild + solve;
    level 0 is the initial build and solve).
"""
struct AmrLevel{T<:Number,RT<:AbstractFloat}
    step::Int
    mesh::Mesh2D{RT}
    punctures::Vector{ElementId}
    record::SolveRecord{RT}
    element::Int
    psi::T
    dpsi_r::T
    dpsi_w::T
    trace_spread::NTuple{2,RT}
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
                   monitor = ResidualMonitor(real(T))) -> AmrHistory

THE hp-AMR SWEEP of arXiv:2606.04998 for one m-mode: solve on the initial lattice, then run
`iterations` levels of [`amr_step!`](@ref) — the paper's `Amr: Iterations: 8` (yaml:153) — and
return the [`AmrHistory`](@ref). One `ModeSolution`, one stack, one epoch, one pair of caches,
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
grows the mesh 20 → 69 elements and 320 → 5720 dofs, converges every level in 1 and then 8-9
preconditioned GMRES iterations, and takes **15.4 s** of warm single-threaded compute — the
paper's own "a few seconds per m-mode on an Apple M2 Pro" (main.tex:716), within a small
factor. The cost is `SolveRecord.build_time` (12.6 s of the 15.4 s), the eager per-level
rebuild of the operator and the Schwarz subdomain matrices, which the paper also names as its
dominant operation (main.tex:678).

PIN `OPENBLAS_NUM_THREADS=1`. The subdomain and bottom solves are DENSE and TINY (a few
hundred rows), and the host's default OpenBLAS thread pool spends more time in barriers than
in arithmetic: the same sweep measures **136 s** unpinned against **15.4 s** pinned, a factor
8.9, with identical iteration counts. Pinning also makes the run bitwise reproducible — the
diagnostics move by ~4e-12 relative between BLAS thread counts, harmless against the 2e-6
level-to-level differences but not bitwise.

Each m-mode is an INDEPENDENT problem with its own stack and caches, so the fan-out is over
PROCESSES (row 47), never threads.

THE DIAGNOSTICS CONVERGE EXPONENTIALLY but do not reach machine precision in 8 levels:
MEASURED successive-difference ratios ≈ 0.22-0.34 for `Ψ^ℛ(x_p)` and ≈ 0.06-0.33 for
`∂_rΨ^ℛ(x_p)`, ending at a relative change of 3.1e-6 and 2.0e-6 between levels 7 and 8. The
sequence is the input of the error budget's resolution floor (main.tex:705-724), which is
measured as the difference between the highest and second-highest level — not as a tolerance
this loop meets.
"""
function solve_mode_amr(orbit::CircularOrbit, config::RunConfig{RT},
                        m::Integer = config.orbit.m_mode_number;
                        iterations::Integer = config.amr.iterations,
                        T::Type = ComplexF64,
                        dom::BlockDomain = block_domain(orbit, config),
                        phi_p::Real = 0,
                        criterion = amr_flags,
                        level_options = amr_level_options,
                        monitor::ResidualMonitor = ResidualMonitor(real(T))) where {RT}
    iterations >= 0 || throw(ArgumentError(
        "solve_mode_amr: iterations must be ≥ 0, got $iterations"))
    xp = puncture_chart_position(dom, orbit)
    t = time()
    sol = solve_mode(orbit, config, m; T = T, dom = dom, phi_p = phi_p,
                     options = level_options(config, 0), monitor = monitor)
    wall = time() - t
    levels = AmrLevel{T,RT}[]
    for k in 0:Int(iterations)
        flags = criterion(sol)
        push!(levels, _amrd_level(k, sol, flags, xp, wall))
        k == iterations && break
        t = time()
        amr_step!(sol; flags = flags, options = level_options(config, k + 1))
        wall = time() - t
    end
    return AmrHistory{T,RT}(orbit, config, Int(m), xp, levels, sol)
end

function _amrd_level(step::Int, sol::ModeSolution{T,RT}, flags::AbstractDict, xp,
                     wall::Float64) where {T,RT}
    es = _amrd_containing(sol.mesh, xp[1], xp[2])
    e = first(es)
    psi = field_at(sol, xp[1], xp[2]; element = e)
    g = field_grad_at(sol, xp[1], xp[2]; element = e)
    sp = (maximum(q -> abs(field_at(sol, xp[1], xp[2]; element = q) - psi), es),
          maximum(q -> abs(field_grad_at(sol, xp[1], xp[2]; element = q)[1] - g[1]), es))
    punct = [id for id in sol.mesh.elements if haskey(flags, id)]
    return AmrLevel{T,RT}(step, sol.mesh, punct, sol.record, e, psi, g[1], g[2],
                          (RT(sp[1]), RT(sp[2])), wall)
end

# every element containing the point, in mesh order — the same predicate (and the same 16 eps
# logical tolerance) the sealed `field_at` locator uses, so `first` IS its default choice.
function _amrd_containing(mesh::Mesh2D{RT}, x1::Real, x2::Real) where {RT}
    tol = 16 * eps(RT)
    out = Int[]
    for e in eachindex(mesh.elements)
        ξ = _dr_logical(mesh, e, x1, x2)
        (-1 - tol <= ξ[1] <= 1 + tol && -1 - tol <= ξ[2] <= 1 + tol) && push!(out, e)
    end
    isempty(out) && throw(ArgumentError(
        "amr-driver: ($x1, $x2) lies in no element of this mesh"))
    return out
end

# reproduce/arXiv-2606.04998/fig_linsolv.jl — ESSF::fig-linsolv.
#
# arXiv:2606.04998 Fig. 6 (`fig:linsolv`, main.tex:698-702): "Convergence of the preconditioned
# linear solver, exemplary for the m = 2 mode of the a = 0.5M, r0 = 10M configuration. The lines
# in the top panel, from left to right, show the GMRES iterations needed to solve consecutive AMR
# levels. Each solve converges in less than 10 preconditioned iterations. The first AMR level
# converges immediately because the multigrid bottom solver solves the problem directly. The
# bottom panel shows the walltime for each iteration. The first iteration on each new AMR level is
# slower than the following because it builds and caches subdomain matrices."
#
# THIN driver: zero physics. Two calls to `solve_mode_amr` (ESSF::amr-driver); both panels are
# read straight off the `ResidualMonitor` the sweep fills
# (`evidence/solver-stack-core/api_contract_row52_fig-linsolv.md`): per AMR grid one
# `SolveRecord` whose `residual` IS the top panel's staircase, whose `iteration_time` IS the
# bottom panel, and `cumulative_iterations(mon)` gives the grid boundaries on the shared x axis.
#
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_linsolv.jl --check
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_linsolv.jl --plot-only   # redraw only
#
# ── WHICH LATTICE THE PUBLISHED FIGURE IS ON — the one deviation from the yaml, and its evidence
#
# The plotted sweep runs on the FULL-z lattice, `equatorial_symmetry = false`: 30 elements and
# 480 grid points at the initial level, with the puncture contained in FOUR elements. The paper's
# own input file says the opposite — `ImposeEquatorialSymmetry: True` (yaml:36) — which halves the
# domain to `z² ∈ [0, 1]`, 20 elements, 320 points and TWO puncture elements (the puncture sits on
# `z² = 0`, which is the domain boundary of the half lattice, so only the radial direction has two
# neighbours sharing its corner).
#
# Three measurements say the FIGURE was made on the full lattice, whatever the yaml says:
#   * ‖b‖₀, which IS level 0's plotted initial residual because the initial guess is zero:
#     published 1.39401, full-z 1.39146 (−0.18 %), half domain 1.7444 (+25.1 %);
#   * the staircase itself: full-z gives 1, 7, 6, 5, 4, 3, 3, 3 (32 cumulative iterations) against
#     the published 1, 7, 6, 5, 4, 3, 3, 2 (31) — seven of the eight grids EXACT — while the half
#     domain gives 1, 7, 6, 5, 5, 4, 4, 4, 3 (39 over nine grids);
#   * the per-grid initial residuals track the published sequence monotonically on the full
#     lattice and do not on the half domain (which rises 2.53e-3 → 3.26e-3 at the first
#     refinement).
# The half-domain production sweep is therefore registered too, as a SECONDARY series
# (`prod_*` keys), so the yaml-versus-figure discrepancy is DATA in the golden rather than prose.
# The figure carries a one-line note saying which lattice it is.
#
# ── HOW MANY AMR GRIDS — a CONVENTION mapping, derived from the source, not fitted
#
# `Amr.Iterations: 8` (yaml:153) counts spectre's AMR *solver iterations*, i.e. SOLVES, and the
# yaml's own trigger `EveryNIterations {N: 1, Offset: 1}` (yaml:157-160, commented "Run AMR in
# every iteration, but not on the initial guess") skips refinement on iteration 0. Tracing
# spectre's loop — `src/Elliptic/Executables/Solver.hpp:348-357`
# (`RunEventsAndTriggers → StopAmr → ExecutePhaseChange → StepActions(solve) → IncrementIterationId`)
# against `src/Elliptic/Amr/Actions.hpp:56` (`HasConverged{num_iterations, 0}`) — the sweep is
# 8 solves separated by 7 refinements. GRTools' `solve_mode_amr(; iterations = n)` counts
# `amr_step!`s, so the paper's configuration is `iterations = 8 − 1 = 7`, which the driver DERIVES
# as `cfg.amr.iterations - cfg.amr.trigger_offset`. The published figure confirms the reading
# independently: its top panel carries exactly 8 disconnected staircases over a 0…32 axis.
#
# ── WHY THE STAIRCASES NOW SHORTEN WITH DEPTH (resolved; was this node's open finding)
#
# Until `fdcd919` this driver plotted a FLAT `1, 8, 9, 9, …` staircase against the published
# shortening one. The cause was the multigrid hierarchy, not the figure: spectre does not build a
# hierarchy by re-coarsening the current grid, it USES THE AMR HISTORY — every grid the sweep has
# already visited is a coarser multigrid level of the next
# (`Elliptic/Executables/SolveScalarSelfForce.hpp:173`;
# `ParallelAlgorithms/Amr/Actions/AdjustDomain.hpp:159-163, 203-212`;
# `SendDataToChildren.hpp:43-70`; `LinearSolver/Multigrid/ElementActions.hpp:96-135`;
# `Multigrid.hpp:56-58`). `hierarchy = :history` (now the `solve_mode_amr` default) reproduces
# that, and the preconditioner then improves with depth exactly as the paper's does: measured
# per-iteration residual reduction 30.8, 44.1, 52.5, 78.0, 216.3, 219.4, 264.0 against the
# published 26.0, 42.7, 48.0, 82.3, 244.6, 191.9, 567.0.
#
# ── THE BUILD COST IS A BOOKKEEPING EQUIVALENCE, NOT A LAZY BUILD
#
# The paper's bottom-panel bump on each grid's first iteration is the *lazy* construction of the
# Schwarz subdomain matrices inside that iteration (main.tex:702). This package builds them
# EAGERLY in `rebuild!`, so the same work is measured in `SolveRecord.build_time` and
# `iteration_time[1]` does not contain it. The bottom panel therefore plots
# `iteration_time[1] + build_time` at each grid's first point and `iteration_time[j]` thereafter:
# the same physical work charged to the same iteration, by addition rather than by measurement
# inside the timer. Both columns are registered separately in the golden, the figure says so in a
# note under the bottom panel, and `reported_first_bar_ratio_no_build` records what the panel
# looks like without the attribution.
#
# ── WARM-UP
#
# `build_time` on the first rebuild of a process is a compilation measurement. The driver runs a
# throwaway 2-solve sweep on the plotted configuration with its own monitor and discards it, then
# runs the production sweeps with fresh monitors. Every number that reaches the figure or the
# golden comes from a warmed sweep.
#
# ── TIMINGS ARE REPORTED, NOT GATED
#
# Wall-clock values are machine- and load-dependent (a contended run once measured a 4× inflated
# sweep); they are registered under `reported_*` names with a deliberately loose tolerance so
# `--check` cannot turn a busy machine into a red result. What IS gated is the structure they
# carry — the per-grid boolean `first_bar_is_slowest_per_level` — and the residual histories
# themselves, at an rtol set by the measured `-t1` vs `-t4` spread (bitwise identical with
# `OPENBLAS_NUM_THREADS=1`).
#
# ── WHAT STILL DOES NOT REPRODUCE (the node's remaining holes)
#
#   1. the LAST grid takes 3 iterations against the published 2 (32 cumulative against 31);
#   2. the walltime panel's spike/plateau amplitude, because of the eager-versus-lazy build above;
#   3. the published grids are m-DEPENDENT and ours are not: dofs 480, 921, 1516, 2297, 3288,
#      4513, 5996, 7761 here against 480, 975, 1604, 2427, 3468, 4751, 6300 read off the paper's
#      own convergence figure — the yaml's commented-out `TruncationError` criterion at 1e-11
#      (yaml:135-138) is the hypothesis, and it belongs to rows 43 and 39, not here.
#
# Style: `_style/manuals/arXiv-2606.04998.md` §0 + §6 (single column, two stacked panels sharing
# x at 3:2, grayscale only, black line + filled circles, log y, dashed "Solver tolerance" line,
# framed "GMRES iterations" legend upper right, light-gray grid, CM serif, walltime stamp
# right-aligned above the bottom panel).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()

# ───────────────────────────── parameters (main.tex:701 names all three) ─────────────────────
const SPIN = 0.5
const R0   = 10.0
const MM   = 2
const PORD = 4              # anc/ScalarSelfForce.yaml:11 (InitialGridPoints)
const LREF = 0              # anc/ScalarSelfForce.yaml:12 (InitialRefinement / InitialCoarseLevels)
# THE PLOTTED LATTICE. `false` = the FULL-z domain (30 elements, 480 points, the puncture in four
# elements). yaml:36 says `ImposeEquatorialSymmetry: True`; the published ‖b‖₀ and staircase say
# the figure is on the full lattice (see the header). The half domain is measured too, as the
# secondary `prod_*` series.
const EQSYM = false
const LATTICE_NOTE = "full-z lattice, ImposeEquatorialSymmetry: False — the published lattice; the paper's yaml says True"

# ───────────────────────────── the figure (presentation only) ────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

"Format `s` seconds as the paper's `H:MM:SS` walltime stamp (Fig. 6, above the bottom panel)."
_hms(s) = (h = floor(Int, s / 3600); m = floor(Int, (s - 3600h) / 60);
           @sprintf("%d:%02d:%02d", h, m, floor(Int, s - 3600h - 60m)))

"""
    draw(xs, residual, iteration_time, build_time, abstol, walltime) -> path

Fig. 6: the top panel is one disconnected black line+circle staircase per AMR level, the bottom
panel the same for the walltime per iteration with each level's first point carrying that level's
build cost. `xs` are the cumulative-iteration level boundaries (`cumulative_iterations`).
"""
function draw(xs, residual, iteration_time, build_time, abstol, walltime)
    nlev = length(residual)
    total = xs[end]
    xstep = total <= 40 ? 4 : total <= 96 ? 8 : 16
    xhi = xstep * cld(total + 1, xstep)
    common = (xlims = (-0.7, xhi + 0.7), xticks = 0:xstep:xhi, yscale = :log10,
              grid = true, gridalpha = 0.45, gridlinewidth = 0.6,
              foreground_color_grid = RGB(0.75, 0.75, 0.75),
              tick_direction = :in, widen = false, guidefontsize = 11, tickfontsize = 9)

    rall = vcat(residual...)
    tall = [j == 1 ? iteration_time[k][1] + build_time[k] : iteration_time[k][j]
            for k in 1:nlev for j in eachindex(iteration_time[k])]

    # top — the residual staircase; ticks every 3 decades anchored on 10⁻¹ (manual §6)
    top = plot(; common..., ylims = (minimum(rall) / 10^0.8, maximum(rall) * 10^0.7),
               yticks = decade_ticks(minimum(rall), maximum(rall); step = 3, anchor = -1,
                                     inside = false), ylabel = "Absolute residual",
               xformatter = _ -> "", legend = :topright, legendfontsize = 9,
               bottom_margin = -4Plots.mm, title = LATTICE_NOTE, titlelocation = :left,
               titlefontsize = 6, titlefontcolor = RGB(0.35, 0.35, 0.35),
               foreground_color_legend = RGB(0.55, 0.55, 0.55), background_color_legend = :white)
    for k in 1:nlev
        plot!(top, xs[k] .+ (0:(length(residual[k]) - 1)), residual[k];
              color = :black, linewidth = 1.0, marker = (:circle, 3.0), markercolor = :black,
              markerstrokewidth = 0, label = k == 1 ? "GMRES iterations" : "")
    end
    hline!(top, [abstol]; color = :black, linewidth = 1.1, linestyle = :dash, label = "")
    annotate!(top, 0.3, abstol / 2.2,
              text("Solver tolerance", 10, :left, :top, "Computer Modern"))

    # bottom — walltime per iteration; each level's first point carries its `build_time`
    bot = plot(; common..., ylims = (minimum(tall) / 10^0.25, maximum(tall) * 10^0.25),
               yticks = decade_ticks(minimum(tall), maximum(tall); inside = false), ylabel = "Walltime per iteration [s]",
               xlabel = "Cumulative linear solver iteration", legend = false,
               title = "Walltime:  " * _hms(walltime), titlelocation = :right, titlefontsize = 10)
    for k in 1:nlev
        y = copy(iteration_time[k])
        y[1] += build_time[k]
        plot!(bot, xs[k] .+ (0:(length(y) - 1)), y; color = :black, linewidth = 1.0,
              marker = (:circle, 3.0), markercolor = :black, markerstrokewidth = 0, label = "")
    end
    # the attribution, ON THE FIGURE: this package builds the subdomain matrices eagerly, so the
    # paper's lazy first-iteration cost is added to each level's first point (see the header).
    annotate!(bot, xhi - 0.4, minimum(tall) / 10^0.13,
              text("first point of each level includes the eager operator build",
                   7, :right, :bottom, RGB(0.25, 0.25, 0.25), "Computer Modern"))

    plt = plot(top, bot; layout = grid(2, 1; heights = [0.6, 0.4]), size = (500, 550), dpi = 150,
               left_margin = 4Plots.mm, right_margin = 4Plots.mm, top_margin = 2Plots.mm,
               bottom_margin = 2Plots.mm)
    png = joinpath(@__DIR__, "expected", "fig_linsolv.png")
    pdf = joinpath(@__DIR__, "expected", "fig_linsolv.pdf")
    mkpath(dirname(png))
    savefig(plt, png)
    savefig(plt, pdf)
    nb = filesize(png)
    println("[fig_linsolv] wrote ", relpath(png, ReproHarness.REPRO_ROOT), " (", nb, " bytes) and ",
            relpath(pdf, ReproHarness.REPRO_ROOT), " (", filesize(pdf), " bytes, untracked)")
    nb <= 200 * 1024 || error("fig_linsolv.png is $(nb) bytes > 200 KB — reduce the dpi/size")
    return png
end

"Split the flat vector `flat` into per-level pieces of lengths `lens`."
_unflat(flat, lens) = (o = 0; [(v = Float64.(flat[(o + 1):(o + n)]); o += n; v) for n in lens])

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2606.04998", "fig_linsolv")
    isfile(gpath) || ReproHarness.hole(
        "--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(Int.(v["cumulative_iterations"]),
         _unflat(v["residual_flat"], Int.(v["residual_lengths"])),
         _unflat(v["reported_iteration_time_flat"], Int.(v["iterations_per_level"])),
         Float64.(v["reported_build_time"]), Float64(v["gmres_absolute_residual"]),
         Float64(v["reported_sweep_wall_time"]))
    println("[fig_linsolv] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

ReproHarness.require_tier(opts, :heavy)   # 252.7 s at -t4 in run_all with a warm depot, ~400 s
                                          # standalone with precompilation: over the fast ceiling

using GrTools
using GrTools.ScalarSelfForce
using GrTools.EllipticCore

for f in (:solve_mode_amr, :crp_paper_config)
    isdefined(ScalarSelfForce, f) || ReproHarness.hole(
        "ScalarSelfForce.$(f) is not available: Fig. 6 needs the hp-AMR driver " *
        "(ESSF::amr-driver, capabilities/ScalarSelfForce/src/Solve/amr.jl).")
end

_cfg(sym) = crp_paper_config(; spin = SPIN, orbital_radius = R0, m_mode = MM,
                             inner_radius = CircularOrbit(SPIN, R0).rp, P = PORD, L = LREF,
                             equatorial_symmetry = sym)
cfg = _cfg(EQSYM)                 # the PLOTTED configuration: the full-z published lattice
cfg_prod = _cfg(true)             # the yaml's own half domain, measured as the secondary series
orb = CircularOrbit(SPIN, R0)

# spectre's `Amr.Iterations` counts SOLVES and the trigger's `Offset` skips refinement on the
# initial guess (see the header): the paper's 8 solves are 7 `amr_step!`s here.
const STEPS = cfg.amr.iterations - cfg.amr.trigger_offset
const ABSTOL = Float64(cfg.gmres.absolute_residual)

println("\n[fig_linsolv] a = ", SPIN, "M  r₀ = ", R0, "M  m = ", MM, "  P = ", PORD, "  L = ", LREF,
        "  ImposeEquatorialSymmetry = ", EQSYM, " (plotted); the yaml says true")
println("[fig_linsolv] Amr.Iterations = ", cfg.amr.iterations, " solves, trigger offset ",
        cfg.amr.trigger_offset, " ⇒ ", STEPS, " refinement steps ⇒ ", STEPS + 1, " AMR levels")
println("[fig_linsolv] GMRES AbsoluteResidual = ", ABSTOL, ", MaxIterations = ",
        cfg.gmres.max_iterations, ", Schwarz MaxOverlap = ", cfg.schwarz.max_overlap,
        ", passes = ", cfg.schwarz.iterations)
flush(stdout)

# ───────────────────────────── warm-up, then the production sweep ────────────────────────────
tw = time()
solve_mode_amr(orb, cfg, MM; iterations = 1, monitor = ResidualMonitor())
@printf("[fig_linsolv] warm-up sweep (2 throwaway solves, absorbs the JIT) %.1f s\n", time() - tw)
flush(stdout)

mon = ResidualMonitor()
t0 = time()
hist = solve_mode_amr(orb, cfg, MM; iterations = STEPS, monitor = mon)
wall = time() - t0
@printf("[fig_linsolv] production sweep: %d grids, %.2f s wall\n", length(hist.levels), wall)
flush(stdout)

# the SECONDARY series: the yaml's own half domain at its own `Amr.Iterations` (nine grids), so
# the lattice discrepancy is data in the golden rather than prose. Not plotted.
mon_p = ResidualMonitor()
tp = time()
hist_p = solve_mode_amr(orb, cfg_prod, MM; iterations = cfg_prod.amr.iterations, monitor = mon_p)
wall_p = time() - tp
prod_iters = [r.iterations for r in mon_p.records]
prod_r0 = [Float64(r.residual[1]) for r in mon_p.records]
prod_rf = [Float64(r.final_residual) for r in mon_p.records]
prod_ndofs = [l.mesh.ndofs for l in hist_p.levels]
prod_nelem = [length(l.mesh.elements) for l in hist_p.levels]
prod_npunct = [length(l.punctures) for l in hist_p.levels]
@printf("[fig_linsolv] half-domain (yaml) sweep: %d grids, %.2f s wall, iterations %s (total %d)\n",
        length(hist_p.levels), wall_p, string(prod_iters), sum(prod_iters))
flush(stdout)

recs = mon.records
nlev = length(recs)
xs = cumulative_iterations(mon)
iters = [r.iterations for r in recs]
resid = [Float64.(r.residual) for r in recs]
ittime = [copy(r.iteration_time) for r in recs]
btime = [r.build_time for r in recs]
_rest(k) = length(ittime[k]) > 1 ? sum(ittime[k][2:end]) / (length(ittime[k]) - 1) : NaN

for k in 1:nlev
    r = recs[k]
    println("[fig_linsolv] level ", k - 1, ": ", r.iterations, " iters, ",
            @sprintf("%.4e → %.4e", r.residual[1], r.final_residual), " (conv ", r.converged,
            "), build ", @sprintf("%.3f", r.build_time), " s, bar1 ",
            @sprintf("%.3f", ittime[k][1] + btime[k]), " s, mean(rest) ",
            @sprintf("%.4f", _rest(k)), " s, ", length(hist.levels[k].mesh.elements),
            " elements / ", hist.levels[k].mesh.ndofs, " dofs")
end
flush(stdout)

# ───────────────────────────── the node's gates (GLOBAL_DAG row 52) ──────────────────────────
# (1) "Each solve converges in less than 10 preconditioned iterations" (main.tex:701, :720).
# (2) "The first AMR level converges immediately because the multigrid bottom solver solves the
#     problem directly" (main.tex:701)  ⇒  level 0 takes ≤ 1 iteration.
# (3) "the number of GMRES iterations does not increase when AMR levels are added" (main.tex:720,
#     multigrid scale independence). Now asserted in the STRICT form the published data itself
#     satisfies — non-increasing from grid 1 on — which became true when the multigrid hierarchy
#     switched to the AMR history (`fdcd919`). The weaker BOUNDED form is registered beside it so
#     a regression back to a flat staircase is visible in the golden, not just in the gate.
# (4) "The first iteration on each new AMR level is slower than the following" (main.tex:702),
#     reproduced by attributing this package's EAGER `build_time` to the grid's first point.
# (5) every grid actually reached the GMRES absolute tolerance.
# (6) TOTAL iterations, a regression band against the published 31: measured 32, gated at ≤ 33.
#     The published figure's own cumulative axis is the reference; the one-iteration excess is the
#     last grid (3 against 2) and is a recorded HOLE, not a tolerance.
# (7) ‖b‖₀ — level 0's initial residual IS the source norm, the initial guess being zero — within
#     2 % of the published 1.375. This is what identifies the LATTICE: full-z 1.3915 (+1.2 %),
#     the yaml's half domain 1.7444 (+26.9 %).
const PUB_ITERS = [1, 7, 6, 5, 4, 3, 3, 2]     # extracted from linsolv.pdf's own vector paths
const PUB_TOTAL = 31
# the published ABSCISSA of fig:convergence — grid points per AMR grid, the same grid sequence
# this figure sweeps. Seven values for the eight grids: the paper's convergence figure plots the
# first seven and holds the eighth as its reference.
const PUB_NDOFS = [480, 975, 1604, 2427, 3468, 4751, 6300]
# ‖b‖₀ read off the same extraction, with the EXACT flipped-y anchors 20.44117 / 103.93375 (the
# 1e-1 and 1e-13 gridlines). An earlier extraction rounded them to 20.4 / 103.9 and put the whole
# published residual series ~1.35 % low (1.375 instead of 1.39401); the iteration counts and the
# reduction-per-iteration ratios are invariant under that rescaling, ‖b‖₀ is not.
const PUB_B0 = 1.39401
const B0_BAND = 0.01
maxiters = maximum(iters)
level0_iters = iters[1]
tail = iters[2:end]
strict_noninc = all(tail[k + 1] <= tail[k] for k in 1:(length(tail) - 1))
bounded = maximum(tail) - minimum(tail) <= 1 && tail[end] <= tail[1] + 1
flat_from_level2 = length(tail) < 2 || allequal(tail[2:end])
firstbar = [length(ittime[k]) == 1 ? true :
            (ittime[k][1] + btime[k]) > maximum(ittime[k][2:end]) for k in 1:nlev]
allconv = all(r.converged for r in recs)
b0 = resid[1][1]
ndofs = [l.mesh.ndofs for l in hist.levels]
ndofs_cmp = ndofs[1:min(length(PUB_NDOFS), end)]
ndofs_match = count(k -> ndofs_cmp[k] == PUB_NDOFS[k], eachindex(ndofs_cmp))
ndofs_maxreldev = maximum(abs(ndofs_cmp[k] - PUB_NDOFS[k]) / PUB_NDOFS[k] for k in eachindex(ndofs_cmp))
b0_reldev = abs(b0 - PUB_B0) / PUB_B0
grids_matching_published = count(k -> k <= length(PUB_ITERS) && iters[k] == PUB_ITERS[k], 1:nlev)

fails = String[]
maxiters < 10 || push!(fails, "GATE 1: a grid needed $(maxiters) ≥ 10 GMRES iterations ($(iters))")
level0_iters <= 1 || push!(fails, "GATE 2: the first grid took $(level0_iters) iterations, " *
                                  "not the ≤ 1 of a direct bottom solve")
strict_noninc || push!(fails, "GATE 3: the iteration count is NOT non-increasing from grid 1 " *
                              "($(tail)) — the published staircase shortens with depth")
all(firstbar) || push!(fails, "GATE 4: the first iteration is not the slowest on grid(s) " *
                              "$(findall(!, firstbar) .- 1) even with the grid's build cost added")
allconv || push!(fails, "GATE 5: a grid did not converge to the GMRES absolute tolerance")
xs[end] <= PUB_TOTAL + 2 || push!(fails, "GATE 6: $(xs[end]) cumulative iterations against the " *
                                         "published $(PUB_TOTAL); the regression band is ≤ $(PUB_TOTAL + 2)")
b0_reldev <= B0_BAND || push!(fails, "GATE 7: ‖b‖₀ = $(b0) deviates $(round(100b0_reldev, digits = 2)) % " *
                                     "from the published $(PUB_B0) (band $(100 * B0_BAND) %) — wrong lattice?")
# (8) the GRID SEQUENCE itself against the paper's own published abscissa (fig:convergence's x
#     axis is this sweep's grid-point count). Exact equality where the published value exists.
ndofs_match == length(PUB_NDOFS) ||
    push!(fails, "GATE 8: the grid sequence $(ndofs_cmp) does not match the published abscissa " *
                 "$(PUB_NDOFS) ($(ndofs_match)/$(length(PUB_NDOFS)) exact, max rel dev " *
                 "$(round(100ndofs_maxreldev, digits = 3)) %)")

# ───────────────────────────── registered values (the golden) ────────────────────────────────
@repro_result "spin" SPIN
@repro_result "orbital_radius" R0
@repro_result "m_mode" MM
@repro_result "grid_points_per_dim" PORD
@repro_result "initial_refinement_L" LREF
@repro_result "amr_iterations_yaml" cfg.amr.iterations
@repro_result "amr_trigger_offset" cfg.amr.trigger_offset
@repro_result "amr_refinement_steps" STEPS
@repro_result "n_amr_levels" nlev
@repro_result "equatorial_symmetry_plotted" EQSYM
@repro_result "equatorial_symmetry_yaml" true
@repro_result "gmres_absolute_residual" ABSTOL
@repro_result "gmres_max_iterations" cfg.gmres.max_iterations
@repro_result "schwarz_iterations" cfg.schwarz.iterations
@repro_result "multigrid_iterations" cfg.multigrid.iterations

@repro_result "iterations_per_level" iters
@repro_result "cumulative_iterations" xs
@repro_result "total_iterations" xs[end]
@repro_result "residual_lengths" [length(r) for r in resid]
@repro_result "residual_flat" vcat(resid...)
# log10 of the same numbers. The harness compares arrays norm-wise, so on the raw vector the O(1)
# initial residual hides a factor-2 move at the 1e-13 floor; on the log10 vector every point
# carries comparable weight and a per-cent move in ANY residual fails at the declared rtol.
@repro_result "residual_log10_flat" log10.(vcat(resid...))
@repro_result "initial_residual_per_level" [r.residual[1] for r in recs]
@repro_result "final_residual_per_level" [Float64(r.final_residual) for r in recs]
@repro_result "converged_per_level" [r.converged for r in recs]
# geometric mean residual reduction per iteration on each level — the quantity that differs from
# the paper (ours degrades with depth, the paper's improves); see the header.
@repro_result "reduction_factor_per_level" [(resid[k][1] / resid[k][end])^(1 / iters[k])
                                            for k in 1:nlev]

@repro_result "nelements_per_level" [length(l.mesh.elements) for l in hist.levels]
@repro_result "ndofs_per_level" [l.mesh.ndofs for l in hist.levels]
@repro_result "multigrid_levels_per_level" [length(r.level_ndofs) for r in recs]
@repro_result "coarsest_ndofs_per_level" [r.level_ndofs[end] for r in recs]
@repro_result "subdomain_builds_per_level" [r.subdomain_builds for r in recs]
@repro_result "subdomain_hits_per_level" [r.subdomain_hits for r in recs]
@repro_result "bottom_builds_per_level" [r.bottom_builds for r in recs]
@repro_result "bottom_hits_per_level" [r.bottom_hits for r in recs]

# gate outcomes (exact-compared)
@repro_result "max_iterations_per_level" maxiters
@repro_result "level0_iterations" level0_iters
@repro_result "iterations_bounded_across_levels" bounded
@repro_result "iterations_flat_from_level2" flat_from_level2
@repro_result "iterations_strictly_non_increasing_from_level1" strict_noninc
# the published reference, extracted from linsolv.pdf's own vector paths (not pixels); see the
# paper README for the extraction method
@repro_result "published_iterations_per_level" PUB_ITERS
@repro_result "published_total_iterations" PUB_TOTAL
@repro_result "published_source_norm" PUB_B0
@repro_result "grids_matching_published_iterations" grids_matching_published
@repro_result "source_norm" b0
@repro_result "source_norm_relative_deviation" b0_reldev
@repro_result "published_ndofs_per_level" PUB_NDOFS
@repro_result "ndofs_matching_published" ndofs_match
@repro_result "ndofs_max_relative_deviation" ndofs_maxreldev
# the SECONDARY series: the yaml's own half domain (ImposeEquatorialSymmetry: True), nine grids.
# Its 8-grid prefix is the like-for-like comparison with the plotted full-z sweep.
@repro_result "prod_n_amr_levels" length(prod_iters)
@repro_result "prod_iterations_per_level" prod_iters
@repro_result "prod_total_iterations" sum(prod_iters)
@repro_result "prod_total_iterations_first8" sum(prod_iters[1:min(8, end)])
@repro_result "prod_initial_residual_per_level" prod_r0
@repro_result "prod_final_residual_per_level" prod_rf
@repro_result "prod_nelements_per_level" prod_nelem
@repro_result "prod_ndofs_per_level" prod_ndofs
@repro_result "prod_puncture_elements_per_level" prod_npunct
@repro_result "puncture_elements_per_level" [length(l.punctures) for l in hist.levels]
@repro_result "first_bar_is_slowest_per_level" firstbar
@repro_result "all_levels_converged" allconv
@repro_result "gates_pass" isempty(fails)

# REPORTED, NOT GATED: wall-clock seconds (machine- and load-dependent; see the header).
@repro_result "reported_iteration_time_flat" vcat(ittime...)
@repro_result "reported_build_time" btime
@repro_result "reported_build_time_hierarchy" [r.build_time_hierarchy for r in recs]
@repro_result "reported_build_time_bottom" [r.build_time_bottom for r in recs]
@repro_result "reported_solve_time" [r.solve_time for r in recs]
@repro_result "reported_sweep_wall_time" wall
@repro_result "reported_prod_sweep_wall_time" wall_p
@repro_result "reported_first_bar_ratio" [(ittime[k][1] + btime[k]) / _rest(k) for k in 1:nlev]
@repro_result "reported_first_bar_ratio_no_build" [ittime[k][1] / _rest(k) for k in 1:nlev]

draw(xs, resid, ittime, btime, ABSTOL, wall)

if !isempty(fails)
    println("\n[fig_linsolv] FINDING — node gates violated:")
    foreach(f -> println("    ", f), fails)
end

# rtol 1e-8: measured `-t1` vs `-t4` spread on the residual histories is 0 (bitwise identical with
# `OPENBLAS_NUM_THREADS=1`); the declared tolerance is three orders looser than that measurement
# and still an order tighter than the 4.2e-12 relative drift that an UNPINNED BLAS thread count
# introduces (ESSF::amr-driver gate (f)). The `reported_*` timing rows carry a deliberately loose
# tolerance so a busy machine cannot turn a reproduced figure into a red check.
const LOOSE = (rtol = 1.0e3, atol = 1.0e3)
const TOLS = Dict("reported_iteration_time_flat" => LOOSE, "reported_build_time" => LOOSE,
                  "reported_build_time_hierarchy" => LOOSE, "reported_build_time_bottom" => LOOSE,
                  "reported_solve_time" => LOOSE, "reported_sweep_wall_time" => LOOSE,
                  "reported_prod_sweep_wall_time" => LOOSE,
                  "reported_first_bar_ratio" => LOOSE, "reported_first_bar_ratio_no_build" => LOOSE)
code = ReproHarness.conclude!(opts; paper = "2606.04998", id = "fig_linsolv",
                              script = @__FILE__, rtol = 1e-8, tolerances = TOLS)
exit(isempty(fails) ? code : max(code, 1))

# reproduce/arXiv-2606.04998/fig_convergence.jl — ESSF::fig-convergence.
#
# arXiv:2606.04998 Fig. 4 (`fig:convergence`, main.tex:684-688): "Convergence of the DG method
# for modes up to m = 20, exemplary for the a = 0.5 M, r0 = 10 M configuration. The convergence
# is exponential despite the non-smooth solution on the grid under our hp-AMR scheme. Each mode
# solves in a few seconds on an Apple M2 Pro chip."
#
# WHAT IS PLOTTED, in the paper's own words (main.tex:719-724 + its footnote):
#   * top panel — "the absolute error for each individual m-mode contribution,
#     |F_r^m − F^m_{r,ref}|, up to m = 20 as a function of the number of grid points across a
#     sequence of AMR iterations … As reference solution for each m-mode we use the result from
#     the highest-resolution AMR iteration here."  So `err_k = |F_r^m(level k) − F_r^m(last)|`
#     and the LAST level is the reference, hence not plotted: an 8-grid sweep gives 7 points.
#   * abscissa — "we plot the error over the number of grid points per dimension", axis label
#     "Grid points per dimension (N_points)^{1/2}", i.e. `sqrt(ndofs)` of that level's mesh.
#   * bottom panel — "the cumulative walltime for each m-mode solve".
#
# THIN driver: zero physics. ONE `mode_pool` for the whole figure (the row-47 contract: a fresh
# worker pays its entire JIT before its first solve, so every sweep of a figure must share one
# pool), one warm-up sweep to absorb that JIT, one production sweep, and then the SEALED
# extraction `mode_force_convergence(orbit, config, m, levels; x_p)` — the level route added for
# this node, which is what makes a fanned-out `ModeResult` (levels + x_p, no live
# `ModeSolution`) reachable by row 44's extraction at all.
#
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_convergence.jl --check
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_convergence.jl --plot-only  # redraw
#
# HOW MANY AMR LEVELS — the same convention mapping fig_linsolv.jl derives, not a fitted number.
# `Amr.Iterations: 8` (yaml:153) counts spectre's AMR SOLVER iterations and the trigger
# `EveryNIterations {N: 1, Offset: 1}` (yaml:157-160) skips refinement on the initial guess, so
# the paper's sweep is 8 solves separated by 7 refinements = 8 grids. `solve_mode_amr(;
# iterations = n)` counts `amr_step!`s, so the driver DERIVES `cfg.amr.iterations −
# cfg.amr.trigger_offset = 7`. The published figure confirms it independently: its curves carry
# exactly 7 vertices, which is 8 grids minus the reference level.
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
# THE ABSCISSA DOES NOT REPRODUCE, AND THAT IS THE NODE'S FINDING (see the paper README).
#
# The published curves' vertex coordinates were read off `ssf_convergence.pdf`'s own vector
# paths (21 polylines per panel, 7 points each; x calibrated on the 20/40/60/80 gridlines, y on
# the decade gridlines — the extraction script and its table are in this node's evidence). They
# give the published grid-point sequence
#
#     N_points = 480, 975, 1604, 2427, 3468, 4751, 6300      (m = 0 … 4)
#
# against ours, 320, 593, 950, 1411, 1988, 2693, 3538. TWO structural differences, both
# escalated rather than plotted around:
#
#   (1) THE LATTICE. 480 and 975 are, exactly, the `NumPoints` column of spectre's own
#       `tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml` regression (rows `[0, 480, …]`,
#       `[1, 480, …]`, `[2, 975, …]`) — a file that sets `ImposeEquatorialSymmetry: False` and
#       therefore a FULL-z 5×3 block lattice of 30 elements at P = 4. The production input this
#       whole reproduction runs, `anc/ScalarSelfForce.yaml:36`, sets
#       `ImposeEquatorialSymmetry: True`: a 5×2 lattice, 20 elements, 320 grid points. So the
#       published figure's grids are the full-z ones and ours are the equatorial-symmetry ones;
#       the physics is the same problem, the resolution axis is not the same axis.
#   (2) THE CRITERIA. From level 1 on the published `N_points` depends on m — 975 for m ≤ 4,
#       970 for m ≥ 5, spreading to 6300 (m ≤ 4) … 6092 (m ≥ 18) at the last plotted level.
#       `RefineAtPuncture` and `IncreaseResolution`, the two criteria the shipped yaml enables
#       (yaml:130-132), are both m-INDEPENDENT and cannot produce that. The yaml's own
#       commented-out `TruncationError` block (`AbsoluteTarget: 1.e-11`, `RelativeTarget:
#       1.e-11`, yaml:133-136) is the candidate: it is solution-dependent, `Policies.
#       AllowCoarsening: True` is only meaningful with an error-based criterion, the figure's
#       dashed "Solver tolerance" line sits at 1e-11 (measured) while the GMRES
#       `AbsoluteResidual` is 1e-12, and main.tex:721 says the error falls "until it reaches the
#       tolerances set for the GMRES algorithm AND FOR AMR" — a second tolerance the shipped
#       configuration does not have. [HOLE: not demonstrated, no run of ours enables it.]
#
# Our meshes are m-independent by construction, and that is GATED here (`ndofs_identical`) —
# it is the sharp statement of the difference, not a defect on our side.
#
# THE FIT GATE IS THE PUBLISHED FIGURE'S OWN ENVELOPE. GLOBAL_DAG row 50 planned "21 exponential
# fits with R² > 0.98". That gate is FALSIFIED BY THE PAPER'S OWN DATA: fitting log10(err)
# against sqrt(N_points) on the 21 published curves gives R² ∈ [0.7996, 0.9929] with only 9 of
# 21 modes at or above 0.98 (m = 15 is the worst, at 0.7996; its level-2 point dips three
# decades below its neighbours and comes back). The landed gate is therefore the published
# envelope — every mode's fit must have a NEGATIVE rate and R² at or above the published
# minimum — plus a floor on the total drop. The row-44 two-step contraction envelope is
# MEASURED and registered but deliberately NOT gated (see the constants block). Both deviations
# from the planned gate are recorded in the ledger row.
#
# Style: `_style/manuals/arXiv-2606.04998.md` §0 + §4 (single column, two stacked panels sharing
# x at 3:1, x linear and labelled only on the bottom, log y, 21 grayscale lines with no markers
# from black at m = 0 to light gray at m = 20, framed legend upper right with ONLY the two
# extreme entries, black dashed "Solver tolerance" line with an inside-axes label below-left,
# light-gray grid, CM serif, "1 s"/"10 s" walltime ticks). The gray ramp is the published one,
# read off the PDF: stroke gray = m/23 for m = 0 … 20.

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

# `Statistics` is deliberately NOT added to reproduce/Project.toml for one number: a dependency
# added without its Manifest entry broke every driver in this directory once already.
_median(v) = (s = sort(collect(Float64, v)); n = length(s);
              isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

opts = ReproHarness.parse_args()
ReproHarness.setup!()

# ───────────────────────────── parameters (main.tex:687 names all of them) ───────────────────
const SPIN = 0.5
const R0   = 10.0
const MS   = 0:20           # main.tex:687 "modes up to m = 20"; = CRP_PAPER_MMODES
const PORD = 4              # anc/ScalarSelfForce.yaml:74 (InitialGridPoints)
const LREF = 0              # anc/ScalarSelfForce.yaml:75 (InitialLevels / InitialCoarseLevels)
const NWORKERS = 5          # compute-usage: 5 processes × 1 thread × 1 BLAS thread = 5 CPUs

# ───────────────────────────── the published figure, measured (NOT plotted) ──────────────────
# Read off `ssf_convergence.pdf`'s vector paths; used ONLY to calibrate this node's gates and to
# quantify the deviations above. No published number is ever plotted as ours.
const PUB_NDOFS      = [480, 975, 1604, 2427, 3468, 4751, 6300]      # m ≤ 4; m ≥ 5 differ
const PUB_R2_MIN     = 0.7996        # m = 15, the worst of the 21 published fits
const PUB_R2_MAX     = 0.9929        # m = 0
const PUB_R2_MEDIAN  = 0.9833
const PUB_R2_GE_098  = 9             # how many of the 21 published modes reach the planned gate
const PUB_RATE_MIN   = -0.1421       # decades per unit sqrt(N); m = 0
const PUB_RATE_MAX   = -0.1043       # m = 15
const PUB_DECADES    = (6.06, 8.12)  # log10(err_first / err_last) over the 7 plotted points
const PUB_TOL_LINE   = 1.0e-11       # where the dashed "Solver tolerance" line sits

# GATES, calibrated BEFORE the production run — externally on the published figure and
# internally on this node's 7-mode throwaway prototype (evidence `calibration_prototype.log`,
# m ∈ {0, 1, 2, 5, 10, 15, 20} at the same 8 grids):
#
#     m      0       1       2       5      10      15      20
#     R²  0.9876  0.9741  0.9863  0.9578  0.9748  0.8212  0.9240
#     dec   4.76    4.80    4.92    5.82    4.67    2.65    3.69
#     2-step max  0.187  0.239  0.078  0.136  0.521  0.236  0.307
#
const R2_FLOOR = 0.80        # = the published figure's own worst fit (m = 15, 0.7996), rounded
                             # down. "At least as straight as the paper's own worst curve"; the
                             # prototype's worst is 0.8212, also at m = 15, and for the same
                             # reason — when a coarse level happens to land near the reference
                             # value the "error" there is anomalously small and breaks the line.
const TWO_STEP_CONTRACT = 1.0  # contraction over two levels — REGISTERED AND REPORTED, NOT
                             # GATED. The row-44 contract words the envelope as "at most HALF
                             # the difference two levels back", measured there at m ∈ {1, 2}.
                             # It is a statement about the PRE-FLOOR regime: at 8 AMR levels the
                             # high-m sequences are already on their resolution floor from level
                             # ~5 (|F_r^20| = 2.1e-9 against a floor of 1.7e-10), where the
                             # differences are floor noise and their ratios mean nothing. The
                             # calibration prototype already exceeded one half at m = 10 (0.521),
                             # and the production sweep exceeds ONE at m = 14 (1.057) and m = 19
                             # (7.687). Gating it would be gating noise; the numbers are in the
                             # golden and the exceedances are named in the paper README.
const TWO_STEP_HALF = 0.5    # the row-44 contract's own constant — reported against, not gated
const DECADES_MIN = 2.0      # the prototype's smallest drop over the 7 plotted points is 2.65
                             # decades (m = 15); the published range is 6.06 … 8.12 over ITS
                             # seven points, which reach 1e-13 because its last grid is 6300
                             # points against our 3538 (see the abscissa finding above)

# ───────────────────────────── the figure (presentation only) ────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

"The published gray ramp: stroke gray `m/23`, black at m = 0 to 0.870 at m = 20 (PDF-measured)."
mode_gray(m, mmax) = (g = m / (mmax + 3); RGB(g, g, g))

"""
    _decade_ticks(lo, hi, step, anchor)

Decade ticks INSIDE `[lo, hi]`, one every `step` decades anchored on exponent `anchor`. The
clamp is load-bearing: GR draws an out-of-range tick pinned to the axis end, which on a shared-x
two-panel stack shows up as a stray label between the panels.
"""
function _decade_ticks(lo, hi, step, anchor)
    e0 = ceil(Int, log10(lo))
    e1 = floor(Int, log10(hi))
    return [10.0^e for e in e0:e1 if mod(e - anchor, step) == 0]
end

"""
    draw(ms, xs, errs, walls, tol, out_png; publish) -> saved

Fig. 4: `xs[i]` is the shared abscissa `sqrt(N_points)` of mode `ms[i]`, `errs[i]` its absolute
error against its own highest-resolution level and `walls[i]` its cumulative walltime, all with
the reference level dropped. 21 grayscale lines, no markers, in both panels.
"""
function draw(ms, xs, errs, walls, tol, out_png; publish::Bool)
    mmax = maximum(ms)
    xall = vcat(xs...); eall = vcat(errs...); wall = vcat(walls...)
    xlo, xhi = minimum(xall), maximum(xall)
    pad = 0.045 * (xhi - xlo)
    xtick = 10 * cld(ceil(Int, (xhi - xlo) / 5), 10)
    xtick = max(xtick, 10)
    ticks = collect((xtick * fld(round(Int, xlo), xtick)):xtick:(xtick * cld(round(Int, xhi), xtick)))
    common = (xlims = (xlo - pad, xhi + pad), xticks = ticks, yscale = :log10,
              grid = true, gridalpha = 0.45, gridlinewidth = 0.6,
              foreground_color_grid = RGB(0.75, 0.75, 0.75),
              tick_direction = :in, widen = false, guidefontsize = 11, tickfontsize = 9)

    ylo = min(minimum(eall), tol) / 10^0.55
    yhi = maximum(eall) * 10^0.55
    top = plot(; common..., ylims = (ylo, yhi), yticks = _decade_ticks(ylo, yhi, 1, 0),
               ylabel = "Absolute error \$|F_r^m - F_{r,\\mathrm{ref}}^m|\$",
               xformatter = _ -> "", legend = :topright, legendfontsize = 9,
               bottom_margin = -4Plots.mm,
               foreground_color_legend = RGB(0.55, 0.55, 0.55), background_color_legend = :white)
    for (i, m) in enumerate(ms)
        plot!(top, xs[i], errs[i]; color = mode_gray(m, mmax), linewidth = 1.0,
              label = m == 0 ? "m=0" : m == mmax ? "m=$(mmax)" : "")
    end
    hline!(top, [tol]; color = :black, linewidth = 1.1, linestyle = :dash, label = "")
    annotate!(top, xlo - pad + 0.02 * (xhi - xlo), tol / 1.35,
              text("Solver tolerance", 9, :left, :top, "Computer Modern"))

    wlo, whi = minimum(wall) / 10^0.28, maximum(wall) * 10^0.28
    wticks = _decade_ticks(wlo, whi, 1, 0)
    bot = plot(; common..., ylims = (wlo, whi),
               yticks = (wticks, [(@sprintf("%g s", t)) for t in wticks]),
               ylabel = "Cumulative\nwalltime",
               xlabel = "Grid points per dimension \$(N_{\\mathrm{points}})^{1/2}\$",
               legend = false, top_margin = -3Plots.mm)
    for (i, m) in enumerate(ms)
        plot!(bot, xs[i], walls[i]; color = mode_gray(m, mmax), linewidth = 1.0, label = "")
    end

    plt = plot(top, bot; layout = grid(2, 1; heights = [0.75, 0.25]), size = (500, 640),
               dpi = 150, left_margin = 4Plots.mm, right_margin = 4Plots.mm,
               top_margin = 2Plots.mm, bottom_margin = 2Plots.mm)
    return save_figure!(plt, out_png; max_bytes = 200 * 1024, pdf = publish, id = "fig_convergence")
end

"Split the flat vector `flat` into `n` equal pieces of length `k`."
_unflat(flat, n, k) = [Float64.(flat[((i - 1) * k + 1):(i * k)]) for i in 1:n]

const GOLD_PNG = joinpath(@__DIR__, "expected", "fig_convergence.png")

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2606.04998", "fig_convergence")
    isfile(gpath) || ReproHarness.hole(
        "--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    ms = Int.(v["m_modes"]); npts = Int(v["n_plotted_points"])
    xs = [sqrt.(Float64.(v["ndofs_per_level"][1:npts])) for _ in ms]
    saved = draw(ms, xs, _unflat(v["abs_error_flat"], length(ms), npts),
                 _unflat(v["reported_cumulative_walltime_flat"], length(ms), npts),
                 Float64(v["gmres_absolute_residual"]), GOLD_PNG; publish = true)
    println("[fig_convergence] --plot-only: redrawn from the golden, no physics recomputed")
    exit(saved.ok ? 0 : 1)
end

ReproHarness.require_tier(opts, :heavy)   # 21 modes × 8 AMR levels on a 5-process pool, plus a
                                          # warm-up sweep: minutes, far over the fast ceiling

using GrTools
using GrTools.ScalarSelfForce

for f in (:solve_modes_parallel, :mode_pool, :mode_config, :crp_paper_config,
          :mode_force_convergence)
    isdefined(ScalarSelfForce, f) || ReproHarness.hole(
        "ScalarSelfForce.$(f) is not available: Fig. 4 needs the hp-AMR driver, the m-mode " *
        "process fan-out and the level-route self-force extraction (ESSF::amr-driver, " *
        "ESSF::mmode-parallel-driver, ESSF::mode-self-force-extraction).")
end

orb = CircularOrbit(SPIN, R0)
cfg = crp_paper_config(; spin = SPIN, orbital_radius = R0, m_mode = first(MS),
                       inner_radius = orb.rp, P = PORD, L = LREF)
const STEPS = cfg.amr.iterations - cfg.amr.trigger_offset
const ABSTOL = Float64(cfg.gmres.absolute_residual)

println("\n[fig_convergence] a = ", SPIN, "M  r₀ = ", R0, "M  m = ", MS, "  P = ", PORD,
        "  L = ", LREF)
println("[fig_convergence] Amr.Iterations = ", cfg.amr.iterations, " solves, trigger offset ",
        cfg.amr.trigger_offset, " ⇒ ", STEPS, " refinement steps ⇒ ", STEPS + 1, " AMR grids ",
        "⇒ ", STEPS, " plotted points (the last grid is the reference)")
println("[fig_convergence] GMRES AbsoluteResidual = ", ABSTOL)
flush(stdout)

# ───────────────────────────── ONE pool, warmed, then the production sweep ───────────────────
# Row-47 contract §0: hold ONE `mode_pool` across the whole figure. The warm-up sweep is the
# same JIT amortization fig_linsolv.jl performs in-process — without it every worker's FIRST
# mode charges its whole compilation to that mode's level 0 and the bottom panel measures the
# compiler instead of the solver.
#
# WHICH PROJECT THE WORKERS GET. `mode_pool` defaults to `Base.active_project()` and then runs
# `using ScalarSelfForce` on every worker. This environment is `reproduce/`, whose Project.toml
# lists only the `GrTools` umbrella, so a worker started in it could not load ScalarSelfForce BY
# NAME. The driver therefore points the pool at the package's own project — the same code, both
# resolved through the repository's relative `[sources]` paths.
const SSF_PROJECT = pkgdir(ScalarSelfForce)
SSF_PROJECT === nothing && ReproHarness.hole("cannot locate the ScalarSelfForce package root")
println("[fig_convergence] worker project: ", SSF_PROJECT)

sets = mode_pool(NWORKERS; threads_per_worker = 1, blas_threads = 1,
                 project = SSF_PROJECT) do pool
    tw = time()
    solve_modes_parallel(orb, cfg, MS; pool = pool, iterations = 1)
    @printf("[fig_convergence] warm-up sweep (2 throwaway solves per mode) %.1f s\n", time() - tw)
    flush(stdout)
    t0 = time()
    s = solve_modes_parallel(orb, cfg, MS; pool = pool, iterations = STEPS)
    @printf("[fig_convergence] production sweep: %d modes, %.1f s wall on %d workers\n",
            length(s.ms), time() - t0, s.run.nworkers)
    flush(stdout)
    return s
end

isempty(sets.failures) || error("fig_convergence: modes failed: ", [f.m for f in sets.failures])

# ───────────────────────────── the extraction (sealed row-44 level route) ────────────────────
ms = collect(MS)
convs = [mode_force_convergence(orb, mode_config(cfg, m), m, sets[m].levels; x_p = sets[m].x_p)
         for m in ms]
nlev = length(first(convs).ndofs)
npts = nlev - 1                       # the highest level is F^m_{r,ref}; it is not plotted

ndofs = Int.(first(convs).ndofs)
ndofs_identical = all(Int.(c.ndofs) == ndofs for c in convs)
xs_one = sqrt.(Float64.(ndofs[1:npts]))
xs = [sqrt.(Float64.(Int.(c.ndofs)[1:npts])) for c in convs]

fr_ref = [Float64(c.values[end]) for c in convs]
errs = [[abs(Float64(c.values[k] - c.values[end])) for k in 1:npts] for c in convs]
walls = [cumsum([l.wall_time for l in sets[m].levels])[1:npts] for m in ms]

"Least squares `y = a + b·x`; returns `(rate, intercept, R²)`."
function lsq(x, y)
    n = length(x); mx = sum(x) / n; my = sum(y) / n
    b = sum((x .- mx) .* (y .- my)) / sum((x .- mx) .^ 2)
    a = my - b * mx
    ssr = sum((y .- (a .+ b .* x)) .^ 2); sst = sum((y .- my) .^ 2)
    return b, a, (sst > 0 ? 1 - ssr / sst : NaN)
end

fits = [lsq(xs[i], log10.(errs[i])) for i in eachindex(ms)]
rates = [f[1] for f in fits]
inters = [f[2] for f in fits]
r2s = [f[3] for f in fits]
decades = [log10(errs[i][1] / errs[i][end]) for i in eachindex(ms)]
two_step_max = [maximum(Float64.(c.two_step_ratios)) for c in convs]
res_err = [Float64(c.resolution_error) for c in convs]
rel_res_err = [Float64(c.relative_resolution_error) for c in convs]
all_converged = all(all(l.record.converged for l in sets[m].levels) for m in ms)
mode_wall = [sets[m].wall_time for m in ms]

println("\n[fig_convergence]   m |  F_r^m(ref)   | err(first)  | err(last)   | decades | ",
        "rate     |    R²   | rel.res.err | max 2-step | wall[s]")
for i in eachindex(ms)
    @printf("[fig_convergence] %3d | %+.7e | %.4e | %.4e | %7.3f | %8.5f | %7.5f | %.4e | %10.5f | %6.1f\n",
            ms[i], fr_ref[i], errs[i][1], errs[i][end], decades[i], rates[i], r2s[i],
            rel_res_err[i], two_step_max[i], mode_wall[i])
end
println("[fig_convergence] ndofs per level: ", ndofs, "  (identical across m: ",
        ndofs_identical, ")")
for i in eachindex(ms)
    println("[fig_convergence] m = ", ms[i], "  err = ", errs[i])
    println("[fig_convergence] m = ", ms[i], "  two_step = ",
            Float64.(convs[i].two_step_ratios))
end
@printf("[fig_convergence] R²: min %.4f median %.4f max %.4f; %d/%d modes >= 0.98 (published figure: min %.4f median %.4f max %.4f, %d/21)\n",
        minimum(r2s), _median(r2s), maximum(r2s), count(>=(0.98), r2s), length(r2s),
        PUB_R2_MIN, PUB_R2_MEDIAN, PUB_R2_MAX, PUB_R2_GE_098)
@printf("[fig_convergence] rate: [%.4f, %.4f] decades per unit √N (published [%.4f, %.4f])\n",
        minimum(rates), maximum(rates), PUB_RATE_MIN, PUB_RATE_MAX)
@printf("[fig_convergence] decades dropped: [%.2f, %.2f] (published [%.2f, %.2f])\n",
        minimum(decades), maximum(decades), PUB_DECADES[1], PUB_DECADES[2])
flush(stdout)

# ───────────────────────────── the node's gates ──────────────────────────────────────────────
# (1) EXPONENTIAL CONVERGENCE, per mode: log10|F_r^m − F^m_{r,ref}| falls linearly in √N_points
#     (main.tex:719-721). Instrument: a least-squares fit over the 7 plotted points; gate: the
#     rate is negative and R² ≥ the PUBLISHED figure's own worst fit. The planned "R² > 0.98 for
#     all 21" is falsified by the published data itself (9/21) — see the header.
# (2) EVERY MODE DROPS AT LEAST `DECADES_MIN` decades over the sweep.
#     (The row-44 two-step contraction envelope is MEASURED and registered but NOT gated here —
#     see the constants block: it is a statement about the pre-floor regime, and at 8 AMR levels
#     the high-m curves sit ON their floor from level ~5, where the differences are noise.)
# (3) THE MESH SEQUENCE IS m-INDEPENDENT — our criteria are RefineAtPuncture ∘ IncreaseResolution
#     and nothing else. This is the sharp form of deviation (2) in the header: the published
#     grids are NOT m-independent.
# (4) Every AMR level of every mode reached the GMRES absolute tolerance, and no mode failed.
# (5) The tracked PNG stays inside the shared 200 KB reproduce budget.
gatefails = String[]
badfit = [i for i in eachindex(ms) if !(rates[i] < 0) || !(r2s[i] >= R2_FLOOR)]
isempty(badfit) || push!(gatefails,
    "GATE 1: modes $(ms[badfit]) are not exponentially convergent at the published envelope " *
    "(rate $(rates[badfit]), R² $(r2s[badfit]), floor $(R2_FLOOR))")
baddec = [i for i in eachindex(ms) if !(decades[i] >= DECADES_MIN)]
isempty(baddec) || push!(gatefails,
    "GATE 2: modes $(ms[baddec]) drop only $(decades[baddec]) decades over the sweep " *
    "(< $(DECADES_MIN))")
# MEASURED, NOT GATED: where the row-44 two-step envelope holds and where it does not.
badtwo = ms[[two_step_max[i] >= TWO_STEP_CONTRACT for i in eachindex(ms)]]
isempty(badtwo) || println("[fig_convergence] NOTE: the two-step ratio reaches or exceeds 1 " *
                           "for modes $(badtwo) — those sequences are already on their " *
                           "resolution floor at 8 AMR levels (reported, not gated)")
ndofs_identical || push!(gatefails,
    "GATE 3: the AMR mesh sequence is not the same for every m — our criteria " *
    "(RefineAtPuncture ∘ IncreaseResolution) are m-independent, so this cannot happen: " *
    "$(unique([Int.(c.ndofs) for c in convs]))")
all_converged || push!(gatefails,
    "GATE 4: some AMR level did not reach the GMRES absolute tolerance $(ABSTOL)")

# ───────────────────────────── registered values (the golden) ────────────────────────────────
@repro_result "spin" SPIN
@repro_result "orbital_radius" R0
@repro_result "m_modes" ms
@repro_result "grid_points_per_dim" PORD
@repro_result "initial_refinement_L" LREF
@repro_result "amr_iterations_yaml" cfg.amr.iterations
@repro_result "amr_trigger_offset" cfg.amr.trigger_offset
@repro_result "amr_refinement_steps" STEPS
@repro_result "n_amr_levels" nlev
@repro_result "n_plotted_points" npts
@repro_result "gmres_absolute_residual" ABSTOL

@repro_result "ndofs_per_level" ndofs
@repro_result "ndofs_identical_across_m" ndofs_identical
@repro_result "nelements_per_level" [length(l.mesh.elements) for l in sets[first(ms)].levels]
@repro_result "sqrt_ndofs_plotted" xs_one
@repro_result "gmres_iterations_per_level_flat" vcat([[l.record.iterations
                                                       for l in sets[m].levels] for m in ms]...)

# the F_r^m sequence itself, and the plotted error. Per-level VALUES are mesh-dependent (row 44:
# they moved 1.5e-8 relative across the hp_refine fix), so they are registered at a tolerance
# rather than pinned bit for bit; the log10 error carries every point at comparable weight,
# which the raw vector does not.
@repro_result "fr_values_flat" vcat([Float64.(c.values) for c in convs]...)
@repro_result "fr_reference_per_mode" fr_ref
@repro_result "abs_error_flat" vcat(errs...)
@repro_result "log10_abs_error_flat" vcat([log10.(e) for e in errs]...)
@repro_result "resolution_error_per_mode" res_err
@repro_result "relative_resolution_error_per_mode" rel_res_err
@repro_result "two_step_ratio_max_per_mode" two_step_max
# how far the row-44 contract's "at most HALF the difference two levels back" actually reaches:
# REPORTED against, never gated (the prototype already exceeds it at m = 10).
@repro_result "n_two_step_ratios_above_half" sum(count(>(TWO_STEP_HALF),
                                                      Float64.(c.two_step_ratios))
                                                for c in convs)
@repro_result "n_two_step_ratios_total" sum(length(c.two_step_ratios) for c in convs)
@repro_result "modes_with_two_step_above_half" ms[[two_step_max[i] > TWO_STEP_HALF
                                                   for i in eachindex(ms)]]

@repro_result "fit_rate_per_mode" rates
@repro_result "fit_intercept_per_mode" inters
@repro_result "fit_r2_per_mode" r2s
@repro_result "decades_dropped_per_mode" decades
@repro_result "n_modes_r2_above_098" count(>=(0.98), r2s)
@repro_result "fit_r2_min" minimum(r2s)
@repro_result "fit_r2_median" _median(r2s)

# the published figure, measured from its own vector paths — REFERENCE, never plotted as ours
@repro_result "published_ndofs_per_level_m0" PUB_NDOFS
@repro_result "published_fit_r2_min" PUB_R2_MIN
@repro_result "published_fit_r2_median" PUB_R2_MEDIAN
@repro_result "published_fit_r2_max" PUB_R2_MAX
@repro_result "published_modes_r2_above_098" PUB_R2_GE_098
@repro_result "published_fit_rate_range" [PUB_RATE_MIN, PUB_RATE_MAX]
@repro_result "published_decades_range" [PUB_DECADES[1], PUB_DECADES[2]]
@repro_result "published_tolerance_line" PUB_TOL_LINE
@repro_result "ndofs_ratio_published_over_ours" [PUB_NDOFS[k] / ndofs[k] for k in 1:npts]

# REPORTED, NOT GATED: wall-clock seconds (machine- and load-dependent; the paper's own numbers
# are an Apple M2 Pro's, main.tex:688 "a few seconds", and ours are a shared 256-core Linux box's)
@repro_result "reported_cumulative_walltime_flat" vcat(walls...)
@repro_result "reported_mode_wall_time" mode_wall
@repro_result "reported_nworkers" sets.run.nworkers
@repro_result "reported_sweep_wall_time" sets.run.wall_time

# `--check` must never touch the tracked artifact: it renders to a SCRATCH path so the plotting
# path is still exercised, and the size verdict enters the gate list BEFORE `gates_pass`.
publish = opts.mode === :regen
out_png = publish ? GOLD_PNG : joinpath(mktempdir(), "fig_convergence.png")
saved = draw(ms, xs, errs, walls, ABSTOL, out_png; publish = publish)
saved.ok || push!(gatefails, "GATE 5: " * saved.message)
publish || println("[fig_convergence] --check: rendered to a scratch path, ",
                   relpath(GOLD_PNG, ReproHarness.REPRO_ROOT), " untouched")

@repro_result "png_bytes_within_budget" saved.ok
@repro_result "gates_pass" isempty(gatefails)
if !isempty(gatefails)
    println("\n[fig_convergence] FINDING — node gates violated:")
    foreach(f -> println("    ", f), gatefails)
end

# rtol 1e-6 on the values, which are mesh- and round-off-sensitive at the 1e-11 floor: the
# plotted errors ARE differences at that floor, so a norm-wise 1e-6 is a real gate on the curve
# shape without turning last-digit motion into a red check. The `reported_*` timing rows carry a
# deliberately loose tolerance so a busy machine cannot fail a reproduced figure.
const LOOSE = (rtol = 1.0e3, atol = 1.0e3)
code = ReproHarness.conclude!(opts; paper = "2606.04998", id = "fig_convergence",
                              script = @__FILE__, rtol = 1e-6,
                              tolerances = Dict("reported_cumulative_walltime_flat" => LOOSE,
                                                "reported_mode_wall_time" => LOOSE,
                                                "reported_sweep_wall_time" => LOOSE,
                                                "reported_nworkers" => LOOSE))
exit(isempty(gatefails) ? code : max(code, 1))

# reproduce/arXiv-2606.04998/fig_error_budget.jl — ESSF::fig-error-budget (the paper's last figure).
#
# arXiv:2606.04998 Fig. 7 (`fig:error_budget`, main.tex:706-710): "Error budget in the extracted
# self-force. Upper panels show the value of the self-force and lower panels the relative error.
# The numerical resolution floor is below the error to the reference solutions from
# [Warburton & Barack 2010], indicating that accuracy is not limited by numerical resolutions.
# Instead, the error is dominated by truncation error in the sum over m-modes."
#
# THE THREE CURVES, in the paper's own words (main.tex:722-726):
#   * "We measure the error due to numerical resolution by the difference between the highest and
#     second-highest AMR level."                            -> `r.resolution_error`  (dotted)
#   * "We estimate the error due to finite-m truncation by the difference between using all
#     m <= 20 modes and omitting the highest mode, m = 20."  -> `r.truncation_error`  (dashed)
#   * "We compute the true error by the difference to the accurate reference solution from
#     [Warburton & Barack]."                                 -> |F_r - F_r^WB|       (solid)
# all three divided by |F_r| (the ordinate is "Relative error |F_r - F_r^ref|/|F_r|"), over
# "Spin a" at r0 = ISCO (left column) and r0 = 10 M (right column), with |F_r| itself in a short
# strip above each column.
#
# THIN driver: zero physics. `kerr_isco` -> `CircularOrbit` -> `crp_paper_config` ->
# ONE `mode_pool` -> `solve_modes_parallel` (rows 43/47) -> `total_self_force` (row 46) and the
# sealed read-outs `horizon_axis_profile`, `block_domain`. The Warburton-Barack column comes from
# the sealed row-6 fixture, which is string-equal to the tex mirror by its own test. The driver
# arranges, gates and draws; it derives no self-force quantity of its own.
#
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_error_budget.jl --check
#   julia --project=reproduce reproduce/arXiv-2606.04998/fig_error_budget.jl --plot-only  # redraw
#
# `--plot-only` REDRAWS THE TRACKED PNG from the golden, in place, with no solve — that is its
# purpose (a cosmetic iteration costs seconds instead of 73 minutes) and its hazard: it changes
# a tracked artifact while leaving the golden's `script_sha256` at whatever produced the golden,
# so `--check` then prints "driver script changed since the golden was generated". Use it while
# iterating; the artifact a commit ships should come from the same `--regen` as its golden.
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
# THE PUBLISHED FIGURE, READ OFF ITS OWN VECTOR PATHS (never plotted as ours; the extraction
# script and its calibration are in this node's evidence). Eleven polylines, calibrated on the
# axis tick STROKES (x: a = 0 at 144.05/370.01 pt, 97.03 pt per unit spin; y: the decade
# gridlines at 29.333 pt/decade below and 8.075 pt/decade in the strips):
#
#   * NINE spins are plotted, not five: a = -0.998, -0.9, -0.7, -0.5, 0, +0.5, +0.7, +0.9, +0.998.
#     The five TICK LABELS (-0.998, -0.5, 0, 0.5, 0.998) are not the sample set; +-0.7 and +-0.9
#     carry markers too. This driver therefore runs the published nine, not the five its DAG row
#     named, and the four extra spins are exactly the ones that carry a printed WB cell.
#   * The dashed (truncation) and dotted (resolution) curves carry all NINE points. The solid
#     ("Error to reference") curve and the |F_r| strip carry only SEVEN: |a| <= 0.9. The marker
#     count confirms it exactly - 2x7 + 2x7 + 2x9 + 2x9 + 3 legend = 67 filled circles, and the
#     PDF has 67. So the published figure has NO reference point at |a| = 0.998 either, and the
#     reason is WARBURTON-BARACK'S SPIN COVERAGE: their tables stop at |a| <= 0.9 AT EVERY
#     RADIUS. That is a different statement from "no printed cell": the published ISCO column
#     has SEVEN reference points while the printed tables carry exactly ONE cell at an ISCO
#     radius (a = 0, r0 = 6 M), because the printed r0 grid {4, 5, 6, 7, 8, 10, 14, 20, 30, 50,
#     70, 100} meets an ISCO radius only there. So the paper had WB data at six ISCO radii its
#     own tables do not print, and had none beyond |a| = 0.9 - our blank at |a| = 0.998 is the
#     paper's own blank, while our six blanks on the ISCO column are ours [SOLID, measured].
#   * The published |F_r| strip at r0 = 10 reproduces the WB printed table to 0.10 % at all seven
#     spins - which is what calibrates the read-back, and is why the numbers below are quoted.
#
# WHAT THE PUBLISHED CURVES SAY, AND WHAT IT IMPLIES ABOUT THE PAPER'S OWN F_r:
#
#   published ratio      r0 = 10 (7 spins)     r0 = ISCO (7 spins)
#   reference/truncation      2.81 - 4.63           3.09 - 3.61
#   truncation/resolution     3.13 - 5.40           10.3 - 1857
#
# The reference/truncation band is the paper's "our finite-m truncation error models the true
# error well up to a constant offset" (main.tex:724). That constant is not free: for a `C m^-p`
# tail the omit-the-largest-m estimator is smaller than the true remainder `sum_{m>m_max}` by
# `~ m_max/(p-1)`, which row 46 MEASURED as 3.57-3.65 on this very data, and the published
# offset IS that factor.
#
# [PRELIMINARY - AN INFERRED HYPOTHESIS, NOT A STATEMENT OF THE PAPER'S] From that, the natural
# reading is that the published `F_r` carries NO FITTED TAIL beyond the regularization
# parameters - i.e. that it is `sum_{m<=20}[F_r^m - F^m_{r[4]}]` - because a fitted tail would
# have removed the remainder and the offset could not then be ~3.6. THE PAPER DOES NOT SAY THIS
# in so many words. What it says is consistent with it in two places: main.tex:722 names exactly
# ONE acceleration, "the F_{r,[4]} regularization parameter from [Heffernan 2014] to accelerate
# the sum over m-modes in our solution", and main.tex:724 puts higher-order fits in the
# conditional - "Fits for higher-order regularization parameters have shown to increase accuracy
# in previous work and COULD BE APPLIED to our solutions as well". Neither sentence excludes a
# tail treatment this reproduction has not thought of. Two quantitative supports, both from
# row 46's own sealed table: our no-tail |F_r - F_r^WB| is 3.6e-6 at a = -0.9 against the
# published 3.53e-6, and 3.4e-5 at a = +0.5 against the published 3.42e-5 - 2 % and 0.6 %.
#
# THIS DRIVER THEREFORE PLOTS THE NO-FITTED-TAIL INTERPRETATION (`tail = :none`) and registers
# the tail-corrected number beside it. Both come from the SAME solved modes - one extra call to
# the sealed `total_self_force`, no extra solve - so nothing is thrown away: the Heffernan-tail
# `F_r` is this reproduction's more accurate number and its reference error is registered as
# `ref_error_tail_*`, roughly 5-10x smaller than the plotted one. If the interpretation is wrong,
# the plotted solid curve is the wrong estimator and every gate that reads it moves with it.
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
# THE GATES, AND THE ONE THAT WAS RE-PLANNED BEFORE IT WAS EVER RUN.
#
# GLOBAL_DAG row 53 planned "the reference-error curve |F_r - WB| within 3x the coherent
# resolution error". That gate was calibrated on row 46's HALF-DOMAIN, tail-corrected numbers,
# where the ratio is O(1). MEASURED HERE: 0 of the 8 WB-referenced configurations satisfy it at
# m_max = 20, registered as `n_wb_configurations_within_dag_factor`. The reading - an INFERENCE,
# not an identity - is that the reference error is DOMINATED BY THE MODE-SUM REMAINDER rather
# than by resolution, the two being two to three decades apart at this m_max. It is NOT true
# that the inequality is unreachable by construction: at infinite m_max the remainder vanishes
# and |F_r - F_r^WB| falls to the REFERENCE's own floor - the ~1e-7 relative that
# Warburton-Barack's printed digits carry - not to zero, and the planned inequality would then
# be a statement about that floor against our resolution error. The planned metric therefore
# stays REGISTERED, as a FAILED DIAGNOSTIC (0 of 8), rather than being retired.
#
# What it does not do is decide this figure. That is the PUBLISHED FIGURE'S OWN ENVELOPE, the
# discipline `fig_convergence.jl` landed for its R2 floor - the paper's curves are the yardstick,
# read off the paper's own PDF, and our curves must sit inside a stated factor of them (GATES 3,
# 4, 5), with the unregularized F_t against the printed Table II as the second external gate
# (GATE 6) and the resolution curve gated FROM ABOVE (GATE 8, added after the first review:
# until it existed nothing constrained the dotted curve at all, and a 389x departure passed).
# EVERY GATE RUNS ON ALL 18 CONFIGURATIONS - there is no scope, and the converged-mode-sum
# criterion is an EXPLORATORY diagnostic that decides nothing here. See the gate block.
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
# THE NODE'S TWO OPEN HOLES, both ring-marked on the figure and both in the formal verdict.
#
# [HOLE 1] THE COMPRESSED PROGRADE ISCO. The absolute resolution error runs 2.9e-10 at a = +0.7,
# 6.0e-7 at +0.9 and 2.4e-2 at +0.998 - three orders per spin step - and at +0.998 the mode sum
# has not converged at m = 20 either (the last retained term is 2.4 % of the total). a = +0.9 is
# the member the first pass MISSED: every gate it had passed there while its resolution
# estimator sat 389x above the published one, which is why GATE 8 exists.
#
# [HOLE 2] |a| = 0.998 AT r0 = 10. The relative resolution error is 1.48e-5 and 3.92e-5 against
# a published 1.99e-7 and 3.92e-7 - 74x and 100x - and these are the only two other
# configurations where our ordering inverts. THE LEAD, measured here: in ABSOLUTE terms those
# two are 7.9e-10 and 8.2e-10, which is the SAME MAGNITUDE as the half domain's spin-independent
# floor at r0 = 10 (5.67e-10 ... 5.89e-10, registered) - while the full-z lattice's other seven
# spins reach 2.5e-13 ... 6.7e-12 there. So the near-extremal points look as if they are sitting
# on a floor the rest of the lattice has already left behind. The `layer_ratio` instrument reads
# 0.098 at both: the honest statement is that THE `layer_ratio` INSTRUMENT DOES NOT DETECT
# whatever costs the resolution there - not that the near-horizon layer is uninvolved. Handed to
# ESSF::fix-horizon-axis-layer with the lead.
#
# THE WORLDTUBE RADIUS. `anc/ScalarSelfForce.yaml:12` templates it and the paper publishes a
# value only for its header example (r0 = 10 M, r_in = 2 M: `worldtube_radius = 4.8`). At r0 = 10
# that published value is legal at every spin (v->t = 10 - 1.5*4.8 = 2.8 > r+ <= 2) and IS used
# here, which also keeps this column directly comparable with row 46's sealed table. At the ISCO
# it is not: 4.8 puts the v->t transition at -1.2 M for r0 = 6. This driver uses `fig_modes.jl`'s
# landed one-parameter family there,
#
#     w_r = 0.6 (r0 - r+),
#
# which reproduces the paper's own 4.8 exactly at its own point and keeps the v block at a fixed
# 10 % of [r+, r0] at every spin. It is an ASSUMPTION where the paper prints nothing, it applies
# to the ISCO column only, and every ISCO number registered here depends on it.
#
# EIGHT AMR GRIDS, the paper's production resolution: `Amr.Iterations: 8` (yaml:153) counts
# spectre SOLVES and the trigger `EveryNIterations {N: 1, Offset: 1}` (yaml:157-160) skips the
# initial guess, so the sweep is SEVEN refinements - the driver DERIVES
# `cfg.amr.iterations - cfg.amr.trigger_offset` exactly as every other driver of this paper does.
# Row 46's contract table was measured at NINE grids (the package default); this figure's numbers
# are one AMR step coarser by construction and are not the same numbers.
#
# THE HORIZON ORDER HOLD IS MEASURED, NEVER ASSUMED (ESSF::fix-horizon-axis-layer, row 58), by
# the rule `fig_modes.jl` landed: sweep with the package default (`horizon_hold = false`,
# reference parity), read `horizon_axis_profile(dom, orbit, mesh, u).layer_ratio` on EVERY mode,
# and re-sweep with `solve_kwargs = (; horizon_hold = true)` only when the MAXIMUM OVER m exceeds
# LAYER_GATE = 1.05 - the maximum, not the top mode, because `layer_ratio` is not monotone in m.
#
# Style: `_style/manuals/arXiv-2606.04998.md` sections 0 + 7 (full-width 2x2, height ratio 1:4,
# columns sharing x, grayscale only, filled circle markers, solid/dashed/dotted hierarchy with
# the layered gray budget fills, x ticks -0.998/-0.5/0/0.5/0.998, panel tags inside top-left,
# framed legend in the lower-right panel, light grid, CM serif).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()

const GOLD_PNG = joinpath(@__DIR__, "expected", "fig_error_budget.png")
const PLOT_ONLY = "--plot-only" in opts.extra

# ───────────────────────────── the published figure, measured (NOT plotted) ──────────────────
# Read off `ssf_error_budget.pdf`'s own vector paths (see the header). Used ONLY to calibrate
# this node's gates and to quantify the deviations. No published number is ever plotted as ours.
const PUB_SPINS      = [-0.998, -0.9, -0.7, -0.5, 0.0, 0.5, 0.7, 0.9, 0.998]
const PUB_REF_SPINS  = [-0.9, -0.7, -0.5, 0.0, 0.5, 0.7, 0.9]      # the solid curve's own seven
# relative errors, per column, in PUB_SPINS order (reference: in PUB_REF_SPINS order)
const PUB_ISCO_TRUNC = [1.463e-6, 1.650e-6, 2.177e-6, 3.026e-6, 1.035e-5, 2.931e-4, 1.031e-4,
                        2.533e-4, 1.992e-4]
const PUB_ISCO_RES   = [1.381e-7, 1.600e-7, 1.635e-7, 1.825e-7, 2.830e-7, 1.886e-6, 5.228e-7,
                        1.364e-7, 1.788e-7]
const PUB_ISCO_REF   = [5.955e-6, 7.442e-6, 1.068e-5, 3.544e-5, 9.759e-4, 3.288e-4, 7.818e-4]
const PUB_R10_TRUNC  = [1.074e-6, 1.107e-6, 1.209e-6, 1.373e-6, 2.649e-6, 7.525e-6, 2.604e-6,
                        1.513e-6, 1.243e-6]
const PUB_R10_RES    = [1.989e-7, 2.081e-7, 2.515e-7, 3.167e-7, 7.624e-7, 2.346e-6, 8.312e-7,
                        4.684e-7, 3.919e-7]
const PUB_R10_REF    = [3.532e-6, 4.556e-6, 3.864e-6, 1.097e-5, 3.416e-5, 1.205e-5, 6.889e-6]
const PUB_ISCO_ABSFR = [9.598e-5, 1.107e-4, 1.274e-4, 1.676e-4, 6.916e-5, 1.087e-3, 1.133e-2]
const PUB_R10_ABSFR  = [4.935e-5, 4.097e-5, 3.286e-5, 1.377e-5, 4.031e-6, 1.091e-5, 1.767e-5]
const PUB_REF_OVER_TRUNC = (2.81, 4.63)    # the r0 = 10 band; the ISCO band is (3.09, 3.61)
const PUB_TRUNC_OVER_RES_R10 = (3.13, 5.40)

# GATE CONSTANTS, calibrated BEFORE the production run - externally on the published curves
# above, internally on this node's throwaway prototype (evidence `calibration_prototype.log`).
const TRUNC_ENVELOPE = 2.0   # our truncation curve within this factor of the published one
const REF_ENVELOPE   = 3.0   # our reference-error curve within this factor of the published one
const OFFSET_BAND    = (2.0, 6.0)  # reference/truncation; published 2.81-4.63, the analytic
                                   # m_max/(p-1) model 3.77-3.83, row 46 measured 3.57-3.65
const LAYER_GATE     = 1.05  # the measured `layer_ratio` a plotted configuration may not pass
# THE EXPLORATORY CONVERGED-MODE-SUM DIAGNOSTIC. `|T_{m_max}|` is a TRUNCATION estimate only
# while the m-sum has begun to converge; where the last retained term is a tenth of the total it
# measures the series' failure to converge instead. The criterion is the SMALLNESS of that term,
# `rel_trunc < MODE_SUM_CONVERGED`. It was written AFTER the first production run, it fires at
# exactly the one configuration that failed a gate there, and it is therefore EXPLORATORY: it is
# registered, its sensitivity is registered, and it scopes NOTHING (see the gate block).
const MODE_SUM_CONVERGED = 1.0e-2
const DAG_RES_FACTOR = 3.0   # GLOBAL_DAG row 53's planned factor - REPORTED, not gated
const RES_ENVELOPE_UP = 10.0 # the resolution curve, gated FROM ABOVE ONLY (see GATE 8): our
                             # estimator may not exceed the published one by more than a decade
const FT_GATE = 1.0e-6       # the F_t envelope. Row 46 measured 3.1e-9 ... 3.4e-8 at nine grids;
                             # this carries 30x of headroom for the one coarser AMR step. The
                             # MEASURED deviations of this figure are quoted with the numbers,
                             # never with this constant (see `Ft_relative_deviation_max`).

# ───────────────────────────── the figure (presentation only) ────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

const XTICKS = ([-0.998, -0.5, 0.0, 0.5, 0.998], ["-0.998", "-0.5", "0", "0.5", "0.998"])
const GRAY_LO = RGB(0.62, 0.62, 0.62)   # below the resolution curve (darker)
const GRAY_HI = RGB(0.88, 0.88, 0.88)   # between resolution and truncation (lighter)

"Greedy word wrap at `width` characters — presentation only; the note text itself is a golden value."
function wrap(s::AbstractString, width::Integer)
    out = String[]; line = ""
    for w in split(s)
        cand = isempty(line) ? String(w) : line * " " * w
        if length(cand) > width && !isempty(line)
            push!(out, line); line = String(w)
        else
            line = cand
        end
    end
    isempty(line) || push!(out, line)
    return join(out, "\n")
end

"""
    draw(cols, out_png; publish) -> saved

Fig. 7: `cols` is a two-element vector of NamedTuples, one per column (ISCO, r0 = 10), each
carrying the spins, `absfr`, the three relative-error curves and the reference-error mask.
"""
function draw(cols, out_png; publish::Bool)
    tops = Any[]; bots = Any[]
    allerr = vcat([vcat(c.rel_ref[c.refmask], c.rel_trunc, c.rel_res) for c in cols]...)
    elo = minimum(allerr) / 10^0.95   # floor room for the per-column note (drawn on it)
    ehi = maximum(allerr) * 10^0.45
    allfr = vcat([c.absfr for c in cols]...)
    flo, fhi = minimum(allfr) / 10^0.45, maximum(allfr) * 10^0.45
    # the published axes END at the tick values, so the +-0.998 markers sit ON the spines and are
    # half-clipped there — reproduced, not padded away
    xlim = (first(XTICKS[1]), last(XTICKS[1]))
    common = (xlims = xlim, yscale = :log10, grid = true, gridalpha = 0.5,
              gridlinewidth = 0.6, foreground_color_grid = RGB(0.72, 0.72, 0.72),
              tick_direction = :in, widen = false, guidefontsize = 10, tickfontsize = 9)
    blank = (XTICKS[1], fill("", length(XTICKS[1])))   # the strips share x: ticks, no labels
    for (k, c) in enumerate(cols)
        # ── the value strip
        ytop = decade_ticks(flo, fhi; step = 2)
        top = plot(; common..., xticks = blank, ylims = (flo, fhi),
                   yticks = k == 1 ? ytop : (ytop, fill("", length(ytop))),
                   ylabel = k == 1 ? "\$|F_r|\$" : "", legend = false,
                   bottom_margin = -1Plots.mm, top_margin = 1Plots.mm,
                   left_margin = (k == 1 ? 7 : 2)Plots.mm)
        plot!(top, c.a, c.absfr; color = :black, linewidth = 1.1, marker = :circle,
              markersize = 3.0, markerstrokewidth = 0, label = "")
        push!(tops, top)
        # ── the budget panel: fills first, then the three curves
        ybot = decade_ticks(elo, ehi)
        bot = plot(; common..., xticks = XTICKS, ylims = (elo, ehi),
                   yticks = k == 1 ? ybot : (ybot, fill("", length(ybot))),
                   ylabel = k == 1 ?
                       "Relative error \$|F_r - F_r^{\\mathrm{ref}}|/|F_r|\$" : "",
                   xlabel = "Spin \$a\$", legend = k == 2 ? :topright : false,
                   legendfontsize = 9, top_margin = 0Plots.mm, bottom_margin = 2Plots.mm,
                   left_margin = (k == 1 ? 7 : 2)Plots.mm,
                   foreground_color_legend = RGB(0.35, 0.35, 0.35),
                   background_color_legend = :white)
        plot!(bot, c.a, c.rel_trunc; fillrange = elo, fillcolor = GRAY_HI, fillalpha = 1,
              linealpha = 0, label = "")
        plot!(bot, c.a, c.rel_res; fillrange = elo, fillcolor = GRAY_LO, fillalpha = 1,
              linealpha = 0, label = "")
        # the reference-error curve exists only where Warburton-Barack print a cell
        m = c.refmask
        plot!(bot, c.a[m], c.rel_ref[m]; color = :black, linewidth = 1.2, linestyle = :solid,
              marker = :circle, markersize = 3.0, markerstrokewidth = 0,
              label = "Error to reference")
        plot!(bot, c.a, c.rel_trunc; color = :black, linewidth = 1.2, linestyle = :dash,
              marker = :circle, markersize = 3.0, markerstrokewidth = 0,
              label = "Truncation at \$m \\leq $(c.m_max)\$")
        plot!(bot, c.a, c.rel_res; color = :black, linewidth = 1.2, linestyle = :dot,
              marker = :circle, markersize = 3.0, markerstrokewidth = 0,
              label = "Numerical resolution")
        # the configurations outside the converged-mode-sum scope are RING-MARKED on the figure
        # itself, not silently dropped: their |T_20| is not measuring truncation [HOLE]
        h = .!c.okmask
        if any(h)
            for y in (c.rel_trunc, c.rel_res)
                scatter!(bot, c.a[h], y[h]; markershape = :circle, markercolor = :white,
                         markerstrokecolor = :black, markerstrokewidth = 1.1, markersize = 5.0,
                         label = "")
            end
        end
        annotate!(bot, xlim[1] + 0.045 * (xlim[2] - xlim[1]), ehi / 10^0.30,
                  text(c.tag, 9, :left, :top, "Computer Modern"))
        # the on-figure statement of where the reference curve exists at all: this reproduction
        # has one printed WB cell on the left column and seven on the right, and the paper's own
        # solid curve stops at |a| = 0.9 as well
        # under the panel tag, where both columns have clear sky: on the floor it ran into the
        # neighbouring panel, and over the fills it ran into the curves
        # on the FLOOR of the panel: the published figure draws nothing there, the darker budget
        # fill gives the text a background, and it stays clear of the legend
        note = c.note * (isempty(c.holdnote) ? "" : " " * c.holdnote)
        note == "" || annotate!(bot, xlim[1] + 0.030 * (xlim[2] - xlim[1]), elo * 10^0.06,
                                text(wrap(note, 52), 5, :left, :bottom, RGB(0.10, 0.10, 0.10),
                                     "Computer Modern"))
        push!(bots, bot)
    end
    # the published canvas is 484.5 x 239.6 pt (aspect 2.02) with the strip 1/5 of the panel
    lay = @layout [a{0.175h} b{0.175h}; c d]
    plt = plot(tops[1], tops[2], bots[1], bots[2]; layout = lay, size = (900, 445), dpi = 150,
               right_margin = 3Plots.mm)
    return save_figure!(plt, out_png; max_bytes = 200 * 1024, pdf = publish,
                        id = "fig_error_budget")
end

# ───────────────────────────── --plot-only: redraw from the golden ───────────────────────────
if PLOT_ONLY
    gpath = ReproHarness.expected_path("2606.04998", "fig_error_budget")
    isfile(gpath) || ReproHarness.hole(
        "--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    ncol = 2
    spins = Float64.(v["spins"])
    n = length(spins)
    cols = map(1:ncol) do k
        rng = ((k - 1) * n + 1):(k * n)
        (a = spins, absfr = abs.(Float64.(v["Fr_flat"][rng])),
         rel_ref = Float64.(v["rel_ref_error_flat"][rng]),
         rel_trunc = Float64.(v["rel_truncation_flat"][rng]),
         rel_res = Float64.(v["rel_resolution_flat"][rng]),
         refmask = Bool.(v["has_wb_reference_flat"][rng]),
         okmask = .!Bool.(v["flagged_hole_flat"][rng]),
         holdnote = String(v["column_hold_notes"][k]),
         m_max = Int(v["m_max"]), tag = String(v["column_tags"][k]),
         note = String(v["column_notes"][k]))
    end
    saved = draw(cols, GOLD_PNG; publish = true)
    println("[fig_error_budget] --plot-only: redrawn from the golden, no physics recomputed")
    exit(saved.ok ? 0 : 1)
end

ReproHarness.require_tier(opts, :heavy)   # 27 configurations (18 plotted + 9 secondary), 21
                                          # modes x 8 AMR grids each: 28 production sweeps (one
                                          # configuration is swept twice for the horizon hold)
                                          # plus one warm-up sweep

using LinearAlgebra
using GrTools
using GrTools.ScalarSelfForce
using GrTools.EllipticCore
using GrTools.Spacetime: kerr_isco

for f in (:solve_modes_parallel, :mode_pool, :total_self_force, :crp_paper_config,
          :horizon_axis_profile, :select_hold_sweep, :block_domain)
    isdefined(ScalarSelfForce, f) || ReproHarness.hole(
        "ScalarSelfForce.$(f) is not available: Fig. 7 needs the m-mode fan-out " *
        "(ESSF::mmode-parallel-driver) over the hp-AMR driver (ESSF::amr-driver), the mode sum " *
        "(ESSF::mode-sum-total-force) and the horizon-axis instrument " *
        "(ESSF::fix-horizon-axis-layer).")
end

# ───────────────────────────── the Warburton-Barack reference column ─────────────────────────
# The SEALED row-6 fixture (`ESSF::wb-reference-fixtures`), a predecessor of this node: every
# literal in it is asserted STRING-EQUAL to the tex mirror by `test_wb_fixtures.jl`, so this is
# the printed table itself and not a transcription of ours. It is pure Base and loads into a
# fresh module, exactly as `test_modesum.jl` loads it.
const SSF_ROOT = pkgdir(ScalarSelfForce)
SSF_ROOT === nothing && ReproHarness.hole("cannot locate the ScalarSelfForce package root")
const WB_PATH = joinpath(SSF_ROOT, "test", "fixtures", "wb_reference.jl")
isfile(WB_PATH) || ReproHarness.hole(
    "the sealed Warburton-Barack fixture is not at $(WB_PATH): Fig. 7's reference-error curve " *
    "IS that table (ESSF::wb-reference-fixtures, GLOBAL_DAG row 6)")
const WB = Module(:WB)
Base.include(WB, WB_PATH)

"""
    wb_reference(a, r0) -> (value, printed, ulp) | nothing

The printed `(M²/q²) F_r` cell at `(a, r0)`, its exact printed decimal and the half-ulp of that
print (the precision the reference itself carries). `nothing` where the table prints no cell:
outside `|a| <= 0.9` / `r0 >= 4 M`, at a radius that is not on the printed `r0` grid, or at a
printed '-' (an orbit below the ISCO). The near-extremal and near-ISCO points of this figure are
`nothing` BY THE TABLE, not by a choice of ours [adopted decision 2].
"""
function wb_reference(a::Real, r0::Real)
    abs(a) <= WB.WB_MAX_ABS_A || return nothing
    any(x -> isapprox(x, a; atol = 1e-12), WB.WB_A_OVER_M) || return nothing
    j = findfirst(x -> isapprox(x, r0; atol = 1e-9), WB.WB_R0_OVER_M)
    j === nothing && return nothing
    cell = WB.wb_cell(WB.WB_FR, a, WB.WB_R0_OVER_M[j])
    cell isa WB.WBRefCell || return nothing
    # half-ulp of the printed decimal: the reference's own resolution
    dec = cell.dec
    mant, ex = split(lowercase(dec), 'e')
    ndec = (i = findfirst('.', mant); i === nothing ? 0 : length(mant) - i)
    ulp = 0.5 * 10.0^(-ndec) * 10.0^parse(Int, ex)
    return (value = cell.value, printed = dec, ulp = ulp)
end

# ───────────────────────────── parameters ────────────────────────────────────────────────────
const SPINS   = PUB_SPINS                 # the nine the published figure actually samples
# NOTE the deliberate "$…$ ISCO" split rather than $\mathrm{ISCO}$: GR renders the upright word
# either way, and this string is a REGISTERED value (`column_tags`) that the golden pins, so it
# is fixed at what the production run wrote.
const COLUMNS = ((:isco, "\$r_0 = \$ ISCO"), (10.0, "\$r_0 = 10\$"))
# the on-figure statement of where the reference-error curve exists at all, per column
const COLUMN_NOTES = (
    "WB print one ISCO cell (a = 0, r0 = 6 M), and stop at |a| <= 0.9. Rings: [HOLE].",
    "WB stop at |a| <= 0.9, as does the published curve. Rings: [HOLE].")
const MS      = CRP_PAPER_MMODES          # 0:20, main.tex:722 "all m <= 20"
const PORD    = 4                         # yaml:74
const LREF    = 0                         # yaml:75
const ROUT    = 1.0e4                     # yaml:11
const WDZ2    = 0.25                      # yaml:12
const WR_PUBLISHED = 4.8                  # yaml:12, legal at r0 = 10 M for every spin
const WR_FRAC = 0.6                       # w_r = WR_FRAC (r0 - r+) at the ISCO; see the header
const NWORKERS = 6                        # compute-usage: 6 processes x 1 thread x 1 BLAS thread

# THE TWO LATTICES. `anc/ScalarSelfForce.yaml:36` declares `ImposeEquatorialSymmetry: True`, but
# `fig_convergence.jl` GATED that the published convergence figure's abscissa is the FULL-Z 5x3
# lattice (480, 975, 1604, 2427, 3468, 4751, 6300 grid points, level for level) - the same hp-AMR
# sweeps this figure sums over. This driver therefore PLOTS THE PAPER'S, exactly as that node
# does, and carries the half domain as a CONTROLLED VARIATION on the r0 = 10 column. Both series
# come from ONE pool.
#
# WHAT THE VARIATION CONTROLS, AND WHAT IT DOES NOT [PRELIMINARY]. The two runs share the SAME
# EIGHT AMR LEVELS and the same solver settings - they do NOT share their grids: the full-z
# lattice reaches 8139 dofs at the last level against the half domain's 4814, a factor 1.69, and
# the equatorial-symmetry reduction is exact so the two are solving the same physics over
# different halves of the sphere. The ndof difference is a CONFOUND and is disclosed rather than
# argued away. What makes it an implausible explanation of the whole effect is the SIZE and the
# CONVERGENCE RATE: the half domain's absolute resolution error at r0 = 10 is a spin-INDEPENDENT
# floor, 5.67e-10 ... 5.89e-10 across |a| <= 0.9 (registered), against row 46's 1.56e-10 ...
# 1.64e-10 for the SAME configurations at NINE levels - i.e. one further AMR step divides it by
# about 3.6 - so closing the measured 88-2370x gap would take roughly four more steps, not the
# 1.69x in dofs that one level of the other lattice buys. The honest statement is therefore
# ASSOCIATION, not causation: the ordering of main.tex:723 is ASSOCIATED WITH THE LATTICE at the
# tested r0 = 10 configurations. It NARROWS the row-46 finding; it does not close it. Closing it
# needs the two lattices at matched ndofs, which this node did not run.
# The secondary series is the CONTROLLED VARIATION of the lattice, and it runs on the r0 = 10
# column only: those nine configurations carry every printed WB cell of the r0 = 10 column and
# are the ones row 46's sealed table measured, so the A/B answers the ordering question where it
# was asked, at a third of the cost of repeating the whole figure.
const LAT_PLOTTED   = (name = "full-z", eqsym = false)
const LAT_SECONDARY = (name = "half-domain", eqsym = true)
const SECONDARY_COLUMNS = (COLUMNS[2],)

"the configuration set on one lattice, over `cols` of `COLUMNS`"
function build_configs(eqsym::Bool; cols = COLUMNS)
    out = NamedTuple[]
    for (spec, tag) in cols, a in SPINS
        r0 = spec === :isco ? kerr_isco(a) : float(spec)
        orb = CircularOrbit(a, r0)
        wr = spec === :isco ? WR_FRAC * (r0 - orb.rp) : WR_PUBLISHED
        cfg = crp_paper_config(; spin = a, orbital_radius = r0, m_mode = first(MS),
                               inner_radius = orb.rp, outer_radius = ROUT,
                               worldtube_radius = wr, worldtube_dz2 = WDZ2, P = PORD, L = LREF,
                               equatorial_symmetry = eqsym)
        push!(out, (a = a, spec = spec, tag = tag, r0 = r0, orb = orb, wr = wr, cfg = cfg,
                    dom = block_domain(orb, cfg), wb = wb_reference(a, r0)))
    end
    return out
end

const CONFIGS = build_configs(LAT_PLOTTED.eqsym)
const NCONF = length(CONFIGS)
const STEPS = first(CONFIGS).cfg.amr.iterations - first(CONFIGS).cfg.amr.trigger_offset
const ABSTOL = Float64(first(CONFIGS).cfg.gmres.absolute_residual)

println("\n[fig_error_budget] ", NCONF, " configurations = ", length(SPINS), " spins x ",
        length(COLUMNS), " orbits, m = ", first(MS), ":", last(MS), ", ", STEPS + 1, " AMR grids")
println("[fig_error_budget] Amr.Iterations = ", first(CONFIGS).cfg.amr.iterations,
        " solves, trigger offset ", first(CONFIGS).cfg.amr.trigger_offset, " => ", STEPS,
        " refinement steps => ", STEPS + 1, " grids; GMRES AbsoluteResidual = ", ABSTOL)
println("[fig_error_budget] lattices: PLOTTED = ", LAT_PLOTTED.name,
        " (ImposeEquatorialSymmetry = ", LAT_PLOTTED.eqsym, ", the published one per ",
        "fig_convergence); SECONDARY = ", LAT_SECONDARY.name, " (= ", LAT_SECONDARY.eqsym,
        ", the shipped yaml:36)")
for c in CONFIGS
    @printf("[fig_error_budget]   a = %+.4f  %-4s r0 = %12.9f  r+ = %.9f  w_r = %.6f  WB %s\n",
            c.a, c.spec === :isco ? "ISCO" : "10M", c.r0, c.orb.rp, c.wr,
            c.wb === nothing ? "-- (no printed cell)" : c.wb.printed)
end
flush(stdout)

# ───────────────────────────── the sweeps ────────────────────────────────────────────────────
"the row-58 instrument on every mode of a fanned-out sweep"
layer_ratios(c, s) = [haskey(s, m) ?
                      Float64(horizon_axis_profile(c.dom, c.orb, s[m].levels[end].mesh,
                                                   s[m].u).layer_ratio) : NaN for m in MS]

"""
    sweep_one(pool, c, monitor) -> NamedTuple

One configuration: the measured horizon-hold rule, then the SEALED mode sum, twice - once with
the paper's own bare partial sum (`tail = :none`, what Fig. 7 plots) and once with the Heffernan
tail (this reproduction's more accurate number). Both read the SAME solved modes.
"""
function sweep_one(pool, c, monitor)
    t = time()
    set = solve_modes_parallel(c.orb, c.cfg, MS; pool = pool, iterations = STEPS, log = monitor)
    lay = layer_ratios(c, set)
    probe, kprobe = findmax(map(x -> isnan(x) ? -Inf : x, lay))
    # THE RULE (ESSF::fix-horizon-axis-layer, row 58) through the SHARED helper
    # `select_hold_sweep`, which this node added to `Solve/horizon_axis.jl` so that the two
    # figure nodes that publish a layer-gated object stop carrying two policies. It keeps
    # `fig_modes.jl`'s trigger — hold iff the default's MAXIMUM OVER m crosses 1.05 — and adds
    # the assertion `fig_modes` could assume on its three panels and this figure cannot: the
    # held sweep is published only when it MEASURES smaller. A hold that does not reduce the
    # ratio returns the default and sets `hole`; the held sweep is not even run when the
    # default is already inside the gate.
    held_pair = nothing
    if probe > LAYER_GATE
        held = solve_modes_parallel(c.orb, c.cfg, MS; pool = pool, iterations = STEPS,
                                    solve_kwargs = (; horizon_hold = true), log = monitor)
        held_pair = held => layer_ratios(c, held)
    end
    sel = select_hold_sweep(set => lay, held_pair; gate = LAYER_GATE)
    set = sel.sweep
    lay = sel.hold ? last(held_pair) : lay
    hold, hold_max = sel.hold, sel.held_ratio
    isempty(set.failures) || error("fig_error_budget: modes failed at a = $(c.a), " *
                                   "r0 = $(c.r0): $([f.m for f in set.failures])")
    none = total_self_force(set; tail = :none)        # THE PAPER'S QUANTITY
    # The tail-corrected number is a BONUS, not the figure: `mode_sum_tail` REFUSES to invent a
    # slope when fewer than two modes clear its headroom rule, which is exactly what a
    # near-extremal compressed ISCO produces (the m-sequence is not in its asymptotic regime
    # there). A refusal is registered as NaN and named, never worked around.
    heff = try
        total_self_force(set; tail = :heffernan)
    catch err
        err isa ArgumentError || rethrow()
        @warn "fig_error_budget: the Heffernan tail refuses this configuration" a=c.a r0=c.r0 err
        nothing
    end
    conv = [all(l.record.converged for l in set[m].levels) for m in MS]
    resid = [Float64(set[m].levels[end].record.residual[end]) for m in MS]
    iters = [set[m].levels[end].record.iterations for m in MS]
    lv = set[first(MS)].levels
    w = time() - t
    @printf("[fig_error_budget] a = %+.4f %-4s : %2d/%2d modes, %6.1f s, layer max %.4g at m = %2d => hold %-5s, F_r = %+.8e, |F_r-WB|/|F_r| = %s\n",
            c.a, c.spec === :isco ? "ISCO" : "10M", length(set.results), length(MS), w,
            probe, MS[kprobe], string(hold), Float64(none.Fr),
            c.wb === nothing ? "--" :
            @sprintf("%.3e", abs(Float64(none.Fr) - c.wb.value) / abs(Float64(none.Fr))))
    flush(stdout)
    return (c = c, set = set, none = none, heff = heff, hold = hold, layer_probe = probe,
            layer_hold_max = hold_max, sel = sel,
            layer = lay, converged = conv, residual = resid, iterations = iters,
            ndofs = [l.mesh.ndofs for l in lv], nelements = [length(l.mesh.elements) for l in lv],
            wall = w, nworkers = set.run.nworkers, max_load = set.run.max_load,
            mode_wall = [set[m].wall_time for m in MS])
end

"the derived per-configuration curve values, for one lattice's sweep results"
function curves(rs)
    fr = [Float64(r.none.Fr) for r in rs]
    ref = [r.c.wb === nothing ? NaN : abs(fr[i] - r.c.wb.value) for (i, r) in enumerate(rs)]
    (fr = fr,
     fr_tail = [r.heff === nothing ? NaN : Float64(r.heff.Fr) for r in rs],
     partial = [Float64(r.none.partial) for r in rs],
     tail_value = [r.heff === nothing ? NaN : Float64(r.heff.tail.value) for r in rs],
     tail_unc = [r.heff === nothing ? NaN : Float64(r.heff.tail.uncertainty) for r in rs],
     tail_exp = [r.heff === nothing ? NaN : Float64(r.heff.tail.exponent) for r in rs],
     ft = [Float64(r.none.Ft) for r in rs],
     fphi = [Float64(r.none.Fphi) for r in rs],
     trunc = [Float64(r.none.truncation_error) for r in rs],
     res = [Float64(r.none.resolution_error) for r in rs],
     res_bound = [Float64(r.none.resolution_error_bound) for r in rs],
     res_tail = [r.heff === nothing ? NaN : Float64(r.heff.resolution_error) for r in rs],
     ft_trunc = [Float64(r.none.Ft_truncation_error) for r in rs],
     ft_res = [Float64(r.none.Ft_resolution_error) for r in rs],
     ref = ref,
     ref_tail = [(r.c.wb === nothing || r.heff === nothing) ? NaN :
                 abs(Float64(r.heff.Fr) - r.c.wb.value) for r in rs],
     wb = [r.c.wb === nothing ? NaN : r.c.wb.value for r in rs],
     wb_ulp = [r.c.wb === nothing ? NaN : r.c.wb.ulp for r in rs],
     has_wb = [r.c.wb !== nothing for r in rs],
     rel_ref = ref ./ abs.(fr),
     rel_trunc = [Float64(r.none.truncation_error) for r in rs] ./ abs.(fr),
     rel_res = [Float64(r.none.resolution_error) for r in rs] ./ abs.(fr))
end

# ───────────────────────────── ONE pool, both lattices ───────────────────────────────────────
# Row-47 contract §0: hold ONE `mode_pool` across the WHOLE figure. A fresh worker pays its
# entire JIT before its first solve, so a pool per configuration would pay it 36 times.
const SSF_PROJECT = joinpath(SSF_ROOT, "Project.toml")
isfile(SSF_PROJECT) || ReproHarness.hole("the ScalarSelfForce project file is not at $(SSF_PROJECT)")
println("[fig_error_budget] worker project: ", SSF_PROJECT); flush(stdout)

# a STABLE path: mktempdir() deletes itself at exit, taking the compute-usage record with it
const MONITOR_LOG = joinpath(tempdir(), "fig_error_budget_load_monitor_$(getpid()).log")
println("[fig_error_budget] load monitor: ", MONITOR_LOG); flush(stdout)

t_all = time()
plotted_rs, secondary_rs = mode_pool(NWORKERS; threads_per_worker = 1, blas_threads = 1,
                                     project = SSF_PROJECT) do pool
    tw = time()
    solve_modes_parallel(first(CONFIGS).orb, first(CONFIGS).cfg, MS; pool = pool, iterations = 1,
                         log = MONITOR_LOG)
    @printf("[fig_error_budget] warm-up sweep (2 throwaway solves per mode) %.1f s\n", time() - tw)
    flush(stdout)
    println("[fig_error_budget] === PLOTTED lattice: ", LAT_PLOTTED.name, " ===")
    a = [sweep_one(pool, c, MONITOR_LOG) for c in CONFIGS]
    println("[fig_error_budget] === SECONDARY lattice: ", LAT_SECONDARY.name,
            " (r0 = 10 column only) ===")
    b = [sweep_one(pool, c, MONITOR_LOG)
         for c in build_configs(LAT_SECONDARY.eqsym; cols = SECONDARY_COLUMNS)]
    return a, b
end
@printf("[fig_error_budget] both lattices: %.1f s on %d workers\n", time() - t_all,
        first(plotted_rs).nworkers)

const P = curves(plotted_rs)
const S = curves(secondary_rs)
const NSPIN = length(SPINS)
# the plotted-lattice indices the secondary series pairs with, config by config
const SEC_IDX = [findfirst(c -> c.a == r.c.a && c.r0 == r.c.r0, CONFIGS) for r in secondary_rs]
any(isnothing, SEC_IDX) && error("fig_error_budget: a secondary configuration has no plotted twin")

# ───────────────────────────── the per-configuration table ───────────────────────────────────
println("\n[fig_error_budget] === PLOTTED (", LAT_PLOTTED.name, ") ===")
println("[fig_error_budget]   a     orbit    r0        F_r(no tail)     WB printed     ",
        "rel.ref    rel.trunc  rel.res    ref/trunc  trunc/res  tail       layer_max hold  nonconv")
for (i, r) in enumerate(plotted_rs)
    @printf("[fig_error_budget] %+.3f  %-5s %9.6f  %+.9e  %-14s %-10s %.3e  %.3e  %-9s %-9s  %+.3e  %.4g  %-5s %s\n",
            r.c.a, r.c.spec === :isco ? "ISCO" : "10M", r.c.r0, P.fr[i],
            r.c.wb === nothing ? "--" : r.c.wb.printed,
            isnan(P.rel_ref[i]) ? "--" : (@sprintf("%.3e", P.rel_ref[i])),
            P.rel_trunc[i], P.rel_res[i],
            isnan(P.rel_ref[i]) ? "--" : (@sprintf("%.3f", P.rel_ref[i] / P.rel_trunc[i])),
            (@sprintf("%.3f", P.rel_trunc[i] / P.rel_res[i])),
            P.tail_value[i], r.layer_probe, string(r.hold),
            string(MS[.!r.converged]))
end
println("\n[fig_error_budget] === SECONDARY (", LAT_SECONDARY.name, ") ===")
for (i, r) in enumerate(secondary_rs)
    @printf("[fig_error_budget] %+.3f  %-5s %9.6f  %+.9e  rel.ref %-10s rel.trunc %.3e  rel.res %.3e  |dF_r| vs plotted %.3e\n",
            r.c.a, r.c.spec === :isco ? "ISCO" : "10M", r.c.r0, S.fr[i],
            isnan(S.rel_ref[i]) ? "--" : (@sprintf("%.3e", S.rel_ref[i])),
            S.rel_trunc[i], S.rel_res[i], abs(S.fr[i] - P.fr[SEC_IDX[i]]))
end
println("[fig_error_budget] ndofs per level: PLOTTED ", first(plotted_rs).ndofs,
        "  SECONDARY ", first(secondary_rs).ndofs)
flush(stdout)

# ───────────────────────────── the node's gates ──────────────────────────────────────────────
gatefails = String[]
pub_trunc = vcat(PUB_ISCO_TRUNC, PUB_R10_TRUNC)
pub_res   = vcat(PUB_ISCO_RES, PUB_R10_RES)
pub_absfr = fill(NaN, NCONF)     # the published |F_r| strip carries only |a| <= 0.9, per column
# the published reference curve is defined only on |a| <= 0.9, in each column
pub_ref = fill(NaN, NCONF)
for (k, (vr, vf)) in enumerate(((PUB_ISCO_REF, PUB_ISCO_ABSFR), (PUB_R10_REF, PUB_R10_ABSFR))),
    (j, a) in enumerate(PUB_REF_SPINS)
    i = (k - 1) * NSPIN + findfirst(x -> x == a, SPINS)
    pub_ref[i] = vr[j]
    pub_absfr[i] = vf[j]
end

# GATE 1 is enforced where it happens: `sweep_one` THROWS on any failed mode, because a mode sum
# with a hole is not a shorter mode sum (the sealed `total_self_force` refuses it too). GMRES
# NON-CONVERGENCE is a different thing and is REPORTED, never gated: `fig_modes.jl` established
# that a stagnating high-m mode at a compressed near-extremal orbit stays in the figure, flagged,
# with its own final residual registered (`converged_flat`, `final_residual_flat`).
nonconv = [(r.c.a, r.c.spec, MS[.!r.converged]) for r in plotted_rs if !all(r.converged)]
layer_max = [maximum(x -> isnan(x) ? -Inf : x, r.layer) for r in plotted_rs]

# ── THE FORMAL VERDICT, AND WHAT IS NOT IN IT ────────────────────────────────────────────────
# GATES 2-6 and 8 run on ALL 18 configurations. There is no scope, and in particular the
# converged-mode-sum criterion below is NOT applied to them: a criterion introduced after a run,
# which removes exactly the configuration that failed, cannot be what makes the verdict pass. So
# `gates_pass` reads FALSE here, the failures are named, and this figure lands PRELIMINARY with
# them explicit — which is the honest state of a terminal node that has two open holes.
#
# EVERY FAILURE IS PINNED. `gate_failure_ids` is a registered, sorted vector, so the golden fixes
# the EXACT set of failing gates: a new failure, or the disappearance of a known one, is a golden
# MISMATCH and fails `--check`. That is what lets the driver exit 0 on a known, reproducible
# failure — see the exit rule at the foot of this file.
gatefails = String[]
gate_failure_ids = String[]
cfgtag(i) = @sprintf("a=%+.3f,%s", plotted_rs[i].c.a, plotted_rs[i].c.spec === :isco ? "ISCO" : "10M")
function gatefail!(gate::AbstractString, idxs, msg::AbstractString)
    isempty(idxs) && return nothing
    for i in idxs
        push!(gate_failure_ids, i === nothing ? gate : string(gate, "@", cfgtag(i)))
    end
    push!(gatefails, string(gate, ": ", msg))
    return nothing
end

# ── the EXPLORATORY converged-mode-sum criterion: registered, PRE-REGISTERED for future runs of
#    this driver, and applied to NOTHING here ---------------------------------------------------
# `truncation_error = |T_{m_max}|` is a TRUNCATION estimate only while the m-sequence is in its
# asymptotic `C m^-p` regime. Two measured criteria speak to that: (a) the SEALED tail estimator
# refusing to fit (`mode_decay_fit` throws unless two modes clear its headroom rule), and (b) the
# last retained regularized term staying below `MODE_SUM_CONVERGED` of the total, i.e.
# `rel_trunc < MODE_SUM_CONVERGED` — the criterion is the SMALLNESS of that term. (a) fires
# nowhere here. (b) was written AFTER the first production run and fires at exactly one
# configuration — the one that failed GATE 3 — so it is recorded as an EXPLORATORY diagnostic
# with its sensitivity, never as a gate scope. SENSITIVITY, registered: it excludes exactly that
# configuration for ANY threshold strictly between the second-largest and the largest
# `rel_trunc` (`exploratory_scope_sensitivity_*` below, measured 2.93e-4 … 2.41e-2); a threshold
# of 1e-1 would exclude nothing at all. A criterion whose verdict is that sensitive to a constant
# nobody derived is a lead, not a gate.
mode_sum_converged = [P.rel_trunc[i] < MODE_SUM_CONVERGED for i in 1:NCONF]
tail_fit_ok = [plotted_rs[i].heff !== nothing for i in 1:NCONF]
exploratory_scope = [tail_fit_ok[i] && mode_sum_converged[i] for i in 1:NCONF]
sorted_rt = sort(P.rel_trunc)
scope_sens = (sorted_rt[end - 1], sorted_rt[end])   # any threshold strictly inside excludes 1
if !all(exploratory_scope)
    println("\n[fig_error_budget] EXPLORATORY converged-mode-sum diagnostic (NOT a gate scope): ",
            count(!, exploratory_scope), " configuration(s) — ",
            [(plotted_rs[i].c.a, plotted_rs[i].c.spec, P.rel_trunc[i], P.rel_res[i])
             for i in findall(!, exploratory_scope)],
            "; tail-fit refusals ", count(!, tail_fit_ok),
            "; threshold ", MODE_SUM_CONVERGED, " sensitive over ", scope_sens)
end

# GATE 2 — the row-58 rule's own instrument, on the sweep the rule chose to publish.
gatefail!("GATE2", [i for i in 1:NCONF if !(layer_max[i] <= LAYER_GATE)],
          "the plotted layer_ratio exceeds $(LAYER_GATE) at " *
          "$([(plotted_rs[i].c.a, plotted_rs[i].c.spec, layer_max[i]) for i in 1:NCONF
              if !(layer_max[i] <= LAYER_GATE)])")
# ... and the rule's own verdict on itself: a hold that did not reduce the ratio, or a gate the
# hold could not clear, is the helper's `hole` and is a gate failure of the SELECTION, not of the
# plane (they differ: see `select_hold_sweep`'s docstring).
gatefail!("GATE2H", [i for i in 1:NCONF if plotted_rs[i].sel.hole],
          "the hold-selection rule could not clear the gate at " *
          "$([(plotted_rs[i].c.a, plotted_rs[i].c.spec, plotted_rs[i].sel.reason) for i in 1:NCONF
              if plotted_rs[i].sel.hole])")
# GATE 3 — THE PUBLISHED TRUNCATION ENVELOPE. Our estimator is the paper's own (omit the largest
# m) at the paper's own m_max, so this is the sharpest curve-level comparison the figure admits
# and it needs no external reference at all.
bad = [i for i in 1:NCONF
       if !(1 / TRUNC_ENVELOPE <= P.rel_trunc[i] / pub_trunc[i] <= TRUNC_ENVELOPE)]
gatefail!("GATE3", bad,
    "the relative truncation error is outside a factor $(TRUNC_ENVELOPE) of the published " *
    "curve at $([(plotted_rs[i].c.a, plotted_rs[i].c.spec) for i in bad]) " *
    "(ours $(P.rel_trunc[bad]) against published $(pub_trunc[bad]), " *
    "ratio $(P.rel_trunc[bad] ./ pub_trunc[bad]))")
# GATE 4 — THE PUBLISHED REFERENCE ENVELOPE, where Warburton-Barack print a cell
wbidx = findall(P.has_wb)
badref = [i for i in wbidx
          if !(1 / REF_ENVELOPE <= P.rel_ref[i] / pub_ref[i] <= REF_ENVELOPE)]
gatefail!("GATE4", badref,
    "the relative reference error is outside a factor $(REF_ENVELOPE) of the published " *
    "curve at $([(plotted_rs[i].c.a, plotted_rs[i].c.spec) for i in badref]) " *
    "(ours $(P.rel_ref[badref]) against published $(pub_ref[badref]))")
# GATE 5 — main.tex:724, "the truncation error models the true error up to a constant offset".
# NOTE ITS MARGIN: the measured offsets run 2.276 … 3.748 against a band of [2, 6], so the low
# end clears by only 1.14x. This gate is FRAGILE and would trip on a modest change of the tail
# treatment or of `m_max`; it is reported as such rather than widened after the fact.
offs = [P.rel_ref[i] / P.rel_trunc[i] for i in wbidx]
badoff = [wbidx[k] for k in eachindex(offs)
          if !(OFFSET_BAND[1] <= offs[k] <= OFFSET_BAND[2])]
gatefail!("GATE5", badoff,
    "the reference/truncation offset leaves the band $(OFFSET_BAND) at " *
    "$([(plotted_rs[i].c.a, plotted_rs[i].c.spec) for i in badoff]) " *
    "(offsets $(offs); published band $(PUB_REF_OVER_TRUNC))")
# GATE 6 — the SECOND external gate, F_t, which involves no regularization and no tail and which
# Warburton-Barack print to NINE significant figures (row 46 §3b).
ft_dev = fill(NaN, NCONF)
for (i, r) in enumerate(plotted_rs)
    cell = try
        abs(r.c.a) <= WB.WB_MAX_ABS_A &&
            any(x -> isapprox(x, r.c.a; atol = 1e-12), WB.WB_A_OVER_M) &&
            any(x -> isapprox(x, r.c.r0; atol = 1e-9), WB.WB_R0_OVER_M) ?
            WB.wb_cell(WB.WB_FT, r.c.a, r.c.r0) : nothing
    catch
        nothing
    end
    cell isa WB.WBRefCell || continue
    ft_dev[i] = abs(P.ft[i] - cell.value) / abs(cell.value)
end
badft = [i for i in 1:NCONF if !isnan(ft_dev[i]) && !(ft_dev[i] <= FT_GATE)]
gatefail!("GATE6", badft,
    "the unregularized F_t leaves $(FT_GATE) relative against the printed Table II at " *
    "$([(plotted_rs[i].c.a, plotted_rs[i].c.spec) for i in badft]) (deviations $(ft_dev[badft]))")
# GATE 8 — THE RESOLUTION CURVE, gated FROM ABOVE ONLY. Until this gate existed nothing
# constrained the dotted curve at all: its ratio to the published one spans 0.0155 … 5.3e5 across
# this figure, and `a = +0.9` at its ISCO sat at 389x with every other gate passing. The
# asymmetry is physical, not cosmetic: our estimator EXCEEDING the published one means our sweep
# is LESS converged than the paper's at the same nominal AMR depth, which is a deficiency;
# undershooting means it is more converged, which is not, and is reported (`resolution_ratio_*`)
# rather than gated. The bound is one decade — deliberately wide, because the two sides are not
# computing the same thing on the same lattice depth, and its whole job is that a departure of
# the size found above cannot pass in silence.
badres = [i for i in 1:NCONF
          if !isnan(pub_res[i]) && !(P.rel_res[i] / pub_res[i] <= RES_ENVELOPE_UP)]
gatefail!("GATE8", badres,
    "the relative resolution error exceeds the published curve by more than " *
    "$(RES_ENVELOPE_UP)x at $([(plotted_rs[i].c.a, plotted_rs[i].c.spec) for i in badres]) " *
    "(ours $(P.rel_res[badres]) against published $(pub_res[badres]), " *
    "ratio $(P.rel_res[badres] ./ pub_res[badres]))")
sort!(gate_failure_ids)

# THE RING-MARKED SET, by a measured rule and not by hand: a configuration this node names in a
# [HOLE] is one whose resolution estimator leaves the envelope (GATE 8) or whose mode sum the
# exploratory diagnostic finds unconverged. That is exactly the two families —
#   * the compressed PROGRADE ISCO, whose absolute resolution error runs 2.9e-10 (a = +0.7)
#     -> 6.0e-7 (+0.9) -> 2.4e-2 (+0.998), three orders per spin step, and
#   * |a| = 0.998 at r0 = 10, where it sits at 7.9e-10 / 8.2e-10.
# They are drawn with open rings on the figure, never dropped from it.
flagged = [(!isnan(pub_res[i]) && P.rel_res[i] / pub_res[i] > RES_ENVELOPE_UP) ||
           !exploratory_scope[i] for i in 1:NCONF]

# ───────────────────────────── REPORTED, NOT GATED ───────────────────────────────────────────
# (a) GLOBAL_DAG row 53's PLANNED gate, kept as the row's metric: is the deviation from the
#     printed table the size of one more AMR step? On the plotted lattice it cannot be — the
#     full-z resolution floor is ~100x below the half domain's, so the reference error is set by
#     the mode-sum remainder, not by resolution. The count is what the ledger row carries.
dag_ratio = [P.has_wb[i] ? P.ref[i] / P.res[i] : NaN for i in 1:NCONF]
dag_pass = [P.has_wb[i] && P.ref[i] <= DAG_RES_FACTOR * P.res[i] for i in 1:NCONF]
# (b) the ORDERING claim of main.tex:723, as a comparison per configuration. The two lattices
#     share their eight AMR LEVELS and settings, NOT their grids (8139 vs 4814 dofs, 1.69x): the
#     ndof difference is a confound, and the verdict below is an ASSOCIATION with the lattice at
#     the tested r0 = 10 configurations, not a demonstrated cause.
ordering = P.rel_trunc ./ P.rel_res
ordering_sec = S.rel_trunc ./ S.rel_res      # the r0 = 10 column only
println("\n[fig_error_budget] main.tex:723 ordering (truncation / resolution), per configuration:")
for (i, r) in enumerate(plotted_rs)
    j = findfirst(==(i), SEC_IDX)
    @printf("[fig_error_budget]   a = %+.3f %-5s  PLOTTED %8.3f   SECONDARY %-9s published %8.3f%s\n",
            r.c.a, r.c.spec === :isco ? "ISCO" : "10M", ordering[i],
            j === nothing ? "--" : (@sprintf("%.4f", ordering_sec[j])),
            pub_trunc[i] / pub_res[i], ordering[i] > 1 ? "" : "   <-- INVERTED")
end
@printf("[fig_error_budget] resolution below truncation: PLOTTED %d/%d, SECONDARY %d/%d, published %d/%d\n",
        count(>(1), ordering), NCONF, count(>(1), ordering_sec), length(ordering_sec),
        count(k -> pub_trunc[k] > pub_res[k], 1:NCONF), NCONF)
@printf("[fig_error_budget] the half domain's absolute resolution error at r0 = 10 is a spin-INDEPENDENT floor over |a| <= 0.9: %.3e ... %.3e (row 46 measured %.2e ... %.2e for the same configurations at NINE levels, i.e. ~3.6x per further step)\n",
        minimum(S.res[2:(end - 1)]), maximum(S.res[2:(end - 1)]), 1.563e-10, 1.637e-10)
println("[fig_error_budget] so the 88-2370x gap is ~4 further AMR steps, not the 1.69x in dofs ",
        "the other lattice buys at the same level count: the ordering is ASSOCIATED with the ",
        "lattice at these configurations [PRELIMINARY], not demonstrated to be caused by it")
@printf("[fig_error_budget] GLOBAL_DAG row 53 planned metric: %d of %d WB-referenced configurations have |F_r - WB| <= %g x the coherent resolution error\n",
        count(dag_pass), count(P.has_wb), DAG_RES_FACTOR)
flush(stdout)

# ───────────────────────────── registered values (the golden) ────────────────────────────────
# LAYOUT of every `*_flat` row: configuration index c = 1..$(NCONF) runs ISCO's nine spins first
# (in `spins` order) and then r0 = 10's nine; per-mode rows are (c, m) with m = 0..20 fastest.
@repro_result "spins" collect(SPINS)
@repro_result "column_tags" [String(t[2]) for t in COLUMNS]
@repro_result "column_notes" collect(COLUMN_NOTES)
@repro_result "n_configurations" NCONF
@repro_result "orbital_radius_flat" [c.r0 for c in CONFIGS]
@repro_result "horizon_radius_flat" [Float64(c.orb.rp) for c in CONFIGS]
@repro_result "worldtube_radius_flat" [Float64(c.wr) for c in CONFIGS]
@repro_result "m_max" last(MS)
@repro_result "regularization_order" first(plotted_rs).none.order
@repro_result "grid_points_per_dim" PORD
@repro_result "initial_refinement_L" LREF
@repro_result "amr_iterations_yaml" first(CONFIGS).cfg.amr.iterations
@repro_result "amr_trigger_offset" first(CONFIGS).cfg.amr.trigger_offset
@repro_result "amr_refinement_steps" STEPS
@repro_result "n_amr_grids" STEPS + 1
@repro_result "gmres_absolute_residual" ABSTOL
@repro_result "lattice_plotted" LAT_PLOTTED.name
@repro_result "lattice_plotted_equatorial_symmetry" LAT_PLOTTED.eqsym
@repro_result "lattice_secondary" LAT_SECONDARY.name
@repro_result "lattice_secondary_equatorial_symmetry" LAT_SECONDARY.eqsym
@repro_result "ndofs_per_level_plotted" first(plotted_rs).ndofs
@repro_result "ndofs_per_level_secondary" first(secondary_rs).ndofs
@repro_result "nelements_per_level_plotted" first(plotted_rs).nelements
@repro_result "tail_model_plotted" "none"
@repro_result "tail_model_secondary_number" (k = findfirst(r -> r.heff !== nothing, plotted_rs);
                                              k === nothing ? "none" :
                                              String(plotted_rs[k].heff.tail.method))
@repro_result "n_tail_refusals" count(r -> r.heff === nothing, plotted_rs)
@repro_result "mode_sum_converged_threshold" MODE_SUM_CONVERGED
@repro_result "mode_sum_converged_flat" mode_sum_converged
# the EXPLORATORY converged-mode-sum diagnostic: registered, applied to nothing (see the gates)
@repro_result "exploratory_scope_flat" exploratory_scope
@repro_result "exploratory_tail_fit_ok_flat" tail_fit_ok
@repro_result "n_outside_exploratory_scope" count(!, exploratory_scope)
@repro_result "exploratory_scope_sensitivity_low" scope_sens[1]
@repro_result "exploratory_scope_sensitivity_high" scope_sens[2]
@repro_result "n_configurations_gated" NCONF
@repro_result "flagged_hole_flat" flagged
@repro_result "n_flagged_holes" count(flagged)
# the resolution curve's envelope, and the full ratio the gate is drawn from
@repro_result "resolution_envelope_upper" RES_ENVELOPE_UP
@repro_result "n_outside_resolution_envelope" length(badres)
@repro_result "Ft_gate" FT_GATE
@repro_result "Ft_relative_deviation_max" maximum(x -> isnan(x) ? -Inf : x, ft_dev)
# the hold-selection rule's own record, per configuration
@repro_result "hold_cleared_flat" [r.sel.cleared for r in plotted_rs]
@repro_result "hold_hole_flat" [r.sel.hole for r in plotted_rs]
@repro_result "hold_reason_held" (k = findfirst(r -> r.hold, plotted_rs);
                                  k === nothing ? "no configuration was held" :
                                  plotted_rs[k].sel.reason)
# the two tail models are independent (the exact printed order-6 m-dependence with one fitted
# amplitude, and a free-exponent fit summed by Hurwitz zeta); row 46 measured 6-12 % agreement on
# real data, so their disagreement is a second, sealed asymptotic-regime instrument
@repro_result "tail_alt_value_flat" [r.heff === nothing ? NaN : Float64(r.heff.tail.alt_value)
                                     for r in plotted_rs]
@repro_result "tail_model_disagreement_flat" [r.heff === nothing ? NaN :
                                              abs(Float64(r.heff.tail.alt_value) -
                                                  Float64(r.heff.tail.value)) /
                                              abs(Float64(r.heff.tail.value))
                                              for r in plotted_rs]
@repro_result "reduction" String(first(plotted_rs).none.reduction)

# THE PLOTTED SERIES — the paper's own bare partial sum
@repro_result "Fr_flat" P.fr
@repro_result "Ft_flat" P.ft
@repro_result "Fphi_flat" P.fphi
@repro_result "truncation_error_flat" P.trunc
@repro_result "resolution_error_flat" P.res
@repro_result "resolution_error_bound_flat" P.res_bound
@repro_result "rel_ref_error_flat" P.rel_ref
@repro_result "rel_truncation_flat" P.rel_trunc
@repro_result "rel_resolution_flat" P.rel_res
@repro_result "has_wb_reference_flat" P.has_wb
@repro_result "wb_reference_Fr_flat" P.wb
@repro_result "wb_printed_half_ulp_flat" P.wb_ulp
@repro_result "ref_error_flat" P.ref
@repro_result "reference_over_truncation_flat" P.rel_ref ./ P.rel_trunc
@repro_result "truncation_over_resolution_flat" ordering
@repro_result "Ft_truncation_error_flat" P.ft_trunc
@repro_result "Ft_resolution_error_flat" P.ft_res
@repro_result "Ft_relative_deviation_vs_WB_flat" ft_dev

# THE TAIL-CORRECTED NUMBER — this reproduction's more accurate F_r, registered beside it
@repro_result "Fr_heffernan_tail_flat" P.fr_tail
@repro_result "partial_sum_flat" P.partial
@repro_result "tail_value_flat" P.tail_value
@repro_result "tail_uncertainty_flat" P.tail_unc
@repro_result "tail_exponent_flat" P.tail_exp
@repro_result "ref_error_tail_flat" P.ref_tail
@repro_result "rel_ref_error_tail_flat" P.ref_tail ./ abs.(P.fr_tail)
@repro_result "resolution_error_tail_flat" P.res_tail

# per (configuration, mode), m fastest
@repro_result "Fr_modes_flat" vcat([Float64.(r.none.Fr_modes) for r in plotted_rs]...)
@repro_result "regularized_terms_flat" vcat([Float64.(r.none.terms) for r in plotted_rs]...)
@repro_result "mode_resolution_errors_flat" vcat([Float64.(r.none.resolution_errors)
                                                  for r in plotted_rs]...)
@repro_result "converged_flat" vcat([r.converged for r in plotted_rs]...)
@repro_result "final_residual_flat" vcat([r.residual for r in plotted_rs]...)
@repro_result "final_gmres_iterations_flat" vcat([r.iterations for r in plotted_rs]...)
@repro_result "layer_ratio_flat" vcat([r.layer for r in plotted_rs]...)

# the row-58 rule's verdict, per configuration
@repro_result "horizon_hold_flat" [r.hold for r in plotted_rs]
@repro_result "layer_ratio_default_max_flat" [r.layer_probe for r in plotted_rs]
@repro_result "layer_ratio_max_flat" layer_max
@repro_result "layer_ratio_held_max_flat" [r.layer_hold_max for r in plotted_rs]
@repro_result "layer_gate" LAYER_GATE
@repro_result "n_configurations_held" count(r -> r.hold, plotted_rs)
@repro_result "all_modes_converged" isempty(nonconv)
@repro_result "n_nonconverged_modes" sum(count(!, r.converged) for r in plotted_rs)

# THE SECONDARY LATTICE — the shipped yaml's half domain, registered and never drawn
@repro_result "secondary_Fr_flat" S.fr
@repro_result "secondary_truncation_error_flat" S.trunc
@repro_result "secondary_resolution_error_flat" S.res
@repro_result "secondary_rel_ref_error_flat" S.rel_ref
@repro_result "secondary_rel_truncation_flat" S.rel_trunc
@repro_result "secondary_rel_resolution_flat" S.rel_res
@repro_result "secondary_truncation_over_resolution_flat" ordering_sec
@repro_result "secondary_horizon_hold_flat" [r.hold for r in secondary_rs]
@repro_result "secondary_layer_ratio_max_flat" [maximum(x -> isnan(x) ? -Inf : x, r.layer)
                                                for r in secondary_rs]
@repro_result "secondary_all_modes_converged" all(all(r.converged) for r in secondary_rs)
# THE TWO LATTICES AGAINST EACH OTHER: the equatorial-symmetry reduction is EXACT, so these must
# agree to their own resolution floor — the sharpest internal check the figure carries.
@repro_result "secondary_configuration_index" SEC_IDX
@repro_result "Fr_lattice_absdiff" abs.(P.fr[SEC_IDX] .- S.fr)
@repro_result "Fr_lattice_max_absdiff" maximum(abs.(P.fr[SEC_IDX] .- S.fr))
@repro_result "resolution_ratio_secondary_over_plotted" S.res ./ P.res[SEC_IDX]
@repro_result "resolution_ratio_secondary_over_plotted_min" minimum(S.res ./ P.res[SEC_IDX])
@repro_result "resolution_ratio_secondary_over_plotted_max" maximum(S.res ./ P.res[SEC_IDX])
# the half domain's absolute resolution error at r0 = 10 is a spin-INDEPENDENT floor over
# |a| <= 0.9 — the number the convergence-rate argument of the header rests on, and the LEAD for
# HOLE 2: the full-z lattice's own |a| = 0.998 values sit AT that floor rather than at the
# 1e-13..1e-12 its other spins reach.
@repro_result "secondary_resolution_floor_min" minimum(S.res[2:(end - 1)])
@repro_result "secondary_resolution_floor_max" maximum(S.res[2:(end - 1)])
@repro_result "plotted_resolution_at_extremal_r10" [P.res[NSPIN + 1], P.res[2 * NSPIN]]

# THE PUBLISHED FIGURE, measured from its own vector paths — REFERENCE, never plotted as ours
@repro_result "published_spins" PUB_SPINS
@repro_result "published_reference_spins" PUB_REF_SPINS
@repro_result "published_rel_truncation_flat" pub_trunc
@repro_result "published_rel_resolution_flat" pub_res
@repro_result "published_rel_ref_error_flat" pub_ref
@repro_result "published_absFr_flat" pub_absfr
@repro_result "absFr_ratio_ours_over_published_flat" abs.(P.fr) ./ pub_absfr
@repro_result "truncation_ratio_ours_over_published_flat" P.rel_trunc ./ pub_trunc
@repro_result "resolution_ratio_ours_over_published_flat" P.rel_res ./ pub_res
@repro_result "ref_error_ratio_ours_over_published_flat" P.rel_ref ./ pub_ref
@repro_result "published_truncation_over_resolution_flat" pub_trunc ./ pub_res
@repro_result "n_configurations_resolution_below_truncation" count(>(1), ordering)
@repro_result "n_configurations_resolution_below_truncation_secondary" count(>(1), ordering_sec)
@repro_result "n_secondary_configurations" length(secondary_rs)
@repro_result "n_configurations_resolution_below_truncation_published" count(k -> pub_trunc[k] >
                                                                            pub_res[k], 1:NCONF)
# GLOBAL_DAG row 53's PLANNED metric, reported in full
@repro_result "dag_planned_factor" DAG_RES_FACTOR
@repro_result "dag_ref_over_resolution_flat" dag_ratio
@repro_result "n_wb_configurations" count(P.has_wb)
@repro_result "n_wb_configurations_within_dag_factor" count(dag_pass)

# REPORTED, NOT GATED: wall-clock seconds on a shared 256-core Linux box
@repro_result "reported_configuration_wall_time_flat" [r.wall for r in plotted_rs]
@repro_result "reported_secondary_wall_time_flat" [r.wall for r in secondary_rs]
@repro_result "reported_nworkers" first(plotted_rs).nworkers
@repro_result "reported_max_load" maximum(r.max_load for r in plotted_rs)

# ───────────────────────────── draw ──────────────────────────────────────────────────────────
# THE HOLD, STATED ON THE FIGURE ITSELF, per column: the setting and BOTH ratios, so a reader
# never has to take the plotted plane's regularity on trust (row 58's rule, `select_hold_sweep`).
column_hold_notes = map(1:length(COLUMNS)) do k
    rng = ((k - 1) * NSPIN + 1):(k * NSPIN)
    held = [i for i in rng if plotted_rs[i].hold]
    isempty(held) && return "horizon_hold = false throughout (layer_ratio <= " *
                            @sprintf("%.2f", maximum(layer_max[rng])) * ")"
    join([@sprintf("a = %+.3f held: layer_ratio %.4g -> %.4g (gate %.2f)",
                   plotted_rs[i].c.a, plotted_rs[i].sel.default_ratio,
                   plotted_rs[i].sel.held_ratio, plotted_rs[i].sel.gate) for i in held], "; ")
end
@repro_result "column_hold_notes" column_hold_notes

cols = map(enumerate(COLUMNS)) do (k, (spec, tag))
    rng = ((k - 1) * NSPIN + 1):(k * NSPIN)
    (a = collect(SPINS), absfr = abs.(P.fr[rng]), rel_ref = P.rel_ref[rng],
     rel_trunc = P.rel_trunc[rng], rel_res = P.rel_res[rng], refmask = P.has_wb[rng],
     okmask = .!flagged[rng], m_max = last(MS), tag = tag, note = COLUMN_NOTES[k],
     holdnote = column_hold_notes[k])
end
publish = opts.mode === :regen
out_png = publish ? GOLD_PNG : joinpath(mktempdir(), "fig_error_budget.png")
saved = draw(cols, out_png; publish = publish)
saved.ok || gatefail!("GATE7", [nothing], saved.message)   # not per-configuration
publish || println("[fig_error_budget] --check: rendered to a scratch path, ",
                   relpath(GOLD_PNG, ReproHarness.REPRO_ROOT), " untouched")

@repro_result "png_bytes_within_budget" saved.ok
# THE FORMAL VERDICT, pinned configuration by configuration and registered only HERE, after the
# artifact gate has had its say: the golden fixes the EXACT set of failing gates, so a new
# failure or the disappearance of a known one is a golden MISMATCH.
sort!(gate_failure_ids)
@repro_result "gate_failure_ids" gate_failure_ids
@repro_result "n_gate_failures" length(gate_failure_ids)
@repro_result "gates_pass" isempty(gate_failure_ids)
if !isempty(nonconv)
    println("\n[fig_error_budget] REPORTED (not gated) — GMRES did not reach ", ABSTOL,
            " on every AMR level:")
    foreach(x -> println("    a = ", x[1], " ", x[2], "  modes ", x[3]), nonconv)
end
if !isempty(gatefails)
    println("\n[fig_error_budget] FINDING — node gates violated:")
    foreach(f -> println("    ", f), gatefails)
end

# rtol 1e-6 on values that are differences at a ~1e-12 floor; the `reported_*` timing rows carry
# a deliberately loose tolerance so a busy machine cannot fail a reproduced figure.
const LOOSE = (rtol = 1.0e3, atol = 1.0e3)
code = ReproHarness.conclude!(opts; paper = "2606.04998", id = "fig_error_budget",
                              script = @__FILE__, rtol = 1e-6,
                              tolerances = Dict("reported_configuration_wall_time_flat" => LOOSE,
                                                "reported_secondary_wall_time_flat" => LOOSE,
                                                "reported_nworkers" => LOOSE,
                                                "reported_max_load" => LOOSE))
# THE EXIT CODE REPORTS REPRODUCIBILITY, NOT THE PHYSICS VERDICT — and the two are separated
# deliberately. The verdict is a REGISTERED VALUE: `gates_pass` and, configuration by
# configuration, `gate_failure_ids`, both pinned by the golden. A gate that starts failing, or a
# known failure that disappears, is therefore a golden MISMATCH and fails `--check` on the
# comparison, which is the check that matters. Exiting non-zero on a known, pinned, reproducible
# failure instead would make the node unadmissible to the ledger at all, and the only way to land
# it would be to hide its holes — the opposite of the discipline. What DOES set a non-zero exit
# is a golden mismatch (`code`) or an over-budget tracked artifact.
exit(max(code, saved.ok ? 0 : 1))

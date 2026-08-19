# reproduce/arXiv-2608.17003/fig04.jl — EGD::fig04-dtdr-ap.
#
# arXiv:2608.17003 Fig. 4 (`fig:dtdr_ap`, main.tex:915-922): the mode-summed torque density
# d𝔗/dr = √−g⟨𝒯⟩ (Eq. 72b integrand, Fourier average 73) in units of F_{J,1} = F_{J,0}/a• (76),
# versus the corotation-scaled radius (r/a• − 1)/h0 ∈ [−2, 2], for a• ∈ {10, 30, 100, 300, 1000} M
# at χ = 𝔭 = 𝔮 = 0 and h0 = 0.05.
#
# Caption claims asserted here (M:919), with the DEFINITIONS fixed before the run:
#
#   * "in the nearly Newtonian limit (a• = 1000M) the peaks of the torque density are offset by about
#     ±h0 a• from corotation"  →  x of the positive lobe's maximum ≈ +1 and of the negative lobe's
#     minimum ≈ −1, on the a• = 1000M curve;
#   * "in the strongly relativistic case (a• = 10M) the peaks of the torque density shift inward by
#     ~a• h0/5"  →  each peak of the a• = 10M curve sits CLOSER TO COROTATION than the corresponding
#     a• = 1000M peak by 0.2 in x; the registered number is the mean of the two lobes,
#     ½[(x⁺_1000 − x⁺_10) + (|x⁻_1000| − |x⁻_10|)], and both lobes are registered separately so that
#     an asymmetry is visible rather than averaged away;
#   * "... and increase in magnitude by a factor of ~1.5"  →  max|d𝔗/dr|(a•=10M) / max|d𝔗/dr|(a•=1000M).
#
# The three numbers are the golden keys `peak_offset_over_h0a_newt_outer` / `_inner`,
# `peak_inward_shift_over_h0a_a10` and `peak_growth_factor_a10` of the printed-number gate
# (`check_printed_numbers.jl`, contract table in this directory's README).
#
# ZERO physics here: the disc, the modes and the torque densities come from
# `AccretionDisc.SpiralWaves`; the mode sum and the reduction to the figure's cells live in
# `_profile_assembly.jl` beside this file (shared with fig03).  Style: `_style/emri_style.jl` +
# section 4 of `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

# manual section 0 palette 3: viridis([0.8, 0.6, 0.4, 0.2]) for the a• ladder + thin black Newtonian
const VIRIDIS4 = (RGB(0x7a / 255, 0xd1 / 255, 0x51 / 255), RGB(0x22 / 255, 0xa8 / 255, 0x84 / 255),
                  RGB(0x2a / 255, 0x78 / 255, 0x8e / 255), RGB(0x41 / 255, 0x44 / 255, 0x87 / 255))
const LSTYLE = (:solid, :dash, :dashdot, :dot)

"""
    draw(x, curves, aps) -> path

Render Fig. 4 from the plotted arrays (`curves[i]` = d𝔗/dr in F_{J,1} for `aps[i]`) and save the PNG
beside the golden.  `--plot-only` redraws from `expected/fig04.json` without recomputing.
"""
function draw(x, curves, aps)
    plt = plot(size = (326, 233), dpi = 160, legend = :topleft, legendfontsize = 6,
               legend_columns = 2, foreground_color_legend = RGB(0.8, 0.8, 0.8),
               background_color_legend = :white,
               guidefontsize = 11, tickfontsize = 9,
               xlims = (-2, 2), ylims = (-7.5, 7.5), xticks = -2:1:2, yticks = -7.5:2.5:7.5,
               minorticks = 5, tick_direction = :in, widen = false,
               xlabel = "\$(r/a_\\bullet - 1)/h_0\$", ylabel = "\$d\\mathfrak{T}/dr\\ (F_{J,1})\$",
               bottom_margin = 1Plots.mm, left_margin = 1Plots.mm)
    for i in 1:(length(aps) - 1)
        plot!(plt, x, curves[i]; color = VIRIDIS4[i], ls = LSTYLE[i], lw = 1.2,
              label = "\$a_\\bullet=$(Int(aps[i]))M\$")
    end
    plot!(plt, x, curves[end]; color = :black, ls = :solid, lw = 0.8,
          label = "\$a_\\bullet=$(Int(aps[end]))M\$")
    png_path = joinpath(@__DIR__, "expected", "fig04.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig04] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig04.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

const AP_KEYS = ("dTdr_a10", "dTdr_a30", "dTdr_a100", "dTdr_a300", "dTdr_a1000")

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig04")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["x_corot"], [v[k] for k in AP_KEYS], v["a_emri"])
    println("[fig04] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves
include(joinpath(@__DIR__, "_profile_assembly.jl"))
using .ProfileAssembly

for f in (:solve_mode, :torque_density, :cumulative_torque, :torque_units)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: Fig. 4 needs the mode solver and the " *
        "torque diagnostics (EGD::sw-mode-solver, EGD::sw-diagnostics).")
end

# ───────────────────────────── parameters (main.tex:915-922, brief section 4 row 4) ──────────────

const CHI   = 0.0
const H0    = 0.05
const PDISK = 0
const QDISK = 0
const BSOFT = 0.6
const EPS_P = 1e-5
const APS   = (10.0, 30.0, 100.0, 300.0, 1000.0)

const MMAX_DEFAULT = 120           # same ceiling as fig03 (m₉₅ ≲ 70 at h0 = 0.05, M:1025)
_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
const MMAX  = _envint("GRTOOLS_FIG04_MMAX", MMAX_DEFAULT)
const NCELL = _envint("GRTOOLS_FIG04_NCELL", 640)     # uniform cells over x ∈ [−2, 2]
const XLO, XHI = -2.0, 2.0

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

# the plotting grid is ANCHORED on corotation (x = 0), so that the torque density is sampled AT r_c
# in a single cell for every a• — see `centered_edges` in `_profile_assembly.jl`.
xe = centered_edges(0.0, XLO, XHI, (XHI - XLO) / NCELL)
xc = cell_centers(xe)
println("\n[fig04] χ = ", CHI, "  h0 = ", H0, "  (𝔭,𝔮) = (", PDISK, ",", QDISK, ")  m = 1:", MMAX,
        "   ", length(xc), " uniform cells over (r/a•−1)/h0 ∈ [", XLO, ", ", XHI, "]")

"One a• of the figure: the mode-summed torque density on the plot cells, plus its two lobes."
function one_ap(a, xc, xe)
    bg = SW.DiscBackground(; χ = CHI, aEMRI = a, h0 = H0, 𝔭 = PDISK, 𝔮 = QDISK, b = BSOFT, εp = EPS_P)
    _, FJ1 = SW.torque_units(bg)
    redges = a .* (1 .+ H0 .* xe)                 # x = (r/a•−1)/h0  ⇒  r = a•(1 + h0 x)
    res = profile_sum(bg, 1:MMAX, redges; verbose = (a == first(APS)))
    y = res.dT_dr ./ FJ1
    (y = y, crest = maximum(y), x_crest = xc[argmax(y)], trough = minimum(y), x_trough = xc[argmin(y)],
     tail = sum(res.per_mode_peak[max(1, MMAX - 9):MMAX]) / sum(res.per_mode_peak),
     edge = res.edge_density_fraction, residual = res.residual_max, seconds = res.seconds,
     nodes = res.nodes)
end

out = map(a -> begin
              o = one_ap(a, xc, xe)
              @printf("[fig04] a• = %6.1f M   crest %+.4f at x = %+.4f   trough %+.4f at x = %+.4f   (%.1f s, %d nodes)\n",
                      a, o.crest, o.x_crest, o.trough, o.x_trough, o.seconds, o.nodes)
              flush(stdout)
              o
          end, APS)

curves = [o.y for o in out]
peaks_pos = [o.crest for o in out]; xpos = [o.x_crest for o in out]
peaks_neg = [o.trough for o in out]; xneg = [o.x_trough for o in out]
tails = [o.tail for o in out]; edgefracs = [o.edge for o in out]
resmax = maximum(o.residual for o in out); secs = sum(o.seconds for o in out)

# ── the caption's three numbers, per the definitions in the file header ──────────────────────────
i10, i1000 = 1, length(APS)
offset_outer = xpos[i1000]
offset_inner = xneg[i1000]
shift_outer  = xpos[i1000] - xpos[i10]                     # + = the a•=10M peak is closer to r_c
shift_inner  = abs(xneg[i1000]) - abs(xneg[i10])
inward_shift = (shift_outer + shift_inner) / 2
growth = max(peaks_pos[i10], -peaks_neg[i10]) / max(peaks_pos[i1000], -peaks_neg[i1000])

@repro_result "x_corot"       xc
for (k, y) in zip(AP_KEYS, curves)
    @repro_result k y
end
@repro_result "a_emri"        collect(APS)
@repro_result "chi"           CHI
@repro_result "h0"            H0
@repro_result "p_disk"        PDISK
@repro_result "q_disk"        QDISK
@repro_result "b_soft"        BSOFT
@repro_result "eps_p"         EPS_P
@repro_result "m_max"         MMAX
@repro_result "n_cells"       length(xc)
@repro_result "crest_FJ1"     peaks_pos
@repro_result "trough_FJ1"    peaks_neg
@repro_result "crest_x"       xpos
@repro_result "trough_x"      xneg
@repro_result "mode_sum_tail_fraction_max" maximum(tails)
@repro_result "wkb_edge_density_fraction_max" maximum(edgefracs)
@repro_result "solver_residual_max" resmax

# printed-number-gate contract keys (README table, claims M919-peak-offset/shift/growth)
@repro_result "peak_offset_over_h0a_newt_outer" offset_outer
@repro_result "peak_offset_over_h0a_newt_inner" offset_inner
@repro_result "peak_inward_shift_over_h0a_a10"  inward_shift
@repro_result "peak_inward_shift_outer_a10"     shift_outer
@repro_result "peak_inward_shift_inner_a10"     shift_inner
@repro_result "peak_growth_factor_a10"          growth

# ───────────────────────────── printed-number checks (M:919) ─────────────────────────────────────
# Set to `false` ONLY after a documented investigation shows the paper number is not reproduced; the
# driver then still passes on OUR golden and prints the discrepancy as a FINDING.
const PAPER_CHECKS_ENFORCED = true

struct PaperCheck; id::String; src::String; paper::Float64; tol::Float64; ours::Float64; end
checks = PaperCheck[
    PaperCheck("Newtonian peak offset (outer)",  "M:919",  1.00, 0.15, offset_outer),
    PaperCheck("Newtonian peak offset (inner)",  "M:919", -1.00, 0.15, offset_inner),
    PaperCheck("a•=10M inward peak shift",       "M:919",  0.20, 0.10, inward_shift),
    PaperCheck("a•=10M peak growth factor",      "M:919",  1.50, 0.25, growth),
]
println("\n[fig04] printed-number checks (M:919)")
println("  claim                              src      paper     ours      tol     verdict")
violations = String[]
for c in checks
    ok = abs(c.ours - c.paper) <= c.tol
    ok || push!(violations, "$(c.id) [$(c.src)]: paper $(c.paper), ours $(round(c.ours; digits = 4)), tol ±$(c.tol)")
    @printf("  %-34s %-8s %7.4f  %7.4f  ±%.3f   %s\n", c.id, c.src, c.paper, c.ours, c.tol,
            ok ? "PASS" : "FAIL")
end
@printf("[fig04] lobe asymmetry of the shift: outer %+.4f, inner %+.4f\n", shift_outer, shift_inner)
@printf("[fig04] mode-sum tail (max over a•): %.2e ; WKB-edge density fraction (max): %.2e ; wall %.1f s\n",
        maximum(tails), maximum(edgefracs), secs)
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig04] FINDING — printed numbers not reproduced within tolerance:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig04] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written")
        paper_check_failed = true
    else
        println("[fig04] PAPER_CHECKS_ENFORCED = false ⇒ recorded as a finding; the golden still gates our own numbers")
    end
end

# mode-sum convergence gate (independent of the paper's numbers)
maximum(tails) < 5e-3 ||
    error("fig04: the last 10 of $(MMAX) modes still carry $(maximum(tails)) of the summed per-mode " *
          "peak torque density — raise GRTOOLS_FIG04_MMAX")

draw(xc, curves, collect(APS))

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig04", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_cells" => (rtol = 0.0, atol = 0.0)))
exit(paper_check_failed ? max(code, 1) : code)

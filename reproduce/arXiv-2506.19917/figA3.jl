# reproduce/arXiv-2506.19917/figA3.jl — ECC::fd25-figA3.
#
# arXiv:2506.19917 Fig. A3 (`fig:dTdr_convergence_sim`, Sections/app_modal_convergence.tex:44-63,
# source figure src/Figures/dTdr_convergence_sim.pdf, page 21): the Disco resolution test for
# (p, q) = (0.5, 0.0), e = 0.5 — (a) fiducial (black solid) vs double resolution (red dashed) radial
# torque density dT/dr, (b) the residual Δ dT/dr = fiducial − high resolution (eq A2), (c) the
# cumulative torque error ΔT integrated out from r = 0.1 a_p.
#
# SIM ONLY — DIGITISED REFERENCE, NEVER A REPRODUCED RESULT.  Every curve of this figure is a Disco
# hydrodynamics run, and the mission default is that the simulations are not re-run (FD25 style
# manual §9, provenance SIM).  This driver therefore contains ZERO physics and calls no GrTools
# package: it re-plots, in the FD25 style, the series that ECC::fd25-vector-digitise read off the
# vector paths of the published PDF (axes calibrated on the tick marks, worst refit 1.56e-3 pt), and
# registers summary values of the plotted series so that `--check` detects any change of the inputs.
#
# INPUTS (emri-disk mission repo; the sibling of this GRTools checkout by default, or the directory
# named by GRTOOLS_EMRI_DISK).  Absent inputs are a documented HOLE (exit 3), never a fabricated plot:
#   results/emri-eccentric/paper_EMRI-ECC/evidence/fd25-vector-digitise/digitised/
#       figA3_dTdr_convergence_sim.csv   the digitised series (columns panel,series,...,x,y,...)
#       manifest.json                    per-panel axes rectangles and tick calibrations (points)
#   results/emri-eccentric/paper_EMRI-ECC/evidence/fd25-figA3/figA3_reference.json
#       the reference package: units, caption, series counts/ranges, the CSV sha256 it was built on
# The sha256 of all three are registered (exact string compare), and the figure checks require the
# CSV to be the bytes the reference package names and every series count/range to match it.
#
# AXES are measured, not chosen: each panel's data limits are the manifest's axes rectangle mapped
# through that panel's tick calibration (x 0.1 … 2.0; (a) matplotlib's 5 % autoscale margin,
# (b) ±0.15, (c) ±8e-4).  Style: `_style/manuals/arXiv-2506.19917.md` §0 and §8 (A3/A4 stacks).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf
using SHA
using DelimitedFiles

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

const PAPER = "2506.19917"
const ID = "fd25_figA3"
const EMRI_DISK = get(ENV, "GRTOOLS_EMRI_DISK",
                      normpath(joinpath(ReproHarness.REPRO_ROOT, "..", "..", "emri-disk")))
const EVID = joinpath(EMRI_DISK, "results", "emri-eccentric", "paper_EMRI-ECC", "evidence")
const CSV_PATH = joinpath(EVID, "fd25-vector-digitise", "digitised", "figA3_dTdr_convergence_sim.csv")
const MANIFEST_PATH = joinpath(EVID, "fd25-vector-digitise", "digitised", "manifest.json")
const REFERENCE_PATH = joinpath(EVID, "fd25-figA3", "figA3_reference.json")

for p in (CSV_PATH, MANIFEST_PATH, REFERENCE_PATH)
    isfile(p) || ReproHarness.hole("Fig A3 is a digitised SIM reference; input $(p) is absent (set " *
                                   "GRTOOLS_EMRI_DISK to the emri-disk mission repo)")
end

sha_hex(path) = bytes2hex(open(sha256, path))

# ───────────────────────────── inputs (reference data only) ──────────────────────────────────────

const REF = ReproHarness.json_read(read(REFERENCE_PATH, String))
const MAN = only(f for f in ReproHarness.json_read(read(MANIFEST_PATH, String))["figures"]
                 if f["figure"] == "A3")

raw, hdr = readdlm(CSV_PATH, ','; header = true)
const COL = Dict(String(h) => k for (k, h) in enumerate(vec(hdr)))
const KEYS = ("a_dTdr/fiducial", "a_dTdr/double_res", "b_delta/fiducial", "b_delta/zero_line",
              "c_cumulative/fiducial", "c_cumulative/zero_line")

"The digitised series `panel/series` in drawing order (x, y as Float64)."
function series(key)
    panel, name = split(key, "/")
    rows = [k for k in axes(raw, 1) if raw[k, COL["panel"]] == panel && raw[k, COL["series"]] == name]
    order = sortperm(Int.(raw[rows, COL["i"]]))
    (x = Float64.(raw[rows, COL["x"]])[order], y = Float64.(raw[rows, COL["y"]])[order])
end
const S = Dict(k => series(k) for k in KEYS)

"Data limits of a panel: its axes rectangle (pt) through its own tick calibration."
function panel_limits(name)
    p = only(q for q in MAN["panels"] if q["name"] == name)
    x0, y0, x1, y1 = Float64.(p["rect_pt"])
    cal(ax, v) = Float64(p[ax]["fit_a"]) + Float64(p[ax]["fit_b_per_pt"]) * v
    (x = (cal("x", x0), cal("x", x1)), y = extrema((cal("y", y0), cal("y", y1))),
     yticks = Float64.(p["y"]["tick_values"]), xticks = Float64.(p["x"]["tick_values"]))
end
const LIM = Dict(n => panel_limits(n) for n in ("a_dTdr", "b_delta", "c_cumulative"))

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

const RED = RGB(1.0, 0.0, 0.0)                         # A3 double resolution #ff0000 (manual §0)
const ZERO_GREY = RGB(0x80 / 255, 0x80 / 255, 0x80 / 255)  # zero lines #808080 dashed (manual §0)
const GAINSBORO = RGB(0xdc / 255, 0xdc / 255, 0xdc / 255)  # rounded panel-letter boxes (manual §0)
const UNIT_DTDR = "[q_\\mathrm{r}^2 h_\\mathrm{p}^{-3} a_\\mathrm{p}^{-1}]"
const YLABELS = Dict("a_dTdr" => "\$dT/dr\\ $(UNIT_DTDR)\$",
                     "b_delta" => "\$\\Delta dT/dr\\ $(UNIT_DTDR)\$",
                     "c_cumulative" => "\$\\Delta T\\ [q_\\mathrm{r}^2 h_\\mathrm{p}^{-3}]\$")

"Tick labels as LaTeX with a true minus sign and the published number of decimals."
tick_labels(vals, digits) = [v == 0 ? "\$" * (digits == 0 ? "0" : "0." * "0"^digits) * "\$" :
                             "\$" * (v < 0 ? "-" : "") * (digits == 0 ? string(round(Int, abs(v))) :
                             @sprintf("%.*f", digits, abs(v))) * "\$" for v in vals]

"Rounded gainsboro box with an italic panel letter, lower right (axes fractions as in the PDF)."
function letter!(plt, lim, letter)
    (xa, xb), (ya, yb) = lim.x, lim.y
    fx(f) = xa + f * (xb - xa); fy(f) = ya + f * (yb - ya)
    u0, u1, v0, v1, rr = 0.887, 0.986, 0.075, 0.278, 0.022       # box edges; corner radius (x-fraction)
    rv = rr * 368 / (408 / 3)                                     # same radius on the ~3:1 panel
    θ = range(0, π / 2; length = 8)
    corner(cu, cv, a0) = [(cu + rr * cos(a0 + t), cv + rv * sin(a0 + t)) for t in θ]
    pts = vcat(corner(u1 - rr, v1 - rv, 0.0), corner(u0 + rr, v1 - rv, π / 2),
               corner(u0 + rr, v0 + rv, π), corner(u1 - rr, v0 + rv, 3π / 2))
    plot!(plt, Shape([fx(p[1]) for p in pts], [fy(p[2]) for p in pts]);
          fillcolor = GAINSBORO, linecolor = GAINSBORO, linewidth = 0, label = "")
    annotate!(plt, fx((u0 + u1) / 2), fy((v0 + v1) / 2), text("\$($(letter))\$", 10, :black))
end

function panel(name; bottom = false)
    lim = LIM[name]
    digits = name == "a_dTdr" ? 0 : name == "b_delta" ? 1 : 4
    # Ticks on the left and bottom only (manual §0): GR's :box frame mirrors them onto the top and
    # right spines, so the frame is :axes and the top/right spines are drawn as plain lines.
    plt = plot(; xlims = lim.x, ylims = lim.y, widen = false, framestyle = :axes, grid = false,
               tick_direction = :out, xticks = (lim.xticks, bottom ? tick_labels(lim.xticks, 1) : fill("", 4)),
               yticks = (lim.yticks, tick_labels(lim.yticks, digits)), ylabel = YLABELS[name],
               xlabel = bottom ? "\$r/a_\\mathrm{p}\$" : "", guidefontsize = 9, tickfontsize = 9,
               legend = false, top_margin = (name == "a_dTdr" ? 0.5 : 0.0) * Plots.mm,
               bottom_margin = bottom ? 0Plots.mm : -4.7Plots.mm)
    # drawn ON the plot-area edge, so half of the stroke is clipped: width 2 shows as the axes' 1
    plot!(plt, [lim.x[1], lim.x[2], lim.x[2]], [lim.y[2], lim.y[2], lim.y[1]]; color = :black,
          linewidth = 2.0, label = "")
    plt
end

function draw(png)
    pa = panel("a_dTdr")
    plot!(pa, S["a_dTdr/fiducial"].x, S["a_dTdr/fiducial"].y; color = :black, linewidth = 1.3,
          label = "Fiducial resolution")
    plot!(pa, S["a_dTdr/double_res"].x, S["a_dTdr/double_res"].y; color = RED, linestyle = :dash,
          linewidth = 1.3, label = "Double resolution")
    plot!(pa; legend = :bottomleft, legendfontsize = 7, foreground_color_legend = RGB(0.8, 0.8, 0.8),
          background_color_legend = :white)
    letter!(pa, LIM["a_dTdr"], "a")
    pb = panel("b_delta")
    plot!(pb, S["b_delta/fiducial"].x, S["b_delta/fiducial"].y; color = :black, linewidth = 1.3, label = "")
    plot!(pb, S["b_delta/zero_line"].x, S["b_delta/zero_line"].y; color = ZERO_GREY, linestyle = :dash,
          linewidth = 1.6, label = "")
    letter!(pb, LIM["b_delta"], "b")
    pc = panel("c_cumulative"; bottom = true)
    plot!(pc, S["c_cumulative/fiducial"].x, S["c_cumulative/fiducial"].y; color = :black, linewidth = 1.3,
          label = "")
    plot!(pc, S["c_cumulative/zero_line"].x, S["c_cumulative/zero_line"].y; color = ZERO_GREY,
          linestyle = :dash, linewidth = 1.6, label = "")
    letter!(pc, LIM["c_cumulative"], "c")
    plt = plot(pa, pb, pc; layout = grid(3, 1, heights = [0.3326, 0.3326, 0.3348]), size = (368, 408), dpi = 200,
               left_margin = 1.1Plots.mm, right_margin = 2.9Plots.mm, link = :x)
    save_figure!(plt, png; id = ID)
end

# ───────────────────────────── registered summary of the plotted series ──────────────────────────

violations = String[]
check(ok, msg) = (ok || push!(violations, msg); ok)

@repro_result "input_csv_sha256" sha_hex(CSV_PATH)
@repro_result "input_manifest_sha256" sha_hex(MANIFEST_PATH)
@repro_result "input_reference_sha256" sha_hex(REFERENCE_PATH)
@repro_result "reference_status" String(REF["status"])
check(sha_hex(CSV_PATH) == REF["series_file"]["sha256"],
      "the CSV is not the bytes figA3_reference.json was built on ($(REF["series_file"]["sha256"]))")
check(REF["all_gates_pass"] === true, "figA3_reference.json does not report all_gates_pass")

for k in KEYS
    s = S[k]; tag = replace(k, "/" => "_")
    ref = REF["series"][k]
    n = length(s.x); xlo, xhi = extrema(s.x); ylo, yhi = extrema(s.y)
    @repro_result "$(tag)_n" n
    @repro_result "$(tag)_x_range" [xlo, xhi]
    @repro_result "$(tag)_y_range" [ylo, yhi]
    check(n == ref["n"], "$(k): $(n) points, reference package says $(ref["n"])")
    check(isapprox([xlo, xhi], Float64.(ref["x_range"]); rtol = 1e-12) &&
          isapprox([ylo, yhi], Float64.(ref["y_range"]); rtol = 1e-12, atol = 1e-15),
          "$(k): ranges differ from the reference package")
    if !endswith(k, "zero_line")
        @repro_result "$(tag)_x_at_y_min" s.x[argmin(s.y)]
        @repro_result "$(tag)_x_at_y_max" s.x[argmax(s.y)]
    end
end
for (n, lim) in sort!(collect(LIM); by = first)
    @repro_result "axes_$(n)_xlim" collect(lim.x)
    @repro_result "axes_$(n)_ylim" collect(lim.y)
end

# Structure the caption and A:63 describe, measured on the digitised series (no physics):
fa, da = S["a_dTdr/fiducial"], S["a_dTdr/double_res"]
cb = S["c_cumulative/fiducial"]
@repro_result "c_outer_r" cb.x[end]
@repro_result "c_outer_value" cb.y[end]
printed = 7.5e-4                                           # app_modal_convergence.tex:63 "~7.5x10^-4" (magnitude)
@repro_result "c_outer_magnitude_rel_to_printed" abs(abs(cb.y[end]) - printed) / printed
bb = S["b_delta/fiducial"]
@repro_result "b_drawn_integral_to_outer" sum((bb.x[k + 1] - bb.x[k]) * (bb.y[k + 1] + bb.y[k]) / 2
                                              for k in 1:length(bb.x) - 1)
check(issorted(fa.x) && issorted(da.x) && issorted(bb.x) && issorted(cb.x), "a drawn path is not monotone in r")

png = joinpath(@__DIR__, "expected", ID * ".png")
art = draw(png)
check(art.ok, art.message)
@repro_result "png_within_budget" art.ok
@repro_result "figure_checks_pass" isempty(violations)

println("\n[$(ID)] SIM reference only: ", REF["status"])
@printf("[%s] (a) trough %.4f at r = %.4f (fiducial), %.4f at r = %.4f (double)\n",
        ID, minimum(fa.y), fa.x[argmin(fa.y)], minimum(da.y), da.x[argmin(da.y)])
@printf("[%s] (b) spikes %+.4f at r = %.4f and %+.4f at r = %.4f\n", ID, maximum(bb.y), bb.x[argmax(bb.y)],
        minimum(bb.y), bb.x[argmin(bb.y)])
@printf("[%s] (c) hump %+.3e at r = %.4f; outer c(%.4f) = %+.4e vs printed ~7.5e-4 (magnitude)\n", ID,
        maximum(cb.y), cb.x[argmax(cb.y)], cb.x[end], cb.y[end])
println("[$(ID)] figure checks: ", isempty(violations) ? "PASS" : "FAIL")
foreach(v -> println("    ", v), violations)

code = ReproHarness.conclude!(opts; paper = PAPER, id = ID, script = @__FILE__, rtol = 1e-12,
                              tolerances = Dict("$(replace(k, "/" => "_"))_n" => (rtol = 0.0, atol = 0.0)
                                                for k in KEYS))
exit(isempty(violations) ? code : max(code, 1))

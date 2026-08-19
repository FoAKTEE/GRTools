# reproduce/arXiv-2608.17003/fig11.jl — EGD::fig11-convergence.
#
# arXiv:2608.17003 Fig. 11 (`fig:convergence`, appendix.tex:68-73): convergence of the PER-MODE
# torque T_m (Eq. 74 evaluated for one azimuthal mode) with the App.-C grid resolution.  Top panel:
# the absolute error |T_m − T_{m,ref}|; bottom panel: the relative error |T_m − T_{m,ref}|/|T_{m,ref}|,
# both versus m ∈ [1, 200] for the resolution triples (n_L = n_R, n_C) = (30, 3000), (60, 6000),
# (120, 12000) against the reference triple (n_L, n_R, n_C) = (240, 240, 24000) of A:140.
# Disc/EMRI parameters: a• = 25M, χ = 0, 𝔭 = 𝔮 = 0, h0 = 0.05.
#
# UNITS (MEASURED, not assumed).  The published absolute panel carries no unit label.  A pixel read
# of the reference render (evidence `fig11_pixel_compare.txt`) settles it: the published RELATIVE
# errors (bottom panel, unit-free) agree with ours to 10-20 % at m >= 100 for the two coarse
# resolutions, while the published ABSOLUTE errors are smaller than ours by a constant factor
# 18-21 = 1/h0 across every series and every m.  Since rel = abs/|T_m,ref|, that factor is the unit
# of T_m: the paper's per-mode torque axis is F_{J,0} (76a), not the F_{J,0}h0 used on its TOTAL
# torque axes (Figs. 6, 7, 10, 13).  T_m is therefore quoted here in F_{J,0}; the total torque is
# additionally registered in F_{J,0}h0 so it stays comparable with the other figures.
#
# ZERO physics here: the per-mode torques come from `AccretionDisc.SpiralWaves.torque_survey`
# (EGD::sw-survey-parallel), which solves each mode on its own App.-C grid and reduces it with the
# AMF/torque diagnostics.  The grid triple is passed through `grid_kwargs`; the DOMAIN of a mode is
# fixed by (C1)-(C3) and does not depend on (n_L, n_R, n_C), so the four runs differ only in
# resolution — which is what makes the difference a discretisation error.
# Style: `_style/emri_style.jl` + section 11 of `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md`.

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

const BLUE   = RGB(0x0f / 255, 0x4c / 255, 0x81 / 255)   # (30, 3000)
const ORANGE = RGB(0xe2 / 255, 0x58 / 255, 0x3e / 255)   # (60, 6000)
const GREEN  = RGB(0x00 / 255, 0x94 / 255, 0x73 / 255)   # (120, 12000)
const SERIES = ((BLUE, :star5), (ORANGE, :rect), (GREEN, :circle))

"Decade limits containing `vs` with one decade of headroom, and ticks every `step` decades."
function _decades(vs, step)
    lo = floor(Int, log10(minimum(x -> x > 0 ? x : Inf, vs))) - 1
    hi = ceil(Int, log10(maximum(vs))) + 1
    (lims = (10.0^lo, 10.0^hi), ticks = [10.0^e for e in (lo + 1):step:(hi - 1)])
end

"""
    draw(ms, abs_errs, rel_errs, labels) -> path

Render the two stacked log-y panels of Fig. 11 and save the PNG beside the golden.
"""
function draw(ms, abs_errs, rel_errs, labels)
    common = (xlims = (0, 205), xticks = 0:50:200, minorticks = 5, tick_direction = :in,
              yscale = :log10, grid = true, gridalpha = 0.35, gridlinewidth = 0.6,
              foreground_color_grid = RGB(0.75, 0.75, 0.75), widen = false,
              guidefontsize = 10, tickfontsize = 8)
    dtop = _decades(vcat(abs_errs...), 2)          # top panel: ticks every 2 decades (manual s.11)
    dbot = _decades(vcat(rel_errs...), 1)          # bottom panel: one tick per decade
    top = plot(; common..., ylims = dtop.lims, yticks = dtop.ticks,
               ylabel = "\$|T_m - T_{m,\\mathrm{ref}}|\$", xformatter = _ -> "",
               legend = :topright, legendfontsize = 6, legend_columns = 1,
               foreground_color_legend = RGB(0.8, 0.8, 0.8), background_color_legend = :white)
    for (i, (c, mk)) in enumerate(SERIES)
        scatter!(top, ms, abs_errs[i]; color = c, marker = (mk, 2.2), markerstrokewidth = 0,
                 label = labels[i])
    end
    bot = plot(; common..., ylims = dbot.lims, yticks = dbot.ticks, xlabel = "\$m\$",
               ylabel = "\$|T_m-T_{m,\\mathrm{ref}}|/|T_{m,\\mathrm{ref}}|\$", legend = false)
    for (i, (c, mk)) in enumerate(SERIES)
        scatter!(bot, ms, rel_errs[i]; color = c, marker = (mk, 2.2), markerstrokewidth = 0, label = "")
    end
    plt = plot(top, bot; layout = grid(2, 1), size = (320, 359), dpi = 160,
               left_margin = 3Plots.mm, bottom_margin = 1Plots.mm)
    png_path = joinpath(@__DIR__, "expected", "fig11.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig11] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig11.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

const ABS_KEYS = ("abs_err_30_3000", "abs_err_60_6000", "abs_err_120_12000")
const REL_KEYS = ("rel_err_30_3000", "rel_err_60_6000", "rel_err_120_12000")
const LABELS = ("\$n_L=n_R=30,\\ n_C=3000\$", "\$n_L=n_R=60,\\ n_C=6000\$",
                "\$n_L=n_R=120,\\ n_C=12000\$")

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig11")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["m"], [v[k] for k in ABS_KEYS], [v[k] for k in REL_KEYS], LABELS)
    println("[fig11] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves

for f in (:torque_survey, :mode_torques, :torque_units)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: Fig. 11 needs the threaded m-survey " *
        "(EGD::sw-survey-parallel) over the mode solver and diagnostics.")
end

# ───────────────────────────── parameters (appendix.tex:68-73, A:140) ────────────────────────────

const CHI   = 0.0
const AEMRI = 25.0
const H0    = 0.05
const PDISK = 0
const QDISK = 0
const BSOFT = 0.6
const EPS_P = 1e-5

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
const MMAX = _envint("GRTOOLS_FIG11_MMAX", 200)          # the published m range
const TRIPLES = ((30, 3000), (60, 6000), (120, 12000))
const REF_TRIPLE = (240, 24000)                          # (n_L = n_R, n_C) — A:140

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

bg = SW.DiscBackground(; χ = CHI, aEMRI = AEMRI, h0 = H0, 𝔭 = PDISK, 𝔮 = QDISK, b = BSOFT, εp = EPS_P)
FJ0, _ = SW.torque_units(bg)
unit = FJ0                                               # per-mode torque axis: F_{J,0} (see the header)

ms = collect(1:MMAX)
println("\n[fig11] a• = ", AEMRI, " M  χ = ", CHI, "  h0 = ", H0, "  (𝔭,𝔮) = (", PDISK, ",", QDISK, ")")
println("[fig11] m = 1:", MMAX, "   triples ", TRIPLES, "   reference (", REF_TRIPLE[1], ", ",
        REF_TRIPLE[1], ", ", REF_TRIPLE[2], ")")

run_triple(nLR, nC; verbose = false) = begin
    t0 = time()
    s = SW.torque_survey(bg, ms; grid_kwargs = (nL = nLR, nR = nLR, nC = nC), verbose = verbose)
    Tm = SW.mode_torques(s) ./ unit
    dt = time() - t0
    @printf("[fig11] (n_L=n_R=%3d, n_C=%5d)  %6.1f s   %10d nodes   ΣT_m = %.10f F_{J,0} = %.10f F_{J,0}h0\n",
            nLR, nC, dt, sum(mt.nodes for mt in s.modes), sum(Tm), sum(Tm) / H0)
    flush(stdout)
    (Tm = Tm, seconds = dt, residual = maximum(mt.residual for mt in s.modes),
     nodes = sum(mt.nodes for mt in s.modes))
end

ref = run_triple(REF_TRIPLE...; verbose = true)
runs = [run_triple(t...) for t in TRIPLES]

abs_errs = [abs.(r.Tm .- ref.Tm) for r in runs]
rel_errs = [abs.(r.Tm .- ref.Tm) ./ abs.(ref.Tm) for r in runs]

@repro_result "m"            ms
@repro_result "T_m_ref"      ref.Tm
for (k, v) in zip(ABS_KEYS, abs_errs); @repro_result k v; end
for (k, v) in zip(REL_KEYS, rel_errs); @repro_result k v; end
@repro_result "a_emri"       AEMRI
@repro_result "chi"          CHI
@repro_result "h0"           H0
@repro_result "p_disk"       PDISK
@repro_result "q_disk"       QDISK
@repro_result "b_soft"       BSOFT
@repro_result "eps_p"        EPS_P
@repro_result "m_max"        MMAX
@repro_result "torque_unit_F_J0"    unit
@repro_result "T_total_ref_FJ0"     sum(ref.Tm)
@repro_result "T_total_ref_FJ0h0"   sum(ref.Tm) / H0
@repro_result "solver_residual_max" max(ref.residual, maximum(r.residual for r in runs))
@repro_result "nodes_reference"     ref.nodes

# ── convergence summary: the median relative error of each triple, and the ratios between them ───
_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? float(s[(n+1)÷2]) : (s[n÷2] + s[n÷2+1]) / 2)
med = [_median(r) for r in rel_errs]
@repro_result "rel_err_median_30_3000"    med[1]
@repro_result "rel_err_median_60_6000"    med[2]
@repro_result "rel_err_median_120_12000"  med[3]
@repro_result "rel_err_ratio_30_over_60"  med[1] / med[2]
@repro_result "rel_err_ratio_60_over_120" med[2] / med[3]
@repro_result "total_rel_err_30_3000"     abs(sum(runs[1].Tm) - sum(ref.Tm)) / abs(sum(ref.Tm))
@repro_result "total_rel_err_60_6000"     abs(sum(runs[2].Tm) - sum(ref.Tm)) / abs(sum(ref.Tm))
@repro_result "total_rel_err_120_12000"   abs(sum(runs[3].Tm) - sum(ref.Tm)) / abs(sum(ref.Tm))

println("\n[fig11] T_m is quoted in F_{J,0} (76a) — the unit measured from the published render, see the file header")
println("[fig11] convergence summary (median relative error over m = 1:", MMAX, ")")
for (i, (nLR, nC)) in enumerate(TRIPLES)
    @printf("  (n_L=n_R=%3d, n_C=%5d)  median rel err = %.3e   max = %.3e   total-torque rel err = %.3e\n",
            nLR, nC, med[i], maximum(rel_errs[i]),
            abs(sum(runs[i].Tm) - sum(ref.Tm)) / abs(sum(ref.Tm)))
end
@printf("  ratios: (30,3000)/(60,6000) = %.2f   (60,6000)/(120,12000) = %.2f\n",
        med[1] / med[2], med[2] / med[3])

# figure claim (A:68-73): the error must DECREASE monotonically with resolution — this is the whole
# statement of the convergence study.  A violation exits 1.
violations = String[]
med[1] > med[2] > med[3] ||
    push!(violations, "the median relative error is not monotone in resolution: $(med)")
println("\n[fig11] figure check: median relative error decreases with resolution   ",
        med[1] > med[2] > med[3] ? "PASS" : "FAIL")
@repro_result "figure_checks_pass" isempty(violations)

check_failed = false
if !isempty(violations)
    println("\n[fig11] FINDING — figure claim not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw(ms, abs_errs, rel_errs, LABELS)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig11", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m" => (rtol = 0.0, atol = 0.0),
                                                "m_max" => (rtol = 0.0, atol = 0.0)))
exit(check_failed ? max(code, 1) : code)

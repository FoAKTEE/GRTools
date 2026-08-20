# reproduce/arXiv-2608.17003/fig07.jl — EGD::fig07-torque-ap-h0.
#
# arXiv:2608.17003 Fig. 7 (`fig:torque_ap_h0`, main.tex:944-949): the TOTAL torque (74) in units of
# F_{J,0}h0 (76a) against a•/M ∈ [10, 1000] (log x, the same 40 log-spaced points as Fig. 6) for the
# five disc thicknesses h0 ∈ {0.05, 0.075, 0.1, 0.125, 0.15} at χ = 0, in two stacked panels: top
# (𝔭,𝔮) = (1.5,1), bottom (0,0).
#
# Caption claims: "Thicker disks are impacted strongly by the inner regions of the disk leading to some
# oscillatory behavior. Observe that the torque monotonically increases until h0 = 0.125 and then
# decreases for very thick disks h0 = 0.15" (M:947). Together with M:925 ("no such torque reversal, at
# least for disk thicknesses 0.05 ≤ h0 ≤ 0.15") these are THE two printed claims the paper attaches to
# this h0 range, and they are the two golden keys the printed-number gate reads from this driver:
#   `total_torque_all_positive`                      (claim M925-no-reversal)
#   `torque_peak_h0`, `torque_monotone_up_to_peak`   (claim M947-h0-turnover)
# The h0 ordering is read in the PLATEAU of the (𝔭,𝔮) = (0,0) panel (a• ≥ 100M, where the published
# markers are cleanly ordered and the hole is irrelevant); the value at a single a• would be dominated
# by the oscillations the caption itself describes. That choice is recorded in the golden
# (`plateau_a_min`) and printed.
#
# ZERO physics here: `AccretionDisc.SpiralWaves.torque_survey` does the per-mode solve, reduction and
# sum. Style: `_style/emri_style.jl` + section 7 of the VLM manual
# `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo).

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

_rgb(hex) = RGB(((hex >> 16) & 0xff) / 255, ((hex >> 8) & 0xff) / 255, (hex & 0xff) / 255)

# manual section 0, palette 2: plasma sampled at linspace(0.1, 0.9, 5), hollow, with the marker sequence
const THICK = (
    (h0 = 0.050, color = _rgb(0x41049d), marker = :rect,      label = "\$h_0=0.050\$"),
    (h0 = 0.075, color = _rgb(0x8f0da4), marker = :star5,     label = "\$h_0=0.075\$"),
    (h0 = 0.100, color = _rgb(0xcc4778), marker = :circle,    label = "\$h_0=0.100\$"),
    (h0 = 0.125, color = _rgb(0xf2844b), marker = :rtriangle, label = "\$h_0=0.125\$"),
    (h0 = 0.150, color = _rgb(0xfcce25), marker = :hexagon,   label = "\$h_0=0.150\$"),
)
# `ylim` = the y-limits MEASURED off the reference PDF at 200 dpi (major-tick rows 45.5/107.5 = 2.00/1.75
# in the top panel spanning rows 1-231, and 273.5/371.0 = 1.2/1.0 in the bottom panel spanning rows
# 234-464), NOT the manual's eyeball estimate.
const PANELS = ((𝔭 = 1.5, 𝔮 = 1.0, stamp = "\$\\mathfrak{p}=1.5,\\ \\mathfrak{q}=1.0\$", ylim = (1.252, 2.179),
                 yticks = [1.25, 1.50, 1.75, 2.00]),
                (𝔭 = 0.0, 𝔮 = 0.0, stamp = "\$\\mathfrak{p}=0.0,\\ \\mathfrak{q}=0.0\$", ylim = (0.809, 1.281),
                 yticks = [1.0, 1.2]))

# The paper's markers are `mfc='none'`. In Plots/GR a fully transparent face is drawn BLACK and
# `markeralpha = 0` hides the stroke as well, so the hollow look is produced with a WHITE face on the
# white page (measured: identical except where two markers overlap, where ours occludes instead of
# showing through).
const HOLLOW = :white

"Published y-limits, widened when our data would otherwise be clipped (never silently)."
function _ylim(pan, series)
    lo = minimum(minimum, series); hi = maximum(maximum, series)
    pad = 0.04 * (pan.ylim[2] - pan.ylim[1])
    out = (min(pan.ylim[1], lo - pad), max(pan.ylim[2], hi + pad))
    out == pan.ylim || @warn "fig07: data outside the published y-range $(pan.ylim) → widened to $(out)"
    out
end

"""
    draw(a_over_M, T) -> path

Render Fig. 7 from `T[p][k]` = the torque series of panel `p`, thickness `k`, and save the PNG beside
the golden (`--plot-only` redraws from the golden without recomputing).
"""
function draw(a_over_M, T)
    plts = map(enumerate(PANELS)) do (p, pan)
        yl = _ylim(pan, T[p])
        last_panel = p == length(PANELS)
        plt = plot(xscale = :log10, xlims = (9.5, 1080),
                   xticks = ([10.0, 100.0, 1000.0], last_panel ?
                             ["\$10^{1}\$", "\$10^{2}\$", "\$10^{3}\$"] : ["", "", ""]),
                   ylims = yl, yticks = pan.yticks, minorticks = 5,
                   tick_direction = :in, widen = false, legend = p == 1 ? (0.72, 0.62) : false,
                   legendfontsize = 7, foreground_color_legend = RGB(0.8, 0.8, 0.8),
                   background_color_legend = :white, guidefontsize = 10, tickfontsize = 8,
                   ylabel = "\$T\\ (F_{J_0}h_0)\$",
                   xlabel = last_panel ? "\$a_\\bullet/M\$" : "")
        for (k, th) in enumerate(THICK)
            scatter!(plt, a_over_M, T[p][k]; marker = th.marker, markersize = 4.2,
                     markercolor = HOLLOW, markerstrokecolor = th.color, markerstrokewidth = 1.1,
                     label = th.label)
        end
        annotate!(plt, 10.6, yl[1] + 0.93 * (yl[2] - yl[1]),
                  text(pan.stamp, :left, 8, "Computer Modern"))
        plt
    end
    plt = plot(plts...; layout = grid(2, 1), size = (330, 300), dpi = 160,
               bottom_margin = 0Plots.mm, left_margin = 1Plots.mm, top_margin = 0Plots.mm)
    png_path = joinpath(@__DIR__, "expected", "fig07.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig07] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig07.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig07")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["a_over_M"], [[v["T_$(p)_$(k)"] for k in 1:length(THICK)] for p in 1:length(PANELS)])
    println("[fig07] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves
for f in (:torque_survey, :torque_units, :mode_torques)
    isdefined(SW, f) || ReproHarness.hole("AccretionDisc.SpiralWaves.$(f) is not available (EGD::sw-survey-parallel).")
end

# ───────────────────────────── parameters (main.tex:944-949, brief section 4 row 7) ──────────────

const CHI   = 0.0
const BSOFT = 0.6                  # nominal softening (M:884)
const EPS_P = 1e-5                 # nominal pole displacement (M:735)

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

# m ceiling PER THICKNESS: ⌈K/h0⌉ with K = 4 (MEASURED: at h0 = 0.05 that is m ≤ 80 and truncates T by
# ≈ 0.5 %; m₉₅h₀ ≃ 2.3 there, inside the published 2.0-3.5 band of M:1025, so 4/h0 sits well beyond
# m₉₅ for every thickness). The tail actually reached is measured per configuration and registered.
const KMODE = _envint("GRTOOLS_FIG07_K", 4)
mmax_of(h0) = ceil(Int, KMODE / h0)
const GRIDKW = (nL = _envint("GRTOOLS_FIG07_NL", 30), nR = _envint("GRTOOLS_FIG07_NL", 30),
                nC = _envint("GRTOOLS_FIG07_NC", 3000))
const NA = _envint("GRTOOLS_FIG07_NA", 40)
const AS = 10.0 .^ range(1.0, 3.0; length = NA)
const PLATEAU_A = 100.0            # the a• ≥ 100M plateau used to read the h0 ordering (M:947)

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

"""
    survey_panel(𝔭, 𝔮) -> (T, T_L, tail, res)

The five thickness series of one panel: T(a•)/(F_{J,0}h0) for every a• of `AS` and every h0 of `THICK`.
"""
function survey_panel(𝔭, 𝔮)
    T = [Float64[] for _ in THICK]; TL = [Float64[] for _ in THICK]
    tail = 0.0; res = 0.0
    for (k, th) in enumerate(THICK), a in AS
        bg = SW.DiscBackground(; χ = CHI, aEMRI = a, h0 = th.h0, 𝔭 = 𝔭, 𝔮 = 𝔮, b = BSOFT, εp = EPS_P)
        mmax = mmax_of(th.h0)
        t0 = time()
        sv = SW.torque_survey(bg, 1:mmax; threads = true, grid_kwargs = GRIDKW)
        FJ0, _ = SW.torque_units(bg)
        Tm = SW.mode_torques(sv)
        ntail = max(1, mmax ÷ 20)
        tl = sum(Tm[end-ntail+1:end]) / sv.T_total
        push!(T[k], sv.T_total / (FJ0 * th.h0)); push!(TL[k], sv.T_L_total / (FJ0 * th.h0))
        tail = max(tail, abs(tl)); res = max(res, maximum(mt.residual for mt in sv.modes))
        @printf("  (𝔭,𝔮)=(%+.1f,%.1f) h0=%.3f a•=%8.2f M (m≤%3d): T = %+.6f  T_L = %+.6f  (tail %.1e, %.0f s)\n",
                𝔭, 𝔮, th.h0, a, mmax, T[k][end], TL[k][end], tl, time() - t0)
        flush(stdout)
    end
    (T, TL, tail, res)
end

println("\n[fig07] χ = ", CHI, "   b = ", BSOFT, "   m ≤ ", KMODE, "/h0 = ",
        [mmax_of(th.h0) for th in THICK], "   grid = ", GRIDKW,
        "   a•/M = ", NA, " log-spaced points in [", first(AS), ", ", last(AS), "]")
println("[fig07] compute budget: ", SW.compute_budget().note)

t_start = time()
panels = [survey_panel(pan.𝔭, pan.𝔮) for pan in PANELS]
T  = [p[1] for p in panels]
TL = [p[2] for p in panels]
tailmax = maximum(p[3] for p in panels)
resmax  = maximum(p[4] for p in panels)
@printf("[fig07] %d configurations in %.0f s (%.1f min)\n", length(PANELS) * length(THICK) * NA,
        time() - t_start, (time() - t_start) / 60)

@repro_result "a_over_M" collect(AS)
@repro_result "chi"      CHI
@repro_result "b_soft"   BSOFT
@repro_result "eps_p"    EPS_P
@repro_result "k_mode"   KMODE
@repro_result "n_L"      GRIDKW.nL
@repro_result "n_C"      GRIDKW.nC
@repro_result "plateau_a_min" PLATEAU_A
for (k, th) in enumerate(THICK)
    @repro_result "h0_$(k)"    th.h0
    @repro_result "m_max_$(k)" mmax_of(th.h0)
end
for (p, pan) in enumerate(PANELS)
    @repro_result "p_$(p)" pan.𝔭
    @repro_result "q_$(p)" pan.𝔮
    for k in eachindex(THICK)
        @repro_result "T_$(p)_$(k)"   T[p][k]
        @repro_result "T_L_$(p)_$(k)" TL[p][k]
    end
end
@repro_result "mode_sum_tail_fraction_max" tailmax
@repro_result "solver_residual_max"        resmax

# ───────────────────────────── the paper's printed checks ────────────────────────────────────────
# M:925 (no torque reversal for 0.05 ≤ h0 ≤ 0.15) and M:947 (T increases with h0 up to 0.125, then
# decreases at 0.15). The gate `check_printed_numbers.jl` reads the three keys registered here.

allpos = all(all(>(0), T[p][k]) for p in eachindex(PANELS), k in eachindex(THICK))
@repro_result "total_torque_all_positive" allpos

const IPLATEAU = findall(>=(PLATEAU_A), AS)
plateau_T = [[sum(T[p][k][IPLATEAU]) / length(IPLATEAU) for k in eachindex(THICK)] for p in eachindex(PANELS)]
const IPQ00 = findfirst(pan -> pan.𝔭 == 0 && pan.𝔮 == 0, PANELS)          # the (0,0) panel of M:947
kpeak = argmax(plateau_T[IPQ00])
peak_h0 = THICK[kpeak].h0
monotone = all(plateau_T[IPQ00][k] < plateau_T[IPQ00][k+1] for k in 1:(kpeak - 1))
@repro_result "torque_peak_h0"              peak_h0
@repro_result "torque_monotone_up_to_peak"  monotone
for (p, _) in enumerate(PANELS)
    @repro_result "plateau_T_$(p)" plateau_T[p]
end

violations = String[]
allpos || push!(violations, "M:925 torque reversal: T ≤ 0 somewhere over h0 ∈ [0.05, 0.15], a• ∈ [10, 1000]M")
peak_h0 == 0.125 ||
    push!(violations, "M:947 turnover: the plateau torque of the (0,0) panel peaks at h0 = $(peak_h0), not 0.125")
monotone ||
    push!(violations, "M:947 monotonicity: the plateau torque is not increasing in h0 below the peak: $(plateau_T[IPQ00])")

println("\n[fig07] figure checks")
@printf("  M:925 no torque reversal (T > 0 everywhere)                     %s\n", allpos ? "PASS" : "FAIL")
@printf("  M:947 plateau (a• ≥ %.0f M) torque of the (0,0) panel vs h0: %s\n", PLATEAU_A,
        join([@sprintf("%.3f→%.4f", THICK[k].h0, plateau_T[IPQ00][k]) for k in eachindex(THICK)], "  "))
@printf("        peak at h0 = %.3f (paper 0.125)   monotone below the peak: %s   %s\n", peak_h0, monotone,
        (peak_h0 == 0.125 && monotone) ? "PASS" : "FAIL")
@printf("  (1.5,1) plateau vs h0: %s\n",
        join([@sprintf("%.3f→%.4f", THICK[k].h0, plateau_T[1][k]) for k in eachindex(THICK)], "  "))
@printf("  mode-sum tail (last 5%% of modes)  %.2e\n", tailmax)
@printf("  worst solver residual             %.2e\n", resmax)

check_failed = false
if !isempty(violations)
    println("\n[fig07] FINDING — printed claims not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw(collect(AS), T)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig07", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("k_mode" => (rtol = 0.0, atol = 0.0),
                                                "n_L" => (rtol = 0.0, atol = 0.0),
                                                "n_C" => (rtol = 0.0, atol = 0.0)))
exit(check_failed ? max(code, 1) : code)

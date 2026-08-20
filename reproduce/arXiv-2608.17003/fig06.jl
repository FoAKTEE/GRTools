# reproduce/arXiv-2608.17003/fig06.jl — EGD::fig06-torque-chi-ap.
#
# arXiv:2608.17003 Fig. 6 (`fig:torque_chi_ap`, main.tex:934-939): the TOTAL torque (74) in units of
# F_{J,0}h0 (76a) against the EMRI semi-major axis a•/M ∈ [10, 1000] (log x, 40 log-spaced points —
# counted on the reference PDF: the marker pitch is 13.565 px on a 264.5 px/decade axis, i.e.
# Δlog10 = 2/39), for the three spins χ ∈ {0.9, 0, −0.9}, in two stacked panels: top (𝔭,𝔮) = (1.5,1),
# bottom (0,0); h0 = 0.05 throughout.
#
# Caption claim (M:936-939): "Relativistic effects are important as we approach a• ∼ O(50M); the
# retrograde orbits (χ < 0) face stronger relativistic effects since the ISCO is closer to the EMRI
# orbit compared to the prograde (χ > 0) and the non-spinning (χ = 0) cases." Plus M:925, the
# statement this figure's parameter range backs: there is NO torque reversal (T > 0 everywhere).
# Both are asserted below (a violation exits 1); `total_torque_all_positive` is also the golden key the
# printed-number gate reads for M:925 (it reads fig07's copy — this driver registers its own).
#
# ZERO physics here: `AccretionDisc.SpiralWaves.torque_survey` does the per-mode solve, reduction and
# sum. Style: `_style/emri_style.jl` + section 6 of the VLM manual
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

# manual section 0, palette 4 (Fig. 6 has its own χ triple) — hollow markers
const SPINS = (
    (χ =  0.9, color = _rgb(0xe70234), marker = :rect,   label = "\$\\chi=0.9\$"),
    (χ =  0.0, color = _rgb(0x580772), marker = :star5,  label = "\$\\chi=0.0\$"),
    (χ = -0.9, color = _rgb(0x2584da), marker = :circle, label = "\$\\chi=-0.9\$"),
)
# `ylim` = the y-limits MEASURED off the reference PDF at 200 dpi (frame rows against the major-tick
# rows: top panel 2.0/1.5/1.0 at rows 24.5/106.5/189.5 in a panel spanning rows 1-231; bottom panel
# 0.75/0.50 at 424.5/360.5 in rows 234-464), NOT the manual's eyeball estimate — matplotlib's default
# 5 % margin puts the frame slightly outside the data, and reproducing it is what keeps the panel
# stamp clear of the markers.
const PANELS = ((𝔭 = 1.5, 𝔮 = 1.0, stamp = "\$\\mathfrak{p}=1.5,\\ \\mathfrak{q}=1.0\$", ylim = (0.747, 2.143),
                 yticks = [1.0, 1.5, 2.0]),
                (𝔭 = 0.0, 𝔮 = 0.0, stamp = "\$\\mathfrak{p}=0.0,\\ \\mathfrak{q}=0.0\$", ylim = (0.346, 1.244),
                 yticks = [0.50, 0.75, 1.00]))

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
    out == pan.ylim || @warn "fig06: data outside the published y-range $(pan.ylim) → widened to $(out)"
    out
end

"""
    draw(a_over_M, T) -> path

Render Fig. 6 from `T[p][s]` = the torque series of panel `p`, spin `s`, and save the PNG beside the
golden (`--plot-only` redraws from the golden without recomputing).
"""
function draw(a_over_M, T)
    plts = map(enumerate(PANELS)) do (p, pan)
        yl = _ylim(pan, T[p])
        last_panel = p == length(PANELS)
        plt = plot(xscale = :log10, xlims = (9.5, 1080),
                   xticks = ([10.0, 100.0, 1000.0], last_panel ?
                             ["\$10^{1}\$", "\$10^{2}\$", "\$10^{3}\$"] : ["", "", ""]),
                   ylims = yl, yticks = pan.yticks, minorticks = 5,
                   tick_direction = :in, widen = false, legend = p == 1 ? (0.80, 0.55) : false,
                   legendfontsize = 7, foreground_color_legend = RGB(0.8, 0.8, 0.8),
                   background_color_legend = :white, guidefontsize = 10, tickfontsize = 8,
                   ylabel = "\$T\\ (F_{J_0}h_0)\$",
                   xlabel = last_panel ? "\$a_\\bullet/M\$" : "")
        for (s, sp) in enumerate(SPINS)
            scatter!(plt, a_over_M, T[p][s]; marker = sp.marker, markersize = 4.4,
                     markercolor = HOLLOW, markerstrokecolor = sp.color, markerstrokewidth = 1.1,
                     label = sp.label)
        end
        annotate!(plt, 10.6, yl[1] + 0.93 * (yl[2] - yl[1]),
                  text(pan.stamp, :left, 8, "Computer Modern"))
        plt
    end
    plt = plot(plts...; layout = grid(2, 1), size = (330, 300), dpi = 160,
               bottom_margin = 0Plots.mm, left_margin = 1Plots.mm, top_margin = 0Plots.mm)
    png_path = joinpath(@__DIR__, "expected", "fig06.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig06] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig06.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig06")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["a_over_M"], [[v["T_$(p)_$(s)"] for s in 1:length(SPINS)] for p in 1:length(PANELS)])
    println("[fig06] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves
for f in (:torque_survey, :torque_units, :mode_torques)
    isdefined(SW, f) || ReproHarness.hole("AccretionDisc.SpiralWaves.$(f) is not available (EGD::sw-survey-parallel).")
end

# ───────────────────────────── parameters (main.tex:934-939, brief section 4 row 6) ──────────────

const H0    = 0.05
const BSOFT = 0.6                  # nominal softening (M:884)
const EPS_P = 1e-5                 # nominal pole displacement (M:735)

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

# m ceiling: ⌈4/h0⌉ = 80 (MEASURED truncation ≈ 0.5 % of T at h0 = 0.05 — see fig05.jl and the node
# report); the tail actually reached is measured per configuration and registered.
const MMAX = _envint("GRTOOLS_FIG06_MMAX", ceil(Int, 4 / H0))
const GRIDKW = (nL = _envint("GRTOOLS_FIG06_NL", 30), nR = _envint("GRTOOLS_FIG06_NL", 30),
                nC = _envint("GRTOOLS_FIG06_NC", 3000))
const NA = _envint("GRTOOLS_FIG06_NA", 40)
const AS = 10.0 .^ range(1.0, 3.0; length = NA)

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

"""
    survey_panel(𝔭, 𝔮) -> (T, T_L, tail, res)

The three spin series of one panel: T(a•)/(F_{J,0}h0) for every a• of `AS` and every χ of `SPINS`.
"""
function survey_panel(𝔭, 𝔮)
    T = [Float64[] for _ in SPINS]; TL = [Float64[] for _ in SPINS]
    tail = 0.0; res = 0.0
    for (s, sp) in enumerate(SPINS), a in AS
        bg = SW.DiscBackground(; χ = sp.χ, aEMRI = a, h0 = H0, 𝔭 = 𝔭, 𝔮 = 𝔮, b = BSOFT, εp = EPS_P)
        t0 = time()
        sv = SW.torque_survey(bg, 1:MMAX; threads = true, grid_kwargs = GRIDKW)
        FJ0, _ = SW.torque_units(bg)
        Tm = SW.mode_torques(sv)
        ntail = max(1, MMAX ÷ 20)
        tl = sum(Tm[end-ntail+1:end]) / sv.T_total
        push!(T[s], sv.T_total / (FJ0 * H0)); push!(TL[s], sv.T_L_total / (FJ0 * H0))
        tail = max(tail, abs(tl)); res = max(res, maximum(mt.residual for mt in sv.modes))
        @printf("  (𝔭,𝔮)=(%+.1f,%.1f) χ=%+.1f a•=%8.2f M : T = %+.6f  T_L = %+.6f  (tail %.1e, %.0f s)\n",
                𝔭, 𝔮, sp.χ, a, T[s][end], TL[s][end], tl, time() - t0)
        flush(stdout)
    end
    (T, TL, tail, res)
end

println("\n[fig06] h0 = ", H0, "   b = ", BSOFT, "   m = 1:", MMAX, "   grid = ", GRIDKW,
        "   a•/M = ", NA, " log-spaced points in [", first(AS), ", ", last(AS), "]")
println("[fig06] compute budget: ", SW.compute_budget().note)

t_start = time()
panels = [survey_panel(pan.𝔭, pan.𝔮) for pan in PANELS]
T  = [p[1] for p in panels]
TL = [p[2] for p in panels]
tailmax = maximum(p[3] for p in panels)
resmax  = maximum(p[4] for p in panels)
@printf("[fig06] %d configurations in %.0f s (%.1f min)\n", length(PANELS) * length(SPINS) * NA,
        time() - t_start, (time() - t_start) / 60)

@repro_result "a_over_M" collect(AS)
@repro_result "h0"       H0
@repro_result "b_soft"   BSOFT
@repro_result "eps_p"    EPS_P
@repro_result "m_max"    MMAX
@repro_result "n_L"      GRIDKW.nL
@repro_result "n_C"      GRIDKW.nC
for (s, sp) in enumerate(SPINS)
    @repro_result "chi_$(s)" sp.χ
end
for (p, pan) in enumerate(PANELS)
    @repro_result "p_$(p)" pan.𝔭
    @repro_result "q_$(p)" pan.𝔮
    for s in eachindex(SPINS)
        @repro_result "T_$(p)_$(s)"   T[p][s]
        @repro_result "T_L_$(p)_$(s)" TL[p][s]
    end
end
@repro_result "mode_sum_tail_fraction_max" tailmax
@repro_result "solver_residual_max"        resmax

# ───────────────────────────── the figure's own printed checks ───────────────────────────────────
# (i)  M:925 — no torque reversal anywhere in the surveyed range (this is also the fig07 golden key
#      the printed-number gate reads; registered here under the same name for this figure's own range);
# (ii) M:936-939 — the retrograde disc is the one relativity suppresses most: at the smallest a• the
#      ordering is T(χ=−0.9) < T(χ=0) < T(χ=0.9), and all three converge to a common plateau at large
#      a• (where the hole is irrelevant): the spread at a• = 1000M is ≲ 2 %.

allpos = all(all(>(0), T[p][s]) for p in eachindex(PANELS), s in eachindex(SPINS))
@repro_result "total_torque_all_positive" allpos

i_small = 1; i_big = NA
order_small = [T[p][3][i_small] < T[p][2][i_small] < T[p][1][i_small] for p in eachindex(PANELS)]
spread_big = [(maximum(T[p][s][i_big] for s in eachindex(SPINS)) -
               minimum(T[p][s][i_big] for s in eachindex(SPINS))) /
              maximum(T[p][s][i_big] for s in eachindex(SPINS)) for p in eachindex(PANELS)]
@repro_result "spin_ordering_at_a10"  all(order_small)
@repro_result "plateau_spread_a1000"  maximum(spread_big)
for (p, pan) in enumerate(PANELS)
    @repro_result "plateau_T_$(p)" T[p][2][i_big]         # the χ = 0 plateau value at a• = 1000M
end

violations = String[]
allpos || push!(violations, "M:925 torque reversal: T ≤ 0 somewhere in a• ∈ [10, 1000]M")
all(order_small) ||
    push!(violations, "M:936 spin ordering at a• = $(AS[i_small])M violated: T = $([[T[p][s][i_small] for s in eachindex(SPINS)] for p in eachindex(PANELS)])")
maximum(spread_big) <= 0.05 ||
    push!(violations, "the three spins do not converge at a• = 1000M: spread $(maximum(spread_big))")

println("\n[fig06] figure checks")
@printf("  M:925 no torque reversal (T > 0 everywhere)                    %s\n", allpos ? "PASS" : "FAIL")
@printf("  M:936 T(χ=−0.9) < T(χ=0) < T(χ=0.9) at a• = %.1f M            %s\n", AS[i_small],
        all(order_small) ? "PASS" : "FAIL")
@printf("  spin spread at a• = %.0f M                        %.3f  (≤ 0.05)  %s\n", AS[i_big],
        maximum(spread_big), maximum(spread_big) <= 0.05 ? "PASS" : "FAIL")
for (p, pan) in enumerate(PANELS)
    @printf("  (𝔭,𝔮) = (%+.1f,%.1f): plateau T(χ=0, a•=1000M) = %.4f  F_J0h0;  T range [%.4f, %.4f]\n",
            pan.𝔭, pan.𝔮, T[p][2][i_big], minimum(minimum, T[p]), maximum(maximum, T[p]))
end
@printf("  mode-sum tail (last 5%% of modes)  %.2e\n", tailmax)
@printf("  worst solver residual             %.2e\n", resmax)

check_failed = false
if !isempty(violations)
    println("\n[fig06] FINDING — figure claims not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw(collect(AS), T)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig06", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_L" => (rtol = 0.0, atol = 0.0),
                                                "n_C" => (rtol = 0.0, atol = 0.0)))
exit(check_failed ? max(code, 1) : code)

# reproduce/arXiv-2608.17003/fig05.jl — EGD::fig05-tgr-tnewt.
#
# arXiv:2608.17003 Fig. 5 (`fig:TGR_vs_TNewt`, main.tex:927-932): the ratio of the RELATIVISTIC total
# torque (74) to the NEWTONIAN one, as a function of the SMBH spin χ ∈ [−0.9, 0.9] (10 points), for the
# five disc profiles (𝔭,𝔮) ∈ {(−1.5,0), (−1.5,1), (0,0), (1.5,0), (1.5,1)} at h0 = 0.05, a• = 12M.
# Caption claim: "disks with negative (positive) surface density slope 𝔭 > 0 (𝔭 < 0) can decrease
# (increase) the relativistic torque compared to the Newtonian calculation" — i.e. the ratio is < 1 for
# every 𝔭 ≥ 0 series and > 1 for every 𝔭 = −1.5 series. That is the figure's assertable content and it
# is checked below (a violation exits 1).
#
# WHAT "T_Newt" IS (the definition this driver commits to, recorded because the paper does not spell it
# out).  T_Newt is the torque of the FULLY NEWTONIAN linear problem with the same disc parameters:
# `AccretionDisc.SpiralWaves` `model = :newtonian` (src/SpiralWaves/newtonian.jl, EGD::sw-newtonian-limit)
#   • background   Ω_K = r^{−3/2}, Ω•,N = a•^{−3/2}, Σ ∝ r^{−𝔭} (12a), c_s ∝ r^{−𝔮/2} (12c), Γ₁ = 1 + c_s²,
#                  rotation from the NEWTONIAN radial balance Ω² = 1/r³ − (𝔭+𝔮)c_s²/r², κ² = 4Ω² + 2rΩΩ′;
#   • master eqs   the printed Newtonian limit (50a)-(50g);
#   • forcing      the Plummer-softened Newtonian potential Φ = −m•/√(ε_soft² + r² + a•² − 2ra•cosΨ),
#                  ε_soft = a•h0b (59), cosine-projected by (60);
#   • AMF/torque   (72a) with √−g → r, u_φ → r²Ω, ε+p → Σ and the COVARIANT δû_φ = r × (50a)
#                  (ERRATUM E-M1), and the Goldreich–Tremaine torque density −4πmrΦ̂ Im δΣ̂ = the M/r → 0
#                  limit of (72b);
#   • domain       the same App.-C machinery with the fiducial (C1) inner edge r_ISCO(χ = 0) = 6M.
# It is therefore χ-INDEPENDENT (verified bitwise in the package tests): the entire χ-dependence of the
# plotted ratio comes from the numerator. That reading is FORCED by the published figure itself — the
# measured T_GR(χ=+0.9)/T_GR(χ=−0.9) at these parameters is 1.578 / 2.422 / 0.838 for (𝔭,𝔮) = (0,0) /
# (1.5,0) / (−1.5,0), while the published ratio-of-ratios is 1.617 / 2.67 / 0.843; a χ-dependent
# denominator could not reproduce that.
#
# ZERO physics here: the disc, the modes, the two models and the torque sums all come from
# `AccretionDisc.SpiralWaves`. Style: `_style/emri_style.jl` + section 5 of the VLM manual
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

# manual section 0, palette 1 (the "Pantone-like" categorical cycle) with the section-5 pairing
const SERIES = (
    (𝔭 = -1.5, 𝔮 = 0.0, color = _rgb(0x0f4c81), marker = :circle,    label = "\$\\mathfrak{p}=-1.5,\\mathfrak{q}=0.0\$"),
    (𝔭 = -1.5, 𝔮 = 1.0, color = _rgb(0xe2583e), marker = :star5,     label = "\$\\mathfrak{p}=-1.5,\\mathfrak{q}=1.0\$"),
    (𝔭 =  0.0, 𝔮 = 0.0, color = _rgb(0x009473), marker = :diamond,   label = "\$\\mathfrak{p}=0.0,\\mathfrak{q}=0.0\$"),
    (𝔭 =  1.5, 𝔮 = 0.0, color = _rgb(0xbf1932), marker = :dtriangle, label = "\$\\mathfrak{p}=1.5,\\mathfrak{q}=0.0\$"),
    (𝔭 =  1.5, 𝔮 = 1.0, color = _rgb(0x5f4b8b), marker = :rect,      label = "\$\\mathfrak{p}=1.5,\\mathfrak{q}=1.0\$"),
)

"""
    draw(chis, ratios) -> path

Render Fig. 5 from the plotted arrays (`ratios[k]` = the ratio series of `SERIES[k]`) and save the PNG
beside the golden. Split out so `--plot-only` can redraw from `expected/fig05.json` after a style change
without repeating the surveys.
"""
function draw(chis, ratios)
    plt = plot(size = (330, 240), dpi = 160, legend = :bottomright, legendfontsize = 7,
               guidefontsize = 11, tickfontsize = 9, foreground_color_legend = RGB(0.8, 0.8, 0.8),
               background_color_legend = :white,
               xlims = (-1.0, 1.0), ylims = (0.15, 1.6), xticks = -0.5:0.5:0.5,
               yticks = 0.25:0.25:1.50, minorticks = 5, tick_direction = :in, widen = false,
               xlabel = "\$\\chi\$", ylabel = "\$T_{\\mathrm{GR}}/T_{\\mathrm{Newt}}\$",
               bottom_margin = 1Plots.mm, left_margin = 1Plots.mm)
    for (k, s) in enumerate(SERIES)
        plot!(plt, chis, ratios[k]; seriestype = :scatter, marker = s.marker, markersize = 4.6,
              markercolor = s.color, markerstrokecolor = s.color, markerstrokewidth = 0.3,
              label = s.label)
    end
    png_path = joinpath(@__DIR__, "expected", "fig05.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig05] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig05.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig05")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["chi"], [v["ratio_$(k)"] for k in 1:length(SERIES)])
    println("[fig05] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves

for f in (:torque_survey, :torque_units, :mode_torques)
    isdefined(SW, f) || ReproHarness.hole("AccretionDisc.SpiralWaves.$(f) is not available (EGD::sw-survey-parallel).")
end
isdefined(SW, :newtonian_solve_mode) ||
    ReproHarness.hole("AccretionDisc.SpiralWaves.newtonian_solve_mode is not available: Fig. 5 needs the " *
                      "Newtonian-limit pipeline (EGD::sw-newtonian-limit) for its denominator.")

# ───────────────────────────── parameters (main.tex:927-932, brief section 4 row 5) ──────────────

const H0    = 0.05
const AEMRI = 12.0
const BSOFT = 0.6                  # nominal softening (M:884)
const EPS_P = 1e-5                 # nominal pole displacement (M:735)

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

# m ceiling of the mode sum. MEASURED convergence at these parameters (χ = 0, a• = 12M, h0 = 0.05,
# (0,0), grid (30,30,3000)): the cumulative torque reaches 92.6 % of its m → ∞ value at m = 40, 99.54 %
# at 80, 99.89 % at 100 and 99.98 % at 120, and m₉₅ = 46 (m₉₅h₀ = 2.3, inside the published 2.0-3.5 band
# of M:1025). ⌈4/h₀⌉ = 80 therefore truncates T by ≈ 0.5 %, well below the ≈ 3 % scatter of the published
# markers; the truncation actually reached is MEASURED per configuration below (`trunc_frac_max`) and
# registered, never assumed. (M:105 quotes ~3/h₀ modes for ~5 % accuracy — that would be too coarse here.)
const MMAX = _envint("GRTOOLS_FIG05_MMAX", ceil(Int, 4 / H0))

# grid resolution: MEASURED to give |ΔT_m|/T_m ≤ 5e-4 against the paper's reference (240,240,24000)
# (A:140) for m = 1 … 120 at these parameters — see the node report.
const GRIDKW = (nL = _envint("GRTOOLS_FIG05_NL", 30), nR = _envint("GRTOOLS_FIG05_NL", 30),
                nC = _envint("GRTOOLS_FIG05_NC", 3000))

# the published χ sampling: 10 points, −0.9 … 0.9 in steps of 0.2 (counted on the reference figure)
const NCHI = _envint("GRTOOLS_FIG05_NCHI", 10)
const CHIS = collect(range(-0.9, 0.9; length = NCHI))

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

"Total torque of one configuration in units of F_{J,0}h0, plus the tail fraction of the mode sum."
function total_torque(; χ, 𝔭, 𝔮, model)
    bg = SW.DiscBackground(; χ = χ, aEMRI = AEMRI, h0 = H0, 𝔭 = 𝔭, 𝔮 = 𝔮, b = BSOFT, εp = EPS_P)
    s = SW.torque_survey(bg, 1:MMAX; threads = true, model = model, grid_kwargs = GRIDKW)
    FJ0, _ = SW.torque_units(bg)
    Tm = SW.mode_torques(s)
    ntail = max(1, MMAX ÷ 20)
    (T = s.T_total / (FJ0 * H0), T_L = s.T_L_total / (FJ0 * H0),
     tail = sum(Tm[end-ntail+1:end]) / s.T_total,          # last 5 % of the modes
     res = maximum(mt.residual for mt in s.modes), seconds = s.seconds)
end

println("\n[fig05] a• = ", AEMRI, " M   h0 = ", H0, "   b = ", BSOFT, "   m = 1:", MMAX,
        "   grid = ", GRIDKW, "   χ = ", CHIS)
println("[fig05] compute budget: ", SW.compute_budget().note)

"""
    run_surveys() -> (ratios, Tgr, Tnewt, tailmax, resmax)

Every configuration of the figure: for each (𝔭,𝔮) series the ONE Newtonian torque and the ten
relativistic ones, plus the worst mode-sum tail fraction and solver residual over all of them.
"""
function run_surveys()
    ratios = Vector{Vector{Float64}}(undef, length(SERIES))
    Tgr    = Vector{Vector{Float64}}(undef, length(SERIES))
    Tnewt  = Vector{Float64}(undef, length(SERIES))
    tailmax = 0.0
    resmax  = 0.0
    for (k, s) in enumerate(SERIES)
        n = total_torque(; χ = 0.0, 𝔭 = s.𝔭, 𝔮 = s.𝔮, model = :newtonian)
        Tnewt[k] = n.T
        tailmax = max(tailmax, abs(n.tail)); resmax = max(resmax, n.res)
        @printf("[fig05] (𝔭,𝔮) = (%+.1f,%.1f)  T_Newt = %.6f F_J0h0  (T_L = %.4f, tail %.1e, %.0f s)\n",
                s.𝔭, s.𝔮, n.T, n.T_L, n.tail, n.seconds)
        flush(stdout)
        Tgr[k] = Float64[]
        for χ in CHIS
            g = total_torque(; χ = χ, 𝔭 = s.𝔭, 𝔮 = s.𝔮, model = :relativistic)
            push!(Tgr[k], g.T)
            tailmax = max(tailmax, abs(g.tail)); resmax = max(resmax, g.res)
            @printf("          χ = %+.2f  T_GR = %.6f  ratio = %.4f   (tail %.1e, res %.1e, %.0f s)\n",
                    χ, g.T, g.T / n.T, g.tail, g.res, g.seconds)
            flush(stdout)
        end
        ratios[k] = Tgr[k] ./ Tnewt[k]
    end
    (ratios, Tgr, Tnewt, tailmax, resmax)
end

t_start = time()
ratios, Tgr, Tnewt, tailmax, resmax = run_surveys()
@printf("[fig05] %d configurations in %.0f s (%.1f min)\n", (1 + NCHI) * length(SERIES), time() - t_start,
        (time() - t_start) / 60)

@repro_result "chi"        CHIS
@repro_result "a_emri"     AEMRI
@repro_result "h0"         H0
@repro_result "b_soft"     BSOFT
@repro_result "eps_p"      EPS_P
@repro_result "m_max"      MMAX
@repro_result "n_L"        GRIDKW.nL
@repro_result "n_C"        GRIDKW.nC
for (k, s) in enumerate(SERIES)
    @repro_result "p_$(k)"      s.𝔭
    @repro_result "q_$(k)"      s.𝔮
    @repro_result "T_newt_$(k)" Tnewt[k]
    @repro_result "T_gr_$(k)"   Tgr[k]
    @repro_result "ratio_$(k)"  ratios[k]
end
@repro_result "mode_sum_tail_fraction_max" tailmax
@repro_result "solver_residual_max"        resmax

# ───────────────────────────── the figure's own printed checks ───────────────────────────────────
# Fig. 5 prints no number of its own; its CLAIM (M:930-932) is the SIGN of the relativistic correction:
# ratio > 1 for the 𝔭 < 0 discs and ratio < 1 for the 𝔭 ≥ 0 ones, at every surveyed spin.

violations = String[]
for (k, s) in enumerate(SERIES)
    if s.𝔭 < 0
        all(>(1), ratios[k]) ||
            push!(violations, "(𝔭,𝔮) = ($(s.𝔭),$(s.𝔮)): T_GR/T_Newt ≤ 1 at χ = $(CHIS[findall(<=(1), ratios[k])])")
    else
        all(<(1), ratios[k]) ||
            push!(violations, "(𝔭,𝔮) = ($(s.𝔭),$(s.𝔮)): T_GR/T_Newt ≥ 1 at χ = $(CHIS[findall(>=(1), ratios[k])])")
    end
end
@repro_result "claim_p_negative_enhances" all(k -> all(>(1), ratios[k]), findall(s -> s.𝔭 < 0, collect(SERIES)))
@repro_result "claim_p_positive_suppresses" all(k -> all(<(1), ratios[k]), findall(s -> s.𝔭 ≥ 0, collect(SERIES)))

println("\n[fig05] figure checks (M:930-932)")
for (k, s) in enumerate(SERIES)
    ok = s.𝔭 < 0 ? all(>(1), ratios[k]) : all(<(1), ratios[k])
    @printf("  (𝔭,𝔮) = (%+.1f,%.1f)  ratio ∈ [%.3f, %.3f]  %s 1   %s\n", s.𝔭, s.𝔮,
            minimum(ratios[k]), maximum(ratios[k]), s.𝔭 < 0 ? ">" : "<", ok ? "PASS" : "FAIL")
end
@printf("  mode-sum tail (last 5%% of modes)  %.2e\n", tailmax)
@printf("  worst solver residual             %.2e\n", resmax)

# FIGURE-READ comparison (NOT a gate, and NOT stored in the golden): the marker values read off the
# published PDF at 200 dpi in this session, printed so the deviation is visible in the run log. These
# are hand-read positions with a ±0.015 reading error; they are never used to tune anything.
const FIGREAD = (
    [1.370, 1.455, 1.440, 1.450, 1.400, 1.330, 1.280, 1.240, 1.200, 1.155],
    [1.500, 1.515, 1.520, 1.470, 1.490, 1.465, 1.410, 1.385, 1.350, 1.315],
    [0.575, 0.680, 0.735, 0.800, 0.825, 0.845, 0.855, 0.875, 0.905, 0.930],
    [0.330, 0.435, 0.505, 0.585, 0.655, 0.695, 0.730, 0.775, 0.830, 0.880],
    [0.555, 0.625, 0.685, 0.695, 0.775, 0.800, 0.830, 0.845, 0.880, 0.925],
)
function report_figread(ratios)
    println("\n[fig05] deviation from the published markers (figure-read, ±0.015; diagnostic only)")
    worst = 0.0
    for (k, s) in enumerate(SERIES)
        d = ratios[k] .- FIGREAD[k]
        worst = max(worst, maximum(abs, d))
        @printf("  (%+.1f,%.1f)  mean %+.3f  max |Δ| %.3f  (Δ = %s)\n", s.𝔭, s.𝔮,
                sum(d) / length(d), maximum(abs, d), join([@sprintf("%+.2f", x) for x in d], " "))
    end
    @printf("  worst |Δ| over all 50 markers: %.3f\n", worst)
    worst
end
NCHI == 10 && report_figread(ratios)

check_failed = false
if !isempty(violations)
    println("\n[fig05] FINDING — figure claims not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw(CHIS, ratios)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig05", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_L" => (rtol = 0.0, atol = 0.0),
                                                "n_C" => (rtol = 0.0, atol = 0.0)))
exit(check_failed ? max(code, 1) : code)

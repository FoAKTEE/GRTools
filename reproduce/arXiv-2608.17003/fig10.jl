# reproduce/arXiv-2608.17003/fig10.jl — EGD::fig10-per-mode-reversal.
#
# arXiv:2608.17003 Fig. 10 (`fig:LR_torque_and_reversal`, main.tex:1007-1012): the PER-MODE Lindblad
# torque T_{L,m} (Eq. (75) evaluated for a single azimuthal mode) in units of F_{J,0}h₀ (76a), versus
# m, in Schwarzschild spacetime (χ = 0) at a• = 15M, 𝔮 = 0, b = 0.6, for the four combinations
# (𝔭, h₀) ∈ {(5, 0.03), (5, 0.05), (0, 0.03), (0, 0.05)}.
#
#   T_{L,m} = [F_J(r_out) − F_J(r_in)]_m       — 2608.17003 (75); GrTools:
#             `SpiralWaves.{torque_survey, mode_lindblad_torques}` (the survey solves each mode on its
#             own Appendix-C adaptive grid and reduces it with the AMF diagnostics (72a)).
#   unit    F_{J,0} h₀ with F_{J,0} = ε₀a•⁴Ω•²q²/h₀³  — (76a); GrTools: `SpiralWaves.torque_units`.
#
# The physics claim (M:1014) is the LOW-m SIGN: for 𝔭 = 5 the modal Lindblad torque is negative for
# m ≲ 10 (h₀ = 0.05) and m ≲ 15 (h₀ = 0.03), while the TOTAL torque stays positive because the high-m
# modes accumulate.  The totals are therefore summed over the FULL converged mode range 1:⌈4/h₀⌉
# (m₉₅h₀ ≲ 3.5, M:1025, so 4/h₀ clears the cutoff), not over the plotted window 4:40 — the plotted
# window follows the reference figure, whose first marker sits at m = 4.
#
# GOLDEN REUSE (checked, does not apply): no other driver of this directory carries per-mode Lindblad
# torques at a• = 15M with 𝔭 = 5 or h₀ = 0.03 (fig09 is h₀ = 0.05, 𝔭 ∈ {0}, a• ∈ [12,100]M and stores
# only the mode SUMS; fig08/fig13/fig14 store cutoffs / totals).  No number here is recomputed
# elsewhere.
#
# ZERO physics here: every number comes from those package calls.  Style: reproduce/_style/
# emri_style.jl + the VLM manual section 10 (progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── parameters (main.tex:1007-1012, brief row 10) ─────────────────────

const CHI   = 0.0                    # Schwarzschild
const AEMRI = 15.0                   # a• (M)
const QDISK = 0.0                    # 𝔮
const BSOFT = 0.6                    # nominal softening (M:884)
const M_PLOT = 4:40                  # plotted window (manual section 10: first marker at m = 4)
# (𝔭, h₀) in the LEGEND order of the reference figure (column-major: blue, orange | purple, brown)
const CONFIGS = ((5.0, 0.03), (5.0, 0.05), (0.0, 0.03), (0.0, 0.05))
const KEYS    = ("p5_h003", "p5_h005", "p0_h003", "p0_h005")

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
_envflt(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
# mode-sum ceiling for the TOTALS: ⌈K/h₀⌉ with K = 4 (m₉₅h₀ ≲ 3.5, M:1025)
const K_MMAX  = _envflt("GRTOOLS_FIG10_K", 4.0)
mmax_of(h0) = ceil(Int, K_MMAX / h0)

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

# manual section 0, palette 1 ("Pantone-like" categorical list), assignment from section 10
const BLUE   = RGB(0x0f / 255, 0x4c / 255, 0x81 / 255)   # (𝔭=5, h₀=0.03) ▽
const ORANGE = RGB(0xe2 / 255, 0x58 / 255, 0x3e / 255)   # (𝔭=5, h₀=0.05) △
const PURPLE = RGB(0x5f / 255, 0x4b / 255, 0x8b / 255)   # (𝔭=0, h₀=0.03) ■
const BROWN  = RGB(0xa4 / 255, 0x78 / 255, 0x64 / 255)   # (𝔭=0, h₀=0.05) ◆
const SERIES = ((BLUE, :dtriangle, "\$\\mathfrak{p}=5, h_0=0.03\$"),
                (ORANGE, :utriangle, "\$\\mathfrak{p}=5, h_0=0.05\$"),
                (PURPLE, :rect, "\$\\mathfrak{p}=0, h_0=0.03\$"),
                (BROWN, :diamond, "\$\\mathfrak{p}=0, h_0=0.05\$"))

"""
    draw(ms, curves) -> path

Render Fig. 10 from the four plotted arrays (`curves[k]` = T_{L,m}/(F_{J,0}h₀) over `ms`) and save the
PNG beside the golden.  Split out so that `--plot-only` can redraw from `expected/fig10.json` after a
style change WITHOUT repeating the mode sums (the numbers are the golden's, never recomputed).
"""
function draw(ms, curves)
    plt = plot(size = (470, 340), dpi = 160, legend = :top, legend_columns = 2, legendfontsize = 8,
               foreground_color_legend = nothing, background_color_legend = nothing,
               guidefontsize = 12, tickfontsize = 10,
               xlims = (0, 40), ylims = (-0.135, 0.150), xticks = 0:10:40, yticks = -0.10:0.05:0.10,
               xlabel = "\$m\$", ylabel = "\$T_{\\mathrm{L},m}\\ (F_{J,0}h_0)\$",
               minorticks = 5, tick_direction = :in, widen = false, bottom_margin = 2Plots.mm)
    hline!(plt, [0.0]; color = :black, lw = 2.0, label = "")          # manual: black zero line
    # GR fills a 2-column legend row-major; SERIES is in the reference's column-major legend order
    # (blue, orange | purple, brown), so the series are emitted as 1,3,2,4 to render that layout.
    for i in (1, 3, 2, 4)
        col, mk, lab = SERIES[i]
        scatter!(plt, ms, curves[i]; marker = (mk, 4.5), color = col, markerstrokewidth = 0,
                 label = lab)
    end
    png_path = joinpath(@__DIR__, "expected", "fig10.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig10] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig10.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

# `--plot-only` redraws the figure from the committed golden (style iteration without recomputing).
if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig10")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["m_plot"], [v["TL_" * k] for k in KEYS])
    println("[fig10] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves

for f in (:torque_survey, :mode_lindblad_torques, :mode_torques, :torque_units)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: the per-mode Lindblad torque of " *
        "arXiv:2608.17003 (75) needs the mode solver, the AMF/torque diagnostics and the survey " *
        "(EGD::sw-mode-solver, EGD::sw-diagnostics, EGD::sw-survey-parallel).")
end

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

"""
    run_config(𝔭, h0) -> (; ms, TL, T, TL_plot, unit, total_TL, total_T, m_last_neg, sign_changes, seconds)

One series of the figure: solve m = 1…⌈4/h₀⌉ on the paper's reference grid triple (240,240,24000)
(A:140, the package default), reduce to the per-mode Lindblad torques T_{L,m} (75) and the per-mode
total torques T_m (74), and normalise by F_{J,0}h₀ (76a).
"""
function run_config(𝔭, h0)
    t0 = time()
    bg = SW.DiscBackground(; χ = CHI, aEMRI = AEMRI, h0 = h0, 𝔭 = 𝔭, 𝔮 = QDISK, b = BSOFT)
    mmax = mmax_of(h0)
    surv = SW.torque_survey(bg, 1:mmax)
    unit = SW.torque_units(bg)[1] * bg.h0                       # F_{J,0} h₀
    TL = SW.mode_lindblad_torques(surv) ./ unit
    Tm = SW.mode_torques(surv) ./ unit
    neg = findall(<(0), TL)
    m_last_neg = isempty(neg) ? 0 : last(neg)
    sign_changes = count(i -> sign(TL[i]) != sign(TL[i+1]), 1:(length(TL) - 1))
    tail = sum(@view Tm[max(1, mmax - 9):mmax]) / sum(Tm)
    (; ms = collect(1:mmax), TL, Tm, unit, mmax,
       TL_plot = TL[M_PLOT], total_TL = sum(TL), total_T = sum(Tm), tail,
       TL_low = TL[1:3], TL_from4 = sum(@view TL[4:end]),
       m_last_neg, sign_changes, seconds = time() - t0, budget = surv.budget)
end

println("\n[fig10] χ = ", CHI, ", a• = ", AEMRI, " M, 𝔮 = ", QDISK, ", b = ", BSOFT,
        "; plotted m = ", first(M_PLOT), ":", last(M_PLOT), "; totals summed to ⌈", K_MMAX, "/h₀⌉")
flush(stdout)

res = Dict{String,Any}()
for (k, (𝔭, h0)) in zip(KEYS, CONFIGS)
    r = run_config(𝔭, h0)
    res[k] = r
    @printf("[fig10] 𝔭=%.0f h0=%.2f  m=1:%d  ΣT=%+.5f  ΣT_L=%+.5f  (Σ_{m≥4}T_L=%+.5f; T_L(1,2,3)=%+.4f,%+.4f,%+.4f)  T_L(4)=%+.5f  T_L(40)=%+.5f  last m with T_L<0 = %d  sign changes = %d  tail10=%.2e  %.1f s\n",
            𝔭, h0, r.mmax, r.total_T, r.total_TL, r.TL_from4, r.TL_low[1], r.TL_low[2], r.TL_low[3],
            r.TL[4], r.TL[40], r.m_last_neg, r.sign_changes, r.tail, r.seconds)
    flush(stdout)
end

# ───────────────────────────── registered results ────────────────────────────────────────────────

@repro_result "chi"            CHI
@repro_result "a_emri"         AEMRI
@repro_result "q_disk"         QDISK
@repro_result "b_soft"         BSOFT
@repro_result "m_plot"         collect(M_PLOT)
for k in KEYS
    register!("TL_" * k, res[k].TL_plot)                        # the plotted arrays (F_{J,0}h₀)
end
for k in KEYS
    register!("m_max_summed_" * k, res[k].mmax)
    register!("total_lindblad_" * k, res[k].total_TL)          # Σ_{m=1}^{m_max} T_{L,m}
    register!("lindblad_from_m4_" * k, res[k].TL_from4)        # Σ_{m≥4}: the PLOTTED mode range
    register!("total_torque_" * k, res[k].total_T)             # Σ_{m=1}^{m_max} T_m  = T of (74)
    register!("TL_low_m_" * k, res[k].TL_low)                  # m = 1,2,3 — below the plotted window
end
@repro_result "tail_fraction_max" maximum(abs(res[k].tail) for k in KEYS)

# ── the printed-number-gate contract keys (reproduce/arXiv-2608.17003/README.md) ──────────────────
@repro_result "m_last_negative_p5_h005"  res["p5_h005"].m_last_neg
@repro_result "m_last_negative_p5_h003"  res["p5_h003"].m_last_neg
# M:1014 "The total torque for 𝔭 = 5 is ultimately positive".  "The total torque" is T of Eq. (74)
# throughout arXiv:2608.17003 — it is the quantity the paper plots as "T (F_{J,0}h_0)" in Figs. 5, 6,
# 7 and 13 and the one Sec. III quotes numerically; the per-mode Lindblad torque of THIS figure is
# always written T_{L,m}.  The gate key therefore asserts Σ_m T_m > 0 at both thicknesses.
@repro_result "total_torque_p5_positive"    (res["p5_h005"].total_T > 0 && res["p5_h003"].total_T > 0)
# FINDING (registered, never hidden): the Lindblad-only mode sum Σ_m T_{L,m} is NEGATIVE for 𝔭 = 5 —
# the m = 1,2,3 modes, which the published figure does not plot, are strongly negative; restricted to
# the paper's own plotted range m ≥ 4 the Lindblad sum is positive.  Both are registered above.
@repro_result "total_lindblad_p5_positive"  (res["p5_h005"].total_TL > 0 && res["p5_h003"].total_TL > 0)
@repro_result "lindblad_from_m4_p5_positive" (res["p5_h005"].TL_from4 > 0 && res["p5_h003"].TL_from4 > 0)
# the 𝔭 = 0 series is stated to be positive mode by mode (M:1014, caption)
@repro_result "p0_all_modes_positive"    (all(>(0), res["p0_h005"].TL) && all(>(0), res["p0_h003"].TL))

# mode-sum convergence of the TOTALS: the last 10 of ⌈4/h₀⌉ modes must carry a negligible share
let tf = maximum(abs(res[k].tail) for k in KEYS)
    tf < 5e-2 || error("fig10: the top-10 modes still carry $(tf) of the total torque at m ≤ ⌈$(K_MMAX)/h₀⌉ " *
                       "— the mode sum is not converged; raise GRTOOLS_FIG10_K")
    @printf("[fig10] mode-sum convergence: top-10 modes carry at most %.2e of the total torque\n", tf)
end

# ───────────────────────────── printed-number checks (brief section 5, M:1014) ───────────────────
const PAPER_CHECKS_ENFORCED = true

struct PaperCheck; id::String; src::String; paper::Float64; tol::Float64; ours::Float64; end
checks = PaperCheck[
    PaperCheck("last m with T_L<0, 𝔭=5 h0=0.05", "M:1014", 10.0, 3.0, float(res["p5_h005"].m_last_neg)),
    PaperCheck("last m with T_L<0, 𝔭=5 h0=0.03", "M:1014", 15.0, 4.0, float(res["p5_h003"].m_last_neg)),
]
println("\n[fig10] printed-number checks")
println("  claim                                src      paper     ours      tol     verdict")
violations = String[]
for c in checks
    ok = abs(c.ours - c.paper) <= c.tol
    ok || push!(violations, "$(c.id) [$(c.src)]: paper $(c.paper), ours $(c.ours), tol ±$(c.tol)")
    @printf("  %-36s %-8s %7.3f  %7.3f  ±%.2f   %s\n", c.id, c.src, c.paper, c.ours, c.tol,
            ok ? "PASS" : "FAIL")
end
for (lab, ok) in (("total torque 𝔭=5 positive (M:1014)", res["p5_h005"].total_T > 0 && res["p5_h003"].total_T > 0),
                  ("every 𝔭=0 mode positive (M:1014)", all(>(0), res["p0_h005"].TL) && all(>(0), res["p0_h003"].TL)))
    ok || push!(violations, lab * ": FALSE")
    @printf("  %-36s %-8s %7s  %7s  %6s   %s\n", lab, "M:1014", "true", string(ok), "exact", ok ? "PASS" : "FAIL")
end
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig10] FINDING — printed numbers not reproduced within tolerance:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig10] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written")
        println("        or compared, so a long run is never thrown away on a paper mismatch")
        paper_check_failed = true
    else
        println("[fig10] PAPER_CHECKS_ENFORCED = false ⇒ recorded as a finding; the golden still gates our own numbers")
    end
end

println("\n[fig10] FINDING (not a paper claim — recorded for the ledger): the Lindblad-only mode sum")
for k in ("p5_h005", "p5_h003")
    @printf("   %s: Σ_{m≥1}T_L = %+.5f  Σ_{m≥4}T_L = %+.5f  T_L(1,2,3) = %+.4f, %+.4f, %+.4f  (F_{J,0}h0)\n",
            k, res[k].total_TL, res[k].TL_from4, res[k].TL_low[1], res[k].TL_low[2], res[k].TL_low[3])
end
println("   The m ≤ 3 modes are outside the published figure's range (its first marker is m = 4).")

draw(collect(M_PLOT), [res[k].TL_plot for k in KEYS])

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig10", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_last_negative_p5_h005" => (rtol = 0.0, atol = 0.0),
                                                "m_last_negative_p5_h003" => (rtol = 0.0, atol = 0.0)))
exit(paper_check_failed ? max(code, 1) : code)

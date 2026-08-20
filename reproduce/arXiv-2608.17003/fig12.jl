# reproduce/arXiv-2608.17003/fig12.jl — EGD::fig12-pole-errors.
#
# arXiv:2608.17003 Fig. 12 (`fig:pole_disp_errors`, appendix.tex:131-136): convergence of the TOTAL
# torque T (Eq. 74, summed over m) with the Landau pole-displacement parameter ε_p of (62)-(63).
# |T − T_ref|/T_ref versus the disc thickness h0 ∈ {0.15, 0.13, 0.11, 0.09, 0.07, 0.05} for
# ε_p ∈ {10⁻³, 10⁻⁴, 10⁻⁵}, against the reference solution at ε_p = 10⁻⁶.  Parameters:
# 𝔭 = 𝔮 = χ = 0, a• = 100M; grid n_L = n_R = 240, n_C = 2400 (A:162).
#
# Printed claim asserted here (A:162): "the relative errors scale as O(ε_p) for all values of h0" ⇒
# the ratio of the ε_p = 10⁻³ error to the ε_p = 10⁻⁴ error (and 10⁻⁴ to 10⁻⁵) is ≈ 10 at every h0.
# The registered numbers are the MEDIAN over h0 of those two ratios — the golden keys
# `relerr_ratio_1e3_over_1e4` and `relerr_ratio_1e4_over_1e5` of the printed-number gate.
#
# ZERO physics here: every torque is `AccretionDisc.SpiralWaves.torque_survey(...).T_total`, the
# fixed-order mode sum of (74) over the same m range for every (h0, ε_p), so the mode-sum truncation
# cancels in the difference and what is measured is the pole-displacement error alone.
# Style: `_style/emri_style.jl` + section 12 of `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md`.
#
# ── RESOLUTION FINDING (measured; why this driver does NOT use the grid printed at A:162) ─────────
# The Landau prescription (62)-(63) displaces the corotation pole by ε(r_c) = ε_p Ω•(r_c), i.e. by
#     Δr_pole = ε_p Ω•(r_c) / (m |dΩ/dr|_{r_c})
# in radius.  The (C4)-(C5) corotation block has σ_wc = 10⁻³r_c and n_C points, so its central
# spacing is Δr_C = σ_wc s_c /[sinh(s_c)] · 2/(n_C−1) — INDEPENDENT of m, while Δr_pole ∝ 1/(m ε_p⁻¹).
# The displaced pole is therefore unresolved above m ≈ ε_p Ω•/( |dΩ/dr| Δr_C ).
# MEASURED at a• = 100M, h0 = 0.15 (evidence `probe_pole.log`, single-mode torques):
#   n_C =  2400 → Δr_C = 1.222e-5 ; at ε_p = 10⁻⁶ the pole is resolved only for m ≲ 5
#                 error ratios 1e-3/1e-4, 1e-4/1e-5 = (9.99, 10.98) at m=1 → (0.86, 0.99) at m=30
#   n_C = 24000 → Δr_C = 1.222e-6 ; ratios (9.99, 10.99) m=1, (10.03, 10.98) m=3,
#                 (10.11, 11.03) m=10, (10.28, 14.02) m=30  — clean O(ε_p) at every m tested.
# So the O(ε_p) claim of A:162 IS reproduced, but only on a corotation block that resolves the
# ε_p = 10⁻⁶ reference pole.  On the grid printed at A:162 (n_L = n_R = 240, n_C = 2400) a mode sum
# to m = 120 saturates at the discretisation floor (measured: relative errors 6.8e-4, 1.2e-4, 2.0e-4
# at h0 = 0.15 — non-monotone in ε_p).  This driver therefore runs on the paper's REFERENCE triple
# n_C = 24000 (A:140) and truncates the mode sum at the paper's own cutoff m_max = ⌈3/h0⌉ (M:105),
# both quantities taken from the paper; the departure from A:162 is the finding, not a tuning knob,
# and `GRTOOLS_FIG12_NC` / `GRTOOLS_FIG12_MMAX` reproduce the saturated A:162 result on demand.

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

const BLUE   = RGB(0x0f / 255, 0x4c / 255, 0x81 / 255)   # ε_p = 1e-03
const ORANGE = RGB(0xe2 / 255, 0x58 / 255, 0x3e / 255)   # ε_p = 1e-04
const GREEN  = RGB(0x00 / 255, 0x94 / 255, 0x73 / 255)   # ε_p = 1e-05
const SERIES = ((BLUE, :star5, 6.5), (ORANGE, :rect, 4.5), (GREEN, :circle, 5.0))
const EPS_LABELS = ("\$\\epsilon_p=1e-03\$", "\$\\epsilon_p=1e-04\$", "\$\\epsilon_p=1e-05\$")

"""
    draw(h0s, rel_errs) -> path

Render Fig. 12 (reversed linear x, log y, three marker-only series) and save the PNG beside the
golden.  `--plot-only` redraws from `expected/fig12.json` without recomputing.
"""
function draw(h0s, rel_errs)
    plt = plot(size = (316, 232), dpi = 160, legend = :top, legendfontsize = 6,
               legend_columns = 3, foreground_color_legend = RGB(0.8, 0.8, 0.8),
               background_color_legend = :white, guidefontsize = 11, tickfontsize = 9,
               xflip = true, xlims = (0.04, 0.16), xticks = collect(h0s),
               yscale = :log10, ylims = (6e-6, 3.2e-3), yticks = [1e-5, 1e-4, 1e-3],
               minorticks = 5, tick_direction = :in, widen = false,
               xlabel = "\$h_0\$", ylabel = "\$|T-T_\\mathrm{ref}|/T_\\mathrm{ref}\$",
               bottom_margin = 1Plots.mm, left_margin = 1Plots.mm)
    for (i, (c, mk, ms)) in enumerate(SERIES)
        scatter!(plt, h0s, rel_errs[i]; color = c, marker = (mk, ms), markerstrokewidth = 0,
                 label = EPS_LABELS[i])
    end
    png_path = joinpath(@__DIR__, "expected", "fig12.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig12] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig12.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

const REL_KEYS = ("rel_err_1e3", "rel_err_1e4", "rel_err_1e5")

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig12")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["h0"], [v[k] for k in REL_KEYS])
    println("[fig12] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves

for f in (:torque_survey, :torque_units)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: Fig. 12 needs the threaded m-survey " *
        "(EGD::sw-survey-parallel) over the mode solver and diagnostics.")
end

# ───────────────────────────── parameters (appendix.tex:131-136, A:162) ──────────────────────────

const CHI   = 0.0
const AEMRI = 100.0
const PDISK = 0
const QDISK = 0
const BSOFT = 0.6
const H0S   = (0.15, 0.13, 0.11, 0.09, 0.07, 0.05)
const EPSS  = (1e-3, 1e-4, 1e-5)
const EPS_REF = 1e-6
_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
# Corotation-block resolution: the paper's REFERENCE triple (A:140), not the n_C = 2400 of A:162 —
# see the RESOLUTION FINDING in the file header.  `GRTOOLS_FIG12_NC=2400` restores A:162.
const NC    = _envint("GRTOOLS_FIG12_NC", 24000)
const GRID  = (nL = 240, nR = 240, nC = NC)
# Mode ceiling: the paper's own cutoff m_max ≈ 3/h (M:105 footnote), evaluated per h0 (the SAME
# ceiling for every ε_p at a given h0, so the mode-sum truncation cancels in |T − T_ref|).
# `GRTOOLS_FIG12_MMAX=<n>` overrides it with one fixed ceiling for every h0.
const MMAX_FIXED = _envint("GRTOOLS_FIG12_MMAX", 0)
m_ceiling(h0) = MMAX_FIXED > 0 ? MMAX_FIXED : ceil(Int, 3 / h0)

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

println("\n[fig12] a• = ", AEMRI, " M  χ = ", CHI, "  (𝔭,𝔮) = (", PDISK, ",", QDISK, ")")
println("[fig12] grid (n_L, n_R, n_C) = (", GRID.nL, ", ", GRID.nR, ", ", GRID.nC, ")   h0 = ", H0S,
        "   ε_p = ", EPSS, " vs reference ", EPS_REF)
println("[fig12] mode ceiling per h0 (m_max ≈ 3/h, M:105): ", [m_ceiling(h) for h in H0S])

run_case(h0, εp) = begin
    bg = SW.DiscBackground(; χ = CHI, aEMRI = AEMRI, h0 = h0, 𝔭 = PDISK, 𝔮 = QDISK, b = BSOFT, εp = εp)
    FJ0, _ = SW.torque_units(bg)
    s = SW.torque_survey(bg, 1:m_ceiling(h0); grid_kwargs = GRID)
    (T = s.T_total / (FJ0 * h0), seconds = s.seconds,
     residual = maximum(mt.residual for mt in s.modes))
end

"One disc thickness: the reference torque at ε_p = 1e-6 and the three displaced ones."
function one_h0(h0)
    r = run_case(h0, EPS_REF)
    eps = [run_case(h0, εp) for εp in EPSS]
    line = @sprintf("[fig12] h0 = %.2f   T(ε_p=1e-06) = %.10f (F_{J,0}h0)", h0, r.T)
    for (j, εp) in enumerate(EPSS)
        line *= @sprintf("   |Δ|/T(%.0e) = %.3e", εp, abs(eps[j].T - r.T) / abs(r.T))
    end
    println(line); flush(stdout)
    (T_ref = r.T, T = [e.T for e in eps], rel = [abs(e.T - r.T) / abs(r.T) for e in eps],
     residual = max(r.residual, maximum(e.residual for e in eps)),
     seconds = r.seconds + sum(e.seconds for e in eps))
end

cases = map(one_h0, H0S)
T_ref = [c.T_ref for c in cases]
rel   = [[c.rel[j] for c in cases] for j in eachindex(EPSS)]
T_eps = [[c.T[j]   for c in cases] for j in eachindex(EPSS)]
resmax = maximum(c.residual for c in cases)
wall   = sum(c.seconds for c in cases)

@repro_result "h0"          collect(H0S)
@repro_result "eps_p"       collect(EPSS)
@repro_result "eps_p_ref"   EPS_REF
for (k, v) in zip(REL_KEYS, rel); @repro_result k v; end
@repro_result "T_ref_FJ0h0" T_ref
@repro_result "T_1e3_FJ0h0" T_eps[1]
@repro_result "T_1e4_FJ0h0" T_eps[2]
@repro_result "T_1e5_FJ0h0" T_eps[3]
@repro_result "a_emri"      AEMRI
@repro_result "chi"         CHI
@repro_result "p_disk"      PDISK
@repro_result "q_disk"      QDISK
@repro_result "b_soft"      BSOFT
@repro_result "m_max"       [m_ceiling(h) for h in H0S]
@repro_result "n_L"         GRID.nL
@repro_result "n_C"         GRID.nC
@repro_result "solver_residual_max" resmax

# printed-number-gate contract keys (README table, claim A162-pole-scaling)
_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? float(s[(n+1)÷2]) : (s[n÷2] + s[n÷2+1]) / 2)
ratio_34 = _median(rel[1] ./ rel[2])
ratio_45 = _median(rel[2] ./ rel[3])
@repro_result "relerr_ratio_1e3_over_1e4" ratio_34
@repro_result "relerr_ratio_1e4_over_1e5" ratio_45
@repro_result "relerr_ratio_1e3_over_1e4_min" minimum(rel[1] ./ rel[2])
@repro_result "relerr_ratio_1e3_over_1e4_max" maximum(rel[1] ./ rel[2])
@repro_result "relerr_ratio_1e4_over_1e5_min" minimum(rel[2] ./ rel[3])
@repro_result "relerr_ratio_1e4_over_1e5_max" maximum(rel[2] ./ rel[3])

# ───────────────────────────── printed-number checks (A:162) ─────────────────────────────────────
const PAPER_CHECKS_ENFORCED = true

struct PaperCheck; id::String; src::String; paper::Float64; tol::Float64; ours::Float64; end
checks = PaperCheck[
    PaperCheck("O(ε_p): 1e-03/1e-04 error ratio", "A:162", 10.0, 20.0 / 3, ratio_34),
    PaperCheck("O(ε_p): 1e-04/1e-05 error ratio", "A:162", 10.0, 20.0 / 3, ratio_45),
]
println("\n[fig12] printed-number checks (A:162 — the errors scale as O(ε_p), factor-3 window)")
println("  claim                              src      paper     ours      tol     verdict")
violations = String[]
for c in checks
    ok = abs(c.ours - c.paper) <= c.tol
    ok || push!(violations, "$(c.id) [$(c.src)]: paper $(c.paper), ours $(round(c.ours; digits = 3)), tol ±$(c.tol)")
    @printf("  %-34s %-8s %7.3f  %7.3f  ±%.3f   %s\n", c.id, c.src, c.paper, c.ours, c.tol,
            ok ? "PASS" : "FAIL")
end
@printf("[fig12] per-h0 ratio ranges: 1e-3/1e-4 ∈ [%.2f, %.2f]   1e-4/1e-5 ∈ [%.2f, %.2f]\n",
        minimum(rel[1] ./ rel[2]), maximum(rel[1] ./ rel[2]),
        minimum(rel[2] ./ rel[3]), maximum(rel[2] ./ rel[3]))
@printf("[fig12] wall %.1f s (%d surveys, m ≤ 3/h0), max solver residual %.2e\n",
        wall, length(H0S) * (length(EPSS) + 1), resmax)
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig12] FINDING — printed numbers not reproduced within tolerance:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig12] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written")
        paper_check_failed = true
    else
        println("[fig12] PAPER_CHECKS_ENFORCED = false ⇒ recorded as a finding; the golden still gates our own numbers")
    end
end

draw(collect(H0S), rel)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig12", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_L" => (rtol = 0.0, atol = 0.0),
                                                "n_C" => (rtol = 0.0, atol = 0.0)))
exit(paper_check_failed ? max(code, 1) : code)

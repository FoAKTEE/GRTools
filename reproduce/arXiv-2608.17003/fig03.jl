# reproduce/arXiv-2608.17003/fig03.jl — EGD::fig03-balance.
#
# arXiv:2608.17003 Fig. 3 (`fig:angular_momentum_balance`, main.tex:903-908): the angular-momentum
# balance law (71) of the paper's linear theory, drawn as the two mode-summed radial densities
#
#   blue  d𝔗/dr  = √−g⟨𝒯⟩(r)                       Eq. (72b) integrand, Fourier average (73)
#   red   dF_J/dr, F_J = √−g[u_φ⟨δp δu^r⟩ + (ε+p)⟨δu^r δu_φ⟩]   Eq. (72a)
#
# both in units of F_{J,1} = F_{J,0}/a• (76b), versus r/a• ∈ [0.4, 3] for (𝔭,𝔮) = (0,0), χ = 0,
# h0 = 0.05, a• = 15M (r_ISCO/a• = 0.4).  They coincide everywhere except at corotation, where the
# EMRI DEPOSITS angular momentum in the disc: (71) then fails by exactly Δ_c𝔗 − Δ_cF_J, the
# corotation torque of (75).  That is the whole content of the figure and of its inset.
#
# ZERO physics here.  The disc, the modes, the fluxes and the torque densities all come from
# `AccretionDisc.SpiralWaves` (EGD::sw-background/grid/master/source/mode-solver/diagnostics); the
# mode sum and the reduction to the figure's own cells live in `_profile_assembly.jl` beside this
# file (shared with fig04, so it is written once).  Style: `_style/emri_style.jl` + section 3 of the
# VLM manual `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo).
#
# PLOTTING-GRID NOTE (presentation; recorded because a plotted number depends on it).  MEASURED here
# (probe at m ≤ 20, cells of 10⁻⁵ a•): across corotation F_J carries a STEP of 0.0266 F_{J,0} plus a
# positive bump narrower than 10⁻⁵ a•, while 𝔗 stays smooth.  A cell average of dF_J/dr over a cell
# of width Δ centred on r_c therefore reads (step + what is left of the bump inside Δ)/Δ: the plotted
# HEIGHT of the corotation spike is a property of the PLOT sampling (in the published figure exactly
# as here — the paper's ≈45 corresponds to Δ ≈ 6×10⁻⁴a•), while the resolution-INDEPENDENT content is
# the deposited angular momentum Δ_cF_J − Δ_c𝔗, which by (75) must equal T_L − T.  That identity is
# what this driver GATES; the spike height is registered as a presentation number only.
#
# The grid is anchored so that r_c is a cell CENTRE (otherwise the split of the near-delta between
# two cells is an accident of the offset) and is refined to Δ = 5×10⁻⁴ over r/a• ∈ [0.9, 1.1] so that
# the published inset (r/a• ∈ [0.99, 1.01]) is drawn with ≈40 cells; the smooth tails keep Δ = 2×10⁻³.

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

const BLUE   = RGB(0x0f / 255, 0x4c / 255, 0x81 / 255)   # manual section 0, palette 1, series 1
const ORANGE = RGB(0xe2 / 255, 0x58 / 255, 0x3e / 255)   # series 2

"""
    draw(x, dTdr, dFdr, x_isco, inset_ymax) -> path

Render Fig. 3 from the three plotted arrays and save the PNG beside the golden.  Split out so that
`--plot-only` can redraw from `expected/fig03.json` after a style change WITHOUT repeating the mode
sum (the numbers are then the golden's, never recomputed).
"""
function draw(x, dTdr, dFdr, x_isco, inset_ymax)
    base = 10.0 ^ floor(log10(max(inset_ymax, 1.0)))          # ~3 round inset ticks, as in the paper
    ytop = ceil(inset_ymax / base) >= 4 ? 2 * base : base
    ytick = filter(t -> t <= inset_ymax, collect(0:ytop:(10 * ytop)))
    plt = plot(size = (326, 233), dpi = 160, legend = :topright, legendfontsize = 8,
               guidefontsize = 10, tickfontsize = 8, foreground_color_legend = nothing,
               xlims = (0.33, 3.0), ylims = (-7, 8), xticks = 0.5:0.5:3.0, yticks = -5.0:2.5:7.5,
               minorticks = 5, tick_direction = :in, widen = false,
               xlabel = "\$r/a_\\bullet\$", ylabel = "Torque Balance",
               bottom_margin = 1Plots.mm, left_margin = 1Plots.mm)
    plot!(plt, x, dTdr; color = BLUE, lw = 2.0, label = "\$d\\mathfrak{T}/dr\\ (F_{J,1})\$")
    plot!(plt, x, dFdr; color = ORANGE, lw = 1.6, ls = :dash,
          label = "\$dF_J/dr\\ (F_{J,1})\$")
    vline!(plt, [x_isco]; color = :black, ls = :dashdot, lw = 1.1, label = "ISCO")

    # inset: the corotation neighbourhood (manual section 3 — lower-right quadrant, x 0.99…1.01)
    sel = findall(t -> 0.989 <= t <= 1.011, x)
    plot!(plt, x[sel], dTdr[sel]; inset = (1, bbox(0.545, 0.135, 0.345, 0.335, :bottom, :left)),
          subplot = 2, color = BLUE, lw = 1.4, label = "", framestyle = :box,
          xticks = ([0.99, 1.0, 1.01], ["0.99", "1.00", "1.01"]),
          yticks = ytick, tickfontsize = 6, tick_direction = :in, minorticks = 2,
          xlims = (0.989, 1.011), ylims = (-0.07 * inset_ymax, 1.10 * inset_ymax),
          grid = false, widen = false)
    plot!(plt, x[sel], dFdr[sel]; subplot = 2, color = ORANGE, lw = 1.2, ls = :dash, label = "")

    png_path = joinpath(@__DIR__, "expected", "fig03.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig03] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig03.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

# `--plot-only` redraws the figure from the committed golden (style iteration without recomputing).
if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig03")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["r_over_a"], v["dTdr_FJ1"], v["dFdr_FJ1"], v["r_isco_over_a"], v["corotation_spike_height_FJ1"])
    println("[fig03] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves
include(joinpath(@__DIR__, "_profile_assembly.jl"))
using .ProfileAssembly

for f in (:solve_mode, :amf, :torque_density, :cumulative_torque, :torque_units, :corotation_split)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: Fig. 3 needs the mode solver and the " *
        "AMF/torque diagnostics (EGD::sw-mode-solver, EGD::sw-diagnostics).")
end

# ───────────────────────────── parameters (main.tex:903-908, brief section 4 row 3) ──────────────

const CHI   = 0.0
const AEMRI = 15.0
const H0    = 0.05
const PDISK = 0                    # 𝔭
const QDISK = 0                    # 𝔮
const BSOFT = 0.6                  # nominal softening (M:884)
const EPS_P = 1e-5                 # nominal pole displacement (M:735)

# m ceiling of the mode sum.  m₉₅h₀ ≃ 2.0–3.5 (M:1025) ⇒ m₉₅ ≲ 70 at h0 = 0.05, and ~3/h ≈ 60 modes
# already give ~5 % of the total torque (M:105 footnote); 120 leaves a factor ≈ 2 of headroom.  The
# per-mode peak of the last modes is measured below and gates the choice.
const MMAX_DEFAULT = 120
_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
const MMAX = _envint("GRTOOLS_FIG03_MMAX", MMAX_DEFAULT)

# the figure's own sampling (see the header note): Δ_in over the corotation neighbourhood so that the
# inset gets ≈40 cells, Δ_out over the smooth tails, anchored so that r_c/a• = 1 is a cell centre
const XLO, XHI = 0.4, 3.0
const XMID_LO, XMID_HI = 0.9, 1.1
const DX_IN  = 5e-4
const DX_OUT = 2e-3

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

bg = SW.DiscBackground(; χ = CHI, aEMRI = AEMRI, h0 = H0, 𝔭 = PDISK, 𝔮 = QDISK, b = BSOFT, εp = EPS_P)
FJ0, FJ1 = SW.torque_units(bg)
x_isco = SW.r_isco(bg) / AEMRI

# the plotting grid is ANCHORED on corotation (r_c/a• = 1 for every m, since ω_d = m(Ω• − Ω) and
# Ω(r_c) = Ω•), refined there and coarse on the smooth tails — see `_profile_assembly.jl`.
let mid = centered_edges(1.0, XMID_LO, XMID_HI, DX_IN)
    global edges = vcat(graded_edges(((XLO, first(mid), DX_OUT),))[1:end-1], mid,
                        graded_edges(((last(mid), XHI, DX_OUT),))[2:end])
end
xc = cell_centers(edges)
println("\n[fig03] a• = ", AEMRI, " M   χ = ", CHI, "   h0 = ", H0, "   (𝔭,𝔮) = (", PDISK, ",", QDISK, ")")
println("[fig03] m = 1:", MMAX, "   plot cells = ", length(xc), " over r/a• ∈ [", XLO, ", ", XHI,
        "] (Δ = ", DX_IN, " on [", XMID_LO, ", ", XMID_HI, "], ", DX_OUT, " outside, anchored on r_c)")
println("[fig03] F_{J,0} = ", FJ0, "   F_{J,1} = ", FJ1, "   r_ISCO/a• = ", x_isco)

res = profile_sum(bg, 1:MMAX, edges .* AEMRI)
dT = res.dT_dr ./ FJ1
dF = res.dF_dr ./ FJ1
@printf("[fig03] %d modes in %.1f s (%d grid nodes, max solver residual %.2e, WKB-edge density fraction %.2e)\n",
        MMAX, res.seconds, res.nodes, res.residual_max, res.edge_density_fraction)

# extrema of the antisymmetric double peak around corotation, and the corotation cell
i_lo = argmin(dT); i_hi = argmax(dT)
i_spike = argmax(abs.(dF .- dT))
cell_w = edges[i_spike+1] - edges[i_spike]
spike = dF[i_spike]

# balance law away from corotation: the figure's claim, in the figure's own cells
far = findall(t -> abs(t - 1) > 0.02, xc)
near = findall(t -> abs(t - 1) <= 0.02, xc)
scale = maximum(abs, dT)
rel_far = maximum(abs.(dF[far] .- dT[far])) / scale

# the resolution-INDEPENDENT content of the spike: the angular momentum deposited at corotation,
# ∫(dF_J/dr − d𝔗/dr)dr over the corotation neighbourhood, in F_{J,0}.  By (75) it must equal
# T_L − T, which `profile_sum` returns from the package's own `corotation_split` — an independent
# cross-check of the figure against the diagnostics that produced it.
cellw = [edges[j+1] - edges[j] for j in eachindex(xc)]
deposition = sum((dF[near] .- dT[near]) .* cellw[near])                       # F_{J,0}
deposition_from_split = (res.T_L_total - res.T_total) / FJ0
deposition_far = sum((dF[far] .- dT[far]) .* cellw[far])
dep_rel = abs(deposition - deposition_from_split) / abs(deposition_from_split)

# mode-sum convergence: share of the last 10 modes in the per-mode peak torque density
tail = sum(res.per_mode_peak[max(1, MMAX - 9):MMAX]) / sum(res.per_mode_peak)

@printf("[fig03] trough %.4f at r/a• = %.4f   crest %.4f at r/a• = %.4f\n",
        dT[i_lo], xc[i_lo], dT[i_hi], xc[i_hi])
@printf("[fig03] corotation cell r/a• = %.5f (width %.2e): dF_J/dr = %.4g, d𝔗/dr = %.4g\n",
        xc[i_spike], cell_w, spike, dT[i_spike])
@printf("[fig03] corotation deposition ∫(dF_J/dr − d𝔗/dr)dr = %.8f F_J0 ; from (75) T_L − T = %.8f F_J0 (rel %.2e)\n",
        deposition, deposition_from_split, dep_rel)
@printf("[fig03] balance |dF_J/dr − d𝔗/dr|/max|d𝔗/dr| away from corotation (|r/a•−1| > 0.02): %.3e\n", rel_far)
@printf("[fig03] mode-sum tail: last 10 of %d modes carry %.2e of the summed per-mode peak\n", MMAX, tail)

@repro_result "r_over_a"        xc
@repro_result "dTdr_FJ1"        dT
@repro_result "dFdr_FJ1"        dF
@repro_result "a_emri"          AEMRI
@repro_result "chi"             CHI
@repro_result "h0"              H0
@repro_result "p_disk"          PDISK
@repro_result "q_disk"          QDISK
@repro_result "b_soft"          BSOFT
@repro_result "eps_p"           EPS_P
@repro_result "m_max"           MMAX
@repro_result "n_cells"         length(xc)
@repro_result "F_J0"            FJ0
@repro_result "F_J1"            FJ1
@repro_result "r_isco_over_a"   x_isco
@repro_result "trough_FJ1"      dT[i_lo]
@repro_result "trough_r_over_a" xc[i_lo]
@repro_result "crest_FJ1"       dT[i_hi]
@repro_result "crest_r_over_a"  xc[i_hi]
@repro_result "corotation_spike_r_over_a"    xc[i_spike]
@repro_result "corotation_spike_height_FJ1" spike
@repro_result "plot_cell_width_r_over_a"     cell_w
# resolution-INDEPENDENT invariant behind the spike: the angular momentum deposited at corotation,
# Δ_cF_J in units of F_{J,0} (height × cell width × a•, then /a• for F_{J,1} → F_{J,0})
@repro_result "corotation_deposition_FJ0"    deposition
@repro_result "corotation_deposition_from_split_FJ0" deposition_from_split
@repro_result "corotation_deposition_rel_diff"       dep_rel
@repro_result "deposition_outside_corotation_FJ0"    deposition_far
@repro_result "balance_rel_far"              rel_far
@repro_result "T_total_FJ0h0"                res.T_total / (FJ0 * H0)
@repro_result "T_L_total_FJ0h0"              res.T_L_total / (FJ0 * H0)
@repro_result "mode_sum_tail_fraction"       tail
@repro_result "solver_residual_max"          res.residual_max
@repro_result "wkb_edge_density_fraction"    res.edge_density_fraction

# ───────────────────────────── the figure's own printed checks ───────────────────────────────────
# Fig. 3 prints no number of its own; its CLAIM is the balance law (71) away from corotation, which
# the lane task fixes at rel ≤ 1e-4, plus the qualitative statement that the deposition is confined
# to corotation.  Both are asserted here; a violation exits 1.
const BALANCE_TOL = 1e-4
violations = String[]
rel_far <= BALANCE_TOL ||
    push!(violations, "balance away from corotation: $(rel_far) > $(BALANCE_TOL)")
abs(xc[i_spike] - 1) <= 2 * cell_w ||
    push!(violations, "the AMF/torque discrepancy peaks at r/a• = $(xc[i_spike]), not at corotation")
abs(spike) > scale ||
    push!(violations, "no corotation spike: |dF_J/dr| there is only $(abs(spike)/scale)× the profile scale")
dep_rel <= 1e-3 ||
    push!(violations, "the plotted deposition $(deposition) F_J0 disagrees with T_L − T = " *
                      "$(deposition_from_split) F_J0 by $(dep_rel) (relative)")
abs(deposition_far) <= 0.01 * abs(deposition) ||
    push!(violations, "the deposition is not confined to corotation: $(deposition_far) F_J0 of " *
                      "$(deposition) F_J0 falls outside |r/a•−1| ≤ 0.02")
tail < 5e-3 ||
    push!(violations, "mode sum not converged: the last 10 modes carry $(tail) of the summed peak")
# the implicit-midpoint solve must be CONVERGED; its residual is roundoff-level (~7e-9) and is gated
# here as a BOUND rather than against a stored 12-digit golden value (see the tolerances below)
const RESIDUAL_BOUND = 1e-6
res.residual_max <= RESIDUAL_BOUND ||
    push!(violations, "the implicit-midpoint solve did not converge: max residual " *
                      "$(res.residual_max) > $(RESIDUAL_BOUND)")

println("\n[fig03] figure checks")
@printf("  balance away from corotation   %.3e  (≤ %.0e)   %s\n", rel_far, BALANCE_TOL,
        rel_far <= BALANCE_TOL ? "PASS" : "FAIL")
@printf("  discrepancy peak at corotation r/a• = %.5f (cell %.1e)          %s\n", xc[i_spike], cell_w,
        abs(xc[i_spike] - 1) <= 2 * cell_w ? "PASS" : "FAIL")
@printf("  spike/profile-scale            %.1f   (> 1)                   %s\n", abs(spike) / scale,
        abs(spike) > scale ? "PASS" : "FAIL")
@printf("  deposition = T_L − T           %.2e  (rel ≤ 1e-3)            %s\n", dep_rel,
        dep_rel <= 1e-3 ? "PASS" : "FAIL")
@printf("  deposition confined to r≈r_c   %.2e  (≤ 1%% of the total)     %s\n",
        abs(deposition_far) / abs(deposition),
        abs(deposition_far) <= 0.01 * abs(deposition) ? "PASS" : "FAIL")
@printf("  mode-sum tail                  %.2e  (< 5e-3)                %s\n", tail,
        tail < 5e-3 ? "PASS" : "FAIL")
@printf("  implicit-midpoint residual     %.2e  (≤ %.0e)              %s\n", res.residual_max,
        RESIDUAL_BOUND, res.residual_max <= RESIDUAL_BOUND ? "PASS" : "FAIL")
@repro_result "figure_checks_pass" isempty(violations)

check_failed = false
if !isempty(violations)
    println("\n[fig03] FINDING — figure claims not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw(xc, dT, dF, x_isco, abs(spike))

# Per-key tolerances.  `m_max`, `n_cells` are exact integers.  The other three are NOISE-DERIVED
# quantities, not predictions, and are gated PHYSICALLY by the driver assertions above (`dep_rel ≤
# 1e-3`, `|deposition_far| ≤ 1 % of the deposition`, `residual_max ≤ RESIDUAL_BOUND`); the golden
# comparison on them is a determinism check, so its tolerance must reflect their own conditioning:
#   • `deposition_outside_corotation_FJ0` — a residue after the corotation window is removed (1e-6,
#     established when the golden was first written);
#   • `corotation_deposition_rel_diff` — (Δ_cF_J − [T_L−T])/Δ_cF_J, a ratio of two 0.0266 numbers
#     that agree to 2.7e-5, so a relative perturbation δ of each can move it by up to δ/2.7e-5;
#   • `solver_residual_max` — the max implicit-midpoint residual, ~7e-9, i.e. ROUNDOFF of the solve.
# EGD::jd-swap measured both against the route-R3 master entries: the physical keys moved by ≤ 6.7e-9
# (all inside rtol 1e-8) while these two moved by 2.5e-7 and 1.9e-2 — see the node evidence. Gating
# roundoff against a stored 12-digit value is what failed, not the physics. [FINDING]
code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig03", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_cells" => (rtol = 0.0, atol = 0.0),
                                                "deposition_outside_corotation_FJ0" =>
                                                    (rtol = 1e-6, atol = 1e-12),
                                                "corotation_deposition_rel_diff" =>
                                                    (rtol = 1e-2, atol = 1e-9),
                                                "solver_residual_max" =>
                                                    (rtol = 1e-1, atol = 1e-11)))
exit(check_failed ? max(code, 1) : code)

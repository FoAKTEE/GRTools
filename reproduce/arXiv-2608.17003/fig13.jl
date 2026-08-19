# reproduce/arXiv-2608.17003/fig13.jl — EGD::fig13-torque-vs-b.
#
# arXiv:2608.17003 Fig. 13 (`fig:variation_with_respect_to_b`, appendix.tex:143-149): the TOTAL torque
# T (74) in units of F_{J,0}h₀ (76a) versus the Plummer softening parameter b ∈ {0.3,…,0.8} (59), for
# the disc thicknesses h₀ ∈ {0.05, 0.075, 0.1, 0.125, 0.15}, with 𝔭 = 𝔮 = χ = 0 and a• = 100M.
#
#   T      = 𝔗(r_out) − 𝔗(r_in) summed over m            — 2608.17003 (74); GrTools:
#            `SpiralWaves.torque_survey(...).T_total` (per-mode Appendix-C grids + the (72b) torque
#            density; the m-sum is a fixed-order serial reduction, bitwise reproducible).
#   b      = the softening of (59), ε_soft = a• h₀ b     — GrTools: `DiscBackground(; b = …)`.
#   unit   F_{J,0}h₀, F_{J,0} = ε₀a•⁴Ω•²q²/h₀³ (76a)     — GrTools: `SpiralWaves.torque_units`.
#
# The appendix claim (A:164) is MONOTONE DECREASE with b; the wave-0 brief §5 reads ≈1.28 → 0.97
# F_{J,0}h₀ for h₀ = 0.05 off the published figure (figure-read provenance, looser tolerance).
#
# GOLDEN REUSE (checked, does not apply): fig07 also sweeps h₀ at 𝔭 = 𝔮 = χ = 0 but at the NOMINAL
# softening b = 0.6 and versus a•, and fig12 sweeps ε_p at a• = 100M with a fixed b; neither stores a
# b-sweep.  fig14 shares the b grid but a• = 30M and stores m₉₅, not the torque.
#
# ZERO physics here: every number comes from those package calls.  Style: reproduce/_style/
# emri_style.jl + the VLM manual section 13 (progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── parameters (appendix.tex:143-149, brief row 13) ───────────────────

const CHI   = 0.0
const AEMRI = 100.0                  # a• (M)
const PDISK = 0.0
const QDISK = 0.0
const BS    = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
const H0S   = (0.05, 0.075, 0.1, 0.125, 0.15)

_envflt(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
# mode-sum ceiling ⌈K/h₀⌉ with K = 4 (m₉₅h₀ ≲ 3.5, M:1025, so 4/h₀ clears the cutoff)
const K_MMAX = _envflt("GRTOOLS_FIG13_K", 4.0)
mmax_of(h0) = ceil(Int, K_MMAX / h0)
# GRID (recorded, VERIFIED, not assumed): 30 configurations × ⌈4/h₀⌉ modes at the paper's reference
# triple (240,240,24000) (A:140) costs ≈ 1 h at -t8, and the mission budget must also pay for the
# `--check` re-run.  The surveys therefore use (60,60,6000), one of the paper's OWN convergence-run
# resolutions (A:139).  `RES_CHECK` below re-runs two points at the reference triple inside this
# driver and requires the total torque to agree to 1e-3 relative (a pre-check measured 6e-6).
const GRID     = (nL = _envint("GRTOOLS_FIG13_NL", 60), nR = _envint("GRTOOLS_FIG13_NR", 60),
                  nC = _envint("GRTOOLS_FIG13_NC", 6000))
const REF_GRID = (nL = 240, nR = 240, nC = 24000)

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

# manual section 0, palette 2: plasma(linspace(0.1,0.9,5)) for the h₀ family, hollow markers
const PLASMA5 = (RGB(0x41 / 255, 0x04 / 255, 0x9d / 255),   # h₀ = 0.050  square
                 RGB(0x8f / 255, 0x0d / 255, 0xa4 / 255),   # h₀ = 0.075  star
                 RGB(0xcc / 255, 0x47 / 255, 0x78 / 255),   # h₀ = 0.100  circle
                 RGB(0xf2 / 255, 0x84 / 255, 0x4b / 255),   # h₀ = 0.125  triangle-right
                 RGB(0xfc / 255, 0xce / 255, 0x25 / 255))   # h₀ = 0.150  hexagon
const MARKERS = (:rect, :star5, :circle, :rtriangle, :hexagon)

"""
    draw(bs, curves) -> path

Render Fig. 13 from the five plotted arrays (`curves[i]` = T/(F_{J,0}h₀) over `bs` for `H0S[i]`) and
save the PNG beside the golden.  Split out so that `--plot-only` can redraw from `expected/fig13.json`
after a style change WITHOUT repeating the mode sums.
"""
function draw(bs, curves)
    plt = plot(size = (470, 340), dpi = 160, legend = :topright, legend_columns = 2, legendfontsize = 8,
               foreground_color_legend = RGB(0.8, 0.8, 0.8), background_color_legend = :white,
               guidefontsize = 12, tickfontsize = 10,
               xlims = (0.28, 0.82), ylims = (0.95, 1.40), xticks = 0.3:0.1:0.8, yticks = 1.0:0.1:1.3,
               xlabel = "\$b\$", ylabel = "\$T\\ (F_{J,0}h_0)\$",
               minorticks = 5, tick_direction = :in, widen = false, bottom_margin = 2Plots.mm)
    # The reference markers are matplotlib `mfc='none'` (transparent fill).  GR ignores the alpha of a
    # marker fill colour, so the fill is painted WHITE and the DRAW ORDER is chosen so that no series
    # hides another: the two coincident purple series (h₀ = 0.050 square, 0.075 star) are drawn
    # largest-first, which reproduces the reference's "star inside a square".
    for (i, ms) in ((5, 6.5), (4, 6.5), (3, 6.5), (1, 7.5), (2, 5.6))
        scatter!(plt, bs, curves[i]; marker = (MARKERS[i], ms), markercolor = :white,
                 markerstrokecolor = PLASMA5[i], markerstrokewidth = 1.6, label = "")
    end
    # legend entries are DUMMY series (no data) so that the legend order is independent of the draw
    # order: GR fills a 2-column legend row-major, so this order renders the reference's column-major
    # descending layout (0.150, 0.125, 0.100 | 0.075, 0.050).
    for i in (5, 2, 4, 1, 3)
        scatter!(plt, [NaN], [NaN]; marker = (MARKERS[i], 6.5), markercolor = :white,
                 markerstrokecolor = PLASMA5[i], markerstrokewidth = 1.6,
                 label = "\$h_0=" * @sprintf("%.3f", H0S[i]) * "\$")
    end
    png_path = joinpath(@__DIR__, "expected", "fig13.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig13] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig13.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig13")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["b_grid"], [v["T_h0_" * replace(@sprintf("%.3f", h), "." => "")] for h in H0S])
    println("[fig13] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves

for f in (:torque_survey, :mode_torques, :torque_units)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: the total torque of arXiv:2608.17003 (74) " *
        "needs the mode solver, the torque diagnostics and the survey (EGD::sw-mode-solver, " *
        "EGD::sw-diagnostics, EGD::sw-survey-parallel).")
end

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

key(h0) = "T_h0_" * replace(@sprintf("%.3f", h0), "." => "")

"""
    neglected_tail(Tm) -> fraction

Estimate of the mode-sum truncation error: the share of the total torque carried by the modes BEYOND
the ceiling, extrapolated geometrically from the last two computed modes
(ρ = |T_{m_max}|/|T_{m_max−1}| ⇒ Σ_{m>m_max}T_m ≈ |T_{m_max}| ρ/(1−ρ)).  The naive "share of the last
ten modes" is NOT a truncation measure when m_max is small (at h₀ = 0.15, m_max = 27, ten modes are
more than a third of the whole sum), so both are registered but only this one is gated.
"""
function neglected_tail(Tm::AbstractVector)
    n = length(Tm); n < 3 && return NaN
    tot = sum(Tm); iszero(tot) && return NaN
    ρ = abs(Tm[n]) / max(abs(Tm[n-1]), eps(float(one(eltype(Tm)))))
    ρ = min(ρ, 0.999)
    abs(Tm[n]) * ρ / (1 - ρ) / abs(tot)
end


"""
    run_point(h0, b) -> (; T, T_L, tail, mmax, seconds)

One marker of the figure: the total torque (74) and Lindblad torque (75) summed over m = 1…⌈4/h₀⌉ at
softening `b`, normalised by F_{J,0}h₀ (76a).  `grid` is forwarded to `build_grid` ((60,60,6000) by
default, `REF_GRID` for the in-driver resolution cross-check).
"""
function run_point(h0, b; grid = GRID)
    t0 = time()
    bg = SW.DiscBackground(; χ = CHI, aEMRI = AEMRI, h0 = h0, 𝔭 = PDISK, 𝔮 = QDISK, b = b)
    mmax = mmax_of(h0)
    surv = SW.torque_survey(bg, 1:mmax; grid_kwargs = grid)
    unit = SW.torque_units(bg)[1] * bg.h0
    Tm = SW.mode_torques(surv)
    (; T = surv.T_total / unit, T_L = surv.T_L_total / unit, mmax,
       tail = sum(@view Tm[max(1, mmax - 9):mmax]) / surv.T_total, trunc = neglected_tail(Tm),
       seconds = time() - t0)
end

println("\n[fig13] χ = ", CHI, ", a• = ", AEMRI, " M, (𝔭,𝔮) = (", PDISK, ",", QDISK,
        "); b = ", BS, "; h₀ = ", H0S, "; m = 1:⌈", K_MMAX, "/h₀⌉; grid = ", GRID)
flush(stdout)

curves = [Float64[] for _ in H0S]
curvesL = [Float64[] for _ in H0S]
tails = Float64[]; truncs = Float64[]
for (i, h0) in enumerate(H0S), b in BS
    p = run_point(h0, b)
    push!(curves[i], p.T); push!(curvesL[i], p.T_L); push!(tails, p.tail); push!(truncs, p.trunc)
    @printf("[fig13] h0=%.3f b=%.2f  m=1:%3d  T=%.5f  T_L=%.5f  (F_{J,0}h0)  tail10=%.2e  trunc=%.2e  %.1f s\n",
            h0, b, p.mmax, p.T, p.T_L, p.tail, p.trunc, p.seconds)
    flush(stdout)
end

# resolution cross-check: the same torque at the paper's reference triple (240,240,24000)
const RES_CHECK = ((0.05, 0.3), (0.15, 0.8))
res_red = Float64[]; res_ref = Float64[]
for (h0, b) in RES_CHECK
    a = run_point(h0, b)
    r = run_point(h0, b; grid = REF_GRID)
    push!(res_red, a.T); push!(res_ref, r.T)
    @printf("[fig13] resolution check h0=%.3f b=%.2f: T = %.6f at %s vs %.6f at %s (rel %.2e)\n",
            h0, b, a.T, string(GRID), r.T, string(REF_GRID), abs(a.T - r.T) / abs(r.T))
    flush(stdout)
end
resolution_ok = all(abs(res_red[i] - res_ref[i]) <= 1e-3 * abs(res_ref[i]) for i in eachindex(res_ref))

# ───────────────────────────── registered results ────────────────────────────────────────────────

@repro_result "chi"      CHI
@repro_result "a_emri"   AEMRI
@repro_result "p_disk"   PDISK
@repro_result "q_disk"   QDISK
@repro_result "b_grid"   collect(BS)
@repro_result "h0_grid"  collect(H0S)
@repro_result "m_max_over_h0" K_MMAX
for (i, h0) in enumerate(H0S)
    register!(key(h0), curves[i])
    register!("TL_h0_" * replace(@sprintf("%.3f", h0), "." => ""), curvesL[i])
end
@repro_result "tail_fraction_max" maximum(abs, tails)
@repro_result "truncation_estimate_max" maximum(abs, truncs)
@repro_result "all_decreasing" all(all(diff(c) .< 0) for c in curves)
@repro_result "grid_nL"   GRID.nL
@repro_result "grid_nC"   GRID.nC
@repro_result "T_reduced_resolution" res_red
@repro_result "T_ref_resolution"     res_ref
@repro_result "resolution_check_pass" resolution_ok

# ── the printed-number-gate contract keys (reproduce/arXiv-2608.17003/README.md) ──────────────────
@repro_result "torque_h005_b03" curves[1][1]
@repro_result "torque_h005_b08" curves[1][end]

@printf("[fig13] mode-sum convergence: top-10 modes carry at most %.2e of the total torque; the geometric estimate of the NEGLECTED sum beyond m ≤ ⌈%.1f/h₀⌉ is at most %.2e\n",
        maximum(abs, tails), K_MMAX, maximum(abs, truncs))

# ───────────────────────────── printed-number checks (brief section 5, A:164) ────────────────────
const PAPER_CHECKS_ENFORCED = true

struct PaperCheck; id::String; src::String; paper::Float64; tol::Float64; ours::Float64; end
checks = PaperCheck[
    PaperCheck("T at h0=0.05, b=0.3 ≈ 1.28", "A:164", 1.28, 0.12, curves[1][1]),
    PaperCheck("T at h0=0.05, b=0.8 ≈ 0.97", "A:164", 0.97, 0.12, curves[1][end]),
]
println("\n[fig13] printed-number checks")
println("  claim                                src      paper     ours      tol     verdict")
violations = String[]
for c in checks
    ok = abs(c.ours - c.paper) <= c.tol
    ok || push!(violations, "$(c.id) [$(c.src)]: paper $(c.paper), ours $(round(c.ours; digits=4)), tol ±$(c.tol)")
    @printf("  %-36s %-8s %7.4f  %7.4f  ±%.3f   %s\n", c.id, c.src, c.paper, c.ours, c.tol,
            ok ? "PASS" : "FAIL")
end
let ok = all(all(diff(c) .< 0) for c in curves)
    ok || push!(violations, "monotone decrease of T with b (A:164): FALSE")
    @printf("  %-36s %-8s %7s  %7s  %6s   %s\n", "T decreases with b (A:164)", "A:164", "true",
            string(ok), "exact", ok ? "PASS" : "FAIL")
end
let tf = maximum(abs, truncs), ok = tf < 5e-2
    ok || push!(violations, "mode-sum truncation: the geometric estimate of the neglected sum is $(tf) > 5e-2; raise GRTOOLS_FIG13_K")
    @printf("  %-36s %-8s %7s  %7.1e  %6s   %s\n", "neglected mode-sum tail", "M:105", "<5e-2", tf,
            "5e-2", ok ? "PASS" : "FAIL")
end
let ok = resolution_ok
    ok || push!(violations, "resolution check: T at $(GRID) = $(res_red) vs $(res_ref) at $(REF_GRID) (rtol 1e-3)")
    @printf("  %-36s %-8s %7s  %7s  %6s   %s\n", "T unchanged at the reference grid", "A:140", "true",
            string(ok), "1e-3", ok ? "PASS" : "FAIL")
end
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig13] FINDING — printed numbers not reproduced within tolerance:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig13] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written or compared")
        paper_check_failed = true
    end
end

draw(collect(BS), curves)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig13", script = @__FILE__, rtol = 1e-8)
exit(paper_check_failed ? max(code, 1) : code)

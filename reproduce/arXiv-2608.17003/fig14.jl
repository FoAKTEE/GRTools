# reproduce/arXiv-2608.17003/fig14.jl — EGD::fig14-m95-vs-b.
#
# arXiv:2608.17003 Fig. 14 (`fig:m95_vs_b`, appendix.tex:153-160): the torque-cutoff parameter m₉₅
# scaled by h₀ versus the Plummer softening parameter b ∈ {0.3,…,0.8} (59), for 𝔭 = 𝔮 = χ = 0,
# a• = 30M, h₀ = 0.05 (a single series).
#
#   m₉₅  = the SMALLEST azimuthal mode number at which the cumulative torque reaches 95 % of the
#          final torque obtained after summing all computed modes (M:990, verbatim); GrTools:
#          `SpiralWaves.m_cut(mode_torques(survey), 0.95)` over the per-mode TOTAL torques (74) —
#          "the total cumulative torque" of the Fig.-8 caption (M:958).
#   b    = the softening of (59), ε_soft = a• h₀ b     — GrTools: `DiscBackground(; b = …)`.
#
# The appendix claim (A:164, brief §5) is a MONOTONE DECREASE of m₉₅h₀ from ≈2.9 (b = 0.3) to ≈1.95
# (b = 0.8): a smaller softening puts more power in the high-m forcing modes.
#
# m_max SENSITIVITY (recorded, not hidden): m₉₅ is defined against "the final torque obtained after
# summing all computed modes", so it inherits a weak dependence on the mode ceiling.  The driver sums
# to m ≤ ⌈K/h₀⌉ with K = 4 (m₉₅h₀ ≲ 3.5, M:1025) and RE-RUNS the b = 0.3 point at K = 6, registering
# both so the sensitivity is visible in the golden.
#
# GOLDEN REUSE (checked, does not apply): fig08 stores m₉₅h₀ samples at the NOMINAL b = 0.6 over
# (h₀, χ, a•) — never a b sweep — and its a• grid is {15,30,120}M; fig13 shares this b grid but at
# a• = 100M and stores the torque, not the cutoff.  No number here is recomputed elsewhere.
#
# ZERO physics here.  Style: reproduce/_style/emri_style.jl + the VLM manual section 14
# (progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── parameters (appendix.tex:153-160, brief row 14) ───────────────────

const CHI   = 0.0
const AEMRI = 30.0                   # a• (M)
const PDISK = 0.0
const QDISK = 0.0
const H0    = 0.05
const BS    = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

_envflt(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
# PRIMARY ceiling ⌈6/h₀⌉ since wave 6 (decision D2, resolved in DECISIONS_W6.md); the pre-registered
# ⌈4/h₀⌉ stays REGISTERED beside it at the b = 0.3 point as `m95_h0_Ksens_b03`.
const K_MMAX = _envflt("GRTOOLS_FIG14_K", 6.0)      # mode ceiling ⌈K/h₀⌉  (PRIMARY)
const K_SENS = _envflt("GRTOOLS_FIG14_KSENS", 4.0)  # sensitivity ceiling for the b = 0.3 point
# GRID (recorded, VERIFIED, not assumed): the surveys run at (60,60,6000), one of the paper's own
# convergence-run resolutions (A:139), because the reference triple (240,240,24000) (A:140) costs
# ≈ 3× the mission budget for `--regen` plus `--check`.  `RES_CHECK` re-runs the b = 0.3 and b = 0.8
# points at the reference triple inside this driver and requires the SAME m₉₅.
const GRID     = (nL = _envint("GRTOOLS_FIG14_NL", 60), nR = _envint("GRTOOLS_FIG14_NR", 60),
                  nC = _envint("GRTOOLS_FIG14_NC", 6000))
const REF_GRID = (nL = 240, nR = 240, nC = 24000)

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

const BLUE = RGB(0x0f / 255, 0x4c / 255, 0x81 / 255)   # manual section 0, palette 1, series 1

"""
    draw(bs, m95h0) -> path

Render Fig. 14 (one blue filled-circle series, no legend, no connecting line) and save the PNG beside
the golden.  Split out so `--plot-only` can redraw from `expected/fig14.json` without recomputing.
"""
function draw(bs, m95h0)
    lo, hi = minimum(m95h0), maximum(m95h0)
    pad = max(0.05, 0.08 * (hi - lo))
    plt = plot(size = (470, 340), dpi = 160, legend = false,
               guidefontsize = 12, tickfontsize = 10,
               xlims = (0.28, 0.82), ylims = (lo - pad, hi + pad), xticks = 0.3:0.1:0.8,
               xlabel = "\$b\$", ylabel = "\$m_{95}\\,h_0\$",
               minorticks = 5, tick_direction = :in, widen = false, bottom_margin = 2Plots.mm)
    scatter!(plt, bs, m95h0; marker = (:circle, 7), color = BLUE, markerstrokewidth = 0, label = "")
    png_path = joinpath(@__DIR__, "expected", "fig14.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig14] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig14.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig14")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["b_grid"], v["m95_h0"])
    println("[fig14] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves

for f in (:torque_survey, :mode_torques, :m_cut)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: the m₉₅ cutoff of arXiv:2608.17003 (M:990) " *
        "needs the mode solver, the torque diagnostics and the survey (EGD::sw-mode-solver, " *
        "EGD::sw-diagnostics, EGD::sw-survey-parallel).")
end

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

"""
    neglected_tail(Tm) -> fraction

Estimate of the mode-sum truncation error: the share of the total torque carried by the modes BEYOND
the ceiling, extrapolated geometrically from the last two computed modes
(ρ = |T_{m_max}|/|T_{m_max−1}| ⇒ Σ_{m>m_max}T_m ≈ |T_{m_max}| ρ/(1−ρ)).  The naive "share of the last
ten modes" is NOT a truncation measure when m_max is small (at h₀ = 0.15 the ceiling is 27 modes, so
ten of them are more than a third of the whole sum), so both are registered but only this one is
gated.
"""
function neglected_tail(Tm::AbstractVector)
    n = length(Tm); n < 3 && return NaN
    tot = sum(Tm); iszero(tot) && return NaN
    ρ = abs(Tm[n]) / max(abs(Tm[n-1]), eps(float(one(eltype(Tm)))))
    ρ = min(ρ, 0.999)
    abs(Tm[n]) * ρ / (1 - ρ) / abs(tot)
end

"""
    run_point(b; K = K_MMAX) -> (; m95, m95L, T, tail, mmax, seconds)

One marker of the figure: m₉₅ from the cumulative per-mode TOTAL torque (74) over m = 1…⌈K/h₀⌉
(`m95L` = the same cutoff on the Lindblad torques (75), registered as a diagnostic).
"""
function run_point(b; K = K_MMAX, grid = GRID)
    t0 = time()
    bg = SW.DiscBackground(; χ = CHI, aEMRI = AEMRI, h0 = H0, 𝔭 = PDISK, 𝔮 = QDISK, b = b)
    mmax = ceil(Int, K / H0)
    surv = SW.torque_survey(bg, 1:mmax; grid_kwargs = grid)
    Tm = SW.mode_torques(surv); TLm = SW.mode_lindblad_torques(surv)
    (; m95 = SW.m_cut(Tm, 0.95), m95L = SW.m_cut(TLm, 0.95), T = surv.T_total, mmax,
       tail = sum(@view Tm[max(1, mmax - 9):mmax]) / surv.T_total, trunc = neglected_tail(Tm),
       seconds = time() - t0)
end

println("\n[fig14] χ = ", CHI, ", a• = ", AEMRI, " M, (𝔭,𝔮) = (", PDISK, ",", QDISK, "), h₀ = ", H0,
        "; b = ", BS, "; m = 1:⌈", K_MMAX, "/h₀⌉ = 1:", ceil(Int, K_MMAX / H0), "; grid = ", GRID)
flush(stdout)

m95s = Int[]; m95Ls = Int[]; tails = Float64[]; truncs = Float64[]
for b in BS
    p = run_point(b)
    push!(m95s, p.m95); push!(m95Ls, p.m95L); push!(tails, p.tail); push!(truncs, p.trunc)
    @printf("[fig14] b=%.2f  m=1:%d  m95=%3d  m95h0=%.3f  (m95 on T_L = %3d)  tail10=%.2e  trunc=%.2e  %.1f s\n",
            b, p.mmax, p.m95, p.m95 * H0, p.m95L, p.tail, p.trunc, p.seconds)
    flush(stdout)
end

# resolution cross-check at the paper's reference triple (240,240,24000)
res_red = Int[]; res_ref = Int[]
for b in (BS[1], BS[end])
    a = run_point(b); r = run_point(b; grid = REF_GRID)
    push!(res_red, a.m95); push!(res_ref, r.m95)
    @printf("[fig14] resolution check b=%.2f: m95 = %d at %s vs %d at %s\n",
            b, a.m95, string(GRID), r.m95, string(REF_GRID)); flush(stdout)
end
resolution_ok = res_red == res_ref

sens = run_point(BS[1]; K = K_SENS)                  # mode-ceiling sensitivity at the largest m₉₅
@printf("[fig14] sensitivity: b=%.2f at m ≤ %d gives m95 = %d (m95h0 = %.3f) vs %d at m ≤ %d\n",
        BS[1], sens.mmax, sens.m95, sens.m95 * H0, m95s[1], ceil(Int, K_MMAX / H0))

m95h0 = m95s .* H0

# ───────────────────────────── registered results ────────────────────────────────────────────────

@repro_result "chi"      CHI
@repro_result "a_emri"   AEMRI
@repro_result "p_disk"   PDISK
@repro_result "q_disk"   QDISK
@repro_result "h0"       H0
@repro_result "b_grid"   collect(BS)
@repro_result "m_max_over_h0" K_MMAX
@repro_result "m95"      m95s
@repro_result "m95_h0"   m95h0
@repro_result "m95_lindblad" m95Ls
@repro_result "tail_fraction_max" maximum(abs, tails)
@repro_result "truncation_estimate_max" maximum(abs, truncs)
@repro_result "grid_nL"   GRID.nL
@repro_result "grid_nC"   GRID.nC
@repro_result "m95_reduced_resolution" res_red
@repro_result "m95_ref_resolution"     res_ref
@repro_result "resolution_check_pass"  resolution_ok
@repro_result "m95_h0_Ksens_b03" sens.m95 * H0
@repro_result "m95_monotone_decreasing" all(diff(m95h0) .<= 0)

# ── the printed-number-gate contract keys (reproduce/arXiv-2608.17003/README.md) ──────────────────
@repro_result "m95h0_b03" m95h0[1]
@repro_result "m95h0_b08" m95h0[end]

# ───────────────────────────── printed-number checks (brief section 5, A:164) ────────────────────
const PAPER_CHECKS_ENFORCED = true

struct PaperCheck; id::String; src::String; paper::Float64; tol::Float64; ours::Float64; end
checks = PaperCheck[
    PaperCheck("m95 h0 at b=0.3 ≈ 2.9",  "A:164", 2.90, 0.25, m95h0[1]),
    PaperCheck("m95 h0 at b=0.8 ≈ 1.95", "A:164", 1.95, 0.25, m95h0[end]),
]
println("\n[fig14] printed-number checks")
println("  claim                                src      paper     ours      tol     verdict")
violations = String[]
for c in checks
    ok = abs(c.ours - c.paper) <= c.tol
    ok || push!(violations, "$(c.id) [$(c.src)]: paper $(c.paper), ours $(round(c.ours; digits=4)), tol ±$(c.tol)")
    @printf("  %-36s %-8s %7.4f  %7.4f  ±%.3f   %s\n", c.id, c.src, c.paper, c.ours, c.tol,
            ok ? "PASS" : "FAIL")
end
let ok = all(diff(m95h0) .<= 0) && m95h0[end] < m95h0[1]
    ok || push!(violations, "m95 h0 decreases with b (A:164): FALSE")
    @printf("  %-36s %-8s %7s  %7s  %6s   %s\n", "m95 h0 decreases with b", "A:164", "true",
            string(ok), "exact", ok ? "PASS" : "FAIL")
end
let tf = maximum(abs, truncs), ok = tf < 5e-2
    ok || push!(violations, "mode-sum truncation: the geometric estimate of the neglected sum is $(tf) > 5e-2")
    @printf("  %-36s %-8s %7s  %7.1e  %6s   %s\n", "neglected mode-sum tail", "M:105", "<5e-2", tf,
            "5e-2", ok ? "PASS" : "FAIL")
end
let ok = resolution_ok
    ok || push!(violations, "resolution check: m95 at $(GRID) = $(res_red) ≠ $(res_ref) at $(REF_GRID)")
    @printf("  %-36s %-8s %7s  %7s  %6s   %s\n", "m95 identical at the reference grid", "A:140",
            string(res_ref), string(res_red), "exact", ok ? "PASS" : "FAIL")
end
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig14] FINDING — printed numbers not reproduced within tolerance:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig14] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written or compared")
        paper_check_failed = true
    end
end

draw(collect(BS), m95h0)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig14", script = @__FILE__, rtol = 1e-8,
                              tolerances = Dict("m95" => (rtol = 0.0, atol = 0.0),
                                                "m95_lindblad" => (rtol = 0.0, atol = 0.0),
                                                "m95_reduced_resolution" => (rtol = 0.0, atol = 0.0),
                                                "m95_ref_resolution" => (rtol = 0.0, atol = 0.0)))
exit(paper_check_failed ? max(code, 1) : code)

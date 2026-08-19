# reproduce/arXiv-2608.17003/fig09.jl — EGD::fig09-num-vs-analytic.
#
# arXiv:2608.17003 Fig. 9 (`fig:num_vs_analytical`, main.tex:965-970): the ratio of the ANALYTIC
# pressureless Lindblad torque to the Lindblad torque of the paper's own linear theory, in
# Schwarzschild spacetime, versus the EMRI semi-major axis a• ∈ [12, 100] M, for (𝔭,𝔮) = (0,0)
# and (0,1) at h0 = 0.05.
#
#   numerator   T_L^analytical = m_max² ⟨Ė'⟩^{(1)}_{Lr,+} / Ω(a•)   — arXiv:2608.17003 (78)
#                                                                     = arXiv:2510.03564 (24),
#               built from Eqs. (18b) + (20) of arXiv:2510.03564 exactly as main.tex:984 prescribes;
#               GrTools: `AccretionDisc.Modes.hkr_lindblad_torque`.
#   denominator T_L^num = Σ_m [F_J(r_out) − F_J(r_in)]_m               — arXiv:2608.17003 (75);
#               GrTools: `AccretionDisc.SpiralWaves.{solve_mode, lindblad_torque}`.
#   cutoff      m_max = m_45, the smallest m at which the NUMERICAL cumulative Lindblad torque
#               reaches 45 % of its converged value (main.tex:988-1000); the paper quotes
#               m_45 ~ 0.35/h0 (80).  GrTools: `SpiralWaves.m_cut(Tm, 0.45)`.
#
# GOLDEN REUSE (checked, does not apply): fig06 also runs (𝔭,𝔮)=(0,0), χ=0, h0=0.05, but on a
# LOG a• grid over [10,1000]M and it registers the TOTAL torque (74), not the per-mode Lindblad
# torques.  m₄₅ needs the per-mode cumulative sum, and a per-mode array for every a• would blow the
# 200 KB golden budget, so fig09 computes its own mode sums.  No other golden of this directory
# carries T_{L,m}.
#
# ZERO physics here: every number comes from those package calls.  Style: reproduce/_style/
# emri_style.jl + the VLM manual section 9 (progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md).

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

const BLUE   = RGB(0x0f / 255, 0x4c / 255, 0x81 / 255)   # (0,0) — manual section 9
const ORANGE = RGB(0xe2 / 255, 0x58 / 255, 0x3e / 255)   # (0,1)

"""
    draw(xs, y00, y01) -> path

Render Fig. 9 from three plotted arrays and save the PNG beside the golden. Split out so that
`--plot-only` can redraw the figure from `expected/fig09.json` after a style change WITHOUT
repeating the multi-hour mode sums (the numbers are the golden's, never recomputed).
"""
function draw(xs, y00, y01)

    # geometry/fonts tuned by the VLM style probe against the 160-dpi render of the reference PDF
    plt = plot(size = (460, 330), dpi = 160, legend = :topright, legendfontsize = 8,
               guidefontsize = 11, tickfontsize = 9,
               xlims = (8, 104), ylims = (0.855, 1.535), xticks = 20:20:100, yticks = 0.9:0.1:1.5,
               xlabel = "\$a_\\bullet\\ (M)\$", ylabel = "\$T_\\mathrm{L}\$ (analytical/num)",
               minorticks = 5, tick_direction = :in, widen = false)
    scatter!(plt, xs, y00; marker = (:diamond, 4.5), color = BLUE, markerstrokewidth = 0,
             label = "\$(\\mathfrak{p},\\mathfrak{q})=(0,0)\$")
    scatter!(plt, xs, y01; marker = (:cross, 6), color = ORANGE, markerstrokewidth = 3.2,
             markerstrokecolor = ORANGE, label = "\$(\\mathfrak{p},\\mathfrak{q})=(0,1)\$")

png_path = joinpath(@__DIR__, "expected", "fig09.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig09] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig09.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end


# `--plot-only` redraws the figure from the committed golden (style iteration without recomputing).
if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig09")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw(v["a_emri"], v["ratio_pq00"], v["ratio_pq01"])
    println("[fig09] --plot-only: redrawn from the golden, no physics recomputed")
    exit(0)
end


using GrTools
using GrTools.AccretionDisc
const SW = GrTools.AccretionDisc.SpiralWaves
const MD = GrTools.AccretionDisc.Modes

# The numerical half needs the mode solver + diagnostics + survey
# (EGD::sw-mode-solver, EGD::sw-diagnostics, EGD::sw-survey-parallel).
for f in (:torque_survey, :mode_lindblad_torques, :m_cut, :torque_units)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: the numerical Lindblad torque of " *
        "arXiv:2608.17003 (75) needs the mode solver and the AMF/torque diagnostics " *
        "(EGD::sw-mode-solver, EGD::sw-diagnostics). The analytic half (Modes.hkr_lindblad_torque) " *
        "is landed; re-run this driver once those nodes are in.")
end

# ───────────────────────────── parameters (main.tex:965-970, brief row 9) ────────────────────────

const CHI   = 0.0                  # Schwarzschild
const H0    = 0.05
const BSOFT = 0.6                  # nominal softening (M:884)
const CONFIGS = ((0, 0), (0, 1))   # (𝔭, 𝔮)
const A_MIN, A_MAX = 12.0, 100.0
const N_A_FULL = 19                # ≈ the point count of the published figure (manual section 9)
# m ceiling of the numerical mode sum. The paper needs ~3/h ≈ 60 modes for ~5 % accuracy
# (main.tex:105 footnote) and m₉₅h₀ ≃ 2.0–3.5 (main.tex:1025) ⇒ m₉₅ ≲ 70 at h0 = 0.05; 120 leaves
# a factor ≈ 2 of headroom. Overridable for the reduced run.
const MMAX_DEFAULT = 120

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
const MMAX    = _envint("GRTOOLS_FIG09_MMAX", MMAX_DEFAULT)
const REDUCED = ("--reduced" in opts.extra) || haskey(ENV, "GRTOOLS_FIG09_REDUCED")
const N_A     = REDUCED ? _envint("GRTOOLS_FIG09_NA", 7) : N_A_FULL
const APS     = collect(range(A_MIN, A_MAX; length = N_A))

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

"""
    torque_pair(𝔭, 𝔮, a) -> (; T_num, T_ana, m45, Tm, seconds)

One point of the figure: the numerical Lindblad torque summed over m = 1…MMAX, the cutoff m₄₅ read
off that same mode sum, and the analytic torque of arXiv:2510.03564 evaluated with m_max = m₄₅.
Both are per q² with ε₀ = 1 and M = 1, i.e. the same units, so the ratio is unit-free.

SIGN CONVENTION. arXiv:2510.03564 Eq. (20) gives ⟨L̇'_z⟩, the secular torque ON THE SECONDARY
(negative for the discs studied here — the EMRI loses angular momentum and migrates inward),
whereas arXiv:2608.17003 Eq. (75) T_L = F_J(r_out) − F_J(r_in) is the angular momentum the spiral
waves carry INTO the disc, i.e. minus the same quantity (positive). The published ratio is
positive (0.86–1.53), so it compares magnitudes: `ratio = −T_ana / T_num`. `signs_opposite` is
registered in the golden so a future convention change is caught rather than silently absorbed.
"""
function torque_pair(𝔭, 𝔮, a)
    t0 = time()
    bg = SW.DiscBackground(; χ = CHI, aEMRI = a, h0 = H0, 𝔭 = 𝔭, 𝔮 = 𝔮, b = BSOFT)
    # the m-loop belongs to the package (threaded, fixed-order reduction, bitwise reproducible)
    surv = SW.torque_survey(bg, 1:MMAX)
    Tm = SW.mode_lindblad_torques(surv)
    T_num = surv.T_L_total
    m45 = SW.m_cut(Tm, 0.45)
    σ = MD.PowerLawSigma(bg.ε0, a, 𝔭)                 # ε = ε0 (r/a•)^{−𝔭}, 2608 (12a) = 2510 (22a)
    T_ana = MD.hkr_lindblad_torque(CHI, a, σ; m_max = m45)
    tail = sum(@view Tm[max(1, MMAX - 9):MMAX]) / T_num     # mode-sum convergence diagnostic
    (; T_num, T_ana, m45, Tm, tail, seconds = time() - t0, budget = surv.budget)
end

println("\n[fig09] a• grid: ", N_A, " points in [", A_MIN, ", ", A_MAX, "] M; m = 1:", MMAX,
        REDUCED ? "  (REDUCED run)" : "")
flush(stdout)

results = Dict{Tuple{Int,Int},NamedTuple}()
for (𝔭, 𝔮) in CONFIGS
    ratios = Float64[]; nums = Float64[]; anas = Float64[]; m45s = Int[]; secs = Float64[]
    tails = Float64[]; opposite = Bool[]
    for a in APS
        p = torque_pair(𝔭, 𝔮, a)
        push!(nums, p.T_num); push!(anas, p.T_ana); push!(m45s, p.m45); push!(secs, p.seconds)
        push!(tails, p.tail); push!(opposite, sign(p.T_ana) == -sign(p.T_num))
        push!(ratios, -p.T_ana / p.T_num)
        @printf("[fig09] (𝔭,𝔮)=(%d,%d)  a•=%7.3f M   T_num=%+.6e  T_ana=%+.6e  ratio=%.4f  m45=%3d  tail=%.2e  %.1f s\n",
                𝔭, 𝔮, a, p.T_num, p.T_ana, -p.T_ana / p.T_num, p.m45, p.tail, p.seconds)
        flush(stdout)   # a multi-hour run must leave a readable trail even if it is interrupted
    end
    results[(𝔭, 𝔮)] = (; ratios, nums, anas, m45s, secs, tails, opposite)
end

r00 = results[(0, 0)]
r01 = results[(0, 1)]

@repro_result "a_emri"            APS
@repro_result "h0"                H0
@repro_result "chi"               CHI
@repro_result "m_max_numerical"   MMAX
@repro_result "ratio_pq00"        r00.ratios
@repro_result "ratio_pq01"        r01.ratios
@repro_result "T_num_pq00"        r00.nums
@repro_result "T_num_pq01"        r01.nums
@repro_result "T_ana_pq00"        r00.anas
@repro_result "T_ana_pq01"        r01.anas
@repro_result "m45_pq00"          r00.m45s
@repro_result "m45_pq01"          r01.m45s
@repro_result "signs_opposite"    all(r00.opposite) && all(r01.opposite)
@repro_result "tail_fraction_max" max(maximum(abs, r00.tails), maximum(abs, r01.tails))

# Mode-sum convergence: share of the Lindblad torque still carried by the top 10 of MMAX modes.
# Reported and REGISTERED, never fatal — an assertion here would throw away a multi-hour run, and
# the number itself belongs in the golden where a reviewer can see it. The threshold below is the
# paper's own accuracy scale: 2608 main.tex:105 states that ~3/h modes give a ~5 % accurate torque,
# so a residual tail of 5e-3 is an order of magnitude tighter than the physics claim it supports.
let tf = max(maximum(abs, r00.tails), maximum(abs, r01.tails))
    @printf("[fig09] mode-sum convergence: top-10 of %d modes carry at most %.2e of the total%s\n",
            MMAX, tf, tf < 5e-3 ? "" : "   ** ABOVE 5e-3: raise GRTOOLS_FIG09_MMAX **")
end

# ── the printed-number-gate contract keys (reproduce/arXiv-2608.17003/README.md) ─────────────────
# m₄₅h₀ is quoted as ~0.35 for (𝔭,𝔮) = (0,0) (main.tex:999); report the median over the a• grid.
_median(v) = (s = sort(collect(v)); n = length(s); isodd(n) ? float(s[(n+1)÷2]) : (s[n÷2] + s[n÷2+1]) / 2)
@repro_result "m45_times_h0"      _median(r00.m45s) * H0
@repro_result "ratio_pq00_a12"    r00.ratios[1]
@repro_result "ratio_pq00_a100"   r00.ratios[end]
@repro_result "ratio_pq01_min"    minimum(r01.ratios)
@repro_result "ratio_pq01_max"    maximum(r01.ratios)

# ───────────────────────────── printed-number checks (brief section 5, M:999/M:1002) ─────────────
# Set to `false` ONLY after a documented investigation shows the paper number is not reproduced;
# the driver then still passes on OUR golden and prints the discrepancy as a FINDING (the node
# status drops to preliminary and the ledger carries the reason).
const PAPER_CHECKS_ENFORCED = true

struct PaperCheck; id::String; src::String; paper::Float64; tol::Float64; ours::Float64; end
checks = PaperCheck[
    PaperCheck("m45 h0 ≈ 0.35",              "M:999",  0.35, 0.10, _median(r00.m45s) * H0),
    PaperCheck("ratio (0,0) at a•=12M ≈ 1.5", "M:1002", 1.50, 0.15, r00.ratios[1]),
    PaperCheck("ratio (0,0) at a•=100M ≈ 0.92","M:1002", 0.92, 0.06, r00.ratios[end]),
    PaperCheck("ratio (0,1) min ≈ 0.88",     "M:1002", 0.88, 0.04, minimum(r01.ratios)),
    PaperCheck("ratio (0,1) max ≈ 0.96",     "M:1002", 0.96, 0.04, maximum(r01.ratios)),
]
println("\n[fig09] printed-number checks")
println("  claim                              src      paper     ours      tol     verdict")
violations = String[]
for c in checks
    ok = abs(c.ours - c.paper) <= c.tol
    ok || push!(violations, "$(c.id) [$(c.src)]: paper $(c.paper), ours $(round(c.ours; digits=4)), tol ±$(c.tol)")
    @printf("  %-34s %-8s %7.4f  %7.4f  ±%.3f   %s\n", c.id, c.src, c.paper, c.ours, c.tol,
            ok ? "PASS" : "FAIL")
end
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig09] FINDING — printed numbers not reproduced within tolerance:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig09] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written")
        println("        or compared, so a multi-hour run is never thrown away on a paper mismatch")
        paper_check_failed = true
    else
        println("[fig09] PAPER_CHECKS_ENFORCED = false ⇒ recorded as a finding; the golden still gates our own numbers")
    end
end

draw(APS, r00.ratios, r01.ratios)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig09", script = @__FILE__,
                             rtol = 1e-8,
                             tolerances = Dict("m45_pq00" => (rtol = 0.0, atol = 0.0),
                                               "m45_pq01" => (rtol = 0.0, atol = 0.0)))
exit(paper_check_failed ? max(code, 1) : code)

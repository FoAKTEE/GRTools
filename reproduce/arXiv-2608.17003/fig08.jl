# reproduce/arXiv-2608.17003/fig08.jl — EGD::fig08-m95-violins.
#
# arXiv:2608.17003 Fig. 8 (`fig:m_95_plot`, main.tex:958-963): violin plots of the torque-cutoff
# parameter m₉₅ scaled by h₀, one group per disc profile (𝔭,𝔮) ∈ {(0,0),(1.5,0),(2,0),(1.5,1),(2,1)}
# and three violins per group for a• ∈ {15, 30, 120}M.  Each violin is the DISTRIBUTION of m₉₅h₀ over
# the scan h₀ ∈ [0.05, 0.15] × χ ∈ {−0.9, 0, +0.9} at the nominal softening b = 0.6.
#
#   m₉₅  = the SMALLEST azimuthal mode number at which the cumulative torque reaches 95 % of the
#          final torque obtained after summing all computed modes (M:990, verbatim); GrTools:
#          `SpiralWaves.m_cut(mode_torques(survey), 0.95)` over the per-mode TOTAL torques (74) —
#          "the total cumulative torque" of the caption (M:958).  The Lindblad-only cutoff is
#          registered alongside as a diagnostic.
#
# Claims gated here: M:1025 (m₉₅h₀ ≃ 2.0–3.5 across the profiles) and M:992 (m₉₅ ∝ 1/h₀, i.e. m₉₅h₀
# approximately universal at fixed (𝔭,𝔮,χ,a•) as h₀ is scanned).
#
# M:1025 CRITERION, CORRECTED IN WAVE 5.  The first version demanded that our band EDGES sit inside
# [1.75,2.25] and [3.25,3.75], i.e. AT the edges of the paper's rounded prose interval "m₉₅h₀ ≃
# 2.0–3.5".  A pixel read-back of the PUBLISHED figure (emri-disk
# results/emri-gasdisk/paper_EMRI-GASDISK/evidence/fig08-m95-violins/REFERENCE_READBACK.md, recorded
# while this survey was still running) measures the published band as [2.208, 3.157]: the published
# figure FAILS that criterion itself (3.157 < 3.25).  The prose interval is a rounded ENVELOPE, not
# a pair of attained extremes.  The check below now mirrors `check_printed_numbers.jl` exactly —
# containment in [2.0,3.5]±0.25, ≥ 55 % coverage of that interval (the published band covers 63.3 %),
# and the ⌈K_SENS/h₀⌉ values within ±0.30 of the published pixels.  See that file's `M1025-m95h0-band`
# comment and reproduce/arXiv-2608.17003/README.md for the full argument and the falsification tests.
#
# ── COMPUTE / RESOLUTION (recorded, verified, not assumed) ───────────────────────────────────────
# This is the heaviest driver of the paper: N_pq × N_a × N_χ × N_h₀ = 5×3×3×N_h₀ mode SURVEYS, each
# summing m = 1…⌈4/h₀⌉ (m₉₅h₀ ≲ 3.5 ⇒ 4/h₀ clears the cutoff).  At the paper's reference grid triple
# (240,240,24000) (A:140) one survey at h₀ = 0.05 costs ≈ 170 s at -t8, i.e. ≳ 13 h for the full
# scan — beyond the mission compute budget.  The surveys are therefore run at the triple
# (60,60,6000), which is one of the paper's OWN convergence-run resolutions (A:139).  That choice is
# VERIFIED INSIDE THIS DRIVER, not assumed: `RES_CHECK` re-runs three spot configurations at
# (240,240,24000) and requires the SAME m₉₅ (a measured pre-check gave m₉₅ = 22 and 58 at both
# resolutions, with the total torque agreeing to 6e-6 relative).
#
# GOLDEN REUSE (checked, does not apply): fig14 stores m₉₅h₀ versus the SOFTENING b at a• = 30M,
# h₀ = 0.05, χ = 0 only; fig09 stores m₄₅ (45 %, Lindblad) at h₀ = 0.05 on an a• grid.  Neither
# carries the (h₀, χ) distribution this figure plots.
#
# ZERO physics here: the KDE/violin geometry below is PRESENTATION ONLY (matplotlib's violinplot is
# reproduced by hand because StatsPlots is not in the pinned reproduce Manifest).  Style:
# reproduce/_style/emri_style.jl + the VLM manual section 8 (progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── parameters (main.tex:958-963, brief row 8) ────────────────────────

const PQ      = ((0.0, 0.0), (1.5, 0.0), (2.0, 0.0), (1.5, 1.0), (2.0, 1.0))
const PQ_KEY  = ("pq00", "pq150", "pq200", "pq151", "pq201")
const PQ_LAB  = ("\$(0,0)\$", "\$(1.5,0)\$", "\$(2,0)\$", "\$(1.5,1)\$", "\$(2,1)\$")
const APS     = (15.0, 30.0, 120.0)
const CHIS    = (-0.9, 0.0, 0.9)
const BSOFT   = 0.6

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
_envflt(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
const N_H0   = _envint("GRTOOLS_FIG08_NH0", 5)          # h₀ samples in [0.05, 0.15]
const H0S    = collect(range(0.05, 0.15; length = N_H0))
const K_MMAX = _envflt("GRTOOLS_FIG08_K", 6.0)          # mode ceiling ⌈K/h₀⌉  (PRIMARY)
# m₉₅ is defined against "the final torque obtained after summing ALL COMPUTED modes" (M:990), so it
# inherits a dependence on that ceiling.  fig14 measured it directly at a• = 30M, b = 0.3, h₀ = 0.05:
# m₉₅ = 54 at ⌈4/h₀⌉ and 58 at ⌈6/h₀⌉ (+7.4 %).  WAVE-6 DECISION D2 (progress/emri-gasdisk/
# DECISIONS_W5.md, resolved in DECISIONS_W6.md): the PRIMARY ceiling is now ⌈6/h₀⌉ — three
# independent measurements (this driver's own K-scan, fig14's b = 0.3 point, fig13's 0.3–1.7 %
# torque deficit) indicate the published sums ran to ≈ 6/h₀.  The pre-registered K = 4 stays
# REGISTERED beside it (`m95_h0_*_Ksens`) at the two configurations that set the band, so the
# ceiling sensitivity remains visible in the golden rather than hidden.
const K_SENS = _envflt("GRTOOLS_FIG08_KSENS", 4.0)
const GRID   = (nL = _envint("GRTOOLS_FIG08_NL", 60), nR = _envint("GRTOOLS_FIG08_NR", 60),
                nC = _envint("GRTOOLS_FIG08_NC", 6000))
const REF_GRID = (nL = 240, nR = 240, nC = 24000)       # the paper's reference triple (A:140)
vkey(ip, ia) = "m95h0_" * PQ_KEY[ip] * "_a" * string(Int(APS[ia]))

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

# manual section 0, palette 5: matplotlib Set2, first three colours
const SET2 = (RGB(0x66 / 255, 0xc2 / 255, 0xa5 / 255),   # a• = 15M   (left violin)
              RGB(0xfc / 255, 0x8d / 255, 0x62 / 255),   # a• = 30M   (middle)
              RGB(0x8d / 255, 0xa0 / 255, 0xcb / 255))   # a• = 120M  (right)
const OFFSETS = (-0.3, 0.0, 0.3)
const VWIDTH  = 0.25

_mean(v) = sum(v) / length(v)
_std(v) = (n = length(v); n < 2 ? 0.0 : sqrt(sum(abs2, v .- _mean(v)) / (n - 1)))

"""
    kde_profile(x; npts = 100) -> (coords, halfwidth)

matplotlib `Axes.violinplot` body geometry reproduced exactly: a 1-D Gaussian KDE with SciPy's
default Scott factor (`n^(-1/5)`, covariance from the unbiased sample variance) evaluated on `npts`
points spanning [min(x), max(x)], rescaled so the widest point of the body equals `VWIDTH`
(half-width `VWIDTH/2` on each side).  Presentation only.
"""
function kde_profile(x::AbstractVector{<:Real}; npts::Int = 100)
    n = length(x)
    lo, hi = float(minimum(x)), float(maximum(x))
    σ = _std(x)
    coords = collect(range(lo, hi; length = npts))
    if σ == 0 || lo == hi                       # degenerate sample: a thin bar at the single value
        return ([lo - 1e-6, lo + 1e-6], [VWIDTH / 2, VWIDTH / 2])
    end
    h = n^(-1 / 5) * σ                          # Scott bandwidth
    vals = [sum(exp(-0.5 * ((c - xi) / h)^2) for xi in x) / (n * h * sqrt(2π)) for c in coords]
    (coords, vals ./ maximum(vals) .* (VWIDTH / 2))
end

"""
    draw(samples) -> path

Render Fig. 8 from the 15 sample vectors (`samples[ip][ia]` = the m₉₅h₀ distribution of group `ip`,
radius `ia`) and save the PNG beside the golden.  `--plot-only` redraws from the golden.
"""
function draw(samples)
    allv = reduce(vcat, [samples[ip][ia] for ip in 1:length(PQ), ia in 1:length(APS)])
    # the reference fixes ylim = (2.00, 3.50) with ticks every 0.25 (manual §8); keep that tick PHASE
    # and only widen the range when our sample falls outside it — data are never clipped.
    ylo, yhi = minimum(allv), maximum(allv)
    ylo = min(2.0, floor(ylo * 4) / 4); yhi = max(3.5, ceil(yhi * 4) / 4)
    plt = plot(size = (470, 350), dpi = 160, legend = :topright, legendfontsize = 8,
               foreground_color_legend = nothing, background_color_legend = nothing,
               guidefontsize = 12, tickfontsize = 10,
               xlims = (-0.6, 4.6), ylims = (ylo, yhi),
               xticks = (0:4, collect(PQ_LAB)), yticks = ylo:0.25:yhi,
               yformatter = y -> @sprintf("%.2f", y),
               xlabel = "\$(\\mathfrak{p},\\mathfrak{q})\$", ylabel = "\$m_{95}h_0\$",
               minorticks = 5, tick_direction = :in, widen = false, bottom_margin = 2Plots.mm)
    # light-grey vertical gridlines at the group centres, BEHIND the violins (manual section 8)
    for c in 0:4
        plot!(plt, [c, c], [ylo, yhi]; color = RGB(0.85, 0.85, 0.85), lw = 0.8, label = "")
    end
    for ip in 1:length(PQ), ia in 1:length(APS)
        x = (ip - 1) + OFFSETS[ia]
        v = samples[ip][ia]
        coords, hw = kde_profile(v)
        xs = vcat(x .- hw, reverse(x .+ hw)); ys = vcat(coords, reverse(coords))
        plot!(plt, Shape(xs, ys); fillcolor = SET2[ia], fillalpha = 0.85, linecolor = :black,
              linewidth = 0.6, label = "")
        lo, hi = minimum(v), maximum(v); mu = _mean(v)
        plot!(plt, [x, x], [lo, hi]; color = :black, lw = 0.8, label = "")           # cbars
        for y in (lo, hi, mu)                                                        # cmins/cmaxes/cmeans
            plot!(plt, [x - VWIDTH / 4, x + VWIDTH / 4], [y, y]; color = :black, lw = 0.8, label = "")
        end
    end
    for ia in 1:length(APS)                                                          # legend patches
        scatter!(plt, [NaN], [NaN]; marker = (:rect, 7), color = SET2[ia], markerstrokewidth = 0,
                 label = "\$a_\\bullet=" * string(Int(APS[ia])) * "M\$")
    end
    annotate!(plt, -0.45, yhi - 0.03 * (yhi - ylo),
              text("distribution over \$h_0,\\chi\$", :left, :top, 9, "Computer Modern"))
    png_path = joinpath(@__DIR__, "expected", "fig08.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig08] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || error("fig08.png is $(n) bytes > 200 KB — reduce the dpi/size")
    png_path
end

if "--plot-only" in opts.extra
    gpath = ReproHarness.expected_path("2608.17003", "fig08")
    isfile(gpath) || ReproHarness.hole("--plot-only needs the golden $(relpath(gpath, ReproHarness.REPRO_ROOT)); run --regen first")
    _, v = ReproHarness.read_expected(gpath)
    draw([[collect(Float64.(v[vkey(ip, ia)])) for ia in 1:length(APS)] for ip in 1:length(PQ)])
    println("[fig08] --plot-only: redrawn from the golden, no physics recomputed")
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
    m95_of(𝔭, 𝔮, ap, χ, h0; grid = GRID) -> (; m95, m95L, T, tail, mmax, seconds)

One sample: solve m = 1…⌈K/h₀⌉ for the disc (𝔭,𝔮,a•,χ,h₀) at b = 0.6 and read the 95 % cutoff off the
cumulative per-mode total torque (74).
"""
function m95_of(𝔭, 𝔮, ap, χ, h0; grid = GRID, K = K_MMAX)
    t0 = time()
    bg = SW.DiscBackground(; χ = χ, aEMRI = ap, h0 = h0, 𝔭 = 𝔭, 𝔮 = 𝔮, b = BSOFT)
    mmax = ceil(Int, K / h0)
    surv = SW.torque_survey(bg, 1:mmax; grid_kwargs = grid)
    Tm = SW.mode_torques(surv); TLm = SW.mode_lindblad_torques(surv)
    (; m95 = SW.m_cut(Tm, 0.95), m95L = SW.m_cut(TLm, 0.95), T = surv.T_total, mmax,
       tail = sum(@view Tm[max(1, mmax - 9):mmax]) / surv.T_total, trunc = neglected_tail(Tm),
       seconds = time() - t0)
end

println("\n[fig08] (𝔭,𝔮) = ", PQ, "; a• = ", APS, " M; χ = ", CHIS, "; b = ", BSOFT)
println("[fig08] h₀ scan (", N_H0, " values) = ", round.(H0S; digits = 4),
        "; m = 1:⌈", K_MMAX, "/h₀⌉; grid = ", GRID)
println("[fig08] ", length(PQ) * length(APS) * length(CHIS) * N_H0, " surveys, ",
        sum(ceil(Int, K_MMAX / h) for h in H0S) * length(PQ) * length(APS) * length(CHIS),
        " mode solves")
flush(stdout)

samples  = [[Float64[] for _ in APS] for _ in PQ]
samplesL = [[Float64[] for _ in APS] for _ in PQ]
cfgs     = NamedTuple[]            # (𝔭,𝔮,a•,χ,h₀) in the order the samples were pushed
vals     = Float64[]
spread   = Float64[]                 # per (𝔭,𝔮,a•,χ): (max−min)/mean of m₉₅h₀ over the h₀ scan
tails    = Float64[]; truncs = Float64[]
t_start  = time()
for ip in 1:length(PQ), ia in 1:length(APS)
    𝔭, 𝔮 = PQ[ip]; ap = APS[ia]
    for χ in CHIS
        row = Float64[]
        for h0 in H0S
            s = m95_of(𝔭, 𝔮, ap, χ, h0)
            push!(samples[ip][ia], s.m95 * h0); push!(samplesL[ip][ia], s.m95L * h0)
            push!(row, s.m95 * h0); push!(tails, s.tail); push!(truncs, s.trunc)
            push!(cfgs, (; 𝔭, 𝔮, ap, χ, h0)); push!(vals, s.m95 * h0)
            @printf("[fig08] (𝔭,𝔮)=(%.1f,%.1f) a•=%5.1f χ=%+.1f h0=%.4f  m=1:%3d  m95=%3d  m95h0=%.3f  tail10=%.1e  %.1f s\n",
                    𝔭, 𝔮, ap, χ, h0, s.mmax, s.m95, s.m95 * h0, s.tail, s.seconds)
            flush(stdout)
        end
        push!(spread, (maximum(row) - minimum(row)) / _mean(row))
    end
    @printf("[fig08] --- group (%.1f,%.1f) a•=%.0f done; elapsed %.1f min\n", 𝔭, 𝔮, ap, (time() - t_start) / 60)
    flush(stdout)
end

# mode-ceiling sensitivity at the two configurations that set the M:1025 band
sens = NamedTuple[]
for i in (argmin(vals), argmax(vals))
    c = cfgs[i]
    a = m95_of(c.𝔭, c.𝔮, c.ap, c.χ, c.h0; K = K_SENS)
    push!(sens, (; c..., m95h0_K = vals[i], m95h0_Ksens = a.m95 * c.h0))
    @printf("[fig08] m_max sensitivity (𝔭,𝔮)=(%.1f,%.1f) a•=%.0f χ=%+.1f h0=%.4f: m95h0 = %.3f at ⌈%.1f/h₀⌉ vs %.3f at ⌈%.1f/h₀⌉\n",
            c.𝔭, c.𝔮, c.ap, c.χ, c.h0, vals[i], K_MMAX, a.m95 * c.h0, K_SENS)
    flush(stdout)
end

# resolution cross-check: the same m₉₅ at the paper's reference triple (240,240,24000)
const RES_CHECK = ((0.0, 0.0, 15.0, 0.0, 0.10), (2.0, 1.0, 15.0, -0.9, 0.05), (1.5, 0.0, 120.0, 0.9, 0.15))
res_red = Int[]; res_ref = Int[]
for (𝔭, 𝔮, ap, χ, h0) in RES_CHECK
    a = m95_of(𝔭, 𝔮, ap, χ, h0)
    b = m95_of(𝔭, 𝔮, ap, χ, h0; grid = REF_GRID)
    push!(res_red, a.m95); push!(res_ref, b.m95)
    @printf("[fig08] resolution check (𝔭,𝔮)=(%.1f,%.1f) a•=%.0f χ=%+.1f h0=%.2f: m95 = %d at %s vs %d at %s\n",
            𝔭, 𝔮, ap, χ, h0, a.m95, string(GRID), b.m95, string(REF_GRID))
    flush(stdout)
end
resolution_ok = res_red == res_ref

allv = reduce(vcat, [samples[ip][ia] for ip in 1:length(PQ), ia in 1:length(APS)])

# ───────────────────────────── registered results ────────────────────────────────────────────────

@repro_result "h0_grid"       H0S
@repro_result "chi_grid"      collect(CHIS)
@repro_result "a_grid"        collect(APS)
@repro_result "b_soft"        BSOFT
@repro_result "m_max_over_h0" K_MMAX
@repro_result "grid_nL"       GRID.nL
@repro_result "grid_nC"       GRID.nC
for ip in 1:length(PQ), ia in 1:length(APS)
    register!(vkey(ip, ia), samples[ip][ia])
end
for ip in 1:length(PQ), ia in 1:length(APS)
    register!("lindblad_" * vkey(ip, ia), samplesL[ip][ia])
end
@repro_result "tail_fraction_max"     maximum(abs, tails)
@repro_result "truncation_estimate_max" maximum(abs, truncs)
@repro_result "m95_ref_resolution"    res_ref
@repro_result "m95_reduced_resolution" res_red
@repro_result "resolution_check_pass" resolution_ok

# ── the printed-number-gate contract keys (reproduce/arXiv-2608.17003/README.md) ──────────────────
@repro_result "m95_h0_min" minimum(allv)
@repro_result "m95_h0_max" maximum(allv)
# M:992 "m₉₅ scales as 1/h₀": at FIXED (𝔭,𝔮,χ,a•) the product m₉₅h₀ must be near-constant as h₀ is
# scanned.  Registered as the WORST relative spread (max−min)/mean over the 45 fixed-parameter scans.
@repro_result "m95_h0_spread_fraction" maximum(spread)
@repro_result "m95_h0_min_Ksens" sens[1].m95h0_Ksens
@repro_result "m95_h0_max_Ksens" sens[2].m95h0_Ksens
@repro_result "m_max_over_h0_sens" K_SENS

@printf("[fig08] mode-sum convergence: top-10 modes carry at most %.2e of the total torque; the geometric estimate of the NEGLECTED sum beyond m ≤ ⌈%.1f/h₀⌉ is at most %.2e\n",
        maximum(abs, tails), K_MMAX, maximum(abs, truncs))

# ───────────────────────────── printed-number checks (brief section 5, M:1025 / M:992) ───────────
const PAPER_CHECKS_ENFORCED = true

println("\n[fig08] printed-number checks")
println("  claim                                src      paper           ours          verdict")
violations = String[]
# The three parts of the corrected M:1025 criterion (identical to check_printed_numbers.jl).
const PUB_BAND = (2.208, 3.157)     # published Fig. 8 pixel read-back (REFERENCE_READBACK.md)
let lo = minimum(allv), hi = maximum(allv)
    loK, hiK = sens[1].m95h0_Ksens, sens[2].m95h0_Ksens
    contained = lo <= hi && 1.75 <= lo <= 3.75 && 1.75 <= hi <= 3.75          # [2.0,3.5] ± 0.25
    cover = max(0.0, min(hi, 3.5) - max(lo, 2.0)) / 1.5                       # published: 0.633
    # (c) compares LIKE WITH LIKE: the published read-back is the band of the published FULL survey,
    # so it is held against OUR full survey at the primary ceiling (K = 6 since wave 6), not against
    # the two-configuration sensitivity re-run.
    pixels = abs(lo - PUB_BAND[1]) <= 0.30 && abs(hi - PUB_BAND[2]) <= 0.30
    ok = contained && cover >= 0.55 && pixels
    ok || push!(violations, "m95h0 band [M:1025]: ours ⌈$(K_MMAX)/h₀⌉ [$(round(lo;digits=3)), $(round(hi;digits=3))] " *
                "(contained=$(contained), coverage=$(round(100cover;digits=1))% of [2.0,3.5], need >=55%; " *
                "vs published pixels [$(PUB_BAND[1]), $(PUB_BAND[2])] within ±0.30 = $(pixels)); " *
                "⌈$(K_SENS)/h₀⌉ sensitivity [$(round(loK;digits=3)), $(round(hiK;digits=3))]")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "m95 h0 band inside 2.0-3.5 ±0.25", "M:1025", "[1.75, 3.75]",
            "[$(round(lo;digits=3)), $(round(hi;digits=3))]", contained ? "PASS" : "FAIL")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "  coverage of [2.0, 3.5]", "M:1025", "≥ 55%",
            "$(round(100cover;digits=1))%", cover >= 0.55 ? "PASS" : "FAIL")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "  ⌈$(K_MMAX)/h₀⌉ band vs published pixels", "M:1025",
            "[$(PUB_BAND[1]), $(PUB_BAND[2])]",
            "[$(round(lo;digits=3)), $(round(hi;digits=3))]", pixels ? "PASS" : "FAIL")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "  ⌈$(K_SENS)/h₀⌉ sensitivity band", "D2",
            "registered", "[$(round(loK;digits=3)), $(round(hiK;digits=3))]", "—")
end
let s = maximum(spread)
    ok = s <= 0.35
    ok || push!(violations, "m95 h0 universality [M:992]: worst relative spread over h0 = $(round(s;digits=4)) > 0.35")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "m95 ∝ 1/h0 (spread over h0)", "M:992", "≤ 0.35",
            string(round(s; digits = 4)), ok ? "PASS" : "FAIL")
end
let tf = maximum(abs, truncs), ok = tf < 5e-2
    ok || push!(violations, "mode-sum truncation: the geometric estimate of the neglected sum is $(tf) > 5e-2")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "neglected mode-sum tail", "M:105", "< 5e-2",
            string(round(tf; sigdigits = 3)), ok ? "PASS" : "FAIL")
end
let ok = resolution_ok
    ok || push!(violations, "resolution check: m95 at $(GRID) = $(res_red) ≠ $(res_ref) at $(REF_GRID)")
    @printf("  %-36s %-8s %-15s %-13s %s\n", "m95 identical at the reference grid", "A:140",
            string(res_ref), string(res_red), ok ? "PASS" : "FAIL")
end
@repro_result "printed_checks_pass" isempty(violations)

paper_check_failed = false
if !isempty(violations)
    println("\n[fig08] FINDING — checks not satisfied:")
    foreach(v -> println("    ", v), violations)
    if PAPER_CHECKS_ENFORCED
        println("[fig08] PAPER_CHECKS_ENFORCED = true ⇒ the driver exits 1 AFTER the golden is written or compared")
        paper_check_failed = true
    end
end

draw(samples)

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig08", script = @__FILE__, rtol = 1e-8,
                              tolerances = Dict("m95_ref_resolution" => (rtol = 0.0, atol = 0.0),
                                                "m95_reduced_resolution" => (rtol = 0.0, atol = 0.0)))
exit(paper_check_failed ? max(code, 1) : code)

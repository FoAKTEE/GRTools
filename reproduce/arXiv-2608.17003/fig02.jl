# reproduce/arXiv-2608.17003/fig02.jl — EGD::fig02-spirals-mmax.
#
# arXiv:2608.17003 Fig. 2 (`fig:spirals_conv`, main.tex:861-866, reference raster
# `ref-paper/arXiv-2608.17003/src/figures/spirals_two_rows.png`, 3226x2001 px): the CONVERGENCE of
# the spiral pattern of a 𝔭 = 2, 𝔮 = 1, h0 = 0.1 disc with the azimuthal-mode ceiling m_max, for an
# EMRI at a• = 15M.
#
#   columns  m_max = 10 | 30 | 120        (the SAME modes, summed to a different ceiling)
#   top      χ = −0.9, full frame |x|, |y| ≤ 40M
#   bottom   χ = +0.9, ZOOM   |x|, |y| ≤ 10M
#
# The caption's claim, and M:893 ("most of the large-scale Fourier ringing subsides by m ≥ 40, though
# some subtle ringing can persist … even to m ≳ 100"), is made quantitative here by the RINGING
# RESIDUAL — the relative L2 distance between the truncated image and the m_max = 120 image over the
# panel's own pixels,
#
#     residual(m) = ‖ δU_{≤m} − δU_{≤120} ‖₂ / ‖ δU_{≤120} ‖₂ ,
#
# measured on a ladder m ∈ {5, 10, 20, 30, 40, 60, 90} for both rows.  The printed-number gate
# (`check_printed_numbers.jl`, claim `M893-ringing`) reads `ringing_residual_m10` and
# `ringing_residual_m30`, which are the LARGE-SCALE (χ = −0.9, |r| ≤ 40) row of that ladder — the
# view M:893 speaks about — and requires residual(30) < residual(10)/3.
#
# ZERO physics here: `AccretionDisc.SpiralWaves.wake_modes` solves m = 1…120 ONCE per row and every
# column is a partial sum of that one cached set of solutions (no mode is solved twice).  Style:
# `_style/emri_style.jl` + `_spiral_panels.jl` beside this file, section 2 of the VLM manual
# `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo).
#
# GOLDEN NOTE (as in fig01.jl): the six maps are far larger than the 200 KB golden budget, so what is
# registered is the reduced content — the ringing ladder, the field extrema, the saturation fraction
# against the ±0.77 colour limits, and the radial cut along the EMRI ray.  `GRTOOLS_FIG02_CACHE=<dir>`
# dumps the raw maps for `--plot-only` style iteration (development aid; never read by --check).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── parameters (main.tex:861-866, brief section 4 row 2) ──────────────

const AEMRI = 15.0
const PDISK = 2.0
const QDISK = 1.0
const H0    = 0.1
const BSOFT = 0.6
const EPS_P = 1e-5

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

const MMAX   = _envint("GRTOOLS_FIG02_MMAX", 120)       # the published converged column
# the three published ceilings and the ringing ladder (clamped to MMAX so the driver can also be
# exercised at a reduced ceiling; at the published MMAX = 120 both lists are exactly as printed)
const MCOLS   = unique([min(m, MMAX) for m in (10, 30, 120)])
const MLADDER = [m for m in (5, 10, 20, 30, 40, 60, 90) if m < MMAX]
const NPIX   = _envint("GRTOOLS_FIG02_NPIX", 601)
# IMAGE MESH (same as fig01.jl): the published rasters are a polar mesh drawn to r = 40M, so the ±40
# top row leaves the corners (r > 40) white while the ±10 zoom is covered everywhere (r ≤ 14.1 < 40).
const MESH_R = 40.0
# App.-C grid resolution: the paper's REFERENCE triple (A:140), registered in the golden
const NL = _envint("GRTOOLS_FIG02_NL", 120)
const NR = _envint("GRTOOLS_FIG02_NR", NL)
const NC = _envint("GRTOOLS_FIG02_NC", 12000)
const NREG   = 96
const NRDIAG = 480                 # radii of the winding curve (log-spaced r_ISCO … frame edge)
const NΨDIAG = 720                 # azimuths of the winding curve (0.5°)
# radial window in which the winding is quoted, per row (the top row is the ±40 frame, the bottom
# row the ±10 zoom): the published raster gives −1.445 / −1.478 / −1.456 turns over 18M → 35M for
# m_max = 10 / 30 / 120 (evidence `reference_arm_measurement.log`), i.e. the winding is a property
# of the DISPERSION RELATION and must not move with the mode ceiling.  This driver asserts that.
const WINDOW = [(18.0, 35.0), (3.0, 9.0)]

# (χ, frame half-width, y ticks, x ticks) per row — the rows do NOT share axes (M:861-866)
const ROWS = [(-0.9, 40.0, [-40, -20, 0, 20, 40], [-20, 0, 20]),
              ( 0.9, 10.0, [-10, -5, 0, 5, 10],   [-5, 0, 5])]

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
include(joinpath(@__DIR__, "_spiral_panels.jl"))
using .SpiralPanels
gr()

_num(x) = isinteger(x) ? string(Int(x)) : string(x)
stamp_of(c) = "\$(\\mathfrak{p},\\mathfrak{q},h_0,\\chi)=($(_num(PDISK)),$(_num(QDISK)),$(_num(H0)),$(_num(c)))\$"

"""
    draw(grids, maps, riscos, rplus) -> path

Render Fig. 2 (2×3 panels with per-ROW frames + one shared RdBu colorbar) and save the PNG beside
the golden.  `grids[row]` is the Cartesian axis of that row, `maps[row][j][…]` its three column maps
(already transposed for `heatmap`).
"""
function draw(grids, maps, riscos, rplus)
    panels = Plots.Plot[]
    for row in 1:2, col in 1:3
        c, lim, yt, xt = ROWS[row]
        g = grids[row]
        push!(panels, wake_panel(g, g, maps[row][col]; lim = lim, r_isco = riscos[row], r_plus = rplus[row],
            xticks = (xt, row == 2 ? string.(xt) : fill("", length(xt))),
            yticks = (yt, col == 1 ? string.(yt) : fill("", length(yt))),
            xlabel = row == 2 ? "\$r\\cos\\phi\$" : "",
            ylabel = col == 1 ? "\$r\\sin\\phi\$" : "",
            title = row == 1 ? "\$m_\\mathrm{max}=$(MCOLS[col])\$" : "",
            stamp = col == 1 ? stamp_of(c) : "", stamp_box = true,
            isco_legend = (row == 1 && col == 1),
            labelsize = 12, ticksize = 9, stampsize = 8))
    end
    cb = colorbar_panel(label = "\$\\delta U\\,M/m_\\bullet\$", ticksize = 9, labelsize = 12)
    plt = plot(panels..., cb; layout = (@layout [grid(2, 3) a{0.075w}]), size = (1240, 790), dpi = 100,
               fontfamily = "Computer Modern")
    for k in 1:6
        row = k <= 3 ? 1 : 2
        col = ((k - 1) % 3) + 1
        plot!(plt[k];
              left_margin   = col == 1 ?  7Plots.mm : -7Plots.mm,
              right_margin  = -4Plots.mm,
              top_margin    = row == 1 ?  1Plots.mm :  0Plots.mm,
              bottom_margin = row == 2 ?  4Plots.mm : -2Plots.mm)
    end
    plot!(plt[7]; left_margin = 2Plots.mm, right_margin = 3Plots.mm,
          top_margin = 1Plots.mm, bottom_margin = 3Plots.mm)
    png_path = joinpath(@__DIR__, "expected", "fig02.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig02] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || @warn "fig02.png is $(n) bytes > 200 KB — reduce the figure size"
    png_path
end

# ── raw-map cache (development aid; never read by --check/--regen) ────────────────────────────────
cache_dir() = get(ENV, "GRTOOLS_FIG02_CACHE", "")
cache_file(row, col) = joinpath(cache_dir(), "fig02_r$(row)c$(col).bin")
function cache_write(row, col, M)
    isempty(cache_dir()) && return
    mkpath(cache_dir())
    open(cache_file(row, col), "w") do io
        write(io, Int64(size(M, 1))); write(io, Int64(size(M, 2))); write(io, M)
    end
end
function cache_read(row, col)
    open(cache_file(row, col), "r") do io
        n1 = read(io, Int64); n2 = read(io, Int64)
        M = Matrix{Float64}(undef, n1, n2); read!(io, M); M
    end
end

axis_of(row) = collect(range(-ROWS[row][2], ROWS[row][2]; length = NPIX))

if "--plot-only" in opts.extra
    isempty(cache_dir()) && ReproHarness.hole("--plot-only needs GRTOOLS_FIG02_CACHE=<dir> from a previous full run")
    using GrTools
    SWp = GrTools.AccretionDisc.SpiralWaves
    riscos = [SWp.r_isco(SWp.DiscBackground(; χ = ROWS[r][1], aEMRI = AEMRI, h0 = H0,
                                            𝔭 = PDISK, 𝔮 = QDISK)) for r in 1:2]
    rplus = [1 + sqrt(1 - ROWS[r][1]^2) for r in 1:2]
    draw([axis_of(r) for r in 1:2],
         [[permutedims(cache_read(r, c)) for c in 1:3] for r in 1:2], riscos, rplus)
    println("[fig02] --plot-only: redrawn from the cached maps, no physics recomputed")
    exit(0)
end

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves
using GrTools.Spacetime: horizons

for f in (:wake_modes, :synthesize, :synthesize_cartesian, :arm_phase)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: Fig. 2 needs the mode solver and the " *
        "image synthesis (EGD::sw-mode-solver, EGD::sw-image-synthesis).")
end

println("\n[fig02] a• = ", AEMRI, " M   (𝔭,𝔮,h0) = (", PDISK, ",", QDISK, ",", H0, ")   m = 1:", MMAX)
println("[fig02] columns m_max = ", MCOLS, "; ringing ladder ", MLADDER, "; ", NPIX, "×", NPIX, " px per panel")
@printf("[fig02] compute budget: %d thread(s), %d CPUs, load %s\n",
        Threads.nthreads(), Sys.CPU_THREADS, string(round.(Sys.loadavg(); digits = 2)))

"Relative L2 distance between two maps over the pixels where BOTH are finite."
function residual(A, B)
    m = isfinite.(A) .& isfinite.(B)
    any(m) || return NaN
    sqrt(sum(abs2, A[m] .- B[m]) / sum(abs2, B[m]))
end

"Solve one row once and reduce it to the three column maps, the ringing ladder and the field stats."
function row_result(row)
    c, lim, _, _ = ROWS[row]
    bg = SW.DiscBackground(; χ = c, aEMRI = AEMRI, h0 = H0, 𝔭 = PDISK, 𝔮 = QDISK, b = BSOFT, εp = EPS_P)
    r_isco = SW.r_isco(bg)
    t0 = time()
    sols = SW.wake_modes(bg, 1:MMAX; verbose = false,          # ONE solve of every mode
                         grid_kwargs = (nL = NL, nR = NR, nC = NC))
    t_solve = time() - t0
    @assert [s.m for s in sols] == collect(1:MMAX)

    xs = axis_of(row)
    t1 = time()
    # every ceiling is a PARTIAL SUM of the same cached solutions
    map_at(mm) = SW.synthesize_cartesian(view(sols, 1:mm), xs, xs; r_lo = r_isco, r_hi = MESH_R)
    ref = map_at(MMAX)
    cols = [mm == MMAX ? ref : map_at(mm) for mm in MCOLS]
    ladder = [residual(map_at(mm), ref) for mm in MLADDER]
    # winding curve of each COLUMN: the crest azimuth ψ(r), from the same cached solutions
    rq_w = exp.(collect(range(log(1.0005 * r_isco), log(lim); length = NRDIAG)))
    Ψq = collect(range(0, 2π; length = NΨDIAG + 1)[1:NΨDIAG])
    psis = [SW.arm_phase(SW.synthesize(view(sols, 1:mm), rq_w, Ψq), Ψq) for mm in MCOLS]
    wa, wb = WINDOW[row]
    iw(x) = argmin(abs.(rq_w .- x))
    turns_col = [(ψ[iw(wb)] - ψ[iw(wa)]) / 2π for ψ in psis]
    t_img = time() - t1

    stats = map(cols) do Z
        v = Z[isfinite.(Z)]
        (dU_max = maximum(v), dU_min = minimum(v), dU_rms = sqrt(sum(abs2, v) / length(v)),
         frac_sat = count(x -> abs(x) > SpiralPanels.CBLIM, v) / length(v))
    end
    # radial cut along the EMRI ray Ψ = 0, out to the panel frame, at the three ceilings
    rq = collect(range(1.0005 * r_isco, lim; length = NREG))
    cut = [SW.synthesize(view(sols, 1:mm), rq, [0.0])[:, 1] for mm in MCOLS]

    decw = round.(Int, range(1, NRDIAG; length = NREG))
    (r_isco = r_isco, r_plus = horizons(c)[1], cols = cols, ladder = ladder, stats = stats,
     rq_w = rq_w[decw], psis = [ψ[decw] for ψ in psis], turns_col = turns_col, window = (wa, wb),
     rq = rq, cut = cut, residual = maximum(s -> s.residual, sols),
     nodes = maximum(s -> length(s.r), sols), t_solve = t_solve, t_img = t_img)
end

rows = Vector{Any}(undef, 2)
for row in 1:2
    rows[row] = row_result(row)
    R = rows[row]
    @printf("[fig02] row %d (χ = %+.1f, |r| ≤ %.0f): %.0f s solve + %.0f s image, max nodes %d, resid %.1e\n",
            row, ROWS[row][1], ROWS[row][2], R.t_solve, R.t_img, R.nodes, R.residual)
    for (j, mm) in enumerate(MCOLS)
        s = R.stats[j]
        @printf("        m_max = %3d: δU ∈ [%.4f, %.4f], rms %.4f, |δU| > %.2f on %.3f%% of the panel\n",
                mm, s.dU_min, s.dU_max, s.dU_rms, SpiralPanels.CBLIM, 100 * s.frac_sat)
        cache_write(row, j, R.cols[j])
    end
    @printf("        winding turns over r ∈ [%.0f, %.0f] M at m_max = %s: %s\n",
            R.window[1], R.window[2], string(MCOLS), join([@sprintf("%+.3f", t) for t in R.turns_col], ", "))
    print("        ringing residual vs m_max = 120:")
    for (j, mm) in enumerate(MLADDER)
        @printf("  %d→%.3f", mm, R.ladder[j])
    end
    println()
    flush(stdout)
end

# ───────────────────────────── registered results ────────────────────────────────────────────────

@repro_result "a_emri"    AEMRI
@repro_result "p_disk"    PDISK
@repro_result "q_disk"    QDISK
@repro_result "h0"        H0
@repro_result "b_soft"    BSOFT
@repro_result "eps_p"     EPS_P
@repro_result "m_max"     MMAX
@repro_result "m_columns" MCOLS
@repro_result "m_ladder"  MLADDER
@repro_result "n_pix"     NPIX
@repro_result "mesh_r_max" MESH_R
@repro_result "grid_nL"   NL
@repro_result "grid_nR"   NR
@repro_result "grid_nC"   NC
@repro_result "solver_residual_max" maximum(R.residual for R in rows)

const ROWTAG = ["chineg", "chipos"]
for row in 1:2
    R = rows[row]
    pre = "row_$(ROWTAG[row])_"
    @repro_result pre * "chi"          ROWS[row][1]
    @repro_result pre * "frame_half_M" ROWS[row][2]
    @repro_result pre * "r_isco"       R.r_isco
    @repro_result pre * "r_plus"       R.r_plus
    @repro_result pre * "ringing_ladder" R.ladder
    @repro_result pre * "cut_r"        R.rq
    @repro_result pre * "arm_r"        R.rq_w
    @repro_result pre * "window_lo"    R.window[1]
    @repro_result pre * "window_hi"    R.window[2]
    for (j, mm) in enumerate(MCOLS)
        s = R.stats[j]
        @repro_result pre * "m$(mm)_dU_max"   s.dU_max
        @repro_result pre * "m$(mm)_dU_min"   s.dU_min
        @repro_result pre * "m$(mm)_dU_rms"   s.dU_rms
        @repro_result pre * "m$(mm)_frac_saturated" s.frac_sat
        @repro_result pre * "m$(mm)_cut_psi0" R.cut[j]
        @repro_result pre * "m$(mm)_arm_phase" R.psis[j]
        @repro_result pre * "m$(mm)_turns"     R.turns_col[j]
    end
end

# the printed-number gate's keys (claim `M893-ringing`): the LARGE-SCALE row, χ = −0.9, |r| ≤ 40
const i10 = findfirst(==(10), MLADDER)
const i30 = findfirst(==(30), MLADDER)
@repro_result "ringing_residual_m10" rows[1].ladder[i10]
@repro_result "ringing_residual_m30" rows[1].ladder[i30]
@repro_result "ringing_residual_m10_zoom" rows[2].ladder[i10]
@repro_result "ringing_residual_m30_zoom" rows[2].ladder[i30]

# ───────────────────────────── the figure's own printed checks ───────────────────────────────────
# M:893 (the gate's `M893-ringing` claim, factor 3) plus monotone convergence of the ladder.

violations = String[]
r10, r30 = rows[1].ladder[i10], rows[1].ladder[i30]
r30 < r10 / 3 ||
    push!(violations, "M:893 large-scale ringing does not subside: residual(30) = $(r30) ≥ residual(10)/3 = $(r10 / 3)")
for row in 1:2
    L = rows[row].ladder
    issorted(L; rev = true) ||
        push!(violations, "row $row: the ringing residual is not monotone in m_max: $(L)")
    L[end] < L[1] / 10 ||
        push!(violations, "row $row: residual(90) = $(L[end]) is not ≪ residual(5) = $(L[1])")
end
# the winding must NOT move with the mode ceiling (the published raster: −1.445/−1.478/−1.456)
for row in 1:2
    tc = rows[row].turns_col
    all(t -> t < 0, tc) ||
        push!(violations, "row $row: the arm does not trail at every m_max (turns = $(tc))")
    maximum(abs.(tc .- tc[end])) <= 0.15 ||
        push!(violations, "row $row: the winding moves with m_max: turns = $(tc) over r ∈ $(rows[row].window)")
end
i40 = findfirst(==(40), MLADDER)
rows[1].ladder[i40] < 0.10 ||
    push!(violations, "M:893 \"ringing subsides by m ≥ 40\": residual(40) = $(rows[1].ladder[i40]) ≥ 0.10")

println("\n[fig02] printed-number checks (M:893)")
@printf("  residual(10) = %.4f, residual(30) = %.4f   →   residual(30) < residual(10)/3 = %.4f   %s\n",
        r10, r30, r10 / 3, r30 < r10 / 3 ? "PASS" : "FAIL")
@printf("  residual(40) = %.4f  (< 0.10, \"ringing subsides by m ≥ 40\")                       %s\n",
        rows[1].ladder[i40], rows[1].ladder[i40] < 0.10 ? "PASS" : "FAIL")
for row in 1:2
    @printf("  row %d ladder monotone decreasing: %s\n", row,
            issorted(rows[row].ladder; rev = true) ? "PASS" : "FAIL")
    @printf("  row %d winding independent of m_max over r ∈ [%.0f, %.0f] M: %s   %s\n", row,
            rows[row].window[1], rows[row].window[2],
            join([@sprintf("%+.3f", t) for t in rows[row].turns_col], ", "),
            maximum(abs.(rows[row].turns_col .- rows[row].turns_col[end])) <= 0.15 ? "PASS" : "FAIL")
end

@repro_result "printed_checks_pass" isempty(violations)

check_failed = false
if !isempty(violations)
    println("\n[fig02] FINDING — claims not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw([axis_of(r) for r in 1:2],
     [[permutedims(rows[r].cols[c]) for c in 1:3] for r in 1:2],
     [R.r_isco for R in rows], [R.r_plus for R in rows])

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig02", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_pix" => (rtol = 0.0, atol = 0.0),
                                                "grid_nL" => (rtol = 0.0, atol = 0.0),
                                                "grid_nR" => (rtol = 0.0, atol = 0.0),
                                                "grid_nC" => (rtol = 0.0, atol = 0.0)))
exit(check_failed ? max(code, 1) : code)

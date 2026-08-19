# reproduce/arXiv-2608.17003/fig01.jl — EGD::fig01-spirals.
#
# arXiv:2608.17003 Fig. 1 (`fig:spirals`, main.tex:853-858, reference raster
# `ref-paper/arXiv-2608.17003/src/figures/spiral_grid.png`, 3276x1949 px): the linear response of six
# accretion discs to an EMRI at a• = 15M, drawn as the enthalpy perturbation δU M/m• on the
# equatorial plane (r cos φ, r sin φ) ∈ [−40, 40]M.
#
#   columns  (𝔭, 𝔮, h0) = (0, 0, 0.15) | (0, 0, 0.1) | (2, 1, 0.1)
#   rows     χ = −0.9 (top, r_ISCO = 8.717M) | χ = +0.9 (bottom, r_ISCO = 2.321M)
#
# The caption makes three falsifiable statements, and this driver measures all three:
#   (i)   "The spiral arms are truncated at the ISCO"          → no field inside r_ISCO;
#   (ii)  "and wind more tightly in thinner disks"             → |turns| larger at h0 = 0.1 than 0.15;
#   (iii) "In disks with ∂_r c_s < 0 [i.e. 𝔮 = 1], the inner spiral arm winds more loosely, while
#          the outer spiral arm winds more tightly"            → |turns_inner| smaller AND
#                                                                |turns_outer| larger at 𝔮 = 1.
# The winding is measured, not eyeballed: `SpiralWaves.arm_phase` returns the unwrapped azimuth of
# the crest as a function of radius, and `turns = Δψ/2π` over a fixed radial window.
#
# ZERO physics here.  The disc, the per-m solutions, the mode-sum image and the winding curve all
# come from `AccretionDisc.SpiralWaves` (EGD::sw-background/grid/master/source/mode-solver/
# image-synthesis).  Style: `_style/emri_style.jl` + `_spiral_panels.jl` beside this file, following
# sections 0-1 of the VLM manual `progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo).
#
# GOLDEN NOTE.  The six 601×601 δU maps are ~17 MB of Float64 and cannot live in a ≤ 200 KB golden.
# What is registered instead is the REDUCED, physics-bearing content of each panel: the winding curve
# ψ(r) (the quantity the caption is about), the radial cut along the EMRI ray Ψ = 0, the field
# extrema and their location, and the saturation fraction against the published ±0.77 colour limits.
# The PNG beside the golden is the VLM artefact; `GRTOOLS_FIG01_CACHE=<dir>` additionally dumps the
# raw maps so `--plot-only` can iterate the STYLE without repeating the mode sum (development aid
# only — `--check` and `--regen` never read it).

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness
using Printf

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :heavy)

# ───────────────────────────── parameters (main.tex:853-858, brief section 4 row 1) ──────────────

const AEMRI = 15.0                 # a• = 15M
const BSOFT = 0.6                  # nominal softening (M:884)
const EPS_P = 1e-5                 # nominal pole displacement (M:735)

_envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

# m ceiling of the mode sum.  M:893: "most of the large-scale Fourier ringing subsides by m ≥ 40,
# though some subtle ringing can persist … even to m ≳ 100"; Fig. 2 shows m_max = 120 as the
# converged column, so 120 is the published choice for these images.
const MMAX = _envint("GRTOOLS_FIG01_MMAX", 120)

# App.-C grid resolution.  DEFAULT = (120, 120, 12000), one step below the paper's torque-convergence
# REFERENCE triple (240, 240, 24000) of A:140.  This is a MEASURED choice, not a guess: the probe in
# `results/emri-gasdisk/paper_EMRI-GASDISK/evidence/fig01-spirals/grid_resolution_probe.log` (emri-disk)
# synthesises the hardest of these panels at four triples and finds the (120, 12000) image to differ
# from the reference-triple image by L2-rel 5.17e-4 and at most 1.07e-3 of the peak — i.e. 0.08 % of
# the published ±0.77 colour scale, two orders of magnitude below the ~5 % with which an amplitude can
# be read off the published raster at all.  The triple is registered in the golden (grid_nL/nR/nC), and
# the env knobs below let the probe (and a re-check at the reference triple) use the SAME driver.
const NL = _envint("GRTOOLS_FIG01_NL", 120)
const NR = _envint("GRTOOLS_FIG01_NR", NL)
const NC = _envint("GRTOOLS_FIG01_NC", 12000)

const FRAME = 40.0                 # the published half-width of every panel, in M
# IMAGE MESH.  The published rasters are a POLAR mesh drawn out to r = 40M: the ±40 panels show a
# disc that touches the frame at the four edge midpoints and leaves the CORNERS (r > 40) white, and
# nothing is drawn inside r_ISCO.  The mode solutions live on the fiducial domain (C1) out to
# 5a• = 75M; only the drawn — and therefore measured — region is r ∈ [r_ISCO, 40M].
const MESH_R = FRAME
const NPIX  = _envint("GRTOOLS_FIG01_NPIX", 601)     # Cartesian pixels per axis (0.133 M/px)
const NRDIAG = 480                 # radii of the winding curve (log-spaced r_ISCO … mesh edge)
const NΨDIAG = 720                 # azimuths of the winding curve (0.5°)
const NREG   = 96                  # decimation of the registered curves

# (𝔭, 𝔮, h0, χ) of the six panels, in reading order
const CASES = [(0.0, 0.0, 0.15, -0.9), (0.0, 0.0, 0.10, -0.9), (2.0, 1.0, 0.10, -0.9),
               (0.0, 0.0, 0.15,  0.9), (0.0, 0.0, 0.10,  0.9), (2.0, 1.0, 0.10,  0.9)]

# radial windows of the winding measurement, in units of a• (the inner one starts above the ISCO so
# the crest is tracked in the wave zone, not in the truncation layer)
const WIN_OUT = (1.20, FRAME / AEMRI)
win_in(r_isco) = (max(1.15 * r_isco / AEMRI, 0.20), 0.85)

# ───────────────────────────── the figure (presentation only) ────────────────────────────────────

using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
include(joinpath(@__DIR__, "_spiral_panels.jl"))
using .SpiralPanels
gr()

_num(x) = isinteger(x) ? string(Int(x)) : string(x)
stamp_of(p, q, h, c) = "\$(\\mathfrak{p},\\mathfrak{q},h_0,\\chi)=($(_num(p)),$(_num(q)),$(_num(h)),$(_num(c)))\$"

"""
    draw(xs, ys, maps, riscos, rplus) -> path

Render Fig. 1 (2×3 panels + one shared RdBu colorbar) and save the PNG beside the golden.
`maps[k][j, i]` is δU M/m• at (xs[i], ys[j]) of panel k, already transposed for `heatmap`.
"""
function draw(xs, ys, maps, riscos, rplus)
    panels = Plots.Plot[]
    for (k, (p, q, h, c)) in enumerate(CASES)
        row = k <= 3 ? 1 : 2
        col = ((k - 1) % 3) + 1
        push!(panels, wake_panel(xs, ys, maps[k]; lim = FRAME, r_isco = riscos[k], r_plus = rplus[k],
            xticks = ([-20, 0, 20], row == 2 ? ["-20", "0", "20"] : ["", "", ""]),
            yticks = ([-20, 0, 20], col == 1 ? ["-20", "0", "20"] : ["", "", ""]),
            xlabel = row == 2 ? "\$r\\cos(\\phi)\$" : "",
            ylabel = col == 1 ? "\$r\\sin(\\phi)\$" : "",
            stamp = stamp_of(p, q, h, c), isco_legend = (k == 1),
            labelsize = 12, ticksize = 9, stampsize = 8))
    end
    cb = colorbar_panel(label = "\$\\delta U\\,\\frac{M}{m_\\bullet}\$", ticksize = 9, labelsize = 12)
    plt = plot(panels..., cb; layout = (@layout [grid(2, 3) a{0.075w}]), size = (1240, 790), dpi = 100,
               fontfamily = "Computer Modern")
    # tight panel packing (manual section 1: wspace/hspace ~ 0.04) — per-subplot margins, since the
    # empty tick labels of the interior panels would otherwise reserve their full width/height
    for k in 1:6
        row = k <= 3 ? 1 : 2
        col = ((k - 1) % 3) + 1
        plot!(plt[k];
              left_margin  = col == 1 ?  6Plots.mm : -7Plots.mm,
              right_margin = -4Plots.mm,
              top_margin   = row == 1 ?  1Plots.mm : -7Plots.mm,
              bottom_margin = row == 2 ? 4Plots.mm : -4Plots.mm)
    end
    plot!(plt[7]; left_margin = 2Plots.mm, right_margin = 3Plots.mm,
          top_margin = 1Plots.mm, bottom_margin = 3Plots.mm)
    png_path = joinpath(@__DIR__, "expected", "fig01.png")
    mkpath(dirname(png_path))
    savefig(plt, png_path)
    n = filesize(png_path)
    println("[fig01] wrote ", relpath(png_path, ReproHarness.REPRO_ROOT), " (", n, " bytes)")
    n <= 200 * 1024 || @warn "fig01.png is $(n) bytes > 200 KB — reduce the figure size or NPIX"
    png_path
end

# ── raw-map cache (development aid for style iteration; never read by --check/--regen) ────────────
cache_dir() = get(ENV, "GRTOOLS_FIG01_CACHE", "")
cache_file(k) = joinpath(cache_dir(), "fig01_panel$(k).bin")
function cache_write(k, M)
    isempty(cache_dir()) && return
    mkpath(cache_dir())
    open(cache_file(k), "w") do io
        write(io, Int64(size(M, 1))); write(io, Int64(size(M, 2))); write(io, M)
    end
end
function cache_read(k)
    open(cache_file(k), "r") do io
        n1 = read(io, Int64); n2 = read(io, Int64)
        M = Matrix{Float64}(undef, n1, n2); read!(io, M); M
    end
end

const XS = collect(range(-FRAME, FRAME; length = NPIX))
const YS = copy(XS)

if "--plot-only" in opts.extra
    isempty(cache_dir()) && ReproHarness.hole("--plot-only needs GRTOOLS_FIG01_CACHE=<dir> from a previous full run")
    using GrTools
    SWp = GrTools.AccretionDisc.SpiralWaves
    riscos = [SWp.r_isco(SWp.DiscBackground(; χ = c, aEMRI = AEMRI, h0 = h, 𝔭 = p, 𝔮 = q))
              for (p, q, h, c) in CASES]
    rplus = [1 + sqrt(1 - c^2) for (_, _, _, c) in CASES]
    draw(XS, YS, [cache_read(k) for k in 1:6], riscos, rplus)
    println("[fig01] --plot-only: redrawn from the cached maps, no physics recomputed")
    exit(0)
end

# ───────────────────────────── the computation (package calls only) ──────────────────────────────

using GrTools
const SW = GrTools.AccretionDisc.SpiralWaves
using GrTools.Spacetime: horizons

for f in (:wake_modes, :synthesize, :synthesize_cartesian, :arm_phase)
    isdefined(SW, f) || ReproHarness.hole(
        "AccretionDisc.SpiralWaves.$(f) is not available: Fig. 1 needs the mode solver and the " *
        "image synthesis (EGD::sw-mode-solver, EGD::sw-image-synthesis).")
end

println("\n[fig01] a• = ", AEMRI, " M   m = 1:", MMAX, "   b = ", BSOFT, "   ε_p = ", EPS_P)
println("[fig01] image ", NPIX, "×", NPIX, " px over [-", FRAME, ", ", FRAME, "]M, polar mesh r ≤ ", MESH_R, "M; winding curve on ",
        NRDIAG, " radii × ", NΨDIAG, " azimuths; App.-C grid (n_L,n_R,n_C) = (", NL, ",", NR, ",", NC, ")")
@printf("[fig01] compute budget: %d thread(s), %d CPUs, load %s\n",
        Threads.nthreads(), Sys.CPU_THREADS, string(round.(Sys.loadavg(); digits = 2)))

"Reduce one panel to the numbers the golden keeps (and the map the figure draws)."
function panel(p, q, h, c)
    bg = SW.DiscBackground(; χ = c, aEMRI = AEMRI, h0 = h, 𝔭 = p, 𝔮 = q, b = BSOFT, εp = EPS_P)
    r_isco = SW.r_isco(bg)
    t0 = time()
    sols = SW.wake_modes(bg, 1:MMAX; verbose = false,
                         grid_kwargs = (nL = NL, nR = NR, nC = NC))
    t_solve = time() - t0
    resid = maximum(s -> s.residual, sols)
    nodes = maximum(s -> length(s.r), sols)

    t1 = time()
    Z = SW.synthesize_cartesian(sols, XS, YS; r_lo = r_isco, r_hi = MESH_R)        # (nx, ny)
    # winding curve on a log-spaced polar grid (the inner windings of the χ = +0.9 row are tight)
    rq = exp.(collect(range(log(1.0005 * r_isco), log(MESH_R); length = NRDIAG)))
    Ψq = collect(range(0, 2π; length = NΨDIAG + 1)[1:NΨDIAG])
    Fp = SW.synthesize(sols, rq, Ψq)
    ψ = SW.arm_phase(Fp, Ψq)
    t_img = time() - t1

    turns(lo, hi) = begin
        i1 = argmin(abs.(rq .- lo * AEMRI)); i2 = argmin(abs.(rq .- hi * AEMRI))
        (ψ[i2] - ψ[i1]) / 2π
    end
    wi = win_in(r_isco)
    t_out = turns(WIN_OUT...)
    t_in  = turns(wi...)

    fin = isfinite.(Z)
    vals = Z[fin]
    dU_max, dU_min = maximum(vals), minimum(vals)
    dU_rms = sqrt(sum(abs2, vals) / length(vals))
    frac_sat = count(v -> abs(v) > SpiralPanels.CBLIM, vals) / length(vals)
    # the ISCO truncation: nothing is drawn inside r_ISCO
    n_inside = count(k -> isfinite(Z[k]) && hypot(XS[k[1]], YS[k[2]]) < r_isco, CartesianIndices(Z))

    # WAKE peak: the strongest response in the CO-ORBITAL annulus |r − a•| ≤ 0.2a•.  (The GLOBAL
    # maximum is registered too but is not at the EMRI for χ = +0.9, where the tight windings just
    # outside the 2.32M ISCO are stronger than the wake — visible in the published bottom row.)
    coorb = [k for k in CartesianIndices(Z) if isfinite(Z[k]) &&
             abs(hypot(XS[k[1]], YS[k[2]]) - AEMRI) <= 0.2 * AEMRI]
    kw = coorb[argmax([abs(Z[k]) for k in coorb])]
    kg = findall(fin)[argmax(abs.(vals))]

    # radial structure on the polar diagnostic grid: the cut along the EMRI ray Ψ = 0, and the
    # azimuthal amplitude √⟨δU²⟩_Ψ(r) whose 1 %-of-peak crossing is the VISIBLE disc edge
    # the VISIBLE disc edge: the outermost radius whose azimuthal amplitude still shows on the
    # published ±0.77 colour scale, taken as 2 % of that scale (≈ 1/8 of the smallest labelled tick
    # 0.154); it saturates at the mesh edge when the wake is still visible there, as in the paper.
    amp = [sqrt(sum(abs2, view(Fp, i, :)) / NΨDIAG) for i in 1:NRDIAG]
    iedge = findlast(a -> a > 0.02 * SpiralPanels.CBLIM, amp)
    r_edge = iedge === nothing ? NaN : rq[iedge]
    dec = round.(Int, range(1, NRDIAG; length = NREG))
    cut0 = [Fp[i, 1] for i in dec]

    (r_isco = r_isco, r_plus = horizons(c)[1], Z = Z, rq = rq[dec], psi = ψ[dec], cut0 = cut0,
     amp = amp[dec], r_edge = r_edge,
     turns_out = t_out, turns_in = t_in, win_in = wi,
     dU_max = dU_max, dU_min = dU_min, dU_rms = dU_rms, frac_sat = frac_sat,
     peak_x = XS[kw[1]], peak_y = YS[kw[2]], peak_val = Z[kw],
     gpeak_x = XS[kg[1]], gpeak_y = YS[kg[2]], gpeak_val = Z[kg],
     n_inside = n_inside, residual = resid, nodes = nodes,
     seconds = t_solve + t_img, t_solve = t_solve, t_img = t_img)
end

panels = Vector{Any}(undef, 6)
for (k, (p, q, h, c)) in enumerate(CASES)
    panels[k] = panel(p, q, h, c)
    P = panels[k]
    @printf("[fig01] panel %d (𝔭,𝔮,h0,χ)=(%s,%s,%s,%s): %.0f s (%.0f solve + %.0f image), max nodes %d, resid %.1e\n",
            k, _num(p), _num(q), _num(h), _num(c), P.seconds, P.t_solve, P.t_img, P.nodes, P.residual)
    @printf("        δU ∈ [%.4f, %.4f], rms %.4f, |δU|>%.2f on %.3f%% of the disc; visible edge r = %.1f M\n",
            P.dU_min, P.dU_max, P.dU_rms, SpiralPanels.CBLIM, 100 * P.frac_sat, P.r_edge)
    @printf("        co-orbital wake peak %+.4f at (%.2f, %.2f); global peak %+.4f at (%.2f, %.2f)\n",
            P.peak_val, P.peak_x, P.peak_y, P.gpeak_val, P.gpeak_x, P.gpeak_y)
    @printf("        winding: %.3f turns over r/a• ∈ [%.2f, %.2f] (outer), %.3f turns over [%.2f, %.2f] (inner)\n",
            P.turns_out, WIN_OUT[1], WIN_OUT[2], P.turns_in, P.win_in[1], P.win_in[2])
    cache_write(k, P.Z)
    flush(stdout)
end

# ───────────────────────────── registered results ────────────────────────────────────────────────

@repro_result "a_emri"       AEMRI
@repro_result "m_max"        MMAX
@repro_result "b_soft"       BSOFT
@repro_result "eps_p"        EPS_P
@repro_result "n_pix"        NPIX
@repro_result "grid_nL"      NL
@repro_result "grid_nR"      NR
@repro_result "grid_nC"      NC
@repro_result "frame_half_M" FRAME
@repro_result "mesh_r_max"   MESH_R
@repro_result "cb_limit"     SpiralPanels.CBLIM
@repro_result "arm_r_top"    panels[1].rq
@repro_result "arm_r_bottom" panels[4].rq
@repro_result "solver_residual_max" maximum(P.residual for P in panels)

for (k, (p, q, h, c)) in enumerate(CASES)
    P = panels[k]
    pre = "p$(k)_"
    @repro_result pre * "chi"          c
    @repro_result pre * "p_disk"       p
    @repro_result pre * "q_disk"       q
    @repro_result pre * "h0"           h
    @repro_result pre * "r_isco"       P.r_isco
    @repro_result pre * "r_plus"       P.r_plus
    @repro_result pre * "dU_max"       P.dU_max
    @repro_result pre * "dU_min"       P.dU_min
    @repro_result pre * "dU_rms"       P.dU_rms
    @repro_result pre * "frac_saturated" P.frac_sat
    @repro_result pre * "wake_peak_x"     P.peak_x
    @repro_result pre * "wake_peak_y"     P.peak_y
    @repro_result pre * "wake_peak_value" P.peak_val
    @repro_result pre * "global_peak_x"     P.gpeak_x
    @repro_result pre * "global_peak_y"     P.gpeak_y
    @repro_result pre * "global_peak_value" P.gpeak_val
    @repro_result pre * "r_visible_edge"  P.r_edge
    @repro_result pre * "turns_outer"  P.turns_out
    @repro_result pre * "turns_inner"  P.turns_in
    @repro_result pre * "arm_phase"    P.psi
    @repro_result pre * "cut_psi0"     P.cut0
    @repro_result pre * "amp_rms_psi"  P.amp
    @repro_result pre * "n_inside_isco" P.n_inside
end

# ───────────────────────────── the figure's own checks (caption claims) ──────────────────────────
# Fig. 1 prints no number; its CLAIMS are the three caption statements plus the wake location, and
# every one of them is asserted here.  A violation prints a FINDING and exits 1.

violations = String[]
for (k, (p, q, h, c)) in enumerate(CASES)
    P = panels[k]
    hypot(P.peak_x - AEMRI, P.peak_y) <= 2.0 ||
        push!(violations, "panel $k: the co-orbital |δU| peak sits at ($(P.peak_x), $(P.peak_y)), not at the EMRI ($(AEMRI), 0)")
    P.n_inside == 0 ||
        push!(violations, "panel $k: $(P.n_inside) pixels inside r_ISCO carry field (the arm must be truncated there)")
    P.turns_out < 0 ||
        push!(violations, "panel $k: the outer arm is not trailing (turns_outer = $(P.turns_out) ≥ 0)")
    P.turns_in < 0 ||
        push!(violations, "panel $k: the inner arm is not trailing (turns_inner = $(P.turns_in) ≥ 0)")
    P.dU_max >= 0.154 ||
        push!(violations, "panel $k: max δU = $(P.dU_max) never reaches the first labelled colorbar tick (0.154)")
end
# (ii) thinner disc ⇒ tighter winding: column 2 (h0 = 0.10) vs column 1 (h0 = 0.15), both rows
for (a, b, row) in ((1, 2, "χ = −0.9"), (4, 5, "χ = +0.9"))
    abs(panels[b].turns_out) > abs(panels[a].turns_out) ||
        push!(violations, "$row: h0 = 0.10 does not wind more tightly than h0 = 0.15 " *
                          "(|turns_outer| $(abs(panels[b].turns_out)) ≤ $(abs(panels[a].turns_out)))")
end
# (iii) ∂_r c_s < 0 (𝔮 = 1): inner arm LOOSER, outer arm TIGHTER — column 3 vs column 2, both rows
for (a, b, row) in ((2, 3, "χ = −0.9"), (5, 6, "χ = +0.9"))
    abs(panels[b].turns_in) < abs(panels[a].turns_in) ||
        push!(violations, "$row: 𝔮 = 1 inner arm is not looser (|turns_inner| $(abs(panels[b].turns_in)) ≥ $(abs(panels[a].turns_in)))")
    abs(panels[b].turns_out) > abs(panels[a].turns_out) ||
        push!(violations, "$row: 𝔮 = 1 outer arm is not tighter (|turns_outer| $(abs(panels[b].turns_out)) ≤ $(abs(panels[a].turns_out)))")
end

println("\n[fig01] caption checks")
for (k, _) in enumerate(CASES)
    P = panels[k]
    @printf("  panel %d  wake peak at (%6.2f,%6.2f)  inside-ISCO %3d  turns out %+7.3f  in %+7.3f  max δU %6.3f  edge %5.1f\n",
            k, P.peak_x, P.peak_y, P.n_inside, P.turns_out, P.turns_in, P.dU_max, P.r_edge)
end
@printf("  thinner ⇒ tighter (χ=-0.9): |%.3f| > |%.3f|   %s\n", panels[2].turns_out, panels[1].turns_out,
        abs(panels[2].turns_out) > abs(panels[1].turns_out) ? "PASS" : "FAIL")
@printf("  thinner ⇒ tighter (χ=+0.9): |%.3f| > |%.3f|   %s\n", panels[5].turns_out, panels[4].turns_out,
        abs(panels[5].turns_out) > abs(panels[4].turns_out) ? "PASS" : "FAIL")
@printf("  𝔮=1 inner looser  (χ=-0.9): |%.3f| < |%.3f|   %s\n", panels[3].turns_in, panels[2].turns_in,
        abs(panels[3].turns_in) < abs(panels[2].turns_in) ? "PASS" : "FAIL")
@printf("  𝔮=1 outer tighter (χ=-0.9): |%.3f| > |%.3f|   %s\n", panels[3].turns_out, panels[2].turns_out,
        abs(panels[3].turns_out) > abs(panels[2].turns_out) ? "PASS" : "FAIL")
@printf("  𝔮=1 inner looser  (χ=+0.9): |%.3f| < |%.3f|   %s\n", panels[6].turns_in, panels[5].turns_in,
        abs(panels[6].turns_in) < abs(panels[5].turns_in) ? "PASS" : "FAIL")
@printf("  𝔮=1 outer tighter (χ=+0.9): |%.3f| > |%.3f|   %s\n", panels[6].turns_out, panels[5].turns_out,
        abs(panels[6].turns_out) > abs(panels[5].turns_out) ? "PASS" : "FAIL")

@repro_result "caption_checks_pass" isempty(violations)

check_failed = false
if !isempty(violations)
    println("\n[fig01] FINDING — caption claims not reproduced:")
    foreach(v -> println("    ", v), violations)
    check_failed = true
end

draw(XS, YS, [permutedims(P.Z) for P in panels],
     [P.r_isco for P in panels], [P.r_plus for P in panels])

code = ReproHarness.conclude!(opts; paper = "2608.17003", id = "fig01", script = @__FILE__,
                              rtol = 1e-8,
                              tolerances = Dict("m_max" => (rtol = 0.0, atol = 0.0),
                                                "n_pix" => (rtol = 0.0, atol = 0.0),
                                                "grid_nL" => (rtol = 0.0, atol = 0.0),
                                                "grid_nR" => (rtol = 0.0, atol = 0.0),
                                                "grid_nC" => (rtol = 0.0, atol = 0.0),
                                                "p1_n_inside_isco" => (rtol = 0.0, atol = 0.0)))
exit(check_failed ? max(code, 1) : code)

# reproduce/arXiv-2608.17003/_published_grids.jl — the PUBLISHED abscissa grid of Figs. 6 and 7.
#
# Figs. 6 (`chi_rp.pdf`) and 7 (`rp_h0.pdf`) share one a•-grid: fig07's h₀ = 0.05 series of the (0,0)
# panel must stay BITWISE identical to fig06's χ = 0 series, so the grid lives here and both drivers
# include it — it is never written twice.
#
# WHY IT IS NOT `10 .^ range(1, 3; length = 40)`.  That was this directory's original reading (from a
# marker-pitch estimate of 13.565 px on a 264.5 px/decade axis).  It is WRONG at the small-radius end
# and it produced a spurious "5–15 % strong-field deficit": our a• = 10 M value was being compared
# with the published leftmost marker, which is at a• = 11 M.  The measurement that corrected it is
# `results/emri-gasdisk/paper_EMRI-GASDISK/evidence/sw-grid-gr-refinement/{REPORT.md,ref_markers.py,
# ref_fig06_markers.log}` (emri-disk repo), re-confirmed here at 400 dpi:
#
#   * calibration from the axis's own decade ticks — 514.0 px/decade at 400 dpi, a• = 10 at col 177.5,
#     the eight minor ticks per decade fitting to ≤ 0.38 px;
#   * the three ISOLATED markers of the steep χ = −0.9 series sit at a• = 10.986, 12.016, 13.025 M,
#     a pitch of 10.0 and 9.0 px where Δa• = 1 M predicts 10.16 and 9.29 px and a 40-point log grid
#     would predict a constant 13.18 px  ⇒  the small-radius abscissae are the INTEGERS 11, 12, 13 M;
#   * the twelve isolated markers of the χ = 0 series between a• = 248 M and 941 M are geometric with
#     Δlog₁₀ a• = 0.05270 ± 0.0003 (excluding 2/39 = 0.05128 at 6σ), and the last marker centre —
#     half-clipped by the right spine — sits at a• ≈ 1063 M, not 1000 M;
#   * no marker exists between 13.03 M and ≈ 14.6 M (the 13.03 M marker is isolated, and two 13-px
#     markers 9 px apart would have merged), so the integer run stops at 13.
#
# [PRELIMINARY] The exact published CONSTRUCTION is not recoverable at raster precision: 40 points,
# an integer low end and a geometric tail pinned at ≈ 1063 M over-determine nothing between 13 M and
# 15 M.  The grid below reproduces the three DIRECTLY MEASURED small-radius abscissae EXACTLY (that is
# where the curve is steep — dlnT/dln a• ≈ 2.8 at χ = −0.9 — and where the spurious deficit lived),
# reproduces the measured last abscissa, and matches the measured tail ratio to 2.4 %.  Beyond
# a• ≳ 16 M the reproduction already agrees with the published markers to 1–3 % and T(a•) is flat, so
# the residual tail ambiguity is immaterial there.
const PUBLISHED_AP_GRID =
    vcat([11.0, 12.0, 13.0], 10.0 .^ range(log10(15.0), log10(1063.4); length = 37))

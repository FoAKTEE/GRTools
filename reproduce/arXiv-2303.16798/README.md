# arXiv:2303.16798 — linear-in-spin fluxes from generic Kerr orbits

Reproduced with GrTools `Perturbation.SpinFluxes` (modes, amplitudes, fluxes) on top of
`Spacetime.SpinningSecondary` (the linear-in-spin orbit-correction engine).
Legacy sources: `emri-disk/results/gw-flux/figures/2303.16798/repro_*.jl` and the published-table
testset `emri-disk/packages/TeukolskySpinFluxes/test/runtests.jl:731-760` (both read-only).

| published object | driver | result id | tier | status | wall (-t4) |
|---|---|---|---|---|---|
| Table I (`tab:error_generic`) row 1 | `check_tableI.jl` | `tableI` | fast | solid | 117 s |
| Table I, all six rows | `check_tableI_full.jl` | `tableI_full` | heavy | **future** (no golden) | est. 20–25 min |
| Table II (`tab:spherical`), 12 rows | `check_tableII.jl` | `tableII` | fast | solid | 88 s |
| Table III (`tab:generic`), 24 rows | `check_tableIII.jl` | `tableIII` | fast | solid | 94 s |
| `fig:fourier_coefficients` | `fig_fourier_coefficients.jl` | `fig_fourier_coefficients` | fast | solid | 47 s |
| `fig:fluxes_spherical` (FSl{1,2,3}0k) | `fig_FSl0k_spherical.jl` | `fig_FSl0k_spherical` | fast | solid | 82 s |
| `fig:comparison_spherical` | `fig_comparison_spherical.jl` | `fig_comparison_spherical` | heavy | **future** (no golden) | 79 min, not completed |
| `fig:comparison_generic` (error_eccentric ×2) | `fig_error_eccentric.jl` | `fig_error_eccentric` | fast | solid | 220 s |

The three published tables also live as a GrTools test file,
`capabilities/Perturbation/test/test_published_tables_2303_16798.jl` (commit `8bb3a57`), so the
reproduction is asserted by the package suite and not only by this tree.

`[FUTURE]` figures of this paper (documented, not reproduced): `error_spherical_prograde` /
`_retrograde` and `trajectory_difference` — both are differences against an external, un-vendored
time-domain Teukolsky / full-MPD integrator, so the plotted quantity itself is unavailable.
`CSlmnk_k` / `CSlmnk_n` and the generic `FSl{1,2,3}nk` landscapes are not yet ported (legacy
runtimes 900 s+ each).

## Against the legacy packages

Recomputing the same quantities with `emri-disk/packages/TeukolskySpinFluxes` (max relative
difference over all registered entries):

| quantity | max rel diff | verdict |
|---|---|---|
| Table III `C^+_{S,lmnk}` (24 rows) | 1.2e-12 | SAME |
| Table III `C^-_{S,lmnk}` (24 rows) | 1.7e-12 | SAME |
| Table II `F^E_{S,m}` (12 rows) | 3.5e-12 | SAME |
| Table II `F^{Jz}_{S,m}` (12 rows) | 4.8e-12 | SAME |

i.e. agreement at the Float64 re-association level — the GrTools port is the same computation, and
the residual differences against the *published* numbers below are the paper's own truncation, not a
port discrepancy.

## Measured differences vs the published numbers (findings, not tuned)

* **Table I.** Row 1 lands 0.30 % from the published −2.8259e-6 at the legacy truncation
  (l ≤ 3, |n| ≤ 2, |k| ≤ 2). Rows 2–6 are **not converged** at that truncation — measured
  2.7 % / 3.8 % / 17 % / 6.3 % / 96 % (the last at e = 0.5) — because |n| ≤ 2 misses the radial
  harmonics that carry the flux at larger eccentricity. They are therefore *not registered* as a
  reproduction here; the e-scaled `check_tableI_full.jl` covers them and is honestly marked
  `future` until it has been run.
* **Table II.** Worst |ΔF^E|/|F^E| = 2.5 % at the I = 120°, m = 3 row; worst |ΔF^Jz| = 0.72 %; the
  other rows are well under 1 %. **Paper typo carried, not silently fixed:** the I = 150°, m = 2 row
  prints F^Jz = 1.658e-5 where the value is 1.658e-4 — the golden records both the printed and the
  corrected number.
* **Table III.** Worst |ΔC^+|/|C^+| = 0.28 %, worst |ΔC^-| = 0.040 % over all 24 rows. Sign
  convention: C_GrTools = −C_paper (the flux |C|² is convention-independent).
* **`fig:fluxes_spherical`.** The published claim is reproduced mode by mode: the global peak of
  log10|F^E_{S,lm0k}| sits at (l = max(2,m), k = l − m) for m = 1, 2, 3, and the k = l − m ridge is
  the per-l argmax in 6/8, 7/8 and 8/8 of the l-columns. A VLM comparison of the regenerated
  `FSl20k.png` against the legacy render finds them indistinguishable (same `TemperatureMap`
  colorbar, same peak cell, same warm ridge, same deep-blue lobes).
* **`fig:comparison_generic` (error_eccentric).** The O(I²) message reproduces exactly: log-log
  slopes 1.9957–1.9967 (l = 2, m = 2) and 1.9615–1.9788 (l = 5, m = 4) over I = 0.5°…8°. **Known
  gap, carried unchanged from the legacy reproduction:** the published panels spread the five
  n-curves over ≈ 22×, whereas the validated flux gives 1.42× (l2m2) and 1.49× (l5m4). The inclined
  per-mode flux is oracle-matched against the BHPToolkit Mathematica source at (2,2,0,0) and
  (2,2,1,0), so the spread comes from an undocumented per-n figure-construction step, not from the
  flux; it is registered as the measured `per_n_spread_at_8deg_*` rather than reproduced.
* **`fig:comparison_spherical`.** Driver complete, **not run to completion**: the 3570-mode scan
  measured 4716 s (79 min) at −t4 with one straggler job still running when it was killed — the MST
  radial sums are slow at p = 15 with |k| ≤ 8 — against a 30-minute per-driver budget. It has no golden and `--check` exits 3 (HOLE); status `future` in the manifest.
* **`fig:fourier_coefficients`.** Peak at n = ±1 ∝ e (residual 1.5e-3) and the truncation floor
  reproduce; only the e = 0.6 floor (4.5e-12) is resolvable in Float64. The e = 0.2 / 0.4 floors sit
  at 1.7e4–3.3e6 × ε·peak — the conditioning floor of the Float64 least-squares solve — whereas the
  published figure reaches ~1e-20 (extended precision). The driver registers the measured ratios and
  makes **no** monotone-floor claim.

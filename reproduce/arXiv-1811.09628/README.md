# arXiv:1811.09628 — Miranda & Rafikov, *Multiple spiral arms in protoplanetary disks: linear theory* (MR19)

## Migration status, honestly

Before this migration NONE of the 16 figures could be regenerated inside GrTools: the model they are
built on — the Miranda–Rafikov **Newtonian, locally-cooled, planet-forced** linear disc — had no
GrTools home. `AccretionDisc.Modes` held only the RELATIVISTIC branch of the legacy `DiscModes`
package (`RelativisticDisc`, arXiv:2601.19123) and `AccretionDisc.SpiralWaves` the Page–Thorne
Lagrangian pipeline of arXiv:2608.17003 with its own Newtonian limit — a different model with a
different forcing, no cooling and a different numerical scheme.

This wave promoted the missing physics into the package, tests first
(`capabilities/AccretionDisc/test/Modes/test_newtonian_planet.jl`, gates N1–N8 against a one-shot
capture of the legacy solver):

| legacy source | GrTools home |
|---|---|
| `DiscModes/newtonian_disc.jl` (background, WKB, softened Laplace coefficients) | `AccretionDisc.Modes` — `NewtonianDisc`, `laplace_b_soft{,_deriv,_2deriv}` (unexported profiles) |
| `DiscModes/master_equation_cool_newt.jl` (3.4 k generated lines, MR20 Eqs. 9, 24–26) | `AccretionDisc.Modes.master_equation_cool` |
| `DiscModes/solver_newt.jl` (MR19 shooting method) | `AccretionDisc.Modes` — `solve_m_mode`, `solve_modes`, `vars_from_enthalpy`, `loggrid` |
| `DiscModes/integrator.jl` (adaptive DP5 with grid output) | `NumericsCore.integrate_grid` + `dp5_step_buffered!` (the 4th vendored ODE copy, folded into the single family) |
| figure-local `figures/1811.09628/fast_solver.jl` (cached-potential throughput path) | `AccretionDisc.Modes.PotentialCache` + the `solve_m_mode(d, m, x, cache)` method |
| figure-local `figures/1811.09628/homog_solver.jl` (passive-wake problem, §passive) | `AccretionDisc.Modes.passive_mode`, `passive_sigma` |

The legacy `data/` directory (797 MB, 147 files) was NEVER tracked, so nothing could be baked; every
driver here regenerates from package functions.

## Per-figure record — what regenerates today

`cost` is the measured/estimated wall time at `-t4` for the FULL figure; a forced mode solve on a
20 000-point grid costs 3–8 min per `m` (the softened-potential quadrature dominates), a passive
(homogeneous) mode a fraction of a second.

| # | published object | driver | tier | status | why |
|---|---|---|---|---|---|
| 1 | `mode_profiles` δΣ_m/Σ, m = 2,5,10 | `fig01_mode_profiles.jl` | heavy | **regenerates — SAME as legacy** | 3 forced solves in parallel, 676 s at -t4 |
| 2 | `sigma2d` 4-panel polar (headline) | — | heavy | future | needs the m = 1…120 forced mode bank + 2-D synthesis: ≈ 120 × 5 min / 4 threads ≈ 2.5 h. Command: `solve_modes(d, 2:120, loggrid(d, 0.05, 5, 6000))` then the φ-synthesis |
| 3 | `sigma2d_cartesian` | — | heavy | future | the fig-2 bank + a polar→Cartesian resample |
| 4 | `wake_profiles` 2×2 | — | heavy | future | the fig-2 bank at 4 radii. RECORDED FINDING: the legacy render used the cooling value β = 0.1 and so did NOT reproduce the published curve collapse; the adiabatic β ≳ 10³ limit is the correct regime (same fix as figs 1/2/5) |
| 5 | `wake_profiles_gamma` | — | heavy | future | 4 γ values × the fig-2 bank ≈ 10 h |
| 6 | `wake_profiles_sim` | — | — | **FUTURE (external)** | the solid curves are an external FARGO3D hydro run (N_r = 4505, N_φ = 6144); only the dashed linear-theory curves are in scope |
| 7 | `arm_profiles_h` | — | heavy | future | 3 aspect ratios × the fig-2 bank, plus an arm-extraction heuristic that the paper does not publish |
| 8 | `arm_profiles_q` | — | heavy | future | as fig 7, over the mass ratio |
| 9 | `secondary_amplitude` | — | heavy | future | parameter scan + the same unpublished extraction |
| 10 | `sigma2d_homogeneous` | — | fast | future | `passive_mode` covers the (+ik) outgoing branch of `homog_solver.jl`; the legacy fig-10 script used its own INWARD (−ik) branch, which needs its own gate before a driver can claim it |
| 11 | `wake_profiles_homogeneous` | — | fast | future | m = 2…120 passive modes (cheap) + the Gaussian outer-edge projection; the legacy script carries a third passive-solver variant, to be reconciled with `passive_mode` first |
| 12 | `secondary_amplitude_homogeneous` | — | fast | future | theory curves are cheap; the overlaid POINTS are published reference values — the authors' amplitude-ID driver is unpublished |
| 13 | `mode_phases` | — | heavy | future | m = 2…30 global modes; the legacy run capped at m = 20 as compute-prohibitive |
| 14 | `dphi_phi_theory` | — | heavy | future | the theory curves are cheap and were reproduced exactly by the legacy run; the published NUMERICAL (dotted) overlay needs the full 2-D synthesis + arm extraction |
| 15 | `phase_refine` | `fig15_phase_refine.jl` | heavy | **preliminary — does NOT match yet** | one m = 10 forced solve plus one MR App. A secant step (≈ 380 s at -t1; the caption needs no mode bank). The golden pins our own output for regression only: reference acceptance 17/21, VLM mismatch, unattributed (see "Figure 15" below) |
| 16 | `m_max` convergence | — | heavy | future | one `solve_modes(d, 1:80, x)` (≈ 1.7 h) then m-sum truncation at 10/20/40/80 |

**Summary: 2 of 16 have a driver. Fig. 1 regenerates end-to-end and matches; fig. 15 regenerates our
own computation but does not yet reproduce the published figure (preliminary). 13 are compute-bound
(the physics is now in the package, the wall time is the obstacle) and 1 (fig 6) has a genuine
external dependency.** None of the 15 in-scope figures is blocked on missing physics any more —
that is what the migration wave changed.

## Figure 1 — measured agreement with the legacy run

The driver reproduces the legacy DiscModes run on the same 20 000-point grid to every digit the
legacy record printed (`results/emri-disk/figures/1811.09628/repro_fig1.spec.md`, "Rule 3 — VLM
COMPARISON"; those values are registered in the golden as `legacy_*` beside the computed ones):

| quantity (unit (M_p/M_th)^{-1}) | m = 2 | m = 5 | m = 10 |
|---|---|---|---|
| interior envelope 0.3 < r < 0.9, computed | 0.0362796 | 0.0566626 | 0.0378439 |
| … legacy | 0.0363 | 0.0567 | 0.0378 |
| global max, computed | 0.269520 | 0.274813 | 0.181617 |
| … legacy | 0.2695 | 0.2748 | 0.1816 |
| radius of the global max | 5 | 5 | 5 |

The published amplitude ordering is recovered: m = 5 carries the largest interior wave train
(m* ≈ (2h_p)^{-1} = 5), m = 10 is the most tightly wound, and every panel grows outward to the grid
boundary — the signature the adiabatic regime restores.

**VLM comparison** of `fig01_mode_profiles.png` against the legacy `repro_fig1.png`: **MATCH** —
identical 3-panel stack, identical axes and ticks, identical curve shapes, wavelengths and envelope
magnitudes in every panel, Re black / Im red, corotation line at r = 1 with the m = 2 spike. The
published figure's own remaining difference is unchanged from the legacy assessment: it used ≈ 1e5
grid points, so its exterior red fill is denser than a 2e4-point render.

## Figure 15 — preliminary: the golden is a regression pin, not a reproduction

`fig15_phase_refine.jl` calls `AccretionDisc.Modes.mr19_phase_refinement`: one forced m = 10 solve
on the 60 000-point grid of MR19 §2 (h_p = 0.1 in the paper's adiabatic sound-speed convention,
p = q = 1, γ = 7/5, softening 0.06, β = 10⁸ as the no-cooling limit), followed by ONE secant step of
MR19 App. A on (C_in, C_out), ported from the authors' `main.py` (`phase_error_diagnostic` and the
`phgrad_steps` loop) into `Modes/PhaseGradient.jl` and exposed as `solve_m_mode(...; phgrad_steps = 1)`.
The default `phgrad_steps = 0` leaves the solver bitwise unchanged.

**The golden `expected/fig15_phase_refine.json` stores OUR computed curves and diagnostics so that a
later change to the solver is detected. It does not certify agreement with the paper.** The
preregistered comparison with the published vector figure (emri-disk
`results/emri-eccentric/paper_EMRI-ECC/evidence/mr19-fig15/`, `TOLERANCES.json` fixed before the run)
passes 17 of 21 checks. The four failures:

| check | measured | bound |
|---|---|---|
| before panel, max point-to-polyline distance (computed→published / published→computed) | 2.56 / 2.70 pt | 1.5 pt |
| after panel, max pointwise error for r ≤ 0.7 r_p and r ≥ 1.1 r_p (attained at r ≈ 1.18 r_p, on the steep flank right of r_p) | 3.5 × 10⁻³ | 5 × 10⁻⁴ |
| outer boundary X_after / X_before | 0.281 | 0.2 |

The VLM comparison is a **mismatch**: panels, axes, labels, the off-frame feature at r_p, the dip near
1.25 r_p and the suppression of the inner oscillation agree, but the inner-disc oscillation of panel
(a) is phase-shifted (its extrema sit 1–2 % inward of the published ones, giving an extra lobe beside
the trough at r ≈ 0.73). The cause is unattributed; the reviewer's divergence table between this port,
the authors' code and the paper text is in the emri-disk evidence (`review/findings.json`).

## Physics regime — the finding that governs every forced figure

MR19 (2019) is the **adiabatic, no-cooling** linear theory. The authors' shipped code
(`planet-disk-modes-cooling`) is the MR20 cooling extension, and its committed `parameters.py` sets
β = Ω t_c = 0.1, which damps the outward-propagating waves and destroys the multi-arm structure the
paper is about. The adiabatic limit is β → ∞; β = 80 already realises it (residual damping ~1e-4).
`NewtonianDisc` keeps β as a field so both regimes are reachable, and gate N6 of the package test
asserts the damping is measurable.

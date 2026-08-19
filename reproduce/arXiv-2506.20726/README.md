# arXiv:2506.20726 — spherical inspirals of a spinning secondary

Reproduced with GrTools `Spacetime` (`kerr_isso`, `SpinningSecondary`) and `Perturbation.SpinFluxes`.
Legacy sources: `emri-disk/results/gw-flux/figures/2506.20726/repro_*.jl` (read-only).

| legacy script | driver | result id | tier | status | wall (-t4) | vs legacy |
|---|---|---|---|---|---|---|
| `repro_rISSOg.jl` + `repro_drISSO.jl` | `fig_rISSO.jl` | `fig_rISSO` | fast | solid | 19 s | r_ISSO DIFFERENT by design (≤ 1.2e-5, finding 1); δr SAME (bitwise) |
| `repro_deltaJr.jl` (+ `spin_actions.jl`) | `fig_deltaJr.jl` | `fig_deltaJr` | fast | solid | 17 s | SAME (bitwise, 4 × 40 points) |
| `repro_trajectory_rz.jl` | `fig_trajectory_rz.jl` | `fig_trajectory_rz` | fast | solid | 21 s | SAME (bitwise, 400 points) |
| `packages/SpinningOrbit` testset (T1) | `check_isso_spin_shift.jl` | `check_isso_spin_shift` | fast | solid | 9 s | legacy root-finder unusable at a ≤ 1e-3 (see below) |

`[FUTURE]` figures of this paper (documented, not reproduced): `Example_non_monotonic_mode2`,
`adiabatic_evolution`, `evolution_delta_p/_delta_x`, the four `evolution_delta_Phi_*` panels,
`characteristic_strain_*` and the three `mismatches_faith_*` — all need the un-ported 0PA inspiral
driver or the waveform/LISA pipeline.

## Promotions made for these drivers (one home, no duplicates)

* `Spacetime.kerr_isso` / `Spacetime.isso_radial_curvature` — commit `0307b53`. The ISSO root-finder
  was figure-local in `repro_rISSOg.jl` / `repro_drISSO.jl`.
* `Spacetime.SpinningSecondary.spin_linear_orbital_params` (generic `pex1CTilde`, equatorial AND
  inclined) and `spin_actions_pex` — commit `b799875`. `δJ_r|_{p,e,x}` lived only in the legacy
  `SpinningOrbit` package plus the untracked figure helper `results/gw-flux/figures/spin_actions.jl`.

## Findings vs the legacy scripts

1. **ISSO root-finder (`rISSO`) — DIFFERENT, by construction.** The legacy figure scripts bisected
   on a *finite-difference* `R''(r₀)` with `h = 1e-4`; `Spacetime.kerr_isso` bisects on the
   **closed-form** `R''`, which removes both the O(h²) truncation bias and the O(ε/h²) round-off.
   Max relative difference over the published 200-point x grid, as a function of the legacy step:

   | legacy h | a = 0.3 | a = 0.6 | a = 0.9 |
   |---|---|---|---|
   | 1e-3 | 4.2e-8 | 5.1e-8 | 8.4e-8 |
   | 1e-4 (the legacy value) | 3.7e-6 | 7.7e-6 | 1.2e-5 |
   | 1e-5 | 6.3e-4 | 6.6e-4 | 9.6e-4 |

   The difference is smallest at the largest step and blows up as h shrinks — the signature of the
   FD second difference losing digits to cancellation, with the closed form sitting at the h → 0
   limit. So this is the legacy step error, not a disagreement about the ISSO.
2. **The legacy root-finder is unusable in the Schwarzschild limit** — at a ≤ 1e-3 its radius lands
   ~7e-6 *inside* the true ISSO, and the legacy `KerrGeoFrequencies` then hits `sqrt(<0)`. The
   `check_isso_spin_shift` limit (−√(8/3)) is therefore reachable with `kerr_isso` and not with the
   legacy in-script bisection.
3. **Everything else is BITWISE identical to the legacy packages** (max rel = 0.0): `δJ_r` over the
   4 × 40 published (orbit, e) points, the `(rS, zS)` trajectory correction over all 400 `w_z`
   points, and `δr_ISSO` at the sampled x = −1, ≈0, +1 for all three spins when evaluated at the
   *same* radius. The promotions are faithful ports, not re-derivations.

# arXiv:2603.13482 — a new approach to EMRIs with a spinning secondary

Reproduced with GrTools `Spacetime.SpinningSecondary` (four spin gauges) and
`Perturbation.SpinFluxes`. Legacy sources:
`emri-disk/results/gw-flux/figures/2603.13482/repro_*.jl` (read-only).

| published object | driver | result id | tier | status | wall (-t4) |
|---|---|---|---|---|---|
| `fig:delta_Omega` (4 spin gauges) | `fig_delta_Omega.jl` | `fig_delta_Omega` | fast | solid | 32 s |

Against the legacy packages: the C̃ and J gauge curves are **bitwise identical** (max rel = 0.0 over
3 × 12 sampled `p − p_s`), and the separatrix agrees to 1.6e-11 — the legacy figure script bisected
its own `r₃ → r₂` root predicate to that tolerance, `Spacetime.separatrix` is the package root.

| `fig:Clmnk_k` green "Analytical" curves | `fig_FE1_Cplus1_220k.jl` | `fig_FE1_Cplus1_220k` | fast | solid | 136 s |

`[FUTURE]` figures of this paper (documented, not reproduced): `DeltaFE1_220k` (the plotted quantity
IS the difference against the un-vendored Drummond–Hughes numeric solver), `Delta_Phi1_phi` and
`Delta_Phi1_z` (un-ported 0PA inspiral driver). `isofreq_a` / `isofreq_x` are not yet ported: the
legacy reproduction needs an adaptive BigFloat scan to ~2400 bits to reach the fold at
`p_iso − p_s ≈ 1e-21`, which is a heavy-tier driver in its own right.

## Promotions made for these drivers (one home, no duplicates)

* `Spacetime.SpinningSecondary.spin_frequencies_action_gauge` (`Ω¹|_J`) — commit `b799875`.
* `Spacetime.SpinningSecondary.spin_frequencies_constants_gauge` (`Ω¹|_C`) and
  `spin_frequencies_turning_point_gauge` (`Ω¹|_{p,e,x}`) — commit `a49ffea`. The legacy figure script
  built these two gauges from *finite-difference* Jacobians of `KerrGeoFrequencies` /
  `KerrGeoConstantsOfMotion` inline; the package versions use the analytic Jacobians (validated
  against the FD construction to rtol 1e-6, i.e. to the FD's own accuracy).

## Measured numbers (`fig:delta_Omega`, a = 0.95, e = 0.5, x = 0.5)

Separatrix `p_s = 4.0369921887351845` (`Spacetime.separatrix`; the legacy script bisected its own).
Over `p − p_s ∈ [1e-4, 10]` the maximum |Ω¹| per gauge is

| gauge | max &#124;Ω¹&#124; | near/far growth |
|---|---|---|
| C̃ (shifted constants) | 4.92e-2 | 10.3 |
| p (turning points) | 3.51e+1 | 6.90e4 |
| C (constants of motion) | 1.18e+6 | 1.58e9 |
| J (actions) | 2.01e+5 | 2.25e7 |

which is the published message quantitatively: the C̃ gauge stays **bounded** while C and J blow up
by five to six decades and the turning-point gauge grows mildly. The C̃ near/far ratio of 10.3 is the
mild bump the published panel also shows at `p − p_s ≈ 1e-4…1e-3`, not a divergence.

## `fig:Clmnk_k` — the Analytical curves

Both green curves are the central finite difference in the spin parameter (step 2e-6) of the
**pex-gauge** `TeukolskySpinModeAnalytical`; the CTilde-locked entry point is a different gauge and
lands 8–25× off. Measured at k = 0: |F^{E,1}_{220k}| = 1.806e-7 (published ≈ 1e-7) and
|C^{+,1}_{220k}| = 1.928e-5 (published ≈ 3e-5); both peak at k = 0 and fall by 21.4 and 7.3 decades
to k = +10. The published Drummond–Hughes overlays at k^s_max = 8 / 16 come from an external
un-vendored solver and are deliberately absent.

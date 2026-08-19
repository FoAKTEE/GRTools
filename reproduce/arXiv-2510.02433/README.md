# arXiv:2510.02433 — Lindblad torques on EMRIs in relativistic thin discs

Resonance locations, the Novikov–Thorne disc structure, the pressure-cutoff scalings, and the total
Lindblad torque with its power-law fits.

| legacy script (`results/emri-disk/figures/2510.02433/…`) | driver | result id | tier | status | vs legacy |
|---|---|---|---|---|---|
| `repro_LRlocation.jl` | `figLRlocation.jl` | `figLRlocation` | fast | solid | SAME (same root condition, same frequencies) |
| `repro_discProfile.jl`, `lindblad_torque_lib.jl` | `figDiscProfile.jl` | `figDiscProfile` | fast | solid | SAME to 1e-11 (characterization test in `AccretionDisc/test/Modes/test_thindisc.jl`) |
| `verify_RT_satscale.jl` | `check_satscale.jl` | `check_satscale` | fast | solid | SAME limit; small-m offsets DIFFER because the package uses the RELATIVISTIC κ (see below) |
| `repro_TorqueStrengths.jl`, `lindblad_torque_lib.jl` | `figTorqueStrengths.jl` | `figTorqueStrengths` | fast (4m43s) | preliminary | SAME to 1e-6 for the un-enhanced kernel; the published slopes are steeper (below) |
| `repro_TorqueStrengths.jl` + `relenh_grid.jl` | `figTorqueStrengths_enhanced.jl` | `figTorqueStrengths_enhanced` | heavy (~35 min) | preliminary | **DIFFERENT**: the self-consistent enhancement is < 1 and flattens the slope; the legacy grid (1.9–9.6) double-counted the geometry factor |

The figure-local library `lindblad_torque_lib.jl` (160 lines: Novikov–Thorne Σ and h, the
Shakura–Sunyaev overlay, the pressure cutoff, the signed per-mode strength, the total-torque
assembly, the power-law fit) had **no package home** and was PROMOTED into
`AccretionDisc.Modes.ThinDisc` + `AccretionDisc.Modes.Forcing` before any driver was written
(commit "promote the Novikov-Thorne thin disc and the Hirata Newtonian resonance kernel"); every
number below is a package call.

## Reproduced

* **fig:discProfile** — Σ(r) and h(r) for a/M = 0 and 0.9 with the SS overlay. Landmarks:
  Σ_min = 1.212e5 g/cm² at r = 13.41M (a=0) and 2.823e4 at 4.77M (a=0.9); h_max = 0.0563 at 18.72M
  (a=0) and 0.1417 at 6.86M (a=0.9); the published normalization anchor h(a=0, 1000M) = 0.0033 is
  exact by construction and the large-r log-slopes are 3/2 and −1 (the SS scalings).
* **fig:LRlocations** — |r_LR − r_s| for m = 2, 3, 4, 8 at a/M = 0 and 0.9.
* **saturation scalings** — n_r = 3/2 − 2Σ_p = 9/2 for the radiation-pressure inner disc;
  |r_LR − r_s| → 2H/3 (rel 1.5e-3 at m = 5000); the cutoff falls below ½ at ξ = mh = 0.392.
* **fig:TorqueStrengths / tab:FitParameters ("All" rows)** — see the finding below.

## Findings (not tuned away)

1. **The un-enhanced kernel under-predicts the NT slopes by ~0.2–0.4, and the SS row is exact — and
   the relativistic amplitude correction does NOT close the gap.**
   Fitted n_r (this driver) vs published: a=0 **3.19** vs 3.57; a=0.3 **3.31** vs 3.66;
   a=0.6 **3.49** vs 3.77; a=0.9 **3.70** vs 3.92; **SS 4.387 vs 4.39 (rel 6e-4)**.
   The Shakura–Sunyaev disc is the case where Hirata's NEWTONIAN overlap amplitude is the correct
   one, and it lands on the published slope; every relativistic-disc row is short.
   The natural candidate was the relativistic AMPLITUDE enhancement E(m,r) = |N_Teukolsky /
   N_Newtonian| that `lindblad_strength` does not carry, so `figTorqueStrengths_enhanced.jl`
   (heavy tier) measures it directly from `lindblad_torque` / `lindblad_strength` on a
   9 × 6 × 2 grid at a = 0. **It measures E < 1 everywhere, falling with m** (0.86 / 0.67 / 0.48 /
   0.32 at m = 2 / 5 / 8 / 10, r₀ = 30M, ℓ-converged), so applying it *flattens* the slope instead
   of steepening it. The legacy "enhancement" of 1.9–9.6 was computed against the PURE-Newtonian
   N (fixed ς_m^∓ = ((m∓1)/m)^{2/3}), which double-counts the relativistic resonance-geometry
   shift that the hybrid kernel already carries. **The published slope steepening is therefore not
   explained by the Teukolsky amplitude — it remains an open obligation, not a tuning knob.**
   Measured enhanced fits (converged grid, `lmax = max(16, m+8)`): n_r = **2.54 / 2.70 / 2.88 /
   3.13** for a = 0 / 0.3 / 0.6 / 0.9 — i.e. the gap to the published 3.57 / 3.66 / 3.77 / 3.92
   *widens* from −0.38 to −1.03 at a = 0. (ℓ-truncation matters and was measured independently: at
   lmax = 12 the ratio is 52 % low at m = 10 and 32 % low at m = 8; E is converged by lmax = 16 for
   m ≤ 10, and the converged grid reproduces the probe, E(m=10, r₀=30M) = 0.318 vs 0.3175.)
2. **The torque amplitude A is not predicted in absolute units.** The paper's 10⁻⁹/q² axis is set by
   one global unit factor; here κ_U = 2.281e-18 is anchored to the published a=0 amplitude
   A = 2.33e-9, exactly as the legacy analysis did. The SPIN ORDERING is a prediction:
   A(0.3)/A(0) = 0.696 (published 0.725), A(0.6)/A(0) = 0.426 (0.498), A(0.9)/A(0) = 0.230 (0.309).
3. **Pressure-modified resonance at small m.** `check_satscale` uses the package
   `lindblad_locations_pressure`, which keeps the RELATIVISTIC κ_r; the legacy check used the
   Newtonian κ = Ω. At m = 5 and r_s = 50M this shifts the offset from 6.47M (legacy) to 6.13M.
   The saturated limit — the statement being verified — agrees to 5e-4 at m = 5000.
4. **SS reversal radius is NaN** on the 50-point grid registered here (no sign change between 6.5M
   and 100M at this resolution); the NT reversal radii are 11.82 / 10.36 / 8.11 / 5.02 M for
   a = 0 / 0.3 / 0.6 / 0.9.

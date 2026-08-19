# arXiv:1010.0759 — Hirata, *Resonant recoil in extreme mass ratio binary black hole mergers II*

Relativistic Lindblad resonances of a circular perturber in a Schwarzschild/Kerr disc: resonance
locations, the resonant torque strength `N`, and the amplitude-convention reconciliation that ties
the paper's `c₁₃`, `α`, `Z^out/Z^down` to the GrTools Teukolsky engine.

| legacy script (emri-disk `results/emri-disk/figures/…`) | driver | result id | tier | status | vs legacy |
|---|---|---|---|---|---|
| `rt_n1_torque_strength.jl`, `1010.0759/verify_RT_A2_Sm.jl` | `check_RT_N1.jl` | `check_RT_N1` | fast | solid | SAME (closed form; ≤1.2 % of the published −2.36/−7.50/0.19/8.62) |
| `c13_alpha_reconcile.jl` | `check_c13_alpha.jl` | `check_c13_alpha` | fast | solid | **DIFFERENT (corrected)**: the package α now uses Γ = 2√(1−a²)/(r₊²+a²); see below |
| `1010.0759/repro_tab_ilrS.jl`, `verify_tab_ilrS.jl` | `tabIlrS.jl` | `tabIlrS` | fast (3m18s) | solid | SAME to 1e-4 on locations and ILR strengths |
| `1010.0759/repro_ilrplot.jl` | `figIlrplot.jl` | `figIlrplot` | fast (3m11s) | solid | SAME (same kernel; 12-point grid instead of 24) |
| `1010.0759/repro_kerrm1_locations.jl`, `verify_kerr_m1ILR_loc.jl` | `figKerrm1.jl` | `figKerrm1` | fast | solid | SAME (anchors 39.828 / 46.410 / 43.249) |
| `1010.0759/repro_lkerr_locations.jl` | `figLkerr.jl` | `figLkerr` | fast | solid | SAME |

## What is reproduced, and how well

* **Locations** (tab:ilrS, fig:kerrm1, fig:lkerr, fig:ilrplot left): the roots of
  `m(Ω_φ − Ω_φ(r₀)) = ±κ` on the exact Kerr equatorial geodesic frequencies. Max relative deviation
  from the published 4-significant-figure table: **5.6e-4** — i.e. the table's own rounding.
  Published Kerr anchors r₀ = 250: 39.83 (a=+0.9, paper 39.8), 46.41 (a=−0.9, paper 46.4).
* **ILR strengths** (tab:ilrS, fig:ilrplot right): the full Teukolsky overlap
  `P^{(m)} = −ω⁻² Σ_ℓ Z^down_ℓm (c₁₃ Z^out*_{1,ℓm} + α Z^down*_{1,ℓm})` → `S^{(m)}` → `N`.
  Ratio to the published column: **m=1 ILR 0.99–1.00 for r₀ ≥ 12M**, **m=2 ILR 0.93–0.97**.
  The two outliers (r₀ = 8, 9 for m=1) sit where the ILR is within 1M of the ISCO.
* **Newtonian limit** (RT-N1): −2.356 / −7.504 / +0.1922 / +8.615 against the published
  −2.36 / −7.50 / 0.19 / 8.62, from the closed form `lindblad_strength_newtonian`.

## Findings (not tuned away)

1. **Horizon factor α corrected.** The legacy scripts and the legacy `DiscModes` used
   Γ = 4√(1−a²)/(r₊²+a²), which is Hirata's AS-PRINTED bracket and misses the eq:cplex1 closure
   `|c₁₃|² + α = (2ω)⁶|ℵ|²/|C|²` by up to **9.1e-5**. The Teukolsky–Press Γ = 2√(1−a²)/(r₊²+a²)
   closes it to **≤6e-13** (this is the gate in `check_c13_alpha.jl`). `AccretionDisc.c13_alpha`
   was corrected; `lindblad_torque` moves by only 2.0e-8 because α ≪ |c₁₃|².
2. **OLR strengths are NOT registered.** eq:P-OLR is a near-cancellation of two large terms; in
   Float64 the legacy run reproduces the published OLR column only to factors 1.5–260
   (`results/emri-disk/figures/1010.0759/repro_tab_ilrS.out`). Registering those numbers would
   register noise. Reproducing them needs the BigFloat overlap — an open obligation.
3. **fig:ilrplot m=3 asymptote.** `N(m=3, r₀=250) = −6.41` vs the published asymptotic −7.50 (14 %).
   The ℓ-sum is truncated at ℓ ≤ 24; the m=3 tail converges more slowly than m=1,2.

## Regenerating

```bash
julia --project=reproduce reproduce/run_all.jl --check --paper 1010.0759
julia --project=reproduce -t4 reproduce/arXiv-1010.0759/tabIlrS.jl --regen   # rewrite one golden
```

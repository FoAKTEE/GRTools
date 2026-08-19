# arXiv:1010.0758 — Hirata, *Resonant recoil in extreme mass ratio binary black hole mergers I*

Paper I is a DERIVATION paper: nine of its ten reproduced results are symbolic certificates, kept
here as LaTeX under `derivations/` (copied verbatim from the mission repository
`results/emri-disk/paper/1010.0758/`; the canonical, evolving copies live there). The single
quantitative external validation — the Keplerian limit reproducing Goldreich & Tremaine (1980)
Eq. (18) — is the numeric driver.

| derivation certificate | content | numeric driver |
|---|---|---|
| `RT1_C_coefficients.tex` | the C-coefficients, κ and Z of the resonance expansion | — (checked in `arXiv-1010.0759/check_c13_alpha.jl` conventions) |
| `RT2_perturbation_hamiltonian.tex` | the perturbation Hamiltonian H₁ = h^{ab}p_a p_b/(2p^t) | — |
| `RT3_resonant_integrals.tex` | resonant integrals / driven-SHO reduction | — |
| `RT4_master_torque_collapse.tex`, `RT4_torque_density_chain.tex` | the master torque T ∝ \|S^{(m)}\|² and the torque-density chain | `check_RT8_GT80.jl` (Keplerian limit) |
| `RT5_explicit_Sm.tex` | the explicit S^{(m)} (a factor-2 typo located) | `arXiv-1010.0759/check_RT_N1.jl` (closed form) |
| `RT6_gauge_invariance.tex` | gauge invariance of S^{(m)} (16 cases ≤ 6e-16) | — |
| `RT7_epicyclic_waveform.tex` | epicyclic waveform forms | — |
| `RT8_GT80_newtonian_check.tex` | GT80 Eq. (18) transcription + the 32/81 reduction | **`check_RT8_GT80.jl`** |
| `RT9_disc_heating_flux.tex` | disc-heating flux F and u^r | — |

| driver | result id | tier | status | vs legacy |
|---|---|---|---|---|
| `check_RT8_GT80.jl` | `check_RT8_GT80` | fast | solid | SAME (the 32/81 coefficient and the ς_m^∓ limits are re-derived from the package frequencies) |

## What `check_RT8_GT80.jl` checks

* the GT80 prefactor `κ²Ω²/(2³A⁴)` computed from the package `disc_OmegaK`/`kappar2` and the Oort
  constant `A = (r/2)dΩ/dr` → **0.3950594** at r = 10⁶M vs the exact **32/81 = 0.3950617**
  (rel 6.0e-6 — the residual is the O(M/r) relativistic correction, which is the point);
* κ/Ω → 1 and A/Ω → −3/4 in the same limit (the Keplerian Oort values);
* the Lindblad ratios ς_m^∓ = r₁/r₀ → ((m∓1)/m)^{2/3} to **3.2e-6** (ILR) and **5.1e-7** (OLR),
  from the same relativistic root solve — this is the RT4/RT5 tie to the Newtonian closed form;
* the modified-Bessel moment `[2K₀(2/3)+K₁(2/3)]² = 6.048740` is carried as a LITERATURE-GROUNDED
  constant (certificate `derivations/RT8_GT80_newtonian_check.tex`, mpmath 25 dps): no GrTools
  package owns modified Bessel functions, so it is not recomputed here. That is the one
  unchecked external step of this driver.

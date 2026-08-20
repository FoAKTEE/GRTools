"""
    TimeDomain

Multi-domain nodal discontinuous-Galerkin **time-domain** solver for the distributionally
forced s = 0 (scalar) Teukolsky equation on a Kerr background — the physics stack of
DG-Teulolsky (arXiv:2307.01349 "Paper 1" method; arXiv:2503.11523 "Paper 2"
superconvergence/self-force), ported as a Perturbation submodule (DG-P3).

Contents (mirrors the reference repo's physics split):
  - `hyperboloidal.jl` — hyperboloidal-layer compactification `Layer`, `Ωfun`, `Hfun`, …
    (Paper 1 §II.D, Paper 2 §II.A).
  - `background.jl` — scalar-Teukolsky mode potential `potential_V`, layer-coupling
    `coupling_f`, imaginary friction `friction_mu` (Paper 1 §II). The tortoise coordinate,
    its inverse, horizons and Δ DELEGATE to `Spacetime` (M = 1 canon, exact M-scaling shims).
  - `mode1d.jl` — single-mode Schwarzschild evolution `Mode1DProblem`/`solve_mode1d` +
    Price-tail fitter `fit_tail` (Paper 1 §IV.B).
  - `flux1d.jl` — Schwarzschild circular-orbit flux/self-force `Flux1DProblem`/`solve_flux`/
    `selfforce_ft` with the characteristic-projected delta-flux injection (Paper 1 §III.C,
    Paper 2 §IV.B); `orbital_params` and `Ylm_equator` DELEGATE to `Spacetime`/`Harmonics`.
  - `kerrmodes.jl` — Kerr m = 0 coupled-mode tails `KerrTailProblem`/`solve_kerr_tail`
    (Paper 1 §II.E–F, F1.4).
  - `kerrflux.jl` — Kerr coupled m ≠ 0 circular-orbit flux `KerrFluxProblem`/
    `solve_kerr_flux` with eigen-decomposed source injection (Paper 1 §IV.C Kerr, T1.3/F1.7).
  - `spinweight.jl` / `rtchart.jl` / `spinmodes.jl` — DG-P4 generic-spin (s = ±1, ±2) VACUUM
    m-mode ℓ-coupled sector in the Harms–Bernuzzi–Brügmann RT horizon-penetrating
    hyperboloidal chart (arXiv:1301.1591; proven-regular coefficients at scri for s ≠ 0 —
    the repo layer chart diverges there, dossier DG_SPIN_FORMULATION.md §3.4). Spin-weighted
    coupling coefficients per Barack–Giudice (arXiv:1702.04204 Eqs. (32)–(35));
    `SpinModesProblem`/`solve_spin_modes`. The s = 0 layer-chart path above stays bit-untouched.
  - `spinsource.jl` — DG-P5 SOURCED circular-equatorial point-particle sector at a = 0
    (s ∈ {0, −2, +2}): closed-form δ/δ′/δ″ jump conditions (Barack–Giudice arXiv:1702.04204
    App. A + the GHP source of the in-package FD oracle) enforced as jump-conditioned
    numerical fluxes at the particle interface (arXiv:2211.14390 pattern);
    `SpinSourcedProblem`/`solve_spin_sourced`/`spin_sourced_flux`, gated against the
    frequency-domain oracle `TeukolskyPointParticleModeCircular`.
  - `spinsource_kerr.jl` — DG-S2 SOURCED circular-equatorial sector at KERR a ≠ 0
    (s ∈ {0, −2, +2}): closes the dossier §8.3 spheroidal→spherical re-projection [FUTURE]
    of DG-P5 via Harmonics' SphericalExpansion b-coefficients (per-spheroidal jumps b-summed
    onto the spherical band; Kerr RT chart push-forward DERIVED-CK) LAYERED on the a = 0
    machinery (additive dispatch only); `SpinSourcedKerrProblem`/`spin_sourced_flux_kerr`,
    gated against the FD oracle at a ∈ {0.5, ±0.9} + the SKH published tables
    (arXiv:gr-qc/0703028). Eccentric/inclined Kerr sourced TD remain out of scope
    (monochromatic jump ansatz); |s| = 1 sourced remains [FUTURE].
  - `hphc.jl` — DG-D3 HPHC (Horizon-Penetrating Hyperboloidal-Compactified) **spectral-collocation**
    VACUUM Teukolsky evolution backend (spin-`s` field about Kerr, `a ≠ 0`): port of
    TeukEvolution.jl's `linear_field` path (Ripley, JLRipley314/TeukEvolution.jl). The
    2nd-order-in-time master equation is reduced to `(f,p) = (ψ, ∂_Tψ)` and RK4-evolved on the
    compactified radial grid `R = cl²/r` × the Gauss-Legendre `y = −cosϑ` spin-weighted collocation
    (same RLGP-N chart as Move B). Coefficients `A_pp,A_pq,B_pp,B_pq,B_pf,pre` are the verbatim
    Evolution.jl closed forms (`pre` = the PAPER-VERIFIED inverse `∂_T²ψ` coefficient, RLGP P1);
    FD/KO from `NumericsCore.fd_uniform` (D1), angular collocation from `Harmonics.SphereCollocation`
    (D2), integration from `NumericsCore.rk4!`. `HPHCTeukolskyProblem`/`solve_hphc`/`solve_hphc_mmodes`
    (embarrassingly `m`-parallel). Characterized end-to-end vs the ORIGINAL TeukEvolution run
    (test/fixtures/td_hphc_upstream.jl, a=0.7 s=−2 m=2). This is the distinct HPHC **evolution**
    engine — complementary to the DG **layer-chart** s=0 stack above; no numerics duplicated.
  - `hphc_reconstruction.jl` — DG-D4 (Move D) HPHC VACUUM metric-**reconstruction** hierarchy: port of
    TeukEvolution.jl's `runtype="reconstruction"` path (Ripley/Loutrel et al. arXiv:2010.00162 §II.D +
    arXiv:2008.11770 §V), DEFECT-FIXED and LAYERED on the D3 ψ₄ backend above. The ordered vacuum
    radiation-gauge transport chain ψ₄ → {ψ₃,λ} → {ψ₂,h_m̄m̄} → π → h_lm̄ → μh_ll is evolved as one
    18-field ±m-coupled state by ONE `NumericsCore.rk4!`. Grid GHP applicators
    (`recon_{edth,edth_prime,thorn,thorn_prime}!`) are thin loops SOURCED from `Spacetime.NPGHP`
    coefficient fields (no closed form re-hardcoded). Defect fixes vs upstream: P1 (the ψ₄ momentum
    `pre` normalization), P2 (one standard RK4 — no malformed stage-4 weights / `muhll` index typo),
    P3 (the angular killtop applied in reconstruction mode), quirk-1 (þ′ implemented fresh, not the
    dead `thorn_prime_arr` body). `HPHCReconProblem`/`solve_hphc_recon`/`hphc_recon_residual3`;
    the ψ₄ sector matches the D3 linear_field backend by construction (test/fixtures/td_hphc_recon_upstream.jl).

The problem-agnostic DG core (LGL `RefElement`, `Mesh`, `build_mesh`, `rk4!`) is REUSED from
`NumericsCore.dg` (already the characterized port of the repo's DGCore/TimeStepping); no DG
numerics are duplicated here. Characterized end-to-end against baked one-shot outputs of the
ORIGINAL repo (test/fixtures/td_*.jl): field samples, tail exponents, fluxes and self-force
at near-machine agreement.
"""
module TimeDomain

using LinearAlgebra: mul!, dot, eigen, inv, norm, I as Id
# DG core + RK4 (NumericsCore.dg = the characterized DGCore/TimeStepping port; commit faa4916)
using NumericsCore: RefElement, Mesh, build_mesh, element_of, lagrange_interp, node_coords,
                    rk4!, cfl_dt, bary_weights, bary_row
# Uniform-grid FD stencils + KO dissipation (NumericsCore.fd_uniform = the TeukEvolution Radial port,
# D1 commit 6aaf6d82) — the HPHC spectral-collocation evolution backend (hphc.jl).
using NumericsCore: set_d1!, set_d2!, ko_filter!
# Canonical Kerr background + circular-equatorial geodesics (M = 1 convention; commit 8be5f32).
# (tortoise_kerr_inverse is NOT imported: the deep-horizon floor convention differs — see
# TimeDomain/background.jl header for the measured evidence.)
using Spacetime: tortoise_kerr, geodesic_energy, geodesic_angular_momentum, kerr_frequencies
# NP/GHP coefficient kernel (Spacetime.NPGHP, Move B): the six HPHC NP_0 scalars, the eleven GHP
# pre_* operator coefficient fields, and the {p,q}=ghp_weights(s,b) integer weights — the coefficient
# surface the D4 reconstruction backend (hphc_reconstruction.jl) layers its grid GHP applicators on.
using Spacetime: NPGHP
# Rácz–Tóth PURE chart maps re-homed to Spacetime.RTChart (F5 Part B, MINIMAL move). RE-IMPORT +
# RE-EXPORT the SAME objects (Move-A Lorenz shim): every `rt_r_of_R`/`rt_R_of_r`/`rt_Rplus`/
# `rt_F_of_R` call site in rtchart.jl/spinsource.jl/spinmodes.jl/spinsource_kerr.jl resolves to the
# Spacetime binding (`===`), so the locked time-domain physics is bit-untouched. The RT
# master-equation COEFFICIENT fields (rt_C*, rt_A*) stay defined below (application-side geometry).
using Spacetime.RTChart: rt_r_of_R, rt_R_of_r, rt_Rplus, rt_F_of_R
# Scalar spherical harmonic (equatorial value path)
using Harmonics: Ylm_real
# Gauss-Legendre spin-weighted angular collocation (Harmonics.SphereCollocation = the TeukEvolution
# Sphere port, D2 commit f8bda83c): grids, swal values, Laplacian / killtop matrices + angular apply
# — the HPHC evolution backend's angular sector (hphc.jl).
using Harmonics: swal, Y_vals, cos_vals, sin_vals, swal_laplacian_matrix, swal_killtop_matrix,
                 swal_raising_matrix, swal_lowering_matrix, angular_matrix_mult!, inner_product
# Spin-weighted harmonic values at aω = 0 (≡ exact ₛYₗₘ; DG-P5 sourced sector) and the
# SphericalExpansion b-coefficients ₛS_{ℓ̂mω} = Σ_ℓ b_ℓ ₛYₗₘ (DG-S2 Kerr sourced sector —
# the spheroidal→spherical source re-projection, dossier §8.3 closure)
using Harmonics: SpinWeightedSpheroidalHarmonicS, swsh_derivative, expansion_coefficients
# Amplitude → flux layer (BHPToolkit normalization; AmplitudesFluxes/fluxes.jl) for the
# DG-P5 scri/horizon flux extraction — defined in the parent module before this include.
import ..energy_flux

include("TimeDomain/hyperboloidal.jl")
include("TimeDomain/background.jl")
include("TimeDomain/common.jl")
include("TimeDomain/mode1d.jl")
include("TimeDomain/flux1d.jl")
include("TimeDomain/kerrmodes.jl")
include("TimeDomain/kerrflux.jl")
include("TimeDomain/spinweight.jl")
include("TimeDomain/rtchart.jl")
include("TimeDomain/spinmodes.jl")
include("TimeDomain/tailfits.jl")
include("TimeDomain/spinsource.jl")
include("TimeDomain/spinsource_kerr.jl")
include("TimeDomain/hphc.jl")
include("TimeDomain/hphc_reconstruction.jl")
include("TimeDomain/hphc_qnm_id.jl")
include("TimeDomain/hphc_scd_source.jl")
include("TimeDomain/hphc_scd_solve.jl")

export Layer, Ωfun, Ωprime, denomΩ, Hfun, omega_factor, rstar_of_ρ
export potential_V, coupling_f, friction_mu
export Mode1DProblem, solve_mode1d, ρ_of_rstar, fit_tail
export Flux1DProblem, solve_flux, selfforce_ft, Ylm_equator, orbital_params
export cminus, mode_coupling, mode_list, KerrTailProblem, solve_kerr_tail, fit_tail_window
export KerrFluxProblem, solve_kerr_flux
export spin_cminus, spin_czero, spin_cplus, spin_cos_coupling, spin_mode_coupling,
       spin_mode_list
export rt_Rplus, rt_r_of_R, SpinModesProblem, solve_spin_modes, spin_cfl_dt
export fit_tail_envelope, fit_tail_adaptive, fit_tail_sliding, fit_convergence_order,
       spin_excise_policy, TAIL_SLIDING_STARTS
export rt_R_of_r, spin_source_jumps, bg_jumps_s2, SpinSourcedProblem, solve_spin_sourced,
       spin_sourced_flux
export spin_source_jumps_kerr, rt_H_of_R_kerr, SpinSourcedKerrProblem, spin_sourced_flux_kerr
export HPHCTeukolskyProblem, solve_hphc, solve_hphc_mmodes
export hphc_A_pp, hphc_A_pq, hphc_B_pp, hphc_B_pq, hphc_B_pf, hphc_pre
export HPHCReconProblem, HPHCReconGHP, solve_hphc_recon, hphc_recon_residual3
export hphc_qnm_seed, hphc_qnm_amplitude, hphc_qnm_tower_mode, hphc_qnm_initialdata,
       solve_hphc_recon_qnm
export SCDSourceOps, SCDBundle, scd_bundle, scd_zero_bundle, scd_sd!, scd_st!,
       scd_source!, scd_source_sector
export scd_source_weight, scd_ramp, solve_scd_driven, scd_project_lm, scd_fit_qqnm, scd_qqnm_ratio

end # module TimeDomain

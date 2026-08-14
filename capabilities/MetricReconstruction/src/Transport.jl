"""
    MetricReconstruction.Transport

Li–Yunes traceful radiation-gauge (traceful ORG) metric reconstruction for GENERIC (non-vacuum)
type-D sources — the third gauge-disjoint submodule of `MetricReconstruction` (peer of `Radiation`
and `Lorenz`), porting arXiv:2605.11080v1 (Li & Yunes). Physics: reconstruct the linear metric
perturbation `h_ab^{(1)}` from a single radiative Weyl scalar `Ψ4^{(1)}` and the stress-energy
tensor `T_μν^{(1)}`, via the trace-transport head (the innovation — `h_mm̄^{(1)}` direct from
`T_nn^{(1)}`, no vacuum/RLGP analogue) followed by the seven generalized NP transport equations
(Fig. 1 hierarchy).

Layout (dossier §4):
  - `radial_profile.jl`  — the HOLE-2 distributional radial-profile algebra (`Laurent`,
                           `ShellProfile`): the closed-form value types the static-shell chain uses.
  - `operators.jl`       — the weighted NP operators Δ_{[a,b,c,d]}/D/δ/δ̄ (LY eq:derivative-notation)
                           as an application layer over `Spacetime.NPGHP` coefficients.
  - `sources.jl`         — T_μν^{(1)} → (Φ_ij, Λ) + the Bianchi sources R_g, R_h (HOLE-6).
  - `trace_branch.jl`    — A1 mu_transport_final + A2 trace_transport (the spin-0 head).
  - `transport_chain.jl` — B1…B7 (the seven generalized transport eqs) + the S_μ auxiliary + the
                           Fig.-1 DAG driver `transport_solve`.
  - `completion.jl`      — the traceful-ORG assembly h_ab→h_μν, the exact-solution target, ξ_t gauge.

DECISIVE gate: the paper's own Schwarzschild static thin-shell worked example (Supplemental
Material) — the reconstructed `h_μν^{(1)}` equals the linearized exact Frauendiener metric after the
explicit ξ_t gauge transform (test/transport_shell_example.jl). Every equation cites its LY label;
the whole shell chain is closed-form and re-runs at BigFloat-256. Vacuum cross-gate: the transport
operators' weight quartets vs the D4 assembler algebra (test/transport_vacuum_crossgate.jl).

Background NP/GHP quantities are SOURCED from `Spacetime.NPGHP` (never re-hardcoded; the a=0 Laurent
table is gated field-by-field against `np_spincoeff(0,r,θ)`). {T}-generic (BigFloat 256…2048-bit
inherited).

**[AXIOM-E] inputs** (physical choices not derivable from the transport ODEs — documented defaults):
outside-Schwarzschild branch amplitudes (`c0⁺=c1⁺=0`, `Cp=0`); the mass-shift BCs (`Ψ2,reg⁺=−M/r³`,
`[Ψ2,reg]=−M/R³`); continuity of `h_mm̄`, `h_ll`; the free residual-gauge constants `C_t^±`.
"""
module Transport

using Spacetime

include("Transport/radial_profile.jl")
include("Transport/operators.jl")
include("Transport/sources.jl")
include("Transport/trace_branch.jl")
include("Transport/transport_chain.jl")
include("Transport/completion.jl")
include("Transport/smooth_source.jl")   # F4: generic SMOOTH static-source path via NumericsCore DP5
include("Transport/dipole_completion.jl")  # F6-B: static ℓ=1 axial angular-momentum (Lense–Thirring)

# distributional radial-profile algebra (HOLE-2)
export Laurent, laurent, lzero, lscale, lderiv, lintegrate, leval, lapprox
export ShellProfile, shellprofile, pdelta, pscale, psmooth, psmooth_vd, pderiv, pjump, pcont, peval
# operators
export WeightedOp, Δw, Dw, δw, δbarw, weight_combo, schw_f, schw_np_laurent, delta_static
# sources
export shell_eos, shell_T_np, np_ricci_from_T, shell_ricci_profiles, Rg_shell, Rh_shell
# trace branch
export solve_hmm_shell, mu_from_trace, gamma_from_trace
# transport chain
export CHAIN_DELTA_OPS, vanishing_shell, solve_psi2_shell, smu_cancel_shell
export HllField, hll_eval, solve_hll_shell, ShellSource, SmoothSource, transport_solve
# completion / assembly / gauge
export recon_htt, recon_hrr, recon_htr, recon_hφφ_over_sin2
export exact_htt, exact_hrr, exact_hφφ_over_sin2, gauge_ode_rhs, gauge_xi_t
# F4 smooth-source DP5 path
export SmoothSolution, smooth_norm, smooth_sigma, smooth_sigma_p, smooth_pt, smooth_ricci
export smooth_rhs, smooth_bc, smooth_solve_targets, conservation_residual
export smooth_mu, smooth_hmm, smooth_psi2, smooth_hll, smooth_gamma
export smooth_htt, smooth_hrr, smooth_hφφ_over_sin2
# F6-B dipole (Lense–Thirring) completion
export DipoleCompletion, dipole_psi3_radial, dipole_psi3, dipole_pi_radial, dipole_pi
export dipole_hlmbar_radial, dipole_hlmbar, dipole_psi2_radial, dipole_psi2
export dipole_lambda, dipole_hmbarmbar, dipole_hmmbar, dipole_hll
export dipole_assemble_org, dipole_gauge_c, dipole_gauge_hrphi, dipole_h

end # module Transport

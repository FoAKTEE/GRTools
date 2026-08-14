"""
    SelfForce.LorenzField

Lorenz-gauge BLS first-order gravitational field solver (legacy `H1Lorenz`). CONSUMES
`MetricReconstruction.Lorenz.LorenzRadial`'s radial RHS (`even_ode!`/`odd_ode!`/`monopole_ode!`/
`odd_dipole_ode!`) — it does not own them — assembles the R1..R10 sources and the in/out
boundary-condition series, integrates with the NumericsCore vendored DP5+GBS, builds the gauge
field, and WRITES the h1 corpus via `Data.write_mode_h5`. Type-generic `{T}`.

The Schwarzschild geometry (lapse/tortoise) and circular-orbit constants of motion are CONSUMED from
`Spacetime`; the normalized associated Legendre factors of the BLS sources from `NumericsCore`.
"""
module LorenzField

import Spacetime
using NumericsCore: assoc_legendre_P, legendre_P_thetaderiv

# ── Circular-orbit parameters + frequency-domain BLS source coefficients R1..R10 ─────────────────
include("LorenzField/Orbit.jl")

# ── Analytic static first-order fields: l=m=0 closed-form static monopole ─────────────────────────
include("LorenzField/Static.jl")

# ── Homogeneous in/out boundary-condition series per parity sector (seeds the radial integrator) ──
include("LorenzField/Boundary.jl")

# ── Coupled-set descriptors + parity selection + γ↔(l,m) enumeration (leaf solve-plan layer) ──────
include("LorenzField/CoupledSets.jl")

# ── VoP/two-point-shoot integrator → gauge assembly → NModeData → h1-corpus writer (CLOSES it) ────
include("LorenzField/Solver.jl")

export OrbitalParams, orbit_energy, orbit_ang_mom, orbit_phi_freq
export CoupledSet, CSETS, select_coupled_set, gamma_to_lm, lm_to_gamma, mode_list
export integrate_field_equations, calculate_scaling_coefficients
export NModeData, compute_mode, compute_modes, compute_one_mode, write_h1_mode, h1_schema
export R1_source, R3_source, R5_source, R6_source, R7_source, R9_source, R10_source
export rs_of_r, r_of_rs
export construct_static_monopole
export even_out_bcs, even_in_bcs, even_dipole_out_bcs, even_dipole_in_bcs
export monopole_out_bcs, monopole_in_bcs, odd_out_bcs, odd_in_bcs
export odd_dipole_out_bcs, odd_dipole_in_bcs, even_static_out_bcs, even_static_in_bcs

end # module LorenzField

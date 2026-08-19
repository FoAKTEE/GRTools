"""
    Spacetime

Kerr background-geometry + orbital-mechanics FOUNDATION for GrTools — the Tier-2 spine
root whose dependency closure is EXACTLY `{LinearAlgebra, NumericsCore}` (no Roots, no
HDF5/Symbolics/Wigner/ODE). Every capability above (Perturbation, MetricReconstruction,
SelfForce, EMRIWaveforms, ...) reaches Kerr geometry/geodesics through here, so this
package carries no heavy dependency by construction.

Owns the physics that is *geometry of the Kerr background and motion within it* (per design
§4.1). Type-generic `{T}` throughout so BigFloat precision is inherited, never hard-Float64'd.

Submodules (built incrementally, tests-first):
  - `Geometry`     — BL/In/Out metric + inverse, Kinnersley tetrad, Christoffels + r/θ
                     derivs (incl. the a=0/M-only Schwarzschild method from Punctures),
                     bg-scalar bundle Δ/Σ/ζ/K/Q, tortoise + inverse, horizons, lapse/
                     potential, circular-orbit E/Lz/Ωφ. Coords parametrized as METHODS
                     (BL/In/Out + `:eq` restriction), NOT flattened forks.
  - `Geodesics`    — E/Lz/Q constants of motion, radial/polar roots, Fujita-Hikida Mino +
                     BL frequencies Ω_r/Ω_θ/Ω_φ, closed-form trajectories, proper time.
  - `Separatrix`   — ONE canonical `separatrix(a,e,x)` + `kerr_isco(a)` (the most-copied
                     primitives) over NumericsCore's vendored A42/brent bracketing solver.
                     The FEW `get_separatrix` / KerrBackwards `separatrix` become thin aliases.
  - `SpinningSecondary` — linear-in-spin orbital actions/frequencies + spin-curvature force
                     terms (from SpinningOrbit).
  - `PNKinematics` — numeric Kepler Ω / r0_of_Ω / pn_x (the NUMERIC kinematic helpers and the
                     collapsed single `kerr_isco`/PN-ISCO fork — NOT the symbolic PN series).

Legacy packages absorbed (kept GREEN as the characterization oracle): KerrGeodesics,
SpinningOrbit, plus the geometry forks from MetricReconstruction / DiscModes / KerrLorenzCirc
/ Punctures and the separatrix/ISCO forks from KerrBackwards / FastEMRIWaveforms / WaSABI /
DiscEMRIForce.
"""
module Spacetime

using LinearAlgebra
using NumericsCore

# Submodule files are added as each primitive is TDD-ported and grounded against the legacy
# oracle. Each include + its exports are uncommented as the module lands.

include("Geodesics/constants.jl")
export constants_of_motion, geodesic_energy, geodesic_angular_momentum, carter_constant

include("Geometry/geometry.jl")
export kerr_metric, kerr_metric_inv, kinnersley_tetrad, bg_scalars, horizons, omega_plus
export delta_kerr, sigma_kerr, zeta_kerr, kfun, qfun
export tortoise_kerr, tortoise_kerr_inverse   # canonical Kerr tortoise + inverse (a=0 ⇒ tortoise_schw)

include("Geodesics/frequencies.jl")
export kerr_frequencies, radial_roots, polar_roots, proper_frequency_factor

include("Geodesics/trajectory.jl")
export geodesic_trajectory, geodesic_trajectory_radial, geodesic_trajectory_polar
export geodesic_trajectory_equatorial, geodesic_proper_time

include("Separatrix/separatrix.jl")
export separatrix, kerr_isco, get_separatrix

include("Geometry/christoffel.jl")
export christoffel, christoffel_dr, christoffel_dtheta

include("Geometry/schwarzschild_orbit.jl")
export schwarzschild_lapse, schwarzschild_dlapse, bls_potential, tortoise_schw, tortoise_schw_inverse
export schw_orbit_energy, schw_orbit_ang_mom, schw_orbit_phi_freq

include("PNKinematics/pn_kinematics.jl")
export kepler_omega, r0_of_omega, pn_x, kepler_Ω, r0_of_Ω

# Circular equatorial u^t/u^φ closed forms + the Kerr azimuthal chart offset φ̄−φ and its exact
# r-derivative a/Δ (elliptic-SSF W1-C, arXiv:2606.04998 main.tex:192-196, 223-226). Additive.
include("Geodesics/circular_kinematics.jl")
export kerr_circular_ut_uphi, kerr_varphi_offset, kerr_varphi_offset_dr

# Page–Thorne-form potentials (ν, ψ, μ, ω̃) of the Kerr EQUATORIAL line element + analytic
# r-derivatives + the named circular dt/dτ accessor (EGD spiral-wave pipeline, arXiv:2608.17003
# eq (4) + A1a–e, main.tex:673). Over delta_kerr / kerr_circular_ut_uphi. Additive.
include("Geometry/equatorial_potentials.jl")
export kerr_equatorial_potentials, kerr_equatorial_potentials_dr, kerr_circular_dtdtau

# Linear-in-spin secondary corrections, namespaced `Spacetime.SpinningSecondary.*` to avoid the
# root-package `SpinningOrbit` module-name collision (design §4.1). Not re-exported to top level.
include("SpinningSecondary/spinning_secondary.jl")

# NP/GHP coefficient kernel (NPGHP campaign Move B), namespaced `Spacetime.NPGHP.*` — its NP/GHP
# symbols (ρ, τ, β, μ, pre_edth, …) would collide flat, so it is a nested module (SpinningSecondary
# precedent) and only the module NAME is exported, not its symbols.
include("NPGHP/NPGHP.jl")
export NPGHP

# Rácz–Tóth horizon-penetrating hyperboloidal PURE chart maps (F5 Part B; re-homed from
# Perturbation.TimeDomain.rtchart, MINIMAL move — only the 4 pure maps, no s/ℓs/coefficients),
# namespaced `Spacetime.RTChart.*`; re-imported+re-exported by Perturbation.TimeDomain (same objects).
include("Geometry/rt_charts.jl")
export RTChart

end # module Spacetime

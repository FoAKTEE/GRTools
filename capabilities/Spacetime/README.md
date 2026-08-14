# Spacetime

Kerr background-geometry + orbital-mechanics **foundation** for GrTool (Plan B, design
§4.1). This is the Tier-2 spine root: its entire dependency closure is `{LinearAlgebra,
NumericsCore}`. Every physics capability above it (Perturbation, MetricReconstruction,
SelfForce, EMRIWaveforms, PostNewtonian, Hybrid, AccretionDisc) reaches Kerr geometry,
geodesics, separatrix and ISCO through this single package — so the spine root carries
**no heavy dependency** (no `Roots`, no HDF5/Symbolics/Wigner/external-ODE).

## What physics this capability owns

The geometry of the Kerr background and the mechanics of geodesic + spinning-secondary
motion within it:

| Submodule | Owns |
|-----------|------|
| `Geometry` | BL/In/Out metric + inverse, Kinnersley tetrad, Christoffels (+ r/θ derivs), bg-scalar bundle Δ/Σ/ζ/K/Q, tortoise + inverse, horizons, lapse/potential, circular-orbit E/Lz/Ωφ. Coords are **methods** (BL/In/Out + `:eq` restriction), and the a=0/M-only Schwarzschild Christoffel is an M-only method of the canonical `christoffel`. |
| `Geodesics` | E/Lz/Q constants of motion, radial/polar roots, Fujita–Hikida Mino + BL frequencies Ω_r/Ω_θ/Ω_φ, closed-form trajectories, proper time. |
| `Separatrix` | ONE canonical `separatrix(a,e,x)` + `kerr_isco(a)` over NumericsCore's vendored A42/brent bracketing root-finder (no `Roots`). |
| `SpinningSecondary` | Linear-in-spin orbital actions/frequencies + spin-curvature force terms. |
| `PNKinematics` | Numeric Kepler Ω / `r0_of_Ω` / `pn_x` and the collapsed single ISCO/PN-ISCO fork (the *numeric* kinematic helpers — not the symbolic PN series). |
| `NPGHP` | NP/GHP **coefficient kernel** (NPGHP campaign Move B, nested `Spacetime.NPGHP.*`): the full Kerr Kinnersley NP spin-coefficient table `np_spincoeff(a,r,θ)` + GHP-primed view `ghp_spincoeff`, GHP `{p,q}=(s+b,−s+b)` weights, and the TeukEvolution HPHC-cl chart (`R=cl²/r`) with its 6 background NP_0 scalars + 11 GHP operator pre-coefficients. Owns **coefficients** (pointwise closed-form geometry); operator *application* (þ/þ′/ð/ð′) stays consumer-side. |

Everything is type-generic `{T}`, so BigFloat precision (256…2048-bit) is inherited from
the caller via `NumericsCore.precision`.

## Legacy packages absorbed

Consolidated here, with their original packages kept GREEN as the characterization oracle:

- **KerrGeodesics** — constants of motion, frequencies, radial/polar roots, trajectories (→ `Geodesics`).
- **SpinningOrbit** — linear-in-spin orbit corrections + spin deviations (→ `SpinningSecondary`).
- **MetricReconstruction / DiscModes / KerrLorenzCirc** geometry forks — `kerr_metric`/`kinnersley`/`christoffel`/bg-scalars canonicalized as coords-parametrized methods (→ `Geometry`).
- **Punctures** `worldline_self_field` `_christoffel(r,θ,M)` — a=0/M-only Schwarzschild method of the canonical `christoffel` (→ `Geometry`).
- **KerrBackwards `separatrix` / FastEMRIWaveforms `get_separatrix` / WaSABI `kerr_isco` / DiscModes `_risco` / DiscEMRIForce `_risco_disc`** — collapsed to one `separatrix` + one `kerr_isco`, the legacy names become thin aliases (→ `Separatrix` / `PNKinematics`).

## Public API (built incrementally, tests-first)

`kerr_metric`, `kerr_metric_inv`, `kinnersley`, `christoffel`, `bg_scalars`, `tortoise`,
`horizons`, `circular_orbit` · `constants_of_motion`, `kerr_frequencies`, `radial_roots`,
`polar_roots`, `kerr_trajectory`, `kerr_proper_time` · `separatrix`, `kerr_isco`,
`get_separatrix` · `spin_orbit`, `spin_frequencies`, `spin_deviations` · `kepler_omega`,
`r0_of_omega`, `pn_x` · `NPGHP.{np_spincoeff, ghp_spincoeff, ghp_weights, np_spincoeff_hphc,
hphc_ghp_coeffs, hphc_R_of_r, hphc_r_of_R, hphc_horizon_R}` (nested module, not re-exported flat).

(Scaffold stage: submodules are stubs; each lands grounded against the legacy oracle or an
analytic identity before its exports go live.)

## Dependencies

- **Runtime:** `LinearAlgebra`, `NumericsCore` (Carlson/elliptic/Jacobi for Geodesics;
  precision policy; the vendored A42/brent bracketing root-finder for Separatrix). No
  external deps — spine-root invariant.
- **Test-characterization oracles (light, live):** `KerrGeodesics` (empty `[deps]`) and
  `SpinningOrbit` (depends only on KerrGeodesics). No heavy package is ever a live test-dep
  here, keeping the cold-compile fast.

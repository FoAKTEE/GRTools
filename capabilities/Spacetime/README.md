# Spacetime

Kerr background-geometry + orbital-mechanics **foundation** for GrTools. This is the
Tier-2 spine root: its entire dependency closure is `{LinearAlgebra, NumericsCore}`. Every
physics capability above it (Perturbation, MetricReconstruction, SelfForce, EMRIWaveforms,
PostNewtonian, Hybrid, AccretionDisc) reaches Kerr geometry, geodesics, separatrix and ISCO
through this single package — so the spine root carries **no heavy dependency** (no `Roots`,
no HDF5/Symbolics/Wigner/external-ODE).

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
| `NPGHP` | NP/GHP **coefficient kernel** (nested `Spacetime.NPGHP.*`): the full Kerr Kinnersley NP spin-coefficient table `np_spincoeff(a,r,θ)` + GHP-primed view `ghp_spincoeff`, GHP `{p,q}=(s+b,−s+b)` weights, and the TeukEvolution HPHC-cl chart (`R=cl²/r`) with its 6 background NP_0 scalars + 11 GHP operator pre-coefficients. Owns **coefficients** (pointwise closed-form geometry); operator *application* (þ/þ′/ð/ð′) stays consumer-side. |

Everything is type-generic `{T}`, so BigFloat precision (256…2048-bit) is inherited from
the caller via `NumericsCore.precision`.

## Provenance

Consolidated here from the predecessor codes, which served as the characterization oracle during
construction (their outputs are baked into `test/fixtures/`):

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

## Equatorial eccentric harmonics

`EquatorialOrbit(a,p,e)` stores a stable bound orbit with positive Lz, pericentre
at q=0, and M=1. A negative spin describes retrograde motion relative to the hole.
`worldline_rates(orbit,q)` returns signed `Ur=dr/dλ`, the periodic t/φ deviations,
the coordinate-time measure `dt_dlambda`, and proper-time four-velocity components.
`orbit_frequency(orbit,m,k)` is mΩφ+kΩr.

```julia
using Spacetime
setprecision(BigFloat, 256) do
    orbit = EquatorialOrbit(big"0.9", big"10", big"0.3")
    H(q, state) = state.r / orbit.p + im * state.Ur / orbit.p
    coefficient = orbit_harmonic(H, orbit, 2, -1; N=160, paired=true)
end
```

The callback supplies the azimuthal coefficient H_m, with exp(−imψ) convention.
The projector integrates H_m exp(i[kq+ωΔt−mΔφ]) with the dt/dλ measure and
normalizes by 2πΓ. `paired=true` evaluates both q and 2π−q explicitly; it does
not assume a real or even callback. Nodes are generated in the orbit's precision.
Check convergence by increasing N for the orbit and harmonic of interest.

`turning_point_jacobian` and `orbital_params_jacobian` live in Geodesics;
`SpinningSecondary.spin_orbital_params_jacobian` is the compatibility wrapper.
The constants differentiated are (E,Lz,K), where K=Q+(Lz−aE)². The legacy
Jacobian keeps its original arithmetic and signed-spin constants convention:
physical inverse checks use a≥0 with signed x. A consumer of negative-a element
maps must normalize that chart explicitly. Circular turning-point degeneracy
(e=0) is outside this Jacobian's domain.

The equatorial trajectory fixture is independently generated from Darwin-angle
quadrature; [its provenance](test/fixtures/spacetime_trajectory_PROVENANCE.md)
contains the equations and reproduction command. C-ORBIT/C-HARM checks live
in `test/test_eccentric_harmonics.jl`; C-MAP chain checks run in the same
Spacetime test environment in `test/test_turning_point_jacobian.jl`.

## Dependencies

- **Runtime:** `LinearAlgebra`, `NumericsCore` (Carlson/elliptic/Jacobi for Geodesics;
  precision policy; the vendored A42/brent bracketing root-finder for Separatrix). No
  external deps — spine-root invariant.
- **Characterization oracles:** the predecessor `KerrGeodesics` and `SpinningOrbit` codes,
  whose outputs on the test grids are baked into `test/fixtures/`. No heavy package is ever a
  live test-dep here, keeping the cold-compile fast.

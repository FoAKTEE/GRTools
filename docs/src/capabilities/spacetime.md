# Spacetime

Kerr background-geometry + orbital-mechanics **foundation** of GrTool (design §4.1) — the
Tier-2 *spine root*. Its entire dependency closure is `{LinearAlgebra, NumericsCore}`; every
capability above it (Perturbation, MetricReconstruction, SelfForce, EMRIWaveforms,
PostNewtonian, Hybrid, AccretionDisc) reaches Kerr geometry, geodesics, separatrix and ISCO
through this single package, so the spine root carries **no heavy dependency** by
construction.

The package owns the physics that is *geometry of the Kerr background and motion within it*:
metric/tetrad/connection, bound-geodesic constants of motion and Fujita–Hikida frequencies,
closed-form trajectories, the separatrix/ISCO boundary of bound motion, linear-in-spin
(Mathisson–Papapetrou–Dixon) secondary corrections, and the numeric post-Newtonian kinematic
helpers (Kepler ``\Omega``, ``r_0(\Omega)``, ``x = (M\Omega)^{2/3}``). Everything is
type-generic `{T}`, so BigFloat precision (256…2048-bit) is inherited from the caller via
`NumericsCore.precision` — nothing is hard-Float64'd except where noted below.

## Submodule map

Source is organized by physics area under `src/`. Note the namespace reality: the first four
areas are `include`d **directly into the flat `Spacetime` namespace** (they are source
directories, not nested Julia modules); only `SpinningSecondary` is a real nested module.

| Source area | Owns | Namespace |
|---|---|---|
| `Geometry/` | BL Kerr metric + analytic inverse, Kinnersley tetrad, Christoffels ``\Gamma^\lambda_{\mu\nu}`` (+ the a=0/M-only Schwarzschild method and its r/θ derivatives), bg-scalar bundle ``\Delta,\Sigma,\zeta,K,Q``, horizons + ``\Omega_+``, the canonical Kerr tortoise coordinate `tortoise_kerr` + `tortoise_kerr_inverse`, Schwarzschild lapse/potential/tortoise + circular-orbit ``E, L_z, \Omega_\varphi`` | `Spacetime.*` |
| `Geodesics/` | ``E/L_z/Q`` constants of motion, radial/polar roots, Fujita–Hikida Mino + BL frequencies ``\Omega_r,\Omega_\theta,\Omega_\varphi``, closed-form trajectories, proper time | `Spacetime.*` |
| `Separatrix/` | ONE canonical `separatrix(a,e,x)` + `kerr_isco(a)` over NumericsCore's vendored A42/brent bracketing root-finder (no `Roots`); `get_separatrix` is a thin alias | `Spacetime.*` |
| `PNKinematics/` | numeric `kepler_omega` / `r0_of_omega` / `pn_x` (+ the `kepler_Ω`/`r0_of_Ω` aliases) — the *numeric* kinematic helpers, **not** the symbolic PN series (those live in PostNewtonian) | `Spacetime.*` |
| `SpinningSecondary/` | linear-in-spin orbital constants/frequencies, spin-curvature worldline deviations, the numeric 2-torus force/orbit-correction engine (`spin_orbit_correction` / `spin_orbit_numeric`) | `Spacetime.SpinningSecondary.*` (namespaced, **not** re-exported to top level, to avoid the legacy `SpinningOrbit` module-name collision — design §4.1) |

Coordinate charts are **methods**, not flattened forks: `kerr_metric(a,r,θ; coords=:BL/:In/:Out)`
plus the `:eq` restriction, and the Schwarzschild connection coexists with the Kerr one by
arity — `christoffel(a,r,θ; coords=:BL)` (full Kerr) vs `christoffel(r,θ; M=1)` (M-only,
rehomed from Punctures).

## Dependencies and poison-axis budget

Runtime `[deps]` are exactly `LinearAlgebra` + `NumericsCore` (Carlson/elliptic/Jacobi for
Geodesics, precision policy, and the vendored A42 bracketing solver for Separatrix — `Roots`
was dropped from the spine, design §4.1 Rev-2).

This budget is a **binding CI gate**, not an aspiration: `test/hygiene.jl` (§7.7 gates 1+2)
asserts `allow = ["LinearAlgebra", "NumericsCore"]` against `Project.toml` and forbids all
four poison axes plus the heavy tail — `Roots`, `HDF5`/`JLD2`, `Symbolics`,
`OrdinaryDiffEq*`, `CUDA`/`AMDGPU`/`KernelAbstractions`, `WignerSymbols`, `SpecialFunctions`,
`HypergeometricFunctions`, `NLsolve`, `SciMLBase`, `QuadGK` — and additionally checks after
`using Spacetime` that none of them appear in `Base.loaded_modules` (closure-load gate).
`TestKit` (+ its `Data` spine) is a **test-only** dependency powering those gates; since
Phase 7-B the legacy `KerrGeodesics`/`SpinningOrbit` test-oracles are baked into
`test/fixtures/` and no `packages/*` entry remains in the test target.

## GPU / data needs

None / none. The package is closed-form scalars, 4×4 linear algebra and per-orbit
root-finds; parallelism is caller-side CPU broadcast over ``(a,p,e,x)``. Zero data files —
this is load-bearing for the spine-root property. The one Float64-tier corner is the
`SpinningSecondary` DFT/least-squares engine (LAPACK `pinv`, faithful to the oracle); it is
serial, so no §8 BLAS guard is required.

## Legacy provenance (drawsFrom)

Functions moved here by physics home, regardless of source package (design §6.1). The legacy
packages were kept green as characterization oracles during the port; their outputs are now
frozen fixtures.

| Legacy package / file / function(s) | Landed as |
|---|---|
| KerrGeodesics `kerrgeo_constants.jl` (`KerrGeoEnergy`, `AngularMomentum`, `CarterConstant`, `ConstantsOfMotion`) | `constants_of_motion`, `geodesic_energy`, `geodesic_angular_momentum`, `carter_constant` |
| KerrGeodesics `kerrgeo_frequencies.jl` (`KerrGeoRadialRoots`, `KerrGeoFrequencies`, `ProperFrequencyFactor`) | `kerr_frequencies`, `radial_roots`, `polar_roots`, `proper_frequency_factor` |
| KerrGeodesics `kerrgeo_trajectory.jl` (`KerrGeoTrajectory` + Equatorial/Radial/Polar, `KerrGeoProperTime`) | `geodesic_trajectory*`, `geodesic_proper_time` |
| KerrBackwards `utility.jl` `separatrix`/`_sep_root` (was the `Roots` carrier) | `separatrix` (canonical, over NumericsCore A42) |
| FastEMRIWaveforms `separatrix.jl` `get_separatrix`/`_brentq` | `get_separatrix` (thin alias of `separatrix`) |
| WaSABI `kerr_isco` + DiscModes `_risco` + DiscEMRIForce `_risco_disc` + FEW kerrflux isco (4 forks) | one `kerr_isco` |
| MetricReconstruction `background.jl` + DiscModes `kerr_geometry.jl` + KerrLorenzCirc `operators.jl` (metric/tetrad forks) | `kerr_metric`, `kerr_metric_inv`, `kinnersley_tetrad` (coords-parametrized methods) |
| MetricReconstruction `background.jl`/KerrLorenzCirc (``\Delta,\Sigma,\zeta,K,Q`` bundles) | `bg_scalars`, `delta_kerr`, `sigma_kerr`, `zeta_kerr`, `kfun`, `qfun` |
| MetricReconstruction `{einstein,weyl_scalars}.jl` `christoffel_BL` + DiscModes | `christoffel(a,r,θ)` (Kerr BL) |
| Punctures `worldline_self_field.jl` `_christoffel(r,θ,M)` (+dr, +dθ) | `christoffel(r,θ; M)`, `christoffel_dr`, `christoffel_dtheta` (a=0/M-only method — geometry is Spacetime physics, not NumericsCore; §4.5 Rev-2) |
| MetricReconstruction `weyl_amplitudes.jl` `_horizons`, `_Omega_plus` | `horizons`, `omega_plus` |
| H1Lorenz `orbit.jl` (`rs_of_r`, `r_of_rs`, `_f`, `_df`, `_Vl`, `orbit_*`) + DiscModes `lorenz_metric.jl` | `schwarzschild_lapse/dlapse`, `bls_potential`, `tortoise_schw(_inverse)`, `schw_orbit_energy/ang_mom/phi_freq` |
| SpinningOrbit `linear_parts.jl`, `spin_orbit.jl`, `spin_deviations.jl`, `numeric_correction*.jl`, `orbit_correction.jl` | `Spacetime.SpinningSecondary` (whole module) |
| WaSABI `pn_kinematics.jl` (`kepler_Ω`, `r0_of_Ω`, `pn_x` + aliases) | `kepler_omega`, `r0_of_omega`, `pn_x` (+ `kepler_Ω`, `r0_of_Ω`) |

## Spinning secondary: the S5-1 spin-curvature slice

The S5-1 feature port brought the full linear-in-spin **force machinery** for generic
(inclined) bound orbits into `Spacetime.SpinningSecondary` — beyond the earlier
constants/frequencies slice:

- the 4 spin-curvature worldline deviations ``\delta t_3, \delta r_3, \delta z_3,
  \delta\varphi_3`` (the ``\mathcal{O}(s)`` position shifts off the geodesic, evaluated on
  the geodesic trajectory);
- the 47 algebraic force/coefficient terms + Darwin-angle maps + 2-torus
  ``(\lambda_r,\lambda_\theta)`` DFT `Correction` engine (`_numeric_correction`): the
  spin-curvature force is Fourier-decomposed over the 2-torus, then the orbit correction is
  solved linearly via an M-matrix least-squares (`pinv`) solve;
- the 11 `_oc_*` four-velocity/metric helpers and the public assembly:
  [`Spacetime.SpinningSecondary.spin_orbit_correction`](@ref) (spin-independent correction
  data at ``(a,p,e,x)`` with resolution `(nmax, kmax)`, `kmax` even) and
  [`Spacetime.SpinningSecondary.spin_orbit_numeric`](@ref) (the thin frequency assembly
  ``\{\Gamma,\Upsilon_r,\Upsilon_\theta,\Upsilon_\varphi\} =
  \text{geo} + s\,\cdot\text{correction}``, then BL frequencies ``\Upsilon/\Gamma``).

The algebraic terms are pure arithmetic and type-generic `{T}`; only the DFT/least-squares
solve is Float64-tier (faithful to the oracle). Scope boundary (not a missing oracle): the
closed-form **inclined** `spin_orbit` branch (``x^2 \neq 1``) and the **spherical**
(``e = 0``) engine remain in the legacy oracle only — but the numeric `Correction` engine
here *is* the generic inclined force machinery, evaluated at ``x^2 < 1``. `spin_orbit`
itself throws for ``x^2 \neq 1``.

## Public API

Top-level `Spacetime` exports (flat namespace), then the namespaced
`Spacetime.SpinningSecondary` exports. Internal `_`-prefixed helpers are intentionally
omitted.

```@autodocs
Modules = [Spacetime, Spacetime.SpinningSecondary]
Order = [:module, :type, :constant, :function]
Private = false
```

## What changed vs legacy

- **One separatrix, one ISCO.** `separatrix`/`get_separatrix` (KerrBackwards/FEW) collapse
  to a single implementation; the four `kerr_isco`/`_risco` forks collapse to one. The
  external `Roots` dependency was dropped — the A42 (Alefeld–Potra–Shi) bracketing solver is
  vendored in `NumericsCore.rootfind`, characterization-gated at frozen tolerance against
  the legacy Roots-based oracle.
- **Geometry canonicalized as methods.** The MetricReconstruction/DiscModes/KerrLorenzCirc
  metric/tetrad/Christoffel forks became coords-parametrized methods of one function each;
  Punctures' Schwarzschild `_christoffel(r,θ,M)` became the M-only method of the canonical
  `christoffel` (not a separate name, and not in NumericsCore).
- **The Kerr tortoise coordinate is canonical here** (DG-P2, commit `8be5f32` — closes the
  Phase-7D "no canonical `tortoise_kerr`" leftover): `tortoise_kerr(r, a)` /
  `tortoise_kerr_inverse(rs, a)` in `Geometry/`, M = 1 convention, reusing
  `horizons`/`delta_kerr`; the ``a = 0`` branch delegates to `tortoise_schw` **bit-exactly**
  (`==`-gated). The inverse is a Newton solve on the ``r_*`` residual with an
  `eps(T)`-scaled horizon clamp, so it is BigFloat-inherited (round-trips gated to 1e-12 at
  Float64 and 1e-40 at BigFloat-256, ``a \in \{0, 0.5, 0.9, 0.998\}``, plus 4 baked DG-repo
  characterization values at rtol 1e-14). `Perturbation.TimeDomain` consumes `tortoise_kerr`
  from here (its *inverse* stays local by documented deviation — a different deep-horizon
  floor convention), while the `:GSN` backend keeps a local map textually faithful to its
  upstream export; both are gated against the same anchors.
- **Naming.** Legacy `KerrGeo*` CamelCase verbs became snake_case (`KerrGeoEnergy` →
  `geodesic_energy`, etc.); KerrBackwards' `constants_of_motion` survives as the delegating
  wrapper name. `SpinningOrbit.KerrSpinOrbit*` became
  `Spacetime.SpinningSecondary.spin_orbit*` — qualified access only, never re-exported flat.
- **Precision.** Legacy Float64 hard-coding removed: all closed-form/algebraic paths are
  `{T}`-generic and BigFloat-inherit; the sole Float64-tier remnant is the SpinningSecondary
  DFT/`pinv` engine, kept faithful to the oracle.
- **Tests stand alone.** Since Phase 7-B the characterization oracles are baked fixtures
  under `test/fixtures/`; the test target has no live `packages/*` dependency.

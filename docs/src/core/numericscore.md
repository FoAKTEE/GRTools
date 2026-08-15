# NumericsCore (`GrToolsBase`)

The **dimensionless math + infrastructure sink** of GrTools — tier C0, the unique root
of the dependency DAG (design doc §3.2/§5.1). Its only runtime dependency is
`LinearAlgebra`, and a binding CI dep-budget gate keeps it that way. It carries **no
physics geometry**: Kerr metric/Christoffel/geodesics are `Spacetime` physics, and even the
Schwarzschild `M`-only `_christoffel` from Punctures landed in `Spacetime.Geometry`, not
here.

Everything in this package is type-generic `{T}`, so BigFloat precision (256/512/1024/2048
bits) is *inherited*, never hard-`Float64`'d on a generic path. Physicist mental map:

- **Special functions** — Carlson symmetric integrals ``R_F, R_C, R_D, R_J``; complete and
  incomplete Legendre elliptic ``K(m), E(m), \Pi(n,m)``; Jacobi amplitude and
  ``\mathrm{sn}, \mathrm{cn}, \mathrm{dn}`` (real and complex arguments).
- **Legendre / angular seed** — normalized associated-Legendre ``P_\ell^m`` with
  θ-derivative, plus the half-angle term-list machinery behind the spin-weighted spherical
  ``{}_sY_{\ell m}`` seed that `Harmonics` builds on.
- **Interpolation library** — not-a-knot cubic/bicubic/tricubic, GSL-style natural bicubic,
  barycentric and log-Chebyshev, Chebyshev-spectral differentiation, Fornberg
  finite-difference weights, cubic Hermite, and WL-faithful `Interp1/2/3`. All with a
  coefficient/eval split (coefficients built once, evaluation batchable — the GPU-eval seam).
- **Mode indexing** — the ``(l,m,n,k)`` `ModeIndexTable`, the flat ``\gamma = l(l+1)/2 + m``
  counter, ladder factors ``\hat\lambda, \hat\Lambda``, Barack–Lousto ``a_{il}``.
- **Precision policy** — `bc_tol`/`rk_tol`/`fd_step` selectors that thread the element type
  (Float64 keeps the validated legacy values `1e-14`/`1e-12`; BigFloat scales as
  ``2^{10-\mathrm{bits}}`` / ``2^{20-\mathrm{bits}}``), `bigprec_for_mag`, `with_precision`.
- **Units & conventions** — split-named constants (`Gpc_in_m` vs `Gpc_in_s`, `MTSUN_SI`,
  `YRSID_SI`, ...); `p_to_y`, ``\epsilon \leftrightarrow \omega``, ``q \leftrightarrow a``,
  `imagunit`.
- **DG spectral-element core** (`dg.jl`, the DG-Teulolsky port) — 1D nodal
  discontinuous-Galerkin building blocks in Hesthaven–Warburton conventions:
  `RefElement(N; T)` (LGL nodes/weights, differentiation matrix `D`, mass matrix `M`),
  multi-domain `Mesh`/`build_mesh` with forced interfaces, `element_of`/`node_coords`,
  barycentric `lagrange_interp`, and the `rk4!`/`cfl_dt` time-stepping helpers. The
  original's `FastGaussQuadrature` dependency is **eliminated** — the LGL rule is computed
  self-contained (Newton on the roots of ``P'_N`` from Chebyshev–Gauss–Lobatto seeds), which
  also makes the operators type-generic: `RefElement(N; T=BigFloat)` yields genuine BigFloat
  spectral operators (~1e-77 residuals at 256 bits), a capability the Float64-only original
  lacked. Problem-agnostic — no Teukolsky physics enters; `Perturbation.TimeDomain` is the
  consumer.
- **Infra** — `blas_guard` (scoped `BLAS.set_num_threads`), the `AbstractBackend` /
  `CPUBackend` / `select_backend` dispatch seam, vendored adaptive ODE integrators
  (DP5, GBS, DOPR853), continued-fraction primitives (finite / convergent / modified-Lentz),
  a self-contained series Kummer ``{}_1F_1``, root-finders (secant, bisection, Brent, and the
  A42 Alefeld–Potra–Shi bracketing solver), and dep-free `DataSchema` contract structs.

## Submodule map

`NumericsCore` is a single flat module; source files map one-to-one onto the capability list
above:

| Source file | Owns |
|---|---|
| `specialfun_carlson.jl` | Carlson `RF`/`RC`/`RD`/`RJ`, elliptic `ellipK`/`ellipE`/`ellipPi` (+ incomplete), Jacobi `jacobiAM`/`SN`/`CN`/`DN` |
| `precision.jl` | `bc_tol`, `rk_tol`, `fd_step`, `bigprec_for_mag`, `with_precision` |
| `modeindex.jl` | `gamma_to_lm`, `lm_to_gamma`, `mode_list`, `ModeIndexTable`, `a_il`, `lambdahat`/`Lambdahat` |
| `units.jl` | `MTSUN_SI`, `YRSID_SI`, `MRSUN_SI`, `AUsec`, `Gpc_in_m`, `Gpc_in_s`, `GPCINSEC`, `PI`/`PI_2`/`PI_4` |
| `conventions.jl` | `p_to_y`, `eps_to_omega`/`omega_to_eps`, `q_to_a`/`a_to_q`, `imagunit` |
| `blas_guard.jl` | `with_blas_threads`, `blas_single_thread` |
| `rootfind.jl` | `secant`, `bisect`, `brent`, `a42` |
| `continued_fraction.jl` | `cf_finite`, `cf_convergent`, `cf_lentz` |
| `dataschema.jl` | `H1Schema`, `AmpGridSchema`, `FluxGridSchema`, `validate_schema`, `schema_version` |
| `interp.jl` | the full interpolation library (see overview) |
| `ode.jl` | DP5 (`dp5_step!`, `integrate_to_targets!`), GBS (`gragg_midpoint!`, `gbs_*`), and the nested **`NumericsCore.Dopr853`** submodule (`DOPR853`, `take_step`, `integrate`, `dense_eval`) |
| `legendre.jl` | `assoc_legendre_P`, `legendre_P_thetaderiv`, `Ylm_real`, `Ylm_eq`, `spin_weighted_Y` |
| `kummer.jl` | `kummer_1F1` (self-contained Maclaurin series), `kummer_U` (Γ-free asymptotic Tricomi `U`) |
| `dg.jl` | `RefElement`, `Mesh`, `build_mesh`, `element_of`, `node_coords`, `lagrange_interp`, `rk4!`, `cfl_dt` (nodal-DG spectral-element core, self-contained LGL) |
| `backend_iface.jl` | `AbstractBackend`, `CPUBackend`, `select_backend`, `backend_array`, `backend_collect` |

The only genuine nested module is `NumericsCore.Dopr853` — the Dormand–Prince 8(5,3)
trajectory integrator with degree-7 dense output, ported verbatim (Float64 tableau) from
FastEMRIWaveforms; its exports are re-exported at top level.

## Dependencies & poison-axis declaration

Runtime `[deps]`: **`LinearAlgebra` — exactly, and nothing else.**

This is not an aspiration but a CI gate: `test/hygiene.jl` runs
`TestKit.assert_dep_budget` with `allow = ["LinearAlgebra"]` and explicitly **forbids**
the full poison-axis set — HDF5, JLD2, Symbolics, `OrdinaryDiffEq*`, CUDA/AMDGPU/Metal/
oneAPI/GPUArrays/KernelAbstractions, WignerSymbols, Roots, NLsolve, SciMLBase, QuadGK,
ForwardDiff, SpecialFunctions, HypergeometricFunctions — **and every other GrTools
package** (nothing may point *up* from the sink). A second closure-load gate
(`TestKit.assert_not_loaded`) asserts none of the poison packages sneak into the load
closure after `using NumericsCore`.

Consequences physicists rely on:

- depending on NumericsCore inherits **zero** heavy deps — this is what keeps the
  `Spacetime`/`Harmonics` spine roots at closure exactly `{LinearAlgebra, NumericsCore}`;
- the self-contained `kummer_1F1` here is *distinct by decision* from
  `BigSpecial.kummer_hyp` (HypergeometricFunctions-delegated) — the heavy special-function
  stack never leaks onto the spine;
- vendoring A42 here is what let `Spacetime` drop the external `Roots` package.

Test-only: `Test`, `TestKit` (+ its `Data` spine) power the hygiene gates. All legacy
characterization oracles are **baked fixtures** (post Phase-7B) — see below.

## GPU / data needs

None. NumericsCore is CPU-only and file-IO-free:

- `backend_iface.jl` defines the dispatch *seam* only; the GPU backend type and the
  `gpu=true` branch of `select_backend` are injected by the `Kernels` package extension.
  No CUDA/AMDGPU/KernelAbstractions symbol is referenced here.
- `dataschema.jl` holds pure struct shapes + field-invariant validators; the HDF5/JLD2
  reader/writer *bodies* live in `Data` weakdep-extensions. The grid layout is a versioned
  contract decoupled from the storage backend.
- The interpolation coefficient/eval split is deliberate: coefficients are built on the
  (possibly BigFloat) CPU tier; batched Float64 evaluation is the seam `Kernels` consumes.

## Legacy provenance (drawsFrom)

Primitives were consolidated — one canonical form per *same decision* — from thirteen legacy
packages. Each ported function was grounded by characterization against its legacy oracle at
the legacy tolerance, or against an analytic identity (interpolants are exact on the
polynomial degree they interpolate). Where the ported function lives now:

| Legacy package / file | What landed here |
|---|---|
| KerrGeodesics (`_RF`/`_RC`/`_RD`/`_RJ`, `_ellip*`, `_jacobi*`) | Carlson + elliptic + Jacobi family (`specialfun_carlson.jl`), made type-generic |
| SpinWeightedSpheroidalHarmonics (`SpinWeightedSphericalHarmonicY`, `_y_dtheta`; `leaver.jl`) | spin-weighted `spin_weighted_Y` seed; self-contained Kummer ``{}_1F_1`` pattern + continued-fraction primitives (promoted per design H-CORE-STRAT) |
| Teukolsky (`continued_fraction.jl` `_cf_finite`/`_cf_convergent`; `radial_numint.jl` `_hermite`, `_secant`) | `cf_finite`/`cf_convergent` (+ new `cf_lentz`); `cubic_hermite`; `secant` |
| H1Lorenz (`math_funcs.jl`, `driver.jl`, `integrator.jl`, `io_hdf5.jl`) | `assoc_legendre_P` + θ-derivative; `gamma_to_lm`/`lm_to_gamma`/`mode_list`; `bc_tol`/`rk_tol` closed forms; DP5/GBS integrators; `blas_guard`; `H1Schema` field layout |
| SecondOrderRicci (`modes.jl`) | Barack–Lousto `a_il` |
| KerrBackwards (`interpolant.jl`, `gslbicubic.jl`, `constants.jl`, `p_to_y`) | `CubicSplineInterpolant`, `GslBicubic` + `natural_spline_node_derivs`/`gsl_bsearch`, canonical LAL constants, `p_to_y` |
| FastEMRIWaveforms (`ampinterp2d.jl`, `kerrflux.jl`, `dopr853.jl`, `baseclasses` mode table, `_brentq`/separatrix A42) | `BicubicNotAKnot`, `TricubicNotAKnot`, `AmpGridSchema`/`FluxGridSchema` layouts, `Dopr853` submodule, `ModeIndexTable`, `brent`/`a42` |
| WaSABI (`interpolation.jl`) | WL-faithful `Interp1WL`/`Interp2Hermite`/`Interp3Lagrange` |
| DiscModes (`integrator.jl`, `metric_recon.jl`) | DP5 vendoring dedup; the duplicate Condon–Shortley Legendre kernel folded into `_legendre_Plm` |
| DiscEMRIForce (`numerics.jl`) | `fornberg_weights`, `cheb_*`, `bary_*`, `local_interp` |
| KerrLorenzCirc (`radial_kerr.jl`, `operators.jl`) | `LogChebInterp`/`BigLogRInterp`, `fd1c`/`fd2c`; precision-policy plumbing |
| Punctures, SpinningOrbit | absorbed at the primitive level (shared conventions/precision plumbing); their *geometry* went to `Spacetime`, not here |
| DG-Teulolsky (`src/DGCore.jl`, `src/TimeStepping.jl`) — *external research repo, not `packages/*`* | the `dg.jl` spectral-element core; characterized vs baked FastGaussQuadrature-original fixtures `test/fixtures/dg_lgl_ref.jl` (N = 2..16: node/weight rel dev ≤ 5.3e-15, D/M Frobenius ≤ 3.2e-15) + the original testsets ported verbatim (D exactness, mass row-sums, forced-interface mesh, advection superconvergence) |

The FEW-vs-KerrBackwards fork collisions (`CubicSplineInterpolant`, bicubic/tricubic,
`GslBicubic`, `ModeIndexTable`, constants, Carlson) are *same-decision* primitives and were
hoisted to the single copy here rather than namespaced (design §6/C3).

## What changed vs legacy

- **BigFloat inheritance everywhere.** Legacy KerrGeodesics hard-`Float64`'d the Jacobi
  functions; all Carlson/elliptic/Jacobi paths here promote to
  `T = promote_type(...)`. Root-finder default tolerances are `eps(T)`-scaled, so they
  tighten automatically at higher precision.
- **Units forks reconciled, not averaged.** KerrBackwards carries LAL `global.h` values;
  FEW re-derives from `lisaconstants` (differs at ``\sim 10^{-7}`` in `MTSUN`/`MRSUN`/
  `YRSID`). The LAL values are canonical here. The `Gpc` symbol collision (meters in FEW,
  seconds in the C `GPCINSEC`) is resolved by *splitting the axis*: `Gpc_in_m` vs
  `Gpc_in_s` (with `GPCINSEC` kept as a legacy alias).
- **One ODE family instead of four vendored copies** (Teukolsky, H1Lorenz, DiscModes, FEW),
  with step tolerance flowing through `rk_tol(T)`.
- **`Roots` dropped from the spine** — A42 is vendored in `rootfind.jl`; the
  separatrix/ISCO oracle was characterization-gated at frozen tolerance before the legacy
  `Roots` call was swept.
- **`FastGaussQuadrature` never entered** — the DG core's upstream LGL dependency was
  replaced by the self-contained Newton-on-``P'_N`` rule in `dg.jl`, keeping the sink budget
  at exactly `LinearAlgebra` *and* upgrading the spectral operators to BigFloat inheritance.
- **A documented misnomer kept for compatibility:** `legendre_P_thetaderiv` reproduces the
  legacy `H1Lorenz.associated_legendre_P_theta_deriv` byte-for-byte, which actually returns
  ``-\,dP_\ell^m/d(\cos\theta)`` — i.e. ``(1/\sin\theta)\,dP_\ell^m/d\theta``; the genuine
  colatitude derivative is ``\sin\theta`` times this value. See its docstring.
- **Test oracles are baked fixtures.** Post Phase-7B no legacy `packages/*` code is on the
  test load path: Carlson/Legendre/Jacobi and SWSH references are captured-once fixtures
  under `test/fixtures/` (each with a provenance header naming the legacy function and
  arguments that produced it); units are checked against exact legacy literals + analytic
  cross-identities; `bc_tol`/`rk_tol` against the documented H1Lorenz closed forms.

## Public API

```@autodocs
Modules = [NumericsCore, NumericsCore.Dopr853]
Order = [:module, :type, :constant, :function]
```

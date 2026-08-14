# NumericsCore (`GrToolBase`)

The **dimensionless math + infrastructure sink** of GrTool. It is the unique root of the
dependency DAG: its only dependency is `LinearAlgebra`. A CI dep-budget gate keeps
it from re-growing — no heavy deps (HypergeometricFunctions, SpecialFunctions, HDF5, Symbolics,
KernelAbstractions, ODE solvers) and **no physics geometry** (geometry is `Spacetime` physics).

## What it owns

Pure, convention-stable, widely-reused primitives, all type-generic `{T}` (BigFloat-inheriting):

- **Special functions** — Carlson symmetric integrals `RF`/`RC`/`RD`/`RJ`; complete & incomplete
  Legendre elliptic `K`/`E`/`Π`; Jacobi amplitude + sn/cn/dn (real & complex).
- **Legendre seed** — normalized associated-Legendre `P_l^m` + θ-derivative; canonical spin-weighted
  spherical-Y term machinery.
- **Interpolation library** — not-a-knot cubic/bicubic/tricubic, GSL-bicubic, barycentric/log-Cheb,
  Chebyshev-spectral, Fornberg-FD, cubic-Hermite, WL-faithful Interp1/2/3 (coeff/eval split).
- **Mode indexing** — `(l,m,n,k)` `ModeIndexTable`, ladder factors `λ̂`/`Λ̂`, Barack-Lousto `a_il`.
- **Precision policy** — BigFloat bit selection, `bc_tol`/`rk_tol`/`fd_step` selectors, `with_precision`.
- **Units & conventions** — split-named constants (`Gpc_in_m` vs `Gpc_in_s`, `MTSUN_SI`, `YRSID_SI`);
  `p_to_y`, ε↔ω, q↔a, imaginary-unit helpers.
- **Infra** — `blas_guard` (scoped `BLAS.set_num_threads`), `backend_iface`
  (`AbstractBackend`/`CPUBackend`/`select_backend`), vendored adaptive-RK (DP5/GBS/DOPR853),
  continued-fraction primitive, self-contained Kummer 1F1, root-finders
  (secant/bisect/brent + A42 Alefeld-Potra-Shi), dep-free `DataSchema` structs.

## Legacy packages absorbed

Primitives are consolidated (one canonical form per *same decision*) from:
KerrGeodesics, SpinWeightedSpheroidalHarmonics, Teukolsky, H1Lorenz, KerrBackwards,
FastEMRIWaveforms, DiscModes, KerrLorenzCirc, SecondOrderRicci, DiscEMRIForce, Punctures,
WaSABI, SpinningOrbit. Each ported function is grounded by characterization against its legacy
oracle (kept green) at the legacy tolerance, or against an upstream/analytic identity.

## Dependencies

- Runtime: `LinearAlgebra`.
- Test-only (`[extras]`, post Phase-7B): `Test`, `TestKit`, `Data` — the §7.7 hygiene gates.
  NO legacy packages/* remain in the test sandbox: every legacy characterization oracle
  (KerrGeodesics Carlson/elliptic/Jacobi, SWSH spin-weighted-Y, KerrBackwards/FEW units,
  H1Lorenz tolerances, FEW mode table) is BAKED into a captured-once fixture under
  `test/fixtures/`, each with a provenance header naming the legacy fn + args that produced it.

## Status

Scaffold. Primitives are TDD-ported in rounds (Carlson/elliptic, precision, mode-index, units,
conventions, blas_guard first; interpolation, vendored-ODE, continued-fraction, root-finders next).

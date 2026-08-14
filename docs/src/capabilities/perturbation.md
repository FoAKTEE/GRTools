# Perturbation

Tier-3 vacuum/test perturbed-field solver for the Kerr spacetime — the **Teukolsky-radial
engine** (frequency domain) plus a **discontinuous-Galerkin time-domain solver**
(`Perturbation.TimeDomain`), and the **light-dep firewall** of GrTool (design §4.3).

## Overview

`Perturbation` owns the homogeneous solutions of the radial Teukolsky equation and everything
built directly on them for a *vacuum / test* (point-particle) perturbation:

- **Homogeneous radial solutions** ``R^{\mathrm{in}}/R^{\mathrm{up}}`` behind one dispatcher
  ([`TeukolskyRadial`](@ref)) over the built backends: the MST (Mano–Suzuki–Takasugi)
  hypergeometric representation (`:MST`), a static ``\omega = 0`` closed-form ``{}_2F_1`` branch
  (`:Static`), confluent-Heun (`:HeunC`), the Sasaki–Nakamura form (`:SasakiNakamura`, riding
  a *vendored* DP5 integrator — never an external ODE package), and the Generalized
  Sasaki–Nakamura backend (`:GSN`, all spins × In/Up at real ω — see below).
- **Time-domain evolution** — the [`Perturbation.TimeDomain`](@ref perturbation-timedomain)
  nested submodule: a multi-domain nodal-DG hyperboloidal solver for the 1+1 m-mode Teukolsky
  equation (vacuum ``s \in \{0, \pm 1, \pm 2\}``; sourced circular-orbit ``s \in \{0, \pm 2\}``
  at ``a = 0``), riding `NumericsCore.dg` — see its own section below.
- **Renormalized angular momentum** ``\nu`` — the MST continued-fraction root
  ([`RenormalizedAngularMomentum`](@ref)) — and the MST minimal-solution coefficients ``f_n``,
  generated through the BigFloat 256..2048-bit tier.
- **Point-particle amplitudes** ``Z/C`` and the **energy / angular-momentum fluxes** from the
  point-particle source (`TeukolskyPointParticleMode{Circular,Eccentric,Spherical,Generic}`,
  [`energy_flux`](@ref), [`angular_momentum_flux`](@ref)), for circular / eccentric / spherical /
  generic orbits × spin weights.
- **MST scattering phase** ``e^{2i\delta}`` ([`scattering_phase`](@ref)) + the Starobinsky
  constant + the ``B^{\mathrm{refl}}/B^{\mathrm{inc}}`` reflection/transmission machinery — a
  strict dep-subset of the radial tier.
- **Spinning-secondary spin-fluxes** — the complete 12-entrypoint `TeukolskySpinMode*` surface
  (linear-in-spin ``s=-2`` GW flux modes for a spinning small body; see below).

Every primitive is type-generic `{T}`, so BigFloat precision (256/512/1024/2048 bits) is
inherited end-to-end: the working precision flows in through the argument types (e.g.
`T = float(a)` on the static path) and out as `Complex{T}`.

## Source layout (flat namespace + one nested submodule)

`Perturbation` is a flat Julia module with exactly **one nested submodule**: the frequency-domain
directories under `src/` are a physics grouping of the source files, *not* nested Julia modules,
while `TimeDomain.jl` declares the genuine nested module `Perturbation.TimeDomain` (verified:
the only `module` declarations in `src/` are `module Perturbation` and `module TimeDomain`).
Frequency-domain public names live directly in `Perturbation`; the time-domain surface lives in
(and is exported from) `Perturbation.TimeDomain`.

| Source group | Files | What it holds |
| --- | --- | --- |
| `Radial/` | `renorm_angmom.jl`, `mst.jl`, `hypergeom.jl`, `radial.jl`, `carrier.jl`, `heunc.jl`, `teukolsky_radial_static.jl`, `teukolsky_radial.jl`, `radial_sn.jl`, `gsn_{coordinates,transform,potentials,conversions,asymptotics,boundary,solve,backend}.jl` | ``\nu``-solve, MST coefficients, exact ``{}_2F_1``/``U`` MST terms, raw ``R^{\mathrm{in}}/R^{\mathrm{up}}`` + ``r``-derivatives, the callable `TeukolskyRadialFunction` carrier, the `:HeunC`/`:Static`/`:SasakiNakamura`/`:GSN` backends, and the public `TeukolskyRadial` dispatcher |
| `AmplitudesFluxes/` | `amplitudes.jl`, `fluxes.jl`, `pointparticle.jl` | MST asymptotic amplitudes (In/Up trans/inc/ref), Wronskian energy/``L_z`` fluxes from ``Z``, GHP source convolution → ``Z`` + the `TeukolskyPointParticleMode*` entrypoints |
| `Scattering/` | `scattering.jl` | `scattering_phase` (numerical MST + closed-form analytic ``\varepsilon``-series), `starobinsky_constant`, `mst_coefficient`, `mst_convergence`, `scattering_phase_table` |
| `SpinFluxes/` | `circular.jl`, `eccentric.jl`, `inclined.jl`, `fab_basis.jl`, `from_correction.jl`, `spin_mode_num.jl`, `spin_mode_correction.jl`, `spherical*.jl`, `darwin.jl`, `spinning_*_local.jl`, `entrypoints.jl` | the 12 public `TeukolskySpinMode*` entrypoints plus the *local* spinning-orbit kinematics (inclined / numeric-correction / spherical) that `Spacetime.SpinningSecondary`'s equatorial-only slice does not carry |
| `TimeDomain/` (nested module) | `hyperboloidal.jl`, `background.jl`, `common.jl`, `mode1d.jl`, `flux1d.jl`, `kerrmodes.jl`, `kerrflux.jl`, `spinweight.jl`, `rtchart.jl`, `spinmodes.jl`, `spinsource.jl` | the DG time-domain solver: hyperboloidal charts, mode potentials, single-mode / coupled-mode evolution, circular-orbit flux + self-force, generic-spin vacuum sector, sourced ``|s| = 2`` sector |

`_`-prefixed internals stay unexported; the public dispatchers and entrypoints are exported from
the top-level module (or from `Perturbation.TimeDomain` for the TD surface).

## The `:GSN` radial backend (Generalized Sasaki–Nakamura)

`TeukolskyRadial(s, l, m, a, ω; method=:GSN)` is the absorption of the upstream
`GeneralizedSasakiNakamura.jl` v0.8.5 (arXiv:2306.16469) as a native backend (commit
`4339845`) — it **closes the documented `radial_sn` boundary** ("`s = -2` Up only"): the full
``s \in \{0, \pm 1, \pm 2\}`` × {In, Up} homogeneous surface at real ω, self-normalized to unit
transmission via the GSN asymptotic ansatzes (orders 0–8, in-house Taylor arithmetic — no
`TaylorSeries` dependency), with the SN→Teukolsky conversion handled internally and ``\lambda``
taken from `Harmonics` (per-row ≤ 1e-11 vs the upstream SWSH). Redesigns vs upstream: the
`DifferentialEquations`/`AutoVern9` solver is replaced by the **vendored NumericsCore DP5**
(Riccati + linear forms, cross-checked < 1e-9) with **quintic-Hermite FSAL dense output**
(O(h⁶) — measured GSN-function error 1.38e-7 → 1.07e-9 vs the cubic-Hermite prototype);
**8 upstream dependencies are eliminated** (firewall intact — zero new `[deps]`). Everything is
type-generic `{T}`; for the BigFloat tier and the documented D-GBS deviation see
[Precision](../topics/precision.md). Grounded against a 706-row upstream oracle fixture net
(1e-10 trust floor, incl. 48 In/Up unit-transmission amplitudes ≤ 6.2e-10) with the `:MST`
backend as referee and a condition-aware Wronskian gate. **Boundary:** complex ω (QNM contours)
and the inhomogeneous GSN subsystem are quarantined (`@test_throws`-gated) — `[FUTURE]` GSN-3.

## [The TimeDomain submodule (DG time-domain solver)](@id perturbation-timedomain)

`Perturbation.TimeDomain` is the multi-domain nodal discontinuous-Galerkin **time-domain**
Teukolsky solver (the TD counterpart of the frequency-domain `Radial`/`AmplitudesFluxes`
stack), merged from the DG-Teulolsky research repo (arXiv:2307.01349 / arXiv:2503.11523
method) in commits `1745571` (s = 0), `350eea4` (generic-spin vacuum) and `f45dd71` (sourced
``|s| = 2``). The problem-agnostic DG core (LGL `RefElement`, `Mesh`, `rk4!`) is **reused from
`NumericsCore.dg`** — no DG numerics are duplicated here — and the Kerr background delegates to
`Spacetime` (`tortoise_kerr`, geodesics) and `Harmonics` (harmonic values). The firewall holds:
the TD stack added **zero new dependencies**.

**Two hyperboloidal charts, by spin (the chart story):**

- ``s = 0`` — the repo's **hyperboloidal-layer chart** (`hyperboloidal.jl` `Layer`; Paper 1
  §II.D): `Mode1DProblem` (Schwarzschild single-mode tails), `Flux1DProblem` (circular-orbit
  flux + self-force with characteristic-projected δ-flux injection), `KerrTailProblem` /
  `KerrFluxProblem` (Kerr m-mode ℓ-coupled tails and fluxes). This path is characterized
  **bit-exact in construction** against baked one-shot fixtures of the original repo
  (DG-Teulolsky@0a88967): Mode1D tails ≤ 7.2e-14 of peak, Flux1D ``\dot E`` ≤ 1.4e-12 rel,
  KerrModes a = 0.9 ≤ 1.5e-14, KerrFlux ≤ 7.1e-14 rel.
- ``s \ne 0`` — the layer chart's scri coefficients diverge (``\mu_s \sim 2s/r``), so the
  generic-spin sector uses the **Harms–Bernuzzi–Brügmann RT horizon-penetrating hyperboloidal
  chart** (arXiv:1301.1591, proven-regular coefficients for every s; `rtchart.jl`), with the
  Barack–Giudice (arXiv:1702.04204) 1+1 m-mode ℓ-coupled system (`spinweight.jl` closed-form
  ``\cos\theta``/``\sin^2\theta`` couplings ≡ Wigner-CG fixtures ≤ 1e-13) in
  `SpinModesProblem`/`solve_spin_modes`. The locked s = 0 layer-chart path is bit-untouched.

**Spin coverage (gated):**

- **Vacuum ``s \in \{0, \pm 1, \pm 2\}``**, ``a`` generic — Price-tail exponents gated against
  the arXiv:1301.1591 laws (e.g. s = −2, a = 0 scri measured −5.999 vs law −6; coupled a = 0.9
  ℓ-mixed tails), the a = 0, s = −2 operator ≡ the chain-rule Bardeen–Press operator, and
  boundedness for all combos. The s = +2 scri tail *exponent* is skipped with evidence (tail
  below the 7.5e-7 numerical floor); stability is still gated.
- **Sourced ``s \in \{0, \pm 2\}``, circular equatorial, ``a = 0``** (`spinsource.jl`):
  closed-form δ/δ′/δ″ jump conditions (Barack–Giudice App. A) enforced as jump-conditioned
  numerical fluxes at the particle interface (arXiv:2211.14390 pattern);
  `SpinSourcedProblem`/`solve_spin_sourced`/`spin_sourced_flux`. **Decisive gate — TD vs the
  in-package FD oracle** (`TeukolskyPointParticleModeCircular`), s = −2, r₀ = 10M: (2,2)
  ``\dot E^\infty`` 1.3e-3 / ``\dot E^H`` 5.4e-5; (3,3) 2.8e-3 / 6.2e-6; (2,1) 1.5e-2 /
  4.4e-6; s = +2 (2,2) 3.8e-4 — with resolution-ladder convergence asserted.

**`[FUTURE]` boundaries (documented in-source, constructor-enforced where applicable):**
Kerr ``a \ne 0`` *sourced* TD (the spheroidal→spherical source re-projection is not in
print — dossier `DG_SPIN_FORMULATION.md` §8.3; the sourced constructor enforces a = 0);
``|s| = 1`` sourced; m = 0 sourced excluded; the s = +2 tail exponent below the numerical
floor. Vacuum Kerr is *not* restricted.

## Dependencies and the poison-axis budget

Runtime `[deps]` (from `Project.toml`):

- **GrTool:** `NumericsCore` (continued-fraction, interpolation, conventions, precision),
  `BigSpecial` (``\Gamma`` / ``{}_2F_1`` / Pochhammer / Kummer), `Spacetime` (geodesic constants,
  frequencies, roots; `SpinningSecondary` equatorial spin-frequency corrections), `Harmonics`
  (SWSH eigenvalue + ``S``).
- **stdlib:** `LinearAlgebra`, `Printf`.
- **external:** `HypergeometricFunctions`, `SpecialFunctions` — the *sanctioned* light closure
  (the same one legacy `Teukolsky` carried), reaching the user through `BigSpecial`.

**Poison-axis declaration** — this package is the light-dep firewall of GrTool. The binding CI
budget lives in `test/hygiene.jl` (design §7.7 gates 1+2), which asserts TOML-statically that
the runtime `[deps]` are exactly the allow-list above and **forbids**

```
HDF5, JLD2, Symbolics, OrdinaryDiffEq*, CUDA, AMDGPU, KernelAbstractions,
Kernels, Data, WignerSymbols, Roots, NLsolve, SciMLBase
```

and additionally asserts *after* `using Perturbation` that none of these poison axes is loaded
into the session (closure-load gate). "External-ODE-free" rests on the vendored DP5 +
Sasaki–Nakamura radial integrators — these must never be refactored onto OrdinaryDiffEq.
`TestKit` and `Data` appear only in `[extras]`/`[targets].test` (test harness), never in the
runtime closure.

## GPU / data needs

None / none in the solver path — this is the HDF5-free firewall; all outputs are returned in
memory and persisted downstream (via `Data`) by the consumers. The BigFloat MST/continued-fraction
generation tier stays on CPU (BigFloat-CPU / Float64-GPU boundary); heavy outer flux sweeps scale
via `Distributed` pmap + CPU threads ≤ 20 under the BLAS-guard. Any future GPU flux-grid
generation must remain an optional extension, never a hard `Kernels` dep, or the firewall is lost.

## Legacy provenance (drawsFrom)

Functions relocated to this physics home from (design §6 crosswalk rows dest=Perturbation):

| Legacy package | What moved in |
| --- | --- |
| `Teukolsky` | homogeneous radial In/Up (`TeukolskyRadial`, `TeukolskyRadialFunction`, `teukolsky_radial_deriv`), `RenormalizedAngularMomentum` ``\nu``, MST coefficients, amplitudes ``Z/C``, point-particle modes + fluxes |
| `ScatteringPhaseMST` | `starobinsky_constant`, ``B^{\mathrm{refl}}/B^{\mathrm{inc}}`` internals, `scattering_phase` / `scattering_phase_numerical` / `scattering_phase_table`, `mst_coefficient`, `mst_convergence` |
| `TeukolskySpinFluxes` | all **12** exported `TeukolskySpinMode*` spinning-secondary flux entrypoints (S5-2 port), characterized vs baked BHPToolkit-MMA fixtures |

Grounding: test-only characterization against the (light) legacy oracles `Teukolsky` and
`ScatteringPhaseMST` at their frozen tolerances, plus analytic identities (``\nu``
continued-fraction root, MST ``{}_2F_1`` connection, flux Newtonian-quadrupole ``32/5`` limit,
Starobinsky–Chandrasekhar scattering identity). Since Phase 7-B the oracle outputs are **baked
fixtures** (`test/fixtures/leg_*.jl`) — no `packages/*` edge exists at test time.

Two later absorptions came from *outside* the `packages/*` legacy set, grounded the same way
(baked one-shot fixtures + in-package referees): the DG-Teulolsky research repo →
`Perturbation.TimeDomain` (fixtures `test/fixtures/td_*.jl`), and the vendored
`GeneralizedSasakiNakamura.jl` v0.8.5 → the `:GSN` backend (fixtures
`test/fixtures/gsn_oracle.jl`, 706 rows).

## The SpinFluxes surface (12 entrypoints)

The complete spinning-secondary linear-in-spin ``s=-2`` surface, with the BHPToolkit names
preserved:

- `TeukolskySpinModeCircularCorrectionAnalytical` (circular-equatorial) and
  `TeukolskySpinModeEquatorialCorrectionAnalyticalPhases` (eccentric-equatorial) — the original
  circular/eccentric slice;
- the ten S5-2 ports: `TeukolskySpinMode`, `TeukolskySpinModeFromCorrection`,
  `TeukolskySpinModeAnalytical`, `TeukolskySpinModeCorrection`,
  `TeukolskySpinModeCorrectionAnalytical`, `TeukolskySpinModeSpherical`,
  `TeukolskySpinModeSphericalCorrectionNum`, `TeukolskySpinModeSphericalCorrection`,
  `TeukolskySpinModeSphericalCorrectionAnalytical`,
  `TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin`.

Each returns a NamedTuple `(; l, m, n, k, omega, α, Amplitudes = (I, H), S, Fluxes =
Energy/AngularMomentum → (I, H), …)`. The convenience methods materialize the spinning orbit
internally (`KerrSpinOrbit*` kinematics, local to this package); the core methods accept a
pre-built orbit NamedTuple so a whole ``(l,m,n,k)`` sweep reuses one orbit solve.

## What changed vs legacy

- **Entrypoint count corrected: 12, not 13.** Design §4.3 listed "13 `TeukolskySpinMode*`
  entrypoints"; the legacy `TeukolskySpinFluxes` exports exactly 12 (its own runtests asserts
  "12 public API symbols") — the 13 counted a wrapper list with one duplicate. This package
  exports the corrected 12.
- **Static path is now type-generic `{T}` (S3 fix).** The legacy static ``\omega = 0`` branch
  was precision-capped; here the working precision `T = float(a)` flows through the closed-form
  ``{}_2F_1`` internals (`Complex{T}` via `BigSpecial.gamma_big`), so `:Static` inherits BigFloat
  precision like every other backend. Characterized against the baked pre-`{T}` references.
- **`scattering_phase_analytic` is implemented** (design §4.3 flagged it as an unported stub):
  the published `phase{slm}Expansion` low-frequency ``\varepsilon``-series (exact through
  ``\varepsilon^3``) are transcribed as explicit polynomials in
  ``\gamma_E, \pi, \zeta(3), \log`` for the tabulated modes ``(s,l,m) \in
  \{(0,0,0), (0,1,0), (-1,1,0), (-2,2,0)\}`` — with **no** symbolic-series dependency (firewall
  preserved).
- **Local spinning-orbit kinematics.** `Spacetime.SpinningSecondary` carries only the
  equatorial slice (`spin_linear_frequencies`, `spin_mino_time_correction`); the
  inclined / numeric-correction / spherical spinning-orbit kinematics live *here* as unexported
  internals (`KerrSpinOrbit`, `KerrSpinOrbitCorrection`, `KerrSpinOrbitNum`, …) — zero new deps.
- **Namespace collapse.** The legacy three-package split (`Teukolsky`, `ScatteringPhaseMST`,
  `TeukolskySpinFluxes`) becomes one flat `Perturbation` namespace; the Phase-0
  module-name collision around the root spinning-flux packages is resolved by this fold. Function
  names are otherwise preserved 1:1.
- **Test oracles baked (Phase 7-B).** The former live test-oracle path-deps (`Teukolsky`,
  `ScatteringPhaseMST`, `KerrGeodesics`, `SpinWeightedSpheroidalHarmonics`) were removed from
  `[sources]`; their outputs on the exact test grids are baked into `test/fixtures/leg_*.jl` at
  the same tolerances.
- **The `radial_sn` boundary is closed by `:GSN`.** The `:SasakiNakamura` backend's documented
  upstream limit ("`s = -2`, Up only"; everything else `NotYetPorted`) is superseded: `:GSN`
  serves the full ``s \in \{0, \pm 1, \pm 2\}`` × {In, Up} surface at real ω. `radial_sn` stays
  as an independent cross-check backend.
- **A time-domain solver now lives beside the frequency-domain engine.** The
  `Perturbation.TimeDomain` nested module (the one deliberate exception to the flat-namespace
  rule) brings DG TD evolution in with zero new dependencies — the firewall budget is unchanged.

## Public API

```@autodocs
Modules = [Perturbation, Perturbation.TimeDomain]
Order = [:module, :type, :constant, :function]
```

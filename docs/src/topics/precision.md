# [Precision / BigFloat inheritance](@id topic-precision)

How arbitrary precision flows through GrTools: one policy home in
`NumericsCore`, two BigFloat-generic spine roots, a type-generic `{T}` radial chain
from [Perturbation](../capabilities/perturbation.md) through
[MetricReconstruction](../capabilities/metricreconstruction.md) to
[SelfForce](../capabilities/selfforce.md)/[AccretionDisc](../capabilities/accretiondisc.md),
one *named* Float64 ceiling on the angular side, and a Float64/GPU tier on top that
imposes no constraint downward. This page is the **single home** of the topic
— package pages link here rather than restating it.

The governing rule everywhere below: numerical code is type-generic `{T}`, so a
`BigFloat` input under the active `setprecision(BigFloat, bits)` (256/512/1024/2048 in
practice) inherits its full mantissa, while `Float64` input reproduces the validated
legacy path bit-for-bit. A hard `Float64` literal in a tolerance or an accumulator
would silently cap the result at ~16 digits regardless of mantissa — the bug class the
`{T}` discipline exists to prevent.

## [The policy home: `NumericsCore` precision selectors](@id precision-policy)

`NumericsCore/src/precision.jl` canonicalizes the BigFloat plumbing that was scattered
across the legacy SWSH/Teukolsky/H1Lorenz/KerrLorenzCirc codes into precision-inheriting
selectors. Every selector threads the element type `T`:

| Selector | `Float64` | `BigFloat` (`bits = precision(T)`) | Role |
|---|---|---|---|
| `bc_tol(T)` | `1e-14` (validated C value) | `max(2.0^(10-bits), 8*eps(T))` | series / boundary-condition truncation |
| `rk_tol(T)` | `1e-12` (validated C value) | `2.0^(20-bits)` | adaptive-RK step/solve tolerance |
| `fd_step(T; kind)` | `T(1)/T(10)^6` radial, `T(1)/T(10)^4` theta | same expressions (exact in `T`) | finite-difference step for background derivatives |
| `bigprec_for_mag(magpow)` | — | `clamp(128 + 6*ceil(magpow), 160, 320)` bits | mantissa sized to a working magnitude `~10^magpow` |
| `with_precision(f, bits)` | — | `setprecision(f, BigFloat, bits)` do-block (restored on exception) | scoped global-precision runs |

The vendored ODE family (`NumericsCore/src/ode.jl`) consumes `rk_tol(T)` directly: DP5
is type-generic and BigFloat-clean, GBS is the high-order BigFloat path
(`gbs_kmax(T)` returns 0 for `Float64`), DOPR853 is the Float64 FEW trajectory
integrator. The same selector formulas appear as private per-package selectors where a
capability owns its solve loop: `SelfForce.LorenzField` carries `_bc_tol`
(`Boundary.jl`) and `_rk_tol` (`Solver.jl`), and the MetricReconstruction Kerr radial
carries `_bigprec_for_mag` (`KerrRadial.jl`) — identical formulas, kept local to the
consuming solver.

## [The two BigFloat spine roots](@id precision-spine-roots)

Everything high-precision is seeded by two type-generic foundations:

- **Harmonics eigenvalues.** `SpinWeightedSpheroidalEigenvalue(s, l, m, γ)` with
  `method=:Auto` routes a `BigFloat` γ to the type-generic Leaver continued-fraction
  backend — which inherits the active `setprecision` — and a `Float64` γ to the fast
  LAPACK spectral backend (`SpheroidalEigenvalues.jl`). The asymptotic-series backends
  (`:SmallGamma`/`:LargeGamma`) have Float64 coefficient recurrences and therefore carry
  a Float64 accuracy floor even for BigFloat γ — the BigFloat home is the Leaver path.
  `return_residual=true` routes Leaver unconditionally and returns `(λ, |Eq(A)|)`, the
  eigenvalue plus converged CF residual consumed by the BigFloat radial
  precision-convergence gate (the fold of KerrLorenzCirc `swsh_eigenvalue_big`).
- **Spacetime.** Geodesic constants and Fujita–Hikida frequencies
  (`Geodesics/constants.jl`, `kerr_frequencies`) and the spinning-secondary
  closed-form/algebraic corrections (`SpinningSecondary`) are pure arithmetic on
  NumericsCore's complex Carlson elliptic layer — "type-generic `{T}` throughout so
  BigFloat precision is inherited, never hard-Float64'd" (`Spacetime.jl`).

## [The angular asymmetry: `angular_precision_floor()`](@id precision-angular-floor)

The angular solve is deliberately asymmetric, and
`Harmonics.angular_precision_floor([T=BigFloat])` (`Harmonics/src/precision.jl`) is the
public contract naming it:

```julia
angular_precision_floor(T) == (value_tier = eps(Float64), eigenvalue_tier = eps(T))
```

- **eigenvalue tier — floor `eps(T)`.** The Leaver CF eigenvalue path is type-generic
  and tracks the active `setprecision(BigFloat, bits)`.
- **value tier — ceiling `eps(Float64)`.** The harmonic *value* surface
  (`SWSHSFunction` Float64 expansion coefficients, the LAPACK spectral eigensolve, the
  spin-weighted-``Y`` downcast) is Float64-capped: ``S(\theta)`` and
  ``\partial_\theta S(\theta)`` stay at `eps(Float64)` even under BigFloat input.

The cap is accepted by design: precision is driven through the **radial** chain,
not the angular value — KerrLorenzCirc-derived code Float64-casts its angular factors
and relies on `eigenvalue_tier ≤ value_tier`. The function also marks the seam where a
future un-downcast BigFloat angular-value path would lift the floor.

## [The radial inheritance chain](@id precision-radial-chain)

### Perturbation — MST / CF / renormalized ν / static, all `{T}`

- **ν (renormalized angular momentum).** `RenormalizedAngularMomentum` with
  `method=:Auto` sends BigFloat inputs to the precision-inherited `:Monodromy`
  confluent-Heun path (all recurrences and Pochhammer products in `Complex{BigFloat}`,
  run under `setprecision(BigFloat, bits)` with `bits = max(256, precision(BigFloat))`)
  and Float64 inputs to machine-precision `:FindRoot` (`renorm_angmom.jl`).
- **MST coefficients.** `_mst_fn_array` infers its accumulator type from the
  recurrence-ratio type — `ComplexF64` for Float64 inputs, `Complex{BigFloat}` for the
  high-``\ell`` small-``\omega`` extended-precision amplitude path (`mst.jl`); the CF
  tails reuse `NumericsCore.cf_convergent`/`cf_finite`.
- **Static (ω = 0) branch — de-capped to `{T}`.** The closed-form ``{}_2F_1`` path
  (`teukolsky_radial_static.jl`) was ComplexF64-pinned and is now type-generic: the
  working precision `T = float(a)` threads ``\tau``/``\kappa``, the asymptotic
  amplitudes, the norm and the radial value/derivative, so a BigFloat `a` gives
  `Complex{BigFloat}` end-to-end. Γ dispatches per precision —
  `_static_gamma(z::ComplexF64) = SpecialFunctions.gamma(z)` keeps Float64
  **bit-identical** to the pre-`{T}` baked regression literals, every other `Complex{T}`
  routes `BigSpecial.gamma_big`; ``{}_2F_1`` goes through `_safe_2f1` (byte-identical
  upstream `_₂F₁` for Float64, BigSpecial Pfaff/``z\to\infty`` fallback for BigFloat).
- **`:GSN` backend — `{T}`-generic with a documented BigFloat truncation.** The
  Generalized Sasaki–Nakamura backend (`gsn_*.jl`) is type-generic
  end-to-end: the closed-form layers (potentials, conversions, boundary series) are exact
  in `T` (256-vs-512-bit agreement < 1e-60), and the radial solve runs the *same*
  type-generic vendored DP5 mesh at every precision, giving cross-precision solve-tier
  determinism at the 1e-74 class. The **deviation**: for BigFloat input the solve keeps a
  default **1e-12 step truncation** instead of the `2^(20-bits)` `rk_tol(T)` deep
  tolerance — the GBS phase-form (Riccati) tier planned for BigFloat has a *measured*
  non-``h^2`` extrapolation floor near the horizon (at tol = 1e-15 the ladder bottoms out
  and rises again after ~214k RHS calls; at tol = 1e-10 the 256- and 512-bit amplitudes
  agree with each other but sit 5e-6 from Float64 — the floor biases accuracy, not just
  cost; evidence in `gsn_solve.jl`). So BigFloat GSN runs are exact in *arithmetic*
  (precision threading, determinism) under an algorithmic truncation shared across
  precisions; **deep-tolerance BigFloat GSN integration is a documented `[FUTURE]`**.
- The public `TeukolskyRadial(s, l, m, a, ω)` dispatcher composes all of the above, so
  its In/Up solutions are BigFloat-inheritable for both ``\omega \ne 0`` (MST/Monodromy)
  and ``\omega = 0`` (static).

### MetricReconstruction — Kerr and Schwarzschild radial reconstruction

- **`SchwarzschildRadial` (a = 0) — de-capped to `{T}`.** `jump_data_radial(r0, l, m)`
  takes its element type from `r0` via `orbit_background` (`T = typeof(bg.r0)`;
  BigFloat in ⇒ BigFloat out); the upstream chain
  (`orbit_background`/`jumps_pm`/`hj1`/`closure_a0`) is already `{T}`, the reconstruction
  structs are `HomField{T}`/`KappaField{T}`, the variation-of-parameters quadratures are
  precision-aware via `atol`/`rtol`, and the ``\kappa`` chain self-converges 256→512
  bits. Float64 input reproduces the pre-`{T}` path exactly (the conversions are exact).
- **`KerrRadial` (a ≠ 0) — magnitude-triggered BigFloat.** The homogeneous
  ``R^{\rm up}/R^{\rm in}`` legs are Float64 until the sampled `|R^up|` exceeds
  `_BIGMAG = 1e8` (where Float64 loses ≳8 relative digits); the leg is then rebuilt as a
  **coarse-node BigFloat Chebyshev evaluation + Float64 barycentric interpolation** —
  ~26 BigFloat node evals per leg inside `setprecision(BigFloat, bits)` blocks with
  `bits = _bigprec_for_mag(magpow)`, giving full relative precision at near-Float64
  evaluation cost. The dynamic-range problem lives only in *building* the field, not in
  using it.

### Downstream consumers

- **SelfForce.** The owned LorenzField solver/boundary machinery is type-generic `{T}`
  with the precision-inherited `_bc_tol`/`_rk_tol` selectors above; parallel BigFloat
  mode builds use a per-task `setprecision(BigFloat, bits)` for thread safety
  (`Solver.jl`). Only the second-order-Ricci GPU path is Float64 (next section).
- **AccretionDisc.** Every struct is type-generic (`RelativisticDisc{T<:Real}`,
  `FluidVars{T<:Real}`), so disc modes and the disc force evaluate at whatever precision
  the reconstruction chain delivers.

### Neutral and partially capped

- **PostNewtonian is precision-neutral.** It emits `Symbolics.Num` coefficient
  expressions; downstream substitution evaluates at the caller's chosen precision — it
  neither caps nor forces anything.
- **ModifiedGravity is partially de-capped.** `hertz_potential.jl` threads `T` (no
  `float(a)`/`float(ω)` coercion, no ComplexF64 pin), but the source/RHS amplitude
  accumulators (`source_amplitudes.jl`, `rhs_amplitudes.jl`) still initialize as
  `zero(ComplexF64)` — a remaining Float64 cap on that surface.

## [The Float64/GPU tier](@id precision-float64-tier)

The 256..2048-bit generation tier (Harmonics eigenvalues, Perturbation MST/CF,
MetricReconstruction Kerr radial) **never touches the device**. On top of it:

- **Perturbation.TimeDomain is a Float64-validated tier.** The DG TD problem structs are
  `{T}`-parametrized (`Mode1DProblem{T}`, `SpinModesProblem{T}`, …, default `T = Float64`)
  and the underlying `NumericsCore.dg` spectral operators are BigFloat-clean
  (`RefElement(N; T=BigFloat)`), but every TD characterization/flux gate runs at `Float64`
  — TD evolution accuracy is set by DG resolution ``(K, N)`` and the RK4 step, not by the
  mantissa, and no BigFloat TD run is claimed or gated. Treat TD outputs as
  Float64-class; the BigFloat generation tier remains the frequency-domain radial chain
  above.
- **EMRIWaveforms** is Float64/GPU-only: it calls the BigFloat-generic lower tiers *at*
  `Float64` and imposes no precision constraint upward (`EMRIWaveforms.jl`).
- **Hybrid** is a Float64 tier by construction: its model containers are
  `Interp1{Float64, ComplexF64}`/`Interp1{Float64, Float64}` interpolants over Float64
  grids (`Models/`).
- **Kernels**: "Float64 ON DEVICE ONLY" (`Kernels.jl`, binding). The CPU time-domain
  mode-sum kernel itself stays type-generic (BigFloat phases/amplitudes flow through),
  but the FD-SPA path is Float64 by construction (Float64 series constants), and the
  `r2_1..r2_10` second-order-Ricci kernels are pure Float64 array arithmetic. CPU-vs-GPU
  `sincos`/`fma` divergence is absorbed by TestKit's `:gpu` tolerance class
  (`rtol = 1e-4`) — bit-equality is required only CPU-serial vs CPU-threaded (see
  [GPU / Parallel](gpu-parallel.md)).

## [BigFloat and threads](@id precision-threads)

Two operational rules interact with precision (details in
[GPU / Parallel](gpu-parallel.md)): BigFloat workloads scale via processes / modest
thread counts, not wide threading (GC contention makes `-t16` net slower than `-t4`,
measured), and every threaded LAPACK/BLAS region — including the Float64 spectral
eigensolve behind the angular value tier — runs under the
`NumericsCore.blas_single_thread` guard. Parallel BigFloat mode builds set precision
per task (`setprecision` is task-local state), as in `SelfForce.LorenzField.Solver`.

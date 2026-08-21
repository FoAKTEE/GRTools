# `AccretionDisc/derivation` — route R3, the Julia-native derivation

Standalone project that derives the relativistic master system (47)–(49), the eliminated velocities
and the diagnostics integrands of **arXiv:2608.17003** from the paper's own equations, entirely
inside Julia with `Symbolics.jl`.

It is **not** a dependency of `AccretionDisc`: nothing under `capabilities/AccretionDisc/src` may
`import R3Derivation`, and the package dependency budget of `DESIGN.md` §4.10 is unchanged. The
pattern is the one already used by `test/SpiralWaves/fixtures/capture/` and `PostNewtonian/test/` — a
quarantined environment with a committed `Manifest.toml` that resolves **offline** (`Symbolics`
7.36.0, the version the depot already carries).

```
julia --project=capabilities/AccretionDisc/derivation -t4 \
      capabilities/AccretionDisc/derivation/scripts/derive_all.jl 50   # full run + artefacts
julia --project=capabilities/AccretionDisc/derivation -t4 \
      capabilities/AccretionDisc/derivation/test/runtests.jl           # fast suite (R3_NPTS=12)
```

## Why route R3

Today `SpiralWaves/master_generated.jl` is transcribed from the authors' supplementary *Mathematica*
notebook (route **R1**) and `r2_generated.jl` from an independent *Mathematica* derivation (route
**R2**). Route R3 removes every external CAS from the build path. R1 and R2 are **demoted to
oracles**: they gate R3 at `EGD::jd-equivalence-gate`, they never feed the build. Nothing in this
project reads an R1/R2 expression — the independence is structural, not a promise.

## Layout (stages 1-2)

| file | content |
|---|---|
| `src/symbols.jl` | the alphabet and the radial derivation `ddr` (chain rule over `Symbolics.derivative`) |
| `src/geometry.jl` | `christoffel_from_metric(ginv, dg)` — the ONE generic connection helper (promotion candidate for `Spacetime`) |
| `src/jets.jl` | first-order λ-jets carrying the Fourier phase implicitly |
| `src/fluid.jl` | Eq. (3) written **once**, element-type generic — used for the background and for its linearisation |
| `src/background.jl` | line element (4), u^α (5),(8)–(10), hydrostatic balance (11) |
| `src/certify.jl` | the admission modalities: symbolic, exact-rational, pole order |
| `src/series.jl` | truncated Laurent series over ℚ + a direct evaluator |
| `src/stage1.jl` | background certificate suite |
| `src/perturbation.jl` | the linearised system (15b),(15c) and the perturbed normalisation |
| `src/reduce.jl` | the reduction M:604–608 → C, S, δû^t, δû^φ, δε̂ |
| `src/newtonian.jl` | the paper-facing gate: Newtonian limit ≡ Eq. (50) |
| `src/stage2.jl` | back-substitution and shape gates for the master system |
| `src/serialize.jl` | the hand-off artefacts consumed by `EGD::jd-codegen` |

Stage 3 (diagnostics) lands in the next node of the series.

## Design decisions that make it work

* **Generic background.** ν, ψ, μ, ω̃ are free functions of r; the Kerr closed forms (A1a–e) are never
  used here — they enter only when the emitted code is evaluated, through
  `Spacetime.kerr_equatorial_potentials`.
* **No transcendentals.** e^ν, e^ψ, e^μ are the symbols `enu, epsi, emu` and W = 1/√(1−v²) is a symbol
  with the single relation W²(1−v²) = 1. Every expression is then a rational function of the
  alphabet — which is what makes the certificates exact and the emitted code precision-inheriting.
* **Nothing is hand-derived.** The perturbed Christoffels, δq^{μν} and δ(u^ν∇_νu^μ) come out of jet
  arithmetic; each elimination is a `solve_linear`/`lin_subs` with an asserted affineness.
* **Two exact modalities.** `simplify` where it closes; otherwise evaluation at random exact
  `Rational{BigInt}` points with residual **exactly** `0//1`. Rational arithmetic has no rounding, so
  the second is an exact modality, not a numerical one. The conic parametrisation
  v = 2t/(1+t²), W = (1+t²)/(1−t²) removes the one algebraic constraint and makes `simplify_fractions`
  a decision procedure — substituting W → 1/√(1−v²) instead leaves uncollapsed `(1−v²)^{n/2}` powers
  and stalls (recorded as a finding).

## Findings recorded while building it (see the node reports)

1. `Symbolics.substitute` refuses to bind a **complex** value to a real-typed symbol
   (`promote_symtype(complex, …)` asserts `S <: Real`); every amplitude replacement therefore goes
   through `lin_subs`, coefficient extraction, which also asserts linearity.
2. `Symbolics.simplify_fractions` alone is **not** a normal form — it leaves a sum of fractions;
   `expand` first, then `simplify_fractions`, gives the single fraction the pole-order test needs.
3. Substituting W → 1/√(1−v²) leaves uncollapsed `(1−v²)^{n/2}` powers and stalls; the conic
   parametrisation closes the same identity instantly.
4. Taking the Newtonian limit by symbolic substitution + `expand` blows up (>10 GB resident, no result
   in 12 min); evaluating the same tree in ℚ((s)) finishes in seconds.

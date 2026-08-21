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
      capabilities/AccretionDisc/derivation/scripts/emit_sources.jl    # stage 4: emit the shipped sources
julia --project=capabilities/AccretionDisc/derivation -t4 \
      capabilities/AccretionDisc/derivation/test/runtests.jl           # fast suite (R3_NPTS=12)
```

## Why route R3

The master system used to be transcribed from the authors' supplementary *Mathematica* notebook
(route **R1**) and its diagnostics derived in an independent *Mathematica* derivation of our own
(route **R2**). Route R3 removes every external CAS from the build path; R2 is **demoted to an
oracle** that gates R3 at `EGD::jd-equivalence-gate` and never feeds the build
(`test/SpiralWaves/fixtures/oracle/{r2_master_generated,r2_generated}.jl`), and the notebook-derived
R1 oracle was **removed from this repository entirely** on 2026-08-22 by owner decision
(`EGD::jd-notebook-purge`) — its lineage and the full R1-vs-R2-vs-R3 equivalence evidence are
archived in the mission repository (emri-disk `results/emri-gasdisk` `paper_EMRI-GASDISK`,
`evidence/{R1-vs-R2-gate,jd-equivalence-gate,jd-notebook-purge}`). Nothing in this project ever read
an R1/R2 expression — the independence is structural, not a promise.

## Layout

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
| `src/diagnostics.jl` | ξ from (16), U/V/W (32),(33),(42), δ_ξT (40),(41), (72a),(72b), the (69) gate |
| `src/serialize.jl` | the hand-off artefacts of stages 2-3 (`out/R3_*.jls.txt`) |
| `src/codegen.jl` | stage 4: the Julia printer — artefact → `src/SpiralWaves/{master,r2}_derived.jl` |

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
5. `Symbolics.build_function` was evaluated for stage 4 and is NOT used: it emits a correct function
   but with generated argument names and no control over the wrapper (signature, `MasterInputs`
   destructuring, docstring, provenance header), while the emission conventions the shipped files
   need (`n//d` rationals, hoisted `enu^a*epsi^b*emu^c`, no transcendental) are properties of the R3
   ALPHABET, not of the printer.  `string(::Num)` plus a 40-line wrapper emitter does the whole job,
   and the emitted text is certified by an EXACT round-trip (`codegen_selfcheck`): the file is
   re-parsed with `Meta.parse`, interpreted over `Rational{BigInt}` and compared with `exact_eval` of
   the R3 expression — residual exactly `0//1` at every point.

## Stage 4 — what the emitter ships

`master_derived.jl` defines `_master_entries(inp::MasterInputs)` and `r2_derived.jl` defines
`_r2_delta_u_phi`, `_r2_amf`, `_r2_torque_density` — the SAME names and signatures as the route-R1/R2
generated files they replaced, so the swap was a one-line include change at every call site.  Two
consumer conventions are applied at emission time and are documented in each file header:
`circmetriczero` (ĥ_tr = ĥ_rφ = 0 for a circular equatorial orbit — which is what makes every
imaginary part of the master system vanish identically), and the change of variables
ω̃ = Ω − v e^{ν−ψ}, v′ = (Ω′−ω̃′)e^{ψ−ν} − v(ν′−ψ′) of (10), because `MasterInputs` carries Ω and Ω′
where route R3 carries ω̃ and v′.  The emission is deterministic (two runs byte-identical, asserted in
the test suite) and the shipped files are asserted to be exactly what the emitter produces.

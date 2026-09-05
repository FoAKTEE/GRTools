# [Reproducing arXiv:2606.04998](@id topic-reproducing-2606-04998)

Vu, Nishimura, Osburn, Thompson, Kidder, Upton & Wardell, *Self-force calculations with
numerical relativity methods* — the m-mode **elliptic** Kerr scalar self-force. This page is
the index: what is reproduced, where it lives, how to re-check it, and what is still open.

The physics is in two packages and nowhere else. [EllipticCore](../core/ellipticcore.md) is the
physics-free complex-parametric 2D hp-DG elliptic solver stack (reference elements, multi-block
meshes and mortars, the matrix-free operator, boundary-data hooks, additive Schwarz, geometric
multigrid, complex GMRES, hp-AMR); [ScalarSelfForce](../capabilities/scalarselfforce.md) owns
the m-mode formulation (orbits, the first-order flux form, the vtu null slicing, the
Detweiler–Whiting puncture and effective source, the block lattice, the worldtube, the jump and
boundary tables, the solve drivers, the extraction and the mode sum). The drivers under
`reproduce/` are **thin**: zero physics, per the `reproduce/README.md` contract.

## What is reproduced

Six of the paper's seven figures have landed drivers and goldens. Everything below is
registered in `reproduce/MANIFEST.toml`, with the per-figure evidence, gates, VLM verdicts and
holes written up at length in `reproduce/arXiv-2606.04998/README.md` — that file is the
normative record, this page is the map to it.

| paper figure | id | driver | golden | tier |
|---|---|---|---|---|
| Fig. 1 `fig:domain` — the DG domain, the worldtube, the AMR cascade | `fig_domain` | `reproduce/arXiv-2606.04998/fig_domain.jl` | `expected/fig_domain.json` | heavy |
| Fig. 2 `fig:matrix` — DG operator sparsity at the coarsest grid | `fig_matrix` | `reproduce/arXiv-2606.04998/fig_matrix.jl` | `expected/fig_matrix.json` | heavy |
| Fig. 4 `fig:convergence` — the hp-AMR convergence of every ``F_r^m`` | `fig_convergence` | `reproduce/arXiv-2606.04998/fig_convergence.jl` | `expected/fig_convergence.json` | heavy |
| Fig. 5 `fig:modes` — 21 stacked ``\log_{10}\|\Psi_m\|`` planes at three spins | `fig_modes` | `reproduce/arXiv-2606.04998/fig_modes.jl` | `expected/fig_modes.json` | heavy |
| Fig. 6 `fig:linsolv` — the GMRES residual staircase over AMR grids | `fig_linsolv` | `reproduce/arXiv-2606.04998/fig_linsolv.jl` | `expected/fig_linsolv.json` | heavy |
| Fig. 7 `fig:ssf_error_budget` — the ``F_r`` error budget over spin, at ``r_0 = r_\mathrm{ISCO}`` and ``r_0 = 10M`` | `fig_error_budget` | `reproduce/arXiv-2606.04998/fig_error_budget.jl` | `expected/fig_error_budget.json` | heavy |

Beyond the figures, the mode sum itself is gated against an **external** table nothing in this
suite feeds: the summed ``(M^2/q^2)F_r`` at ``r_0 = 10M`` against the printed Warburton–Barack
Table III at ``a/M = 0, \pm0.5, \pm0.9`` (worst relative deviation ``6.4\times10^{-6}``), with
``\sum_m F_t^m`` against their Table II as a second, regularization-free gate. Those numbers,
their error budget and their provenance live on the
[ScalarSelfForce page](../capabilities/scalarselfforce.md).

Fig. 3 `fig:precond`, the preconditioner schematic, is not claimed at all; nothing on this page
should be read as covering it.

### Fig. 7 is `preliminary`, and its failures are pinned rather than scoped away

The error-budget figure is worth reading with its `[[result]]` row open, because it is the one
row here whose registered verdict is a FAILURE. Its three curves per column are main.tex:722's:
the solid reference error against Warburton–Barack, the dashed finite-``m`` truncation estimator
and the dotted numerical-resolution estimator. Reading the published figure's own vector paths
established two facts the raster cannot give: it samples **nine** spins, not the five it labels,
and its solid curve carries only the seven with ``|a| \le 0.9`` — Warburton–Barack's tables stop
there, at every radius.

What reproduces: the published truncation curve to 0.3–3.6 % at 17 of 18 configurations, the
reference curve to within 23 % at all 8 printed comparison points, ``|F_r|`` to ``\le 0.19`` %
across the strip, and ``F_t`` against the printed nine-figure Table II to
``2.91\times10^{-8}``.

What does not, stated as the row states it — **every gate runs on all 18 configurations, so
`gates_pass` is `false` and `gate_failure_ids` names the exact five**: GATE 3 (relative
truncation within a factor 2 of the published curve) fails at ``a = +0.998`` at its ISCO, at
121×, and GATE 8 (the resolution curve bounded a decade from above) fails at four
configurations. GATES 2, 2H, 4, 5, 6 and 7 pass everywhere. `--check` verifies
REPRODUCIBILITY, not the verdict: the verdict is a registered value pinned by the golden, so a
gate that starts failing — or a known failure that quietly disappears — is a golden MISMATCH.

One interpretive choice governs the plotted solid curve and is flagged as such: the published
``F_r`` appears to carry **no fitted tail** (its reference/truncation ratio, 2.81–4.63, is the
``m_\mathrm{max}/(p-1)`` factor between the true remainder and the omit-the-largest-``m``
estimator). The paper does not say so; the driver plots that reading, registers the
Heffernan-tail ``F_r`` beside it from the same solved modes, and says that if the reading is
wrong the plotted solid curve is the wrong estimator.

## Re-checking it

Every row is deterministic and re-checked from its golden by `reproduce/run_all.jl`:

```bash
JULIA_PKG_OFFLINE=true julia --project=reproduce -e 'using Pkg; Pkg.instantiate()'
julia --project=reproduce reproduce/run_all.jl --check --paper 2606.04998
julia --project=reproduce reproduce/run_all.jl --check --only fig_linsolv   # one row
julia --project=reproduce reproduce/run_all.jl --check --only fig_error_budget
julia --project=reproduce reproduce/run_all.jl --regen --only fig_linsolv   # rewrite a golden
```

All six rows are **heavy tier**, so they are outside the fast set CI runs — re-check them
manually. Determinism rules are the suite's own (see [Reproducibility](reproducibility.md)):
`OPENBLAS_NUM_THREADS=1`, results bit-identical at `-t1` and `-t4`, the pinned
`reproduce/Manifest.toml` and Julia 1.12.6. A golden is rewritten only with `--regen`, after
inspecting what moved, and the commit body says why.

The two package suites are the gate underneath all of it:

```bash
JULIA_NUM_THREADS=4 julia --project=core/EllipticCore          -e 'using Pkg; Pkg.test()'
JULIA_NUM_THREADS=4 julia --project=capabilities/ScalarSelfForce -e 'using Pkg; Pkg.test()'
```

## Two lattices, and an erratum

The paper ships its own input file (`anc/ScalarSelfForce.yaml`), and the reproduction follows it
by default — including `ImposeEquatorialSymmetry: True`, the ``z^2`` half domain with
20 elements / 320 grid points at the initial lattice. **The published figures were not made on
that lattice.** Vector extraction of two published figures says so directly: `fig:convergence`'s
abscissa starts at ``N = 480`` for all 21 curves, and `fig:linsolv`'s level-0 residual, 1.394,
is the full-``z`` ``\|b\|_0`` (ours 1.391, 0.18 % out) rather than the half domain's 1.744. So
the figures were run with the symmetry OFF, on the 5×3 full-``z`` lattice — 30 elements /
480 points. That is recorded as an erratum against the shipped input file; the package default
stays `true` (the paper's own file), and `crp_paper_config(; equatorial_symmetry = false)`
builds the published-figure lattice. `fig_linsolv` and `fig_convergence` are frozen on it, with
the production half domain registered beside them as a secondary series so the discrepancy is
data rather than prose.

A second erratum is a print typo with no downstream consequence: the last term of the radial
effective potential (main.tex:190) carries the wrong sign — `dV_eff/dr` of the printed potential
does not vanish on the paper's own circular-orbit ``(E, L)`` (normalized residual ``\sim 0.5``
at all 35 fixture points), while the corrected circularity condition holds to
``2.1\times10^{-15}`` at Float64 and ``8.7\times10^{-78}`` at BigFloat-256. Every other closed
form the paper prints is consistent with the corrected potential.

## Open holes

Stated from the ledgers, not from optimism; each is visible in a `[[result]]` row that is not
`solid`, in a figure annotation, or in the per-figure README section.

- **The published grid sizes depend on ``m``; ours cannot.** With the 2:1 refinement decision
  closed in every dimension, the hp-AMR sweep reproduces the published grid sequence exactly
  for ``m \le 4`` — 480, 975, 1604, 2427, 3468, 4751, 6300 points — but the reference publishes
  970 rather than 975 at grid 1 for ``m \ge 5``. Nothing in `RefineAtPuncture ∘
  IncreaseResolution ∘ 2:1 balance` reads the solution, so none of them can make a grid size
  depend on ``m``; the yaml's commented-out `TruncationError` criterion is the candidate and
  enabling it is the experiment that would settle it.
- **`fig:linsolv`'s last grid.** Our staircase on the published lattice is
  1, 7, 6, 5, 4, 3, 3, 3 (32 cumulative) against the published 1, 7, 6, 5, 4, 3, 3, 2 (31) —
  seven of eight grids exact, the finest one iteration long. The walltime panel's amplitude
  differs for a known reason: our Schwarz subdomain matrices are built eagerly and the
  reference's lazily.
- **GMRES stagnates at ``m = 19, 20`` for ``a = 0.998``, ``r_0 = 2M``** — the published
  panel 3 of `fig:modes`. Those two planes are drawn dashed, named on the figure, and carry
  their own residual; every other mode of every panel converged.
- **The error budget orders the other way from the paper's — but the lattice is now
  implicated.** Our total's numerical-resolution error is a few times LARGER than the
  finite-``m`` truncation estimator, where main.tex:723 reports resolution pushed well below
  truncation. Fig. 7 NARROWS that, without closing it: on the paper's own full-``z`` lattice
  the ordering comes out the paper's way at 15 of 18 configurations (published 18/18), against
  0 of 9 on the shipped half domain at the same AMR settings *[PRELIMINARY]*. The two are not
  the same grids — 8139 against 4814 dofs, a 1.69× confound that is disclosed rather than
  argued away — but size and rate make the confound an implausible explanation of the whole
  effect: the half domain's absolute resolution error at ``r_0 = 10M`` is a spin-INDEPENDENT
  floor of ``5.7``–``5.9\times10^{-10}``, roughly four further AMR steps from closing the
  measured gap, not 1.69× in dofs. Settling it needs the two lattices at MATCHED ndofs, which
  has not been run.
- **The compressed prograde ISCO has three members, not one.** The absolute resolution error
  runs ``2.9\times10^{-10}`` at ``a = +0.7``, ``6.0\times10^{-7}`` at ``+0.9`` and
  ``2.4\times10^{-2}`` at ``+0.998``; at ``+0.998`` the mode sum has not converged at
  ``m = 20`` either (the last retained term is 2.4 % of the total). ``a = +0.9`` is the member
  a first pass missed. Both this family and the two ``|a| = 0.998``, ``r_0 = 10M``
  configurations whose ordering also inverts are RING-MARKED on the figure by a measured rule,
  never dropped — and the layer-ratio instrument reads 0.098 at the latter two, so the honest
  statement is that the instrument does not detect the cause, not that the near-horizon layer
  is uninvolved.
- **The asymptotic ``m^{-6}`` decay of the regularized residual is not discharged.** The
  ``m \in \{2, 5\}`` slope sits on ``-6``; the now-resolved ``\{5, 10\}`` slope does not.
- **The near-horizon layer.** At high ``\nu`` the horizon element's growing polynomial order
  admits a bounded, non-analytic ``(r - r_+)^{i\nu}`` branch that the inert horizon face cannot
  reject. `horizon_hold` cures it, is chosen by a MEASURED rule (the layer ratio, never the
  spin and never ``\nu``), and defaults to **off** for reference parity — the paper's own
  observables move by ``2.6\times10^{-5}`` relative when it is switched on.

## Related pages

- [EllipticCore](../core/ellipticcore.md) — the solver stack the reproduction runs on.
- [ScalarSelfForce](../capabilities/scalarselfforce.md) — the physics, the API contracts and
  the measured tables, all stamped with the commit they were measured at.
- [Reproducibility](reproducibility.md) — environment pinning, determinism rules, the
  `reproduce/` tree contract and the suite-wide list of known holes.
- [Fixtures](fixtures.md) — the per-artifact provenance ledger (the baked `EffectiveSource`,
  Warburton–Barack and runtime-config oracles this reproduction gates against).

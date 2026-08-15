# Wigner

Hand-rolled angular-momentum coupling — the Tier-1 core stratum (C2). This package exists to
**eliminate the `WignerSymbols` package** from every GrTools runtime path: it hosts the
promoted hand-rolled implementations, so the consumers that previously hand-rolled or
`[deps]`-ed WignerSymbols inherit **no** new heavy dependency. It sits *beside* `Harmonics`
and is deliberately kept **off** the in-degree light-numeric spine.

## Capability overview

The physics this package owns:

- **3j / 6j symbols** — [`wigner3j`](@ref), [`wigner6j`](@ref): Racah closed forms, exact in
  `BigInt`/`Rational{BigInt}` combinatorics.
- **Clebsch-Gordan** — [`clebsch_gordan`](@ref): the Racah single-sum CG
  ``\langle j_1 m_1; j_2 m_2 | j_3 m_3\rangle``, the 3j↔CG counterpart. Used downstream for the
  spheroidal→spherical mode mixing in the Kerr-Lorenz reconstruction.
- **3-harmonic angular coupling** — [`coupling`](@ref): couples two spin-weighted harmonics
  ``(l_a,m_a,s_a)``, ``(l_b,m_b,s_b)`` onto a third ``(l_c,m_c,s_c)`` via two 3j symbols and the
  normalisation
  ```math
  (-1)^{m_c+s_c}\sqrt{\frac{(2l_a+1)(2l_b+1)(2l_c+1)}{4\pi}}
  \begin{pmatrix} l_c & l_a & l_b \\ s_c & -s_a & -s_b \end{pmatrix}
  \begin{pmatrix} l_c & l_a & l_b \\ -m_c & m_a & m_b \end{pmatrix}.
  ```
  This is the canonical home of the second-order Ricci source coupling.
- **Reduced Wigner-d / fixed-angle Wigner-D** — [`wignerd`](@ref) (general ``\beta``),
  [`wignerd_halfpi`](@ref) (the ``\beta=\pi/2`` collapse, trig product → ``2^{-j}``), and
  [`wignerD_fixed`](@ref) (the ``D^j_{m_1 m_2}(\pi,\pi/2,\pi/2)`` mode-mixing rotation used by
  the puncture-field machinery).

All primitives are type-generic `{T}` — pass `BigFloat` and the working precision is inherited;
the exact combinatorics run in `BigInt`/rational arithmetic and downcast only at the boundary.
Every function returns `zero(T)` on selection-rule violations (triangle, ``m``-sum,
``|m|>j`` slots) rather than throwing, matching the GSL convention consumers relied on.

## Module map

`Wigner` is a **flat module** — no submodules — assembled from four source files:

| File | Contents |
|---|---|
| `src/wigner3j.jl` | `wigner3j`, `wigner6j` (Racah closed forms) + shared `_fbig`/`_sqrt_rat` combinatorial helpers |
| `src/clebsch_gordan.jl` | `clebsch_gordan` (Racah CG) |
| `src/coupling.jl` | `coupling` (3-harmonic coupling over the in-package `wigner3j`) |
| `src/wignerD.jl` | `wignerd`, `wignerd_halfpi`, `wignerD_fixed` |

## Dependencies and poison-axis declaration

**Runtime `[deps]`: `NumericsCore` only.** This budget is CI-enforced by
`test/hygiene.jl`, which asserts `allow=["NumericsCore"]` and forbids —
statically in the TOML and, for the loadable set, again after `using Wigner`:

- `WignerSymbols` — **the defining poison axis for this package**: runtime-eliminated.
  It survives *only* as a sanctioned **test-only characterization oracle** in
  `[extras]`/`[targets]` (small, light, registered), grounding the hand-rolled 3j/6j/CG over a
  grid of ``(j_1,j_2,j_3,m_1,m_2,m_3)``.
- The other poison axes and heavy numerics: `HDF5`, `JLD2`, `Symbolics`, `OrdinaryDiffEq*`,
  `CUDA`, `AMDGPU`, `KernelAbstractions`, `Kernels`, `Data`, `SpecialFunctions`,
  `HypergeometricFunctions`, `Roots`, `NLsolve`, `SciMLBase`, `QuadGK`.

`TestKit` (+ its `Data` spine) appears in `[extras]` purely to power the hygiene gates.

## GPU / data needs

None. Pure CPU integer/rational combinatorics plus `T`-precision trig and square roots — no GPU
kernels, no runtime data files. The only bundled artifact is the test-side
`test/wignerD_reference.csv` (Wigner-d oracle values).

## Provenance

Promoted via **rewrite**, not copy; the predecessor codes served as the characterization oracle
during construction and were never edited:

| Legacy source | Legacy function | New home |
|---|---|---|
| `KerrLorenzCirc/src/angular.jl` | `cg` | `Wigner.clebsch_gordan` |
| `SecondOrderRicci/src/coupling.jl` | `coupling` (had the only real `WignerSymbols` `[deps]`) | `Wigner.coupling` |
| `Punctures/src/mode_mixing.jl` | `_wignerd_halfpi`, `_wignerD_fixed` | `Wigner.wignerd_halfpi`, `Wigner.wignerD_fixed` |

The 3j/6j Racah closed forms are new canonical implementations grounded on analytic identities
(orthogonality, Regge symmetries, triangle/selection rules, 3j↔CG) plus the test-only
`WignerSymbols.jl` oracle; the Wigner-d rotation matrices are grounded on the closed forms and a
Mathematica `WignerD` reference.

## What changed from the predecessor codes

- **`WignerSymbols` dropped from the runtime graph**: SecondOrderRicci's `coupling` is
  re-implemented over the in-package `wigner3j` (originally GSL `gsl_sf_coupling_3j`), so its
  descendant `SelfForce.SecondOrder` carries no `WignerSymbols` dependency
  (characterization-gated against WignerSymbols output).
- **Hard `Float64` → type-generic `{T}`**: the legacy KerrLorenzCirc `cg` and Punctures
  ``\pi/2`` routines were `Float64`-only; every promoted primitive now takes a leading
  `::Type{T}` and inherits BigFloat precision. Legacy-signature `Float64` convenience methods
  are retained for facade re-export.
- **General-``\beta`` consolidation**: legacy only had the ``\beta=\pi/2`` special case;
  `wignerd` implements the general reduced Wigner-d, with `wignerd_halfpi` recovered as its
  special case.

## Public API

```@autodocs
Modules = [Wigner]
Order = [:module, :type, :constant, :function]
```

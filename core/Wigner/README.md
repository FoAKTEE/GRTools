# Wigner (GrToolWigner) — core stratum C2

Hand-rolled angular-momentum coupling for GrTool. This package exists to **eliminate the
`WignerSymbols` package** from every runtime path (design §5.3 / §7.4 / L4): it hosts the
promoted hand-rolled implementations so the consumers that previously hand-rolled or
`[deps]`-ed Wigner inherit **no** new heavy dependency.

It is a Tier-1 core stratum whose **only** dependency is `NumericsCore`. It sits *beside*
`Harmonics` and is deliberately kept off the in-degree light-numeric spine.

## Physics this capability owns

- **3j / 6j symbols** (`wigner3j`, `wigner6j`) — Racah closed form, BigInt-exact.
- **Clebsch-Gordan** (`clebsch_gordan`) — Racah CG; the 3j↔CG counterpart.
- **3-harmonic angular coupling** (`coupling`) — couples two spin-weighted harmonics onto a
  third via two 3j symbols and the √((2lₐ+1)(2l_b+1)(2l_c+1)/4π) normalisation; the canonical
  home of the second-order Ricci source coupling.
- **Reduced Wigner-d / fixed-angle Wigner-D** (`wignerd`, `wignerD_halfpi`) — general β plus the
  β=π/2 special case and the D(π,π/2,π/2) mode-mixing rotation.

All primitives are type-generic `{T}` (BigFloat precision inherited); the combinatorics run in
`BigInt`/rational arithmetic and downcast only at the boundary.

## Legacy packages absorbed (via REWRITE/promotion, not copy)

| Legacy source | Function | New home |
|---|---|---|
| `KerrLorenzCirc/src/angular.jl` | `cg` | `Wigner.clebsch_gordan` |
| `SecondOrderRicci/src/coupling.jl` | `coupling` (drops `WignerSymbols` dep) | `Wigner.coupling` |
| `Punctures/src/mode_mixing.jl` | `_wignerd_halfpi`, `_wignerD_fixed` | `Wigner.wignerD` |

The legacy packages stay green as the characterization oracle; they are never edited.

## Public API (built TDD in the BUILD rounds)

`wigner3j`, `wigner6j`, `clebsch_gordan`, `coupling`, `wignerd`, `wignerD_halfpi`.

## Dependencies

- Runtime: `NumericsCore` only. **No `WignerSymbols`.**
- Test-only oracle: `WignerSymbols.jl` (small, light, registered) in `[extras]`/`[targets]` to
  ground the hand-rolled 3j/6j/CG over a grid of `(j1,j2,j3,m1,m2,m3)`.

## Grounding

Analytic 3j/6j identities (orthogonality, Regge symmetries, closed forms, triangle/selection
rules, 3j↔CG) plus test-only characterization against `WignerSymbols.jl`, and the Wigner-d
closed forms / Mathematica `WignerD` oracle for the rotation matrices.

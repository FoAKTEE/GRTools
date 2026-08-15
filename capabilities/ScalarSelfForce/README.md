# ScalarSelfForce — elliptic scalar self-force capability

The m-mode **elliptic** Kerr scalar self-force capability: the GrTools
reproduction of arXiv:2606.04998 (Vu et al., *Self-force calculations with
numerical relativity methods*). It owns the physics of the elliptic formulation
and consumes the physics-free hp-DG solver stack from `core/EllipticCore`.

## Physics this capability owns

- **Orbits/** — `CircularOrbit{T}` thin adapters over Spacetime (E, L, Ω, uᵗ,
  u^φ, r_ISCO, r_*, φ̄ offset).
- **MModeOperator/** — the operator pipeline: mode ansatz → t-sliced △̃_m →
  hyperboloidal slicing → vtu junctions → horizon-penetrating first-order flux
  form (Fⁱ, β_m, γⁱ_m) → equatorial-symmetry z² variant → Bayliss–Turkel outer
  radiation boundary condition.
- **Puncture/** — 4th-order Detweiler–Whiting Φ^P (A_ijk tensor basis), m-mode
  integrals over `BigSpecial.legendre_P3`, analytic derivatives.
- **Domain/** — `RunConfig` runtime configuration (`crp_paper_config`,
  `CRP_PAPER_MMODES`) + EllipticCore `Mesh2D` block instantiation.

Planned (each landing tests-first once grounded):
`EffectiveSource/` (S^eff_m direct + near-worldline series switch),
`Worldtube/` (Ψ^P mortar jump data, fixed-source assembly), `Solve/`
(single-mode driver, hp-AMR loop, process-parallel m fan-out), `Extraction/`
(F_r^m, Heffernan regularization, mode-sum + tail).

## Conventions (locked package-wide)

- Geometric units `G = c = 1`, background mass `M = 1`; Kerr spin
  `a ∈ (−1, 1)`, retrograde orbits encoded via `a < 0` with `m ≥ 0` and the
  ×2 real-part mode-sum convention.
- Mode ansatz `Φ = Σ_m Ψ_m(r, z) · sin^{|m|}θ / r · e^{i m (φ̄ − Ω t)}` with
  `z = cosθ` and the horizon-regular azimuthal chart
  `φ̄ = φ + a/(r₊−r₋)·ln((r−r₊)/(r−r₋))`.
- Radial chart: horizon-penetrating Boyer–Lindquist `r` via
  `∂_{r_*} = Δ/(r²+a²) ∂_r`; vtu slicing (v→t→u null/Boyer slabs) with
  piecewise boost `H ∈ {−1, 0, +1}` and junction conditions
  `[∂_{r_*}Ψ] = −imΩΨ`.
- Equatorial symmetry via the z² variant on `z ∈ [0, 1]`.

## Dependencies

Binding dependency budget, enforced by `test/hygiene.jl`: exactly
`{LinearAlgebra, Printf, NumericsCore, BigSpecial, Spacetime, EllipticCore,
Data}`. Special functions come only through BigSpecial; the EffectiveSource C
oracles are baked test fixtures, never live deps (no GSL/ccall channel); GPU
stacks are banned.

## Type policy

Type-generic `{T}` throughout so BigFloat precision is inherited from
NumericsCore/BigSpecial (never hard-Float64'd); the production path is
ComplexF64.

## Testing

`Pkg.test` standalone. Grounding oracles: baked EffectiveSource-C fixtures,
Warburton–Barack reference tables, spectre-verified operator fixtures, and
exact analytic identities; the dependency budget is itself a test
(`test/hygiene.jl`).

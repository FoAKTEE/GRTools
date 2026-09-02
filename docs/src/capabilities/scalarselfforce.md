# ScalarSelfForce

The m-mode **elliptic** Kerr scalar self-force capability: the GrTools
reproduction of arXiv:2606.04998 (Vu et al., *Self-force calculations with
numerical relativity methods*). It owns the physics of the elliptic
formulation and consumes the physics-free hp-DG solver stack from
[EllipticCore](../core/ellipticcore.md).

## Capability overview

- **Orbits** — `CircularOrbit{T}` thin adapters over Spacetime: E, L, Ω, uᵗ,
  u^φ, r_ISCO, tortoise r_*, and the horizon-regular azimuthal offset φ̄.
- **MModeOperator** — the operator pipeline: mode ansatz → t-sliced
  ``\tilde\triangle_m`` → hyperboloidal slicing → vtu junctions →
  horizon-penetrating first-order flux form (``F^i``, ``\beta_m``,
  ``\gamma^i_m``) → equatorial-symmetry ``z^2`` variant → Bayliss–Turkel
  outer radiation boundary condition.
- **Puncture** — 4th-order Detweiler–Whiting ``\Phi^P`` (``A_{ijk}`` tensor
  basis), m-mode integrals over `BigSpecial.legendre_P3`, analytic
  derivatives.
- **Domain** — `RunConfig` runtime configuration (`crp_paper_config`,
  `CRP_PAPER_MMODES`) and EllipticCore `Mesh2D` block instantiation.
- **Extraction** — Heffernan m-mode regularization parameters of the radial
  self-force for circular equatorial Kerr orbits (arXiv:1211.6446), in the
  folded ``m \ge 0`` ×2 real-part convention:
  - `reg_param_Fr(order, orbit)` — the m-independent circular-Kerr
    coefficients ``F_{r[2]}``, ``F_{r[4]}`` (`HEFFERNAN_CIRCULAR_FR_ORDERS`),
    purely algebraic in ``(a, r_0, M)``; the printed ``a \to 0`` Schwarzschild
    limits are given in the docstring.
  - `mmode_reg_factor_over_pi(order, m)` — the printed m-dependence
    ``\pi F^m_{r[2j]} / F_{r[2j]}`` as an exact `Rational{BigInt}`, for every
    published order `HEFFERNAN_MMODE_ORDERS = (2, 4, 6, 8)`;
    `mmode_reg_factor([T,] order, m)` is the same divided by ``\pi``.
  - `mmode_reg_Fr(order, orbit, m)` — the assembled per-mode subtraction term
    ``F^m_{r[order]}``.
  - `mmode_reg_factor_tail_over_pi(order, mstart)` /
    `mmode_reg_factor_tail([T,] order, mstart)` /
    `mmode_reg_tail_Fr(order, orbit, mstart)` — the closed-form analytic tail
    ``\sum_{m \ge m_{\rm start}} F^m_{r[order]}`` that completes a truncated
    mode sum.

Planned, each landing written test-first once it is grounded:
`EffectiveSource` (``S^{\rm eff}_m`` direct + near-worldline series switch),
`Worldtube` (``\Psi^P`` mortar jump data, fixed-source assembly), `Solve`
(single-mode driver, hp-AMR loop, process-parallel m fan-out), and the rest of
`Extraction` (``F_r^m`` extraction from a solved mode, and the mode sum itself).

## Conventions (locked package-wide)

- Geometric units ``G = c = 1``, background mass ``M = 1``; Kerr spin
  ``a \in (-1, 1)``, retrograde orbits encoded via ``a < 0`` with ``m \ge 0``
  and the ×2 real-part mode-sum convention.
- Mode ansatz
  ``\Phi = \sum_m \Psi_m(r, z)\, \sin^{|m|}\theta / r\, e^{i m (\bar\varphi - \Omega t)}``
  with ``z = \cos\theta`` and the horizon-regular azimuthal chart
  ``\bar\varphi = \varphi + \frac{a}{r_+ - r_-}\ln\frac{r - r_+}{r - r_-}``.
- Radial chart: horizon-penetrating Boyer–Lindquist ``r`` via
  ``\partial_{r_*} = \frac{\Delta}{r^2 + a^2}\partial_r``; vtu slicing
  (v→t→u null/Boyer slabs) with piecewise boost ``H \in \{-1, 0, +1\}`` and
  junction conditions ``[\partial_{r_*}\Psi] = -i m \Omega \Psi``.
- Equatorial symmetry via the ``z^2`` variant on ``z \in [0, 1]``.

## Dependencies and poison-axis declaration

Binding dependency budget, enforced by `test/hygiene.jl`: exactly
`{LinearAlgebra, Printf, NumericsCore, BigSpecial, Spacetime, EllipticCore,
Data}`. Special functions come only through BigSpecial (no direct
SpecialFunctions / HypergeometricFunctions); the EffectiveSource C oracles are
baked test fixtures, never live deps (no GSL/ccall channel); GPU stacks are
banned.

## Type policy

Type-generic `{T}` throughout so BigFloat precision is inherited from
NumericsCore/BigSpecial (never hard-Float64'd); the production path is
ComplexF64.

## Public API

```@autodocs
Modules = [ScalarSelfForce]
Order = [:module, :type, :constant, :function]
```

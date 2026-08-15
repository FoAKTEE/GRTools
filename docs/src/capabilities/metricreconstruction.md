# MetricReconstruction

Tier-4 capability (design §4.4): reconstruction of the linear metric perturbation
``h_{\mu\nu}`` of Kerr from perturbed-field scalars. It fuses two legacy packages into
disjoint submodule namespaces so the gauge-dependent concept-name overlaps
(`weyl_scalar`, `to_BLS`) coexist without rename.

## Physics overview

- **Radiation-gauge CCK reconstruction** (`MetricReconstruction.Radiation`): the
  Chrzanowski–Cohen–Kegeles / Hertz-potential procedure in the IRG (``s=-2``) and ORG
  (``s=+2``) gauges. From the Teukolsky ``\psi_4`` amplitudes it builds the Hertz
  potential, the six tetrad-frame `H`-components, the tetrad→coordinate projection to
  ``h_{\mu\nu}``, the linearized-Einstein vacuum residual ``G^{(1)}_{\mu\nu}[h]``, and the
  gauge-invariant linearized ``\psi_0/\psi_4`` regenerated from the reconstructed ``h``
  (the physics gate: ``\psi_4[h_\mathrm{ORG}]`` must equal the Teukolsky ``s=-2`` shape
  ``\rho^4\,R^\mathrm{up}\,S_{-2}``).
- **Lorenz-gauge reconstruction** (`MetricReconstruction.Lorenz`): the Dolan-2024
  ``a\neq 0`` Lorenz-gauge metric reconstruction for a circular equatorial Kerr orbit —
  angular/radial closure (`closure_kerr`/`closure_a0`), App-H jump conditions
  (`jumps_pm`/`hj1`), the radial fields ``P_{-2}, P_{-1}, h, \kappa``
  (`reconstruct_kerr`/`reconstruct`), ten-component tetrad assembly and BLS conversion
  (`tetrad_modes_kerr`/`to_BLS_kerr` and their ``a=0`` counterparts), and the
  gauge-invariant ``\psi_4`` cross-check (`weyl_psi4_lorenz`).
- **`Lorenz.LorenzRadial`**: the CANONICAL home of the ``a=0`` first-order Lorenz-gauge
  radial field-equation RHS (`even_ode!`/`odd_ode!`/`monopole_ode!`/`odd_dipole_ode!`,
  plus `even_dipole_ode!`/`even_static_ode!`) that `SelfForce.LorenzField` and
  `AccretionDisc.Force` consume but do not own (design §5.3 leak fix L3, Rev-2
  physics-justified single home).

Every primitive is type-generic `{T}`, so BigFloat precision (256..2048-bit) is inherited
from the input element type.

## Submodule map

| Module | Contents |
|---|---|
| `MetricReconstruction` | Namespace root; re-exports `Radiation` and `Lorenz`. |
| `MetricReconstruction.Radiation` | CCK/IRG-ORG: `SepField`/`RTJet` separable ``(r,\theta)`` jet algebra, `Dn`/`Ln` GHP ladder operators, per-mode radial/angular jets, six Hertz tetrad components `HIRG_*`/`HORG_*` (+`_jet` forms), Teukolsky–Starobinsky `weyl_amplitudes`, tetrad→coordinate projection, `lin_einstein_residual`, regenerated `weyl_scalar`, `h_coord_modesum`/`einstein_grid` mode-sum surface, `reconstruct`/`psi` entry points. |
| `MetricReconstruction.Lorenz` | Dolan-2024: Chandrasekhar operator scalars + 1D radial jet algebra (`jet_ops`), `OrbitBackground`, jumps (`jumps_pm`/`hj1`), closures (`closure_kerr`/`closure_a0`), radial fields (`reconstruct_kerr`/`reconstruct`, `HomFieldK`/`KappaFieldK` and ``a=0`` twins), S-matrices (`smatrices_kerr`), ten-component assembly (`tetrad_modes_kerr`/`to_BLS_kerr`, `build_recon`/`to_BLS`), `weyl_psi4_lorenz`. |
| `MetricReconstruction.Lorenz.LorenzRadial` | Canonical ``a=0`` Lorenz-gauge radial field-equation RHS (BLS tensor-harmonic systems, tortoise variable), pure arithmetic, `{T}`-generic. |

The two jet algebras (Radiation's 2D separable `RTJet`/`SepField` vs Lorenz's 1D radial
Chandrasekhar towers) are DIFFERENT decisions and are deliberately kept separate (§4.4).

## Dependencies and poison-axis budget

Runtime `[deps]` (Project.toml): GrTools `NumericsCore`, `BigSpecial`, `Wigner`,
`Spacetime`, `Harmonics`, `Perturbation`; external `QuadGK`; stdlib `LinearAlgebra`,
`Printf`. `BigSpecial` is a SANCTIONED Rev-3 amendment on top of the §3.2 tier line
(terminating-``{}_2F_1``/``\Gamma`` primitives used by the Lorenz submodule).

The binding budget is enforced by `test/hygiene.jl` (§7.7 gates, run in CI): the allow
list is exactly the set above, and the package must NOT carry — statically in
`Project.toml` and dynamically after `using MetricReconstruction` — any of `HDF5`, `JLD2`,
`Symbolics`, `OrdinaryDiffEq*`, `CUDA`, `AMDGPU`, `KernelAbstractions`, `Kernels`, `Data`,
`WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase`. This capability is the HDF5-FREE
reconstruction surface of the DAG (design §4.4 "EXCLUDED: HDF5, Symbolics,
KernelAbstractions"); `AccretionDisc` reaches the Lorenz fields through it without
inheriting the SelfForce HDF5/Kernels closure.

## GPU / data needs

None: no GPU path and ZERO runtime I/O. The package is BigFloat-CPU heavy
(log-Chebyshev-barycentric Kerr radial interpolants, VoP/Abel quadrature, banded closure
solves); per-``(l,m)`` work is threaded with lock-guarded caches. Legacy test fixtures are
baked as literals in `test/` — the legacy oracles are not live test dependencies.

## Legacy provenance (drawsFrom, design §6.4)

| Legacy pkg / file / functions | New home |
|---|---|
| `MetricReconstruction` / `api.jl`, `operators.jl`, `modes.jl`, `tetrad_components*.jl`, `weyl_amplitudes.jl`, `projection*.jl`, `einstein.jl`, `weyl_scalars.jl` | `MetricReconstruction.Radiation` |
| `KerrLorenzCirc` / `assembly.jl`, `operators.jl`, `angular.jl`, `background.jl`, `jumps.jl`, `closure.jl`, `radial.jl`, `radial_kerr.jl`, `components*.jl` | `MetricReconstruction.Lorenz` |
| `KerrLorenzCirc` / `angular.jl` / `cg` | promoted to `Wigner.clebsch_gordan` (a thin `cg` adapter is retained in `Lorenz`) |
| `H1Lorenz` / `ode_systems.jl` + `DiscModes` / `lorenz_metric.jl` / `even_ode!`, `odd_ode!`, `monopole_ode!`, `odd_dipole_ode!` | `MetricReconstruction.Lorenz.LorenzRadial` (canonical; fork collapsed) |
| `MetricReconstruction` / `background.jl`, `einstein.jl` / `kerr_metric`, `kinnersley`, `christoffel_BL`, `_dg` | delegated to `Spacetime` geometry |

## What changed vs legacy

- **Namespacing, no rename**: legacy `MetricReconstruction` and `KerrLorenzCirc` shared
  gauge-dependent concept names (`weyl_scalar`, `to_BLS`); they now live as
  `Radiation.weyl_scalar` vs the Lorenz variants (`Lorenz.weyl_psi4_lorenz`,
  `Lorenz.to_BLS`/`to_BLS_kerr`) in disjoint submodules.
- **`monopole_ode!` arity seam reconciled**: H1Lorenz's `(g,rs,y,l,ω)` carried a dead `l`;
  DiscModes' `(g,rs,y,ω)` is the canonical arity-4 primitive, with a trivial arity-5
  method retained for the H1Lorenz consumer convention (bodies verified bit-identical
  before the fork collapse, per the §4.4 gate).
- **Promotions/delegations**: KLC's hand-rolled Clebsch–Gordan is promoted to
  `Wigner.clebsch_gordan`; background geometry (``\Delta,\Sigma,\zeta``, metric,
  Kinnersley legs, Christoffels) delegates to `Spacetime`; Teukolsky radial legs come from
  `Perturbation`; the LorenzRadial background scalars (`f`, `f'`, ``V_\ell``, tortoise
  inverse) delegate to `Spacetime`.
- **Precision**: the whole reconstruction chain is `{T}`-generic and inherits BigFloat
  precision (256..2048-bit) from the input element type; where the Float64 Kerr
  ``R^\mathrm{up}`` magnitude overflows, the radial legs switch to log-Chebyshev-barycentric
  BigFloat interpolants sampled once per mode.

## Public API

```@autodocs
Modules = [MetricReconstruction, MetricReconstruction.Radiation, MetricReconstruction.Lorenz, MetricReconstruction.Lorenz.LorenzRadial]
Order = [:module, :type, :constant, :function]
```

# MetricReconstruction

Tier-4 capability: reconstruction of the linear metric perturbation
``h_{\mu\nu}`` of Kerr from perturbed-field scalars. It fuses two predecessor codes into
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
  `AccretionDisc.Force` consume but do not own (the physics-justified single home).

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
Chandrasekhar towers) are DIFFERENT decisions and are deliberately kept separate.

## Dependencies and poison-axis budget

Runtime `[deps]` (Project.toml): GrTools `NumericsCore`, `BigSpecial`, `Wigner`,
`Spacetime`, `Harmonics`, `Perturbation`; external `QuadGK`; stdlib `LinearAlgebra`,
`Printf`. `BigSpecial` is a sanctioned amendment on top of the dependency tier line
(terminating-``{}_2F_1``/``\Gamma`` primitives used by the Lorenz submodule).

The binding budget is enforced by `test/hygiene.jl` (run in CI): the allow
list is exactly the set above, and the package must NOT carry — statically in
`Project.toml` and dynamically after `using MetricReconstruction` — any of `HDF5`, `JLD2`,
`Symbolics`, `OrdinaryDiffEq*`, `CUDA`, `AMDGPU`, `KernelAbstractions`, `Kernels`, `Data`,
`WignerSymbols`, `Roots`, `NLsolve`, `SciMLBase`. This capability is the HDF5-FREE
reconstruction surface of the dependency graph (HDF5, Symbolics and KernelAbstractions
are excluded by design); `AccretionDisc` reaches the Lorenz fields through it without
inheriting the SelfForce HDF5/Kernels closure.

## GPU / data needs

None: no GPU path and ZERO runtime I/O. The package is BigFloat-CPU heavy
(log-Chebyshev-barycentric Kerr radial interpolants, VoP/Abel quadrature, banded closure
solves); per-``(l,m)`` work is threaded with lock-guarded caches. Legacy test fixtures are
baked as literals in `test/` — the legacy oracles are not live test dependencies.

## Provenance

| Predecessor code / file / functions | New home |
|---|---|
| `MetricReconstruction` / `api.jl`, `operators.jl`, `modes.jl`, `tetrad_components*.jl`, `weyl_amplitudes.jl`, `projection*.jl`, `einstein.jl`, `weyl_scalars.jl` | `MetricReconstruction.Radiation` |
| `KerrLorenzCirc` / `assembly.jl`, `operators.jl`, `angular.jl`, `background.jl`, `jumps.jl`, `closure.jl`, `radial.jl`, `radial_kerr.jl`, `components*.jl` | `MetricReconstruction.Lorenz` |
| `KerrLorenzCirc` / `angular.jl` / `cg` | promoted to `Wigner.clebsch_gordan` (a thin `cg` adapter is retained in `Lorenz`) |
| `H1Lorenz` / `ode_systems.jl` + `DiscModes` / `lorenz_metric.jl` / `even_ode!`, `odd_ode!`, `monopole_ode!`, `odd_dipole_ode!` | `MetricReconstruction.Lorenz.LorenzRadial` (canonical; fork collapsed) |
| `MetricReconstruction` / `background.jl`, `einstein.jl` / `kerr_metric`, `kinnersley`, `christoffel_BL`, `_dg` | delegated to `Spacetime` geometry |

## What changed from the predecessor codes

- **Namespacing, no rename**: the earlier `MetricReconstruction` and `KerrLorenzCirc` research codes shared
  gauge-dependent concept names (`weyl_scalar`, `to_BLS`); they now live as
  `Radiation.weyl_scalar` vs the Lorenz variants (`Lorenz.weyl_psi4_lorenz`,
  `Lorenz.to_BLS`/`to_BLS_kerr`) in disjoint submodules.
- **`monopole_ode!` arity seam reconciled**: H1Lorenz's `(g,rs,y,l,ω)` carried a dead `l`;
  DiscModes' `(g,rs,y,ω)` is the canonical arity-4 primitive, with a trivial arity-5
  method retained for the H1Lorenz consumer convention (bodies verified bit-identical
  before the fork collapse).
- **Promotions/delegations**: KLC's hand-rolled Clebsch–Gordan is promoted to
  `Wigner.clebsch_gordan`; background geometry (``\Delta,\Sigma,\zeta``, metric,
  Kinnersley legs, Christoffels) delegates to `Spacetime`; Teukolsky radial legs come from
  `Perturbation`; the LorenzRadial background scalars (`f`, `f'`, ``V_\ell``, tortoise
  inverse) delegate to `Spacetime`.
- **Precision**: the whole reconstruction chain is `{T}`-generic and inherits BigFloat
  precision (256..2048-bit) from the input element type; where the Float64 Kerr
  ``R^\mathrm{up}`` magnitude overflows, the radial legs switch to log-Chebyshev-barycentric
  BigFloat interpolants sampled once per mode.

## Radial node-count options (`KerrRadialOptions`)

Where the Float64 ``R^\mathrm{up}`` magnitude overflows, `Lorenz.branchK` replaces the affected
legs with log-Chebyshev-barycentric panels sampled once per mode. Five integers set how many nodes
those panels carry — a pure cost/accuracy dial, but one that CHANGES the numbers a run produces:

| field | default | meaning |
|---|---|---|
| `ncheb_floor` | `30` | nodes on the Float64-gate Up panel (`_branchK_up_interp`) |
| `nin_floor` | `18` | In-leg Chebyshev-node floor (18 nodes reach full precision, measured) |
| `nin_factor` | `3` | In-leg node law `N_in = max(nin_floor, nin_factor*l)` (before the √-span bump) |
| `nup_floor` | `70` | BigFloat ``s=0`` Up-leg node floor |
| `nup_factor` | `3` | BigFloat ``s=0`` Up-leg node law `N = max(nup_floor, nup_factor*l)` |

The defaults are the F8-validated production values. Override them per call, or for the process:

```julia
using MetricReconstruction: Lorenz
rec = Lorenz.reconstruct_kerr(0.5, 6.0, 2; opts = Lorenz.KerrRadialOptions(; nin_floor = 24))

old = Lorenz.kerr_radial_options()
Lorenz.set_kerr_radial_options!(; nup_factor = 4)   # keyword form patches the current default
Lorenz.set_kerr_radial_options!(old)                # restore
```

These five knobs were process environment variables (`KLC_NCHEB_FLOOR`, `KLC_NIN_FLOOR`,
`KLC_NIN_FACTOR`, `KLC_NUP_FLOOR`, `KLC_NUP_FACTOR`) through v0.1.0. They no longer exist: a
numerical result that silently depends on the environment of the process that produced it cannot be
reproduced from the tree, so the knobs are now explicit arguments and `test/lorenz_kerr_radial.jl`
gates both the defaults and the absence of any `KLC_*` read in `src/`.

## Public API

```@autodocs
Modules = [MetricReconstruction, MetricReconstruction.Radiation, MetricReconstruction.Lorenz, MetricReconstruction.Lorenz.LorenzRadial]
Order = [:module, :type, :constant, :function]
```

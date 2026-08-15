# PostNewtonian

Tier-4 GrTools capability — the **post-Newtonian series engine** and the **Symbolics quarantine**.
It fuses the two legacy PN leaf packages (`TeukolskyPN` + `PostNewtonianSelfForce`)
into one capability and is a top-of-DAG leaf: nothing depends UP on it, so `Symbolics` (a HARD dep
here and ONLY here) never leaks its ~1.9 s precompile floor onto the numeric/GPU waveform path.

## Physics owned

PN expansions of Teukolsky radial/source modes and of the PN self-force, carried as
`Symbolics.Num`-valued Laurent/Taylor series. The central object is `PNSeries`, the Julia mirror of
Mathematica `SeriesData`:

```math
\mathrm{SeriesData}[x, x_0, \{a_0,\dots,a_k\}, n_{\min}, n_{\max}, d]
\;=\; \sum_{i=0}^{k} a_i\,(x-x_0)^{(n_{\min}+i)/d} \;+\; O\!\left[(x-x_0)^{n_{\max}/d}\right],
```

with ``k = (n_{\max}-n_{\min})/d - 1`` stored coefficients. On top of this the capability builds the
symbolic Mano–Suzuki–Takasugi (MST) machinery — the renormalized angular momentum ``\nu`` solve and
the ``a_n`` coefficients as ``\varepsilon = 2\omega`` series (an ``\varepsilon``-power IS a PN
order) — the ``\eta``-graded radial In/Up solutions, the point-particle source assembly with
asymptotic amplitudes ``Z^\infty/Z^H`` and energy flux (Schwarzschild ``\ell=2,m=2`` recovers the
Newtonian quadrupole ``\tfrac{32}{5}v^{10}``), and the PN self-force series store with its
``y^n \log^{n_L}\! y`` coefficient extraction and resummed series.

The two quarantine anchors are struct FIELD types — `PNSeries.coeffs :: Vector{Symbolics.Num}` and
`ResummedSeriesData.sf :: Symbolics.Num` — which is precisely why `Symbolics` cannot be a weakdep
and must be a hard dep confined to this capability.

## Submodule map

Seven submodules (verified against `src/`); load order `Series → Catalog → MMAInputForm →
SelfForce → MST → RadialMode → SourceMode`.

| Submodule | Owns |
|-----------|------|
| `Series` | the `PNSeries` Num-algebra engine: Laurent/Taylor arithmetic (`+ - * / inv ^`), `Normal` materialization, `polyToSeries`, `series_take/collect/terms/of`, change-of-parameter, power-counting, `expand_log`/`expand_gamma` |
| `Catalog` | PN `.m` discovery + substring browse over the owned SeriesData store (`SERIES_DATA_DIR`, `pn_series_catalog`, `select_pn_series`, `pn_browse`); stdlib `walkdir` only |
| `MMAInputForm` | the unified MMA-InputForm→`Num` symbolic walker (`read_mma_inputform`, `parse_pn_data_file`): `Meta.parse` + tree walk, NO `eval`/world-age; covers the full circular/spherical/eccentric/Kerr alphabet incl. `ResummedSeriesData`, `HoldForm` generating-function heads, `TagBox` box-language, and balanced reduction for the ~35k-term flux giants |
| `SelfForce` | `PostNewtonianData` wrapper + order metadata (`pn_series`, `pn_start_order`, `pn_max_order`), `postnewtonian_coefficient` (the ``y^n \log^{n_L}\! y`` extraction), `ResummedSeriesData` (extends `Series.Normal` on the same method table), `postnewtonian_expansion` load-by-name |
| `MST` | symbolic MST coefficient engine: truncated-``\varepsilon`` Laurent algebra `ES{T}`, the ``\nu``-solve from the pole-cleared characteristic equation, the ``a_n`` ``\varepsilon``-series lifted to `PNSeries` (``\eta``-form and ``\gamma``-FreqRep form), `imagunit` |
| `RadialMode` | ``\eta``-graded symbolic radial Coulomb/MST In/Up solutions (`TeukolskyRadialPN`) + asymptotic mode-amplitude builders, from `MST`'s ``\varepsilon``-series; adds the transcendental series ops (exp/log1p/Pochhammer/Γ-ratio/Hurwitz-ζ/sin π) the Coulomb construction needs |
| `SourceMode` | point-particle source ``\{\delta,\delta',\delta''\}`` coefficient series (exact in ``r_0,\Omega``, keyed by ``(l,m)``), sourced mode assembly, asymptotic energy flux (`TeukolskyPointParticleModePN`) |

The imaginary unit is carried as the Symbolics symbol `imagunit()` (`𝕚`) so a coefficient stays a
single `Num` (a `Complex{Num}` cannot live in `PNSeries`' `Vector{Num}`); substitute `=> im` to
evaluate numerically.

## Dependencies and poison-axis budget

Declared `[deps]` (Project.toml): `NumericsCore`, `Harmonics`, `Data` (GrTools) + `Symbolics`
(HARD, quarantined), `SpecialFunctions`, and stdlib `Printf`/`LinearAlgebra`.

The binding budget is the CI gate in `test/hygiene.jl` (`TestKit.assert_dep_budget` +
`assert_not_loaded`):

- **may carry** (allow-list): `Printf`, `LinearAlgebra`, `NumericsCore`, `Harmonics`, `Data`,
  `Symbolics`, `SpecialFunctions`. (`LinearAlgebra` is flagged there as a TODO discrepancy:
  stdlib, zero-cost, but grep-verified no `src/` usage.)
- **must declare**: `Symbolics` — this capability is its SOLE hard home in GrTools (gate 5,
  quarantine home).
- **may NOT carry** (forbid-list): `HDF5`, `JLD2`, `CUDA`, `AMDGPU`, `KernelAbstractions`,
  `Kernels`, `WignerSymbols`, `OrdinaryDiffEq*`, `Roots`, `NLsolve`, `SciMLBase`. In particular
  PostNewtonian touches `Data` only for `.m`/`.wl` text artifacts, so it must NOT declare `HDF5`
  (Data's HDF5 weakdep-ext stays dormant here; hygiene gate 5).
- the closure-load gate additionally asserts none of the poison packages the Symbolics closure
  cannot legitimately pull are loaded after `using PostNewtonian`.

The quarantine holds because the DAG is one-directional: `Harmonics` supplies the numeric
`SpinWeightedSpheroidalEigenvalueSeries` (returns `Vector{Float64}`, NOT `Num`, so it carries no
Symbolics), and nothing depends up on PostNewtonian after the dead `FEW → TeukolskyPN` edge cut.

## GPU / data needs

- **GPU: none.** This capability never touches `Kernels`/`CUDA` (forbidden by the hygiene gate).
  Its parallelism model is CPU threads/processes over the mode axis.
- **Data: the package-OWNED `SeriesData` store** — 91 Mathematica-InputForm `.m` files (~38 MB,
  `Kerr/` + `Schwarzschild/` trees), a real copy (no symlink into the predecessor code),
  resolved via `Catalog.SERIES_DATA_DIR`. The externalize-to-`Data.artifacts`
  move is deferred. No HDF5/JLD2 grids.
- **Test grounding** needs no heavy live legacy dep: the wolframscript byte-exact oracle on the 91
  `.m` files, captured fixtures (`test/fixtures/pn_series_ops_oracle.json`), and analytic PN
  identities.

## Provenance

The predecessor codes served as the characterization oracle during construction; their outputs
are baked into `test/fixtures/`.

| Legacy source | Landed in |
|---------------|-----------|
| TeukolskyPN / `pnseries.jl` — `PNSeries`, `Normal`, arithmetic, `change_series_parameter`, `power_counting`, `expand_log`/`expand_gamma` | `PostNewtonian.Series` |
| TeukolskyPN / `mst_coefficients.jl` — `MSTCoefficients`, ν-solve, ε-series, `imagunit` | `PostNewtonian.MST` |
| TeukolskyPN / `radial_pn.jl`, `mode_pn*.jl`, `mode_pn_assemble.jl` — `TeukolskyRadialPN`, `TeukolskyPointParticleModePN` | `PostNewtonian.RadialMode` + `PostNewtonian.SourceMode` |
| PostNewtonianSelfForce / `postnewtonian_data.jl`, `coefficients.jl`, `resummed_series.jl` | `PostNewtonian.SelfForce` |
| PostNewtonianSelfForce / `catalog_selector.jl`, `series_data.jl` | `PostNewtonian.Catalog` |
| TeukolskyPN `mode_pn_source.jl` `_eval_mma_ast` + PostNewtonianSelfForce `series_data.jl` `_eval_pn_ast` (duplicated walker forks) | unified into `PostNewtonian.MMAInputForm.read_mma_inputform`, behind a byte-for-byte `Num`-equivalence gate over all 91 legacy `.m` files |

## What changed from the predecessor codes

- **One capability instead of two packages.** `TeukolskyPN` and `PostNewtonianSelfForce` fuse;
  `Symbolics` gets a single hard home and the dead `FEW → TeukolskyPN` import edge is cut, so no
  waveform-path package pays the Symbolics precompile floor.
- **The MMA walker fork is unified.** The two duplicated InputForm walkers (`_eval_mma_ast` /
  `_eval_pn_ast`) collapse into one `MMAInputForm.read_mma_inputform`, verified `Num`-equivalent
  on all 91 legacy `.m` files. It is DISTINCT from `Data.read_mma_inputform`, which reads the
  numeric `.wl` dialect.
- **Finer submodule split than originally sketched.** The original 5-way split
  (`Series·MST·RadialMode·SelfForce·Catalog`) is realized as 7 modules: `MMAInputForm` (the
  symbolic walker) split out of `Catalog` (pure discovery/browse, stdlib-only), and `SourceMode`
  (source + flux assembly) split out of `RadialMode` (radial solutions + amplitudes).
- **`Normal` is one generic.** `Series` owns it (on `PNSeries`); `SelfForce` extends the same
  method table on `ResummedSeriesData` — as both legacy packages already did, now within one
  namespace.
- **Owned data.** The 91-file SeriesData store is bundled inside the capability (real copy), so
  tests and loaders no longer reach into the earlier `PostNewtonianSelfForce` research code.

## Public API

```@autodocs
Modules = [PostNewtonian, PostNewtonian.Series, PostNewtonian.Catalog,
           PostNewtonian.MMAInputForm, PostNewtonian.SelfForce, PostNewtonian.MST,
           PostNewtonian.RadialMode, PostNewtonian.SourceMode]
Order = [:module, :type, :constant, :function]
```

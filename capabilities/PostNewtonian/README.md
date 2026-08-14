# PostNewtonian

Tier-4 GrTool capability — the **post-Newtonian series engine** and the **Symbolics quarantine**
(design §4.6). It fuses the two legacy PN leaf packages into one capability and is a top-of-DAG leaf:
nothing depends UP on it, so `Symbolics` (a HARD dep here and ONLY here) never leaks its ~1.9s
precompile floor onto the numeric/GPU waveform path.

## Physics owned

PN expansions of Teukolsky radial/source modes and PN self-force, carried as `Symbolics.Num`-valued
Laurent/Taylor series. The two quarantine anchors are struct field types — `PNSeries.coeffs ::
Vector{Symbolics.Num}` and `ResummedSeriesData.sf :: Symbolics.Num` — which is precisely why
Symbolics cannot be a weakdep and must be a hard dep confined to this capability.

## Submodules

| Submodule    | Owns | Legacy source |
|--------------|------|---------------|
| `Series`     | `PNSeries` Num-algebra: Laurent/Taylor arithmetic, `Normal` materialization, change-of-parameter, power-counting, `expand_log`/`expand_gamma` | TeukolskyPN/pnseries.jl |
| `MST`        | symbolic MST coefficients, ν-solve, ε-series, `imagunit` | TeukolskyPN/mst_coefficients.jl |
| `RadialMode` | η-graded symbolic radial + amplitude + source δ-coefficient series, mode assembly | TeukolskyPN/{radial_pn,mode_pn*,mode_pn_assemble}.jl |
| `SelfForce`  | `PostNewtonianData` accessor, `postnewtonian_coefficient`, `ResummedSeriesData`, `postnewtonian_expansion` | PostNewtonianSelfForce/{postnewtonian_data,coefficients,resummed_series}.jl |
| `Catalog`    | PN `.m` discovery + unified `read_mma_inputform` MMA-InputForm→`Num` walker | PostNewtonianSelfForce/{catalog_selector,series_data}.jl + TeukolskyPN/mode_pn_source.jl |

The duplicated `_eval_mma_ast` / `_eval_pn_ast` walker forks unify into `Catalog.read_mma_inputform`
behind a byte-for-byte `Num`-equivalence gate over all 91 legacy `.m` files.

## Absorbs (legacy)

- `TeukolskyPN`
- `PostNewtonianSelfForce`

These stay green as the oracle; this capability never edits them.

## Public API

Populated per-module as primitives land (`PNSeries`, `Normal`, `read_mma_inputform`,
`postnewtonian_coefficient`, `ResummedSeriesData`, `postnewtonian_expansion`, ...).

## Dependencies

- GrTool: `NumericsCore`, `Harmonics` (numeric `SpinWeightedSpheroidalEigenvalueSeries` →
  `Vector{Float64}`, carries no Symbolics), `Data` (`.m`/`.wl` artifact readers — **no HDF5**).
- external: `Symbolics` (HARD, quarantined), std `Printf`/`LinearAlgebra`.

## Grounding (TDD)

No heavy live legacy test dep (both PN packages are Symbolics-heavy). Tests ground against the
wolframscript byte-exact oracle on the 91 `.m` SeriesData files, captured fixtures, and analytic PN
identities (e.g. Schwarzschild ℓ=2,m=2 leading flux = Newtonian quadrupole 32/5).

# Data artifacts

*This page is the single home for the GrTool data-artifact story (design doc §10, per the §12
single-source rule): the schema contract, the HDF5/JLD2 extension-trigger rule, the injected data
root, the reader families, and the format separation between the on-disk dialects. The
[Data](../core/data.md) core page and the capability pages link here rather than restating it.*

GrTool packages consume large precomputed corpora — H1Lorenz retarded-field modes, FEW /
KerrBackwards Teukolsky amplitude and flux grids, the ZNAmps Kerr eccentric-equatorial B-spline
coefficients, ROMAN NN weights, PN series catalogs, WaSABI flux tables, ModifiedGravity source-mode
tables. All of it flows through the tier-1 `Data` package (`GrToolData`) under three binding
rules: a **versioned schema contract** decoupled from storage, an **extension-trigger rule** that
keeps HDF5/JLD2 out of the numeric path, and an **injected data root** that forbids path escapes.

## The schema contract: shapes in `NumericsCore`, bodies in `Data` extensions

The grid/inhomogeneous-field layouts are a **versioned contract**, not a side effect of whichever
library wrote the file. The contract structs are *dep-free* — defined in `NumericsCore`
(`core/NumericsCore/src/dataschema.jl`), which pulls in no I/O library at all — and re-exported by
`Data` so consumers reach them through the data spine:

| Struct | Contract | Layout (read literally from the legacy writers) |
|---|---|---|
| `H1Schema{T}` | one H1Lorenz `(l,m)` mode file | radial `grid`; split field tables `inhom_left` (`4nf × (r0index+1)`) and `inhom_right` (`4nf × (gs−r0index)`); asymptotic amplitudes `C_inf`/`C_horiz`; 0-based `r0index`/`ra_index`/`rb_index` |
| `AmpGridSchema{T}` | Teukolsky amplitude grid (FEW `ampinterp2d`) | strictly increasing axes `w` (eccentricity `e`) and `u` (`y(p,e)`), plus the stored `(l,m,n)` mode list; each mode a `length(w) × length(u)` surface |
| `FluxGridSchema{T}` | one Kerr flux-grid region (FEW `kerrflux`) | tricubic axes `u`,`w`,`z` on `[0,1]`, region tag `:A` (close) / `:B` (far) |

Each struct carries `schema_version(s)` (currently `SCHEMA_VERSION = 1`; bumped on any field-layout
change) and a field-invariant validator `validate_schema(s)` that throws
`ArgumentError`/`DimensionMismatch` on violation — e.g. the `H1Schema` column-count identities
`size(inhom_left, 2) == r0index + 1` and `size(inhom_right, 2) == gs − r0index`, or the ≥3
strictly-increasing nodes per interpolation axis (not-a-knot / tricubic minimum).

The I/O **bodies** that fill these shapes live in `Data`'s weakdep extensions
(`core/Data/ext/DataHDF5Ext.jl`, `DataJLD2Ext.jl`), never in `NumericsCore`. `read_h1` de-dups the
four legacy H1 readers (H1Lorenz, SecondOrderRicci, Punctures, DiscEMRIForce) over one `H1Schema`;
`read_ampgrid` is the single loader behind the byte-identical FEW-vs-KerrBackwards amplitude
schema. The HDF5-free parts of the H1 story — the parity-selected Barack–Lousto field lists
`h1_parity_fields(l, m)` and the component→storage-block map `hret_block(i, l)` (block `k`
occupies dataset rows `(k−1)·4 .+ (1:4)` = `[Re h, Im h, Re ∂ᵣh, Im ∂ᵣh]`) — are closed-form
tables in `Data`'s main source, always available.

## The extension-trigger rule (design §5.4, hygiene gate 5)

`Data` declares `HDF5` and `JLD2` only under `[weakdeps]` (see `core/Data/Project.toml`). A Julia
package extension activates **only when the *loading* package declares the trigger in its own
`[deps]`** — activation is *not* transitive through `Data`. In the main source each
HDF5/JLD2-backed function is a generic stub with an untyped fallback that raises an actionable
error ("load HDF5 (`import HDF5`) — Data's HDF5 reader-writer lives in a weakdep extension"); the
extension adds the strictly-more-specific `(AbstractString, …)` method, refining dispatch without a
method overwrite.

Who declares what (verified against the capability `Project.toml` files):

* **HDF5** — `SelfForce`, `EMRIWaveforms`, `AccretionDisc` (the grid/corpus consumers).
* **JLD2** — `Hybrid` (the WaSABI parse-once cache consumer).
* Everyone else (`Spacetime` → `Harmonics` → `Perturbation` → …, `PostNewtonian`) uses `Data`
  with no trigger, so HDF5/JLD2 never enter their dependency closure.

This is CI-enforced, not aspirational: `TestKit.assert_dep_budget` (design §7.7, gates 1+5;
`core/TestKit/src/hygiene.jl`) statically parses each package's own `Project.toml` and asserts the
trigger via its `must_declare` argument — `capabilities/EMRIWaveforms/test/hygiene.jl` runs
`assert_dep_budget(...; must_declare=["HDF5"], ...)`, `capabilities/Hybrid/test/hygiene.jl` uses
`must_declare=["JLD2"]`. The same gate keeps HDF5/JLD2 *out* of `Data`'s hard budget (they appear
only in its `weakdeps_allow`) and rejects any `[sources]` path into `packages/` or `ref-code/`.

## The injected data root: `set_datadir!` / `datadir`

The legacy loaders hard-coded paths that escaped their package directory (the DiscEMRIForce
`H1DIR` `run_output` path; the FEW `.fewdata` `@__DIR__`/`..` walks). `Data` replaces all of them
with one process-global, explicitly injected root (`core/Data/src/datadir.jl`):

```julia
using Data
set_datadir!("/path/to/corpora")   # or ENV["GRTOOLS_DATADIR"] = ...
path = datapath("h1-l2m2.h5")      # joinpath under the active root
```

* `set_datadir!(path)` — set (or, with `nothing`, clear) the override; returns the normalized
  absolute path.
* `datadir()` — the active root: the override if set, else `default_datadir()`, which reads
  `GRTOOLS_DATADIR` and otherwise **throws**. There is deliberately no implicit `ref-code`
  fallback: a missing configuration fails loudly instead of escaping a package.
* `datapath(parts...)` — `joinpath` under `datadir()`; the canonical way every reader resolves an
  on-disk asset.
* `available_modes(; dir=datadir())` — discover the `(l, m)` h1 mode files (`h1-l{l}m{m}.h5`).

For named corpora, `artifact_path(name)` (`core/Data/src/artifacts.jl`) resolves through an
explicit **override → artifact → error** chain: the injected datadir wins first (HPC/offline and
test redirect), else a declared *and provisioned* `Artifacts.toml` entry resolves lazily via the
`Artifacts` stdlib tree-hash, else `datapath` throws loudly. `stage_datadir!(src, dst)` pre-stages
a corpus into a local root and activates it (for sites without artifact download — e.g. the
4.74 GB ZNAmps grid staged from shared scratch), and `ensure_local(name)` asserts presence with
the staging hint.

## Reader families

| Family (source file) | Entry points | Backend | Corpus |
|---|---|---|---|
| `h1io` | `read_h1`, `write_mode_h5`; `h1_parity_fields`, `hret_block` | HDF5 ext | H1Lorenz `(l,m)` mode files `h1-l{l}m{m}.h5` |
| `ampgrid` | `read_ampgrid`, `schwarzecc_amp_modes` | HDF5 ext | Schwarzschild-eccentric Teukolsky amplitude grid (FEW = KerrBackwards) |
| `fluxgrid` | `read_fluxgrid(path, region)` | HDF5 ext | Kerr flux grid, region-streamed `:A`/`:B` |
| `kerrampgrid` | `read_kerrampgrid`, `KerrAmpGrid`, `kerr_mode_table` | HDF5 ext | ZNAmps Kerr-ecc-eq B-spline coefficients (`ZNAmps_l10_m10_n55_DS2Outer.h5`) |
| `romannn` | `read_roman_nn` | HDF5 ext | ROMAN dense-NN weights (`SchwarzschildEccentricInput.hdf5`) |
| `dlm` | `read_dat`, `read_csv`, `SourceMode`, `read_source_table` | `DelimitedFiles` (always on) | `.dat`/`.csv` numeric tables, incl. the 12-column source-mode table |
| `wl` | `read_mma_inputform`, `read_wl_table`, `read_wl_data`, `read_wl_grid2d`, `read_wl_grid3d`, `read_wl_table_complex` | pure string parsing (always on) | Mathematica InputForm `.wl`/`.m` numeric assets (WaSABI) |
| `jld2io` | `wl_to_jld2`, `read_jld2` | JLD2 ext | parse-once-load-later cache of parsed `.data` catalogs |

Container structs are dep-free even for the extension-backed families: `KerrAmpGrid` (coefficient
slabs + shared spin knots `z_knots`) lives in `Data`'s main source; only the `h5open` body is in
the extension. Evaluators that *consume* the loaded slabs stay with the physics — the clamped
B-spline evaluator in `EMRIWaveforms.KerrAmpSpline`, the ROMAN forward pass in
`Kernels.roman_forward` — because they are numerics, not I/O.

## Format separation (design §10)

Four on-disk dialects, four deliberately separate readers — same-looking files are **not** merged:

* **PN `.m` (symbolic)** — the PN catalog `.m` files carry Mathematica
  `SeriesData`/`Association`/`Module` expressions that must become `Symbolics.Num`. Their reader
  is `PostNewtonian.MMAInputForm` (`parse_pn_data_file`, used by `PostNewtonian.SelfForce` to
  build `PostNewtonianData`), behind the Symbolics quarantine. It is a *sibling* of the `Data`
  reader, not part of it: `Data` must stay Symbolics-free on the numeric path, so there is
  deliberately no symbolic walker in `Data`'s `wl.jl`.
* **WaSABI `.wl`/`.m`/`.data` (numeric InputForm)** — `Data.read_mma_inputform` and the
  `read_wl_*` wrappers: a self-contained recursive-descent parser (no `eval`, no world-age) over
  nested `{}` lists, `<|…|>` associations, `Around[c,e]` (reduced to the central value),
  `` `precision `` markers and `*^` exponents; `want_big=true` keeps full file precision as
  `BigFloat`. `Hybrid`'s inspiral/amplitude models read the 1SF/2SF flux and redshift tables this
  way (e.g. `read_wl_table`/`read_wl_table_complex` in `Hybrid/src/Models/inspiral_1PAT1.jl`), and
  `wl_to_jld2` caches a `.data` catalog so it is parsed once and thereafter loaded bit-for-bit.
* **HDF5 corpora** — everything schema-shaped and large (h1 modes, amplitude/flux grids, ZNAmps
  coefficients, ROMAN weights) goes through the `DataHDF5Ext`-backed readers above, gated by the
  extension-trigger rule.
* **12-column `.dat` source tables** — the ModifiedGravity/ModifiedTeukolsky tab-separated
  source-mode table (`Sourcefiletotable`/`mywrite`, Polcar & Witzany 2025, arXiv:2507.15720):
  integers `m,n,k`, mode frequency `ω`, then the four projected source amplitudes
  `Sl, Sn, Sm, Sm†` as interleaved Re/Im pairs. `Data.read_source_table` returns
  `Vector{SourceMode}`; `radius`/`ringmass` are caller-supplied (not stored in the file).
  `ModifiedGravity` (`src/source_tables.jl`) consumes it.

The naming collision is intentional and disambiguated by module: `Data.read_mma_inputform` parses
the *numeric* dialect; `PostNewtonian.MMAInputForm.read_mma_inputform` walks the *symbolic* one.
Pick by what the file contains, not by its extension — WaSABI ships numeric `.m` files too.

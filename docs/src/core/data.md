# Data (`GrToolsData`)

Tier-1 core **data spine** of GrTools. Every on-disk corpus the physics packages consume — H1Lorenz
retarded-field modes, FEW/KerrBackwards Teukolsky amplitude and flux grids, the ZNAmps Kerr
eccentric-equatorial B-spline coefficients, the ROMAN NN weights, PN `.dat` series and Mathematica
`.wl`/`.m` catalogs — is read through this one package. The schema *shapes* (`H1Schema`,
`AmpGridSchema`, `FluxGridSchema`) are dep-free contract structs defined in `NumericsCore.dataschema`
and re-exported here; `Data` owns the I/O **bodies** that fill them, plus a single injected
data-root policy (`set_datadir!` / `datadir` / `artifact_path`) that replaces the legacy hard-coded
`ref-code/.../run_output` path escapes.

The design rule: storage backends are *decoupled* from schema shapes. HDF5 and
JLD2 never appear in `Data`'s hard `[deps]` — their reader-writer bodies live in weakdep
**extensions**, so the numeric path (Spacetime → Harmonics → Perturbation → …) never pays for HDF5,
and HDF5 stays confined to the tier-5/6 grid consumers.

## Module map

`Data` is a single flat module (no submodules); the source files are functional areas. Two package
extensions carry the backend bodies.

| Source | Responsibility |
|---|---|
| `src/datadir.jl` | process-global data root: `set_datadir!`, `default_datadir` (env `GRTOOLS_DATADIR`), `datadir`, `datapath`, `available_modes` (h1 `(l,m)` mode-file discovery) |
| `src/artifacts.jl` | `Artifacts.toml` resolution: `artifact_path`, `ensure_local`, `stage_datadir!` (offline/HPC pre-staging), `set_artifacts_toml!`, `artifacts_toml` |
| `src/h1io.jl` | H1Lorenz `(l,m)` mode read/write over the versioned `H1Schema`: `read_h1`, `write_mode_h5`, plus the HDF5-free closed-form index tables `h1_parity_fields`, `hret_block` |
| `src/ampgrid.jl` | Teukolsky amplitude-grid loader (`AmpGridSchema`): `read_ampgrid`, `schwarzecc_amp_modes` |
| `src/fluxgrid.jl` | Kerr flux-grid loader, region-streamed (`:A` close / `:B` far) over `FluxGridSchema`: `read_fluxgrid` |
| `src/kerrampgrid.jl` | ZNAmps Kerr-ecc-eq B-spline coefficient grid: `KerrAmpGrid` container, `kerr_mode_table`, `read_kerrampgrid` |
| `src/romannn.jl` | ROMAN dense-NN weight loader (`SchwarzschildEccentricInput.hdf5`): `read_roman_nn` |
| `src/dlm.jl` | `.dat`/`.csv` numeric tables via `DelimitedFiles`: `read_dat`, `read_csv`, and the 12-column `SourceMode` table reader `read_source_table` |
| `src/wl.jl` | Mathematica `.wl`/`.m` InputForm readers (BigFloat-exact): `read_mma_inputform`, `read_wl_table`, `read_wl_data`, `read_wl_grid2d`, `read_wl_grid3d`, `read_wl_table_complex` |
| `src/jld2io.jl` | JLD2 parse-once cache path: `wl_to_jld2`, `read_jld2` |
| `ext/DataHDF5Ext.jl` | HDF5 bodies for `read_h1`, `write_mode_h5`, `read_ampgrid`, `read_fluxgrid`, `read_kerrampgrid`, `read_roman_nn` |
| `ext/DataJLD2Ext.jl` | JLD2 bodies for `wl_to_jld2`, `read_jld2` |

### Extension-activation rule (binding)

A Julia extension activates **only when the loading package declares the trigger in its own
`[deps]`**. The HDF5/JLD2-backed functions above are declared in `Data`'s main source as generic
stubs whose untyped fallback raises an actionable "load the trigger" error; the extension supplies
the strictly-more-specific `(AbstractString, …)` method once the trigger is present. Consequences:

- Grid consumers — `SelfForce`, `EMRIWaveforms`, `AccretionDisc` (Force) for HDF5; `Hybrid` for
  JLD2 — must declare `HDF5`/`JLD2` **directly** in their own `[deps]`. It is *not* transitive
  through `Data`.
- `PostNewtonian` and the `.wl`/`.dat` paths use `Data` without any trigger, so HDF5 never enters
  their dependency closure.

## Dependencies and poison-axis declaration

Hard `[deps]` (from `core/Data/Project.toml`): `NumericsCore` (schema structs, `p_to_y`,
`mode_index_table`), `DelimitedFiles`, `Artifacts`, `Printf` (all stdlib except `NumericsCore`).
Weakdeps: `HDF5`, `JLD2` (extension triggers only).

The dependency budget is CI-enforced by `core/Data/test/hygiene.jl`:

- **allow** (closed hard-dep budget): `NumericsCore`, `DelimitedFiles`, `Artifacts`, `Printf`.
- **weakdeps\_allow** (the sanctioned weak axes, pinned exactly): `HDF5`, `JLD2`.
- **forbid** (poison axes that must never appear): `Symbolics`, `OrdinaryDiffEq*`, `CUDA`,
  `AMDGPU`, `KernelAbstractions`, `WignerSymbols`, `Kernels`, `Roots`, `NLsolve`, `SciMLBase`,
  `SpecialFunctions`, `HypergeometricFunctions`, `QuadGK`.
- A closure-load gate additionally asserts that `using Data` does not pull `Symbolics`,
  `OrdinaryDiffEq*`, `CUDA`, `AMDGPU`, `KernelAbstractions`, or `WignerSymbols` into the session.

## GPU / data needs

`Data` itself has **no GPU path** (GPU kernels live in `Kernels`); it is pure I/O. What it needs is
the on-disk corpora, resolved through the injected root — never through hard-coded paths:

- H1Lorenz mode files `h1-l{l}m{m}.h5` (discovered by `available_modes`).
- FEW/KerrBackwards Teukolsky amplitude grids (`read_ampgrid`) and the 4.74 GB Kerr-ecc-eq flux
  grid `KerrEccEqFluxData.h5`, which is deliberately **region-streamed** (`read_fluxgrid(path, :A)`
  / `:B`, with a `downsample` stride for grid-convergence tests) rather than loaded whole.
- `ZNAmps_l10_m10_n55_DS2Outer.h5` B-spline coefficient slabs (`read_kerrampgrid`, per-`(l,m,n)`
  subset selection).
- `SchwarzschildEccentricInput.hdf5` ROMAN network weights (`read_roman_nn`).
- PN `.dat`/`.csv` tables and Mathematica `.wl`/`.m` catalogs (BigFloat-exact InputForm parse;
  `wl_to_jld2` caches the parse so a catalog is parsed once, then loaded).

Root resolution order: `set_datadir!` override → `GRTOOLS_DATADIR` env var → pinned artifact →
loud error. There is deliberately **no** implicit `ref-code` fallback. `stage_datadir!` supports
offline/HPC pre-staging.

## Provenance

`Data` is a de-duplication point: each reader below collapses several byte-compatible legacy
loaders into one canonical body.

| Data function | Replaces (legacy package · function) |
|---|---|
| `read_h1` / `write_mode_h5` | H1Lorenz `io_hdf5` (`read_grid`/`write_mode_h5`); SecondOrderRicci `read_h1`; Punctures `read_hret_mode`; DiscEMRIForce `load_h1_mode` |
| `read_ampgrid` | FastEMRIWaveforms `ampinterp2d` ≡ KerrBackwards amplitude loader (byte-identical schema strings) |
| `read_fluxgrid` | FastEMRIWaveforms `kerrflux` ≡ KerrBackwards flux loader |
| `read_kerrampgrid` / `kerr_mode_table` | `EMRIWaveforms.KerrAmpSpline` grid-I/O half (itself from FEW `KerrEccentricEquatorial`) |
| `read_roman_nn` | `EMRIWaveforms.AmplitudeModels._load_roman_nn` (FEW ROMAN) |
| `read_dat` / `read_source_table` | PN series tables; ModifiedGravity 12-column source table (Polcar & Witzany 2025) |
| `read_mma_inputform` / `read_wl_*` | WaSABI `wl_parser`; TeukolskyPN/PN InputForm walkers |
| `wl_to_jld2` / `read_jld2` | WaSABI `.data` parse-once cache |
| `set_datadir!` / `datadir` / `available_modes` | DiscEMRIForce `available_modes` + the assorted `H1DIR`/`run_output` path constants across the legacy corpus |

## What changed from the predecessor codes

- **One schema, one reader.** Four independent h1 readers (H1Lorenz / SecondOrderRicci /
  Punctures / DiscEMRIForce) collapse onto the single versioned `H1Schema` contract;
  `write_mode_h5` now takes an `H1Schema` (not the legacy `NModeData`/`OrbitalParams`), so the
  schema *is* the storage contract. `read_h1(path; T)` widens the Float64 corpus to the requested
  working precision.
- **HDF5/JLD2 demoted to weakdeps.** The predecessor codes carried HDF5 as a hard dep everywhere; here
  the bodies moved into `DataHDF5Ext`/`DataJLD2Ext`, and only tier-5/6 consumers trigger them.
- **Reader relocation:** the ZNAmps Kerr-amp grid I/O and the ROMAN
  NN-weight loader moved out of `EMRIWaveforms` into `Data` — the container structs and the
  `(l,m,n)` mode table live in main `src/`, the `h5open` bodies in `ext/`. The *numerics* stayed
  behind: the clamped-B-spline evaluator remains in `EMRIWaveforms.KerrAmpSpline`, the GEMM
  forward pass in `Kernels.roman_forward`. `kerr_mode_table` is a thin adapter over
  `NumericsCore.mode_index_table` (the ZNAmps axis is its un-mirrored ``m \ge 0`` prefix).
- **Path injection.** All readers resolve assets through `datapath` under the active `datadir`;
  the legacy hard-coded `ref-code/.../run_output` escapes are gone, and a missing configuration
  fails loudly instead of silently reading a stale tree.
- **Normalized quirks.** The h1 corpus stores grid-index attributes either as scalars (H1Lorenz
  writer) or as length-1 `Int32` vectors (the C self-force code) — `read_h1` normalizes both;
  `read_fluxgrid` undoes HDF5's on-disk dimension reversal with an explicit `permutedims`.

## Public API

The docstrings below are the trigger-free surface (main module). The HDF5/JLD2 method *bodies* and
their extended docstrings live in the `DataHDF5Ext`/`DataJLD2Ext` extension modules, which are not
loaded in this docs build — calling those functions without the trigger raises an error telling you
which package to load.

```@autodocs
Modules = [Data]
Order = [:module, :type, :constant, :function]
```

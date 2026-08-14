# Data (`GrToolData`)

Tier-1 core **data spine** of GrTool. Owns the file-format I/O *bodies* that sit behind the
dep-free `DataSchema` contract structs (`H1Schema` / `AmpGridSchema` / `FluxGridSchema`, defined in
`NumericsCore.dataschema`). The schema *shapes* are decoupled from the storage backend; this package
fills them with HDF5 / JLD2 / DelimitedFiles / Mathematica readers and a single injected data-root
policy.

## What it owns

| Module | Responsibility | De-dups legacy |
|---|---|---|
| `h1io`     | H1Lorenz `(l,m)` mode read/write over the versioned `H1Schema` | H1Lorenz `io_hdf5`, SecondOrderRicci `read_h1`, Punctures `read_hret_mode`, DiscEMRIForce `load_h1_mode` |
| `ampgrid`  | Teukolsky amplitude-grid loader (`AmpGridSchema`) | FEW `ampinterp2d` ⟂ KerrBackwards (byte-identical schema) |
| `fluxgrid` | Kerr flux-grid loader, region-streamed (`FluxGridSchema`) | FEW `kerrflux` ⟂ KerrBackwards |
| `dlm`      | `.dat`/`.csv` numeric tables (`DelimitedFiles`) | PN series, ModifiedGravity 12-col source |
| `wl`       | Mathematica `.wl`/`.m` InputForm readers | WaSABI `wl_parser`, TeukolskyPN/PN InputForm walkers (→ `read_mma_inputform`) |
| `datadir` / `artifacts` | injected `datadir=` + `set_datadir!` + `Artifacts` resolution | kills the legacy `ref-code/.../run_output` path escapes |

## HDF5 / JLD2 are weakdep extensions (design §5.4)

The HDF5 and JLD2 reader-writer bodies live in `ext/DataHDF5Ext.jl` / `ext/DataJLD2Ext.jl`. A Julia
extension activates **only when the loading package declares the trigger in its own `[deps]`**.
So the grid/corpus consumers — `SelfForce`, `EMRIWaveforms`, `AccretionDisc.Force` (HDF5) and
`Hybrid` (JLD2) — declare the trigger directly; the numeric path never pays for HDF5/JLD2.
`DelimitedFiles` and `Artifacts` are normal `[deps]` (always available).

## Public API

`H1Schema`, `AmpGridSchema`, `FluxGridSchema`, `schema_version`, `validate_schema`, `SCHEMA_VERSION`
(re-exported from `NumericsCore`); `read_h1`, `write_mode_h5`, `read_ampgrid`, `read_fluxgrid`,
`read_dat`, `read_mma_inputform`, `read_wl_table`, `wl_to_jld2`; `set_datadir!`, `datadir`,
`artifact_path`.

## Deps

`NumericsCore` (schema structs). External: `DelimitedFiles`, `Artifacts` (hard); `HDF5`, `JLD2`
(weakdep extension triggers).

## Status

Scaffold: module loads, schema contract re-exported, weakdep-ext seams in place. Reader/writer bodies
are ported tests-first (round-trip identity + checksum-verified FEW grids) in the Data build round.

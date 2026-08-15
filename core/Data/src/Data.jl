"""
    Data (GrToolsData)

Tier-1 core data spine for GrTools. Owns the I/O BODIES behind the dep-free
`DataSchema` contract structs (`H1Schema`/`AmpGridSchema`/`FluxGridSchema`, defined in
`NumericsCore.dataschema`):

  * `h1io`     — versioned H1Lorenz `(l,m)` mode reader/writer (de-dups the H1Lorenz /
                 SecondOrderRicci / Punctures / DiscEMRIForce readers over one `H1Schema`).
  * `ampgrid`  — Teukolsky amplitude-grid loader (the byte-identical FEW-vs-KerrBackwards schema).
  * `fluxgrid` — Kerr flux-grid loader (per-region tricubic axes).
  * `kerrampgrid` — ZNAmps Kerr-ecc-eq B-spline-coefficient grid loader (`KerrAmpGrid` +
                 `kerr_mode_table`; moved from `EMRIWaveforms.KerrAmpSpline`, S4-b §7.2/§10).
  * `romannn`  — ROMAN NN-weight loader (`SchwarzschildEccentricInput.hdf5`; moved from
                 `EMRIWaveforms.AmplitudeModels._load_roman_nn`, S4-b §7.2/§10).
  * `dlm`      — `.dat`/`.csv` readers (PN series, ModifiedGravity 12-col tables) via `DelimitedFiles`.
  * `wl`       — Mathematica `.wl`/`.m` InputForm reader (`read_mma_inputform`).
  * `artifacts`/`datadir` — `Artifacts` resolution + an injected `datadir=` / `set_datadir!`
                 that kills the legacy `ref-code/.../run_output` path escapes.

Extension-activation rule (design §5.4): the HDF5/JLD2 reader-writer bodies live in
`ext/DataHDF5Ext.jl` / `ext/DataJLD2Ext.jl` weakdep extensions. They activate ONLY when the
LOADING package declares `HDF5`/`JLD2` in its own `[deps]`, so the numeric path stays HDF5-free.
The generic functions below are declared here and given methods in those extensions; called
without the trigger loaded they raise an actionable error.
"""
module Data

using NumericsCore: H1Schema, AmpGridSchema, FluxGridSchema, schema_version, validate_schema,
    SCHEMA_VERSION, p_to_y, mode_index_table
using DelimitedFiles
using Artifacts

# Re-export the contract structs so consumers reach them through the data spine.
export H1Schema, AmpGridSchema, FluxGridSchema, schema_version, validate_schema, SCHEMA_VERSION

# datadir / artifacts (no heavy deps — always available)
include("datadir.jl")
include("artifacts.jl")

# Generic I/O entry points. Bodies for HDF5/JLD2-backed formats are added by the weakdep
# extensions; the stubs here error with a load-the-trigger hint until then.
include("h1io.jl")
include("ampgrid.jl")
include("fluxgrid.jl")
include("kerrampgrid.jl")
include("romannn.jl")

# DelimitedFiles / hand-rolled readers — always available (no weakdep trigger needed).
include("dlm.jl")
include("wl.jl")

# JLD2 cache path (parse-once-load-later). Bodies in DataJLD2Ext; stubs error until JLD2 is loaded.
include("jld2io.jl")

end # module Data

# HDF5 reader-writer BODIES for the Data spine (design §5.4 weakdep extension). Activates only when
# the loading package declares HDF5 in its own [deps], so the numeric path never pays for HDF5.
# Provides the bodies for `read_h1`/`write_mode_h5` (H1Schema), `read_ampgrid` (AmpGridSchema),
# `read_fluxgrid` (FluxGridSchema), `read_kerrampgrid` (KerrAmpGrid) and `read_roman_nn`;
# these override the load-the-trigger error stubs in Data's main src once HDF5 is present.
module DataHDF5Ext

using Data
using Data: H1Schema, AmpGridSchema, FluxGridSchema, validate_schema, KerrAmpGrid, kerr_mode_table
using HDF5

# The h1 corpus stores the `grid` index attributes either as scalars (H1Lorenz writer) or as
# length-1 Int32 vectors (the C self-force code). Normalize to a plain Int.
_attr1(x) = x isa AbstractArray ? Int(x[1]) : Int(x)

"""
    Data.read_h1(path; T=Float64) -> H1Schema{T}

Read one H1Lorenz `(l,m)` mode HDF5 file into an `H1Schema`. The single canonical reader collapsing
the H1Lorenz/SecondOrderRicci/Punctures/DiscEMRIForce loaders: `grid` + the `r0index`/`raindex`/
`rbindex` (0-based) attrs, the split `inhom_left` (`4nf × (r0index+1)`) / `inhom_right`
(`4nf × (gs−r0index)`) tables, and the `C_inf`/`C_horiz` amplitude vectors. `T` widens the loaded
Float64 corpus to the requested working precision (the field values stay Float64-accurate).
"""
function Data.read_h1(path::AbstractString; T::Type = Float64)
    h5open(path, "r") do f
        g = f["grid"]
        grid = T.(read(g))
        r0index = _attr1(read_attribute(g, "r0index"))
        ra_index = haskey(HDF5.attributes(g), "raindex") ? _attr1(read_attribute(g, "raindex")) : 0
        rb_index = haskey(HDF5.attributes(g), "rbindex") ? _attr1(read_attribute(g, "rbindex")) :
                   length(grid) - 1
        inhom_left = T.(read(f["inhom_left"]))
        inhom_right = T.(read(f["inhom_right"]))
        C_inf = haskey(f, "C_inf") ? T.(read(f["C_inf"])) : T[]
        C_horiz = haskey(f, "C_horiz") ? T.(read(f["C_horiz"])) : T[]
        return H1Schema{T}(grid, inhom_left, inhom_right, C_inf, C_horiz,
            r0index, ra_index, rb_index)
    end
end

"""
    Data.write_mode_h5(path, s::H1Schema) -> path

Write an `H1Schema` to an HDF5 mode file in the versioned oracle layout: `grid` (Float64) with
`raindex`/`r0index`/`rbindex` attrs, `inhom_left`/`inhom_right`, `C_inf`/`C_horiz`. Rewritten to take
an `H1Schema` (not the legacy `NModeData`/`OrbitalParams`), so the schema is the storage contract.
"""
function Data.write_mode_h5(path::AbstractString, s::H1Schema)
    validate_schema(s)
    h5open(path, "w") do f
        f["grid"] = Float64.(s.grid)
        HDF5.attributes(f["grid"])["raindex"] = Int32(s.ra_index)
        HDF5.attributes(f["grid"])["r0index"] = Int32(s.r0index)
        HDF5.attributes(f["grid"])["rbindex"] = Int32(s.rb_index)
        f["inhom_left"] = Float64.(s.inhom_left)
        f["inhom_right"] = Float64.(s.inhom_right)
        f["C_inf"] = Float64.(s.C_inf)
        f["C_horiz"] = Float64.(s.C_horiz)
    end
    return path
end

include("hdf5_ampgrid.jl")
include("hdf5_kerrampgrid.jl")
include("hdf5_romannn.jl")

# Per-region grid metadata: HDF5 group name + the (NU, NW, NZ) tricubic axis counts. Read verbatim
# from FastEMRIWaveforms/src/trajectory/kerrflux.jl (_region_grid dims + KerrEccEqFlux NU?,NW?,NZ?).
const _FLUX_REGION = Dict(
    :A => ("regionA", 129, 65, 65),
    :B => ("regionB", 65, 33, 33),
)

"""
    Data.read_fluxgrid(path, region; downsample=(1,1,1)) -> (FluxGridSchema{Float64}, channels)

Read one Kerr eccentric-equatorial flux-grid region (`:A` close / `:B` far) from the FEW
`KerrEccEqFluxData.h5`. Returns the tricubic `(u,w,z) ∈ [0,1]³` axes as a `FluxGridSchema{Float64}`
(`u,w,z = range(0,1; length=NU/NW/NZ)`) plus `Dict(:Edot, :Ldot)` of the raw `NU×NW×NZ` flux samples.
The single canonical reader de-duping the FEW-vs-KerrBackwards flux loaders: HDF5 reverses the on-disk
`(NU,NW,NZ)` dims, so each channel is `permutedims(read, (3,2,1))`; `downsample=(su,sw,sz)` strides the
axes and channels in lockstep (`[1:su:end, 1:sw:end, 1:sz:end]`) for grid-convergence/testing.
"""
function Data.read_fluxgrid(path::AbstractString, region::Symbol;
                            downsample::NTuple{3,Int} = (1, 1, 1))
    haskey(_FLUX_REGION, region) ||
        throw(ArgumentError("read_fluxgrid: region must be :A or :B, got $region"))
    name, NU, NW, NZ = _FLUX_REGION[region]
    su, sw, sz = downsample
    u = collect(range(0.0, 1.0; length = NU))[1:su:end]
    w = collect(range(0.0, 1.0; length = NW))[1:sw:end]
    z = collect(range(0.0, 1.0; length = NZ))[1:sz:end]
    Eraw, Lraw = h5open(path, "r") do f
        (read(f["$name/Edot"]), read(f["$name/Ldot"]))
    end
    # HDF5 reverses dims (NZ,NW,NU) -> permute back to (NU,NW,NZ), then stride.
    Edot = permutedims(Eraw, (3, 2, 1))[1:su:end, 1:sw:end, 1:sz:end]
    Ldot = permutedims(Lraw, (3, 2, 1))[1:su:end, 1:sw:end, 1:sz:end]
    schema = FluxGridSchema{Float64}(region, u, w, z)
    return schema, Dict(:Edot => Edot, :Ldot => Ldot)
end

end # module DataHDF5Ext

# Dep-free DataSchema contract structs — the infra seam that `GrToolsData` fills with HDF5/JLD2
# bodies. NO HDF5 here: these are pure struct shapes + field-invariant validators, so the grid/inhom
# layout is a VERSIONED CONTRACT decoupled from the storage backend.
#
# Field layouts are read literally from the legacy writers (which are NOT live-dep'd — H1Lorenz pulls
# HDF5, FEW/KerrBackwards are heavy):
#   - H1Schema  ← packages/H1Lorenz/src/io_hdf5.jl:28-58 (write_mode_h5) + assemble.jl:6-15 (NModeData):
#       inhom_left  is (4nf, r0index+1)   — grid columns 0..r0index   (nL = r0i+1)
#       inhom_right is (4nf, gs-r0index)  — right-limit at r0 + grid columns r0index+1..gs-1
#       grid attrs raindex/r0index/rbindex; C_inf/C_horiz amplitude vectors.
#   - AmpGridSchema  ← packages/FastEMRIWaveforms/src/amplitude/ampinterp2d.jl:84-121: (w=e, u=y) axes
#       + the stored (l,m,n) mode list; each mode a (num_w × num_u) surface.
#   - FluxGridSchema ← packages/FastEMRIWaveforms/src/trajectory/kerrflux.jl:99-114: per-region
#       (u,w,z) ∈ [0,1]³ tricubic axes (region :A close / :B far).

"Current DataSchema contract version. Bump on any field-layout change."
const SCHEMA_VERSION = 1

"""
    H1Schema{T}

Versioned contract for one H1Lorenz `(l,m)` mode file: radial `grid`, the split inhomogeneous-field
tables `inhom_left` (`4nf × (r0index+1)`) and `inhom_right` (`4nf × (gs−r0index)`), asymptotic
amplitude vectors `C_inf`/`C_horiz`, and the important grid indices `r0index`/`ra_index`/`rb_index`
(0-based, as the H1 loader stores them). Struct shape only — the HDF5 body lives in `GrToolsData`.
"""
struct H1Schema{T}
    grid::Vector{T}
    inhom_left::Matrix{T}
    inhom_right::Matrix{T}
    C_inf::Vector{T}
    C_horiz::Vector{T}
    r0index::Int
    ra_index::Int
    rb_index::Int
end

"""
    AmpGridSchema{T}

Versioned contract for a Teukolsky amplitude grid (FEW `ampinterp2d`): the interpolation axes
`w` (eccentricity `e`) and `u` (`y(p,e)`), strictly increasing, and the stored `(l,m,n)` mode list.
Each mode's surface is a `length(w) × length(u)` matrix (the body lives in `GrToolsData`).
"""
struct AmpGridSchema{T}
    w::Vector{T}
    u::Vector{T}
    modes::Vector{NTuple{3,Int}}
end

"""
    FluxGridSchema{T}

Versioned contract for one Kerr flux-grid region (FEW `kerrflux`): the tricubic axes `u`,`w`,`z`
(each strictly increasing on `[0,1]`) and the region tag (`:A` close / `:B` far). Each flux channel's
sample array is `length(u) × length(w) × length(z)` (the body lives in `GrToolsData`).
"""
struct FluxGridSchema{T}
    region::Symbol
    u::Vector{T}
    w::Vector{T}
    z::Vector{T}
end

"""`schema_version(s)` — the contract version of a DataSchema struct (currently `SCHEMA_VERSION`)."""
schema_version(::Union{H1Schema,AmpGridSchema,FluxGridSchema}) = SCHEMA_VERSION

_strictly_increasing(v) = length(v) >= 2 && all(diff(v) .> 0)

"""
    validate_schema(s) -> true

Field-invariant validator (throws `ArgumentError`/`DimensionMismatch` on violation). Checks the
literal H1 layout (`r0index` in range, `inhom_left`/`inhom_right` grid-axis lengths = `r0index+1` /
`gs−r0index`, matching row counts, equal-length amplitude vectors) and the interpolation-axis
invariants for the amplitude/flux grids (≥3 strictly-increasing nodes per axis — not-a-knot needs 3).
"""
function validate_schema(s::H1Schema)
    gs = length(s.grid)
    gs >= 1 || throw(ArgumentError("H1Schema: empty grid"))
    0 <= s.r0index < gs || throw(ArgumentError("H1Schema: r0index=$(s.r0index) not in 0:$(gs-1)"))
    0 <= s.ra_index < gs || throw(ArgumentError("H1Schema: ra_index=$(s.ra_index) not in 0:$(gs-1)"))
    0 <= s.rb_index < gs || throw(ArgumentError("H1Schema: rb_index=$(s.rb_index) not in 0:$(gs-1)"))
    size(s.inhom_left, 2) == s.r0index + 1 ||
        throw(DimensionMismatch("H1Schema: inhom_left has $(size(s.inhom_left,2)) cols, expected r0index+1=$(s.r0index+1)"))
    size(s.inhom_right, 2) == gs - s.r0index ||
        throw(DimensionMismatch("H1Schema: inhom_right has $(size(s.inhom_right,2)) cols, expected gs−r0index=$(gs - s.r0index)"))
    size(s.inhom_left, 1) == size(s.inhom_right, 1) ||
        throw(DimensionMismatch("H1Schema: inhom_left/right row counts differ ($(size(s.inhom_left,1)) vs $(size(s.inhom_right,1)))"))
    size(s.inhom_left, 1) % 4 == 0 ||
        throw(DimensionMismatch("H1Schema: inhom row count $(size(s.inhom_left,1)) must be 4·nf"))
    length(s.C_inf) == length(s.C_horiz) ||
        throw(DimensionMismatch("H1Schema: C_inf/C_horiz length mismatch"))
    return true
end

function validate_schema(s::AmpGridSchema)
    length(s.w) >= 3 || throw(ArgumentError("AmpGridSchema: w axis needs ≥3 nodes (not-a-knot)"))
    length(s.u) >= 3 || throw(ArgumentError("AmpGridSchema: u axis needs ≥3 nodes (not-a-knot)"))
    _strictly_increasing(s.w) || throw(ArgumentError("AmpGridSchema: w axis must be strictly increasing"))
    _strictly_increasing(s.u) || throw(ArgumentError("AmpGridSchema: u axis must be strictly increasing"))
    isempty(s.modes) && throw(ArgumentError("AmpGridSchema: empty mode list"))
    return true
end

function validate_schema(s::FluxGridSchema)
    s.region in (:A, :B) || throw(ArgumentError("FluxGridSchema: region must be :A or :B, got $(s.region)"))
    for (nm, ax) in ((:u, s.u), (:w, s.w), (:z, s.z))
        length(ax) >= 3 || throw(ArgumentError("FluxGridSchema: $nm axis needs ≥3 nodes (tricubic)"))
        _strictly_increasing(ax) || throw(ArgumentError("FluxGridSchema: $nm axis must be strictly increasing"))
    end
    return true
end

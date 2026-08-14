# ZNAmps Kerr eccentric-equatorial amplitude-grid loader (S4-b, design §7.2/§10: ALL grid reader
# BODIES live in Data's HDF5 weakdep-ext). Moved from EMRIWaveforms/src/KerrAmpSpline.jl — the
# dep-free container struct + the (l,m,n) grid-index mode table live here; the `h5open` body is in
# ext/hdf5_kerrampgrid.jl. The clamped-B-spline EVALUATOR that consumes these slabs stays in
# EMRIWaveforms.AmplitudeModels.KerrAmpSpline (it is numerics, not I/O).

"""
    KerrAmpGrid{T}

Raw coefficient slabs parsed from `ZNAmps_l10_m10_n55_DS2Outer.h5` for a requested `(l,m,n)` subset.
`coeffsA`/`coeffsB` are `(n_coefflat, 2, n_modes, n_z)` (ReIm axis = 2; each flat coeff per `(mode,z)`
reshapes C-order to `(n_knots−4)²` — the reshape/evaluator lives in `EMRIWaveforms.KerrAmpSpline`).
`z_knots` are the spin-interpolation knots (shared A/B). The HDF5 body lives in `DataHDF5Ext`.
"""
struct KerrAmpGrid{T}
    modes::Vector{NTuple{3,Int}}
    z_knots::Vector{T}
    twA::Vector{T}; tuA::Vector{T}
    twB::Vector{T}; tuB::Vector{T}
    coeffsA::Array{T,4}
    coeffsB::Array{T,4}
    lmax::Int
    mmax::Int
    nmax::Int
    num_total_modes::Int
end

"""
    kerr_mode_table(; lmax=10, nmax=55) -> (modes, index)

The `(l,m,n)` mode list in ZNAmps grid-axis order (FEW `KerrEccentricEquatorial`: `m=0` block first
then `m>0`, each in `l`-outer/`n`-inner order) and the inverse `Dict((l,m,n) => 1-based grid index)`.
Thin adapter over [`NumericsCore.mode_index_table`](@ref): the ZNAmps axis is its un-mirrored
`m>=0` prefix `l_arr[1:num_modes]`.
"""
function kerr_mode_table(; lmax::Int = 10, nmax::Int = 55)
    tbl = mode_index_table(; lmax = lmax, nmax = nmax)
    nm = tbl.num_modes
    modes = NTuple{3,Int}[(tbl.l_arr[i], tbl.m_arr[i], tbl.n_arr[i]) for i in 1:nm]
    index = Dict{NTuple{3,Int},Int}(m => i for (i, m) in enumerate(modes))
    return modes, index
end

"""
    read_kerrampgrid(file; modes) -> KerrAmpGrid

Open the ZNAmps Kerr-ecc-eq amplitude grid and read the `regionA`/`regionB` clamped knot vectors and
the `CoeffsRegion{A,B}` B-spline coefficient slabs for the requested `(l,m,n)` `modes`. Method
provided by `DataHDF5Ext`.
"""
function read_kerrampgrid end

# Untyped fallback (signature `(Any,)`): the HDF5 ext supplies the strictly-more-specific
# `(AbstractString,)` method, refining dispatch without an illegal precompile-time method-overwrite.
read_kerrampgrid(file; kwargs...) = error("read_kerrampgrid: " * _NEED_HDF5)

export KerrAmpGrid, kerr_mode_table, read_kerrampgrid

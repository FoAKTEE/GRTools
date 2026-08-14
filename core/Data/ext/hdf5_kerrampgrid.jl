# ZNAmps Kerr eccentric-equatorial amplitude-grid HDF5 body (part of DataHDF5Ext). Moved VERBATIM
# from EMRIWaveforms/src/KerrAmpSpline.jl `read_kerrampgrid` (S4-b, design §7.2/§10: ALL grid reader
# BODIES live in Data's HDF5 weakdep-ext). Characterized bit-exactly by the pre-change fixture
# test/fixtures/hdf5_readers_ref.jl (baked from the committed EMRIWaveforms loader @ f409fc7).

"""
    Data.read_kerrampgrid(file; modes, lmax=attr, nmax=attr) -> KerrAmpGrid

Open the ZNAmps Kerr-ecc-eq amplitude grid and read the `regionA`/`regionB` clamped knot vectors and the
`CoeffsRegion{A,B}` B-spline coefficient slabs for the requested `(l,m,n)` `modes` (mapped to grid
indices via [`kerr_mode_table`]). Reads are serial (HDF5.jl is not thread-safe); the downstream bicubic
bank build is the threaded stage. Structural invariant: `n_coeff = n_knots − 4` per axis.
"""
function Data.read_kerrampgrid(file::AbstractString; modes)
    mode_list = NTuple{3,Int}[Tuple(Int.(m)) for m in modes]
    h5open(file, "r") do f
        lmax = Int(read_attribute(f, "lmax"))
        mmax = Int(read_attribute(f, "mmax"))
        nmax = Int(read_attribute(f, "nmax"))
        _, index = Data.kerr_mode_table(; lmax = lmax, nmax = nmax)
        midx = [index[m] for m in mode_list]

        gA = f["regionA"]; gB = f["regionB"]
        twA = Float64.(read(gA["w_knots"])); tuA = Float64.(read(gA["u_knots"]))
        twB = Float64.(read(gB["w_knots"])); tuB = Float64.(read(gB["u_knots"]))
        z_knots = Float64.(read(gA["z_knots"]))
        num_total = size(gA["CoeffsRegionA"], 3)   # Julia dims: (coefflat, 2, mode, z)

        dA = gA["CoeffsRegionA"]; dB = gB["CoeffsRegionB"]
        ncA = size(dA, 1); ncB = size(dB, 1); nz = length(z_knots)
        nm = length(mode_list)
        coeffsA = Array{Float64,4}(undef, ncA, 2, nm, nz)
        coeffsB = Array{Float64,4}(undef, ncB, 2, nm, nz)
        for (k, gi) in enumerate(midx)
            coeffsA[:, :, k, :] = dA[:, :, gi, :]   # serial HDF5 slab read
            coeffsB[:, :, k, :] = dB[:, :, gi, :]
        end
        return KerrAmpGrid{Float64}(mode_list, z_knots, twA, tuA, twB, tuB,
                                    coeffsA, coeffsB, lmax, mmax, nmax, num_total)
    end
end

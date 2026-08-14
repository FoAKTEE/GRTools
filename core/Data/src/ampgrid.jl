# Teukolsky amplitude-grid loader. The FEW (`amplitude/ampinterp2d.jl`) and KerrBackwards amplitude
# loaders use byte-identical schema strings → ONE loader over `AmpGridSchema` + the per-mode surface
# matrices. HDF5 body in ext/DataHDF5Ext.jl.

"""
    read_ampgrid(path) -> (AmpGridSchema, surfaces)

Read a Teukolsky amplitude grid: the `(w=e, u=y)` axes + `(l,m,n)` mode list as an `AmpGridSchema`,
and the stacked per-mode `length(w) × length(u)` complex surfaces. Method provided by `DataHDF5Ext`.
"""
function read_ampgrid end

# Untyped fallback (signature `(Any,)`): the HDF5 ext supplies the strictly-more-specific
# `(AbstractString,)` method, refining dispatch without an illegal precompile-time method-overwrite.
read_ampgrid(path; kwargs...) = error("read_ampgrid: " * _NEED_HDF5)

# u-axis coordinate of the FEW Schwarzschild-eccentric amplitude grid: y(p,e) = log(p − 2e − 2.1)
# (few.utils.mappings.schwarzecc.schwarzecc_p_to_y). Canonical single source: NumericsCore.p_to_y;
# the ext rounds to 8 digits like FEW.
const _amp_y = p_to_y

"""
    schwarzecc_amp_modes(; lmax=10, nmax=30) -> Vector{NTuple{3,Int}}

Stored `(l, m, n)` amplitude modes for the Schwarzschild-eccentric grid, in the file's storage order
(`m=0` block first, then `m>0`), restricted to `m ≥ 0`. Closed-form mode-discovery transcription of
`FastEMRIWaveforms.schwarzecc_amp_modes`.
"""
function schwarzecc_amp_modes(; lmax::Int = 10, nmax::Int = 30)
    base = NTuple{3,Int}[]
    for l in 2:lmax, m in 0:l, n in -nmax:nmax
        push!(base, (l, m, n))
    end
    m0 = [t for t in base if t[2] == 0]
    mpos = [t for t in base if t[2] != 0]
    return vcat(m0, mpos)
end

export read_ampgrid, schwarzecc_amp_modes

# H1Lorenz `(l,m)` mode I/O over the versioned `H1Schema` contract. De-dups the legacy readers in
# H1Lorenz (io_hdf5.jl write_mode_h5 / read_grid), SecondOrderRicci (read_h1), Punctures
# (read_hret_mode) and DiscEMRIForce (load_h1_mode) into one HDF5-backed pair. Bodies live in
# ext/DataHDF5Ext.jl; the stubs here error until HDF5 is loaded by the consumer.

const _NEED_HDF5 = "load HDF5 (`import HDF5`) — Data's HDF5 reader-writer lives in a weakdep extension"

"""
    read_h1(path; T=Float64) -> H1Schema{T}

Read one H1Lorenz `(l,m)` mode HDF5 file into an `H1Schema`. Method provided by `DataHDF5Ext`.
"""
function read_h1 end

"""
    write_mode_h5(path, s::H1Schema) -> path

Write an `H1Schema` to an HDF5 mode file in the versioned layout. Method provided by `DataHDF5Ext`.
"""
function write_mode_h5 end

# Untyped fallbacks (signature `(Any,...)`): the HDF5 ext supplies the strictly-more-specific
# `(AbstractString,...)` methods, so loading HDF5 refines dispatch WITHOUT method-overwrite (which is
# illegal during ext precompilation). Without HDF5 these fire and point at the missing trigger.
read_h1(path; kwargs...) = error("read_h1: " * _NEED_HDF5)
write_mode_h5(path, ::H1Schema) = error("write_mode_h5: " * _NEED_HDF5)

# --- pure closed-form index tables (HDF5-free; stay in main src) -------------------------------------

"""
    h1_parity_fields(l, m) -> Vector{Int}

Barack-Lousto basis indices stored (in column order) in the `inhom_left/right` datasets of an
`h1-l{l}m{m}.h5` file, selected by the parity of `l+m`. The `l=0`/`l=1` special cases drop the
higher-l-only fields. Closed-form transcription of `SecondOrderRicci.read_h1`'s `fields` block.
"""
function h1_parity_fields(l::Integer, m::Integer)
    if iseven(l + m)
        l == 1 ? [1, 3, 5, 6, 2, 4] :
        l == 0 ? [1, 3, 6, 2] :
                 [1, 3, 5, 6, 7, 2, 4]
    else
        l == 1 ? [9, 8] : [9, 10, 8]
    end
end

# Physical component i -> storage block k for the H1Lorenz no-ddh layout. Even-parity {1..7} live in
# l+m-even files (R1,R3,R5,R6,R7 coupled; R2,R4 gauge); odd-parity {8,9,10} in l+m-odd files (R9,R10
# coupled; R8 gauge). Dipole (l=1) and monopole (l=0) drop the higher-l-only fields, shifting the
# gauge blocks down — exactly as the Punctures source special-cases them.
const _HRET_BLOCK_EVEN = Dict(1 => 1, 3 => 2, 5 => 3, 6 => 4, 7 => 5, 2 => 6, 4 => 7)
const _HRET_BLOCK_EVEN_DIPOLE = Dict(1 => 1, 3 => 2, 5 => 3, 6 => 4, 2 => 5, 4 => 6)
const _HRET_BLOCK_MONOPOLE = Dict(1 => 1, 3 => 2, 6 => 3, 2 => 4)
const _HRET_BLOCK_ODD = Dict(9 => 1, 10 => 2, 8 => 3)
const _HRET_BLOCK_ODD_DIPOLE = Dict(9 => 1, 8 => 2)

"""
    hret_block(i, l) -> Int

Storage block (1-based) of physical component `i` in an H1Lorenz no-ddh mode file at multipole `l`.
Block `k` occupies dataset rows `(k-1)*4 .+ (1,2,3,4) = [Re h, Im h, Re ∂_r h, Im ∂_r h]`. Throws if
the component is absent at this `l`. Closed-form transcription of `Punctures.hret_block`.
"""
function hret_block(i::Integer, l::Integer)
    d = if i in (8, 9, 10)
        l == 1 ? _HRET_BLOCK_ODD_DIPOLE : _HRET_BLOCK_ODD
    else
        l == 0 ? _HRET_BLOCK_MONOPOLE :
        l == 1 ? _HRET_BLOCK_EVEN_DIPOLE : _HRET_BLOCK_EVEN
    end
    haskey(d, Int(i)) || throw(ArgumentError("hret_block: component i=$i absent at l=$l"))
    return d[Int(i)]
end

export read_h1, write_mode_h5, h1_parity_fields, hret_block

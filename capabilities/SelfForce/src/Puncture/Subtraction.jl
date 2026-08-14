# Regular-field subtraction hR = hret − hS (and ddhR = ddhret − ddhS) + the HDF5 retarded-mode read.
# Port of Punctures/subtraction.jl (h1S-exact-order-eps2.wl L226-316, L485-527, Wardell et al.
# arXiv:2508.00087), REDESIGNED:
#
#   • the physical-component→storage-block map `hret_block` (and the `_HRET_BLOCK_*` parity/dipole/
#     monopole dicts) is CONSUMED from the canonical `Data.hret_block` (AHA — Data already de-dups the
#     H1Lorenz/SecondOrderRicci/Punctures no-ddh block layout; re-exported here identity).
#   • the retarded-mode HDF5 read is CONSUMED from `Data.read_h1` (the single canonical H1Schema
#     reader): `read_hret_mode` decodes the schema's joined `inhom_left ++ inhom_right` block rows
#     `[Re h, Im h, Re ∂_r h, Im ∂_r h]` into the complex (hret, ∂_r hret) mode.
#
# `compute_hR`/`compute_ddhR` are the aligned elementwise differences on the worldtube grid.

import Data

# ── consumed physical-component → storage-block map (re-export identity; collapse legacy dicts) ────
const hret_block = Data.hret_block

"""
    read_hret_mode(path, i, l; T=Float64) -> (grid, hret, dhret)

Read the retarded mode component `i` (and its first r-derivative) from an H1Lorenz no-ddh mode file
at multipole `l`. Consumes `Data.read_h1` to load the `H1Schema`, joins `inhom_left ++ inhom_right`,
and decodes the `hret_block(i,l)` block rows `[Re h, Im h, Re ∂_r h, Im ∂_r h]` into the complex
`hret`, `dhret` vectors over the joined grid.
"""
function read_hret_mode(path::AbstractString, i::Integer, l::Integer; T::Type = Float64)
    s = Data.read_h1(path; T = T)
    full = hcat(s.inhom_left, s.inhom_right)     # (4nf, ngrid_joined)
    base = (hret_block(i, l) - 1) * 4
    hret = complex.(full[base + 1, :], full[base + 2, :])
    dhret = complex.(full[base + 3, :], full[base + 4, :])
    return (s.grid, hret, dhret)
end

"""
    compute_hR(hret, hS)

Regular-field mode `hR = hret − hS` on the worldtube grid (source L485). Both arguments must be
aligned, equal-length vectors (the caller slices the retarded mode to the worldtube before
subtracting).
"""
function compute_hR(hret::AbstractVector, hS::AbstractVector)
    length(hret) == length(hS) ||
        throw(DimensionMismatch("compute_hR: hret ($(length(hret))) vs hS ($(length(hS)))"))
    return hret .- hS
end

"""
    compute_ddhR(ddhret, ddhS)

Regular-field second r-derivative `ddhR = ddhret − ddhS` on the worldtube grid (source L515). Both
arguments must be aligned, equal-length vectors.
"""
function compute_ddhR(ddhret::AbstractVector, ddhS::AbstractVector)
    length(ddhret) == length(ddhS) ||
        throw(DimensionMismatch("compute_ddhR: length mismatch"))
    return ddhret .- ddhS
end

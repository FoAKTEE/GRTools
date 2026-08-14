# 12-column tab-separated `.dat` source-mode table I/O (Sourcefiletotable / mywrite,
# Sourceincoordinates.nb, Polcar & Witzany 2025, arXiv:2507.15720). Read side lives in Data
# (re-exported here); this file adds the write side, using Data's `_mma_sci` for byte-exact ω.

using Data: SourceMode, read_source_table
import Data: _mma_sci

"""
    source_mode_line(mode) -> String

Render one [`SourceMode`](@ref) as a single newline-terminated `.dat` record in `mywrite` column
order: `m, n, k, ω` (8 sig figs), then Re/Im of `Sl, Sn, Sm, Sm†` (full Float64 precision). The ω
column is formatted to 8 sig figs to reproduce the upstream file byte-for-byte.
"""
function source_mode_line(mode::SourceMode)
    cols = String[
        string(mode.m), string(mode.n), string(mode.k),
        _mma_sci(mode.ω; sigfigs=8),
        _mma_sci(real(mode.Sl)),    _mma_sci(imag(mode.Sl)),
        _mma_sci(real(mode.Sn)),    _mma_sci(imag(mode.Sn)),
        _mma_sci(real(mode.Sm)),    _mma_sci(imag(mode.Sm)),
        _mma_sci(real(mode.Smdag)), _mma_sci(imag(mode.Smdag)),
    ]
    return join(cols, '\t') * "\n"
end

"""
    write_source_mode(io, mode)

Append one mode record (see [`source_mode_line`](@ref)) to `io`.
"""
write_source_mode(io::IO, mode::SourceMode) = write(io, source_mode_line(mode))

"""
    write_source_table(io, modes)

Write a collection of `SourceMode`s, one `mywrite`-format record per line.
"""
function write_source_table(io::IO, modes)
    for mode in modes
        write_source_mode(io, mode)
    end
    return nothing
end

export SourceMode, read_source_table, source_mode_line, write_source_mode, write_source_table

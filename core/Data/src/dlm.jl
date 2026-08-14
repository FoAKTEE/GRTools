# `.dat`/`.csv` readers via DelimitedFiles. De-dups the scattered `readdlm` call sites — FEW
# `trajectory/flux.jl` (PN-scaled flux grid), `amplitude/romannet.jl` (AmplitudeVectorNorm) — and the
# ModifiedTeukolsky 12-column source-table reader. Pure stdlib (DelimitedFiles + Printf), no HDF5.
# Kept SEPARATE from the MMA `.wl`/`.m` readers (a different dialect — AHA, design §10).

import Printf

"""
    read_dat(path; kwargs...) -> Matrix

Read a whitespace-delimited numeric table (`.dat`) — a thin, faithful `DelimitedFiles.readdlm`
wrapper (returns a `Matrix{Float64}` when every cell is numeric). The de-duped loader behind the FEW
flux / amplitude-norm tables; column-convention reshaping is the caller's job (it differs per table).
"""
read_dat(path::AbstractString; kwargs...) = readdlm(path; kwargs...)

"""
    read_csv(path, T=Float64; kwargs...) -> Matrix{T}

Read a comma-delimited table (`.csv`) into a `Matrix{T}` — the `','`-delimited sibling of
[`read_dat`](@ref).
"""
read_csv(path::AbstractString, T::Type = Float64; kwargs...) = readdlm(path, ',', T; kwargs...)

# --- ModifiedGravity / ModifiedTeukolsky source-mode `.dat` table -----------------------------------

"""
    SourceMode

One row of the 12-column tab-separated source-mode table (`Sourcefiletotable`/`mywrite`,
Polcar & Witzany 2025, arXiv:2507.15720): integers `m,n,k`, mode frequency `ω`, then the four NP/GHP
projected source amplitudes `Sl,Sn,Sm,Sm†` as interleaved Re/Im pairs. `radius`/`ringmass` are
caller-supplied (not in the file); `Snorm` is 0 on read.
"""
struct SourceMode
    m::Int
    n::Int
    k::Int
    ω::Float64
    radius::Float64
    ringmass::Float64
    Sl::ComplexF64
    Sn::ComplexF64
    Sm::ComplexF64
    Smdag::ComplexF64   # Sm† — the GHP-prime of Sm
    Snorm::Float64
end

# Parse one tab-separated record into a SourceMode. `radius`/`ringmass` come from the caller.
function _parse_source_line(line::AbstractString, radius::Real, ringmass::Real)
    t = split(line, '\t')
    length(t) == 12 ||
        error("source table row needs 12 tab-separated fields, got $(length(t)): $(repr(line))")
    f(i) = parse(Float64, t[i])
    return SourceMode(
        parse(Int, t[1]), parse(Int, t[2]), parse(Int, t[3]),
        f(4), Float64(radius), Float64(ringmass),
        complex(f(5), f(6)),     # Sl
        complex(f(7), f(8)),     # Sn
        complex(f(9), f(10)),    # Sm
        complex(f(11), f(12)),   # Sm†
        0.0,                     # Snorm := 0 on read
    )
end

"""
    read_source_table(io_or_path; radius, ringmass) -> Vector{SourceMode}

Read every mode row from a source `.dat` file (path or `IO`); blank lines are skipped. `radius` and
`ringmass` are stored on each mode (caller-supplied, not in the file).
"""
function read_source_table(io::IO; radius::Real, ringmass::Real)
    modes = SourceMode[]
    for line in eachline(io)
        isempty(strip(line)) && continue
        push!(modes, _parse_source_line(line, radius, ringmass))
    end
    return modes
end

read_source_table(path::AbstractString; radius::Real, ringmass::Real) =
    open(io -> read_source_table(io; radius = radius, ringmass = ringmass), path)

# Format a real in Mathematica `mantissa "e" exponent` style (e.g. "2.4056261e-2"), trailing
# mantissa zeros trimmed. `sigfigs===nothing` => 17 sig figs (round-trips any Float64). The 2nd-user
# canonical sci-notation formatter promoted out of ModifiedTeukolsky (design §6 crosswalk row dlm).
function _mma_sci(x::Real; sigfigs::Union{Int,Nothing} = nothing)
    isfinite(x) || error("cannot format non-finite value $x")
    digits = sigfigs === nothing ? 16 : sigfigs - 1
    s = Printf.format(Printf.Format("%.$(digits)e"), Float64(x))
    mant, expstr = split(s, 'e')
    e = parse(Int, expstr)
    if occursin('.', mant)
        mant = rstrip(mant, '0')
        mant = rstrip(mant, '.')
    end
    return string(mant, "e", e)
end

export read_dat, read_csv, SourceMode, read_source_table

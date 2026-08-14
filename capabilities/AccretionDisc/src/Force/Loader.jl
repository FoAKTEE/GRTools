# LOADER — the in-memory complex BLS adapter `H1Mode{T}` over the on-disk `Data.H1Schema` contract.
# Loads the a=0 Lorenz-gauge Barack-Lousto field corpus (each `h1-l{l}m{m}.h5`) through the single
# canonical `Data.read_h1` reader (the HDF5 weakdep ext, activated because the Force submodule pulls
# HDF5 into the session), then widens the Float64 corpus to `Complex{T}` so the downstream type-generic
# assembly/force machinery carries BigFloat working precision. `H1Mode` is DISTINCT from `H1Schema`
# (design §4.10): the schema is the storage contract (real `inhom_left/right` tables), `H1Mode` the
# split complex `val_*`/`dr_*` fields the assembly consumes. The data root is the injected
# `Data.datadir()` (kills the legacy ref-code `H1DIR` `..`-escape). From DiscEMRIForce/loader.jl.

using Data: read_h1, datadir, datapath
import Data

"""
    H1Mode{T}

A loaded BLS Lorenz-gauge `(l,m)` mode: complex fields `val_right[j][k]` (`j=1..nfields`) and their
r-derivatives `dr_right[j][k]` on the right (`r≥r0`) and left (`r≤r0`) radial grids, with orbital
radius `r0`, asymptotic amplitudes `C_inf`, and `parity` (`:even` ⇒ `nfields=7`, `:odd` ⇒ `nfields=3`).
Low multipoles are zero-padded into these canonical parity layouts, so downstream assembly never
has to infer on-disk special-case column counts.
Element type `T` inherits from the requested precision (the Float64 corpus widened to `Complex{T}`).
"""
struct H1Mode{T<:Real}
    l::Int; m::Int
    r_right::Vector{T}; val_right::Vector{Vector{Complex{T}}}; dr_right::Vector{Vector{Complex{T}}}
    r_left::Vector{T};  val_left::Vector{Vector{Complex{T}}};  dr_left::Vector{Vector{Complex{T}}}
    r0::T; C_inf::Vector{T}; nfields::Int; parity::Symbol
end

# Split a (4nf, n) real BLS table into nf complex value/deriv column-vectors widened to Complex{T}.
function _split_bls(tab::AbstractMatrix, ::Type{T}) where {T}
    size(tab, 1) % 4 == 0 || throw(DimensionMismatch(
        "BLS table row count $(size(tab, 1)) is not divisible by four"))
    nf = size(tab, 1) ÷ 4
    cx(re, im) = Complex{T}.(T.(re), T.(im))
    val = [cx(tab[(j-1)*4+1, :], tab[(j-1)*4+2, :]) for j in 1:nf]
    dr  = [cx(tab[(j-1)*4+3, :], tab[(j-1)*4+4, :]) for j in 1:nf]
    return val, dr, nf
end

const _H1_EVEN_FIELDS = (1, 3, 5, 6, 7, 2, 4)
const _H1_ODD_FIELDS = (9, 10, 8)

# Map the compact storage layout (whose low-l fields are listed by Data.h1_parity_fields) into the
# fixed seven/even or three/odd layout consumed by Assemble.jl. Missing physical components are
# represented by zero vectors of the requested element type.
function _canonical_bls(tab::AbstractMatrix, ::Type{T}, l::Integer, m::Integer) where {T}
    rawval, rawdr, nf = _split_bls(tab, T)
    stored = Data.h1_parity_fields(l, m)
    nf == length(stored) || throw(DimensionMismatch(
        "H1 (l=$l,m=$m) stores $nf fields, expected $(length(stored)) from h1_parity_fields"))
    parity = iseven(l + m) ? :even : :odd
    canonical = parity === :even ? _H1_EVEN_FIELDS : _H1_ODD_FIELDS
    ngrid = size(tab, 2)
    val = [zeros(Complex{T}, ngrid) for _ in canonical]
    dr = [zeros(Complex{T}, ngrid) for _ in canonical]
    slots = Dict(field => j for (j, field) in enumerate(canonical))
    for (jraw, field) in enumerate(stored)
        j = slots[field]
        val[j] = rawval[jraw]
        dr[j] = rawdr[jraw]
    end
    return val, dr, length(canonical), parity
end

"""
    load_h1_mode(l, m; dir=datadir(), T=Float64) -> H1Mode{T}

Load the BLS Lorenz-gauge fields for `(l,m)` through `Data.read_h1` (the canonical H1Schema reader),
widened to `Complex{T}`. `val_right[j]` (`j=1..nfields`) = complex `h^{(j)}_{lm}(r)` on `r_right`
(`r≥r0`); `val_left` likewise (`r≤r0`); each carries its r-derivative (`dr_*`). `r0` = orbital radius.
Pass `T=BigFloat` (with the ambient `setprecision`) to widen for the downstream high-precision
assembly — the field VALUES stay Float64-accurate (the a=0 oracle ceiling).
"""
function load_h1_mode(l::Integer, m::Integer; dir::AbstractString = datadir(), T::Type = Float64)
    s = read_h1(joinpath(dir, "h1-l$(l)m$(m).h5"))               # H1Schema{Float64}
    valR, drR, nf, parity = _canonical_bls(s.inhom_right, T, l, m)
    valL, drL, nfL, parityL = _canonical_bls(s.inhom_left, T, l, m)
    (nfL == nf && parityL === parity) || error("inconsistent left/right H1 parity layout")
    grid = s.grid; r0i = s.r0index                               # 0-based
    rR = T.(grid[r0i+1:end]); rL = T.(grid[1:r0i+1])
    H1Mode{T}(Int(l), Int(m), rR, valR, drR, rL, valL, drL,
              T(grid[r0i+1]), T.(s.C_inf), nf, parity)
end

"""
    load_modes(lrange, m; dir=datadir(), T=Float64) -> Dict{Int,H1Mode{T}}

Load all `(l,m)` modes for `l in lrange` IN PARALLEL across threads (independent file reads;
race-free preallocated output slots, mirroring `DiscModes.solve_modes`). Returns `l => H1Mode`.
Parallel result is identical to the serial loop (each slot is written by exactly one task).
"""
function load_modes(lrange, m::Integer; dir::AbstractString = datadir(), T::Type = Float64)
    ls = collect(lrange)
    out = Vector{H1Mode{T}}(undef, length(ls))
    Threads.@threads for i in eachindex(ls)
        out[i] = load_h1_mode(ls[i], m; dir = dir, T = T)
    end
    Dict(ls[i] => out[i] for i in eachindex(ls))
end

"""
    available_modes(; dir=datadir()) -> Vector{NTuple{2,Int}}

Discover the `(l,m)` h1 mode files in `dir` (names `h1-l{l}m{m}.h5`), sorted. Reuses the canonical
`Data.available_modes` (the datadir-injected scanner — no `H1DIR` escape) — AHA.
"""
available_modes(; dir::AbstractString = datadir()) = Data.available_modes(; dir = dir)

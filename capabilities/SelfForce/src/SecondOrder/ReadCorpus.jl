# Directory-level h1-corpus reader (legacy SecondOrderRicci/read_h1.jl) + the zero-pad field shim.
# Assembles the per-mode HDF5 files h1-l{l}m{m}.h5 in an injected `dir` into the corpus struct the
# r2 kernels consume. `h1_parity_fields` is the canonical Data primitive (re-exported here through
# the `using Data: h1_parity_fields` in SecondOrder.jl — NOT re-owned). `a_il` (trace-reversal
# normalization) comes from NumericsCore via the SecondOrder front door.
#
# This reader unpacks `inhom_left/right` DIRECTLY via HDF5 because two compatible on-disk layouts
# exist.  The captured SecondOrderRicci corpus stores 6 rows per field (Re/Im of h̄, ∂_r h̄, ∂²_r h̄),
# while the production `Data.H1Schema` written by `LorenzField` stores 4 rows per field (Re/Im of
# h̄, ∂_r h̄).  For the 4-row layout the missing radial jet is closed with the independently
# oracle-tested Lorenz field equation `d2hdr2`.  The special l>=2 even-static sector is deliberately
# rejected in that layout until its algebraically reconstructed fields have a validated second-
# derivative path; legacy 6-row files remain the supported route for those modes.

import HDF5

"""
    ZeroPadField{T}(d, z)

Zero-padded, read-only view of an h1-corpus field dictionary
`d::Dict{NTuple{3,Int},Vector{Complex{T}}}` keyed by `(i, l, m)`. Indexing with a key absent from
the corpus returns the single shared zero vector `z` instead of throwing, mirroring the legacy C++
driver's zero-initialized `boost::multi_array` (modes outside the computed `(l, m)` range contribute
0 to the quadratic `r2` sums). The `r2` kernels only READ fields, so aliasing the shared zero
across all missing keys is safe. Supports `getindex` only.
"""
struct ZeroPadField{T}
    d::Dict{NTuple{3,Int},Vector{Complex{T}}}
    z::Vector{Complex{T}}
end
Base.getindex(p::ZeroPadField{T}, k::NTuple{3,Int}) where {T} = get(p.d, k, p.z)

"""
    H1RadialJetCorpus{T}

Typed in-memory first-order Lorenz corpus consumed by the second-order Ricci kernels.  `h`, `dh`,
and `ddh` are the physical (non-trace-reversed) radial jet, keyed by `(i,l,m)`; `r` includes the two
one-sided copies of the worldline point.  Both legacy 6-row files and production 4-row
`Data.H1Schema` files are adapted to this single representation by [`read_h1_corpus`](@ref).
"""
struct H1RadialJetCorpus{T<:AbstractFloat}
    r0::T
    r::Vector{T}
    f::Vector{T}
    fp::Vector{T}
    h::Dict{NTuple{3,Int},Vector{Complex{T}}}
    dh::Dict{NTuple{3,Int},Vector{Complex{T}}}
    ddh::Dict{NTuple{3,Int},Vector{Complex{T}}}
end

_is_unvalidated_even_static(l::Integer, m::Integer) = l >= 2 && m == 0 && iseven(l)

function _four_row_ddbar!(ddbar, hbar, dhbar, fields, l, m, r, f, fp, r0, ::Type{T}) where {T}
    _is_unvalidated_even_static(l, m) && throw(ArgumentError(
        "read_h1_corpus: production 4-row H1Schema mode (l=$l,m=$m) is even-static; " *
        "its reconstructed R6/R7 second radial derivatives are not yet validated. " *
        "Provide a legacy 6-row radial-jet file for this mode."))

    hb = zeros(Complex{T}, 10)
    db = zeros(Complex{T}, 10)
    for j in eachindex(r)
        fill!(hb, zero(Complex{T}))
        fill!(db, zero(Complex{T}))
        for i in fields
            hb[i] = hbar[i][j]
            db[i] = dhbar[i][j]
        end
        for i in fields
            ddbar[i][j] = d2hdr2(T, i, l, m, r[j], f[j], fp[j], r0, hb, db)
        end
    end
    return ddbar
end

"""
    read_h1_corpus([T=Float64,] dir; l_max) -> H1RadialJetCorpus{T}

Ingest the first-order Lorenz-gauge fields written as `h1-l{l}m{m}.h5` files in directory `dir` (the
injected corpus datadir). Type-generic port of `read_h1_corpus`
(`ref-code/SecondOrderRicci/src/read_h1.jl` / `h1.cc` L125-337).

Returns the worldline radius `r0`, the radial grid `r` (with `r0` duplicated for the left/right
solutions), `f = 1-2M/r`, `fp = 2M/r²`, and the *non*-trace-reversed physical field `h` with its
first/second radial derivatives `dh,ddh` — each a `Dict{NTuple{3,Int},Vector{Complex{T}}}` keyed by
`(i,l,m)`. `M = 1`.

Accepted storage layouts and conventions:
- grid: the `"grid"` dataset holds `N-1` points; `r0` sits at index `N_left-1` (0-based, from the
  `inhom_left` column count) and is duplicated so the right solutions start at the repeated `r0`;
- legacy `inhom_left/right` tables store 6 rows per field
  `(Re,Im, Re∂_r,Im∂_r, Re∂²_r,Im∂²_r)`;
- production `Data.H1Schema` tables store 4 rows per field `(Re,Im, Re∂_r,Im∂_r)` and
  `d2hdr2` supplies ∂²_r on each one-sided grid point.  Four-row even-static modes with
  `l>=2,m=0,l` even throw a precise `ArgumentError` rather than applying an unvalidated equation;
- trace-reversal swaps basis indices 3↔6 when applying `a_il` and the `1/r` factors:
  `h = a·h̄/r`, `dh = a(∂_r h̄ - h̄/r)/r`, `ddh = a(∂²_r h̄ - 2(∂_r h̄ - h̄/r)/r)/r`;
- negative-`m`: `h[i][l][-m] = (-1)^m conj(h[i][l][m])` (the `m=0` case conjugates in place).
"""
function read_h1_corpus(::Type{T}, dir::AbstractString; l_max::Integer) where {T<:AbstractFloat}
    M = one(T)
    pat = r"^h1-l([0-9]+)m([0-9]+)\.h5$"
    names = sort([fn for fn in readdir(dir) if occursin(pat, fn)])
    isempty(names) && error("read_h1_corpus: no h1-l{l}m{m}.h5 files in $dir")
    modes = NTuple{2,Int}[]
    for fn in names
        mm = match(pat, fn)
        l = parse(Int, mm[1]); m = parse(Int, mm[2])
        l > l_max && continue
        push!(modes, (l, m))
    end
    isempty(modes) && error("read_h1_corpus: no modes with l ≤ $l_max in $dir")

    # --- grid structure from the first selected mode ---
    first_l, first_m = first(modes)
    first_path = joinpath(dir, "h1-l$(first_l)m$(first_m).h5")
    grid, N_left0 = HDF5.h5open(first_path, "r") do file
        g = read(file["grid"])
        L = read(file["inhom_left"])   # (4· or 6·num_fields, N_left)
        (Vector{T}(g), size(L, 2))
    end
    G = length(grid)
    N = G + 1
    r0i = N_left0 - 1                  # 0-based worldline index
    r = Vector{T}(undef, N)
    @views r[1:G] .= grid
    r0 = r[r0i + 1]
    for ci in (N - 1):-1:(r0i + 1)     # 0-based duplication loop: r[ci]=r[ci-1]
        r[ci + 1] = r[ci]
    end
    f = @. 1 - 2M / r
    fp = @. 2M / r^2

    h   = Dict{NTuple{3,Int},Vector{Complex{T}}}()
    dh  = Dict{NTuple{3,Int},Vector{Complex{T}}}()
    ddh = Dict{NTuple{3,Int},Vector{Complex{T}}}()

    for (l, m) in modes
        path = joinpath(dir, "h1-l$(l)m$(m).h5")
        mode_grid, L, R = HDF5.h5open(path, "r") do file
            (Vector{T}(read(file["grid"])),
             read(file["inhom_left"]), read(file["inhom_right"]))
        end
        length(mode_grid) == length(grid) || throw(DimensionMismatch(
            "read_h1_corpus: radial grid length $(length(mode_grid)) for (l=$l,m=$m) differs from " *
            "first selected mode (l=$first_l,m=$first_m) length $(length(grid))"))
        if mode_grid != grid
            j = findfirst(k -> mode_grid[k] != grid[k], eachindex(grid))
            throw(ArgumentError(
                "read_h1_corpus: radial grid mismatch for (l=$l,m=$m) at index $j: " *
                "got $(mode_grid[j]), expected $(grid[j]) from first selected mode " *
                "(l=$first_l,m=$first_m); all corpus modes must share an identical grid"))
        end
        N_left = size(L, 2); N_right = size(R, 2)
        N_left + N_right == N || throw(DimensionMismatch(
            "read_h1_corpus: N_left+N_right=$(N_left + N_right) ≠ N=$N for (l=$l,m=$m)"))
        N_left == r0i + 1 || throw(DimensionMismatch(
            "read_h1_corpus: worldline split N_left=$N_left differs from corpus N_left=$(r0i + 1) for (l=$l,m=$m)"))
        fields = h1_parity_fields(l, m)
        nf = length(fields)
        nrows = size(L, 1)
        size(R, 1) == nrows || throw(DimensionMismatch(
            "read_h1_corpus: left/right row counts differ ($(size(L,1)) vs $(size(R,1))) for (l=$l,m=$m)"))
        row_stride = if nrows == 6nf
            6
        elseif nrows == 4nf
            4
        else
            throw(DimensionMismatch(
                "read_h1_corpus: row count $nrows must be 4·#fields=$(4nf) (H1Schema) or " *
                "6·#fields=$(6nf) (legacy radial jet) for (l=$l,m=$m)"))
        end

        # raw trace-reversed bar field and its derivatives for this (l,m)
        hbar  = Dict(i => zeros(Complex{T}, N) for i in fields)
        dhbar = Dict(i => zeros(Complex{T}, N) for i in fields)
        ddbar = Dict(i => zeros(Complex{T}, N) for i in fields)
        for (it, i) in enumerate(fields)
            c = row_stride * (it - 1)
            for j in 0:(N_left - 1)
                hbar[i][j + 1]  = complex(T(L[c + 1, j + 1]), T(L[c + 2, j + 1]))
                dhbar[i][j + 1] = complex(T(L[c + 3, j + 1]), T(L[c + 4, j + 1]))
                row_stride == 6 &&
                    (ddbar[i][j + 1] = complex(T(L[c + 5, j + 1]), T(L[c + 6, j + 1])))
            end
            for j in 0:(N_right - 1)
                hbar[i][N_left + j + 1]  = complex(T(R[c + 1, j + 1]), T(R[c + 2, j + 1]))
                dhbar[i][N_left + j + 1] = complex(T(R[c + 3, j + 1]), T(R[c + 4, j + 1]))
                row_stride == 6 &&
                    (ddbar[i][N_left + j + 1] = complex(T(R[c + 5, j + 1]), T(R[c + 6, j + 1])))
            end
        end
        row_stride == 4 && _four_row_ddbar!(ddbar, hbar, dhbar, fields, l, m, r, f, fp, r0, T)

        # physical (non-trace-reversed) field: 3↔6 swap, a_il + 1/r factors
        for i in fields
            ih = i == 3 ? 6 : (i == 6 ? 3 : i)
            a = a_il(T, ih, l)
            key = (ih, l, m)
            hk = zeros(Complex{T}, N); dhk = zeros(Complex{T}, N); ddhk = zeros(Complex{T}, N)
            for j in 1:N
                hb = hbar[i][j]; db = dhbar[i][j]; ddb = ddbar[i][j]; rj = r[j]
                hk[j]   = a * hb / rj
                dhk[j]  = a * (db - hb / rj) / rj
                ddhk[j] = a * (ddb - 2 * (db - hb / rj) / rj) / rj
            end
            h[key] = hk; dh[key] = dhk; ddh[key] = ddhk
        end

        # negative-m (and m=0 in-place conjugation), iterating over the bar index set
        sgn = isodd(m) ? -one(T) : one(T)
        for i in fields
            src = (i, l, m); dst = (i, l, -m)
            h[dst]   = sgn .* conj.(h[src])
            dh[dst]  = sgn .* conj.(dh[src])
            ddh[dst] = sgn .* conj.(ddh[src])
        end
    end

    return H1RadialJetCorpus{T}(r0, r, f, fp, h, dh, ddh)
end

read_h1_corpus(dir::AbstractString; l_max::Integer) = read_h1_corpus(Float64, dir; l_max = l_max)
read_h1_corpus(dir::AbstractString, T::Type; l_max::Integer) = read_h1_corpus(T, dir; l_max = l_max)

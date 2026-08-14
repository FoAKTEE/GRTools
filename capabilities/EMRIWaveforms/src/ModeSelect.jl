"""
    EMRIWaveforms.ModeSelect

Power-thresholded harmonic mode selection (design §6.8: FEW `modeselector.jl`; KerrBackwards
`modeselector.jl`). v0.1 keeps the two convention forks behind the qualified deep namespaces
`ModeSelect.KerrBackwards` and `ModeSelect.FEW` (§4.8 collisions — never flattened).

The canonical power-threshold algorithm (`|A·y|²` per mode, sort descending, keep until the cumulative
power reaches `(1−eps)` of the row total, union over trajectory points, `±m` partner pairing) lives in
[`KerrBackwards.mode_selector`]; the `FEW` fork threads the extra polar index `k` through and adds the
`"all"` / explicit-mode-list branches over it.
"""
module ModeSelect

export select_modes_by_list

# =====================================================================================
# KerrBackwards fork — the canonical power-threshold selector (equatorial (l,m,n)).
# =====================================================================================
module KerrBackwards

export mode_selector

"""
    mode_selector(teuk_modes, ylms, l_arr, m_arr, n_arr; eps=1e-5)
        -> (; teuk_modes, ylms, l_arr, m_arr, n_arr, keep_modes)

Power-threshold mode filter (FEW `ModeSelector.__call__`, unweighted path). For every trajectory point
the per-mode power `|A·y|²` is formed for the `K` base (`m ≥ 0`) modes and their `m < 0` partners (from
`conj(A)`), sorted descending, and modes are kept until the cumulative power reaches `(1−eps)` of the row
total; a mode kept at any point is kept overall, and a `±m` hit keeps the shared base mode. Base modes
may appear in arbitrary order; partner columns are mapped through the actual `m > 0` base indices.

`ylms` has length `K + #{m>0}` (base harmonics then the `m < 0` partners). Returns the kept base columns,
the kept `ylms` in FEW's `ylmkeep = [keep_modes, temp2]` layout (length `2·#kept`), the kept `(l,m,n)`,
and the 1-based `keep_modes`.
"""
function mode_selector(teuk_modes::AbstractMatrix, ylms::AbstractVector,
                       l_arr::AbstractVector{<:Integer}, m_arr::AbstractVector{<:Integer},
                       n_arr::AbstractVector{<:Integer}; eps::Real = 1e-5)
    (isfinite(eps) && 0 <= eps <= 1) ||
        throw(ArgumentError("eps must be finite and lie in [0,1], got $eps"))
    nt, K = size(teuk_modes)
    (length(l_arr) == K && length(m_arr) == K && length(n_arr) == K) ||
        throw(ArgumentError("l_arr/m_arr/n_arr must have length size(teuk_modes,2) = $K"))
    all(m -> m >= 0, m_arr) ||
        throw(ArgumentError("m_arr must describe the nonnegative-m base modes"))
    mpos = findall(m -> m > 0, m_arr)
    num_m_1_up = length(mpos)
    length(ylms) == K + num_m_1_up ||
        throw(ArgumentError("ylms must have length K + (#m>0) = $(K + num_m_1_up)"))

    ncol = K + num_m_1_up
    # Column-to-base map for [base harmonics; negative-m partner harmonics]. This explicit map is
    # the key invariant when zero- and positive-m base modes are interleaved.
    base_for_power_col = vcat(collect(1:K), mpos)
    partner_col = collect(1:K) # m=0 maps to its own ylm in the 2K output layout
    @inbounds for (p, jb) in enumerate(mpos)
        partner_col[jb] = K + p
    end

    power = Matrix{Float64}(undef, nt, ncol)
    @inbounds for i in 1:nt
        for j in 1:K
            power[i, j] = abs2(teuk_modes[i, j] * ylms[j])
        end
        for (p, jb) in enumerate(mpos)
            power[i, K + p] = abs2(conj(teuk_modes[i, jb]) * ylms[K + p])
        end
    end

    keepset = Set{Int}()
    for i in 1:nt
        row = @view power[i, :]
        order = sortperm(row; rev = true)
        thresh = sum(row) * (1 - eps)
        cum = 0.0
        for k in 1:ncol
            keep = (k == 1) || (cum < thresh)
            if keep
                push!(keepset, base_for_power_col[order[k]])
            end
            cum += row[order[k]]
        end
    end

    keep_modes = sort!(collect(keepset))

    ylmkeep = Vector{Int}(undef, 2 * length(keep_modes))
    @inbounds for (q, km) in enumerate(keep_modes)
        ylmkeep[q] = km
        ylmkeep[length(keep_modes) + q] = partner_col[km]
    end

    return (teuk_modes = teuk_modes[:, keep_modes], ylms = ylms[ylmkeep],
            l_arr = l_arr[keep_modes], m_arr = m_arr[keep_modes],
            n_arr = n_arr[keep_modes], keep_modes = keep_modes)
end

end # module KerrBackwards

# =====================================================================================
# select_modes_by_list — shared explicit-list membership (ModeSelect top level).
# =====================================================================================

"""
    select_modes_by_list(l_arr, m_arr, k_arr, n_arr, requested) -> Vector{Int}

1-based indices of the modes in `requested` (a list of `(l, m, k, n)` tuples) within the mode grid;
errors if any requested mode is absent.
"""
function select_modes_by_list(l_arr, m_arr, k_arr, n_arr, requested)
    keep = Int[]
    for tup in requested
        l, m, k, n = tup
        idx = findfirst(i -> l_arr[i] == l && m_arr[i] == m && k_arr[i] == k && n_arr[i] == n,
                        eachindex(l_arr))
        idx === nothing &&
            throw(ArgumentError("requested mode $(tup) not present in the mode grid"))
        push!(keep, idx)
    end
    return keep
end

# =====================================================================================
# FEW fork — adds the polar index k + "all"/explicit-list branches over the KB selector.
# =====================================================================================
module FEW

using ..KerrBackwards: mode_selector as _kb_mode_selector
using ..ModeSelect: select_modes_by_list

export mode_selector

"""
    mode_selector(teuk_modes, ylms, l_arr, m_arr, k_arr, n_arr; eps=1e-5, mode_selection=nothing)

FEW mode selection over a `(l, m, k, n)` grid. `mode_selection=nothing` delegates to the power threshold
([`ModeSelect.KerrBackwards.mode_selector`]) threading `k` through; `"all"` keeps every mode; a vector of
`(l, m, k, n)` tuples keeps exactly those (via [`select_modes_by_list`]).
"""
function mode_selector(teuk_modes::AbstractMatrix, ylms::AbstractVector,
                       l_arr::AbstractVector{<:Integer}, m_arr::AbstractVector{<:Integer},
                       k_arr::AbstractVector{<:Integer}, n_arr::AbstractVector{<:Integer};
                       eps::Real = 1e-5, mode_selection = nothing)
    (isfinite(eps) && 0 <= eps <= 1) ||
        throw(ArgumentError("eps must be finite and lie in [0,1], got $eps"))
    K = size(teuk_modes, 2)
    (length(l_arr) == K && length(m_arr) == K && length(k_arr) == K && length(n_arr) == K) ||
        throw(ArgumentError("l_arr/m_arr/k_arr/n_arr must have length size(teuk_modes,2) = $K"))
    all(m -> m >= 0, m_arr) ||
        throw(ArgumentError("m_arr must describe the nonnegative-m base modes"))
    mpos = findall(>(0), m_arr)
    length(ylms) == K + length(mpos) ||
        throw(ArgumentError("ylms must have length K + (#m>0) = $(K + length(mpos))"))
    if mode_selection === nothing
        kb = _kb_mode_selector(teuk_modes, ylms, l_arr, m_arr, n_arr; eps = eps)
        keep = kb.keep_modes
        return (teuk_modes = kb.teuk_modes, ylms = kb.ylms, l_arr = kb.l_arr,
                m_arr = kb.m_arr, k_arr = collect(k_arr)[keep], n_arr = kb.n_arr,
                keep_modes = keep)
    elseif mode_selection == "all"
        keep = collect(1:K)
    elseif mode_selection isa AbstractVector
        keep = select_modes_by_list(l_arr, m_arr, k_arr, n_arr, mode_selection)
    else
        throw(ArgumentError("mode_selection must be nothing, \"all\", or a vector of (l,m,k,n)"))
    end
    partner_col = collect(1:K)
    @inbounds for (p, jb) in enumerate(mpos)
        partner_col[jb] = K + p
    end
    ylmkeep = vcat(keep, partner_col[keep])
    return (teuk_modes = teuk_modes[:, keep], ylms = ylms[ylmkeep],
            l_arr = collect(l_arr)[keep],
            m_arr = collect(m_arr)[keep], k_arr = collect(k_arr)[keep],
            n_arr = collect(n_arr)[keep], keep_modes = keep)
end

end # module FEW

end # module ModeSelect

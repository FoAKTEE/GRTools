# (l,m,n) mode-index tables + (l,m)<->gamma enumeration + Barack-Lousto a_il + angular ladder factors.
#
# Consolidates the mode-bookkeeping primitives shared across the GSF/eccentric-amplitude layers:
# H1Lorenz's gamma work-counter, KerrBackwards' FEW (l,m,n) mirrored table, SecondOrderRicci's
# tensor-harmonic normalization, and the spin-1/spin-2 ladder factors used by the disc/Lorenz assemblies.

# --- (l,m) <-> gamma enumeration (gamma = l(l+1)/2 + m, m in 0:l) ---

"""`gamma_to_lm(g)` — decode the flat mode counter `gamma = l(l+1)/2 + m` into `(l, m)`."""
function gamma_to_lm(g::Integer)
    l = 0
    while l*(l+1)÷2 <= g; l += 1; end
    l -= 1
    return (l, g - l*(l+1)÷2)
end

"""`lm_to_gamma(l, m)` — encode `(l, m)` as `gamma = l(l+1)/2 + m`."""
lm_to_gamma(l::Integer, m::Integer) = l*(l+1)÷2 + m

"""`mode_list(lmax)` — `(l, m)` pairs for `l in 0:lmax`, `m in 0:l`."""
mode_list(lmax::Integer) = [(l, m) for l in 0:lmax for m in 0:l]

# --- FEW (l,m,n) mirrored mode-index table (SchwarzschildEccentric base class) ---

"""
    ModeIndexTable

The `(l,m,n)` mode grid of the FEW eccentric amplitude/mode-selection layers. Construction
(`lmax=10, nmax=30`): enumerate `l in 2:lmax, m in 0:l, n in -nmax:nmax`; reorder `m==0` first
(`m0sort`); mirror the `m!=0` block with negated `m`. Indices are 1-based.

- `l_arr,m_arr,n_arr`: full mirrored vectors (length `num_modes + num_m_1_up`).
- `num_modes`: the un-mirrored `m>=0` length; `num_m0`/`num_m_1_up`: `m==0`/`m>0` counts.
- `index_map`: `(l,m,n)->i`; `special_index_map`: `(l,m,n)->` its `m>0` partner position.
"""
struct ModeIndexTable
    l_arr::Vector{Int}
    m_arr::Vector{Int}
    n_arr::Vector{Int}
    num_modes::Int
    num_m0::Int
    num_m_1_up::Int
    num_m_zero_up::Int
    index_map::Dict{Tuple{Int,Int,Int},Int}
    special_index_map::Dict{Tuple{Int,Int,Int},Int}
end

"""`mode_index_table(; lmax=10, nmax=30)` — build the FEW `(l,m,n)` mode-index table."""
function mode_index_table(; lmax::Int=10, nmax::Int=30)
    L = Int[]; M = Int[]; N = Int[]
    for l in 2:lmax, m in 0:l, n in -nmax:nmax
        push!(L, l); push!(M, m); push!(N, n)
    end
    num_modes = length(L)

    zero_idx = [i for i in 1:num_modes if M[i] == 0]
    nz_idx   = [i for i in 1:num_modes if M[i] != 0]
    order = vcat(zero_idx, nz_idx)
    l_arr = L[order]; m_arr = M[order]; n_arr = N[order]

    num_m0 = length(zero_idx)
    num_m_zero_up = num_modes
    num_m_1_up = num_m_zero_up - num_m0

    nzmask = [i for i in 1:num_modes if m_arr[i] != 0]
    l_arr = vcat(l_arr, l_arr[nzmask])
    m_arr = vcat(m_arr, -m_arr[nzmask])
    n_arr = vcat(n_arr, n_arr[nzmask])

    index_map = Dict{Tuple{Int,Int,Int},Int}()
    special_index_map = Dict{Tuple{Int,Int,Int},Int}()
    for i in eachindex(l_arr)
        key = (l_arr[i], m_arr[i], n_arr[i])
        index_map[key] = i
        special_index_map[key] = i <= num_modes ? i : i - num_m_1_up
    end

    return ModeIndexTable(l_arr, m_arr, n_arr, num_modes, num_m0, num_m_1_up,
                          num_m_zero_up, index_map, special_index_map)
end

# --- Barack-Lousto tensor-harmonic normalization + angular ladder factors ---

"""`a_il(::Type{T}, i, l)` — Barack-Lousto tensor-harmonic normalization for basis index `i` at
multipole `l`: base `2^(-1/2)`, divided by `sqrt(l(l+1))` for `i in {4,5,8,9}` and by
`sqrt((l-1)l(l+1)(l+2))` for `i in {7,10}`; unscaled for `i in {1,2,3,6}`. (SecondOrderRicci `a_il`.)"""
function a_il(::Type{T}, i::Integer, l::Integer) where {T<:AbstractFloat}
    a = sqrt(one(T)/2)
    if i in (4, 5, 8, 9)
        a /= sqrt(T(l*(l + 1)))
    elseif i in (7, 10)
        a /= sqrt(T((l - 1)*l*(l + 1)*(l + 2)))
    end
    return a
end
a_il(i::Integer, l::Integer) = a_il(Float64, i, l)

"""`lambdahat(::Type{T}, l)` — spin-1 angular ladder factor `sqrt(l(l+1))`."""
lambdahat(::Type{T}, l::Integer) where {T<:AbstractFloat} = sqrt(T(l*(l + 1)))
lambdahat(l::Integer) = lambdahat(Float64, l)

"""`Lambdahat(::Type{T}, l)` — spin-2 angular ladder factor `sqrt((l-1)l(l+1)(l+2))`."""
Lambdahat(::Type{T}, l::Integer) where {T<:AbstractFloat} = sqrt(T((l - 1)*l*(l + 1)*(l + 2)))
Lambdahat(l::Integer) = Lambdahat(Float64, l)

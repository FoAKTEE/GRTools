# Teukolsky amplitude-grid HDF5 body (part of DataHDF5Ext). De-dups the byte-identical FEW-vs-
# KerrBackwards amplitude loaders into one reader over AmpGridSchema, reproducing the numpy C-order
# reshape + u-axis reversal exactly.

"""
    Data.read_ampgrid(path; modes=schwarzecc_amp_modes()) -> (AmpGridSchema, re, im)

Read a Teukolsky amplitude grid (FEW `ampinterp2d` layout): the `grid` dataset `(4,N)` gives the
`e` (=`w`) and `y(p,e)` (=`u`) interpolation axes; each `l{l}m{m}/n{n}k0` dataset `(2,N)` holds the
real/imag samples. Returns the `(w,u,modes)` `AmpGridSchema` plus the per-mode `num_w × num_u`
real/imag surfaces, reshaped with the numpy C-order + u-reversal index
`idx = (i-1)*num_u + (num_u-j) + 1`.
"""
function Data.read_ampgrid(path::AbstractString; modes = Data.schwarzecc_amp_modes())
    mode_list = NTuple{3,Int}[Tuple(Int.(m)) for m in modes]
    re = Vector{Matrix{Float64}}(undef, length(mode_list))
    im = Vector{Matrix{Float64}}(undef, length(mode_list))
    local unique_w, unique_u, num_w, num_u
    h5open(path, "r") do f
        grid = read(f["grid"])            # (4, N): rows a, p, e, xI
        p = grid[2, :]
        e = grid[3, :]
        uvals = round.(Data._amp_y.(p, e); digits = 8)
        unique_w = sort(unique(e))
        unique_u = sort(unique(uvals))
        num_w = length(unique_w)
        num_u = length(unique_u)
        for (k, (l, m, n)) in enumerate(mode_list)
            d = read(f["l$(l)m$(m)/n$(n)k0"])    # (2, N): real, imag
            valr = @view d[1, :]
            vali = @view d[2, :]
            Mr = Matrix{Float64}(undef, num_w, num_u)
            Mi = Matrix{Float64}(undef, num_w, num_u)
            for i in 1:num_w, j in 1:num_u
                idx = (i - 1) * num_u + (num_u - j) + 1
                Mr[i, j] = valr[idx]
                Mi[i, j] = vali[idx]
            end
            re[k] = Mr
            im[k] = Mi
        end
    end
    return AmpGridSchema{Float64}(unique_w, unique_u, mode_list), re, im
end

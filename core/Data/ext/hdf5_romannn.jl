# ROMAN NN-weight HDF5 body (part of DataHDF5Ext). Moved VERBATIM from
# EMRIWaveforms/src/AmplitudeModels.jl `_load_roman_nn` (S4-b, design §7.2/§10: ALL grid reader
# BODIES live in Data's HDF5 weakdep-ext); only the return is named (positionally identical 6-tuple).
# Characterized bit-exactly by the pre-change fixture test/fixtures/hdf5_readers_ref.jl.

"""
    Data.read_roman_nn(path) -> NamedTuple

Read the dense ROMAN network (`SchwarzschildEccentricInput.hdf5`): `w1..wL` / `b1..bL` + the complex
`reduced_basis` and the scalar attributes (`num_teuk_modes`, `break_index`, `transform_factor`).
Returns `(weights, biases, reduced_basis, transform_factor_inv, break_index, num_teuk_modes)`.
"""
function Data.read_roman_nn(nn_file::AbstractString)
    weights = Matrix{Float64}[]
    biases = Vector{Float64}[]
    local reduced_basis, transform_factor_inv, break_index, num_teuk_modes
    h5open(nn_file, "r") do f
        num_teuk_modes = Int(read_attribute(f, "num_teuk_modes"))
        transform_factor = read_attribute(f, "transform_factor")
        break_index = Int(read_attribute(f, "break_index"))
        transform_factor_inv = 1 / transform_factor
        num_layers = 0
        for k in keys(f)
            k == "reduced_basis" && continue
            num_layers = max(num_layers, parse(Int, k[2:end]))
        end
        for i in 1:num_layers
            push!(weights, Float64.(read(f["w$i"])))     # (n_out, n_in)
            push!(biases, Float64.(read(f["b$i"])))
        end
        reduced_basis = ComplexF64.(read(f["reduced_basis"]))  # (num_teuk_modes, break_index)
    end
    return (weights = weights, biases = biases, reduced_basis = reduced_basis,
            transform_factor_inv = transform_factor_inv, break_index = break_index,
            num_teuk_modes = num_teuk_modes)
end

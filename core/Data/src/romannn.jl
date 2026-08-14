# ROMAN NN-weight loader (S4-b, design §7.2/§10: ALL grid reader BODIES live in Data's HDF5
# weakdep-ext). Moved from EMRIWaveforms/src/AmplitudeModels.jl `_load_roman_nn` — the dense ROMAN
# (reduced-order-model artificial-neuron) network stored in `SchwarzschildEccentricInput.hdf5`.
# The `h5open` body is in ext/hdf5_romannn.jl; the GEMM forward pass stays in `Kernels.roman_forward`.

"""
    read_roman_nn(path) -> NamedTuple

Read the dense ROMAN network file (`SchwarzschildEccentricInput.hdf5`): layers `w1..wL` / `b1..bL`
plus a complex `reduced_basis` and the scalar attributes. Returns the NamedTuple
`(weights, biases, reduced_basis, transform_factor_inv, break_index, num_teuk_modes)` — positionally
identical to the legacy `_load_roman_nn` 6-tuple. Julia's column-major read transposes the HDF5
row-major arrays, giving weights as `(n_out, n_in)` and `reduced_basis` as
`(num_teuk_modes, break_index)` — exactly the layout `Kernels.roman_forward` consumes. (Distinct
schema from the `Teuk_amps` grid, so this is NOT `read_ampgrid`.) Method provided by `DataHDF5Ext`.
"""
function read_roman_nn end

# Untyped fallback (signature `(Any,)`): the HDF5 ext supplies the strictly-more-specific
# `(AbstractString,)` method, refining dispatch without an illegal precompile-time method-overwrite.
read_roman_nn(path; kwargs...) = error("read_roman_nn: " * _NEED_HDF5)

export read_roman_nn

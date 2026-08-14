# Kerr flux-grid loader. Per-region `(u,w,z) ∈ [0,1]³` tricubic axes over `FluxGridSchema`; the
# 4.74GB Kerr-ecc-eq grid is lazy/region-streamed via `_read_region`. HDF5 body in ext/DataHDF5Ext.jl.

"""
    read_fluxgrid(path, region::Symbol) -> (FluxGridSchema, channels)

Read one Kerr flux-grid region (`:A` close / `:B` far): the tricubic axes as a `FluxGridSchema` plus
the per-channel `length(u) × length(w) × length(z)` sample arrays. Method provided by `DataHDF5Ext`.
"""
function read_fluxgrid end

# Untyped `(Any,Any)` fallback: the HDF5 ext supplies the strictly-more-specific
# `(AbstractString,Symbol)` method, so loading HDF5 refines dispatch WITHOUT method-overwrite (illegal
# during ext precompilation). Without HDF5 this fires and points at the missing trigger.
read_fluxgrid(path, region; kwargs...) = error("read_fluxgrid: " * _NEED_HDF5)

export read_fluxgrid

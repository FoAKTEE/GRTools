"""
    Hybrid.Registry

Model registry (design §4.9): the upstream 7-plugin name→constructor catalogue
(`ListInspiralModels`/`ListAmplitudeModels`/`ListModels` + `GetInspiralEquations`/`GetAmplitudes`
dispatch), with the 7 inspiral constructors wired (the WaSABI-C / WaSABI-C_v0.9 entries dispatch to
the HybridC `wasabic_model` SF-PN-hybrid loader) and 4 amplitude constructors (the WaSABI-C
amplitude side remains the upstream-undefined contract).
"""
module Registry

using ..Models: forcing_0PASchwarz, forcing_1PAT1, forcing_0PAKerr, model_1PAT1ea, model_1PAT1R,
    amplitudes_0PASchwarz, amplitudes_1PAT1, amplitudes_0PAKerr, amplitudes_1PAT1ea
using ..HybridC: wasabic_model

include("Registry/registry.jl")

export list_inspiral_models, list_amplitude_models, list_models, get_inspiral, get_amplitudes

end # module Registry

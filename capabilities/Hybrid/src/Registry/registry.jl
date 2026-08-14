# Model registry: the name→constructor catalogue mirroring the upstream `ListInspiralModels` /
# `ListAmplitudeModels` / `ListModels` (which discover `FileBaseName /@ FileNames["*.m", ...]`
# from `InspiralModels/` and `AmplitudeModels/`) and the `GetInspiralEquations` / `GetAmplitudes`
# dispatch. The catalogue is the upstream 7-plugin set; the WaSABI-C / WaSABI-C_v0.9 inspiral
# entries are wired to the HybridC `wasabic_model` SF-PN-hybrid constructor (their amplitude
# entries are still the upstream-undefined contract, so `get_amplitudes` for them errors).

# Upstream plugin names (both `InspiralModels/` and `AmplitudeModels/` ship the same 7 plugins).
const _MODEL_CATALOG = sort(["0PAKerr", "0PASchwarz", "1PAT1", "1PAT1R", "1PAT1ea",
                             "WaSABI-C", "WaSABI-C_v0.9"])

# Name → constructor for the models ported so far. Catalogued names absent here are recognised
# but error as "not yet ported". The constructors adapt each plugin's on-disk data layout from a
# single base directory.
# The WaSABI-C SF-PN hybrid needs both the inspiral `sf_data` and the sibling `AmplitudeModels/
# sf_amp_data` (the (2,2) 1SF amplitude feeds the δx correction). The amplitude root is the
# sibling of the inspiral root: `…/InspiralModels/sf_data` → `…/AmplitudeModels/sf_amp_data`.
_wc_amp_dir(sf) = joinpath(dirname(dirname(sf)), "AmplitudeModels", "sf_amp_data")

const _INSPIRAL_CTORS = Dict{String,Function}(
    "0PASchwarz"    => sf -> forcing_0PASchwarz(joinpath(sf, "1SF_Flux", "Schwarz_Circ")),
    "1PAT1"         => sf -> forcing_1PAT1(sf),
    "0PAKerr"       => sf -> forcing_0PAKerr(sf),
    "1PAT1ea"       => sf -> model_1PAT1ea(sf),
    "1PAT1R"        => sf -> model_1PAT1R(sf),
    "WaSABI-C"      => sf -> wasabic_model(sf, _wc_amp_dir(sf)),
    "WaSABI-C_v0.9" => sf -> wasabic_model(sf, _wc_amp_dir(sf)),
)
const _AMPLITUDE_CTORS = Dict{String,Function}(
    "0PASchwarz" => amp -> amplitudes_0PASchwarz(joinpath(amp, "1SF", "Schwarz_Circ")),
    "1PAT1"      => amp -> amplitudes_1PAT1(amp),
    "0PAKerr"    => amp -> amplitudes_0PAKerr(amp),
    "1PAT1ea"    => amp -> amplitudes_1PAT1ea(amp),
)

"""
    list_inspiral_models() -> Vector{String}

The available inspiral-model plugin names (`ListInspiralModels`): `FileBaseName` of the upstream
`InspiralModels/*.m`.
"""
list_inspiral_models() = copy(_MODEL_CATALOG)

"""
    list_amplitude_models() -> Vector{String}

The available amplitude-model plugin names (`ListAmplitudeModels`): `FileBaseName` of the
upstream `AmplitudeModels/*.m`.
"""
list_amplitude_models() = copy(_MODEL_CATALOG)

"""
    list_models() -> Vector{String}

The available BinaryInspiral model names (`ListModels`); upstream aliases `ListInspiralModels`.
"""
list_models() = list_inspiral_models()

"""
    get_inspiral(name, sf_datadir)

Construct the inspiral forcing model `name` from the self-force data root `sf_datadir`
(`InspiralModels/sf_data`). Dispatches to the model constructor (`GetInspiralEquations`
analogue). Errors for an unknown name and for a catalogued-but-not-yet-ported model.
"""
function get_inspiral(name::AbstractString, sf_datadir::AbstractString)
    name in _MODEL_CATALOG || error("Unknown inspiral model $(repr(name)).")
    haskey(_INSPIRAL_CTORS, name) ||
        error("Inspiral model $(repr(name)) is catalogued but not yet ported.")
    return _INSPIRAL_CTORS[name](sf_datadir)
end

"""
    get_amplitudes(name, amp_datadir)

Construct the amplitude model `name` from the amplitude data root `amp_datadir`
(`AmplitudeModels/sf_amp_data`). Dispatches to the model constructor (`GetAmplitudes` analogue).
Errors for an unknown name and for a catalogued-but-not-yet-ported model.
"""
function get_amplitudes(name::AbstractString, amp_datadir::AbstractString)
    name in _MODEL_CATALOG || error("Unknown amplitude model $(repr(name)).")
    haskey(_AMPLITUDE_CTORS, name) ||
        error("Amplitude model $(repr(name)) is catalogued but not yet ported.")
    return _AMPLITUDE_CTORS[name](amp_datadir)
end

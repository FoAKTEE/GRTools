# Model introspection helpers.
# Port of ref-code/WaSABI/Kernel/{Utilities,Inspiral,Waveform}.m introspection surface.

"""
    list_modes(model, amp_datadir) -> Vector{Tuple{Int,Int}}

List available `(l, m)` modes in `model`, built from `amp_datadir`. Port of `ListModes`
(`Keys[GetAmplitudes[model]]`). Errors for an unknown / not-yet-ported model.
"""
list_modes(model::AbstractString, amp_datadir::AbstractString) =
    mode_list(get_amplitudes(model, amp_datadir))

# --- per-model forcing-channel ODE structure ------------------------------------------------
#
# Each inspiral model contributes a forcing-driven ODE system du/dt = F(u)/M. The symbolic
# structure (left-hand sides, algebraic constraints, state parameters) and the rhs evaluator are
# dispatched on the *constructed forcing object's* type — one source of truth shared by both
# `forcing_terms` and `get_inspiral_equations`. The per-model rhs reuses the already-grounded
# `omega_dot` (circular Schwarzschild), `forcing_Fω` (0PAKerr), `_spin_eom!` (1PAT1ea/R spin) and
# `wc_forcing` (WaSABI-C 9-var) primitives, so the assembly is grounded by self-consistency.

const _SCHWARZ_ODELHS = ["dΩ/dt", "dφ/dt", "dν/dt", "dM/dt"]
const _SCHWARZ_AELHS = ["Ω"]
const _SCHWARZ_PARAMETERS = ["r0", "φ", "ν", "M"]
const _KERR_ODELHS = ["dω/dt", "dφ/dt"]
const _KERR_PARAMETERS = ["ω", "φ", "ν", "M", "χ1"]
const _SPIN_ODELHS = ["dΩ/dt", "dφ/dt", "dν/dt", "dM/dt", "dχt1/dt", "dχt2/dt", "dδm/dt"]
const _SPIN_PARAMETERS = ["Ω", "φ", "ν", "M", "χt1", "χt2", "δm"]
const _WC_ODELHS = ["dω/dt", "dφ/dt", "dν/dt", "dM/dt", "dχ1/dt", "dχ2/dt",
                    "dδm/dt", "dδν/dt", "dδχ/dt"]
const _WC_PARAMETERS = ["ω", "φ", "ν", "M", "χ1", "χ2", "δm", "δν", "δχ"]

# (odelhs, aelhs, parameters) for the constructed forcing object.
_forcing_labels(::Union{Forcing0PASchwarz,Forcing1PAT1}) =
    (copy(_SCHWARZ_ODELHS), copy(_SCHWARZ_AELHS), copy(_SCHWARZ_PARAMETERS))
_forcing_labels(::Forcing0PAKerr) = (copy(_KERR_ODELHS), String[], copy(_KERR_PARAMETERS))
_forcing_labels(::Union{Model1PAT1ea,Model1PAT1R}) =
    (copy(_SPIN_ODELHS), String[], copy(_SPIN_PARAMETERS))
_forcing_labels(::WaSABICModel) = (copy(_WC_ODELHS), String[], copy(_WC_PARAMETERS))

_pv(pv, k, d) = haskey(pv, k) ? float(pv[k]) : d

# rhs du/dt evaluators (Dict keyed by the model's ODElhs).
function _eval_forcing(forcing::Union{Forcing0PASchwarz,Forcing1PAT1}, pv::AbstractDict)
    r0 = float(pv["r0"]); ν = float(pv["ν"]); M = float(pv["M"])
    return Dict("dΩ/dt" => omega_dot(forcing, r0, ν, M),
                "dφ/dt" => kepler_Ω(r0) / M,   # algebraic constraint Ω = √(1/r0³)
                "dν/dt" => 0.0, "dM/dt" => 0.0)
end

function _eval_forcing(forcing::Forcing0PAKerr, pv::AbstractDict)
    ω = float(pv["ω"]); ν = float(pv["ν"]); χ1 = float(pv["χ1"]); M = _pv(pv, "M", 1.0)
    return Dict("dω/dt" => forcing_Fω(forcing, ω, ν, χ1) / M, "dφ/dt" => ω / M)
end

function _eval_forcing(forcing::Union{Model1PAT1ea,Model1PAT1R}, pv::AbstractDict)
    u = [float(pv["Ω"]), _pv(pv, "φ", 0.0), float(pv["ν"]), _pv(pv, "M", 1.0),
         float(pv["χt1"]), float(pv["χt2"]), _pv(pv, "δm", 0.0)]
    du = similar(u)
    _spin_eom!(du, u, forcing)
    return Dict(_SPIN_ODELHS[i] => du[i] for i in eachindex(_SPIN_ODELHS))
end

function _eval_forcing(forcing::WaSABICModel, pv::AbstractDict)
    M = _pv(pv, "M", 1.0)
    F = wc_forcing(forcing, float(pv["ω"]), _pv(pv, "φ", 0.0), float(pv["ν"]), 1,
                   float(pv["χ1"]), float(pv["χ2"]),
                   _pv(pv, "δm", 0.0), _pv(pv, "δν", 0.0), _pv(pv, "δχ", 0.0))
    # u' = (Fω/M, ω/M, Fν/M, Fm, Fχ1/M, Fχ2/M, Fδm/M, Fδν/M, Fδχ/M) (WaSABI-C.m §Evolution).
    rates = (F[1] / M, F[2] / M, F[3] / M, F[4], F[5] / M, F[6] / M, F[7] / M, F[8] / M, F[9] / M)
    return Dict(_WC_ODELHS[i] => rates[i] for i in eachindex(_WC_ODELHS))
end

"""
    ForcingTermsModel

The forcing terms / ODE of a model: the constructed forcing object plus the symbolic structure
(`odelhs`, `aelhs`, `parameters`). Build with [`forcing_terms`](@ref); evaluate the forcing-driven
rhs at a point in parameter space by calling it on a `Dict` of parameter values. Port of the
`ForcingTermsModel[X]` head.
"""
struct ForcingTermsModel{F}
    model::String
    forcing::F
    odelhs::Vector{String}
    aelhs::Vector{String}
    parameters::Vector{String}
end

"""
    forcing_terms(model, sf_datadir) -> ForcingTermsModel

Generate the forcing-driven ODE structure of `model` from the self-force data root `sf_datadir`.
Port of `ForcingTerms`, extended beyond the Schwarzschild-circular split to the 0PAKerr (Fω),
spin (1PAT1ea/R) and WaSABI-C (9-channel Fω/Fν/Fm/Fχ) models. Errors only for an unknown name.
"""
function forcing_terms(model::AbstractString, sf_datadir::AbstractString)
    model in list_inspiral_models() || error("Unknown model $(repr(model)).")
    forcing = get_inspiral(model, sf_datadir)
    odelhs, aelhs, params = _forcing_labels(forcing)
    return ForcingTermsModel(String(model), forcing, odelhs, aelhs, params)
end

"""
    (ftm::ForcingTermsModel)(paramvals) -> Dict{String,Float64}

Evaluate the forcing-driven rhs `du/dt = F(u)/M` at `paramvals` (a `Dict` keyed by
`ftm.parameters`), dispatching on the model's forcing object. For the circular-Schwarzschild
models resolves the algebraic constraint `Ω = √(1/r0³)` (`Spacetime.kepler_Ω`). Port of the
`ForcingTermsModel[assoc][paramvals]` evaluator (Utilities.m:165-174).
"""
(ftm::ForcingTermsModel)(paramvals::AbstractDict) = _eval_forcing(ftm.forcing, paramvals)

# --- inspiral equations (filled H-WASABI-STUBS) ---------------------------------------------

"""
    InspiralEquations

The coupled forcing-driven ODE system of an inspiral model — the modern `GetInspiralEquations`
return. A [`ForcingTermsModel`](@ref) alias-by-composition: it *is* the model's forcing-terms ODE
structure (AHA — the inspiral equations are exactly the forcing-driven rhs `du/dt = F(u)/M`).
"""
const InspiralEquations = ForcingTermsModel

"""
    get_inspiral_equations(model, sf_datadir) -> ForcingTermsModel

Fetch the coupled equations describing the inspiral evolution of `model`. The forcing-driven ODE
system `du/dt = F(u)/M` assembled from the model's forcing channels (`forcing_terms`); evaluate
the rhs at a state by calling the result on a parameter `Dict` (or via [`inspiral_evaluate`](@ref)).
Filled functional body for the upstream declaration-only `GetInspiralEquations` (H-WASABI-STUBS):
the rhs is the already-grounded `omega_dot`/`forcing_Fω`/`_spin_eom!`/`wc_forcing` assembly.
"""
get_inspiral_equations(model::AbstractString, sf_datadir::AbstractString) =
    forcing_terms(model, sf_datadir)

"""
    inspiral_evaluate(eqs::ForcingTermsModel, paramvals) -> Dict{String,Float64}
    inspiral_evaluate(solution, t) -> NamedTuple

Evaluate an inspiral model along a trajectory. With a [`ForcingTermsModel`](@ref) /
[`InspiralEquations`](@ref) and a parameter `Dict`, returns the forcing-driven rhs `du/dt` at that
state. With an integrated `solution` and a time `t`, returns the state along the trajectory
(`HybridC.wc_state_at` for a WaSABI-C run, `(; ω, φ)` from the `Drivers` `*_at` accessors for an
N-var run). Filled functional body for the upstream declaration-only `InspiralEvaluate`
(H-WASABI-STUBS).
"""
inspiral_evaluate(eqs::ForcingTermsModel, paramvals::AbstractDict) = eqs(paramvals)
inspiral_evaluate(s::WaSABICSolution, t::Real) = wc_state_at(s, t)
inspiral_evaluate(s::Union{InspiralSolution,SpinInspiralSolution}, t::Real) =
    (ω = omega_at(s, t), φ = phi_at(s, t))

# --- citation metadata ----------------------------------------------------------------------

"""
    cite_model(model)

Return the bib keys to cite when using `model`. Port of `CiteModel`. Reads the verbatim
`"Citation"` string from the shipped `assets/MetadataModels/<model>.m` asset (a Mathematica
association). An unknown `model` errors, mirroring `CiteModel::nomodel`.
"""
function cite_model end

# Directory of the shipped MetadataModels assets (`<|"Citation"->...|>` files).
_metadata_dir() = joinpath(@__DIR__, "..", "..", "assets", "MetadataModels")

"""
    list_metadata_models() -> Vector{String}

Known citation models == `FileBaseName` of `MetadataModels/*.m`. Port of the model list built
inside `CiteModel` (Utilities.m:184-188).
"""
function list_metadata_models()
    dir = _metadata_dir()
    sort([first(splitext(f)) for f in readdir(dir) if endswith(f, ".m")])
end

# Extract the `"Citation"` value from a Mathematica association source, resolving the `\"` / `\\`
# escapes to their literal characters.
function _extract_citation(text::AbstractString)
    m = match(r"\"Citation\"\s*->\s*\"((?:[^\"\\]|\\.)*)\""s, text)
    m === nothing && error("No Citation entry found.")
    raw = m.captures[1]
    return replace(raw, "\\\\" => "\\", "\\\"" => "\"")
end

function cite_model(model::AbstractString)
    model in list_metadata_models() ||
        error("No citation guideline found for model $(repr(model)).")
    return _extract_citation(read(joinpath(_metadata_dir(), model * ".m"), String))
end

# High-level binary-inspiral assembly built from initial conditions.
# Port of ref-code/WaSABI/Kernel/BinaryInspiral.m (WaSABI`BinaryInspiral` context).

# Initial-conditions formats per model (`InitialConditionsFormat`): the parameters a model needs
# to evaluate `binary_inspiral`. Kept here (rather than Introspection) so `binary_inspiral` can
# validate against it without a submodule cycle; `Introspection.initial_conditions` re-exports.
const _ICS_FORMATS = Dict("0PASchwarz" => ["M", "r0", "ν", "φ"],
                          "1PAT1"      => ["M", "r0", "ν", "φ"],
                          "0PAKerr"    => ["ω", "φ", "ν", "M", "χ1"])

"""
    initial_conditions(model) -> Vector{String}

The parameters in `model` that require values to evaluate `binary_inspiral` or `forcing_terms`
(`InitialConditionsFormat`). Port of `InitialConditions`. Errors for an unknown / not-yet-ported
model.
"""
function initial_conditions(model::AbstractString)
    model in list_inspiral_models() || error("Unknown model $(repr(model)).")
    haskey(_ICS_FORMATS, model) ||
        error("Model $(repr(model)) is catalogued but not yet ported.")
    return copy(_ICS_FORMATS[model])
end

"""
    BinaryInspiralModel

A binary inspiral: the integrated trajectory plus the strain-amplitude model, built by
[`binary_inspiral`](@ref). Query the strain modes with [`waveform`](@ref), the orbital
trajectory with [`trajectory`](@ref), and the available modes with [`mode_list`](@ref).
Port of the upstream `BinaryInspiralModel[<|...|>]` head and its `"Waveform"` / `"Trajectory"`
accessors (`BinaryInspiral.m`).
"""
struct BinaryInspiralModel{A}
    model::String
    ics::Dict{String,Float64}
    inspiral::InspiralSolution
    amp::A
    duration::Float64
end

mode_list(bim::BinaryInspiralModel) = mode_list(bim.amp)

"""
    binary_inspiral(ics; model="1PAT1", sf_datadir, amp_datadir, kwargs...) -> BinaryInspiralModel

Generate a [`BinaryInspiralModel`](@ref) from the initial conditions `ics` (a `Dict` whose keys
match `initial_conditions(model)`). Integrates the registry-dispatched forcing model (built from
`sf_datadir`) and pairs it with the strain-amplitude model (built from `amp_datadir`). Extra
`kwargs` (`reltol`, `abstol`, …) forward to [`int_inspiral`](@ref). Errors (mirroring
`BinaryInspiral::ics`) when the supplied keys do not match the model's required set.
"""
function binary_inspiral(ics::AbstractDict; model::AbstractString = "1PAT1",
                         sf_datadir::AbstractString, amp_datadir::AbstractString, kwargs...)
    Set(string.(keys(ics))) == Set(initial_conditions(model)) ||
        error("Invalid initial conditions $(sort(string.(collect(keys(ics))))) for model " *
              "$(repr(model)); expected keys $(sort(initial_conditions(model))).")
    forcing = get_inspiral(model, sf_datadir)
    amp = get_amplitudes(model, amp_datadir)
    icsf = Dict{String,Float64}(string(k) => float(v) for (k, v) in ics)
    sol = _integrate(forcing, icsf; kwargs...)
    return BinaryInspiralModel(String(model), icsf, sol, amp, duration(sol))
end

# Model-specific dispatch from the initial-conditions dict to the right `int_inspiral` driver.
_integrate(forcing::Union{Forcing0PASchwarz,Forcing1PAT1}, ics; kwargs...) =
    int_inspiral(forcing, ics["r0"], ics["ν"], ics["M"]; φ0 = ics["φ"], kwargs...)
_integrate(forcing::Forcing0PAKerr, ics; kwargs...) =
    int_inspiral(forcing, ics["ω"], ics["ν"], ics["M"], ics["χ1"]; φ0 = ics["φ"], kwargs...)

"""
    duration(bim::BinaryInspiralModel) -> Float64

Inspiral duration `tmax` (`Max[..."Domain"]` upstream).
"""
duration(bim::BinaryInspiralModel) = bim.duration

# Strain amplitude h_{l,m} of the model's amplitude object, evaluated at the trajectory frequency
# (Ω for Schwarz, ω for Kerr). The stored frequency-domain amplitude already absorbs the model's
# ν / spin dependence (`BinaryInspiral.m:168`).
_mode_amplitude(amp::Amp0PASchwarz, l, m, Ω, ics) = h0PA_ampΩ(amp, l, m, Ω, ics["ν"])
_mode_amplitude(amp::Amp1PAT1, l, m, Ω, ics) = h1PA_ampΩ(amp, l, m, Ω, ics["ν"])
_mode_amplitude(amp::Amp0PAKerr, l, m, ω, ics) = real_amp22(amp, ω, ics["χ1"], ics["ν"], ics["M"])

# Out-of-domain test mirroring `outsideDomainQ` (BinaryInspiral.m:147): t<0 or t>tmax.
_outside_domain(t, tmax) = t < 0 || t > tmax

"""
    waveform(bim::BinaryInspiralModel, l, m, t) -> ComplexF64

Spin-weighted spherical-harmonic strain mode `h_{l,m}(t) = A_{l,m}(Ω(t)) exp(-i m φ(t))`
evaluated along the inspiral trajectory. Negative `m` uses the equatorial symmetry
`h_{l,-m} = (-1)^l conj(h_{l,m})` (`BinaryInspiral.m:172`). An unavailable `(l, m)` errors
(`BinaryInspiralModel::nomode`); a time outside `[0, duration]` errors
(`BinaryInspiralModel::endinsp`).
"""
function waveform(bim::BinaryInspiralModel, l::Integer, m::Integer, t::Real)
    if m < 0
        return (-1)^l * conj(waveform(bim, l, -m, t))
    end
    (l, m) in mode_list(bim) ||
        error("Mode $((l, m)) not available in model $(repr(bim.model)).")
    _outside_domain(t, bim.duration) &&
        error("Time $t exceeds inspiral duration $(bim.duration).")
    Ω = omega_at(bim.inspiral, t)
    φ = phi_at(bim.inspiral, t)
    return _mode_amplitude(bim.amp, l, m, Ω, bim.ics) * exp(-im * m * φ)
end

"""
    trajectory(bim::BinaryInspiralModel, param, t) -> Float64

Orbital trajectory parameter `param ∈ {"Ω","ω","r0","φ","ν","M","χ1"}` at time `t`. `"Ω"`/`"ω"`,
`"r0"`, `"φ"` read the integrated solution; the rest are the (constant) initial conditions. A
time outside `[0, duration]` errors. Port of the `"Trajectory"` accessor (`BinaryInspiral.m:179`).
"""
function trajectory(bim::BinaryInspiralModel, param::AbstractString, t::Real)
    _outside_domain(t, bim.duration) &&
        error("Time $t exceeds inspiral duration $(bim.duration).")
    if param == "Ω" || param == "ω"   # Schwarz frequency Ω, Kerr frequency ω
        return omega_at(bim.inspiral, t)
    elseif param == "r0"
        return r0_at(bim.inspiral, t)
    elseif param == "φ"
        return phi_at(bim.inspiral, t)
    elseif haskey(bim.ics, param)      # constant parameters (ν, M, χ1)
        return bim.ics[param]
    else
        error("Unknown trajectory parameter $(repr(param)).")
    end
end

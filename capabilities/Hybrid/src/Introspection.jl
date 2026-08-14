"""
    Hybrid.Introspection

Model introspection (design §4.9): `cite_model` (citation metadata from the shipped
`MetadataModels` assets), `list_modes`, `forcing_terms` (the per-model forcing-channel ODE
structure + its rhs evaluator), `initial_conditions`, and the now-filled `get_inspiral_equations`
/ `inspiral_evaluate` (the forcing-driven inspiral ODE system + its evaluation along a trajectory).
The forcing/ODE rhs are assembled from — and grounded by self-consistency against — the already
-grounded `Models` (`omega_dot`/`forcing_Fω`/`_spin_eom!`) and `HybridC` (`wc_forcing`/
`wc_state_at`) primitives.
"""
module Introspection

using ..Registry: get_inspiral, get_amplitudes, list_inspiral_models
using ..Models: omega_dot, forcing_Fω, Forcing0PASchwarz, Forcing1PAT1, Forcing0PAKerr,
    Model1PAT1ea, Model1PAT1R, _spin_eom!
using ..HybridC: WaSABICModel, wc_forcing, WaSABICSolution, wc_state_at
using ..Drivers: InspiralSolution, SpinInspiralSolution, omega_at, phi_at
using ..BinaryWaveform: mode_list, initial_conditions
using Spacetime: kepler_Ω

include("Introspection/introspection.jl")

export cite_model, list_metadata_models, list_modes, forcing_terms, ForcingTermsModel,
    InspiralEquations, initial_conditions, get_inspiral_equations, inspiral_evaluate

end # module Introspection
